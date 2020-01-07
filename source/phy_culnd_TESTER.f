      
C      SUBROUTINE CULND(WTWB,WTWG,WTRUN,WTDRN,PRECNV,PRECLS,
C     &                 EVAP,RLD,TG1)

C      SUBROUTINE CULND(WTWB,WTWG,WTRUN,WTDRN,PRECNV,PRECLS,
C     &                 EVAP,RLD,TSKIN)

      SUBROUTINE CULND(WTWB,WTWG,WTRUN,WTDRN,PRECNV,PRECLS,
     &                 EVAP)
C
C Subroutine CULND
C   Updates the simple 2 bucket land model of Deardorff (1977), 
C   adapted for water isotope tracers.
C   This is an updated version of the scheme in MUGCM (see Noone and
C   Simmonds, 2002), and code has been formatted to be SPEEDY-fied.
C
C   [Notice this also works over the ocean. Could form the basis of a
C   simple ocean mixed layer model, but would need a deep ocean
C   relaxation term, rather than draining to empty.]
C
C
C   THIS VERSION IS NON-FRACTIONATING. Could be (only EVAP to top layer)!
C
C David Noone <dcn@colorado.edu> - Wed Mar 28 16:42:27 MDT 2012
C
C--   Tracers are assumed:
C--   ITR=1     Q	: normal water vapor (mass) mixing ratio

C     Define all variables
C     Input-only
C       PREC			total precipitation
C       EVAP                    total evaporation
C

C ! Need to pass these back through SUFLUX, and there we will calculate the isotope ratio of evap. 
C ! Next look through we'll pass it the isotopic composition of the soil water, and then we'll pass that to the evaporation.
C ! First we need evaporation from SUFLUX.

C     Input/output
C       WTWB                    bulk water		[m^3/m^3]
C       WTWG                    ground water 		[m^3/m^3]
C       RLD			land isotope ratio

C     Output
C       WTRUN                   runoff			[g/m^2 s]
C       WTDRN                   subsurface drainage	[g/m^2 s]
C

C     Other variables
C      WBXS			bulk excess 		[m^3/m^3]
C      WGXS			ground water excess 	[m^3/m^3]


C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude

      include "com_physcon.h"

C     Large-scale condensation constants

      include "com_lsccon.h"

C     Isotope tracer constants

      include "com_isocon.h"

C     Time parameters
      include "com_tsteps.h"

C 
C     External isotope fractionation functions
C

      REAL ALPLIQ, ALPKOC
      REAL D

C
C Input variables (commented out if in common block PHYSVAR)
C
      REAL EVAP(NGP,1,NTR)
      REAL PREC(NGP,NTR)
      REAL PRECLS(NGP,NTR)
      REAL PRECNV(NGP,NTR)
C      REAL TSKIN(NGP)
C      REAL TG1(NGP,NLEV)

C
C Input/output variables (Tracer soil water fields)
C 
      REAL WTWB(NGP,NTR)			
      REAL WTWG(NGP,NTR)
      REAL WTRUN(NGP,NTR)
      REAL WTDRN(NGP,NTR)


C
C Working variables
C
      REAL WBXS(NGP,NTR)
      REAL WGXS(NGP,NTR)
      
      INTEGER ITR, J
      REAL RLD(NGP,NTR)
C
C Parameters
C

      REAL wgmax , wbmax
      REAL d1    , d2
      REAL wtiny 
      REAL tau1  
      REAL rhow  

      parameter(wgmax = 0.4)              ! ground max volumetic content (m^3/m^3)
      parameter(wbmax = 0.32)             ! bulk max volumetic content (m^3/m^3)
      parameter(d1    = 0.05)             ! depth of upper layer (meters)
      parameter(d2    = 0.15)             ! depth of both layers (meters)
      parameter(wtiny = 0.001)            ! trivial amount to retain ratio
      parameter(tau1  = 86400.)           ! upper damping scale (seconds)
      parameter(rhow  = 1.0e+06)          ! density of water (g/m^3)

C 0.0 Establish molecular diffusion rate ratios
C Establish (D'/D, Merlivat 1977): D'/D = 0.9723 for o18 and 0.9755 for dD.

           if (ITR .EQ. 1) then 
                D = 1.0
      
           else if (ITR .EQ. 2) then 
                D = 1.0
      
           else if (ITR .EQ. 3) then 
                D = 0.9755
      
           else if (ITR .EQ. 4) then 
                D = 0.9723
           
           endif
            
           D = D**0.58

C 0) Compute total Precipitation at point, Convert units of precipitation to account for density
C rhow = divide by 1000 to get kg/m^2/s
      
      DO ITR = 2, NTR
        DO J = 1, NGP
          PREC(J,ITR) = PRECNV(J,ITR) + PRECLS(J,ITR)
        ENDDO
      ENDDO

C PRECLS = large-scale precipitation 	[g/(m^2 s)]   	(2-dim)
C PRECNV = convective precipitation 	[g/(m^2 s)]     (2-dim)
C EVAP   = evaporation 			[g/(m^2 s)]     (2-dim)

C

C 1) Prescribe ratio for fluxes based on surface scheme
C  (this is actually quite disatisfactory, as the "wetness" does not
C  affect cflx...so buckets can drain dry, and still have evap)
C  (THIS SHOULD BE MOVED TO SUFLX)


      DO ITR = 2, NTR
        DO J = 1, NGP

            if (wtwg(J,1) .gt. wtiny) then
               RLD(J,ITR) = (wtwg(J,ITR)/wtwg(J,1))
            else if (wtwb(J,1) > 0.99*wtiny) then 
               RLD(J,ITR) = (wtwb(J,ITR)/wtwb(J,1))
            else                                        	! buckets are dry, assign smow?
               RLD(J,ITR) = 1.0  			! SDEE set equal to 1, because multiply by Rld is already done in SUFLUX.
C               RLD(J,ITR) = 0.0
            end if
            
C            write (*,*) "RLD = ", RLD(J,ITR)

C            EVAP(J,1,ITR)  = RLD(J,ITR)*EVAP(J,1,1)		
C            EVAP(J,1,ITR)= RLD(J,ITR)*EVAP(J,1,1)*D/
C     &                     alpliq(ITR,TSKIN(J)) 

C            EVAP(J,1,ITR)= RLD(J,ITR)*EVAP(J,1,1)*D/
C     &                     alpliq(ITR,TG1(J,NLEV)) 
            
C            write(*,*) 'EVAP alphas',alpliq(ITR,TG1(J,NLEV)),RLD

C Then we'll have two evaporations and weight them. EVAP = VEVAP*0.66 * SEVAP*0.33 (veg and soil)
C SDEE (UPDATE NEEDED) CALL ALPLIQ AND ALPKEX HERE???
C Will come from plans and from soils. Only need K-frac if it's coming from soils. Okay for now. 
            
        ENDDO
      ENDDO
C


C
C 2) Update surface reservoirs:  dbucket = P - E - drain - Runoff
C
      DO J = 1, NGP
        dwb = (PREC(J,1)-EVAP(J,1,1))/(d2*rhow)  			! change in water/change time 
        dwg = (PREC(J,1)-EVAP(J,1,1))/(d1*rhow)
!
        dwg = dwg - (wtwg(J,1)-wtwb(J,1))/tau1  		! drain/recharge
!
        wtwg(J,1) = max(wtwg(J,1) + DELT*dwg, 0.)   		
        wgxs(J,1)  = max(wtwg(J,1) - wgmax, 0.)   		! max, 0, make sure water amount is greater than 0.
        wtwg(J,1) = wtwg(J,1) - wgxs(J,1)			! subtract out excess
!
        wtwb(J,1) = max(wtwb(J,1) + DELT*dwb, 0.)
        wbxs(J,1)  = max(wtwb(J,1) - wbmax, 0.)
        wtwb(J,1) = wtwb(J,1) - wbxs(J,1)			! subtract out excess
!
        wtrun(J,1) = d1*rhow*wgxs(J,1)/DELT				!runoff = (depth of layer*excess)/timestep
        wtdrn(J,1) = d2*rhow*wbxs(J,1)/DELT
!
! If ground bucket has dried up completely, add a tiny mass so that
! we can save the long term isotope ratio.
!
        wtwb(J,1) = max(wtwb(J,1), wtiny)

      ENDDO


C
C Apply same water balance, but for tracers
C
      DO ITR = 1, NTR
        DO J = 1, NGP
!
            dwb = (PREC(J,ITR) -EVAP(J,1,ITR))/(d2*rhow)
            dwg = (PREC(J,ITR) -EVAP(J,1,ITR))/(d1*rhow)
            dwg = dwg - (wtwg(J,ITR)-wtwb(J,ITR))/tau1    ! drain/recharge
!
            wtwg(J,ITR) = max(wtwg(J,ITR) + DELT*dwg, 0.)
            wgxs(J,ITR)  = wgxs(J,1)*RLD(J,ITR)
            wtwg(J,ITR) = wtwg(J,ITR) - wgxs(J,ITR)
!
            wtwb(J,ITR) = max(wtwb(J,ITR) + DELT*dwb, 0.)
            wbxs(J,ITR)  = wbxs(J,1)*RLD(J,ITR)
            wtwb(J,ITR) = wtwb(J,ITR) - wbxs(J,ITR)
!
! If total has dried up apply the saved isotope ratio so that we have
! something sensible to evaporate next time
!
            if (wtwb(J,1) .le. wtiny) then
              wtwb(J,ITR) = RLD(J,ITR)*wtwb(J,1)
            end if
!
            wtrun(J,ITR) = d1*rhow*wgxs(J,ITR)/DELT
            wtdrn(J,ITR) = d2*rhow*wbxs(J,ITR)/DELT
        ENDDO
      ENDDO
C
      RETURN
      END       



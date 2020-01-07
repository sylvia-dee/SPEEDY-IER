
      SUBROUTINE SUFLUX (PSA,UA,VA,TA,QA,RH,PHI,
     &                   PHI0,FMASK,TLAND,TSEA,SWAV,SSR,SLRD,
     &                   USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
     &                   TSFC,TSKIN,U0,V0,T0,Q0,LFLUXLAND,RLD,
     &                   RB,SOILW)
C--
C--   SUBROUTINE SUFLUX (PSA,UA,VA,TA,QA,RH,PHI,
C--  &                   PHI0,FMASK,TLAND,TSEA,SWAV,SSR,SLRD,
C--  &                   USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
C--  &                   TSFC,TSKIN,U0,V0,T0,Q0,LFLUXLAND)
C--
C--   Purpose: Compute surface fluxes of momentum, energy and moisture,
C--            and define surface skin temperature from energy balance
C--   Input:   PSA    = norm. surface pressure [p/p0]   (2-dim)
C--            UA     = u-wind                          (3-dim)
C--            VA     = v-wind                          (3-dim)
C--            TA     = temperature                     (3-dim)
C--            QA     = specific humidity [g/kg]        (3-dim)
C--            RH     = relative humidity [0-1]         (3-dim)
C--            PHI    = geopotential                    (3-dim)
C--            PHI0   = surface geopotential            (2-dim)
C--            FMASK  = fractional land-sea mask        (2-dim)
C--            TLAND  = land-surface temperature        (2-dim)
C--            TSEA   =  sea-surface temperature        (2-dim)
C--            SWAV   = soil wetness availability [0-1] (2-dim)
C--            SSR    = sfc sw radiation (net flux)     (2-dim)
C--            SLRD   = sfc lw radiation (downward flux)(2-dim)
C--            LFLUXLAND   = Logical related ti flux-correction
C--   Output:  USTR   = u stress                        (2-dim)
C--            VSTR   = v stress                        (2-dim)
C--            SHF    = sensible heat flux              (2-dim)
C--            EVAP   = evaporation [g/(m^2 s)]         (2-dim)
C--            SLRU   = sfc lw radiation (upward flux)  (2-dim)
C--            HFLUXN = net heat flux into land/sea     (2-dim)           
C--            TSFC   = surface temperature (clim.)     (2-dim)
C--            TSKIN  = skin surface temperature        (2-dim)
C--            U0     = near-surface u-wind             (2-dim)
C--            V0     = near-surface v-wind             (2-dim)
C--            T0     = near-surface air temperature    (2-dim)
C--            Q0     = near-surface sp. humidity [g/kg](2-dim)
C--
C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude

      include "com_physcon.h"

C     Surface flux constants

      include "com_sflcon.h"      

      include "com_radcon.h"

      include "com_isocon.h"

C     Isotope fractionation functions
      REAL ALPLIQ, ALPKOC

      REAL PSA(NGP), UA(NGP,NLEV), VA(NGP,NLEV), TA(NGP,NLEV),
     &     QA(NGP,NLEV,NTR), RH(NGP,NLEV), PHI(NGP,NLEV),
     &     PHI0(NGP), FMASK(NGP), TLAND(NGP), TSEA(NGP), SWAV(NGP),
     &     SSR(NGP), SLRD(NGP)

      REAL USTR(NGP,3), VSTR(NGP,3), SHF(NGP,3), EVAP(NGP,3,NTR),
     &     SLRU(NGP,3), HFLUXN(NGP,2), TSFC(NGP), TSKIN(NGP),
     &     U0(NGP), V0(NGP), T0(NGP), Q0(NGP,NTR), RLD(NGP,NTR),
     &     SEVAP(NGP,3,NTR), VEVAP(NGP,3,NTR), RB(NGP,NTR), SOILW(NGP),
     &     RVAPOR(NGP,NTR), REVAP(NGP,NTR)
									
      REAL T1(NGP,2), T2(NGP,2), QSAT0(NGP,2), 
     &     DENVVS(NGP,0:2), DENS(NGP,NLEV), DSLR(NGP), DTSKIN(NGP)

      LOGICAL LSCASYM, LSCDRAG, LSKINEB, BUCKET, CRAIGG

      LOGICAL LFLUXLAND
      real VMAG

      SAVE DENVVS

      LSCASYM = .true.   ! true : use an asymmetric stability coefficient
      LSCDRAG = .true.   ! true : use stability coef. to compute drag over sea
      LSKINEB = .true.   ! true : redefine skin temp. from energy balance
      BUCKET  = .false.  ! true : use soil moisture ratio from CULND to define soil wetness factor
      CRAIGG  = .true.   ! true: use Craig-Gordon equations to calculate isotope ratio of evaporative flux

      CLAMBDA = 7.        ! Heat conductivity in skin layer
 
      IF ( LFLUXLAND )  THEN

C--   1. Extrapolation of wind, temp, hum. and density to the surface

C     1.1 Wind components
   
      DO J=1,NGP
        U0(J) = FWIND0*UA(J,NLEV)
        V0(J) = FWIND0*VA(J,NLEV)
      ENDDO

C     1.2 Temperature

      GTEMP0 = 1.-FTEMP0
      RCP = 1./CP
      RDPHI0 =-1./(RD*288.*SIGL(NLEV))
      NL1=NLEV-1
C
      DO J=1,NGP
c       Temperature difference between lowest level and sfc
        DT1 = WVI(NLEV,2)*(TA(J,NLEV)-TA(J,NL1))
c       Extrapolated temperature using actual lapse rate
        T1(J,1) = TA(J,NLEV)+DT1
        T1(J,2) = T1(J,1)+PHI0(J)*DT1*RDPHI0
c       Extrapolated temperature using dry-adiab. lapse rate
        T2(J,2) = TA(J,NLEV)+RCP*PHI(J,NLEV)
        T2(J,1) = T2(J,2)-RCP*PHI0(J)
      ENDDO

      DO J=1,NGP
        IF (TA(J,NLEV).GT.TA(J,NL1)) THEN
          T1(J,1) = FTEMP0*T1(J,1)+GTEMP0*T2(J,1)
          T1(J,2) = FTEMP0*T1(J,2)+GTEMP0*T2(J,2)
        ELSE
          T1(J,1) = TA(J,NLEV)
          T1(J,2) = TA(J,NLEV)
        ENDIF
        T0(J) = T1(J,2)+FMASK(J)*(T1(J,1)-T1(J,2))
      ENDDO

C     1.3 Spec. humidity

      GHUM0 = 1.-FHUM0

      CALL SHTORH (-1,NGP,T0,PSA,1.,Q0(1,1),RH(1,NLEV),QSAT0)

      DO ITR=1,NTR
        DO J=1,NGP
          Q0(J,ITR)=FHUM0*Q0(J,ITR)+GHUM0*QA(J,NLEV,ITR)
        ENDDO
      ENDDO

C     1.4 Density * wind speed (including gustiness factor)

      PRD = P0/RD
      VG2 = VGUST*VGUST

      DO J=1,NGP
        DENVVS(J,0)=(PRD*PSA(J)/T0(J))*
     &              SQRT(U0(J)*U0(J)+V0(J)*V0(J)+VG2)

        DENS(J,NLEV)=(PRD*PSA(J)/T0(J))

C sdee addded DENS to get density of bottom layer--need for alphk to calculate v

      ENDDO

C     A) LAND FLUXES

C     1.5 Define effective skin temperature to compensate for
C         non-linearity of heat/moisture fluxes during the daily cycle

      DO JLAT=1,NLAT
	J0=NLON*(JLAT-1)
        SQCLAT=SQRT(CLAT(JLAT))
        DO J=J0+1,J0+NLON
          TSKIN(J)=TLAND(J)+CTDAY*SQCLAT*SSR(J)*PSA(J)
        ENDDO
      ENDDO

C     1.6 Stability correction

      RDTH  = FSTAB/DTHETA
      ASTAB = 1.
      IF (LSCASYM) ASTAB = 0.5   ! to get smaller dS/dT in stable conditions

      DO J=1,NGP
        IF (TSKIN(J).GT.T2(J,1)) THEN
           DTHL=MIN(DTHETA,TSKIN(J)-T2(J,1))
        ELSE
           DTHL=MAX(-DTHETA,ASTAB*(TSKIN(J)-T2(J,1)))
        ENDIF
        DENVVS(J,1)=DENVVS(J,0)*(1.+DTHL*RDTH)
      ENDDO

C--   2. Computation of fluxes over land and sea

C     2.1 Wind stress

      DO J=1,NGP

        CDLDV = CDL*DENVVS(J,0)*FOROG(J)

        USTR(J,1) = -CDLDV*UA(J,NLEV)
        VSTR(J,1) = -CDLDV*VA(J,NLEV)

      ENDDO

C     2.2 Sensible heat flux 

      CHLCP = CHL*CP

      DO J=1,NGP
        SHF(J,1) = CHLCP*DENVVS(J,1)*(TSKIN(J)-T1(J,1))
      ENDDO

C     2.3 Evaporation

      CALL SHTORH (0,NGP,TSKIN,PSA,1.,QDUMMY,RDUMMY,QSAT0(1,1))

      IF (BUCKET) THEN
      DO J=1,NGP
      EVAP(J,1,1)=CHL*DENVVS(J,1)*MAX(0.,SOILW(J)*QSAT0(J,1)-Q0(J,1))
      ENDDO

      ELSE
      DO J=1,NGP
        EVAP(J,1,1)=CHL*DENVVS(J,1)*MAX(0.,SWAV(J)*QSAT0(J,1)-Q0(J,1))
      ENDDO
      ENDIF

C       For tracers we're taking land isotope ratio, 
C       do can simply scale the corrected total. 
C       (Notice, this needs to be better for production runs).

C      DO ITR=2,NTR
C        DO J=1,NGP
C          EVAP(J,1,ITR) = Rlnd(ITR)*EVAP(J,1,1)
C        ENDDO
C      ENDDO

CCCCCCCCCCCCCCCCCCCCC LAND MODEL SECTION 1 START CCCCCCCCCCCCCCCCCCCCCCC


          DO ITR=2,NTR
            DO J = 1,NGP

C Option 1: Craig Gordon Equation 
              IF (CRAIGG) THEN
                IF (FMASK(J) .gt. 0.0) then
                  hum = min((Q0(J,1)/QSAT0(J,1)),0.9999)
                  hum = max(hum,0.0)

                  alpeq= alpliq(ITR,TSKIN(J))
                    if (Q0(J,2) .lt. qtiny) then 
                      Rvap = 1.		
                    else
                      Rvap = Q0(J,ITR)/Q0(J,2)
                    endif
       
C Add Ice Condition + CG + Vegetation

                  IF (TSKIN(J) .gt. 273.16) then
                   REVAP(J,ITR)=(difr(ITR)**enn)*(((RLD(J,ITR)/alpeq)-
     &                          Rvap*hum)/(1-hum))
                    SEVAP(J,1,ITR) = REVAP(J,ITR)*EVAP(J,1,1) ! soil
                    VEVAP(J,1,ITR) = RB(J,ITR)*EVAP(J,1,1)    ! veg

C Calculate isotopic evaporation as a function of transpiration fraction (fracT)  
             
		    EVAP(J,1,ITR) = (1-fracT)*SEVAP(J,1,ITR) + 
     &              fracT*VEVAP(J,1,ITR)
          
                  ELSE IF (TSKIN(J) .le. 273.16) then
                    EVAP(J,1,ITR) = RLD(J,ITR)*EVAP(J,1,1)
                  ENDIF
                ELSE 
                 EVAP(J,1,ITR) = Rlnd(ITR)*EVAP(J,1,1)
                ENDIF

C Option 2: all water diffuses from ground into evap using land model iso ratio.
              ELSE
                IF (FMASK(J) .gt. 0.0) then
                   EVAP(J,1,ITR) = RLD(J,ITR)*EVAP(J,1,1)
                ENDIF
              ENDIF
              
            ENDDO
          ENDDO

CCCCCCCCCCCCCCCCCCCCC LAND MODEL SECTION 1 END CCCCCCCCCCCCCCCCCCCCCCC

C--   3. Surface energy balance

C     3.1. Emission of lw radiation from the surface
C          and net heat fluxes into land and sea surface

      ESBC  = EMISFC*SBC
      ESBC4 = 4.*ESBC

      DO J=1,NGP

        TSK3     = TSKIN(J)**3
        DSLR(J)  = ESBC4*TSK3
        SLRU(J,1)     = ESBC *TSK3*TSKIN(J)
        SLRU(J,2)     = ESBC*TSEA(J)**4
        SLRU(J,3)  = SLRU(J,2)+FMASK(J)*(SLRU(J,1)-SLRU(J,2))

        HFLUXN(J,1) = SSR(J)+SLRD(J)- 
     &                   (SLRU(J,1)+SHF(J,1)+ALHC*EVAP(J,1,1))

      ENDDO

C     3.2 Re-definition of skin temperature from energy balance

      IF ( LSKINEB ) THEN

C       Compute net heat flux including flux into ground
        DO J=1,NGP
          HFLUXN(J,1) = HFLUXN(J,1)-CLAMBDA*(TSKIN(J)-TLAND(J))
          DTSKIN(J)   = TSKIN(J)+1.
        ENDDO

C       Compute d(Evap) for a 1-degree increment of Tskin

        CALL SHTORH (0,NGP,DTSKIN,PSA,1.,QDUMMY,RDUMMY,QSAT0(1,2))
        
        IF (BUCKET) THEN
        DO J=1,NGP
          IF (EVAP(J,1,1).GT.0) THEN
             QSAT0(J,2) = SOILW(J)*(QSAT0(J,2)-QSAT0(J,1))
          ELSE
             QSAT0(J,2) = 0.
          ENDIF
        ENDDO
        
	ELSE
        DO J=1,NGP
          IF (EVAP(J,1,1).GT.0) THEN
             QSAT0(J,2) = SWAV(J)*(QSAT0(J,2)-QSAT0(J,1))
          ELSE
             QSAT0(J,2) = 0.
          ENDIF
        ENDDO
        ENDIF

C       Redefine skin temperature to balance the heat budget 
        DO J=1,NGP
          DHFDT     = CLAMBDA+DSLR(J)+
     &                CHL*DENVVS(J,1)*(CP+ALHC*QSAT0(J,2))
          DTSKIN(J) = HFLUXN(J,1)/DHFDT
          TSKIN(J)  = TSKIN(J)+DTSKIN(J)
        ENDDO

C       Add linear corrections to heat fluxes
        DO J=1,NGP
          SHF(J,1)    = SHF(J,1) +CHLCP*DENVVS(J,1)*DTSKIN(J)
          EVAP(J,1,1)   = EVAP(J,1,1)+
     &                       CHL*DENVVS(J,1)*QSAT0(J,2)*DTSKIN(J)
          SLRU(J,1)     = SLRU(J,1)  + DSLR(J)*DTSKIN(J)
cfk          SLRU(J,1)   = SLRU(J,1)  +FMASK(J)*DSLR(J)*DTSKIN(J)
          HFLUXN(J,1) = CLAMBDA*(TSKIN(J)-TLAND(J))
        ENDDO

C      DO ITR=2,NTR
C        DO J=1,NGP
C          EVAP(J,1,ITR) = Rlnd(ITR)*EVAP(J,1,1)
C        ENDDO
C      ENDDO

CCCCCCCCCCCCCCCCCCCCC LAND MODEL SECTION 2 CCCCCCCCCCCCCCCCCCCCCCC


          DO ITR=2,NTR
            DO J = 1,NGP

C Option 1: Craig Gordon Equation 
              IF (CRAIGG) THEN
                IF (FMASK(J) .gt. 0.0) then
                  hum = min((Q0(J,1)/QSAT0(J,1)),0.9999)
                  hum = max(hum,0.0)

                  alpeq= alpliq(ITR,TSKIN(J))
                    if (Q0(J,2) .lt. qtiny) then 
                      Rvap = 1.		
                    else
                      Rvap = Q0(J,ITR)/Q0(J,2)
                    endif
       
C Add Ice Condition + CG + Vegetation

                  IF (TSKIN(J) .gt. 273.16) then
                   REVAP(J,ITR)=(difr(ITR)**enn)*(((RLD(J,ITR)/alpeq)-
     &                          Rvap*hum)/(1-hum))
                    SEVAP(J,1,ITR) = REVAP(J,ITR)*EVAP(J,1,1) ! soil
                    VEVAP(J,1,ITR) = RB(J,ITR)*EVAP(J,1,1)    ! veg

C Calculate isotopic evaporation as a function of transpiration fraction (fracT)  
             
		    EVAP(J,1,ITR) = (1-fracT)*SEVAP(J,1,ITR) + 
     &              fracT*VEVAP(J,1,ITR)

                  ELSE IF (TSKIN(J) .le. 273.16) then
                    EVAP(J,1,ITR) = RLD(J,ITR)*EVAP(J,1,1)
                  ENDIF
                ELSE 
                 EVAP(J,1,ITR) = Rlnd(ITR)*EVAP(J,1,1)
                ENDIF

C Option 2: all water diffuses from ground into evap using land model iso ratio.
              ELSE
                IF (FMASK(J) .gt. 0.0) then
                   EVAP(J,1,ITR) = RLD(J,ITR)*EVAP(J,1,1)
                ENDIF
              ENDIF
  
            ENDDO
          ENDDO

CCCCCCCCCCCCCCCCCCCCC LAND MODEL SECTION 2 END CCCCCCCCCCCCCCCCCCCCCCC
      ENDIF
      ENDIF

C     B)   SEA FLUXES 

C     1.6 Stability correction

      RDTH  = FSTAB/DTHETA
      ASTAB = 1.
      IF (LSCASYM) ASTAB = 0.5   ! to get smaller dS/dT in stable conditions

      DO J=1,NGP
        IF (TSEA(J).GT.T2(J,2)) THEN
           DTHS=MIN(DTHETA,TSEA(J)-T2(J,2))
        ELSE
           DTHS=MAX(-DTHETA,ASTAB*(TSEA(J)-T2(J,2)))
        ENDIF
        DENVVS(J,2)=DENVVS(J,0)*(1.+DTHS*RDTH)
      ENDDO

C--   2. Computation of fluxes over sea

C     2.1 Wind stress

      K2 = 0
      IF (LSCDRAG) K2 = 2

      DO J=1,NGP

        CDSDV = CDS*DENVVS(J,K2)

        USTR(J,2) = -CDSDV*UA(J,NLEV)
        VSTR(J,2) = -CDSDV*VA(J,NLEV)

      ENDDO

C     2.2 Sensible heat flux 

      CHSCP = CHS*CP

      DO J=1,NGP
        SHF(J,2) = CHSCP*DENVVS(J,2)*(TSEA(J) -T1(J,2))
      ENDDO

C     2.3 Evaporation

      CALL SHTORH (0,NGP,TSEA ,PSA,1.,QDUMMY,RDUMMY,QSAT0(1,2))

      DO J=1,NGP
        EVAP(J,2,1) = CHS*DENVVS(J,2)*(QSAT0(J,2)-Q0(J,1))
      ENDDO
      DO ITR=2,NTR
        DO J=1,NGP
C          VMAG = SQRT(U0(J)*U0(J)+V0(J)*V0(J)+VG2)
          CHSISO = ALPKOC(ITR,U0(J),V0(J),DENS(J,NLEV))*CHS
          QSATISO = Rocn(ITR)*QSAT0(J,2)/ALPLIQ(ITR,TSEA(j))
          EVAP(J,2,ITR) = CHSISO*DENVVS(J,2)*(QSATISO-Q0(J,ITR))
        ENDDO
      ENDDO


C--   3. Surface energy balance

C     3.1. Emission of lw radiation from the surface
C          and net heat fluxes into land and sea surface

      ESBC  = EMISFC*SBC
      ESBC4 = 4.*ESBC

      DO J=1,NGP
        SLRU(J,2)   = ESBC*TSEA(J)**4
        HFLUXN(J,2) = SSR(J)+SLRD(J)-
     &                  (SLRU(J,2)+SHF(J,2)+ALHC*EVAP(J,2,1))
      ENDDO

      IF ( LFLUXLAND )  THEN

C--   3. Weighted average of surface fluxes and temperatures 
C--      according to land-sea mask

      DO J=1,NGP
        USTR(J,3) = USTR(J,2)+FMASK(J)*(USTR(J,1)-USTR(J,2))
        VSTR(J,3) = VSTR(J,2)+FMASK(J)*(VSTR(J,1)-VSTR(J,2))
         SHF(J,3) =  SHF(J,2)+FMASK(J)*( SHF(J,1)- SHF(J,2))
        SLRU(J,3)  = SLRU(J,2)+FMASK(J)*(SLRU(J,1)-SLRU(J,2))
      ENDDO

      DO ITR=1,NTR
        DO J=1,NGP
          EVAP(J,3,ITR) = EVAP(J,2,ITR)+
     &                 FMASK(J)*(EVAP(J,1,ITR)-EVAP(J,2,ITR))
        ENDDO
      ENDDO

      DO J=1,NGP
        TSFC(J)  = TSEA(J)+FMASK(J)*(TLAND(J)-TSEA(J))
        TSKIN(J) = TSEA(J)+FMASK(J)*(TSKIN(J)-TSEA(J))
      ENDDO

      ENDIF

      RETURN
      END

      SUBROUTINE SFLSET (PHI0)
C--
C--   SUBROUTINE SFLSET (PHI0)
C--
C--   Purpose: compute orographic factor for land surface drag
C--   Input:   PHI0   = surface geopotential            (2-dim)
C--            Initialized common blocks: SFLFIX

C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude
      include "com_physcon.h"

C     Surface flux constants
      include "com_sflcon.h"

      REAL PHI0(NGP)

      RHDRAG = 1./(GG*HDRAG)

      DO J=1,NGP
        FOROG(J)=1.+FHDRAG*(1.-EXP(-MAX(PHI0(J),0.)*RHDRAG))
      ENDDO

C--
      RETURN
      END

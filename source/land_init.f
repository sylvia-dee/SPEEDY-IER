      SUBROUTINE LAND_INIT(ISTART,STANOML,PHI0,PHIS0,GAMMA,GRAV) 
C--
C--   SUBROUTINE LAND_INIT(ISTART,STANOML,PHI0,PHIS0,GAMMA,GRAV)
C--
C--   Purpose : Initialization of land module
C--   Input   : ISTART  = Parameter for deciding if restart run
C--             PHI0    = Orography
C--             PHIS0   = Spectrally truncated orography
C--             GAMMA  = Ref. temperature lapse rate (-dT/dz in deg/km) 
C--             GRAV  = Gravity
C--   Output  : STANOML = Init. value of land sfc temp. anom.
C--             Initialized common blocks: LFORFIX, 
C--                                        LFORMON, lheatflx,
C--                                        lsfconst, lsfcanom,
C--                                        LLSMASK, LFORFIX,
C--                                        LAND_ISTEPS, LAND_ISTFOR,
C--                                        LFORCON,land_anom  
C--	
C--	
							
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL)

      include "com_ts_land.h"
      include "com_forcon_land.h" 
      include "com_land.h"

      real STANOML(IX,IL)

      real*4 veg(ix,il), swl1(ix,il), swl2(ix,il)
      real*4 r4inp(ix,il), dummy4

      real PHI0(ix,il), PHIS0(ix,il)
      real GAMMA, GRAV

      logical LPPRES, LCO2

C
      include "cls_instep.h" 
      include "cls_inphys.h" 

cfk#if 1
      do i=1,ix
       do j=1,il
         phi0land(ix,il) = phi0(ix,il)
         phis0land(ix,il) = phis0(ix,il)
       enddo
      enddo
cfk#endif

      iitest=1

c     set threshold for land-sea mask definition
c     (ie minimum fraction of either land or sea)

      thrsh = 0.1


C--   1. Read time-invariant fields (Land-sea masks and albedos)

c     set threshold for land-sea mask definition
c     (ie minimum fraction of either land or sea)

      read (20) dummy4

      if (iitest.ge.1) print*,' read fractional land-sea mask'  

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          fmask1(i,j) = r4inp(i,j)
        enddo
      enddo

      if (iitest.ge.1) print*,' read surface albedo'  

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          alb0(i,j) = 0.01*r4inp(i,j)
        enddo
      enddo

      REWIND 20 

C--   1.1 Compute additional land-sea masks
c
c     fmask0  = 1 - fmask1 = sea fraction
c     fmaskl1 = 1 wherever there is any land, 0 elsewhere
c     fmasko1 = 1 wherever there is any sea,  0 elsewhere
c
c     Note that because some points have both land and ocean
c     it is *not* true that: fmasko1+fmaskl1=1 everywhere

      do i=1,ix
       do j=1,il

         fmask0(i,j)=1.0-fmask1(i,j)

         if (fmask1(i,j).ge.thrsh) then
           fmaskl1(i,j)=1.0
         else
           fmaskl1(i,j)=0.0
           fmask1 (i,j)=0.0
           fmask0 (i,j)=1.0
         endif

         if (fmask0(i,j).ge.thrsh) then
           fmasko1(i,j)=1.0
         else
           fmasko1(i,j)=0.0
           fmask0 (i,j)=0.0
           fmask1 (i,j)=1.0
         endif

       enddo
      enddo


C--   2. Read monthly-mean climatologies of surface fields

c     2.1 Read land-surface temp

      if (iitest.ge.1) print*,' reading land-surface temp.'
  
      do it = 1,12
        read (23) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            stl12(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.eq.1) print*,' checking land-surface temp.'

      CALL FORCHK (fmaskl1,stl12,ix*il,12,0.,400.,273.)

C--  Orographic adaptation of land sfc temperature 

C
C       GAM = 0.001*GAMMA/GRAV
C
       do it = 1,12
C
Cc        DO J=1,NLAT
Cc         DO I=1,NLON
Cc           STL12(I,J,IT)=STL12(I,J,IT)+GAM*(PHI0(I,J)-PHIS0(I,J))
Cc         ENDDO
Cc        ENDDO   
C
        call ftland (stl12(1,1,it),phi0,phis0,fmaskl1)
C
C        if (iitest.gt.1) then
C          do j = 1,il
C            do i = 1,ix
C              r4inp(i,j) = stl12(i,j,it)
C            enddo
C          enddo
C          write (18) ((r4inp(i,j),i=1,ix),j=il,1,-1)
C        endif
C
      enddo   
C

c     2.2 Read snow depth

      if (iitest.ge.1) print*,' reading snow depth'  

      do it = 1,12
        read (24) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            snow12(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking snow depth'

      CALL FORCHK (fmaskl1,snow12,ix*il,12,0.,20000.,0.)

c    2.3 Read soil moisture and compute soil water availability 
c        using vegetation fraction

      if (iitest.ge.1) print*,' reading soil moisture'  

c     read vegetation fraction (in %)
      read (25) ((veg(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          veg(i,j)=max(0.,0.01*veg(i,j))
        enddo
      enddo

      sdep1 = 70.
      idep2 = 3
      sdep2 = idep2*sdep1

      swwil2= sdep2*swwil
      rsw   = 1./(sdep1*swcap+sdep2*(swcap-swwil))

      do it = 1,12
        read (26) ((swl1(i,j),i=1,ix),j=il,1,-1)
        read (26) ((swl2(i,j),i=1,ix),j=il,1,-1)
        read (26) dummy4
        do j = 1,il
          do i = 1,ix
            swroot = idep2*swl2(i,j)
            soilw12(i,j,it) = rsw*(swl1(i,j)+
     &                        veg(i,j)*max(0.,swroot-swwil2))									
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking soil moisture'

      CALL FORCHK (fmaskl1,soilw12,ix*il,12,0.,10.,0.)


C--   3. Read heat fluxes for surface anomaly model

cfk      if (ialst.gt.0) then
cfk
cfk        if (iitest.ge.1) print*,' reading sfc heat fluxes' 
cfk
cfk        irecl=4*ix*il
cfk        irec =0
cfk
cfk        open ( unit=31, file='fort.31', status='old', 
cfk     &         form='unformatted', access='direct', recl=irecl )
cfk
       do it = 1,12

cfk          irec=irec+1
cfk          read (31,rec=irec) r4inp
          do j = 1,il
            do i = 1,ix
cfk              hflxl12(i,j,it) = r4inp(i,j)
              hflxl12(i,j,it) = 0.
            enddo
          enddo
cfk          irec=irec+1
cfk          read (31,rec=irec) r4inp

        enddo
cfk
cfk        if (iitest.ge.1) print*,' checking sfc heat fluxes'
cfk
cfk        CALL FORCHK (fmaskl1,hflxl12,ix*il,12,-1000.,1000.,0.)
cfk
cfk      endif

C--   4. Initialize land anomaly model
C--  
C--   4.1 Set heat capacities and dissipation times for 
C--       land (root-layer)

C     organic soil layer (depth*heat_cap/m)
c      hcapl =  0.6*2.50e+6 
      hcapl =  1.*2.50e+6
C     land-ice layer (depth*heat_cap/m)
      hcapli =  1.8*1.93e+6


C     Dissipation time (days) for land-surface temperature anomalies
c      tdlsta  = 20.
      tdlsta  = 40.

      rhcapl = 86400./hcapl
      rhcapli = 86400./hcapli

      stdisl = tdlsta/(tdlsta+1.)

C--   4.2 Set variable mixed-layer depth over selected regions

      nlat2=nlat/2

c     Continents/land-ice
      do j=1,nlat
        do i=1,nlon
           if (alb0(i,j).lt.0.4) then
             rhcap2l(i,j) = rhcapl
           else
             rhcap2l(i,j) = 0.5*rhcapli
c             rhcap2l(i,j) = rhcapli
           endif
        enddo
      enddo

C--   4.3 Initialize anomaly model variables

      if (ialst.eq.0.or.istart.eq.0) then
        do j=1,nlat
          do i=1,nlon
            stanoml(i,j)=0.
          enddo
        enddo
      endif

        do j=1,nlat
          do i=1,nlon
            stanoml1(i,j)=stanoml(i,j)
          enddo
        enddo

C      do j=1,nlat
C        do i=1,nlon
C          hfintl(i,j)=0.
C        enddo
C      enddo


cfk      close (unit=31)
C--
      RETURN
      END



      SUBROUTINE FORCHK (FMASK,FIELD,NGP,NF,FMIN,FMAX,FSET)
										
C--   Aux. routine FORCHK: Check consistency of sfc fields with land-sea mask 
C--   and set undefined values to a constant (to avoid over/underflow)

      real fmask(ngp), field(ngp,nf)

      do jf = 1,nf

        nfault=0

        do jgp = 1,ngp
          if (fmask(jgp).gt.0.0) then
            if (field(jgp,jf).lt.fmin.or.field(jgp,jf).gt.fmax)
     *          nfault = nfault+1
          else
            field(jgp,jf) = fset
          endif
        enddo

        print *, ' field: ', jf, '   no. of faulty points:', nfault

      enddo

      print *, ' undefined values set to', fset

      RETURN
      END


      SUBROUTINE FTLAND (STL,PHI0,PHIS0,FMASKL)

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL)

      include "com_dyncon0.h" 
      include "com_dyncon1.h" 

      REAL STL(NLON,NLAT), PHI0(NLON,NLAT), PHIS0(NLON,NLAT),
     &     FMASKL(NLON,NLAT)

      REAL STL2(NLON,NLAT)

      NL8 = NLAT/8
      GAM = 0.001*GAMMA/GRAV

      NLAT1 = 1
      NLAT2 = NL8

      DO JBAND=1,8

        SUMT=0.
        SUMW=0.

        DO J=NLAT1,NLAT2
         DO I=1,NLON
           STL(I,J)=STL(I,J)+GAM*PHI0(I,J)
           SUMT=SUMT+GCOS(J)*FMASKL(I,J)*STL(I,J)
           SUMW=SUMW+GCOS(J)*FMASKL(I,J)
         ENDDO
        ENDDO

        SUMT=SUMT/SUMW

        DO J=NLAT1,NLAT2
         DO I=1,NLON
           IF (FMASKL(I,J).EQ.0.) STL(I,J)=SUMT
         ENDDO
        ENDDO
  
        NLAT1=NLAT1+NL8
        NLAT2=NLAT2+NL8

      ENDDO

      ITR=7
      IDTR=(NTRUN-6)/3

      DO JFIL=1,4

        CALL TRUNCG (ITR,STL,STL2)

        DO J=1,NLAT
         DO I=1,NLON
           IF (FMASKL(I,J).EQ.0.) STL(I,J)=STL2(I,J)
         ENDDO
        ENDDO

        ITR=MIN(ITR+IDTR,NTRUN)

      ENDDO

      CALL TRUNCG (ITR,STL,STL2)

      DO J=1,NLAT
       DO I=1,NLON
         STL(I,J)=STL2(I,J)-GAM*PHIS0(I,J)
       ENDDO
      ENDDO       

      RETURN
      END

      SUBROUTINE TRUNCG (ITR,FG1,FG2)

C--   SUBROUTINE TRUNCG (ITR,FG1,FG2)
C--   Purpose : compute a spectrally-filtered grid-point field
C--   Input   : ITR : spectral truncation (triangular)
C--           : FG1 : original grid-point field
C--   Output  : FG2 : filtered grid-point field

      include "atparam.h"

      REAL FG1 (IX,IL), FG2(IX,IL)
      COMPLEX FSP(MX,NX), ZERO 

      ZERO = (0.,0.)

      CALL SPEC (FG1,FSP)

      DO N=1,NX
        DO M=1,MX
          ITWN=ISC*(M-1)+N-1
          IF (ITWN.GT.ITR) FSP(M,N)=ZERO
        ENDDO
      ENDDO

      CALL GRID (FSP,FG2,1)

      RETURN
      END

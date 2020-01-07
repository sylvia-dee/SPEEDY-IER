      SUBROUTINE SEA_INIT(ISTART,ISSTAN0,SSTAN,
     &                    STANOMI,STANOMS,RADANG) 
C--
C--   SUBROUTINE SEA_INIT(ISTART,ISSTAN0,SSTAN,
C--  &                    STANOMI,STANOMS,RADANG)
C--
C--   Purpose : Initialization of sea module 
C--   Input   : ISTART  = Parameter for deciding if restart run
C--           : ISSTAN  = Parameter for sst anomalies 
C--           : RADANG  = Latitude of gridpoint in y-direc.
C--   Output  : SSTAN   = Init. value of comb. sea-ice sfc temp. anom.
C--             STANOMI = Init. value of ice sfc temp. anom.
C--             STANOMS = Init. value of sea sfc temp. anom.  
C--             Initialized common blocks: SFORMON, sheatflx, 
C--                                        ssfconst, ssfcanom, 
C--                                        SLSMASK, SFORFIX,
C--                                        SEA_ISTEPS, SEA_ISTFOR,
C--                                        sea_anom 
C--	
C--	
							
      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL)

      include "com_ts_sea.h"
      include "com_forcon_sea.h" 
      include "com_sea.h"

      real SSTAN(IX,IL), STANOMI(IX,IL), STANOMS(IX,IL)
      real RADANG(IL)

      real*4 r4inp(ix,il), dummy4

      logical LPPRES, LCO2

C
      include "cls_instep.h"
      include "cls_inphys.h" 

      ISSTAN = ISSTAN0   

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

c     2.1 Read SST 

      if (iitest.ge.1) print*,' reading sst' 

      do it = 1,12
        read (21) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            sst12(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking sst'

      CALL FORCHK (fmasko1,sst12,ix*il,12,0.,400.,273.)

c     2.2 Read sea ice

      if (iitest.ge.1) print*,' reading sea ice'  

      do it = 1,12
        read (22) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            oice12(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking sea ice'

      CALL FORCHK (fmasko1,oice12,ix*il,12,0.,1.,0.)

C--   3. Read SST anomalies for initial and preceeding month

      if (isstan.gt.0) then

        if (iitest.ge.1) print*,' reading sst anomalies' 

        do jrec=1,isst0-2
          read (30) dummy4
        enddo

        read (30) ((r4inp(i,j),i=1,ix),j=il,1,-1)

        do j = 1,il
          do i = 1,ix
            sstan2(i,j,1) = r4inp(i,j)
          enddo
        enddo

        if (isst0.gt.1) read (30) ((r4inp(i,j),i=1,ix),j=il,1,-1)

        do j = 1,il
          do i = 1,ix
            sstan2(i,j,2) = r4inp(i,j)
          enddo
        enddo

        if (iitest.ge.1) print*,' checking sst anomalies'

        CALL FORCHK (fmasko1,sstan2,ix*il,2,-50.,50.,273.)

      endif


C--   4. Read heat fluxes for surface anomaly model

cfk      if (iasst.gt.1.or.iaice.gt.0) then
        if (iasst.gt.1) then

        if (iitest.ge.1) print*,' reading sfc heat fluxes' 

        irecl=4*ix*il
cfk        irecl=ix*il
        irec =0

        open ( unit=31, file='fort.31', status='old', 
     &         form='unformatted', access='direct', recl=irecl )

        do it = 1,12

          irec=irec+1
          read (31,rec=irec) r4inp
          irec=irec+1
          read (31,rec=irec) r4inp
          do j = 1,il
            do i = 1,ix
              hflxs12(i,j,it) = r4inp(i,j)
            enddo
          enddo  

        enddo

        if (iitest.ge.1) print*,' checking sfc heat fluxes'

        CALL FORCHK (fmasko1,hflxs12,ix*il,12,-1000.,1000.,0.)

      endif


C--   5. Initialize sea anomaly model
C--  
C--   5.1 Set heat capacities and dissipation times for 
C--      sea (mixed layer) and sea-ice 

C     oceanic mixed layer (depth*heat_cap/m)
      hcaps = 50.0*4.18e+6
C     sea-ice layer (depth*heat_cap/m)
      hcapi =  1.8*1.93e+6
C     Dissipation time (days) for sea-surface temperature anomalies
      tdssta  = 90.
C     Dissipation time (days) for sea-ice temperature anomalies
c      tdiceta = 10.
      tdiceta = 20.

      rhcaps = 86400./hcaps
      rhcapi = 86400./hcapi

      stdiss = tdssta/(tdssta+1.)
      stdisi = tdiceta/(tdiceta+1.)

      flxice = hcapi/(86400.*tdiceta)

C--   5.2 Set variable mixed-layer depth over selected regions

      nlat2=nlat/2

cfk
cfk   mixed-layer everywhere
cfk
      do j=1,nlat
        do i=1,nlon
          rhcap2s(i,j) = rhcaps
        enddo
      enddo

c     Oceans
cfk      do j=1,nlat
cfk        do i=1,nlon
cfk          rhcap2s(i,j) = -1.
cfk        enddo
cfk      enddo

c     Atlantic Ocean
cfk      dlon=360./nlon
cfk      reast=42.
cfk      rwest=300.-dlon
cfk      do j=nlat2+1,nlat
cfk        do i=1,nlon
cfk          rlon=(i-1)*dlon
cfk          if (rlon.lt.reast.or.rlon.gt.rwest) rhcap2s(i,j) = rhcaps
cfk        enddo
cfk        rwest=max(260.,rwest-2*dlon)
cfk      enddo

cfk      do i=1,nlon
cfk        rhcap2s(i,nlat2+1)=rhcap2s(i,nlat2+1)*0.25
cfk        rhcap2s(i,nlat2+2)=rhcap2s(i,nlat2+2)*0.5
cfk      enddo

C--   5.3 Set weight mask for observed SST anomalies

      do j=1,nlat
        do i=1,nlon
          wobsst(i,j) = 0.
        enddo
      enddo

c     Tropics ( weight = sqrt(cos(3*lat)) for abs(lat) < 30 )
c     Indian + Pacific only

      nwest=1+nint(nlon/16.)
      neast=1+nint(nlon*290./360.)
      rad30=asin(0.5)

      do j=1,nlat
        if (abs(radang(j)).lt.rad30) then
          wob1 = sqrt(cos(3.*radang(j)))
          do i=nwest,neast
            if (rhcap2s(i,j).lt.0.) wobsst(i,j) = wob1
          enddo
        endif
      enddo


C--   5.4 Initialize anomaly model variables

      if (iasst.eq.0.or.istart.eq.0) then
        do j=1,nlat
          do i=1,nlon
            stanoms(i,j)=0.
          enddo
        enddo
      endif

        do j=1,nlat
          do i=1,nlon
            stanoms1(i,j)=stanoms(i,j)
          enddo
        enddo

      if (iaice.eq.0.or.istart.eq.0) then
        do j=1,nlat
          do i=1,nlon
            stanomi(i,j)=0.
          enddo
        enddo
      endif

        do j=1,nlat
          do i=1,nlon
            stanomi1(i,j)=stanomi(i,j)
          enddo
        enddo

      do j=1,nlat
        do i=1,nlon
          sstan(i,j)=0.
        enddo
      enddo

C      do j=1,nlat
C        do i=1,nlon
C          hfints(i,j)=0.
C        enddo
C      enddo



      close (unit=31)

C--
      RETURN
      END




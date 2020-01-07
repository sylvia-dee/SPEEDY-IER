      SUBROUTINE SEA2ATM(IMONTH,IDAY,IRST,SST1,SST01,SST01L,OICE1,ALBS,
     &                   SSTAN,STANOMI,STANOMS,HFLXS1)
C--
C--   SUBROUTINE SEA2ATM(IMONTH,IDAY,SST1,SST01,OICE1,ALBS,
C--  &                   SSTAN,STANOMI,STANOMS,HFLXS1) 
C--   
C--   Purpose :	Compute forcing fields from the sea for the 
C--             current date
C--   Input   : IMONTH  = Actual month 
C--             IDAY    = Actual day
C--   Output  : SST1    = Daily value of sst
C--             SST01   = Daily value of climat. sst 
C--             SST01L  = Daily value of climat. sst prev. day 
C--             OICE1   = Daily value of ice 
C--             ALBS    = Daily value of albedo 
C--             SSTAN   = Daily value of sst anomaly
C--             STANOMI = Daily value of ice sfc temp. anom.
C--             STANOMS = Daily value of sea sfc temp. anom.  
C--             HFLXS1  = Daily value of climat. sea heatflux 
C--   Modified common blocks: SFORDAY, sstanom, sheatflx
C--                           
C--
      include "atparam.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )
cfk---modification for coupling start
      PARAMETER (IXP2=IX+2,ILP2=IL+2)
cfk---modification for coupling end

      include "com_ts_sea.h"
      include "com_forcon_sea.h"
      include "com_sea.h"
C
      integer IRST
  
      real SST1(IX,IL),SST01(IX,IL),SST01L(IX,IL),OICE1(IX,IL),
     &     ALBS(IX,IL),SSTAN(IX,IL),STANOMI(IX,IL),STANOMS(IX,IL),
     &     HFLXS1(IX,IL)
cfk---modification for coupling start
      real SSTCZ(IX,IL)
cfk---modification for coupling end
cfk---modification for coupling start
      integer IALAND(IXP2,ILP2)
cfk---modification for coupling end
 

      IF(IRST.EQ.1) THEN

C--   1. Define the interpolation parameters
 
C     1.1 Define the fraction of month for interpolation 
C         (middle of the month for fixed-month integrations) 

      NDAYS = 30
      IDAYH = 1+NDAYS/2

      IF (IDAY.EQ.0) THEN
        PRINT *, ' Start of SEA2ATM routine (init.)'
        IF (ISEASC.EQ.1) THEN
          FDAY = 0.5
        ELSE
          FDAY = 0.5*NDAYS
        ENDIF
      ELSE
        IF (ISEASC.EQ.1) THEN
          FDAY = IDAY-0.5
        ELSE
          GO TO 900
        ENDIF
      ENDIF

      FMON = FDAY/NDAYS
      TYEAR = (IMONTH-1+FMON)/12.

C--   Surface tempreatures  from previous day

      DO J=1,NLAT
        DO I=1,NLON
            SST01l(I,J)=SST01(I,J)
        ENDDO
      ENDDO


C--   2. Interpolate monthly-mean surface fields

C     Non-linear, mean-conserving interpolation of sea temperature
      CALL FORIN5 (NGP,IMONTH,FMON,SST12,SST01)
C     Linear interpolation of ice fields
      CALL FORINT (NGP,IMONTH,FMON,OICE12,OICE1)

C--   3. Interpolate SST anomaly and heat fluxes 

      IF (ISSTAN.GT.0) THEN

        IF (IDAY.EQ.IDAYH) CALL NEWSST

        IF (IDAY.LT.IDAYH) THEN
          CALL FORINT (NGP,2,FMON,SSTAN2,SSTAN1)
        ELSE
          CALL FORINT (NGP,1,FMON,SSTAN2,SSTAN1)
        ENDIF

      ENDIF

      IF (IASST.GT.0.OR.IAICE.GT.0) THEN

        CALL FORIN5 (NGP,IMONTH,FMON,HFLXS12,HFLXS1)

      ENDIF

C--   4. Surface albedo:
C        sea albedo depends linearly on sea-ice fraction. 

      DALB=ALBICE-ALBSEA

      DO J=1,NLAT
        DO I=1,NLON
          ALBS(I,J)=ALBSEA+DALB*OICE1(I,J)
        ENDDO
      ENDDO



C--   5. Add sfc temp. anomaly to climatological sfc temp.

 900  CONTINUE


C--   5.1 Define SST anomaly from prescribed and/or mixed-layer values

      if (iasst.eq.1.or.iasst.eq.3) then

        sstice=273.2
        do j=1,nlat
          do i=1,nlon
            if (sst01(i,j).le.sstice) sstan1(i,j)=0.
          enddo
        enddo

      endif

      if (iasst.eq.1) then

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=sstan1(i,j)
          enddo
        enddo

      else if (iasst.eq.2) then

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=stanoms1(i,j)
          enddo
        enddo

      else if (iasst.eq.3) then

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=stanoms1(i,j)+
     &                 wobsst(i,j)*(sstan1(i,j)-stanoms1(i,j))
          enddo
        enddo

      else if (iasst.eq.4) then

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=stanoms1(i,j)
          enddo
        enddo

      else if (iasst.eq.7) then

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=sstan1(i,j)
          enddo
        enddo

      endif

C--   5.2 Correct SST using sea-ice temperature anomaly if requested

      if (iaice.gt.0) then

cfk---modification for coupling (mentioned by mpking) start
      if (iasst .eq.0) then
        do j=1,nlat
          do i=1,nlon
           sstan(i,j)=0.
          enddo
        enddo
      endif
cfk---modification for coupling (mentioned by mpking) end

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=sstan(i,j)+
     &                 oice1(i,j)*(stanomi1(i,j)-sstan(i,j))
          enddo
        enddo

      endif

C--   5.3 Superimpose anomalies to sea temperature climatology

      do j=1,nlat
        do i=1,nlon
          sst1(i,j)=sst01(i,j)+sstan(i,j)
        enddo
      enddo

cfk---modification for coupling start
      IF(IASST .EQ. 4 .OR. IASST .EQ. 5 .OR. 
     &   IASST .EQ. 6 .OR. IASST .EQ. 7) THEN
cfk
cfk    do not call ocean for initialization
cfk
       IF (IDAY.GT.0) THEN 
        CALL READERA (SSTCZ,IALAND)
cfk---use ocean sst in common atmosphere and ocean sea points 
c      do j=1,nlat
        do j=17,32
c      do j=16,33
         do i=1,nlon
           if ((fmasko1(i,j).eq.1.).and.(ialand(i+1,j+1).eq.1)) then
               sst1(i,j)=sstcz(i,j)
           endif
         enddo
        enddo
       ENDIF
      ENDIF     
cfk---modification for coupling end

C--   5.4 Define anomalies for atmosphere
  
      do j=1,nlat
        do i=1,nlon
          stanomi(i,j)=stanomi1(i,j)
          stanoms(i,j)=stanoms1(i,j)
        enddo
      enddo            

      ELSE IF(IRST.EQ.2) THEN

C--   5.4 Define anomalies for atmosphere
  
      do j=1,nlat
        do i=1,nlon
          stanomi(i,j)=stanomi1(i,j)
          stanoms(i,j)=stanoms1(i,j)
        enddo
      enddo            

      ENDIF

C--
      RETURN
      END



      SUBROUTINE FORINT (NGP,IMON,FMON,FOR12,FOR1)  

C--   Aux. routine FORINT : linear interpolation of monthly-mean forcing

      REAL FOR12(NGP,*), FOR1(NGP)

      IF (FMON.LE.0.5) THEN
        IMON2 = IMON-1
        IF (IMON.EQ.1) IMON2 = 12
        WMON = 0.5-FMON
      ELSE
        IMON2 = IMON+1
        IF (IMON.EQ.12) IMON2 = 1
        WMON = FMON-0.5
      ENDIF

      DO J=1,NGP
        FOR1(J) = FOR12(J,IMON)+WMON*(FOR12(J,IMON2)-FOR12(J,IMON))
      ENDDO
C--
      RETURN
      END

      subroutine FORIN5 (ngp,imon,fmon,for12,for1)

C--   Aux. routine FORIN5 : non-linear, mean-conserving interpolation 
C--                         of monthly-mean forcing fields

      real for12(ngp,12), for1(ngp)

      im2 = imon-2
      im1 = imon-1
      ip1 = imon+1
      ip2 = imon+2

      if (im2.lt.1)  im2 = im2+12
      if (im1.lt.1)  im1 = im1+12
      if (ip1.gt.12) ip1 = ip1-12
      if (ip2.gt.12) ip2 = ip2-12
 
      c0 = 1./12.
      t0 = c0*fmon
      t1 = c0*(1.-fmon)
      t2 = 0.25*fmon*(1-fmon)

      wm2 =        -t1   +t2
      wm1 =  -c0 +8*t1 -6*t2
      w0  = 7*c0      +10*t2     
      wp1 =  -c0 +8*t0 -6*t2
      wp2 =        -t0   +t2 

      do j=1,ngp
        for1(j) = wm2*for12(j,im2)+wm1*for12(j,im1)
     &           + w0*for12(j,imon)
     &           +wp1*for12(j,ip1)+wp2*for12(j,ip2)
      enddo

      return
      end

      SUBROUTINE NEWSST  

C--   Aux. routine NEWSST : update SST anomaly field 

      include "atparam.h"
      include "com_sea.h"


      real*4 r4inp(ix,il)

      do j = 1,il
        do i = 1,ix
          sstan2(i,j,1) = sstan2(i,j,2)
        enddo
      enddo

      read(30,end=100) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          sstan2(i,j,2) = r4inp(i,j)
        enddo
      enddo

      CALL FORCHK (fmasko1,sstan2(1,1,2),ix*il,1,-50.,50.,273.)

      RETURN

 100  continue

      print *, ' WARNING: end-of-file reached on SST anomaly file'
      print *, ' SST anomaly will be kept constant'

C--
      RETURN
      END




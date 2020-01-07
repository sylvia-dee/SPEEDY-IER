      SUBROUTINE LAND2ATM(IMONTH,IDAY,IRST,STL1,STL01,STL01L,SNOW1,
     &                    SOILW1,ALBL,STANOML,HFLXL1)
C--
C--   SUBROUTINE LAND2ATM(IMONTH,IDAY,STL1,STL01,SNOW1,SOILW1,
C--  &                    ALBL,STANOML,HFLXL1) 
C--   
C--   Purpose :	Compute forcing fields from the land model 
C--             for the current date
C--   Input   : IMONTH = Actual month 
C--             IDAY   = Actual day
C--   Output  : STL1   = Daily value of land-surf. temp. 
C--             STL01  = Daily value of climat. land-surf. temp.
C--             STL01L = Daily value of climat. land-surf. temp. prev. day
C--             SNOW1  = Daily value of snow cover 
C--             SOILW1 = Daily value of soil wetness
C--             ALBL   = Daily value of albedo over land
C--             STANOML= Daily value of lst anomaly 
C--             HFLXL1 = Daily value of climat. land heatflux 
C--   Modified common blocks: LFORDAY, sstanom, lheatflx
C--
      include "atparam.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NGP=NLON*NLAT )

      include "com_ts_land.h"
      include "com_forcon_land.h"
      include "com_land.h"
C
      integer IRST

      real STL1(IX,IL),STL01(IX,IL),STL01L(IX,IL),SNOW1(IX,IL),
     &     SOILW1(IX,IL),ALBL(IX,IL),STANOML(IX,IL),HFLXL1(IX,IL)

cfk#if 1
      real GAM, GRAV

      include 'cls_indyns.h'
      GRAV=9.81
cfk#endif


      IF(IRST.EQ.1) THEN


C--   Surface tempreatures  from previous day

      DO J=1,NLAT
        DO I=1,NLON
            STL01l(I,J)=STL01(I,J)
        ENDDO
      ENDDO

C--   1. Define the interpolation parameters

C     1.1 Define the fraction of month for interpolation 
C         (middle of the month for fixed-month integrations) 

      NDAYS = 30
      IDAYH = 1+NDAYS/2

      IF (IDAY.EQ.0) THEN
        PRINT *, ' Start of LAND2ATM routine (init.)'
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

C--   2. Interpolate monthly-mean surface fields

C     Non-linear, mean-conserving interpolation of land temperature
      CALL FORIN5 (NGP,IMONTH,FMON,STL12,STL01)

C     Linear interpolation of water/snow fields
      CALL FORINT (NGP,IMONTH,FMON,SNOW12,SNOW1)
      CALL FORINT (NGP,IMONTH,FMON,SOILW12,SOILW1)

C--   3. Interpolate heat fluxes 

      IF (IALST.GT.0) THEN

        CALL FORIN5 (NGP,IMONTH,FMON,HFLXL12,HFLXL1)

      ENDIF

C--   4. Set maximum soil water availability to 1

      DO J=1,NLAT
        DO I=1,NLON
           SOILW1(I,J)=MIN(1.,SOILW1(I,J))
        ENDDO
      ENDDO

C--   5. Surface albedo:
C         land albedo depends linearly on snow depth (up to the SDALB
C         threshold). 

      RSD=1./SDALB

      DO J=1,NLAT
        DO I=1,NLON
          SNOWC1(I,J)=MIN(1.,RSD*SNOW1(I,J))
          ALBL(I,J)=ALB0(I,J)+MAX(ALBSN-ALB0(I,J),0.0)*SNOWC1(I,J)
        ENDDO
      ENDDO


C--   6. Add sfc temp. anomaly to climatological sfc temp.

 900  CONTINUE

C--   6.1 Superimpose anomalies to land  temperature climatology

      do j=1,nlat
        do i=1,nlon
          stl1(i,j)=stl01(i,j)+stanoml1(i,j)
        enddo
      enddo

cfk#if 1
c      GAM = 0.001*GAMMA/GRAV
c      DO J=1,NLAT
c         DO I=1,NLON
c            STL1(I,J)=STL1(I,J)+GAM*(PHI0land(I,J)-PHIS0land(I,J))
c         ENDDO
c      ENDDO   
cfk#endif

C--   6.2 Define anomalies for atmosphere
  
      do j=1,nlat
        do i=1,nlon
          stanoml(i,j)=stanoml1(i,j)
        enddo
      enddo

      ELSE IF(IRST.EQ.2) THEN 

C--   6.2 Define anomalies for atmosphere
  
      do j=1,nlat
        do i=1,nlon
          stanoml(i,j)=stanoml1(i,j)
        enddo
      enddo

      ENDIF           

C--
      RETURN
      END

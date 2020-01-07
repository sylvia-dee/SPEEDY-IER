      SUBROUTINE FORDATE
C--
C--   SUBROUTINE FORDATE 
C--   
C--   Purpose :	Compute forcing fields for the current date
C--             and correction terms for horiz. diffusion
C--   Modified common blocks: DATE1, FORDAY, HDIFC4, ANOM
C--
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_tsteps.h"
      include "com_date.h"

      include "com_forcon.h"
      include "com_dyncon0.h"
      include "com_physcon.h"
      include "com_radcon.h"

      include "com_hdifcon.h"
      include "com_forcing.h"
      include "com_for_sea.h"
      include "com_for_land.h"
      include "com_lflags.h"


      REAL GAMLAT(NLAT),
     &     CORH(NLON,NLAT), TSFC(NLON,NLAT), TREF(NLON,NLAT),
     &     PSFC(NLON,NLAT), QSFC(NLON,NLAT), QREF(NLON,NLAT)

      iitest = 0

C--   1. Define the interpolation parameters
 
C     1.1 Define the fraction of month for interpolation 
C         (middle of the month for fixed-month integrations) 

      NDAYS = 30
      IDAYH = 1+NDAYS/2

      IF (IDAY.EQ.0) THEN
        PRINT *, ' Start of FORDATE routine (init.)'
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

C--   2. Surface albedo:
C         defined as a weighed average of land and ocean albedos


      DO J=1,NLAT
        DO I=1,NLON
          ALB1(I,J)=FMASK1(I,J)*ALBL(I,J)+FMASK0(I,J)*ALBS(I,J)
        ENDDO
      ENDDO

C--   3. Call flow-independent parts of physical parametrizations

      IF (IDAY.EQ.0) THEN

        CALL RADSET

        CALL SFLSET (PHIS0)

C       Reference CO2 absorptivity
        ABLCO2_ref = ABLCO2


      ENDIF

      CALL SOL_OZ (SOLC,TYEAR)

C     Linear trend of CO2 absorptivity (Del_CO2 : rate of change per year)

      IYEAR_ref = 1950
      xkappa    = 0.005
C      xkappa    = 0.0033

      IF (LCO2) THEN
          ABLCO2 = ABLCO2_ref*EXP(xkappa*(IYEAR+TYEAR-IYEAR_ref))
      ENDIF

C--   4. Temperature correction term for horizontal diffusion

      CALL SETGAM (TYEAR,GAMLAT)

      DO J=1,NLAT
        DO I=1,NLON
          CORH(I,J)=GAMLAT(J)*PHIS0(I,J)
        ENDDO
      ENDDO

      if (iitest.gt.1.and.iday.eq.0) then
         call outest (19,PHIS0)
         call outest (19,CORH)
      endif

      CALL SPEC (CORH,TCORH)

C--   5. Humidity correction term for horizontal diffusion

      DO J=1,NLAT
        PEXP=1./(RD*GAMLAT(J))
        DO I=1,NLON
CAS          TSFC(I,J)=FMASK1(I,J)*STL01(I,J)+FMASK0(I,J)*SST01(I,J)
          TSFC(I,J)=FMASK1(I,J)*STL01(I,J)+FMASK0(I,J)*SST01(I,J)
cfk          TSFC(I,J)=FMASK1(I,J)*STL1(I,J)+FMASK0(I,J)*SST1(I,J)
          TREF(I,J)=TSFC(I,J)+CORH(I,J)
          PSFC(I,J)=(TSFC(I,J)/TREF(I,J))**PEXP
        ENDDO
      ENDDO

      CALL SHTORH (0,NGP,TREF,  1.,-1.,DUMMY,DUMMY,QREF)
      CALL SHTORH (0,NGP,TSFC,PSFC, 1.,DUMMY,DUMMY,QSFC)

      DO J=1,NLAT
        DO I=1,NLON
          CORH(I,J)=REFRH1*(QREF(I,J)-QSFC(I,J))
        ENDDO
      ENDDO

      if (iitest.gt.1.and.iday.eq.0) call outest (19,CORH)

      CALL SPEC (CORH,QCORH)

 900  CONTINUE

C--   6. Flux correction

c
cjk
cfk      IF (IMONTH .EQ. 1 .AND. IDAY .LT. 9) THEN
cfk       print*,' '
cfk       print*,'  trop SST1, before FC = ',sst1(49,24),iday
cfk       print*,'  nord SST1, before FC = ',sst1(44,34),iday
cfk       print*,' '
cfk      ENDIF
c

      IF ( IFLUXCORR .GT. 0 )  CALL FLUXCORR (FMON)

c
cjk
cfk      IF (IMONTH .EQ. 1 .AND. IDAY .LT. 9) THEN
cfk       print*,' '
cfk       print*,'  trop SST1,  after FC = ',sst1(49,24),iday
cfk       print*,'  nord SST1,  after FC = ',sst1(44,34),iday
cfk       print*,' '
cfk       print*,' '
cfk      ENDIF
c


C--   7. Nino34 forcing

      IF ( IOBSNINO .GT. 0 .AND. IFLUXCORR .GT. 0 ) 
     &                   CALL OBS_SSTA (IDAYH,FMON)

C      IF (TYEAR .LT. 5./12. .OR. TYEAR .GT. 9./12. ) THEN
         IF ( IASST .EQ. 6 ) CALL OBS_SSTA (IDAYH,FMON)
C      ENDIF

C--
      RETURN
      END


      SUBROUTINE SETGAM (TYEAR,GAMLAT)  

C--   Aux. routine GAMLAT : compute reference lapse rate 
C--                         as a function of latitude and date

      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_dyncon0.h"
      include "com_physcon.h"

      REAL GAMLAT(NLAT)

      GAMLAT(1) = GAMMA/(1000.*GG)
      DO J=2,NLAT
        GAMLAT(J) = GAMLAT(1)
      ENDDO
C--
      RETURN
      END


      SUBROUTINE OUTEST (iunit,fout)

C--   Aux. routine OUTEST : write one field on a test output file 

      include "atparam.h"

      real*4 r4out(ix,il)

      do j=1,il
        do i=1,ix
          r4out(i,j)=fout(i,j)
        enddo
      enddo

      write (iunit) r4out

C--
      RETURN
      END

      SUBROUTINE FLUXCORR (FMON)
 
C--   Routine for Flux correction
 
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_tsteps.h"
      include "com_date.h"

      include "com_forcing.h"

      include "com_for_sea.h"

C--   1. Routine for Flux correction  

C     1.1 Calculation (accumulation) of mean ocean temperature


      IYEAR1871 = 1871
      NDAYPM = 30 

      IF(IFLUXCORR .GE. 2) THEN
   
C      IF(IYEAR .GT. IYEAR0 .AND. IYEAR .LT. IYEAR0+11) THEN
       IF(IYEAR .GT. IYEAR0 .OR. IFLUXCORR .EQ. 3) THEN
C       IF(IYEAR .GT. IYEAR1871) THEN
         DO J=1,NLAT
          DO I=1,NLON
             TMP(I,J,IMONTH)=TMP(I,J,IMONTH)+SST1(I,J)/NDAYPM 
           ENDDO
         ENDDO   
C
         IF(IMONTH .EQ. 12 .AND. IDAY .EQ. 30) THEN 
           ICOUNT=ICOUNT+1
           COUNT = FLOAT(ICOUNT-1)/FLOAT(ICOUNT)
           RICOUNT=1./FLOAT(ICOUNT)
c jk
          print*,'  absolute, start, relative year, ICOUNT, COUNT, RICOU
     +NT= ',iyear,iyear0,iyear-iyear0,ICOUNT,COUNT,RICOUNT
c
           DO IT=1,12 
             DO J=1,NLAT
               DO I=1,NLON
                  TMOC(I,J,IT)=COUNT*TMOC(I,J,IT)+
     &                         RICOUNT*TMP(I,J,IT) 
                  TMP(I,J,IT)=0.
               ENDDO
             ENDDO   
           ENDDO
          ENDIF
        ENDIF

       ENDIF

C     1.2 Non-linear, mean-conserving interpolation of mean sea temperature and climatology

      CALL FORIN5 (NGP,IMONTH,FMON,TMOC,TMOC1)
      CALL FORIN5 (NGP,IMONTH,FMON,SST12A,SST01A)

C     1.3 Correction of SSTs

      DO J=1,NLAT
        DO I=1,NLON
             SSTM1(I,J)=SST1(I,J) 
        ENDDO
      ENDDO 


cjk
c      IF(IMONTH .EQ. 1 .AND. IDAY .EQ. 1) print*,'FLUXCOR2= ',IFLUXCORR
c      IF(IMONTH .EQ. 1 .AND. IDAY .EQ. 1) print*,'   IYEAR= ',IYEAR
c      IF(IMONTH .EQ. 1 .AND. IDAY .EQ. 1) print*,'  IYEAR0= ',IYEAR0
c

      IF(IFLUXCORR .EQ. 2.AND.IYEAR.LE.IYEAR0+1) THEN
c      IF(IFLUXCORR .GE. 2.AND.IYEAR.LE.IYEAR1871+1) THEN
        DO J=1,NLAT
           DO I=1,NLON
                IF(OICE1(I,J) .EQ. 0.) THEN
                   SST1(I,J)=SST01A(I,J)
                ENDIF
           ENDDO
         ENDDO 
      ELSE
C      IF(IYEAR .GT. IYEAR0+1) THEN
c       IF(IFLUXCORR .EQ. 2) THEN
c         ascal=(IYEAR+TYEAR-IYEAR0)*12./NMONTS
c       ELSE
         ascal=1.
c       ENDIF
C
c         IF(IMONTH .EQ. 1 .AND. IDAY .LE. 3) then
c          print*,'  absolute, start, relative year, ascal= '
c     +,iyear,iyear0,iyear-iyear0,ascal
c         endif
C
         DO J=1,NLAT
           DO I=1,NLON
              IF(OICE1(I,J) .EQ. 0.) THEN
                 SST1(I,J)=SST1(I,J)+ascal*(SST01A(I,J)-TMOC1(I,J)) 
C                SST1(I,J)=SST01A(I,J)+ascal*(SST1(I,J)-TMOC1(I,J))
              ENDIF 
           ENDDO
         ENDDO 

      ENDIF 

C--
      RETURN
      END   

      SUBROUTINE OBS_SSTA (IDAYH,FMON)

C--   Routine for getting a new SSTA and merging with SST
 
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_tsteps.h"
      include "com_date.h"

      include "com_forcing.h"

      include "com_for_sea.h"
      include "com_dyncon1.h"

      real sstan4(ix,il)
      real*4 r4inp(ix,il)
    
C-- Read new observed SSTA     

      IF (IDAY.EQ.IDAYH) THEN

           do j = 1,il
            do i = 1,ix
              sstan3(i,j,1) = sstan3(i,j,2)
            enddo
           enddo

           read (40,end=100) ((r4inp(i,j),i=1,ix),j=il,1,-1)

           do j = 1,il
              do i = 1,ix
               sstan3(i,j,2) = r4inp(i,j)
             enddo
           enddo 

        ENDIF

 100    continue

C-- Interpolate observed SSTA to daily values

       IF (IDAY.LT.IDAYH) THEN
         CALL FORINT (NGP,2,FMON,SSTAN3,SSTAN4)
       ELSE
         CALL FORINT (NGP,1,FMON,SSTAN3,SSTAN4)
       ENDIF

C     1.2 Non-linear, mean-conserving interpolation of sea  climatology

      CALL FORIN5 (NGP,IMONTH,FMON,SST12A,SST01A)

C-- Merge observed SSTA with modeled 

cfk      nwest=1+nint(nlon*190./360.)
cfk      neast=1+nint(nlon*280./360.)
cfk      rad15=asin(0.2588)

cfk      nwest=1+nint(nlon*140./360.)
cfk      neast=1+nint(nlon*280./360.)
cfk      rad30=asin(0.5)

cfk      do j=1,nlat
cfk        if (abs(radang(j)).lt.rad30) then  
cfk          wob = sqrt(cos(6.*radang(j))) 
cfk          wob = sqrt(cos(3.*radang(j)))        
cfk         do i=nwest,neast
cfk             SST1(I,J) = SST1(I,J) + 
cfk     &                   wob*(SSTAN4(I,J)+SST01A(I,J)-SST1(I,J))   
cfk          enddo
cfk        endif
cfk      enddo

c
c    Observed SST everywhere
c
c      do j=1,nlat
c         do i=1,nlon
c            SST1(I,J) = SST1(I,J) + 
c     &                 (SSTAN4(I,J)+SST01A(I,J)-SST1(I,J))
c         enddo
c      enddo 
c

c
c     Observed SST everwhere outside Indian Ocean
c
      neast=1+nint(nlon*30./360.)
      nwest=1+nint(nlon*136./360.)
      rad30=asin(0.5)

      nlat1=15
      nlat2=34
      do j=1,nlat
         if(j .le. nlat1 .or. j .ge. nlat2) then 
          do i=1,nlon
             SST1(I,J) = SST1(I,J) + 
     &                   (SSTAN4(I,J)+SST01A(I,J)-SST1(I,J))
          enddo 
         else        
           do i=1,neast
              SST1(I,J) = SST1(I,J) + 
     &                    (SSTAN4(I,J)+SST01A(I,J)-SST1(I,J))   
           enddo

           do i=nwest,nlon
              SST1(I,J) = SST1(I,J) + 
     &                    (SSTAN4(I,J)+SST01A(I,J)-SST1(I,J))   
           enddo
        endif
      enddo
c
c
c     set SST in North Atlantic Ocean to Climat.
c

c      nlat2=0
c      dlon=360./ix
c      reast=30.
c      rwest=295.-dlon
c      do j=nlat2+1,il
c        do i=1,ix
c          rlon=(i-1)*dlon
c          if (rlon.lt.reast.or.rlon.gt.rwest) SST1(I,J) = SST01A(I,J)
c        enddo
c      enddo





C--
      RETURN
      END  


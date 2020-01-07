
      SUBROUTINE SETCTL_DAILY (IUNIT,NLON,NLAT,NLEV,NTM,NDAYTM,I3D,N3D,
     *                   N2D,RLAT,RLEV,NAME,NORUN,IYEAR0,IMONT0)
C--
C--   Aux. routine SETCTL_DAILY : write descriptor (.ctl) output file 
C--
C     Post-processing arrays (daily means)
      include "par_tmean.h"

C     Output flags
      include "com_lflags_dmout.h"
C
      CHARACTER*80 LINE(10), LN3D(11), LN2D(33)
C      CHARACTER*80 LINE(10), LN3D(11), LN2D(21)
      CHARACTER*4  LMON(12)
      CHARACTER*5  NAME
      CHARACTER*3  NORUN
      CHARACTER*12 CTLNAME
      INTEGER ILEV(30), NCOUNT
      REAL RLAT(NLAT), RLEV(NLEV)
C
C *** 1. Initialization
C
      DATA LMON /'1jan','1feb','1mar','1apr','1may','1jun',
     &           '1jul','1aug','1sep','1oct','1nov','1dec'/

      DATA LN3D/
     &     'U_850      0  99  zonal (u) wind at 850 hPa       [m/s]',
     &     'V_850      0  99  meridional (v) wind at 850 hPa  [m/s]',
     &     'Q_850      0  99  specific humidity at 850 hPa   [g/Kg]',
     &     'GH_500     0  99  geopotential height at 500 hPa    [m]',
     &     'U_200      0  99  zonal (u) wind at 200 hPa       [m/s]',
     &     'V_200      0  99  meridional (v) wind at 200 hPa  [m/s]',
     &     'H2O_850    0  99  H2O specific hum. at 850 hPa   [g/Kg]',
     &     'HDO_850    0  99  HDO specific hum. at 850 hPa   [g/Kg]',
     &     'H218O_850  0  99  H218O specific hum. at 850 hPa [g/Kg]',
     &      2*' '/

      DATA LN2D/
     &     'PRECLS     0  99  large-scale precipitation    [mm/day]',
     &     'PRLH2O     0  99  H2O large-scale prec         [mm/day]',
     &     'PRLHDO     0  99  HDO large-scale prec         [mm/day]',
     &     'PRLH218O   0  99  H218O large-scale prec       [mm/day]',
     &     'PRECNV     0  99  convective precipitation     [mm/day]',
     &     'PRCH2O     0  99  H2O convective prec          [mm/day]',
     &     'PRCHDO     0  99  HDO convective prec          [mm/day]',
     &     'PRCH218O   0  99  H218O convective prec        [mm/day]',
     &     'EVAP       0  99  evaporation                  [mm/day]',
     &     'EVAPH2O    0  99  H2O evaporation              [mm/day]',
     &     'EVAPHDO    0  99  HDO evaporation              [mm/day]',
     &     'EVAPH218O  0  99  H218O evaporation            [mm/day]',
     &     'WTWG       0  99  Groundwater                  [m^3/m^3]',
     &     'WTWGH2O    0  99  H2O Groundwater              [m^3/m^3]',
     &     'WTWGHDO    0  99  HDO Groundwater              [m^3/m^3]',
     &     'WTWGH218O  0  99  H218O Groundwater            [m^3/m^3]',
     &     'WTWB       0  99  Bulk Water                   [m^3/m^3]',
     &     'WTWBH2O    0  99  H2O Bulk Water               [m^3/m^3]',
     &     'WTWBDO     0  99  HDO Bulk Water               [m^3/m^3]',
     &     'WTWBH218O  0  99  H218O Bulk Water             [m^3/m^3]',
     &     'WTRUN      0  99  Runoff                   	   [mm/day]',
     &     'WTRH2O     0  99  H2O Runoff               	   [mm/day]',
     &     'WTRDO      0  99  HDO Runoff                   [mm/day]',
     &     'WTRH218O   0  99  H218O Runoff                 [mm/day]',
     &     'OLR        0  99  outgoing longwave rad.  (uw.) [W/m^2]',
     &     'USTR       0  99  u-stress                (uw.) [N/m^2]',
     &     'VSTR       0  99  v-stress                (uw.) [N/m^2]',
     &     'LSHF       0  99  heat flux into land sfc (dw.) [W/m^2]',
     &     'SSHF       0  99  heat flux into  sea sfc (dw.) [W/m^2]',
     &     'MSLP       0  99  mean-sea-level pressure         [hPa]',
     &     'TEMP0      0  99  near-surface air temperature   [degK]',
     &      2*' '/

      LINE( 1)='DSET   ^attmdxxx_%y4.grd'
      LINE( 2)='TITLE   Daily Means from PE5L run no. xxx'                 
      LINE( 3)='UNDEF   9.999E+19'
      LINE( 4)='OPTIONS SEQUENTIAL TEMPLATE BIG_ENDIAN 365_day_calendar'
      LINE( 5)='XDEF     nnn  LINEAR     0.000     x.xxx'
      LINE( 6)='YDEF     nnn  LEVELS'
      LINE( 7)='ZDEF      nn  LEVELS       950'
      LINE( 8)='TDEF  nnnnnn  LINEAR            1jan1900      nndy'
      LINE( 9)='VARS      nn'
      LINE(10)='ENDVARS'

      CTLNAME=NAME//NORUN//'.ctl'
      OPEN ( UNIT=IUNIT, FILE=CTLNAME, FORM='FORMATTED' )
C
      C1=90./ASIN(1.)
C
      DO 120 K=1,NLEV
      ILEV(K)=NINT(1000.*RLEV(K))
  120 CONTINUE
C
C *** 2. Insert parameters in strings
C
      LINE(1)( 9:13)= NAME(1:5)
      LINE(1)(14:16)=NORUN(1:3)
      LINE(2)(39:41)=NORUN(1:3)
C
      WRITE (LINE(5)(10:12),'(I3)') NLON
      WRITE (LINE(5)(31:40),'(F10.3)') (360./NLON)
      WRITE (LINE(6)(10:12),'(I3)') NLAT
      WRITE (LINE(7)(10:12),'(I3)') NLEV
C
      WRITE (LINE(8) (7:12),'(I6)') NTM
      WRITE (LINE(8)(37:40),'(I4)') IYEAR0
      LINE(8)(33:36)=LMON(IMONT0)(1:4)
      IF (MOD(NDAYTM,30).EQ.0) THEN
        WRITE (LINE(8)(46:48),'(I3)') NDAYTM/30
        LINE(8)(49:50)='mo'
      ELSE
        WRITE (LINE(8)(46:49),'(I4)') 1460
        LINE(8)(50:51)='mn'
      ENDIF
C
      NCOUNT=0
C      DO 210 N=1,N3D+N2D-1
      DO 210 N=1,N3D+N2D
         IF(LDOUT(N)) NCOUNT = NCOUNT+1
  210 CONTINUE   
      WRITE (LINE(9)(11:12),'(I2)') NCOUNT
C
C *** 3. Write ASCII control file 
C
      DO 310 JLINE=1,5
        WRITE (IUNIT,1000) LINE(JLINE)
  310 CONTINUE
C
      WRITE (IUNIT,1010) LINE(6)(1:20), (C1*RLAT(J),J=1,NLAT)
      WRITE (IUNIT,1020) LINE(7)(1:20), (ILEV(K),K=NLEV,1,-1)

C
      DO 320 JLINE=8,9
        WRITE (IUNIT,1000) LINE(JLINE)
  320 CONTINUE
      
      DO 330 JLINE=1,N3D
        IF(LDOUT(JLINE)) THEN
           WRITE (IUNIT,1000) LN3D(JLINE)
        ENDIF
  330 CONTINUE

C      DO 340 JLINE=1,N2D-1
      DO 340 JLINE=1,N2D
        IF(LDOUT(N3D+JLINE)) THEN
           WRITE (IUNIT,1000) LN2D(JLINE)
        ENDIF
  340 CONTINUE
C
      WRITE (IUNIT,1000) LINE(10)
C
      CLOSE ( UNIT=IUNIT )
C
 1000 FORMAT (A80)
 1010 FORMAT (A20,6F10.3/(8F10.3))
 1020 FORMAT (A20,10I6)
C
      RETURN
      END 

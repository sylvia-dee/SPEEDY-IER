
      SUBROUTINE SETCTL (IUNIT,NLON,NLAT,NLEV,NTM,NDAYTM,I3D,N3D,N2D,
     *                   N2D_D,RLAT,RLEV,NAME,NORUN,IYEAR0,IMONT0)
C--
C--   Aux. routine SETCTL : write descriptor (.ctl) output file 
C--
C      CHARACTER*80 LINE(10), LN3D(33), LN2D(31), LN2D_D(19)
      CHARACTER*80 LINE(10), LN3D(33), LN2D(31), LN2D_D(31)
      CHARACTER*4  LMON(12), NAME
      CHARACTER*3  NORUN
      CHARACTER*11 CTLNAME
      INTEGER ILEV(30)
      REAL RLAT(NLAT), RLEV(NLEV)
C
C *** 1. Initialization
C
      DATA LMON /'1jan','1feb','1mar','1apr','1may','1jun',
     &           '1jul','1aug','1sep','1oct','1nov','1dec'/

      DATA LN3D/
     &     'GH         n  99  geopotential height               [m]',
     &     'TEMP       n  99  abs. temperature               [degK]',
     &     'U          n  99  zonal (u) wind                  [m/s]',
     &     'V          n  99  meridional (v) wind             [m/s]',
     &     'Q          n  99  specific humidity              [g/Kg]',
     &     'RH         n  99  relative humidity                 [%]',
     &     'H2O        n  99  H2O tracer specific humidity   [g/Kg]',
     &     'HDO        n  99  HDO tracer specific humidity   [g/Kg]',
     &     'H218O      n  99  H218O tracer specific humidity [g/Kg]',
     &     'OMEGA      n  99  pressure vertical velocity     [Pa/s]',
     &     'PSI        n  99  streamfunction           [10^6 m^2/s]',
     &     'CHI        n  99  velocity potential       [10^6 m^2/s]',
     
     &     'VARGH      n  99  variance of geop. height        [m^2]',
     &     'VART       n  99  variance of temperature      [degK^2]',
     &     'VARU       n  99  variance of u-wind             [J/Kg]',
     &     'VARV       n  99  variance of v-wind             [J/Kg]',
     &     "COVUV      n  99  u'v' covariance (trans.)       [J/Kg]",
     &     "COVVT      n  99  v'T' covariance (trans.)   [degK m/s]",
 
     &     'DTLSC      n  99  dT/dt by large-scale cond. [degK/day]',
     &     'DTCNV      n  99  dT/dt by convection        [degK/day]',
     &     'DTRSW      n  99  dT/dt by shortwave rad.    [degK/day]',
     &     'DTRLW      n  99  dT/dt by longwave  rad.    [degK/day]',
     &     'DTPBL      n  99  dT/dt by PBL processes     [degK/day]',
     &      10*' '/

      DATA LN2D/
     &     'SP         0  99  surface pressure                [hPa]',
     &     'MSLP       0  99  mean-sea-level pressure         [hPa]',
     &     'ST         0  99  surface temperature            [degK]',
     &     'SKINT      0  99  skin temperature               [degK]',
     &     'SWAV       0  99  soil wetness availability         [%]',
     &     'ALB        0  99  surface albedo                    [%]',
     &     'U0         0  99  near-surface u-wind             [m/s]',
     &     'V0         0  99  near-surface v-wind             [m/s]',
     &     'TEMP0      0  99  near-surface air temperature   [degK]',
     &     'RH0        0  99  near-surface relative humidity    [%]',
     &     'CLC        0  99  cloud cover (total)               [%]',
     &     'CLTOP      0  99  pressure at cloud top           [hPa]',
     &     'IPTOP      0  99  highest precipitation level index  []',
     &     'SHF        0  99  sensible heat flux      (uw.) [W/m^2]',
     &     'TSR        0  99  top shortwave rad.      (dw.) [W/m^2]',
     &     'SSR        0  99  surface shortwave rad.  (dw.) [W/m^2]',
     &     'SLR        0  99  surface longwave rad.   (uw.) [W/m^2]',
cfk     &     'LSTA       0  99  land-surface temp. anomaly     [degK]',
cfk     &     'SSTA       0  99   sea-surface temp. anomaly     [degK]',
     &     'LST        0  99  land-surface temp.             [degK]',
     &     'SST        0  99   sea-surface temp.             [degK]',
     &     'SSTM       0  99  ocean model sea-surface temp.  [degK]',
     &      11*' '/

      DATA LN2D_D/
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
     &      2*' '/

      LINE( 1)='DSET   ^attmxxx_%y4.grd'
      LINE( 2)='TITLE   Means/variances from PE5L run no. xxx'                 
      LINE( 3)='UNDEF   9.999E+19'
      LINE( 4)='OPTIONS   SEQUENTIAL  TEMPLATE nnnnnnnnnnnnnnnn'
      LINE( 5)='XDEF     nnn  LINEAR     0.000     x.xxx'
      LINE( 6)='YDEF     nnn  LEVELS'
      LINE( 7)='ZDEF      nn  LEVELS       950'
      LINE( 8)='TDEF   nnnnn  LINEAR            1jan1900     nnndy'
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
      LINE(1)( 9:12)= NAME(1:4)
      LINE(1)(13:15)=NORUN(1:3)
      LINE(2)(43:45)=NORUN(1:3)
C
      IF (MOD(NDAYTM,30).EQ.0) THEN
        LINE(4)(32:60)='   big_endian                '
      ELSE
        LINE(4)(32:60)='365_day_calendar   big_endian'
      ENDIF
C
      WRITE (LINE(5)(10:12),'(I3)') NLON
      WRITE (LINE(5)(31:40),'(F10.3)') (360./NLON)
      WRITE (LINE(6)(10:12),'(I3)') NLAT
      WRITE (LINE(7)(10:12),'(I3)') NLEV
C
      WRITE (LINE(8) (9:12),'(I4)') NTM
      WRITE (LINE(8)(37:40),'(I4)') IYEAR0
      LINE(8)(33:36)=LMON(IMONT0)(1:4)
      IF (MOD(NDAYTM,30).EQ.0) THEN
        WRITE (LINE(8)(46:48),'(I3)') NDAYTM/30
        LINE(8)(49:50)='mo'
      ELSE
        WRITE (LINE(8)(46:50),'(I5)') NDAYTM*1460
        LINE(8)(51:52)='mn'
      ENDIF
C
      WRITE (LINE(9)(11:12),'(I2)') N3D+N2D+N2D_D
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
C
      DO 330 JLINE=I3D,I3D+N3D-1
        WRITE (LN3D(JLINE)(10:12),'(I3)') NLEV
        WRITE (IUNIT,1000) LN3D(JLINE)
  330 CONTINUE
C
      DO 340 JLINE=1,N2D
        WRITE (IUNIT,1000) LN2D(JLINE)
  340 CONTINUE

      DO 350 JLINE=1,N2D_D
        WRITE (IUNIT,1000) LN2D_D(JLINE)
  350 CONTINUE

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

C--
C--   /PHYGR1/ : Model variables on gaussian grid (updated in PHYPAR)
C--    UG1    = u-wind
C--    VG1    = v-wind
C--    TG1    = abs. temperature
C--    QG1    = specific humidity (g/kg)
C--    PHIG1  = geopotential
C--    PSLG1  = log. of surface pressure

      COMMON /PHYGR1/ UG1(NGP,NLEV), VG1(NGP,NLEV), TG1(NGP,NLEV),
     &                QG1(NGP,NLEV,NTR), PHIG1(NGP,NLEV), PSLG1(NGP)

C--   
C--   /PHYGR2/ : Diagnosed upper-air variables (updated in PHYPAR)
C--    SE     = dry static energy
C--    RH     = relative humidity
C--    QSAT   = saturation specific humidity (g/kg)

      COMMON /PHYGR2/ SE(NGP,NLEV), RH(NGP,NLEV), QSAT(NGP,NLEV)

C--
C--   /PHYGR3/ : Diagnosed surface variables (updated in PHYPAR)
C--    PSG    = surface pressure
C--    TS     = surface temperature
C--    TSKIN  = skin temperature
C--    U0     = near-surface u-wind
C--    V0     = near-surface v-wind
C--    T0     = near-surface air temperature
C--    Q0     = near-surface specific humidity (g/kg)
C--    CLOUDC = total cloud cover (fraction)
C--    CLTOP  = norm. pressure at cloud top
C--    PRTOP  = top of precipitation (level index)
	

      COMMON /PHYGR3/ PSG(NGP), TS(NGP), TSKIN(NGP),
     &                U0(NGP), V0(NGP), T0(NGP), Q0(NGP,NTR),
     &                CLOUDC(NGP), CLTOP(NGP), PRTOP(NGP)

C--
C--   /PHYTEN/ : Physical param. tendencies (updated in PHYPAR)
C--    TT_CNV  =  temperature tendency due to convection
C--    QT_CNV  = sp. humidity tendency due to convection
C--    TT_LSC  =  temperature tendency due to large-scale condensation
C--    QT_LSC  = sp. humidity tendency due to large-scale condensation
C--    TT_RSW  =  temperature tendency due to short-wave radiation
C--    TT_RLW  =  temperature tendency due to long-wave radiation
C--    UT_PBL  =       u-wind tendency due to PBL and diffusive processes
C--    VT_PBL  =       v-wind tendency due to PBL and diffusive processes
C--    TT_PBL  =  temperature tendency due to PBL and diffusive processes
C--    QT_PBL  = sp. humidity tendency due to PBL and diffusive processes

      COMMON /PHYTEN/ TT_CNV(NGP,NLEV), QT_CNV(NGP,NLEV,NTR),
     &                TT_LSC(NGP,NLEV), QT_LSC(NGP,NLEV,NTR),
     &                TT_RSW(NGP,NLEV), TT_RLW(NGP,NLEV),
     &                UT_PBL(NGP,NLEV), VT_PBL(NGP,NLEV),
     &                TT_PBL(NGP,NLEV), QT_PBL(NGP,NLEV,NTR)

C--
C--   /FLUXES/ : Surface and upper boundary fluxes (updated in PHYPAR)
C--    PRECNV = convective precipitation  [g/(m^2 s)], total
C--    PRECLS = large-scale precipitation [g/(m^2 s)], total
C--    SNOWCV = convective precipitation  [g/(m^2 s)], snow only
C--    SNOWLS = large-scale precipitation [g/(m^2 s)], snow only
C--    CBMF   = cloud-base mass flux 
C--    TSR    = top-of-atm. shortwave radiation (downward)
C--    SSR    = surface shortwave radiation (downward)
C--    SLR    = surface longwave radiation (upward) 
C--    OLR    = outgoing longwave radiation (upward)
C--    USTR   = u-stress (1: land, 2: sea, 3: weighted average)
C--    VSTR   = v-stress (1: land, 2: sea, 3: weighted average)
C--    SHF    = sensible heat flux (1: land, 2: sea, 3: w. average)
C--    EVAP   = evaporation [g/(m^2 s)] (1: land, 2: sea, 3: w. average)
C--    HFLUXN = net heat flux into the surface (1: land, 2: sea)
C--    WTWB   = bulk water 	[m^3/m^3], land model
C--    WTWG   = ground water 	[m^3/m^3], land model
C--    WTRUN  = runoff 	            [g/(m^2 s)], land model
C--    WTDRN  = subsurface drainage [g/(m^2 s)], land model
C--    RLD    = Isotope ratio of soil water
C--    RB     = Isotope ratio of bulk water (deep bucket)
C--    SOILW  = Soil water content (ratio 0 to 1)
C--    EMP    = Evap > Precip flag [g/(m^2 s)], land model
C--    ALPLND = equilibrium fractionation factor, land model
C--    HUMLND = q/qsat land model passed from suflux
C--    TLND   = TSKIN to pass to land model from suflux

      COMMON /FLUXES/ PRECNV(NGP,NTR), PRECLS(NGP,NTR),
     &                SNOWCV(NGP), SNOWLS(NGP),
     &                CBMF(NGP), TSR(NGP), SSR(NGP), SLR(NGP), OLR(NGP),
     &                USTR(NGP,3), VSTR(NGP,3), 
     &                SHF(NGP,3), EVAP(NGP,3,NTR), HFLUXN(NGP,2),
     &                WTWB(NGP,NTR), WTWG(NGP,NTR), 
     &                WTRUN(NGP,NTR), WTDRN(NGP,NTR), RLD(NGP,NTR),
     &                RB(NGP,NTR),SOILW(NGP),EMP(NGP),ALPLND(NGP,NTR),
     &                HUMLND(NGP,NTR),TLND(NGP)

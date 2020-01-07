      SUBROUTINE PHYPAR (VOR1,DIV1,T1,Q1,PHI1,PSL1,
     &                   UTEND,VTEND,TTEND,QTEND)
C--
C--   SUBROUTINE PHYPAR (VOR1,DIV1,T1,Q1,PHI1,PSL1,
C--  &                   UTEND,VTEND,TTEND,QTEND)
C--
C--   Purpose: compute physical parametrization tendencies for U, V, T, Q 
C--   and add them to dynamical grid-point tendencies
C--   Input-only  arguments:   VOR1   : vorticity (sp)
C--                            DIV1   : divergence (sp)
C--                            T1     : temperature (sp)
C--                            Q1     : specific humidity (sp)
C--                            PHI1   : geopotential (sp)
C--                            PSL1   : log of sfc pressure (sp)
C--   Input-output arguments:  UTEND  : u-wind tendency (gp)
C--                            VTEND  : v-wind tendency (gp)
C--                            TTEND  : temp. tendency (gp)
C--                            QTEND  : spec. hum. tendency (gp)
C--   Modified common blocks:  PHYGR1, PHYGR2, PHYGR3, PHYTEN, FLUXES
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Constants + functions of sigma and latitude
      include "com_physcon.h"

C     Model variables, tendencies and fluxes on gaussian grid
      include "com_physvar.h"

C     Surface forcing fields (time-inv. or functions of seasonal cycle)
      include "com_forcing.h"
      include "com_for_sea.h"
      include "com_for_land.h"

C     Logical flags
      include "com_lflags.h"

C     Time parameters
      include "com_tsteps.h"

C     Water isotope tracer variables/constants
      include "com_isocon.h"

      COMPLEX VOR1(MX,NX,NLEV), DIV1(MX,NX,NLEV), T1(MX,NX,NLEV),
     &          Q1(MX,NX,NLEV,NTR), PHI1(MX,NX,NLEV), PSL1(MX,NX),
     &          UCOS(MX,NX), VCOS(MX,NX)

      REAL UTEND(NGP,NLEV), VTEND(NGP,NLEV), TTEND(NGP,NLEV),
     &     QTEND(NGP,NLEV,NTR)

      INTEGER IPTOP(NGP), ICLTOP(NGP,2), ICNV(NGP)
      REAL    RPS(NGP), ST4S(NGP), SLRU(NGP,3), SLR_D(NGP), SLRU_3(NGP)

      REAL QT_RSC(NGP,NLEV,NTR)

      iitest=0

C--   1. Compute grid-point fields

C     1.1 Convert model spectral variables to grid-point variables

      if (iitest.eq.1) print *, ' 1.1 in PHYPAR'

      DO K=1,NLEV

        CALL UVSPEC (VOR1(1,1,K),DIV1(1,1,K),UCOS,VCOS)
        CALL GRID   (UCOS,UG1(1,K),2)
        CALL GRID   (VCOS,VG1(1,K),2)

      ENDDO

      DO K=1,NLEV
        CALL GRID   (T1(1,1,K),  TG1(1,K),  1)
        CALL GRID   (PHI1(1,1,K),PHIG1(1,K),1)
      ENDDO

      DO ITR=1,NTR
        DO K=1,NLEV
          CALL GRID   (Q1(1,1,K,ITR),  QG1(1,K,ITR),  1)
        ENDDO
      ENDDO

      CALL GRID (PSL1,PSLG1,1)

C     Check that tracers are correct at start of routine
      CALL CHKQ(QG1,qmag_chk,'phypar1')

C     Remove negative humidity values
C     CALL QNEG (QG1)
      DO ITR=1,NTR
        DO K=1,NLEV
          DO J=1,NGP
	    qg1(j,k,ITR)=max(qg1(j,k,ITR),0.)
          ENDDO
        ENDDO
      ENDDO

C     1.2 Compute thermodynamic variables

      if (iitest.eq.1) print *, ' 1.2 in PHYPAR'

      DO J=1,NGP
       PSG(J)=EXP(PSLG1(J))
       RPS(J)=1./PSG(J)
      ENDDO

      DO K=1,NLEV
       DO J=1,NGP
        SE(J,K)=CP*TG1(J,K)+PHIG1(J,K)
       ENDDO
      ENDDO

      DO K=1,NLEV
       CALL SHTORH (1,NGP,TG1(1,K),PSG,SIG(K),QG1(1,K,1),
     &              RH(1,K),QSAT(1,K))
      ENDDO


C--   2. Precipitation 


C     2.1 Deep convection

cfk#if !defined(KNMI)

      CALL CONVMF (PSG,SE,TG1,QG1,QSAT,
     &             IPTOP,CBMF,PRECNV,TT_CNV,QT_CNV)
cfk#else
cfk      CALL CONVMF (PSG,SE,QG1,QSAT,TS,
cfk     &             IPTOP,CBMF,PRECNV,SNOWCV,TT_CNV,QT_CNV)



      DO K=2,NLEV
       DO J=1,NGP
        TT_CNV(J,K)=TT_CNV(J,K)*RPS(J)*GRDSCP(K)
       ENDDO
      ENDDO
      DO ITR=1,NTR
        DO K=2,NLEV
         DO J=1,NGP
          QT_CNV(J,K,ITR)=QT_CNV(J,K,ITR)*RPS(J)*GRDSIG(K)
         ENDDO
       ENDDO
      ENDDO

      DO J=1,NGP
        ICNV(J)=NLEV-IPTOP(J)
      ENDDO


C     Check that tracer tendencies are correct 
      CALL CHKQ(QT_CNV,qtmag_chk,'xCNVtend')

C     2.2 Large-scale condensation


      CALL LSCOND (PSG,TG1,QG1,QSAT,
     &             IPTOP,PRECLS,TT_LSC,QT_LSC)


C      DO ITR=1,NTR
C        DO J=1,NGP
C          if (EVAP(J,1,1) .ne. 0.0 .and. J .eq. 2000) then
C          write(*,*) "afLS_evap", ((EVAP(J,1,ITR)/EVAP(J,1,1))-1)
C     &   	*1000,EVAP(J,1,ITR),RLD(J,ITR),ITR,J
C          endif
C        ENDDO
C      ENDDO
C     Check that tracer tendencies are correct 
      CALL CHKQ(QT_LSC,qtmag_chk,'xLSCtend')

      DO K=2,NLEV
       DO J=1,NGP
        TTEND(J,K)=TTEND(J,K)+TT_CNV(J,K)+TT_LSC(J,K)
       ENDDO
      ENDDO

      DO ITR=1,NTR
       DO K=1,NLEV
        DO J=1,NGP
         QTEND(J,K,ITR)=QTEND(J,K,ITR) 
     &                    +QT_CNV(J,K,ITR)+QT_LSC(J,K,ITR)
        ENDDO
       ENDDO
      ENDDO



C--   3. Radiation (shortwave and longwave) and surface fluxes

C     3.1 Compute shortwave tendencies and initialize lw transmissivity

      if (iitest.eq.1) print *, ' 3.1 in PHYPAR'

C     The sw radiation may be called at selected time steps

      IF (LRADSW) THEN

        CALL CLOUD (QG1,RH,PRECNV,PRECLS,IPTOP,
     &              ICLTOP,CLOUDC)

        DO J=1,NGP
          CLTOP(J)=SIGH(ICLTOP(J,1)-1)*PSG(J)
          PRTOP(J)=float(IPTOP(J))
        ENDDO

        CALL RADSW (PSG,QG1,ALB1,ICLTOP,CLOUDC,
     &              TSR,SSR,TT_RSW)

        DO K=1,NLEV
         DO J=1,NGP
          TT_RSW(J,K)=TT_RSW(J,K)*RPS(J)*GRDSCP(K)
         ENDDO
        ENDDO

      ENDIF

C     3.2 Compute downward longwave fluxes 

      CALL RADLW (-1,TG1,TS,ST4S,
     &            OLR,SLR,TT_RLW)

C     3.3. Compute surface fluxes and land skin temperature


C sdee check QA values at various points

C      DO ITR=2,NTR
C        DO K=1,NLEV
C          DO J = 1, NGP
CC             if (J .eq. 2200) then
C               write(*,*) "QA pre land/suflux=", QG1(J,K,ITR), ITR
C             endif
C          ENDDO
C        ENDDO
C      ENDDO

C     3.3.0 Run Land Model to conserve isotope tracers over land during evaporation.
 

      CALL SUFLUX (PSG,UG1,VG1,TG1,QG1,RH,PHIG1,
     &             PHIS0,FMASK1,STL1,SST1,SOILW1,SSR,SLR,
     &             USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
     &             TS,TSKIN,U0,V0,T0,Q0,.true.,RLD,RB,SOILW)

      CALL CULND(WTWB,WTWG,WTRUN,WTDRN,PRECNV,PRECLS,EVAP,RLD,RB,
     &           FMASK1,SOILW)


C     3.3.1 save some surface fluxes for later

      DO J=1,NGP
        SLR_D(J)=SLR(J)
        SLRU_3(J)=SLRU(J,3)
      ENDDO     
 
cfk#if !defined(KNMI)
cfk      CALL ADDFLX (HFLUXN,USTR,VSTR,U0,V0,SHF, 
cfk     &              EVAP,SSR,SLR,SLRU, 
cfk     &              PRECNV,PRECLS,T0)
cfk#else
cfk      CALL ADDFLX (HFLUXN,USTR,VSTR,U0,V0,SHF, 
cfk     &              EVAP,SSR,SLR,SLRU, 
cfk     &              PRECNV,PRECLS,SNOWCV,SNOWLS,T0)
cfk#endif

C     3.4 Compute upward longwave fluxes, convert them to tendencies 
C         and add shortwave tendencies

      if (iitest.eq.1) print *, ' 3.4 in PHYPAR'

      CALL RADLW (1,TG1,TS,SLRU(1,3),
     &            OLR,SLR,TT_RLW)

      DO K=1,NLEV
       DO J=1,NGP
        TT_RLW(J,K)=TT_RLW(J,K)*RPS(J)*GRDSCP(K)
        TTEND (J,K)=TTEND(J,K)+TT_RSW(J,K)+TT_RLW(J,K)
       ENDDO
      ENDDO

C--   4. PBL interactions with lower troposphere

C     4.1 Vertical diffusion and shallow convection

      CALL VDIFSC (UG1,VG1,SE,RH,QG1,QSAT,PHIG1,ICNV,
     &             UT_PBL,VT_PBL,TT_PBL,QT_PBL)


C     Check that tracer tendencies are correct 
      CALL CHKQ(QT_PBL,qtmag_chk,'xPBLtend')

C     4.2 Add tendencies due to surface fluxes 

      DO J=1,NGP
       UT_PBL(J,NLEV)=UT_PBL(J,NLEV)+USTR(J,3)*RPS(J)*GRDSIG(NLEV)
       VT_PBL(J,NLEV)=VT_PBL(J,NLEV)+VSTR(J,3)*RPS(J)*GRDSIG(NLEV)
       TT_PBL(J,NLEV)=TT_PBL(J,NLEV)+ SHF(J,3)*RPS(J)*GRDSCP(NLEV)
      ENDDO

      DO ITR=1,NTR
       DO J=1,NGP
        QT_PBL(J,NLEV,ITR)=QT_PBL(J,NLEV,ITR)+ 
     &              EVAP(J,3,ITR)*RPS(J)*GRDSIG(NLEV)
       ENDDO
      ENDDO

C     Check that tracer tendencies are correct 
      CALL CHKQ(QT_PBL,qtmag_chk,'xSFXtend')

      DO K=1,NLEV
       DO J=1,NGP
        UTEND(J,K)=UTEND(J,K)+UT_PBL(J,K)
        VTEND(J,K)=VTEND(J,K)+VT_PBL(J,K)
        TTEND(J,K)=TTEND(J,K)+TT_PBL(J,K)
       ENDDO
      ENDDO

      DO ITR=1,NTR
       DO K=1,NLEV
        DO J=1,NGP
         QTEND(J,K,ITR)=QTEND(J,K,ITR)+QT_PBL(J,K,ITR)
        ENDDO
       ENDDO
      ENDDO

C     Check that tracer tendencies are correct 
      CALL CHKQ(QTEND,qtmag_chk,'xQtend')

C--  5.1 (Optionally) Apply rescaling to isotope tracers

      if (lirescl) then
        CALL RESCALR(QG1, QT_RSC)
        DO ITR=2,NTR
         DO K=1,NLEV
          DO J=1,NGP
           QTEND(J,K,ITR)=QTEND(J,K,ITR)+QT_RSC(J,K,ITR)
          ENDDO
         ENDDO
        ENDDO
      endif


C--  6. Recomputing Sea fluxes in case of flux-correction

      IF (IFLUXCORR .GT. 0) THEN 
         stop 'isotopes not set for this... might work?'

         CALL SUFLUX (PSG,UG1,VG1,TG1,QG1,RH,PHIG1,
     &             PHIS0,FMASK1,STL1,SSTM1,SOILW1,SSR,SLR_D,
     &             USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
     &             TS,TSKIN,U0,V0,T0,Q0,.false.,RLD)
      
      ENDIF


      DO J=1,NGP
c        SLR(J)=SLR_D(J)
        SLRU(J,3)=SLRU_3(J)
      ENDDO     


cfk#if !defined(KNMI)
      CALL ADDFLX (HFLUXN,USTR,VSTR,U0,V0,SHF, 
     &              EVAP,SSR,SLR_D,SLRU, 
     &              PRECNV,PRECLS,T0)
cfk#else
cfk      CALL ADDFLX (HFLUXN,USTR,VSTR,U0,V0,SHF, 
cfk     &              EVAP,SSR,SLR,SLRU, 
cfk     &              PRECNV,PRECLS,SNOWCV,SNOWLS,T0)
cfk#endif


C--   5. Random diabatic forcing 

      IF (LRANDF) THEN

C       5.1 Compute zonal-mean cross sections of diabatic forcing

        IF (LRADSW) THEN
          CALL XS_RDF (TT_LSC,TT_CNV,1)
          CALL XS_RDF (TT_RSW,TT_RLW,2)
        ENDIF

C--     5.2 Compute and store 3-D pattern of random diabatic forcing

        DO K=1,NLEV
         DO J=1,NGP
          TT_CNV(J,K)=TT_CNV(J,K)+TT_LSC(J,K)
         ENDDO
        ENDDO

        CALL SETRDF (TT_LSC)

        DO K=1,NLEV
         DO J=1,NGP
          TTEND(J,K)=TTEND(J,K)+TT_LSC(J,K)
         ENDDO
        ENDDO

      ENDIF


C--
      RETURN
      END

      include "phy_setrdf.f"


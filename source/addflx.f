cfk#if !defined(KNMI)
      SUBROUTINE ADDFLX (HFLUXN,USTR,VSTR,U0,V0,SHF, 
     &                   EVAP,SSR,SLR,SLRU, 
     &                   PRECNV,PRECLS,T0)   
cfk#else
cfk      SUBROUTINE ADDFLX (HFLUXN,USTR,VSTR,U0,V0,SHF, 
cfk     &                   EVAP,SSR,SLR,SLRU, 
cfk     &                   PRECNV,PRECLS,SNOWCV,SNOWLS,T0)   
cfk#endif

C--   SUBROUTINE ADDFLX (HFLUXN,USTR,VSTR,U0,V0,SHF, 
C--   &                  EVAP,SSR,SLR,SLRU, 
C--   &                  PRECNV,PRECLS,T0)   
C--
C--   Purpose: Add up fluxes to provide daily averages 
C--   or integrals to be used in sea- and landmodels
C--


      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX)
      include "com_tsteps.h"
      include "com_forcing.h"
      include "com_physcon.h"
      include "com_atm2sea.h"
      include "com_atm2land.h"

      REAL HFLUXN (IX,IL,2)
      REAL USTR(IX,IL,3), VSTR(IX,IL,3), SHF(IX,IL,3), 
     &     EVAP(IX,IL,3), SLRU(IX,IL,3), U0(IX,IL), 
     &     V0(IX,IL), T0(IX,IL), SSR(IX,IL), SLR(IX,IL), 
     &     PRECNV(IX,IL), PRECLS(IX,IL)
cfk#if defined(KNMI)
cfk      REAL SNOWCV(IX,IL), SNOWLS(IX,IL)
cfk
cfk      include "cls_inphys.h"
cfk#endif

C--   1. Integrated net heat flux 

      do j=1,NLAT
        do i=1,NLON
          hfluxn(i,j,1)=hfluxn(i,j,1)*fmask1(i,j)
          hfluxn(i,j,2)=hfluxn(i,j,2)*fmask0(i,j)
        enddo
      enddo

        do j=1,NLAT
          do i=1,NLON
            hfintl(i,j)=hfintl(i,j)+hfluxn(i,j,1)
          enddo
        enddo

        do j=1,NLAT
          do i=1,NLON
            hfints(i,j)=hfints(i,j)+hfluxn(i,j,2)
          enddo
        enddo

C--    2. Averaged fields and fluxes

C        fac=2.501e3
        fac=ALHC
        rsteps=1./nsteps
        do j=1,NLAT
          do i=1,NLON
            ustrlm(i,j)=ustrlm(i,j)+
     &                  ustr(i,j,1)*rsteps
            vstrlm(i,j)=vstrlm(i,j)+
     &                  vstr(i,j,1)*rsteps
            u0lm(i,j)=u0lm(i,j)+
     &                  u0(i,j)*rsteps
            v0lm(i,j)=v0lm(i,j)+
     &                  v0(i,j)*rsteps
            shflm(i,j)=shflm(i,j)+
     &                  shf(i,j,1)*rsteps
            xlhflm(i,j)=xlhflm(i,j)+
     &                  fac*evap(i,j,1)*rsteps
            ssrlm(i,j)=ssrlm(i,j)+
     &                  ssr(i,j)*rsteps
            slrdlm(i,j)=slrdlm(i,j)+
     &                  slr(i,j)*rsteps
            slrulm(i,j)=slrulm(i,j)+
     &                  slru(i,j,1)*rsteps
            alblm(i,j)=alblm(i,j)+
     &                  alb1(i,j)*rsteps
            preclm(i,j)=preclm(i,j)+
     &                  (precnv(i,j)+precls(i,j))*rsteps
cfk#if !defined(KNMI)
            snowrlm(i,j)=snowrlm(i,j)+
     &                  0.*rsteps
cfk#else
cfk            snowrlm(i,j)=snowrlm(i,j)+
cfk     &                  (snowcv(i,j)+snowls(i,j))*rsteps
cfk#endif
          enddo
        enddo

        do j=1,NLAT
          do i=1,NLON
            ustrsm(i,j)=ustrsm(i,j)+
     &                  ustr(i,j,2)*rsteps
            vstrsm(i,j)=vstrsm(i,j)+
     &                  vstr(i,j,2)*rsteps
            u0sm(i,j)=u0sm(i,j)+
     &                  u0(i,j)*rsteps
            v0sm(i,j)=v0sm(i,j)+
     &                  v0(i,j)*rsteps
cfk#if defined(KNMI)
cfk            ustarsm(i,j)=ustarsm(i,j)+
cfk     &           sqrt(u0(i,j)**2+v0(i,j)**2+vgust**2)*rsteps
cfk#endif
            shfsm(i,j)=shfsm(i,j)+
     &                  shf(i,j,2)*rsteps
            xlhfsm(i,j)=xlhfsm(i,j)+
     &                  fac*evap(i,j,2)*rsteps
            ssrsm(i,j)=ssrsm(i,j)+
     &                  ssr(i,j)*rsteps
            slrdsm(i,j)=slrdsm(i,j)+
     &                  slr(i,j)*rsteps
            slrusm(i,j)=slrusm(i,j)+
     &                  slru(i,j,2)*rsteps
            albsm(i,j)=albsm(i,j)+
     &                  alb1(i,j)*rsteps
            precsm(i,j)=precsm(i,j)+
     &                  (precnv(i,j)+precls(i,j))*rsteps
cfk#if !defined(KNMI)
            snowrsm(i,j)=snowrsm(i,j)+
     &                  0.*rsteps
cfk#else
cfk            snowrsm(i,j)=snowrsm(i,j)+
cfk     &                  (snowcv(i,j)+snowls(i,j))*rsteps
cfk#endif
          enddo
        enddo
     


      return
      end



      subroutine zeroflx

      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL )

      include "com_atm2sea.h"
      include "com_atm2land.h"

C--   1. Set flux integral to zero  

        do j=1,nlat
          do i=1,nlon
            hfints(i,j)=0.
          enddo
        enddo

        do j=1,nlat
          do i=1,nlon
            hfintl(i,j)=0.
          enddo
        enddo

C--   2. set averages to zero

        do j=1,nlat
          do i=1,nlon
            ustrlm(i,j)=0. 
            vstrlm(i,j)=0. 
            u0lm(i,j)=0. 
            v0lm(i,j)=0.
            shflm(i,j)=0.
            xlhflm(i,j)=0.
            ssrlm(i,j)=0.
            slrdlm(i,j)=0. 
            slrulm(i,j)=0. 
            alblm(i,j)=0. 
            preclm(i,j)=0. 
            snowrlm(i,j)=0. 
          enddo
        enddo

        do j=1,nlat
          do i=1,nlon
            ustrsm(i,j)=0. 
            vstrsm(i,j)=0. 
            u0sm(i,j)=0. 
            v0sm(i,j)=0.
cfk#if defined(KNMI)
cfk            ustarsm(i,j)=0.
cfk#endif
            shfsm(i,j)=0.
            xlhfsm(i,j)=0.
            ssrsm(i,j)=0.
            slrdsm(i,j)=0. 
            slrusm(i,j)=0. 
            albsm(i,j)=0. 
            precsm(i,j)=0. 
            snowrsm(i,j)=0. 
          enddo
        enddo


      return
      end    




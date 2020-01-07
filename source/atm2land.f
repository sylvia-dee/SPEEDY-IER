      SUBROUTINE ATM2LAND(HFLXL1,HFINTL,USTRLM,
     &                    VSTRLM,U0LM,V0LM,SHFLM,XLHFLM,
     &                    SSRLM,SLRDLM,SLRULM,ALBLM,
     &                    PRECLM,SNOWRLM,STL01,STL01L)
C-- 
C--   SUBROUTINE ATM2LAND(HFLXL1,HFINTL,USTRLM,
C--  &                    VSTRLM,U0LM,V0LM,SHFLM,XLHFLM,
C--  &                    SSRLM,SLRDLM,SLRULM,ALBLM,
C--  &                    PRECLM,SNOWRLM)
C--   
C--   Purpose :	Update sea-surface tenperature anomalies
C--   Input   : HFLXL1= Daily value of climat. land heat flux
C--           : HFINTL= Net heatflux over land
C--           : USTRLM= Ustress over land 
C--           : VSTRLM= Vstress over land 
C--           : U0LM  = Sfc U wind over land 
C--           : V0LM  = Sfc v wind over land 
C--           : SHFLM = Sensible heatflux over land 
C--           : XLHFLM= Latent heat flux uver land
C--           : SSRLM = Sfc solar radiation land 
C--           : SLRDLM= Sfc downward longw. radiat. land 
C--           : SLRULM= Sfc upward longw. radiat. land
C--           : ALBLM = Albedo over land
C--           : PRECLM= Precip over land 
C--           : SNOWRLM= Snowfall rate over land
C--           : STL01= Daily Land ST
C--           : STL01L= Daily Land ST from previous day
C
C--   Modified common blocks: lsfcanom
C--                           
C--  

      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL )

      include "com_ts_land.h"
      include "com_land.h"

      real hflxl1(ix,il)
      real hfintl(ix,il), ustrlm(ix,il), vstrlm (ix,il), 
     &     u0lm(ix,il), v0lm(ix,il), shflm(ix,il),  
     &     xlhflm(ix,il), ssrlm(ix,il), slrdlm(ix,il), 
     &     slrulm(ix,il), alblm(ix,il), preclm(ix,il), 
     &     snowrlm(ix,il), stl01(ix,il), stl01l(ix,il)
      real hfanom(nlon,nlat)


C--   1. Compute sfc temp. anomalies from flux anomalies

      rnstep = 1./nsteps

      if (ialst.gt.0) then

        do j=1,nlat
          do i=1,nlon
            if (fmaskl1(i,j).gt.0.) then
              fsnow=1.-0.5*snowc1(i,j)
C              hfanom(i,j)=fsnow*(hfintl(i,j)*rnstep-hflxl1(i,j))
              hfanom(i,j)=fsnow*(hfintl(i,j)*rnstep)
              stanoml1(i,j)=stdisl*
     &                      (stanoml1(i,j)+hfanom(i,j)*rhcap2l(i,j))
     &                      + (stl01l(i,j)-stl01(i,j))
            endif
          enddo
        enddo
      endif

C--   3. Set flux integral to zero for next step 

        do j=1,nlat
          do i=1,nlon
            hfintl(i,j)=0.
          enddo
        enddo

      return
      end








 
      SUBROUTINE SEND_TO_COUPLER 
C--
C--   SUBROUTINE SEND_TO_COUPLER
C--
C--   Purpose : coupler sends fields 
C--             sea and atmosphere
C--             updated common blocks: sea_sfcanom, land_sfcanom
C--

      include "atparam.h"
C
      include "com_for_sea.h"
      include "com_for_land.h"
      include "com_atm2sea.h"
      include "com_atm2land.h"

C--  Send fields to  sea-module

      CALL ATM2SEA(OICE1,HFLXS1,HFINTS,USTRSM,
     &             VSTRSM,U0SM,V0SM,SHFSM,XLHFSM,SSRSM,SLRDSM, 
     &             SLRUSM,ALBSM,PRECSM,SNOWRSM,SST01,SST01L)  

C--  Send fields to  land-module

      CALL ATM2LAND(HFLXL1,HFINTL,USTRLM,VSTRLM,
     &              U0LM,V0LM,SHFLM,XLHFLM,SSRLM,SLRDLM, 
     &              SLRULM,ALBLM,PRECLM,SNOWRLM,STL01,STL01L)       


      RETURN
      END

 
      SUBROUTINE GET_FROM_COUPLER(IMONTH,IDAY,IRST)
C--
C--   SUBROUTINE GET_FROM_COUPLER(IMONTH,IDAY)
C--
C--   Purpose : coupler provides fields 
C--             to atmosphere 
C--   Input :   IMONTH, IDAY = specification of date
C--             updated common blocks: SEA_FORDAY, LAND_FORDAY,
C--                                    sea_sfcanom, land_sfcanom
C--             IRST = switch for anomaly definition only in 
C--                    case of writing restart file 

      include "atparam.h"
C
      include "com_for_sea.h"
      include "com_for_land.h"
C   
      integer IRST

C--  Get fields from  sea-module

      CALL SEA2ATM(IMONTH,IDAY,IRST,SST1,SST01,SST01L,OICE1,ALBS,SSTAN,
     &             STANOMI,STANOMS,HFLXS1) 

C--  Get fields from  land-module


      CALL LAND2ATM(IMONTH,IDAY,IRST,STL1,STL01,STL01L,SNOW1,SOILW1,
     &              ALBL,STANOML,HFLXL1)  


      RETURN
      END



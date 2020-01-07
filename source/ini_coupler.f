 
      SUBROUTINE INI_COUPLER(ISTART,ISSTAN)
C--
C--   SUBROUTINE INI_COUPLER(ISTART,ISSTAN)
C--
C--   Purpose : Initialization of coupler
C--   Input   : ISTART  = Parameter for deciding if restart run
C--           : ISSTAN  = Parameter for sst anomalies
C--             initial. common blocks: SEA_FORDAY, LAND_FORDAY,
C--                                     sea_sfcanom, land_sfcanom
C--
      include "atparam.h"
      include "atparam1.h"
C
      include "com_dyncon1.h"
      include "com_dyncon0.h"
      include "com_forcing.h"
      include "com_for_sea.h"
      include "com_for_land.h"

C--  Initialize sea-module

      CALL SEA_INIT(ISTART,ISSTAN,SSTAN,STANOMI,STANOMS,RADANG) 

C--  Initialize land-module

      CALL LAND_INIT(ISTART,STANOML,PHI0,PHIS0,GAMMA,GRAV)

      RETURN
      END


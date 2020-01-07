      SUBROUTINE DELTA2d (ITR,QISO,QTOT)
C--
C--   SUBROUTINE DELTA2d (ITR,QISO,QTOT)
C--
C--   Purpose: Converts isotope mixing ratio to delta value
C--    
C--   Input-    arguments:     QTOT  : reference species
C--                            
C--   Input-output arguments:  QISO  : isotope spechum/delta
C--
C--   David Noone <dcn@colorado.edu> - Tue Feb 22 08:04:06 MST 2011
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      include "com_isocon.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      REAL QISO(NGP), QTOT(NGP)

C-- 1. Convert to delta values
      DO J=1,NGP
        TRRAT = QISO(J)/(QTOT(J)+QTINY)
        QISO(J) = (TRRAT/Rstd(ITR) - 1.)*1000.
      ENDDO

      RETURN
      END


      SUBROUTINE DELTA3d (ITR,QISO,QTOT)
C--
C--   SUBROUTINE DELTA3d (ITR,QISO,QTOT)
C--
C--   Purpose: Converts isotope mixing ratio to delta value
C--    
C--   Input-    arguments:     QTOT  : reference species
C--                            
C--   Input-output arguments:  QISO  : isotope spechum/delta
C--
C--   David Noone <dcn@colorado.edu> - Tue Feb 22 08:04:17 MST 2011
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      include "com_isocon.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      INTEGER ITR
      REAL QISO(NGP,NLEV), QTOT(NGP,NLEV)

C-- 1. Convert to delta values
      DO K=1,NLEV
        DO J=1,NGP
          TRRAT = QISO(J,K)/(QTOT(J,K)+QTINY)
          QISO(J,K) = (TRRAT/Rstd(ITR) - 1.)*1000.
        ENDDO
      ENDDO

      RETURN
      END


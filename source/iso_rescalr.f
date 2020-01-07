      SUBROUTINE RESCALR (QG1,QTEND)
C--
C--   SUBROUTINE RESCALR (QG1,QTEND)
C--
C--   Purpose: Rescales the tracer ratios so that they 
C--       do not drift with respect to the total tracer.
C--       This is essentially a hack. Better not to use it if possible.
C--    
C--   Input-only   arguments:   
C--                            QG1    : specific humidity (gp)
C--   Input-output arguments:  
C--                            QTEND  : spec. hum. tendency (gp)
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      
C     Time parameters
      include "com_tsteps.h"

C     Water isotope tracer variables/constants
      include "com_isocon.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      REAL QG1(NGP,NLEV,NTR), QTEND(NGP,NLEV,NTR)

C-- 1. Initialize for total water
      DO K=1,NLEV
        DO J=1,NGP
          QTEND(J,K,1) = 0.
        ENDDO
      ENDDO

C-- 2. Add a correction to scale tracer mixing ratio with isotope ratio
      DO ITR=2,NTR
        DO K=1,NLEV
          DO J=1,NGP
            if (QG1(j,k,2) .gt. qtiny) then 
              TRRAT = QG1(J,K,ITR)/QG1(J,K,2)
            else
C              TRRAT = 1.0
              TRRAT = 0.0
            endif
            QTEND(J,K,ITR) = (QG1(J,K,1)*TRRAT - QG1(J,K,ITR))/delt2
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END


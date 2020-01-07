      SUBROUTINE CHKQ (QG1,etol,sub)
C--
C--   SUBROUTINE CHKQ (QG1)
C--
C--   Purpose: checks the values of H2O=Q.
C--    
C--   Input-    arguments:     QG1  : specific humidity and tracers
C--                            ETOL : Tolerance for error
C--                            
C--   David Noone <dcn@colorado.edu> - Tue Mar 20 23:04:52 MDT 2012
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      include "com_isocon.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      INTEGER ITR
      REAL QG1(NGP,NLEV,NTR)
      CHARACTER*(*) sub
      REAL ETOL,ERR
      REAL QDEN
      INTEGER NBAD

      Logical lchk

C-- 0. Quick exit, bypasses checking
      if (.not. lqchk) return
     

C-- 1. Do the check
      NBAD = 0

      DO K=1,NLEV
        DO J=1,NGP
          
          qden = max(abs(QG1(J,K,2)),abs(QG1(J,K,1)))
          if (qden .le. qtiny) then 
            err = 0.
          else
            err = abs(QG1(j,K,2)-QG1(J,K,1)) 
          endif

          if (err .gt. etol) then
             write(*,*) 'CHKQ FAILED (',sub,'): ',
     &                    J,K,ERR,QG1(J,K,2),QG1(J,K,1)
             NBAD = NBAD + 1
          endif
        ENDDO
      ENDDO

CC      if (NBAD .eq. 0) write(*,*) 'QCHK (',sub,') All good.'
      if (NBAD .gt. 0) write(*,*) 'QCHK (',sub,') NBAD=',nbad
      if (NBAD .gt. 0) stop'ABORT IN QCHK'

      RETURN
      END


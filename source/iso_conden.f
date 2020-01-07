      SUBROUTINE ISO_CONDEN (NPTS,ITR,EPS,TK,QOLD,DQ,TROLD,DTR)
C--
C--   SUBROUTINE ISO_CONDEN (NPTS,ITR,EPS,TK,QOLD,DQ,TROLD,DTR)
C--
C--   Purpose: Predics changes in isotopic mixing ratio
C--        from total water changes, with isotopic fractionation
C--      (This is just a Rayleigh case, and a place holder
C--       for the Noone 2011 cloud model, which uses eps)
C--    
C--   Input-    arguments:     ITR   : tracer index
C--                            EPS   : precipitation efficiency
C--                            TK    : temperature in kelvin
C--                            QOLD  : initial vapor mass
C--                            DQ    : condensation dendency
C--                            TROLD : tracer vapor mass
C--                            
C--   Input-output arguments:  DTR   : tracer tendency
C--
C--   David Noone <dcn@colorado.edu> - Tue Feb 22 08:03:40 MST 2011
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

C     Time parameters
      include "com_tsteps.h"

C     Water isotope tracer variables
      include "com_isocon.h"

C     External isotope fractionation functions
      REAL ALPLIQ, ALPICE, ALPKSS

      INTEGER ITR,NPTS
      REAL EPS
      REAL TK(NPTS)
      REAL QOLD(NPTS) , DQ(NPTS)
      REAL TROLD(NPTS), DTR(NPTS)

      INTEGER J

      REAL TRRAT
      REAL FRAIN
      REAL ALPEQ, ALPHA
      REAL QNEW,TRNEW

C-- 1. Do a Raylegh fractionation

      DO J=1,NPTS

C       Check for some unusual, but numerically possible
C       bad cases, and do something simpler.


        QNEW  = QOLD(J) + DELT2*DQ(J)

        if (QNEW .gt. QOLD(J) .or.
     &      QNEW .lt. qtiny   .or.
     &      QOLD(J) .lt. qtiny      ) then
       
           if (abs(QOLD(J)) .lt. QTINY) then
              TRRAT = 0.
           else
              TRRAT = TROLD(J)/QOLD(J)
           endif
           DTR(J) = TRRAT*DQ(J)

        else
          if (TK(J) .lt. TFREEZE) then
            alpeq = alpice(ITR,TK(J))
            if (TK(J) .lt. TKINETIC) then
              alpha = alpkss(ITR,TK(J),alpeq)
C              alpha=1.0
            else
              alpha = alpeq
            endif

          else
            alpha = alpliq(ITR,TK(J))
          endif

          FRAIN = QNEW/max(QOLD(J),qtiny)
          TRNEW = TROLD(J)*FRAIN**alpha
          DTR(J) = (TRNEW-TROLD(J))/DELT2

        ENDIF


      ENDDO

      RETURN
      END


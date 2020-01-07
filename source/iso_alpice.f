      FUNCTION ALPICE (ITR,TK)
C--
C--   SUBROUTINE ALPICE (ITR,TK)
C--
C--   Purpose: Evaluates the equilibrium fractionation for
C--       vapour/ice phase change.
C--    
C--   Input- arguments:     ITR   : tracer index
C--                   :     TK    : Temperature in Kelvin
C--                            
C--  David Noone <dcn@colorado.edu> - Sat Nov 12 17:21:22 MST 2011
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      include "com_isocon.h"

C-- Emperical factors for fractionation function

      REAL A_D  , B_D  , C_D  
      parameter(a_D   = 48888., b_D   = -203.10, c_D   = 0.2133)

      REAL A_18O, B_18O, C_18O
      parameter(a_18O = 8312.3, b_18O = -49.192, c_18O = 0.0831)


C
      integer ITR
      real TK
      REAL ALPICE

C-- 0. Check for lousy celcius

      if (TK .lt. 150.) then
          write(*,*) 'ALPICE: TK < 150 (not Kelvin?):',tk
      endif



C-- 1. Do it on a species by species basis

      if (ITR .eq. ixq .or. ITR .eq. ixh2o) then	! Q, or H2O tracer
         alpice = 1.0

      else if (ITR .eq. ixhdo) then			! HDO
CC         alpice = 1.1
         alpice = a_D/(tk*tk) + b_D/tk + c_D
         alpice = exp(alpice)


      else if (ITR .eq. ixh218o) then			! H218O
CC         alpice = 1.01
         alpice = a_18O/(tk*tk) + b_18O/tk + c_18O
         alpice = exp(alpice)

      else
         write(*,*) 'ISO_ALPLIQ: Tracer unknown. ITR=',ITR
         stop'TERMINATED - this is an error'
      endif
    

      RETURN
      END


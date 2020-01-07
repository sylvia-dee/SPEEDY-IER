      FUNCTION ALPLIQ (ITR,TK)
CC      FUNCTION ALPLIQ (ITR)
C--
C--   SUBROUTINE ALPLIQ (ITR,TK)
C--
C--   Purpose: Evaluates the equilibrium fractionation for
C--       vapour/liquid phase change.
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
      parameter(a_D   = 24.844e+3, b_D   = -76.248, c_D   =  52.612e-3)

      REAL A_18O, B_18O, C_18O
      parameter(a_18O = 1.137e+3 , b_18O = -0.4156, c_18O = -2.0667e-3)


C
      integer ITR
      real TK
      REAL ALPLIQ

C-- 0. Check for lousy celcius

      if (TK .lt. 150.) then
          write(*,*) 'ALPLIQ: TK < 150 (not Kelvin?):',tk
      endif


C-- 1. Do it on a species by species basis

      if (ITR .eq. ixq .or. ITR .eq. ixh2o) then	! Q, or H2O tracer
         alpliq = 1.0

      else if (ITR .eq. ixhdo) then			! HDO
C        alpliq = 1.1
        alpliq = a_D/(tk*tk) + b_D/tk + c_D
        alpliq = exp(alpliq)


CC        write(*,*) ITR,tk,a_d,b_d, c_d, alpliq

      else if (ITR .eq. ixh218o) then			! H218O
CC        alpliq = 1.01
        alpliq = a_18O/(tk*tk) + b_18O/tk + c_18O
        alpliq = exp(alpliq)


CC        write(*,*) ITR,tk,a_18o,b_18o, c_18o, alpliq

      else
         write(*,*) 'ISO_ALPLIQ: Tracer unknown. ITR=',ITR
         stop'TERMINATED - this is an error'
      endif
    

      RETURN
      END


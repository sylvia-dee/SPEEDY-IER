
      SUBROUTINE INISO

C--
C--   SUBROUTINE INISO
C--
C--   Purpose : Initialize isotopic constants
C--
C--   initialized common blocks: ISOCON
C--   
C--   Tracers are assumed (set in COM_isocon.h):
C--   ITR=1     Q	: normal water vapor (mass) mixing ratio
C--   ITR=2     H2O	: replicates watercycle, using isotope numerics
C--   ITR=3     HDO	: Variable is (mw_h2o/mw_hdo)*2*Rsmow_DH
C--   ITR=4     H218O	: Variable is (mw_h2o/mw_h218o)*Rsmow_18/16
C--
C--   David Noone <dcn@colorado.edu> - Tue Feb 22 08:04:36 MST 2011
C

      include "atparam.h"
      include "atparam1.h"

      include "com_isocon.h"
 
C--   2. Definition of constants

C     2.1 Standard isotope ratio. this is arbitrary in the model
C      so is best taken as one.
      Rstd(1) = 1.0		
      Rstd(2) = 1.0		
      Rstd(3) = 1.0		
      Rstd(4) = 1.0		

C     2.2 Surface ocean is pretty close to SMOW. Slight enrighment OK.
      Rocn(1) = Rstd(1)
      Rocn(2) = Rstd(2) 
      Rocn(3) = Rstd(3)
      Rocn(4) = Rstd(4)
C      Rocn(3) = Rstd(3)*0.50
C      Rocn(4) = Rstd(4)*0.25

C     2.3 Use a "typical" value for land. 
C      This should be from a soil model.
      Rlnd(1) = Rstd(1)
      Rlnd(2) = Rstd(2) 
      Rlnd(3) = Rstd(3)*(1. - 80./1000.)
      Rlnd(4) = Rstd(4)*(1. - 10./1000.)
C      Rlnd(3) = Rstd(3)*0.5
C      Rlnd(4) = Rstd(4)*0.25

C     2.4 Ratio of diffusivity in air (Merlivat and Jouzel 1979)
      difr(1) = 1.
      difr(2) = 1.
      difr(3) = 0.9755
      difr(4) = 0.9723

C--
      RETURN
      END


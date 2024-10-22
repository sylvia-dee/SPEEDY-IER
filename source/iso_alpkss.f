      FUNCTION ALPKSS (ITR,TK,ALPEQ)
C--
C--   SUBROUTINE ALPKSS (ITR,TK, ALPEQ)
C--
C--   Purpose: Evaluates the kinetic fractionation
C--      for ice deposition under supersaturated conditions.
C--    
C--   Input- arguments:     ITR   : tracer index
C--                   :     TK    : temperature (K)
C--                   :     ALPEQ : equilibrium fractionation
C--                             
C--    David Noone <dcn@colorado.edu> - Sat Nov 12 19:15:06 MST 2011
C--
C Modified Sylvia Dee <sdee@usc.edu> - Wed Apr 4 2012
C Purpose: Calculate kinetic fractionation factor alpha-k as for isotopic fractionation during ice crystal formation in clouds.

C--   Tracers are assumed:
C--   ITR=1     Q	: normal water vapor (mass) mixing ratio
C--   ITR=2     H2O	: replicates watercycle, using isotope numerics
C--   ITR=3     HDO	: Variable is (mw_h2o/mw_hdo)*2*Rsmow_DH
C--   ITR=4     H218O	: Variable is (mw_h2o/mw_h218o)*Rsmow_18/16


C  Define all variables

C     k   = kinetic fractionation factor (isotopic species of water occurring when water is evaporated)
C     D   = molecular diffusivity of water vapor in air, D = 2.44E-5  (Merlivat 1977)
C     Si  = saturation over ice
C     S   = saturation ratio with respect to ice at air temperature
C     as  = alpha s, the equilibrium fractionation coefficient with respect to the solid phase
C     Dr  = (D/D')



C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

C     Water isotope tracer variables
      include "com_isocon.h"
      
C
      integer ITR
      REAL ALPEQ, TK

C  T here will be T2D, the absolute temperature at any point in space.

      real k, Dr, Si

C  1.0: If temperature is above freezing, there will be no ice formation and therefore no kinetic fractionation. Alpha = 1.


      if (TK .ge. 273.16) then
           alphkss = 1.0
           
      else if (TK .lt. 273.16) then

C  1. Establish D/D' value, Dr (D'/D, Merlivat 1977): D'/D = 0.9723 for o18 and 0.9755 for dD.

           if (ITR .EQ. 1) then 
                Dr = 1
      
           else if (ITR .EQ. 2) then 
                Dr = 1
      
           else if (ITR .EQ. 3) then 
                Dr = (1/0.9755)
      
           else if (ITR .EQ. 4) then 
                Dr = (1/0.9723)
           
           endif
      


C  2. Calculate Si (equation based on treatment in IsoCAM, Ciasis and Jouzel 1994)


C Need if statement for temperature here? if it is not less than 0, should not be using this function.
C Equation based on Bergeron-Findesein hypothesis -- uses Celsius so must convert.


           Si = 1.02-(0.0038*(TK-273.16))
           Si = max(Si,1.0)

C           write(*,*) 'Si=', Si
C           write(*,*) 'T=', T
      
C  3. Grab alpha s from alpht function, temperature should be below freezing.
      
      
C           as = alpht(ITR, T) unnecessary because we're using alpheq now instead.
      
C  4. Calculate alpkss for clouds

      
           alpkss = alpeq*(Si/((alpeq*(Dr)*(Si-1))+1))


C           write(*,*) 'alph_ice_cloud = ', alpkss, 'T=', TK
      
      endif
      
      
      
      return
      end
      
      


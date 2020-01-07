C      FUNCTION ALPKOC (ITR,VMAG) (David's Version original text)
C--
C--   SUBROUTINE ALPKOC (ITR,VMAG)
C--
C--   Purpose: Evaluates the kinetic fractionation at the ocean
C--      curface during evaporation. Result is the quantity
C--       1-kmol as per Merlivat and Jouzel (1979)
C--    
C--   Input- arguments:     ITR   : tracer index
C--                   :     VMAG  : 10 meter wind speed
C--                            
C--  David Noone <dcn@colorado.edu> - Sat Nov 12 17:21:22 MST 2011
C--  Modified Sylvia Dee <sdee@usc.edu> - Tues Apr 3 2012


C Sylvia Dee <sdee@usc.edu> - Sun Jul 26 2011


C SUBROUTINE ISO_KINETIC
C Purpose: Calculate kinetic fractionation factor alpha-k as for isotopic fractionation during evaportation, over ocean. 


C--   Tracers are assumed:
C--   ITR=1     Q	: normal water vapor (mass) mixing ratio
C--   ITR=2     H2O	: replicates watercycle, using isotope numerics
C--   ITR=3     HDO	: Variable is (mw_h2o/mw_hdo)*2*Rsmow_DH
C--   ITR=4     H218O	: Variable is (mw_h2o/mw_h218o)*Rsmow_18/16

C     NOTE FROM NIK: alpha= 1-k (alpha should be slightly less than 1).  Then you just multiply the isotope evaporation by alpha (David may have already done that in the code).


C  Define all variables

C     k   = kinetic fractionation factor (isotopic species of water occurring when water is evaporated)
C     D   = molecular diffusivity of water vapor in air, D = 2.44E-5  (Merlivat 1977)
C     ed  = relative difference between molecular diffusivities in air for different isotope species
C     PtPm = (Pt/Pm) describes smooth and rough water-atm interfaces
C     v   = kinematic air viscosity = vmu/density
C     X   = 0.4, the Von Karman Constant
C     z_0 = surface roughness
C     z   = height (m)
C     ubar = sqrt(U0^2+V0^2)
C     us  = u* friction velocity
C     Re  = Surface Roughness Reynolds number = (u*zo)/v
C     n   = 2/3 (smooth surfaces) or 1/2 (rough surfaces)
C     g = gravity, 9.81 m/s^2
C     z = height above sea level = 10 meters

C Variables Brought in from other parts of the code:

C     DENS = density of bottom layer, used to calculate v
C     U0 = near-surface U wind (m/s)
C     V0 = near-surface V wind (m/s)


      real function ALPKOC(ITR,U0,V0)
      implicit none
      
      integer ITR
      real k
      real U0, V0, wind

      wind = sqrt(U0**2 + V0**2)

      if (ITR .eq. 2) then       ! normal water
          alpkoc=1.0
      else if (ITR .eq. 3) then  ! HDO
        if (wind .ge. 7.0) then 
            alpkoc=1-(0.000285*wind+0.00082)*0.88
        else                 
            alpkoc=1-0.006*0.88
        endif
      else if (ITR .eq. 4) then  ! H218O
         if (wind .ge. 7.0) then 
            alpkoc=1-(0.000285*wind+0.00082)
         else                   
            alpkoc=1-0.006
         endif
      
      endif

c
      return
      end

        

CCCCCCCCCCCCCCCCCCCCCC SDEE OLD VERSION-- CCCCCCCCCCCCCCCCCCCCCCCCC
C DXS was too low over oceans. Modified to Match ISOGSM 

C     Function Subroutine iso_kinetic to calculate alphk, alpha for kinetic fractionation


C      real function ALPKOC(ITR,U0,V0,DENS)
C      implicit none
      
C      integer ITR
C      real k, D, ed, PtPm, v, X, z_0, z
C      real U0, V0, ubar, us, Re, n, g, DENS

C 1. Define Constants

C     ed values for O18 and Deuterium are known:

C      if (ITR .EQ. 3) then
C          ed = 0.0251
C      else if (ITR .EQ. 4) then  
C          ed = 0.0285
C      else if (ITR .EQ. 2) then 
C          ed = 0.0
C      endif
      
C     X = Von Karman constant      
C      X = 0.4
      
C      g=9.81
C      z=10
      
      
C 2. Determine us (u* friction velocity), v (kinematic air viscoity), D (molecular diffusivity of water vapor in air)

C C v=vmu/(density of bottom layer).  vmu is the dynamic viscosity of air, which is assumed to be 1.7E-5 
C C D: using values reported in Merlivat, 1977.

      
C      ubar = sqrt(U0**2 + V0**2)
       
C add correction for '0' windspeed. can change later if incorrect. ASK DAVID.

C      if (ubar .EQ. 0) then 
C          us = sqrt(1.5E-3 * 0.1)
C      else
C          us   = sqrt(1.5E-3) * ubar
C      endif

C      write (*,*)'U0, V0, ubar =', U0, V0, us

C      v = 1.7E-5 / DENS
C      D = 2.44E-5 

C      write (*,*) v, D
      
C     D units = (m^2) / s      
      
C 3. Calculate z_o using Charnock's equation

C      z_0 = (us**2)/(81.1*g)
      
    
    
C 4. Determine smooth vs. rough surface suing Re, calculate Re

C      Re = (us**3)/(v*81.1*g)
C      write (*,*) Re
      
C 5. Establish IF statement based on surface roughness and calculate (Pt/Pm)

     
C      if (Re .LT. 1) then
            
C            n=(2/3)
      
C            PtPm = (1/X*log(z/z_0)-5)/(7.3*Re**(1/4)*(v/D)**n)
      
C      else if (Re .GE. 1) then
            
C            n=0.5
      
C            PtPm = (1/X*log((us*z)/(30*v)))/(13.6*(v/D)**n)
            
C      endif
      
      
C 6. Calculate k, alphk

C      k = ((1+ed)**n - 1)/((1+ed)**n+PtPm)

C      if (k .gt. 0) then
C           write(*,*) 'k = ', k
C      endif

C      alpkoc = (1-k)
C      alpkoc =1.0


C      write(*,*) 'alphk = ', alphk
   

C      return 
C      end
     


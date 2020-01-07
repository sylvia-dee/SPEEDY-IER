
C--  /ISOCON/ : Constants for isotope tracers
C--   DIFR   = ratio of diffusivity in air
C--   Rstd   = Standard isotope ratio (arbitrary, dont use SMOW)
C--   Rocn   = Ocean surface isotope ratio (~Rstd)
C--   Rlnd   = Land surface isotope ratio (should be from soil)
C--   ix??   = Tracer index for each species of interest

      REAL DIFR(NTR)
      real RSTD(NTR), ROCN(NTR), RLND(NTR)

      common /ISOCON/ difr, Rstd, Rocn, Rlnd

C--  Set the tracer indices (Use CAM terminolgy)


      parameter(ixq     = 1)
      parameter(ixh2o   = 2)
      parameter(ixhdo   = 3)
      parameter(ixh218o = 4)


C-- TMELT  = temperature above whichh falling condensate is liquid (Kelvin)
      real Tmelt
      parameter (Tmelt = 273.16)

C-- TFREEZE  = temperature below which condensate is frozen (Kelvin)
      real Tfreeze
      parameter (Tfreeze = 261.16)

C-- TKINETIC = temperature below which kinetic effects are applied (Kelvin)
      real Tkinetic
      parameter (Tkinetic = 253.16)

C-- FEQ  = fractional equilibration for falling rain (large-scale and convective)
      real FEQLS, FEQCN
      parameter (FEQLS = 1.0)
      parameter (FEQCN = 1.0)

C-- Parameters for supersaturation function for ice deposition
C   (Note: es/ei - 1 has a slope of -0.0125/K so h = 150% at -40C)

      real ssat0, ssat1
      parameter(ssat0 = 1.0, ssat1=0.004)

C-- Emperical exponent for turbulent/molecular kinetic fractionation
      real enn
      parameter (enn=0.58)

C-- Parameters for effective humidity for falling rain drops. Fraction ambient.(phi-Bony 2009).
      real fheff
      parameter (fheff = 0.925)

C-- Land Model Parameters: Fraction of Transpiration (fracT), Fiddle Factor to fix E>P (fid)
      real fracT
      parameter (fracT = 0.75)
C      real fid
C      parameter (fid = 0.99)
      

C-- Some quantities useful for controlling numerical precision
C     trivially small water mixing ratio
      real qtiny
      parameter (qtiny=1.e-20)

C     trivially small precipitation amount
      real ptiny
      parameter (ptiny = 1.e-9)

C     A number small than but indistringuisable from one.
      real omeps
      parameter (omeps= 1.0-1.e-7)

C     magnitudes to pass checks
      logical lqchk
      real qmag_chk,qtmag_chk
      parameter (lqchk = .false.)
      parameter (qmag_chk  = 1.e-4)	! units g/kg
      parameter (qtmag_chk = 1.e-7)	! units (g/kg)/sec

C-- LISODOUT    = flag to output isotope delta values, rather than mass
      logical lidout
      parameter (lidout = .false.)
C      parameter (lidout = .true.)

C-- LISORESCAL  = flag to rescale tracer H2O to prognostic Q to
C--              avoid numerical drift. Preseve ratios.
      logical lirescl
      parameter (lirescl=.false.)
C      parameter (lirescl=.true.)


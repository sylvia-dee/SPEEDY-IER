C-
C--   /LLSMASK/ land-sea masks (initial. in LAND_INIT)
      common /LLSMASK/ fmask1(ix,il),  fmask0(ix,il),
     .                fmasko1(ix,il), fmaskl1(ix,il)
C--									
C--   /LFORFIX/ Time invariant forcing fields 
C--            (initial. in LAND_INIT)
      common /LFORFIX/ alb0(ix,il)
C--
C--   /LFORMON/ Monthly-mean forcing fields (initial. in LAND_INIT)
      common /LFORMON/ stl12(ix,il,12), snow12(ix,il,12),
     .                 soilw12(ix,il,12)
C--
C--   /LFORDAY/ Daily forcing fields (updated in LAND2ATM)
      common /LFORDAY/ snowc1(ix,il)
C-- 
C--  /land_anom/ Daily anomaly variables 
C--             (updated in ATM2LAND; initial. in LAND_INIT)
      common /land_anom/ stanoml1(ix,il) 
C--
C--   /lheatflx/ Monthly mean and daily forcing fields 
C--     (initial. in LAND_INIT; updated in LAND2ATM)
      common /lsheatflx/ hflxl12(ix,il,12)
C--
C--   /lsfconst/ Constants for sea anomaly model 
C--              (initial. in LAND_INIT)
      common /lsfconst/ rhcapl, stdisl,
     &                  rhcap2l(ix,il)

      common /phi0l/ phi0land(ix,il), phis0land(ix,il)


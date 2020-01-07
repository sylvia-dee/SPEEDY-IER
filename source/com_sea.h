C-
C--   /SLSMASK/ land-sea masks (initial. in ISEA_INIT)
      common /SLSMASK/ fmask1(ix,il),  fmask0(ix,il),
     .                fmasko1(ix,il), fmaskl1(ix,il)
C--									
C--   /SFORFIX/ Time invariant forcing fields 
C--            (initial. in SEA_INIT)
      common /SFORFIX/ alb0(ix,il)
C--
C--   /SFORMON/ Monthly-mean forcing fields (initial. in SEA_INIT)
      common /SFORMON/ sst12(ix,il,12),
     .                oice12(ix,il,12)
C--
C--   /sstanom/ Daily forcing anomaly fields (updated in SEA2ATM)
      common /sstanom/ sstan2(ix,il,2),
     &                 sstan1(ix,il)

C-- 
C--  /sea_anom/ Daily anomaly variables 
C--                (updated in ATM2SEA; initial. in SEA_INIT)
      common /sea_anom/ stanoms1(ix,il), 
     .                  stanomi1(ix,il)

C--
C--   /sheatflx/ Monthly mean and daily forcing fields 
C--     (initial. in SEA_INIT; updated in SEA2ATM)
      common /sheatflx/ hflxs12(ix,il,12)
C--
C--   /ssfconst/ Constants for sea anomaly model
C--              (initial. in SEA_INIT)
      common /ssfconst/ rhcaps, rhcapi, stdiss, stdisi, flxice,
     &                  rhcap2s(ix,il), rhcap2i(ix,il), 
     &                  wobsst(ix,il)


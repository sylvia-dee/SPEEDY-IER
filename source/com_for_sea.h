C--
C--   /SEA_FORDAY/ Daily forcing fields 
C--                (updated in GET_FROM_COUPLER)
      common /SEA_FORDAY/ sst1(ix,il), sst01(ix,il), sst01l(ix,il),
     . oice1(ix,il), albs(ix,il), hflxs1(ix,il) 
C-- 
C--  /sea_sfcanom/ Daily anomaly variables 
C--                (updated in GET_FROM_COUPLER)
      common /sea_sfcanom/ sstan(ix,il), stanoms(ix,il), 
     .                     stanomi(ix,il) 


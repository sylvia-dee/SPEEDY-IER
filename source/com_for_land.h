C--
C--   /LAND_FORDAY/ Daily forcing fields 
C--                 (updated in GET_FROM_COUPLER)
      common /LAND_FORDAY/ stl1(ix,il), stl01(ix,il), stl01l(ix,il),
     . snow1(ix,il),soilw1(ix,il), albl(ix,il), hflxl1(ix,il) 
C-- 
C--  /land_sfcanom/ Daily anomaly variables 
C--             (updated in GET_FROM_COUPLER)
      common /land_sfcanom/ stanoml(ix,il) 



C-- 
C--  /lsfcflux/ Daily integrated surface heatflux(initial. 
C--             in ZEROFLX; updated in ADDFLX, used in ATM2LAND)
      common /lsfcflux/ hfintl(ix,il), ustrlm(ix,il), vstrlm (ix,il), 
     &                  u0lm(ix,il), v0lm(ix,il), shflm(ix,il),  
     &                  xlhflm(ix,il), ssrlm(ix,il), slrdlm(ix,il), 
     &                  slrulm(ix,il), alblm(ix,il), preclm(ix,il), 
     &                  snowrlm(ix,il) 
                       


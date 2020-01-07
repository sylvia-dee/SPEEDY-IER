C-- 
C--  /ssfcflux/ Daily integrated surface heatflux (initial. 
C--             in ZEROFLX; updated in ADDFLX, used in ATM2SEA)      
      common /ssfcflux/ hfints(ix,il), ustrsm(ix,il), vstrsm (ix,il),
     &                  u0sm(ix,il), v0sm(ix,il), ustarsm(ix,il),
     &                  shfsm(ix,il), xlhfsm(ix,il), ssrsm(ix,il),
     &                  slrdsm(ix,il), slrusm(ix,il), albsm(ix,il),
     &                  precsm(ix,il), snowrsm(ix,il)
                      

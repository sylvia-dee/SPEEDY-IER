      FUNCTION ALPKEX (ITR,ALPEQ,heff)
C--
C--   SUBROUTINE ALPKEX (ITR,alpeq,heff)
C--
C--   Purpose: Evaluates the kinetic fractionation
C--      for exchange as rain drops fal
C--    
C--   Input- arguments:     ITR   : tracer index
C--                   :     ALPEQ : equilibrium fractionation
C--                   :     heff  : humidity (< 1)
C--                             
C--    David Noone <dcn@colorado.edu> - Sat Nov 12 19:15:06 MST 2011
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      include "com_isocon.h"
C
      integer ITR
      real heff, check
      REAL ALPEQ

      REAL DDN
        
        
C-- 1. Do it on a species by species basis

      if (ITR .eq. ixq .or. ITR .eq. ixh2o) then	! Q, or H2O tracer
         alpkex = 1.0

      else if (ITR .eq. ixhdo .or.ITR .eq. ixh218o) then  ! H218O or HDO

C 2. Compute gamma, effective fractionation factor as per Stewart (1975) and Bony et al. (2009). Also used in CAM:
    
         ddn = difr(ITR)**enn
C         check = ddn + alpeq*(heff-1.0)
C         if (check .ne. 0.0) then 
C         alpkex = alpeq*ddn*heff/(ddn + alpeq*(heff-1.0))
	 alpkex = (alpeq*ddn*heff)/(ddn + alpeq*(1.0-heff))
C         else 
C           alpkex = alpeq
C         endif
      else
          write(*,*) 'ISO_ALPKEX: Tracer unknown. ITR=',ITR
          stop'TERMINATED - this is an error'
      endif
       

      RETURN
      END



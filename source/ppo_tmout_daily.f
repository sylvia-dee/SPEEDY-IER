      SUBROUTINE TMOUT_DAILY (IMODE)
C--
C--   SUBROUTINE TMOUT_DAILY (IMODE)
C--
C--   Purpose : write daily-means into output files
C--   Input :   IMODE = 0 initialize time-mean arrays to 0
C--             IMODE > 0 write time-means and reset arrays to 0
C--   Modified common blocks : TMSAVE_D,TMSAVE
C--
C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Parameters for post-processing arrays
      include "par_tmean.h"

C     Post-processing arrays (daily means)
      include "com_tmean.h"
      include "com_tmean_daily.h"
      include "com_lflags_dmout.h"

C     Time stepping constants
      include "com_tsteps.h"

C     Physical constants
      include "com_physcon.h"

C     Water isotope tracer variables
      include "com_isocon.h"

      real*4 R4OUT(ngp)

      iitest=0
      if (iitest.eq.1) print *, 'inside TMOUT_DAILY'

      if (imode.eq.0) go to 700


C--   1. Divide the accumulate fields to get the means

      fmean=real(nstppr)/real(nsteps)

      do n=1,ns3d_d
       do j=1,ngp
          SAVE3D_D(j,n)=SAVE3D_D(j,n)*fmean
       enddo
      enddo

      fmean=1./real(nsteps)

      do n=1,ns2d_d-ns2d2
       factor=fmean*FACT2D_D(n)
       do j=1,ngp
         SAVE2D_D(j,n)=SAVE2D_D(j,n)*factor
       enddo
      enddo

      fmean=real(nstppr)/real(nsteps)

      do n=ns2d_d-ns2d2+1,ns2d_d
       factor=fmean*FACT2D_D(n)
       do j=1,ngp
         SAVE2D_D(j,n)=SAVE2D_D(j,n)*factor
       enddo
      enddo

C--   1.1 (Optionally) Convert tracers to isotope delta values
      if (lidout) then
        call delta2d(3,SAVE2D_D(1, 3),SAVE2D_D(1, 2))
        call delta2d(4,SAVE2D_D(1, 4),SAVE2D_D(1, 2))

        call delta2d(3,SAVE2D_D(1, 7),SAVE2D_D(1, 6))
        call delta2d(4,SAVE2D_D(1, 8),SAVE2D_D(1, 6))

        call delta2d(3,SAVE2D_D(1,11),SAVE2D_D(1,10))
        call delta2d(4,SAVE2D_D(1,12),SAVE2D_D(1,10))
C
CC        call delta2d(2,SAVE2D_D(1, 2),SAVE2D_D(1, 1))
CC        call delta2d(2,SAVE2D_D(1, 6),SAVE2D_D(1, 5))
CC        call delta2d(2,SAVE2D_D(1,10),SAVE2D_D(1, 9))
      endif



C--   2. Write daily-mean output file 

      do n=1,ns3d_d
       if(LDOUT(n)) then
         do j=1,ngp
           R4OUT(j)=SAVE3D_D(j,n)
         enddo
         write (17) R4OUT
       endif
      enddo
 
      do n=1,ns2d_d			
       if(LDOUT(ns3d_d+n-1)) then
         do j=1,ngp
           R4OUT(j)=SAVE2D_D(j,n)
         enddo
         write (17) R4OUT
       endif
      enddo


C--   3. Reset daily arrays to zero for the next daily-mean

  700 continue

      if (iitest.eq.1) print*,' reset to zero'

      do n=1,ns3d_d
       do j=1,ngp
	 SAVE3D_D(j,n)=0.
       enddo
      enddo

      do n=1,ns2d_d
       do j=1,ngp
         SAVE2D_D(j,n)=0.
       enddo
      enddo

C--   4. define logical array for daily-mean output

      DO NDOUT=1,NS3D_D+NS2D_D
         LDOUT(NDOUT) = .false.
      ENDDO 

      if (IDOUT.eq.1) then
         LDOUT(4)  = .true.
         LDOUT(7)  = .true.
         LDOUT(14) = .true.
         LDOUT(15) = .true.
      endif

      if (IDOUT.eq.2) then
C         DO NDOUT=1,NS3D_D+NS2D_D-1
         DO NDOUT=1,NS3D_D+NS2D_D
            LDOUT(NDOUT) = .true.
         ENDDO 
      endif

C--   user defined 

      if (IDOUT.eq.3) then 
         LDOUT(14)  = .true.
      endif

C--------------------------------------------------------------
C--
      if (iitest.eq.1) print *, 'end of TMOUT_DAILY'

      RETURN
      END

      SUBROUTINE TMOUT (IMODE)
C--
C--   SUBROUTINE TMOUT (IMODE)
C--
C--   Purpose : write time-means and variances into output files
C--   Input :   IMODE = 0 initialize time-mean arrays to 0
C--             IMODE > 0 write time-means and reset arrays to 0
C--   Modified common blocks : TMSAVE and TMSAVE_D
C--
C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Parameters for post-processing arrays
      include "par_tmean.h"

C     Post-processing arrays (time means)
      include "com_tmean.h"
      include "com_tmean_daily.h"

C     Time stepping constants
      include "com_tsteps.h"

C     Physical constants
      include "com_physcon.h"

C     Forcing fields
      include "com_forcing.h"

C     Isotope standard ratios
      include "com_isocon.h"


C     Fields used to compute omega, psi and chi
      complex VORSP(mx,nx), DIVSP(mx,nx), PSISP(mx,nx), CHISP(mx,nx)
      equivalence (PSISP,CHISP)

      real    DIV3D(ngp,nlev)
      real*4  R4OUT(ngp)
      real*4  R4OUT2(ix,il)

      iitest=0
      if (iitest.eq.1) print *, 'inside TMOUT'

      if (imode.eq.0) go to 700


C--   1. Divide the accumulate fields to get the means

      fmean=real(nstppr)/real(nstout)
      fmean_l=1./real(nstout)

      do n=1,ns3d
       do k=1,kx
        do j=1,ngp
          SAVE3D(j,k,n)=SAVE3D(j,k,n)*fmean
        enddo
       enddo
      enddo

      do n=1,ns2d
       factor=fmean*FACT2D(n)
       do j=1,ngp
         SAVE2D(j,n)=SAVE2D(j,n)*factor
       enddo
      enddo

       do n=1,ns2d_d-ns2d2
        factor=fmean_l*FACT2D_D(n)
        do j=1,ngp
          SAVE2D_L(j,n)=SAVE2D_L(j,n)*factor
        enddo
       enddo



C--   2. Compute omega, psi and chi on p surfaces from wind 

      do k=1,kx
        call VDSPEC (SAVE3D(1,k,3),SAVE3D(1,k,4),
     &               VORSP,DIVSP,2)
        if (ix.eq.iy*4) then
          call TRUNCT (VORSP)
          call TRUNCT (DIVSP)
        endif
        call INVLAP (VORSP,PSISP)
        call GRID (PSISP,SAVE3D(1,k,11),1)
        call INVLAP (DIVSP,CHISP)
        call GRID (CHISP,SAVE3D(1,k,12),1)
        call GRID (DIVSP,DIV3D(1,k),1)
      enddo

      dpr2=0.5*pout(1)*p0
      do j=1,ngp
        SAVE3D(j,1,10)=-DIV3D(j,1)*dpr2
      enddo

      do k=2,kx
        dpr2=0.5*(POUT(k)-POUT(k-1))*p0
        do j=1,ngp
          SAVE3D(j,k,10)=SAVE3D(j,k-1,10)-
     &                  (DIV3D(j,k)+DIV3D(j,k-1))*dpr2
        enddo
      enddo

      do k=1,kx
        do j=1,ngp
          SAVE3D(j,k,11)=SAVE3D(j,k,11)*1.e-6
          SAVE3D(j,k,12)=SAVE3D(j,k,12)*1.e-6
        enddo
      enddo

C--   2.1 (Optionally) convert tracers to isotope delta values
      if (lidout) then
        call delta3d(3,SAVE3D(1,1,8),SAVE3D(1,1,7))
        call delta3d(4,SAVE3D(1,1,9),SAVE3D(1,1,7))
C
CC        call delta3d(2,SAVE3D(1,1,7),SAVE3D(1,1,5))


        call delta2d(3,SAVE2D_L(1, 3),SAVE2D_L(1, 2))
        call delta2d(4,SAVE2D_L(1, 4),SAVE2D_L(1, 2))

        call delta2d(3,SAVE2D_L(1, 7),SAVE2D_L(1, 6))
        call delta2d(4,SAVE2D_L(1, 8),SAVE2D_L(1, 6))

        call delta2d(3,SAVE2D_L(1,11),SAVE2D_L(1,10))
        call delta2d(4,SAVE2D_L(1,12),SAVE2D_L(1,10))
C
CC        call delta2d(2,SAVE2D_L(1, 2),SAVE2D_L(1, 1))
CC        call delta2d(2,SAVE2D_L(1, 6),SAVE2D_L(1, 5))
CC        call delta2d(2,SAVE2D_L(1,10),SAVE2D_L(1, 9))
      endif

C--   3. Write time-mean output file including 3-d and 2-d fields
 
      do n=1,ns3d1
       do k=kx,1,-1
        do j=1,ngp
          R4OUT(j)=SAVE3D(j,k,n)
        enddo
        write (11) R4OUT
       enddo
      enddo

      do n=1,ns2d
       do j=1,ngp
         R4OUT(j)=SAVE2D(j,n)
       enddo
       write (11) R4OUT
      enddo
  
      do n=1,ns2d_d-ns2d2
       do j=1,ngp
         R4OUT(j)=SAVE2D_L(j,n)
       enddo
       write (11) R4OUT
      enddo


C     ----------------------------------------------------------------

      if (ns3d2.gt.0) then

C--   4. Compute variances and covariances

        do n=1,4
         nv=n+ns3d1
         do k=1,kx
          do j=1,ngp
            SAVE3D(j,k,nv)=SAVE3D(j,k,nv)-SAVE3D(j,k,n)**2
          enddo
         enddo
        enddo

        nuv=ns3d1+5
        nvt=ns3d1+6
        do k=1,kx
         do j=1,ngp
           SAVE3D(j,k,nuv)=SAVE3D(j,k,nuv)-SAVE3D(j,k,3)*SAVE3D(j,k,4)
           SAVE3D(j,k,nvt)=SAVE3D(j,k,nvt)-SAVE3D(j,k,2)*SAVE3D(j,k,4)
         enddo
        enddo

C--   5. Write 2-nd order moments 

        do n=ns3d1+1,ns3d1+ns3d2
         do k=kx,1,-1
          do j=1,ngp
           R4OUT(j)=SAVE3D(j,k,n)
          enddo
          write (13) R4OUT
         enddo
        enddo

      endif

C     ----------------------------------------------------------------

      if (ns3d3.gt.0) then

C--   6. Write diabatic forcing fields (in degK/day)

        do n=ns3d1+ns3d2+1,ns3d
         do k=kx,1,-1
          do j=1,ngp
           R4OUT(j)=SAVE3D(j,k,n)*86400.
          enddo
          write (15) R4OUT
         enddo
        enddo

      endif

C     ----------------------------------------------------------------

C--   7. Reset arrays to zero for the next time-mean

  700 continue

      if (iitest.eq.1) print*,' reset to zero'

      do n=1,ns3d
       do k=1,kx
        do j=1,ngp
          SAVE3D(j,k,n)=0.
        enddo
       enddo
      enddo

      do n=1,ns2d
       do j=1,ngp
         SAVE2D(j,n)=0.
       enddo
      enddo
 
      do n=1,ns2d_d-ns2d2
       do j=1,ngp
         SAVE2D_L(j,n)=0.
       enddo
      enddo

C jk test
C      IF (IFLUXCORR .GE. 1) THEN

      IF (IFLUXCORR .GE. 2) THEN

          rewind 51
C          rewind 50 
          do it=1,12
           do j=1,il
             do i=1,ix
              IF(FMASKO1(i,j) .EQ. 1.) THEN
               R4OUT2(i,j)=TMOC(i,j,it)
              ELSE
               R4OUT2(i,j)=9.999E+19
              ENDIF
             enddo
           enddo
           write (51) ((R4OUT2(i,j),i=1,ix),j=il,1,-1)
C           write (50) ((R4OUT2(i,j),i=1,ix),j=il,1,-1)
          enddo

      ENDIF
        
C--
      if (iitest.eq.1) print *, 'end of TMOUT'

      RETURN
      END

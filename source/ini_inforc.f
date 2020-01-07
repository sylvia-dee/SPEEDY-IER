 
      SUBROUTINE INFORC (GRAV)
C--
C--   SUBROUTINE INFORC (GRAV)
C--
C--   Purpose : Read forcing (boundary condition) fields 
C--   Input :   GRAV = gravity accel.
C--             Initialized common blocks: LSMASK, FORFIX

      
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL)

      include "com_tsteps.h" 
      include "com_forcing.h"    
      include "com_anomfor.h"

      real*4 r4inp(ix,il)

      iitest=1

c     set threshold for land-sea mask definition
c     (ie minimum fraction of either land or sea)

      thrsh = 0.1

C--   1. Read time-invariant fields (orography, land-sea mask, sfc albedo)

      if (iitest.ge.1) print*,' read orography' 

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          phi0(i,j) = grav*r4inp(i,j)
        enddo
      enddo

      call truncg (ntrun,phi0,phis0)

cfk#if defined(KNMI)
cfk      do j = 1,il
cfk        do i = 1,ix
cfk        if (fmasko1(i,j) == 1) phis0(i,j)=0.0
cfk        enddo
cfk      enddo
cfk#endif
 
      if (iitest.ge.1) print*,' read fractional land-sea mask'  

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          fmask1(i,j) = r4inp(i,j)
        enddo
      enddo

      if (iitest.ge.1) print*,' read surface albedo'  

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          alb0(i,j) = 0.01*r4inp(i,j)
        enddo
      enddo

      REWIND 20

C--   2. Compute additional land-sea masks
c
c     fmask0  = 1 - fmask1 = sea fraction
c     fmaskl1 = 1 wherever there is any land, 0 elsewhere
c     fmasko1 = 1 wherever there is any sea,  0 elsewhere
c
c     Note that because some points have both land and ocean
c     it is *not* true that: fmasko1+fmaskl1=1 everywhere

      do i=1,ix
       do j=1,il

         fmask0(i,j)=1.0-fmask1(i,j)

         if (fmask1(i,j).ge.thrsh) then
           fmaskl1(i,j)=1.0
         else
           fmaskl1(i,j)=0.0
           fmask1 (i,j)=0.0
           fmask0 (i,j)=1.0
         endif

         if (fmask0(i,j).ge.thrsh) then
           fmasko1(i,j)=1.0
         else
           fmasko1(i,j)=0.0
           fmask0 (i,j)=0.0
           fmask1 (i,j)=1.0
         endif

       enddo
      enddo

C--   3. Initialize Flux correction

C--  3.1 Read SST-climatology

      IF( IFLUXCORR .GT. 0  .OR. IASST .EQ. 6)  THEN 

       do it = 1,12
         read (21) ((r4inp(i,j),i=1,ix),j=il,1,-1)
         do j = 1,il
           do i = 1,ix
             sst12a(i,j,it) = r4inp(i,j)
           enddo
         enddo
       enddo

      CALL FORCHK (fmasko1,sst12a,ix*il,12,0.,400.,273.)

      REWIND 21

      ENDIF
C
C--  3.2 Initialize temperatures

      DO IT=1,12 
         DO J=1,NLAT
           DO I=1,NLON
             TMOC(I,J,IT) = 0.
             TMP(I,J,IT) = 0.
           ENDDO
         ENDDO   
      ENDDO
 
      DO J=1,NLAT
         DO I=1,NLON
           TMOC1(I,J) = 0. 
           SSTM1(I,J) = 0.
           sstan3(I,J,1) = 0.
           sstan3(I,J,2) = 0.
         ENDDO
      ENDDO 

       IF(IFLUXCORR .EQ. 1 .OR. IFLUXCORR .EQ. 3) THEN

         do it = 1,12
           read (50) ((r4inp(i,j),i=1,ix),j=il,1,-1)
           do j = 1,il
             do i = 1,ix
              tmoc(i,j,it) = r4inp(i,j)
             enddo
           enddo
         enddo

       CALL FORCHK (fmasko1,tmoc,ix*il,12,0.,400.,273.)

       ENDIF
  
C--   3.2  Read observed SST anomalies for initial and preceeding month

      IF ( (IOBSNINO .GT. 0 .AND. IFLUXCORR .GT. 0) 
     &     .OR. IASST .EQ. 6) THEN

        do jrec=1,isst0-2
          read (40) dummy4
        enddo

        read (40) ((r4inp(i,j),i=1,ix),j=il,1,-1)

        do j = 1,il
          do i = 1,ix
            sstan3(i,j,1) = r4inp(i,j)
          enddo
        enddo

        if (isst0.gt.1) read (40) ((r4inp(i,j),i=1,ix),j=il,1,-1)

        do j = 1,il
          do i = 1,ix
            sstan3(i,j,2) = r4inp(i,j)
          enddo
        enddo


        CALL FORCHK (fmasko1,sstan3,ix*il,2,-50.,50.,273.)

      ENDIF
C--
      RETURN
      END

cfk      SUBROUTINE TRUNCG (ITR,FG1,FG2)
cfk
cfkC--   SUBROUTINE TRUNCG (ITR,FG1,FG2)
cfkC--   Purpose : compute a spectrally-filtered grid-point field
cfkC--   Input   : ITR : spectral truncation (triangular)
cfkC--           : FG1 : original grid-point field
cfkC--   Output  : FG2 : filtered grid-point field
cfk
cfk      include "atparam.h"
cfk
cfk      REAL FG1 (IX,IL), FG2(IX,IL)
cfk      COMPLEX FSP(MX,NX), ZERO 
cfk
cfk      ZERO = (0.,0.)
cfk
cfk      CALL SPEC (FG1,FSP)
cfk
cfk      DO N=1,NX
cfk        DO M=1,MX
cfk         ITWN=ISC*(M-1)+N-1
cfk       IF (ITWN.GT.ITR) FSP(M,N)=ZERO
cfk      ENDDO
cfk      ENDDO
cfk
cfk      CALL GRID (FSP,FG2,1)
cfk
cfk      RETURN
cfk      END

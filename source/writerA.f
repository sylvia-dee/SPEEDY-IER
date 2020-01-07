



	subroutine writerA (taux,tauy,heatf,saltf,solrf)
        integer imax,jmax
cfk---modification for speedy start
c#include "aparam.h"
        include "atparam.h"
cfk---modification for speedy end

c
c       This is the writer subroutine
c       3 Things need to be done:
c                1) Check if allowed to write , (unit=nf1,istat=1)
c                2) write data (unit=nf2)
c                3) Document that the data has been writen(unit=nf3,istat=1)
c
cfk---modification for speedy start
cfk	parameter (imax=ima,jmax=jma-2)
	parameter (ima=ix,jma=il+2,imap2=ix+2)
cfk      parameter (imax=ima,jmax=jma-2)
        parameter (imax=ima,jmax=jma-2)
cfk---modification for speedy end
        dimension taux(ima,jmax),tauy(ima,jmax),heatf(ima,jmax)
        dimension saltf(ima,jmax),solrf(ima,jmax)
        character*80 fname,fname1,fname2,fname3
        logical first
        integer itime,istat
cfk#ifndef socket
        integer ierrno,stat,statb(12),isize
cfk#endif
        external sleep
cfk#ifndef socket
        external stat
cfk#endif
        data first/.true./
        save first,fname1,fname2,fname3

        if (first) then
           read(200,*)fname
           close(200)
           fname1=trim(fname) // '/fort.104'
           fname2=trim(fname) // '/fort.105'
           fname3=trim(fname) // '/fort.106'
           first=.false.
        endif

        nf1=104 ; nf2=105 ; nf3=106


        itime=5
c
	do
cfk#ifdef socket
cfk          write(6,*)'call sksvr writer'
cfk          call sksrv(iret)
cfk          write(6,*)'exit sksvr writer'
cfk          if(iret.eq.0)then
cfk             istat=1
cfk          else
cfk             istat=0
cfk          endif
cfk#else
clink
          open(nf1,file=fname1)
          read(nf1,*)istat
          close(nf1)

c         isize=8+8
c         isize=12  
c         istat=0
c         do
c            ierrno=stat(trim(fname2),statb)
c            write(6,*)ierrno,statb(8),fname2
c            if((ierrno.eq.0).and.(statb(8).eq.isize))then
c               istat=1
c               exit
c            endif
c            call sleep(2)
c         enddo
cfk#endif
          if (istat.ne.1) then 
             call sleep(itime)
cfk             write(6,*)'not ready to write'
                print*, 'not ready to write'
          else
c
c	Write flux data
c
cfk#ifdef socket
cfk             write(6,*)'skwrt '
cfk             call skwrt(taux,8,imax*jmax,iret)
cfk             call skwrt(tauy,8,imax*jmax,iret)
cfk             call skwrt(heatf,8,imax*jmax,iret)
cfk             call skwrt(saltf,8,imax*jmax,iret)
cfk             call skwrt(solrf,8,imax*jmax,iret)
cfk             write(6,*)'skwrt ',iret
cfk             call skclo()
c
c    write fluxes to disk for restart capability
c
clink
cfk             open(nf2,file=fname2,form='unformatted')
cfk     &            convert='big_endian')
cfk             write(nf2)taux
cfk             write(nf2)tauy
cfk             write(nf2)heatf
cfk             write(nf2)saltf
cfk             write(nf2)solrf
cfk             close(nf2)
c
c 
cfk#else
c
clink
             open(nf2,file=fname2,form='unformatted')
cfk     &            convert='big_endian')
             write(nf2)taux
             write(nf2)tauy
             write(nf2)heatf
             write(nf2)saltf
             write(nf2)solrf
             close(nf2)
c            ierrno=stat(trim(fname2),statb)
c            write(6,*)'finish write',ierrno,statb(8),fname2
cfk             write(6,*)'finish write'
             print*, 'finish write' 
c
c	Document that flux data has been written
c
             open(nf3,file=fname3)
             write(nf3,*)istat
             close(nf3)
c
c	Update nf1 so that data is "not ready to write"
c
             istat=0
             open(nf1,file=fname1)
             write(nf1,*)istat
             close(nf1)
cfk#endif
             exit
          endif
        enddo
        return
	end

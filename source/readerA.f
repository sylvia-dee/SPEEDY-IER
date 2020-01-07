



	subroutine readerA (sstatm,aland)
        integer imax,jmax
cfk---modification for speedy start
c#include "aparam.h"
      include "atparam.h"
cfk---modification for speedy end
c
c       This is the reader subroutine
c       3 Things need to be done:
c           1) Check if the data is ready, (unit=nf1,istat=1)
c           2) Read data (unit=nf2)
c           3) Document that the data has been read(unit=nf3,istat=1)
c
cfk---modification for speedy start
cfk	parameter (imax=ima,jmax=jma-2)
	parameter (ima=ix,jma=il+2,imap2=ix+2)
cfk      parameter (imax=ima,jmax=jma-2)
        parameter (imax=ima,jmax=jma-2)
cfk---modification for speedy end
        dimension sstatm(imax,jmax)
        integer itime,istat,aland(imap2,jma)
cfk#ifndef socket
        integer ierrno,stat,statb(12),isize
cfk#endif
        character*80 fname,fname1,fname2,fname3
        logical first
        external sleep
cfk#ifndef socket
        external stat
cfk#endif
        data first/.true./
        save first,fname1,fname2,fname3

        if (first) then
           read(200,*)fname
	   close(200)
           fname1=trim(fname) // '/fort.101'
           fname2=trim(fname) // '/fort.102'
           fname3=trim(fname) // '/fort.103'
           first=.false.
        endif

        nf1=101 ; nf2=102 ; nf3=103
       

        itime=5
        itime2=15
c
	do
cfk#ifdef socket
cfk          write(6,*)'call sksvr reader'
cfk          call sksrv(iret)
cfk          write(6,*)'exit sksvr reader',iret
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

c         isize=4*imap2*jma+8+8*imax*jmax+8 
c         istat=0
c         do
c           ierrno=stat(trim(fname2),statb)
c           write(6,*)ierrno,isize,statb(8),fname2
c           if((ierrno.eq.0).and.(statb(8).eq.isize))then
c             istat=1
c             exit
c           endif
c           call sleep(2)
c         enddo
          
cfk#endif
c          print*, istat 
          if (istat.ne.1) then 
             call sleep(itime)
cfk             write(6,*)'data not available'
             print*, 'data not available'
          else
c
c	Read Data
c
cfk#ifdef socket
cfk             write(6,*)'call skrd'
cfk             call skrd(sstatm,8,imax*jmax,iret)
cfk             write(6,*)'exit skrd ',iret
cfk             call skrd(aland,4,imap2*jma,iret)
cfk             call skclo()
cfk             write(6,*)'skclo'
cfk#else
clink
             open(nf2,file=fname2,form='unformatted')
cfk     &            convert='big_endian')
             read(nf2)sstatm
             read(nf2)aland
c            rewind(nf2)
c            write(nf2)float(nf2)
             close(nf2)
c            ierrno=stat(trim(fname2),statb)
c            write(6,*)'finish read',ierrno,statb(8),fname2
cfk             write(6,*)'finish read'
             print*, 'finish read'
c
c	Document Data has been read
c
             open(nf3,file=fname3)
             write(nf3,*)istat
             close(nf3)
c
c       Reset nf1 so that next day's data is "not available"
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

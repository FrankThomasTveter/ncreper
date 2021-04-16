      subroutine crefelt(filename,iunit,itype,ltime,itime,
     +                   icode,lspec,ispec,lopt,iopt,ierror)
c
c  NAME:
c     crefelt
c
c  PURPOSE:
c     Create a DNMI felt file.
c     For simple creation of felt files from a program without
c     having to run the standard program NYFELT.
c     If the file exists, it is removed first.
c     The file is closed before exit.
c
c  SYNOPSIS:
c     subroutine crefelt(filename,iunit,itype,ltime,itime,
c    +                   icode,lspec,ispec,lopt,iopt,ierror)
c     integer       iunit,itype,ltime,lopt,icode,lspec,ierror
c     integer       itime(ltime),iopt(lopt)
c     integer*2     ispec(lspec)
c     character*(*) filename
c
c  INPUT:
c     filename  -  the felt file name
c     iunit     -  file unit no. for felt file
c     itype     -  felt file type,
c                      itype=999 : standard felt file
c                      itype=998 : archive  felt file
c                      itype=997 : cyclic archive felt file
c     ltime     -  length of the itime array,
c                      ltime=5 if itype=999
c                      ltime=8 if itype=998,997
c     itime     -  date,time (utc)
c                  for standard felt file (itype=999):
c                      itime(1) : year
c                      itime(2) : month
c                      itime(3) : day
c                      itime(4) : hour
c                      itime(5) : minute
c                  for (cyclic) archive felt file (itype=998,997):
c                      itime(1) : first year
c                      itime(2) : first month
c                      itime(3) : first day
c                      itime(4) : first hour
c                      itime(5) : first minute
c                      itime(6) : no. of time steps
c                      itime(7) : time step type (in itime(8))
c                                   1=years  2=months  3=days
c                                   4=hours  5=minutes
c                      itime(8) : time step
c     icode     -  format code for the ispec array (see ispec)
c     lspec     -  length of the ispec array
c     ispec     -  felt file index specifications,
c                  icode=0:
c                     then mostly 'unspecified' indecies are created,
c                     all with the same producer and grid, lspec=3 and
c                        ispec(1) : producer no. (88=DNMI)
c                        ispec(2) : grid no.
c                        ispec(3) : no. of indecies to be created
c                                   with the above specifications
c                                   and date/time(s) in itime.
c                  icode=1:
c                     consider ispec an array of type ispec(8,nindex),
c                     then lspec=8*nindex, where nindex is the number
c                     of field indecies to be created, then
c                        ispec(1,n) : producer no. (88=DNMI)
c                        ispec(2,n) : grid no.
c                        ispec(3,n) : data type (1=analysis,...)
c                        ispec(4,n) : forecast hour (0 if analysis)
c                        ispec(5,n) : vertical coordinate
c                        ispec(6,n) : parameter
c                        ispec(7,n) : level
c                        ispec(8,n) : second level (usually 0)
c                     if the felt file type is (cyclic) archive,
c                     all indecies will be created for each time step.
c                  other icode's and ispec formats may be added later...
c     lopt      -  length of the iopt array
c                  (possibly required by icode, see icode above)
c     iopt      -  options
c                      iopt(1)=0 : data gaps allowed (default)
c                             =1 : data gaps not allowed
c                  other options may be added later...
c
c  OUTPUT:
c     ierror    -  error exit status
c                       0 = no error
c                       1 = some input specification error
c                           (see printed message),
c                           file not created.
c                       2 = some I/O error (see printed message),
c                           file not created.
c
c  RELATED DOCUMENTATION:
c     Basic document on FELT files:
c          FILE STRUKTUR FOR "SANNTIDS" LAGRING AV GRID-DATA
c               Forskningsavdeling DNMI, oktober 1982
c     Latest update:
c          /usr/local/doc/felt.doc (at DNMI)
c     About the program NYFELT (and felt file features):
c          /usr/local/doc/nyfelt.doc (at DNMI)
c
c  NOTES:
c     - producer numbers are from 1 to 99 (WMO model identifiers),
c       and if more than one used, they must be in increasing order.
c     - the maximum no. of field indecies in one felt file is 32767.
c     - only one producer allowed in (cyclic) archive felt file.
c     - unlike NYFELT, this routine handles minutes.
c     - unlike NYFELT, for cyclic archive felt file the specified
c       date,time will be used as the first in the file.
c
c-----------------------------------------------------------------------
c  DNMI/FoU  26.01.1997  Anstein Foss
c-----------------------------------------------------------------------
c
      implicit none
c
c..input/output
c
      integer       iunit,itype,ltime,lopt,icode,lspec,ierror
      integer       itime(ltime),iopt(lopt)
      integer*2     ispec(lspec)
      character*(*) filename
c
c..local
c
      integer*2 idrec1(1024),idrec2(2,512),indx(16,64)
c
      integer   ierr,nindx,i,j,n,iprod,ios,lrunit,irec,ntime,ni,nt
      integer   itstep,ihstep,mnstep,mn,nh,ny
      integer   itimef(5)
c
c..check input
      ierror=1
      ierr=0
      if(iunit.lt.2) then
        write(6,*) 'CREFELT:'
        write(6,*) 'ERROR in input specification (iunit)'
        return
      end if
      if(itype.lt.997 .or. itype.gt.999) then
        ierr=1
      elseif(itype.eq.999 .and. ltime.ne.5) then
        ierr=1
      elseif(itype.ne.999 .and. ltime.ne.8) then
        ierr=1
      elseif(itype.ne.999) then
c..for standard felt file date,time may be not be a 'real' time,
c..test when (cyclic) archive felt file.
c..note that vtime does not handle minutes and uses itimef(5) as
c..time difference in whole hours (while itimef(5) here is minute).
        mn=itime(5)
        itime(5)=0
        call vtime(itime,ierr)
        itime(5)=mn
        if(mn.lt.0 .or. mn.gt.59) ierr=1
        if(itime(6).lt.1 .or. itime(6).gt.32767) ierr=1
        if(itime(7).lt.1 .or. itime(7).gt.5) ierr=1
        if(itime(8).lt.1) ierr=1
      end if
      if(ierr.ne.0) then
        write(6,*) 'CREFELT:'
        write(6,*) 'ERROR in input specification (itype,ltime,itime)'
        return
      end if
      if(icode.eq.0) then
c..icode=0
        if(lspec.ne.3) then
           write(*,*) 'icode.eq.0.and.lspec .ne. 3'
          ierr=1
        else
          nindx=ispec(3)
          if(nindx.lt.1 .or. nindx.gt.32767) then
             write(*,*) 'icode.eq.0.and.lspec.eq.3.and.'//
     &            '(nindx.lt.1.or.nindx.gt.32767)'
             ierr=1
          end if
          if(ispec(1).lt.1 .or. ispec(1).gt.99) then
             write(*,*) 'icode.eq.0.and.lspec.eq.3.and.'//
     &            '(ispec(1).lt.1.or.ispec(1).gt.99)'
             ierr=1
          end if
        end if
      elseif(icode.eq.1) then
c..icode=1
        nindx=nint(float(lspec)/8.0)
        if(nindx.lt.1 .or. nindx.gt.32767) then
           write(*,*) 'icode.eq.1.and.(nindx.lt.1 .or. nindx.gt.32767)'
          ierr=1
        elseif(lspec.ne.nindx*8) then
           write(*,*) 'icode.eq.1.and.lspec.ne.nindx*8'
          ierr=1
        else
          iprod=1
          do i=1,lspec,8
            if(ispec(i).lt.iprod) then
               write(*,*) ispec
               write(*,*) 'icode.eq.1.and.lspec.ne.nindx*8.and..not.'//
     &              '(nindx.lt.1 .or. nindx.gt.32767).and.'//
     &              'ispec(i).lt.iprod',i,ispec(i),iprod,lspec
               ierr=1
            end if
            iprod=ispec(i)
          end do
        end if
      else
         write(*,*) 'icode.ne.0.and.icode.ne.1'
        ierr=1
      end if
c
      if(ierr.ne.0) then
        write(6,*) 'CREFELT:'
        write(6,*) 'ERROR in input specifications (icode,lspec,ispec)'
        return
      end if
c
      if(lopt.ge.1) then
        if(iopt(1).lt.0 .or. iopt(1).gt.1) ierr=1
      end if
      if(ierr.ne.0) then
        write(6,*) 'CREFELT:'
        write(6,*) 'ERROR in input specifications (lopt,iopt)'
        return
      end if
c
      if(itype.eq.997 .or. itype.eq.998) then
        ntime=itime(6)
        if(ntime*nindx.gt.32767) then
          write(6,*) 'CREFELT:'
          write(6,*) 'ERROR in input specifications (itime,ispec)'
          write(6,*) '      Too many timesteps or field indecies'
          write(6,*) '      The total maximum is 32767'
          return
        end if
      end if
c
      ierror=2
c
c..remove an existing file if it exists (i.e. make a very small file),
c..without doing a system call (with time and memory overhead)
      open(iunit,file=filename,
     +           access='sequential',form='unformatted',
     +           status='old',iostat=ios)
      if(ios.eq.0) then
        i=0
        write(iunit,iostat=ios) i
        close(iunit)
      end if
C     write(*,*) 'A'
c
c..get record length unit (bytes), machine dependant
      call rlunit(lrunit)
c
c..open the felt file
      open(iunit,file=filename,
     +           access='direct',form='unformatted',
     +           recl=2048/lrunit,
     +           status='unknown',iostat=ios,err=900)
C      write(*,*) 'B',lrunit
c
c..initialize (records 1 and 2)
      do i=1,1024
        idrec1(i)=0
      end do
      do i=1,512
        idrec2(1,i)=0
        idrec2(2,i)=0
      end do
c
C      write(*,*) 'C'
      write(iunit,rec=1,iostat=ios,err=900) idrec1
      write(iunit,rec=2,iostat=ios,err=900) idrec2
C      write(*,*) 'D'
c
      irec=2
c
      do n=1,64
        do i=1,16
          indx(i,n)=-32767
        end do
      end do
c
C      write(*,*) 'E'
      ntime=1
      if(itype.eq.997 .or. itype.eq.998) then
        ntime =itime(6)
        itstep=itime(7)
        ihstep=0
        mnstep=0
        if(itstep.eq.3) ihstep=itime(8)*24
        if(itstep.eq.4) ihstep=itime(8)
        if(itstep.eq.5) then
          ihstep=itime(8)/60
          mnstep=itime(8)-ihstep*60
        end if
      end if
c
      do i=1,5
        itimef(i)=itime(i)
      end do
c
C      write(*,*) 'F'
      if(icode.eq.0) then
        iprod=ispec(1)
        nindx=ispec(3)
        idrec2(1,iprod)=irec+1
        idrec2(2,iprod)=1
        idrec1(100+iprod)=nindx*ntime
      elseif(icode.eq.1) then
        nindx=lspec/8
        iprod=0
      end if
c
      ni=0
C      write(*,*) 'G'
c
      do nt=1,ntime
c
        if(icode.eq.0) then
c
          do n=1,nindx
            ni=ni+1
            indx( 1,ni)=ispec(1)
            indx( 2,ni)=ispec(2)
            indx( 3,ni)=itime(1)
            indx( 4,ni)=itimef(2)*100+itimef(3)
            indx( 5,ni)=itimef(4)*100+itimef(5)
            if(ni.eq.64) then
              irec=irec+1
              write(iunit,rec=irec,iostat=ios,err=900) indx
              ni=0
            end if
          end do
c
        elseif(icode.eq.1) then
c
          j=-8
          do n=1,nindx
            j=j+8
            if(iprod.ne.ispec(j+1)) then
              iprod=ispec(j+1)
              idrec2(1,iprod)=irec+1
              idrec2(2,iprod)=ni*16+1
            end if
            idrec1(100+iprod)=idrec1(100+iprod)+1
            ni=ni+1
            indx( 1,ni)=ispec(j+1)
            indx( 2,ni)=ispec(j+2)
            indx( 3,ni)=itimef(1)
            indx( 4,ni)=itimef(2)*100+itimef(3)
            indx( 5,ni)=itimef(4)*100+itimef(5)
            indx( 9,ni)=ispec(j+3)
            indx(10,ni)=ispec(j+4)
            indx(11,ni)=ispec(j+5)
            indx(12,ni)=ispec(j+6)
            indx(13,ni)=ispec(j+7)
            indx(14,ni)=ispec(j+8)
            if(ni.eq.64) then
              irec=irec+1
              write(iunit,rec=irec,iostat=ios,err=900) indx
              ni=0
            end if
          end do
c
        end if
c
        if(nt.lt.ntime) then
c
          if(itstep.eq.1) then
            itimef(1)=itimef(1)+itime(8)
          elseif(itstep.eq.2) then
            itimef(2)=itimef(2)+itime(8)
            ny=(itimef(2)-1)/12
            itimef(1)=itimef(1)+ny
            itimef(2)=itimef(2)-ny*12
          else
c..note that vtime does not handle minutes and uses itimef(5) as
c..time difference in whole hours (while itimef(5) here is minute).
            mn=itimef(5)+mnstep
            nh=mn/60
            mn=mn-nh*60
            itimef(5)=ihstep+nh
            if(itimef(5).ne.0) call vtime(itimef,ierr)
            itimef(5)=mn
          end if
c
        end if
c
      end do
c
C      write(*,*) 'I'
      if(ni.gt.0) then
c..complete last index record
        do n=ni+1,64
          do i=1,16
            indx(i,n)=-1
          end do
        end do
        irec=irec+1
        write(iunit,rec=irec,iostat=ios,err=900) indx
      end if
c
      idrec1( 1)=itype
      idrec1( 5)=itimef(1)
      idrec1( 6)=itimef(2)*100+itimef(3)
      idrec1( 7)=itimef(4)*100+itimef(5)
      idrec1( 8)=irec
      idrec1( 9)=irec
      idrec1(10)=0
      idrec1(11)=nindx*ntime
      if(lopt.ge.1) then
c..data gaps  (0=on 1=off)
        if(iopt(1).eq.1) idrec1(14)=1
      end if
      if(itype.eq.998 .or. itype.eq.997) then
c..(cyclic) archive: first and last date/time.
        idrec1(20)=itime(1)
        idrec1(21)=itime(2)*100+itime(3)
        idrec1(22)=itime(4)*100+itime(5)
        idrec1(23)=itimef(1)
        idrec1(24)=itimef(2)*100+itimef(3)
        idrec1(25)=itimef(4)*100+itimef(5)
        idrec1(26)=ntime
        idrec1(27)=nindx
        idrec1(28)=iprod
        idrec1(29)=itime(7)
        idrec1(30)=itime(8)
      end if
c
C      write(*,*) 'J'
c..update record 2 and 1 (in this order)
      write(iunit,rec=2,iostat=ios,err=900) idrec2
      write(iunit,rec=1,iostat=ios,err=900) idrec1
c
      ierror=0
c
  900 continue
C      write(*,*) 'K',ios
      close(iunit)
c
      return
      end

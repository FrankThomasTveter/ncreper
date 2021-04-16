      subroutine mrfelt(mode,filnam,iunit,in,ipack,lfield,field,fscale,
     *                  ldata,idata,ierror)
c
c  NAME:
c     mrfelt
c
c  PURPOSE:
c     Master read field routine.
c     Open file, read field(s) from Felt file, scale data and close file.
c     Using the standard routine rfelt to read fields.
c
c  SYNOPSIS:
c     subroutine mrfelt(mode,filnam,iunit,in,ipack,lfield,field,fscale,
c    +                  ldata,idata,ierror)
c     integer       mode,iunit,ipack,lfield,ldata,ierror
c     real          fscale
c     real          field(lfield)
c     integer*2     in(16),idata(ldata)
c     character*(*) filnam
c
c  INPUT:
c     mode    -  0 = combination of 1, 2 and 3 below
c                1 = open Felt file
c                2 = read field
c                3 = close Felt file
c		 10,11,12,13 are the same as 0,1,2,3 except that no
c		             error messages are printed
c     filnam  -  file name (if filename(1:1)='*' the file is opened
c			    without 'file=filnam', only with unit no.)
c     iunit   -  file unit no.
c     in      -  parameters describing the field ('innholdsfortegnelse'):
c                  in( 1) = producer no.
c                  in( 2) = grid no.
c                  in( 3) = -32767  or  year
c                  in( 4) = -32767  or  month*100+day
c                  in( 5) = -32767  or  hour*100+minute
c                  in( 9) = data type
c                  in(10) = forecast length (hours)
c                  in(11) = vertical coordinate
c                  in(12) = parameter no.
c                  in(13) = level 1
c                  in(14) = level 2 (usually 0)
c                if at least one of the in(3:5) equals -32767 then
c                date/time is not checked, always set date/time in
c                in(3:5) if the file is an Archive or Cyclic Archive
c                Felt file.
c     ipack   -  0 = return the field in the idata array (no unpacking)
c                1 = return the field in the field array and the
c                    field identification in the first part of
c                    the idata array, the field will be scaled using
c                    scaling exponent in idata(20) and fscale,
c                    not checking undefined/missing values
c                2 = as 1 above, but will now also check
c                    undefined/missing values, i.e. input elements
c                    with the value -32767 will be set to +1.e+35 .
c     lfield  -  maximum length of the field (nx*ny), for ipack=1,2
c     fscale  -  extra scaling of field, for ipack=1,2 (usually 1.)
c     ldata   -  maximum length of identification and field (20+nx*ny),
c		 for mode=1,11 and ldata>1:
c		     returning the first ldata words (max. 1024) in the
c		     idata array from record 1 in the file, and if
c		     mode=1,11 and ldata>1024 also returning ldata-1024
c		     words (max. 1024) from record 2
c
c  OUTPUT:
c     field   -  the field, if ipack=1,2
c     idata   -  identification and field (all ipack values),
c                grid dimensions are found in the identification
c                (i.e. the 20 first elements of idata):
c                       nx=idata(10)
c                       ny=idata(11)
c		 for mode=1,11: see ldata above
c     ierror  -  exit status
c                  0 = no error
c                  1 = some error (see printed message)
c
c  NOTES:
c     +  Note the use of ldata and idata for mode=1,11 if ldata>1,
c	 usually it is enough to return 32 words if used at all.
c     +  If the Felt file are opened elsewhere, then a call with
c        mode=1 is not necessary.
c     +  If the Felt file are closed elsewhere, then a call with
c        mode=3 is not necessary.
c     +  As rfelt really is a simple routine to use, you can keep
c        more than one file open when using mrfelt.
c     +  Use mode=0 only when reading one or a few fields,
c        otherwise this will cause unnecessary overhead.
c     +  Call sequence:
c           call mrfelt(1,filnam,iunit,in,ipack,1,1,1,1,1,ierror)
c           do while <more fields to read>
c             call mrfelt(2,filnam,iunit,in,ipack,lfield,field,fscale,
c        -                  ldata,idata,ierror)
c           end do
c           call mrfelt(3,filnam,iunit,in,ipack,1,1,1,1,1,ierror)
c     +  When the field array is not used (ipack=0):
c           call mrfelt(mode,filnam,iunit,in,0,1,1,1,ldata,idata,ierror)
c
c-----------------------------------------------------------------------
c  DNMI/FoU  01.10.1993  Anstein Foss
c  DNMI/FoU  02.12.1993  Anstein Foss
c  DNMI/FoU  15.03.1995  Anstein Foss ... silent modes
c  DNMI/FoU  14.02.1996  Anstein Foss ... ldata,idata for mode=1,11
c  DNMI/FoU  10.11.2000  Anstein Foss ... automatic byte swap
c  DNMI/FoU  03.03.2003  Anstein Foss ... rfelt ip=-5 error
c-----------------------------------------------------------------------
c
c..input/output
      integer       mode,iunit,ipack,lfield,ldata,ierror
      real          field(lfield),fscale
      integer*2     in(16),idata(ldata)
      character*(*) filnam
c
c..local
      integer   ierr(3),ihelp(6),lrunit,lunit
      real      undef
      logical   swapfile,swap
c
      data lrunit/0/
      data ihelp/6*0/
      data lunit/0/
      data undef/+1.e+35/
c
      ierror=0
c
      imode=mode
      iprint=1
      if(mode.ge.10 .and. mode.le.13) then
	imode=mode-10
	iprint=0
      end if
c
      if(imode.lt.0 .or. imode.gt.3) then
        write(6,*) '**mrfelt** unknown mode: ',mode
        ierror=1
        return
      end if
c
      if(imode.ne.0 .and. imode.ne.1) goto 200
c
      if (lrunit.lt.1) call rlunit(lrunit)
c
      ihelp(1)=0
c
      if(filnam(1:1).ne.'*') then
        open(iunit,file=filnam,
     *             access='direct',form='unformatted',
     *             recl=2048/lrunit,
     *             status='old',iostat=ios)
      else
        open(iunit,access='direct',form='unformatted',
     *             recl=2048/lrunit,
     *             status='old',iostat=ios)
      end if
c
      if (ios.eq.0) swap= swapfile(-iunit)
c
      if(imode.eq.1 .and. ldata.gt.1 .and. ios.eq.0) then
c..return part of record 1 and possibly 2
	lread=min(ldata,1024)
	read(iunit,rec=1,iostat=ios) (idata(i),i=1,lread)
	if(ldata.gt.1024 .and. ios.eq.0) then
	  lread=min(ldata,2048)
	  read(iunit,rec=2,iostat=ios) (idata(i),i=1025,lread)
	end if
	if (ios.eq.0 .and. swap) call bswap2(lread,idata)
      end if
c
      if(ios.ne.0) then
	if(iprint.eq.1) then
          write(6,*) '**mrfelt** open/read error.  iostat= ',ios
          call pioerr(6,ios)
          write(6,*) '** ',filnam
          write(6,*) '** file unit: ',iunit
	end if
	lunit=0
        ierror=1
        return
      end if
c
      if(imode.eq.1) return
c
  200 if(imode.ne.0 .and. imode.ne.2) goto 300
c
      if(iunit.ne.lunit) then
        lunit=iunit
        ihelp(1)=0
      end if
c
      ierr(1)=0
      ierr(2)=0
      ierr(3)=0
      irtime=0
      if(in(3).ne.-32767 .and. in(4).ne.-32767
     +                   .and. in(5).ne.-32767) then
        ierr(2)=999
        irtime=1
      end if
c
      call rfelt(iunit,ip,in,idata,ldata,ierr,ihelp)
c
      if(ip.ne.1 .and. iprint.eq.1) then
c..field not read from file
        write(6,*) '**mrfelt**  field not found'
        write(6,*) '** ',filnam
        write(6,*) '** file unit: ',iunit
        write(6,1020) ip,ierr,in(1),in(2),(in(i),i=9,15)
 1020   format(' **           ip=',i4,'  ierr=',3i9,/,
     *         ' ** in(1,2,9-15): ',9i6)
        if(irtime.eq.1) write(6,1021) (in(i),i=3,5)
 1021   format(' ** in(3-5):      ',3i6)
        if(ip.eq.0) then
          write(6,*) '** no field data'
        elseif(ip.eq.-1) then
          write(6,*) '** i/o error (see ierr)'
          do i=1,3
            if(ierr(i).ne.0) call pioerr(6,ierr(i))
          end do
        elseif(ip.eq.-2) then
          write(6,*) '** field not defined in felt file index'
        elseif(ip.eq.-3) then
          write(6,*) '** idata(ldata) too small.  ldata= ',ldata
        elseif(ip.eq.-4) then
          write(6,*) '** producer not defined: ',in(1)
        elseif(ip.eq.-5) then
          write(6,*) '** corrupt pointer from index to data'
        end if
      end if
      if(ip.ne.1) then
        ierror=1
        return
      end if
c
      if(ipack.eq.1 .or. ipack.eq.2) then
        nx=idata(10)
        ny=idata(11)
        ldat=nx*ny
        if(ldat.gt.lfield) then
	  if(iprint.eq.1) then
            write(6,*) '**mrfelt** field(lfield) too small'
            write(6,*) '**   lfield,nx,ny: ',lfield,nx,ny
	  end if
          ierror=1
          return
        end if
        iscale=idata(20)
        scale=10.**iscale
        if(fscale.ne.0 .and. fscale.ne.1.) scale=scale*fscale
        if(ipack.eq.1) then
          do i=1,ldat
            field(i)=scale*float(idata(20+i))
          end do
        else
          do i=1,ldat
            if(idata(20+i).ne.-32767) then
              field(i)=scale*float(idata(20+i))
            else
              field(i)=undef
            end if
          end do
        end if
      end if
c
  300 if(imode.ne.0 .and. imode.ne.3) return
c
      close(iunit)
c
      return
      end

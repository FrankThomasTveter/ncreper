      subroutine mwfelt(mode,filnam,iunit,ipack,lfield,field,fscale,
     *                  ldata,idata,ierror)
c
c  NAME:
c     mwfelt
c
c  PURPOSE:
c     Master write field routine.
c     Open file, scale data, write field(s) to Felt file and close file.
c     Using the standard routine wfelt to write fields.
c
c  SYNOPSIS:
c     subroutine mwfelt(mode,filnam,iunit,ipack,lfield,field,fscale,
c    +                  ldata,idata,ierror)
c     integer       mode,iunit,ipack,lfield,ldata,ierror
c     real          fscale
c     real          field(lfield)
c     integer*2     idata(ldata)
c     character*(*) filnam
c
c  INPUT:
c     mode    -  0 = combination of 1, 2 and 3 below
c                1 = open Felt file and initialize
c                2 = write field
c                3 = update record 1 and close Felt file
c		 10,11,12,13 are the same as 0,1,2,3 except that no
c		             error messages are printed
c     filnam  -  file name (if filename(1:1)='*' the file is opened
c			    without 'file=filnam', only with unit no.)
c     iunit   -  file unit no.
c     ipack   -  0 = the field is in the idata array (no packing)
c                1 = the field is in the field array and the
c                    field identification is in the first part of
c                    the idata array, the field will be scaled using
c                    scaling exponent in idata(20) and fscale,
c		     if idata(20)=-32767 then the best possible
c		     scaling will be found and used,
c                    not checking undefined/missing values
c                2 = as 1 above, but will now also check
c                    undefined/missing values, i.e. elements with
c                    the value +1.e+35 (-32767 in output field)
c     lfield  -  length of the field (nx*ny), for ipack=1,2
c     field   -  field to be packed, for ipack=1,2
c     fscale  -  extra scaling of field, for ipack=1,2 (usually 1.)
c     ldata   -  maximum length of identification and field (20+nx*ny),
c		 in some cases there also is extra geometry
c		 identification after the field data,
c		 for mode=1,11,3,13 and ldata>1:
c		     returning the first ldata words (max. 1024) in the
c		     idata array from record 1 in the file, and if
c		     mode=1,11 and ldata>1024 also returning ldata-1024
c		     words (max. 1024) from record 2 (for mode=1,11 this
c		     is contents before any update and for mode=3,13
c		     after last update of the file header)
c     idata   -  identification and field if ipack=0,
c                identification and space for output field
c                if ipack=1,2
c
c  OUTPUT:
c     idata   -  for mode=1,11,3,13: see ldata above
c     ierror  -  exit status
c                  0 = no error
c                  1 = some error (see printed message)
c
c  NOTES:
c     +  Note the use of ldata and idata for mode=1,11,3,13 if ldata>1,
c	 usually it is enough to return 32 words if used at all.
c     +  If field exists on a Felt file, it is overwritten.
c     +  One cannot write to a second file without closing the
c        first file.
c     +  Use mode=0 only when writing one or a few fields,
c        otherwise this will cause unnecessary overhead.
c     +  Call sequence:
c           call mwfelt(1,filnam,iunit,ipack,1,1,1,1,1,ierror)
c           do while <more fields to write>
c             call mwfelt(2,filnam,iunit,ipack,lfield,field,fscale,
c        -                  ldata,idata,ierror)
c           end do
c           call mwfelt(3,filnam,iunit,ipack,1,1,1,1,1,ierror)
c     +  When the field array is not used (ipack=0):
c           call mwfelt(mode,filnam,iunit,0,1,1,1,ldata,idata,ierror)
c
c-----------------------------------------------------------------------
c  DNMI/FoU  01.10.1993  Anstein Foss
c  DNMI/FoU  12.03.1994  Anstein Foss
c  DNMI/FoU  15.03.1995  Anstein Foss ... silent modes
c  DNMI/FoU  31.03.1995  Anstein Foss ... automatic scaling
c  DNMI/FoU  09.06.1995  Anstein Foss ... extra ident, gridtype 3,4
c  DNMI/FoU  14.02.1996  Anstein Foss ... ldata,idata for mode=1,11,3,13
c  DNMI/FoU  22.03.2003  Anstein Foss ... automatic byte swap
c-----------------------------------------------------------------------
c
c..record length in integer*2 words
      integer lrec
      parameter (lrec=1024)
c
c..input/output
      integer       mode,iunit,ipack,lfield,ldata,ierror
      real          field(lfield),fscale
      integer*2     idata(ldata)
      character*(*) filnam
c
c..local
      integer*2 idrec1(lrec),idrec2(lrec),idsav1(32)
      integer*2 in(16)
      integer   ierr(3),ihelp(6),lrunit,lunit,nwrite
      real      undef
      character*256 filold
      logical   swapfile,swap
c
      data lrunit/0/
      data ihelp/6*0/
      data lunit/0/
      data nwrite/0/
      data undef/+1.e+35/
      data filold/'***'/
      data swap/.false./
c
      save idrec1, idrec2, idsav1
c
      if (lrunit.lt.1) call rlunit(lrunit)
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
        write(6,*) '**mwfelt** unknown mode: ',mode
        ierror=1
        return
      end if
c
      if(imode.ne.0 .and. imode.ne.1) goto 200
c
      filold=filnam
      lunit=iunit
      nwrite=0
      ihelp(1)=0
      swap=.false.
c
c..check update flag (and possibly correct file header after crash)
      call wcfelt(filnam,iunit,iretur)
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
      if(ios.eq.0) swap= swapfile(-iunit)
c
      if(ios.eq.0) read(iunit,rec=1,iostat=ios) idrec1
      if(ios.eq.0) read(iunit,rec=2,iostat=ios) idrec2
c
      if(ios.eq.0 .and. swap) then
        call bswap2(lrec,idrec1)
        call bswap2(lrec,idrec2)
      end if
c
      if(imode.eq.1 .and. ldata.gt.1 .and. ios.eq.0) then
c..return part of record 1 and possibly 2
	lread=min(ldata,lrec)
	do i=1,lread
	  idata(i)=idrec1(i)
	end do
	if(ldata.gt.lrec) then
	  lread=min(ldata,lrec*2)
	  do i=lrec+1,lread
	    idata(i)=idrec2(i-lrec)
	  end do
	end if
      end if
c
      if(ios.ne.0) then
	if(iprint.eq.1) then
          write(6,*) '**mwfelt** open/read error.  iostat= ',ios
          call pioerr(6,ios)
          write(6,*) '** ',filnam
          write(6,*) '** file unit: ',iunit
	end if
	lunit=0
        filold='****'
        ierror=1
        return
      end if
c
      do i=1,32
        idsav1(i)=idrec1(i)
      end do
c
      if(imode.eq.1) return
c
  200 if(imode.ne.0 .and. imode.ne.2) goto 300
c
      if(filnam.ne.filold .or. iunit.ne.lunit) then
	if(iprint.eq.1) then
        write(6,*) '**mwfelt** you must initialize mwfelt for the file:'
        write(6,*) '** ',filnam
        write(6,*) '** file unit: ',iunit
	end if
        ierror=1
        return
      end if
c
      nx=idata(10)
      ny=idata(11)
      if((ipack.eq.1 .or. ipack.eq.2) .and. lfield.ne.nx*ny) then
	if(iprint.eq.1) then
          write(6,*) '**mwfelt** field size error',
     +                               ' (lfield =/= nx*ny)'
          write(6,*) '**      lfield,nx*ny: ',lfield,nx*ny
          write(6,*) '**             nx,ny: ',nx,ny
	end if
	ierror=1
	return
      end if
      igtype=idata(9)
      lxid=0
      if(igtype.gt.999) then
	i=igtype
	igtype=igtype/1000
	lxid=i-igtype*1000
      end if
      ldat=20+nx*ny+lxid
      if(ldata.lt.ldat) then
	if(iprint.eq.1 .and. lxid.eq.0) then
          write(6,*) '**mwfelt** idata(ldata) too small',
     +                           ' (ldata=20+lfield)'
          write(6,*) '**   ldata,lfield: ',ldata,lfield
	elseif(iprint.eq.1) then
          write(6,*) '**mwfelt** idata(ldata) too small',
     +                           ' (ldata=20+lfield+lxid)'
          write(6,*) '**   ldata,lfield,lxid: ',ldata,lfield,lxid
          write(6,*) '**   (lxid words geometry specification',
     +			    ' after field data)'
	end if
	ierror=1
	return
      end if
c
      if(ipack.eq.1 .or. ipack.eq.2) then
	if(idata(20).eq.-32767) then
c..automatic scaling
	  fmax=0.
          if(ipack.eq.1) then
            do i=1,lfield
              fmax=max(fmax,abs(field(i)))
            end do
          else
            udef=0.9*undef
            do i=1,lfield
              if(field(i).lt.udef) fmax=max(fmax,abs(field(i)))
            end do
          end if
	  if(fmax.gt.0.) then
            if(fscale.ne.0 .and. fscale.ne.1.) fmax=fmax*abs(fscale)
	    iscale=log10(fmax)-4.
	    ifmax=nint(fmax*10.**(-iscale))
	    if(ifmax.lt.3278) then
	      iscale=iscale-1
	      ifmax=nint(fmax*10.**(-iscale))
	    end if
	    if(ifmax.gt.32766) iscale=iscale+1
	    iscale=max(iscale,-30)
	    idata(20)=iscale
	  else
	    idata(20)=0
	  end if
	end if
        iscale=idata(20)
        scale=10.**(-iscale)
        if(fscale.ne.0 .and. fscale.ne.1.) scale=scale*fscale
        if(ipack.eq.1) then
          do i=1,lfield
            idata(20+i)=nint(scale*field(i))
          end do
        else
          udef=0.9*undef
          do i=1,lfield
            if(field(i).lt.udef) then
              idata(20+i)=nint(scale*field(i))
            else
              idata(20+i)=-32767
            end if
          end do
        end if
      end if
c
      in( 1)=idata(1)
      in( 2)=idata(2)
      in( 9)=idata(3)
      in(10)=idata(4)
      in(11)=idata(5)
      in(12)=idata(6)
      in(13)=idata(7)
      in(14)=idata(8)
c.old in(15)=idata(9)
c
      ierr(1)=1
      ierr(2)=0
      ierr(3)=0
c
      iwtime=0
      if(idrec1(1).eq.998 .or. idrec1(1).eq.997) then
c..Archive FELT or Cyclic_Archive FELT file:
c..use date/time
        in( 3)=idata(12)
        in( 4)=idata(13)
        in( 5)=idata(14)
        ierr(2)=999
        if(idrec1(1).eq.997) ierr(2)=997
        iwtime=1
      end if
c
      call wfelt(iunit,ip,in,idrec1,idrec2,ldat,idata,ierr,ihelp)
c
      if(ip.eq.1) then
c
        nwrite=nwrite+1
        if(nwrite.eq.1 .and. imode.ne.0) then
c..waiting as long as possible to touch the file,
c..set update flag
          idrec1(15)=1
          if(swap) call bswap2(lrec,idrec1)
          write(iunit,rec=1,iostat=ios) idrec1
          if(swap) call bswap2(lrec,idrec1)
          if(ios.ne.0) then
	    if(iprint.eq.1) then
            write(6,*) '**mwfelt** write error (record 1).  iostat= ',
     +                             ios
            write(6,*) '** ',filnam
            write(6,*) '** file unit: ',iunit
	    end if
            ierror=1
            return
          end if
        end if
c
      elseif(iprint.eq.1) then
c
c..field not written to file
        write(6,*) '**mwfelt**  field not written'
        write(6,*) '** ',filnam
        write(6,*) '** file unit: ',iunit
        write(6,1020) ip,ierr,in(1),in(2),(in(i),i=9,15)
 1020   format(' **           ip=',i4,'  ierr=',3i9,/,
     *         ' ** in(1,2,9-15): ',9i6)
        if(iwtime.eq.1) write(6,1021) (in(i),i=3,5)
 1021   format(' ** in(3-5):      ',3i6)
        if(ip.eq.0) then
          write(6,*) '** field not defined in felt file index'
        elseif(ip.eq.-1) then
          write(6,*) '** i/o error (see ierr)'
          do i=1,3
            if(ierr(i).ne.0) call pioerr(6,ierr(i))
          end do
        elseif(ip.eq.-2) then
          write(6,*) '** existing field not overwritten'
        elseif(ip.eq.-3) then
          write(6,*) '** existing field with different size'
        elseif(ip.eq.-4) then
          write(6,*) '** producer not defined: ',in(1)
        elseif(ip.eq.-5) then
          write(6,*) '** too small output ident+field: ',ldat
        elseif(ip.eq.-6) then
          write(6,*) '** too big ident+field: ',ldat
        elseif(ip.eq.-7) then
          write(6,*) '** maximum file size exceeded'
        end if
        ierror=1
        return
c
      else
c
        ierror=1
        return
c
      endif
c
  300 if(imode.ne.0 .and. imode.ne.3) return
c
      if(filnam.ne.filold .or. iunit.ne.lunit) then
	if(iprint.eq.1) then
          write(6,*) '**mwfelt** not initialized for the file:'
          write(6,*) '** ',filnam
          write(6,*) '** file unit: ',iunit
	end if
        ierror=1
        return
      end if
c
      if(nwrite.gt.0) then
c
	if(iprint.eq.1) then
          write(6,*) '**mwfelt** record 1 at start and after update:'
          write(6,*) '** ',filnam
          write(6,*) '** file unit: ',iunit
          write(6,fmt='('' **'',12i6)') (idsav1(i),i=1,12)
	end if
c
        call daytim(iyear,month,iday,ihour,minute,isecnd)
        idrec1( 2)=iyear
        idrec1( 3)=month*100+iday
        idrec1( 4)=ihour*100+minute
c
c..reset update_flag
        idrec1(15)=0
c
	if(swap) call bswap2(lrec,idrec1)
c
        write(iunit,rec=1,iostat=ios) idrec1
c
	if(swap) call bswap2(lrec,idrec1)
c
        if(ios.ne.0) then
	  if(iprint.eq.1) then
          write(6,*) '**mwfelt** write error (record 1). iostat=',ios
          call pioerr(6,ios)
          write(6,*) '** ',filnam
          write(6,*) '** file unit: ',iunit
	  end if
          ierror=1
          return
        endif
c
	if(iprint.eq.1) then
          write(6,fmt='('' **'',12i6)') (idrec1(i),i=1,12)
          write(6,*) '**    no. of fields written: ',nwrite
        endif
c
      elseif(iprint.eq.1) then
c
        write(6,*) '**mwfelt** no fields written to file:'
        write(6,*) '** ',filnam
        write(6,*) '** file unit: ',iunit
c
      end if
c
      if(imode.eq.3 .and. ldata.gt.1) then
c..return part of record 1
	lread=min(ldata,lrec)
	do i=1,lread
	  idata(i)=idrec1(i)
	end do
      end if
c
      close(iunit)
c
      filold='****'
      lunit=0
      nwrite=0
c
      return
      end

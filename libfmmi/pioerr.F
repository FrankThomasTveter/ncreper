      subroutine pioerr(iunit,ioerr)
c
c  NAME:
c     pioerr
c
c  PURPOSE:
c     Print i/o error message after open/read/write/close or other
c     i/o failures.
c     Uses the fortran library routine gerror to retrieve the message.
c     Call this routine immediately after an error has occured as the
c     gerror routine simply retrives the last system error !
c
c  SYNOPSIS:
c     subroutine pioerr(iunit,ioerr)
c     integer iunit,ioerr
c
c  INPUT:
c     iunit  -  output unit no. (6 is the screen or 'standard out')
c     ioerr  -  error no. from open/read/write/close or other i/o
c		(iostat=ioerr)
c
c  NOTES:
c     +  Not using the fortran library routine perror as this prints
c        on 'standard error' only (not even 'standard out').
c     +  By using the fortran library routine gerror a lot of coding
c	 is avoided, but with limited functionality as result.
c        It only prints the last system error message and can't always
c	 explain an error that happened some time ago.
c
c---------------------------------------------------------------------
c  DNMI/FoU  03.12.1993  Anstein Foss
c---------------------------------------------------------------------
c
c..input/output
      integer iunit,ioerr
c
c..local
      character*256 errmsg
c
      call gerror(errmsg)
c
      lm=len(errmsg)
      kk=1
      do k=1,lm
	if(errmsg(k:k).ne.' ') kk=k
      end do
c
      ks=70
      k2=0
      do while (k2.lt.kk)
	k1=k2+1
	k2=k2+ks
	if(k2.lt.kk) then
	  k2=k2+1
	  do while (errmsg(k2:k2).ne.' ' .and. k2.gt.k1)
	    k2=k2-1
          end do
	  k2=k2-1
	  if(k2.le.k1) k2=k1+ks-1
	else
	  k2=kk
	end if
	write(iunit,*,iostat=ios) errmsg(k1:k2)
      end do
c
      return
      end

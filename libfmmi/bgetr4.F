      subroutine bgetr4(iunit,irec,buff,lbuff,ibuff,dat,ndat,ierror)
c
c  NAME:
c     bgetr4
c
c  PURPOSE:
c     Buffered file get of real*4 data.
c     The data is copied from an integer*2 buffer without scaling.
c     When all of the buffer is used a new record is read from the file.
c     The file is an unformatted direct access file with record length
c     equal buffer length.
c
c  SYNOPSIS:
c     subroutine bgetr4(iunit,irec,buff,lbuff,ibuff,dat,ndat,ierror)
c     integer   iunit,irec,lbuff,ibuff,ndat,ierror
c     integer*2 buff(lbuff)
c     real*4    dat(ndat)
c
c  INPUT:
c     iunit  - the file unit
c     lbuff  - buffer and file record length
c     ndat   - length of data array
c
c  INPUT/OUTPUT:
c     irec   - last record no. read
c              (initialize to zero before first call to a bget* routine)
c     buff   - buffer
c     ibuff  - last used element in the buffer
c              (initialize to zero before first call to a bget* routine)
c
c  OUTPUT:
c     dat    - data array retrieved from buffer.
c     ierror - error status: 0 = no error
c			     other value is an i/o error
c
c-----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.198x  Anstein Foss
c  DNMI/FoU  06.04.1995  Anstein Foss ... new library version
c-----------------------------------------------------------------------
c
      implicit none
c
      integer   iunit,irec,lbuff,ibuff,ndat,ierror
      integer*2 buff(lbuff)
      real*4    dat(ndat)
c
      integer   i,i1,i2
c
      ierror=0
      i2=0
c
      do while(i2.lt.ndat .and. ierror.eq.0)
	if(ibuff.eq.lbuff .or. irec.eq.0)
     +     call bgetfi(iunit,irec,buff,lbuff,ibuff,ierror)
	if(ierror.eq.0) then
	  i1=i2+1
	  i2=min(i2+lbuff-ibuff,ndat)
	  do i=i1,i2
	    ibuff=ibuff+1
	    dat(i)=buff(ibuff)
	  end do
	end if
      end do
c
      return
      end

      subroutine bgetjp(iunit,irec,buff,lbuff,ibuff,ndat,ierror)
c
c  NAME:
c     bgetjp
c
c  PURPOSE:
c     Jump over (skip) words in the integer*2 buffer. When all of
c     the buffer is skipped a new record is read from the file.
c     The file is an unformatted direct access file with record length
c     equal buffer length.
c
c  SYNOPSIS:
c     subroutine bgetjp(iunit,irec,buff,lbuff,ibuff,ndat,ierror)
c     integer   iunit,irec,lbuff,ibuff,ndat,ierror
c     integer*2 buff(lbuff)
c
c  INPUT:
c     iunit  - the file unit
c     lbuff  - buffer and file record length
c     ndat   - no. of words to skip in the buffer.
c
c  INPUT/OUTPUT:
c     irec   - last record no. read
c              (initialize to zero before first call to a bget* routine)
c     buff   - buffer
c     ibuff  - last used element in the buffer
c              (initialize to zero before first call to a bget* routine)
c
c  OUTPUT:
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
c
      integer   i1,i2
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
	  ibuff=ibuff+(i2-i1+1)
	end if
      end do
c
      return
      end

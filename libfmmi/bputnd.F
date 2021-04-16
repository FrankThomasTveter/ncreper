      subroutine bputnd(iunit,irec,buff,lbuff,ibuff,konvd,ierror)
c
c  NAME:
c     bputnd
c
c  PURPOSE:
c     Buffered file put of integer*2 data.
c     Fill remaining part of buffer with zeros and write it to the file.
c     Usually called when writing the last buffer to file.
c     The file is an unformatted direct access file with record length
c     equal buffer length.
c
c  SYNOPSIS:
c     subroutine bputnd(iunit,irec,buff,lbuff,ibuff,konvd,ierror)
c     integer   iunit,irec,lbuff,ibuff,konvd,ierror
c     integer*2 buff(lbuff)
c
c  INPUT:
c     iunit  - the file unit
c     lbuff  - buffer and file record length
c     konvd  - data conversion: 0 = no conversion
c			        1 = swap bytes (for pc,dec,vax,...)
c
c  INPUT/OUTPUT:
c     irec   - last record no. written
c              (initialize to zero before first call to a bput* routine)
c     buff   - buffer
c     ibuff  - last used element in the buffer,
c	       if ibuff equals zero, nothing is done
c              (initialize to zero before first call to a bput* routine)
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
      integer   iunit,irec,lbuff,ibuff,konvd,ierror
      integer*2 buff(lbuff)
c
      integer   i
c
      ierror=0
c
      if(ibuff.gt.0) then
c
        do i=ibuff+1,lbuff
	  buff(i)=0
        end do
        ibuff=lbuff
c
        call bputfi(iunit,irec,buff,lbuff,ibuff,konvd,ierror)
c
      end if
c
      return
      end

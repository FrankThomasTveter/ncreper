      subroutine bputi4(iunit,irec,buff,lbuff,ibuff,konvd,dat,ndat,
     +			ierror)
c
c  NAME:
c     bputi4
c
c  PURPOSE:
c     Buffered file put of integer*4 data.
c     The data is copied to an integer*2 buffer without scaling.
c     When a buffer is full it is written to the file.
c     The file is an unformatted direct access file with record length
c     equal buffer length.
c
c  SYNOPSIS:
c     subroutine bputi4(iunit,irec,buff,lbuff,ibuff,konvd,dat,ndat,
c    +			ierror)
c     integer   iunit,irec,lbuff,ibuff,konvd,ndat,ierror
c     integer*2 buff(lbuff)
c     integer*4 dat(ndat)
c
c  INPUT:
c     iunit  - the file unit
c     lbuff  - buffer and file record length
c     konvd  - data conversion: 0 = no conversion
c			        1 = swap bytes (for pc,dec,vax,...)
c     dat    - data array for transfer to bufffer
c     ndat   - length of data array
c
c  INPUT/OUTPUT:
c     irec   - last record no. written
c              (initialize to zero before first call to a bput* routine)
c     buff   - buffer
c     ibuff  - last used element in the buffer
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
      integer   iunit,irec,lbuff,ibuff,konvd,ndat,ierror
      integer*2 buff(lbuff)
      integer*4 dat(ndat)
c
      integer   i,i1,i2
c
      ierror=0
      i2=0
c
      do while(i2.lt.ndat .and. ierror.eq.0)
	i1=i2+1
	i2=min(i2+lbuff-ibuff,ndat)
	do i=i1,i2
	  ibuff=ibuff+1
	  buff(ibuff)=dat(i)
	end do
	if(ibuff.eq.lbuff)
     +     call bputfi(iunit,irec,buff,lbuff,ibuff,konvd,ierror)
      end do
c
      return
      end

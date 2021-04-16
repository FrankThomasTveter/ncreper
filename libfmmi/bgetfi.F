      subroutine bgetfi(iunit,irec,buff,lbuff,ibuff,ierror)
c
c  NAME:
c     bgetfi
c
c  PURPOSE:
c     Buffered file get of a integer*2 buffer from a file.
c     The file is an unformatted direct access file with record length
c     equal buffer length.
c
c  SYNOPSIS:
c     subroutine bgetfi(iunit,irec,buff,lbuff,ibuff,ierror)
c     integer   iunit,irec,lbuff,ibuff,ierror
c     integer*2 buff(lbuff)
c
c  INPUT:
c     iunit  - the file unit
c     lbuff  - buffer and file record length
c
c  INPUT/OUTPUT:
c     irec   - last record no. written (increased by 1)
c     ibuff  - last used element in the buffer,
c	       i.e. buffer and file record length,
c	       reset to 0 here.
c
c  OUTPUT:
c     buff   - buffer
c     ierror - error status: 0 = no error
c			     other value is an i/o error
c
c-----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.198x  Anstein Foss
c  DNMI/FoU  06.04.1995  Anstein Foss ... new library version
c  DNMI/FoU  12.06.2003  Anstein Foss ... automatic byteswap
c-----------------------------------------------------------------------
c
      implicit none
c
      integer   iunit,irec,lbuff,ibuff,ierror
      integer*2 buff(lbuff)
c
      integer   ios
c
      logical swapfile,swap
c
      ierror=0
c
      if(irec.eq.0) swap=swapfile(-iunit)
c
      irec=irec+1
      read(iunit,rec=irec,iostat=ios) buff
      if(ios.ne.0) then
	write(*,*) '**bgetfi** read error. file unit: ',iunit
	write(*,*) '**                        record: ',irec
	write(*,*) '**                        iostat: ',ios
	ierror=ios
      elseif(swapfile(iunit)) then
	call bswap2(lbuff,buff)
      end if    
c
      ibuff=0
c
      return
      end

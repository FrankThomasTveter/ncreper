      subroutine bputfi(iunit,irec,buff,lbuff,ibuff,konvd,ierror)
c
c  NAME:
c     bputfi
c
c  PURPOSE:
c     Buffered file put of the integer*2 buffer to the file.
c     The file is an unformatted direct access file with record length
c     equal buffer length.
c
c  SYNOPSIS:
c     subroutine bputfi(iunit,irec,buff,lbuff,ibuff,konvd,ierror)
c     integer   iunit,irec,lbuff,ibuff,konvd,ierror
c     integer*2 buff(ibuff)
c
c  INPUT:
c     iunit  - the file unit
c     buff   - buffer
c     lbuff  - buffer and file record length
c     konvd  - data conversion: 0 = no conversion
c			        1 = swap bytes (for pc,dec,vax,...)
c
c  INPUT/OUTPUT:
c     irec   - last record no. written (increased by 1)
c     ibuff  - last used element in the buffer,
c	       i.e. buffer and file record length,
c	       reset to 0 here.
c
c  OUTPUT:
c     ierror - error status: 0 = no error
c			     other value is an i/o error
c
c-----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.198x  Anstein Foss
c  DNMI/FoU  06.04.1995  Anstein Foss ... new library version
c  DNMI/FoU  29.11.1995  Anstein Foss ... IBM/RS 6000 xl fortran
c-----------------------------------------------------------------------
c
      implicit none
c
      integer   iunit,irec,lbuff,ibuff,konvd,ierror
      integer*2 buff(lbuff)
c
      integer   i,ios
      integer*2 i1,i2,nn8,nnm8,nn255
c
      ierror=0
      nn8=+8
      nnm8=-8
      nn255=255
c
      if(konvd.eq.1) then
c..swap bytes
	do i=1,lbuff
	  i1=ishft(iand(buff(i),nn255),nn8)
	  i2=iand(ishft(buff(i),nnm8),nn255)
	  buff(i)=ior(i1,i2)
	end do
      end if
c
      irec=irec+1
      write(iunit,rec=irec,iostat=ios) buff
      if(ios.ne.0) then
	write(6,*) '**bputfi** write error. file unit: ',iunit
	write(6,*) '**                         record: ',irec
	write(6,*) '**                         iostat: ',ios
	ierror=ios
      end if
c
      ibuff=0
c
      return
      end

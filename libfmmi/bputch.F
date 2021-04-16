      subroutine bputch(iunit,irec,buff,lbuff,ibuff,konvd,konvc,
     +			text,ntext,nchar,ierror)
c
c  NAME:
c     bputi2
c
c  PURPOSE:
c     Buffered file put of character data in text string(s).
c     Two characters are copied to one integer*2 word in the buffer.
c     If the no. of characters are not modulo two, the last
c     character will be lost.
c     When a buffer is full it is written to the file.
c     The file is an unformatted direct access file with record length
c     equal buffer length.
c
c  SYNOPSIS:
c     subroutine bputch(iunit,irec,buff,lbuff,ibuff,konvd,konvc,
c    +			text,ntext,nchar,ierror)
c     integer   iunit,irec,lbuff,ibuff,konvd,konvc,ntext,nchar,ierror
c     integer*2 buff(lbuff)
c     character*(*) text(ntext)
c
c  INPUT:
c     iunit  - the file unit
c     lbuff  - buffer and file record length
c     konvd  - data conversion: 0 = no conversion
c			        1 = swap bytes (for pc,dec,vax,...)
c     konvc  - character conversion: 0 = no conversion
c				     1 = ascii for pc
c				     2 = ebcdic (for ibm)
c			            -1 = ascii for pc, and resetting
c				         byteswapping due to konvd=1,
c					 (for some old pc programs)
c     text   - character string(s) for transfer to bufffer
c     ntext  - no. of character strings
c     nchar  - no. of characters in each text string
c	       (-1 = all characters in each string) 
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
c  DNMI/FoU  29.11.1995  Anstein Foss ... IBM/RS 6000 xl fortran
c-----------------------------------------------------------------------
c
      implicit none
c
      integer   iunit,irec,lbuff,ibuff,konvd,konvc,ntext,nchar,ierror
      integer*2 buff(lbuff)
      character*(*) text(ntext)
c
      integer    mcbuff
      parameter (mcbuff=128)
      character*(mcbuff*2) cbuff
c
      integer   nchr,nchtwo,itext,i,i1,i2,k1,k2,nc
      integer*2 ic1,ic2,nn8
      character*1 chr
c
      ierror=0
      nchr=nchar
      if(nchr.lt.0) nchr=len(text(1))
      nchtwo=nchr/2
      itext=0
      nn8=+8
c
      do while(itext.lt.ntext .and. ierror.eq.0)
	itext=itext+1
        i2=0
c
        do while(i2.lt.nchtwo .and. ierror.eq.0)
c
	  i1=i2+1
	  i2=min(i2+lbuff-ibuff,mcbuff,nchtwo)
	  k1=i1*2-1
	  k2=i2*2
	  nc=0
	  do i=k1,k2
	    nc=nc+1
	    cbuff(nc:nc)=text(itext)(i:i)
	  end do
c
	  if(konvc.eq.1 .or. konvc.eq.-1) then
c..ascii (standard) -> ascii (pc) ......................... pc
	    call ascii1(1,cbuff,nc)
	    if(konvc.eq.-1) then
c..usually to reset swapping performed later due to konvd=1
	      do i=1,nc,2
		chr=cbuff(i:i)
		cbuff(i:i)=cbuff(i+1:i+1)
		cbuff(i+1:i+1)=chr
	      end do
	    end if
	  elseif(konvc.eq.2) then
c..ascii (standard) -> ebcdic (ibm) ....................... ibm
	    call ebcdic(1,cbuff,nc)
	  end if
c
	  do i=1,nc,2
	    ic1=ichar(cbuff(i:i))
	    ic2=ichar(cbuff(i+1:i+1))
	    ibuff=ibuff+1
	    buff(ibuff)=ior(ishft(ic1,nn8),ic2)
	  end do
	  if(ibuff.eq.lbuff)
     +       call bputfi(iunit,irec,buff,lbuff,ibuff,konvd,ierror)
c
	end do
c
      end do
c
      return
      end

      subroutine bgetch(iunit,irec,buff,lbuff,ibuff,text,ntext,nchar,
     +			ierror)
c
c  NAME:
c     bgetch
c
c  PURPOSE:
c     Buffered file get of character data in text string(s).
c     The data is copied from an integer*2 buffer without scaling.
c     Two and two characters are copied from integer*2 words in the
c     buffer. If the no. of characters are not modulo two, the last
c     character will just get a space.
c     When all of the buffer is used a new record is read from the file.
c     The file is an unformatted direct access file with record length
c     equal buffer length.
c
c  SYNOPSIS:
c     subroutine bgetch(iunit,irec,buff,lbuff,ibuff,text,ntext,nchar,
c    +			ierror)
c     integer   iunit,irec,lbuff,ibuff,ntext,nchar,ierror
c     integer*2 buff(lbuff)
c     character*(*) text(ntext)
c
c  INPUT:
c     iunit  - the file unit
c     lbuff  - buffer and file record length
c     ntext  - no. of character strings
c     nchar  - no. of characters in each text string
c	       (-1 = all characters in each string) 
c
c  INPUT/OUTPUT:
c     irec   - last record no. read
c              (initialize to zero before first call to a bget* routine)
c     buff   - buffer
c     ibuff  - last used element in the buffer
c              (initialize to zero before first call to a bget* routine)
c
c  OUTPUT:
c     text   - character string(s) retrieved from buffer.
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
      integer   iunit,irec,lbuff,ibuff,ntext,nchar,ierror
      integer*2 buff(lbuff)
      character*(*) text(ntext)
c
      integer   nchr,nchtwo,itext,i,i1,i2,k1,k2
      integer*2 ic1,ic2,nnm8,nn255
c
      ierror=0
      nchr=nchar
      if(nchr.lt.0) nchr=len(text(1))
      nchtwo=nchr/2
      itext=0
      nnm8=-8
      nn255=255
c
      do while(itext.lt.ntext .and. ierror.eq.0)
	itext=itext+1
	i2=0
c
        do while(i2.lt.nchtwo .and. ierror.eq.0)
c
	  if(ibuff.eq.lbuff .or. irec.eq.0)
     +       call bgetfi(iunit,irec,buff,lbuff,ibuff,ierror)
	  if(ierror.eq.0) then
	    i1=i2+1
	    i2=min(i2+lbuff-ibuff,nchtwo)
	    k1=i1*2-1
	    k2=i2*2
	    do i=k1,k2,2
	      ibuff=ibuff+1
	      ic1=ishft(buff(ibuff),nnm8)
	      ic2= iand(buff(ibuff),nn255)
	      text(itext)(i:i)    =char(ic1)
	      text(itext)(i+1:i+1)=char(ic2)
	    end do
	  end if
c
        end do
c
	if(ierror.eq.0 .and. nchtwo*2.lt.nchr)
     +	   text(itext)(nchr:nchr)=' '
c
      end do
c
      return
      end

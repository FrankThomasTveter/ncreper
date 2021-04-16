      subroutine bpaki2(n1a,n2a,n3a,n4a,n5a,da,
     +			n1b,n2b,n3b,n4b,n5b,db,iab)
c
c  NAME:
c     bpaki2
c
c  PURPOSE:
c     Pack or unpack an integer*2 data array.
c     Packing/unpacking here means decrease/increase dimensions of
c     an array. Useful when the calling routine uses max dimensions
c     while only used parts of the array are stored in a file.
c
c  SYNOPSIS:
c     subroutine bpaki2(n1a,n2a,n3a,n4a,n5a,da,
c    +			n1b,n2b,n3b,n4b,n5b,db,iab)
c     integer   n1a,n2a,n3a,n4a,n5a
c     integer   n1b,n2b,n3b,n4b,n5b,iab
c     integer*2 da(n1a,n2a,n3a,n4a,n5a)
c     integer*2 db(n1b,n2b,n3b,n4b,n5b)
c
c  INPUT:
c     n1a    - first  dimension of the unpacked array da
c     n2a    - second dimension of the unpacked array da
c     n3a    - third  dimension of the unpacked array da (possibly 1)
c     n4a    - fourth dimension of the unpacked array da (possibly 1)
c     n5a    - fifth  dimension of the unpacked array da (possibly 1)
c     n1b    - first  dimension of the packed array db
c     n2b    - second dimension of the packed array db
c     n3b    - third  dimension of the packed array db (possibly 1)
c     n4b    - fourth dimension of the packed array db (possibly 1)
c     n5b    - fifth  dimension of the packed array db (possibly 1)
c     iab    - +1 = packing   (copy contents of da to db)
c	       -1 = unpacking (copy contents of db to da)
c
c  INPUT/OUTPUT:
c     da     - input  array if iab=+1,  output array if iab=-1
c     db     - output array if iab=+1,  input  array if iab=-1
c
c  NOTES:
c     da is always the unpacked array.
c     db is always the   packed array.
c     Requirements (not checked in this routine):
c     n1a >= n1b,  n2a >= n2b,  n3a >= n3b,  n4a >= n4b,  n5a >= n5b.
c
c-----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.198x  Anstein Foss
c  DNMI/FoU  06.04.1995  Anstein Foss ... new library version
c-----------------------------------------------------------------------
c
      implicit none
c
      integer   n1a,n2a,n3a,n4a,n5a
      integer   n1b,n2b,n3b,n4b,n5b,iab
      integer*2 da(n1a,n2a,n3a,n4a,n5a)
      integer*2 db(n1b,n2b,n3b,n4b,n5b)
c
      integer   i1,i2,i3,i4,i5
c
      if(iab.eq.+1) then
c
	do i5=1,n5b
	  do i4=1,n4b
	    do i3=1,n3b
	      do i2=1,n2b
		do i1=1,n1b
		  db(i1,i2,i3,i4,i5)=da(i1,i2,i3,i4,i5)
		end do
	      end do
	    end do
	  end do
	end do
c
      elseif(iab.eq.-1) then
c
	do i5=n5b,1,-1
	  do i4=n4b,1,-1
	    do i3=n3b,1,-1
	      do i2=n2b,1,-1
		do i1=n1b,1,-1
		  da(i1,i2,i3,i4,i5)=db(i1,i2,i3,i4,i5)
		end do
	      end do
	    end do
	  end do
	end do
c
      end if
c
      return
      end


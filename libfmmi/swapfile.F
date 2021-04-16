      logical function swapfile(iunit)
c
c	returns true (swap) or false (not swap)
c
c	input: iunit<0 : reads first word from already open file,
c		         and checks known filetypes
c	       iunit>0 : returns status from a previous check.
c
c	Known filetypes (with access='direct' and integer*2 data):
c	  DNMI Felt file:                 first word = 999,998,997
c         DNMI Observation files:         first word = 1,2,...15
c         DNMI Vertical profile data:     first word = 201
c         DNMI Vertical crossection data: first word = 121
c         DNMI Wave spectrum data:        first word = 251
c         DNMI Meteogram data:            first word = 211 (metdat/posdat)
c         DNMI Wave/sea diagram data:     first word = 212 (metdat/posdat)
c         DNMI misc diagram data resorted:first word = 221 (posdat)
c       Unknown filetypes always of the "not swap" type
c	this is also default for all possible file units.
c
c	It is the programmers responsability to check the file before
c	requesting the swap status.
c
c	A minimal (unsafe) solution when Linux appeared at DNMI
c
c-----------------------------------------------------------------------
c  DNMI/FoU  09.11.2000  Anstein Foss
c  DNMI/FoU  13.01.2003  Anstein Foss  ... +obsfiles 11-15
c  DNMI/FoU  22.03.2003  Anstein Foss  ... +wave spectrum file
c  DNMI/FoU  12.06.2003  Anstein Foss  ... +meteogram data (211,212,221)
c-----------------------------------------------------------------------
c
      implicit none
c
      integer iunit
c
      integer maxunit,maxkey
      parameter (maxunit=256,maxkey=3+15+6)
c
      logical   swap(maxunit)
      logical   first
      integer*2 keys(maxkey,2)
      integer   i,n,iu,ios,keytype
      integer*2 key
c
      save swap,keys
      data first /.true./
      data keys/999,998,997, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     +		201, 121, 251, 211, 212, 221, 24*0/
c
      if (first) then
	do i=1,maxunit
	  swap(i)= .false.
	end do
	do i=1,maxkey
	  keys(i,2)= keys(i,1)
	end do
	call bswap2(maxkey,keys(1,2))
c----------------------------------------------
c	do i=1,maxkey
c	  write(6,*) keys(i,1),keys(i,2)
c	end do
c----------------------------------------------
	first= .false.
      end if
c
      if (iunit.gt.0 .and. iunit.le.maxunit) then
	swapfile= swap(iunit)
      elseif (iunit.lt.0 .and. iunit.ge.-maxunit) then
	iu= -iunit
	swap(iu)= .false.
	read(iu,rec=1,iostat=ios) key
	if (ios.eq.0) then
	  keytype=0
	  do n=1,2
	    do i=1,maxkey
	      if (keys(i,n).eq.key) keytype=n
	    end do
	  end do
	  if (keytype.eq.2) then 
	    swap(iu)= .true.
	  elseif (keytype.eq.0) then
	    write(6,*) 'swapfile: unknown file type, unit=',iu
	  end if
	else
	  write(6,*) 'swapfile: read error, unit=',iu
	end if
	swapfile= swap(iu)
      else
	write(6,*) 'swapfile: bad unit=',iunit
	swapfile= .false.
      end if
c
      return
      end

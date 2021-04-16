      subroutine extendf2(a,b,nx,ny,nxwrap,nxext1,nxext2,nyext1,nyext2)
c
c  extendf2 - extend a(i,j) near the boundaries for
c  hemispheric spherical fields; cyclic in x direction
c  and/or extentions in y-directions
c
      implicit none
c
      integer nx,ny,nxwrap,nxext1,nxext2,nyext1,nyext2
      real    a(nx,ny),b(1-nxext1:nx+nxext2,1-nyext1:ny+nyext2)
c
      integer i,j
c
      do j = 1,ny
      do i = 1,nx
         b(i,j) = a(i,j)
      enddo
      enddo
c
c  west cyclic
      if (nxext1.gt.0 .and. nxwrap.gt.0) then
         do j = 1,ny
	    do i = 1-nxext1,0
               b(i,j) = a(i+nxwrap,j)
	    end do
         enddo
      endif
c
c  east cyclic
      if (nxext2.gt.0 .and. nxwrap.gt.0) then
         do j = 1,ny
	    do i = nx+1,nx+nxext2
               b(i,j) = a(i-nxwrap,j)
	    end do
         enddo
      endif
c
c  "south" extension (may as well be north)
      if (nyext1.gt.0) then
	 do j=1-nyext1,0
            do i = 1-nxext1,nx+nxext2
               b(i,j) = b(i,1)
	    end do
	 end do
      endif
c
c  "north" extension (may as well be south)
      if (nyext2.gt.0) then
	 do j=ny+1,ny+nyext2
            do i = 1-nxext1,nx+nxext2
               b(i,j) = b(i,ny)
	    end do
	 end do
      endif
c
      return
      end

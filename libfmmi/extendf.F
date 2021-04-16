      subroutine extendf(a,b,nx,ny)
c
c  extendf - extend a(i,j) near the boundaries for
c  hemispheric spherical fields; cyclic in x direction
c  and extention in y-direction
c
      integer nx,ny
      real a(nx,ny),b(-2:nx+3,-2:ny+3)
c
      do j = 1,ny
      do i = 1,nx
         b(i,j) = a(i,j)
      enddo
      enddo
c
c  west-east cyclic field
      do j = 1,ny
         b(-2,j) = a(nx-2,j)
         b(-1,j) = a(nx-1,j)
         b(0,j) = a(nx,j)
         b(nx+1,j) = a(1,j)
         b(nx+2,j) = a(2,j)
         b(nx+3,j) = a(3,j)
      enddo
c
c  north-south extension
      do i = -2,nx+3
         b(i,-2) = b(i,1)
         b(i,-1) = b(i,1)
         b(i,0) = b(i,1)
         b(i,ny+1) = b(i,ny)
         b(i,ny+2) = b(i,ny)
         b(i,ny+3) = b(i,ny)
      enddo
c
      return
      end

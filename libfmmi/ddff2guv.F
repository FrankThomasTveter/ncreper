      subroutine ddff2guv(nx,ny,u,v,iturn,igtypr,gr,udef,ierror)
c
c****************************************************************
c
c     ddff2guv - convert true north velocity direction and force
c                to x and y grid components
c
c  purpose:
c
c     convert true north velocity direction and force
c     to grid x and y components in different grids:
c     spherical, spherical rotated , polar stereographic or mercator.
c
c  input/output parameters:
c
c     nx     - no. of grid points in x-direction
c     ny     - no. of grid points in y-direction
c     u      - input:  true north direction (dd, degrees)
c              output: vector (velocity) component in x-direction
c     v      - input:  velocity
c              output: vector (velocity) component in y-direction
c     iturn  - 0 : dd is the 'towards' direction (waves,current etc.)
c              1 : dd is the 'from'    direction (wind in meteorology)
c     udef   - undefined (missing) value
c     igtypr - output grid type
c     gr     - output grid description
c     ierror - output error status, 0=no error
c
c  description of g = gr (for igtype = igtypr):
c
c  for spherical (rotated) grid (igtype=2,3):
c
c     g(1) - western boundary (degrees)
c     g(2) - southern boundary (degrees)
c     g(3) - longitude increment (degrees)
c     g(4) - latitude increment (degrees)
c     g(5) - xcen: longitude position of rotated equator (degrees)
c     g(6) - ycen: latitude  position of rotated equator (degrees)
c            (lamda,theta)=(xcen,ycen) at (lamda',theta')=(0,0),
c            where (lamda,theta) are usual spherical coord. and
c            (lamda',theta') are rotated spherical coord.
c            xcen = ycen = 0 for usual spherical coord.
c
c  for polar stereographic grid (igtype=1,4):
c
c     g(1) - x-position of north pole
c     g(2) - y-position of north pole
c     g(3) - number of grid distances between pole and equator
c     g(4) - rotation angle of the grid (degrees)
c     g(5) - projection latitude (degrees)
c            (60 degrees north for igtype=1)
c     g(6) - 0. (not used)
c
c  for mercator (unrotated) grid, igtype=5:
c
c     g(1) - western boundary (longitude for x=1) (degrees)
c     g(2) - southern boundary (latitude for y=1) (degrees)
c     g(3) - x (longitude) increment (km)
c     g(4) - y (latitude)  increment (km)
c     g(5) - reference (construction) latitude (degrees)
c     g(6) - 0.  (not used)
c
c  for unknown grid types, igtype=0, it is assumed that
c  'true north' is along the (positive) y axis.
c
c  externals:
c
c     uvconvert - turn vector (velocity) components between grids
c
c  history:
c
c     j.e. haugen/dnmi      nov -94 ... grv2grv
c     a.   foss   dnmi   02.02.1995 ... no size limits
c     a.   foss   dnmi   25.08.1995 ... uvconvert
c     a.   foss   dnmi   28.08.1995 ... ddff2guv
c     a.   foss   dnmi   12.06.1996 ... calling uvconvert
c
c****************************************************************
c
      implicit none
c
      integer nx, ny, iturn, igtypr, ierror
      real u(nx*ny), v(nx*ny),
     +     gr(6), udef
c
      integer nmax
c
      parameter ( nmax = 1000 )
c
      real x1(nmax), y1(nmax)
c
      real    gg(6), zpir18, zturn, dd, ff
      integer igtypg, j, j0, j1, j2, nxy, ix, iy
c
      data igtypg/2/
      data gg/1.,1.,1.,1.,0.,0./
c
      zpir18 = 2.0*asin(1.0)/180.
c
      ierror = 0
c
      zturn=+1.
      if (iturn.ne.0) zturn=-1.
c
      do j = 1,nx*ny
         if (u(j).ne.udef .and. v(j).ne.udef) then
            dd = u(j)*zpir18
            ff = v(j)*zturn
            u(j) = ff*sin(dd)
            v(j) = ff*cos(dd)
         else
            u(j) = udef
            v(j) = udef
         endif
      end do
c
      if (igtypr.eq.0) return
c
      do j1 = 1,nx*ny,nmax
c
      j0 = j1 - 1
      j2 = min(j0+nmax,nx*ny)
      nxy = j2 - j0
c
      do j = 1,nxy
         iy = (j0+j+nx-1)/nx
         ix = (j0+j)-(iy-1)*nx
         x1(j) = float(ix)
         y1(j) = float(iy)
      end do
c
      call uvconvert(nxy,x1(1),y1(1),u(j1),v(j1),
     +               igtypg,gg,igtypr,gr,udef,ierror)
      if(ierror.ne.0) return
c
      end do
c
      return
      end

      subroutine pos2pos (igtypa,ga,
     +                    npos,pos,geopos,gxypos,ierror)
c
c****************************************************************
c
c     pos2pos - coordinate conversions as in grd2pos
c
c  purpose:
c
c     coordinate conversions in different grids; spherical,
c     spherical rotated, polar stereographic or mercator.
c
c  input/output parameters:
c
c     igtypa - grid type
c     ga     - grid parameters
c     npos   - no. of positions
c     pos    - input positions,
c		    pos(1,n) : position type, 
c				  1 = x,y position in grid
c			          2 = geographic latitude,longitude
c				      (decimal)
c				 (any other value treated as 1)
c		    pos(2,n) : x or latitude  position
c		    pos(3,n) : y or longitude position
c     geopos- output positions,
c		    geopos(1,n) : geographic latitude
c		    geopos(2,n) : geographic longitude
c     gxypos- output positions,
c		    gxypos(1,n) : x position in grid
c		    gxypos(2,n) : y position in grid
c     ierror - exit status, 0=no error
c
c  description of g = ga (for igtype = igtypa):
c
c  for spherical (rotated) grid, igtype=2,3:
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
c  for polar stereographic grid, igtype=1,4:
c
c     g(1) - x-position of north pole
c     g(2) - y-position of north pole
c     g(3) - number of grid units between pole and equator
c     g(4) - rotation angle of the grid (degrees)
c     g(5) - projection latitude (degrees)
c	     (usually 60 degrees north)
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
c  for unknown grid type, igtype=0:
c     the g array is not used,
c     geographic positions (input) is not allowed.
c
c  externals:
c
c     xyconvert - performs all coordinate conversions
c
c  history:
c
c     j.e. haugen/dnmi   nov -94 ...... grd2grd
c     a.   foss   dnmi   02.02.1995 ... no size limits, interp speedup
c     a.   foss   dnmi   06.02.1995 ... grd2pos
c     a.   foss   dnmi   25.04.1995 ... pos2pos
c     a.   foss   dnmi   12.06.1996 ... new call, using xyconvert,...
c
c****************************************************************
c
      implicit none
c
      integer igtypa, npos, ierror
      real pos(3,npos), geopos(2,npos), gxypos(2,npos), ga(6)
c
      integer nmax
c
      parameter ( nmax = 1000 )
c
      real x1(nmax), y1(nmax)
c
      integer igtypg, ipos, j, n, j0, j2
      real zmvx, zmvy
      real gg(6)
c
      data igtypg/2/
      data gg/1.,1.,1.,1.,0.,0./
c
        zmvx = 0.
        zmvy = 0.
c
	do j = 1,npos
	   geopos(1,j) = pos(2,j)
	   geopos(2,j) = pos(3,j)
	   gxypos(1,j) = pos(2,j)
	   gxypos(2,j) = pos(3,j)
	enddo
c
c convert geographic latitude,longitude coordinates to x,y coordinates
	j = 0
c
	do while (j.lt.npos)
c
	j0 = j
	n  = 0
c
	do while (n.lt.nmax .and. j.lt.npos)
	   j = j + 1
	   ipos = nint(pos(1,j))
	   if (ipos.eq.2) then
	      n = n + 1
              x1(n) = pos(3,j)
              y1(n) = pos(2,j)
	   endif
	enddo
	j2=j
c
	if (n.gt.0) then
c
	  call xyconvert(n,x1(1),y1(1),igtypg,gg,igtypa,ga,ierror)
	  if(ierror.ne.0) then
	    write(6,*) '**POS2POS** XYCONVERT ERROR: ',ierror
	    return
	  end if
c
	  j = j0
	  n = 0
c
	  do while (j.lt.j2)
	     j = j + 1
	     ipos = nint(pos(1,j))
	     if (ipos.eq.2) then
	        n = n + 1
                gxypos(1,j) = x1(n) - zmvx
                gxypos(2,j) = y1(n) - zmvy
	     endif
	  enddo
c
	endif
c
	enddo
c
	if (igtypa.ne.0) then
c
c convert x,y coordinates to geographic coordinates
	j = 0
c
	do while (j.lt.npos)
c
	j0 = j
	n  = 0
c
	do while (n.lt.nmax .and. j.lt.npos)
	   j = j + 1
	   ipos = nint(pos(1,j))
	   if (ipos.eq.1) then
	      n = n + 1
              x1(n) = pos(2,j) + zmvx
              y1(n) = pos(3,j) + zmvy
	   endif
	enddo
	j2=j
c
	if (n.gt.0) then
c
	  call xyconvert(n,x1(1),y1(1),igtypa,ga,igtypg,gg,ierror)
	  if(ierror.ne.0) then
	    write(6,*) '**POS2POS** XYCONVERT ERROR: ',ierror
	    return
	  end if
c
c  check xyconvert output 
	  do j = 1,n
	    if (x1(j).le.-180.) x1(j) = x1(j) + 360.
	    if (x1(j).gt.+180.) x1(j) = x1(j) - 360.
	  end do
c
	  j = j0
	  n = 0
c
	  do while (j.lt.j2)
	     j = j + 1
	     ipos = nint(pos(1,j))
	     if (ipos.eq.1) then
	        n = n + 1
                geopos(1,j) = y1(n)
                geopos(2,j) = x1(n)
	     endif
	  enddo
c
	endif
c
	enddo
c
	endif
c
      return
      end

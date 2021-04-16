      subroutine movegrid(igtype,grid,dx,dy,gridmv,ierror)
c
c  NAME:
c     movegrid
c
c  PURPOSE:
c     Adjust grid parameters to a new origo.
c     May be used when only a part of a field is handled and other
c     routines can't handle a sub area and need the grid parameters.
c     The origo (or rather position x=1,y=1) is moved dx,dy grid units
c     in x and y direction.
c
c  SYNOPSIS:
c     subroutine movegrid(igtype,grid,dx,dy,gridmv,ierror)
c     integer igtype,ierror
c     real    grid(6),dx,dy,gridmv(6)
c
c  INPUT:
c     igtype - grid type, 1 = polarstereographic grid (true at 60 deg. N)
c			  2 = geographic
c			  3 = spherical rotated grid
c			  4 = polarstereographic grid
c			  5 = mercator (unrotated) grid
c     grid   - grid parameters
c                igtype=1,4, polarstereographic grid:
c			   grid(1) = x position of north pole (xp)
c			   grid(2) = y position of north pole (yp)
c			   grid(3) = no. of grid units between
c				     North Pole and Equator
c			   grid(4) = grid rotation angle (degrees),
c				     positive east, negative west
c			   grid(5) = projection latitude
c				     (degrees) standard is 60 (60 deg. N)
c			   grid(6) = 0.  (not used)
c                igtype=2,3, geographic or spherical rotated grid:
c			   grid(1) = western boundary (degrees)
c				     (longitude for x=1)
c			   grid(2) = southern boundary (degrees)
c				     (latitude for y=1)
c			   grid(3) = longitude increment (degrees)
c			   grid(4) = latitude  increment (degrees)
c			   grid(5) = longitude position of rotated equator
c				     (degrees)  (0 if geographic grid)
c			   grid(6) = latitude  position of rotated equator
c				     (degrees)  (0 if geographic grid)
c                igtype=5, mercator (unrotated) grid:
c			   grid(1) = western boundary (degrees)
c				     (longitude for x=1)
c			   grid(2) = southern boundary (degrees)
c				     (latitude for y=1)
c			   grid(3) = x (longitude) increment (km)
c			   grid(4) = y (latitude)  increment (km)
c			   grid(5) = reference (construction) latitude
c                                    (degrees)
c			   grid(6) = 0.  (not used)
c     dx     - no. of gridpoints move in x direction
c     dy     - no. of gridpoints move in y direction
c
c  OUTPUT:
c     gridmv - grid parameters after move
c	       (same description as for input grid)
c     ierror - error status: 0 = no error
c			     1 = unknown grid type (igtype)
c
c-----------------------------------------------------------------------
c  DNMI/FoU  04.06.1996  Anstein Foss
c  DNMI/FoU  18.05.2001  Anstein Foss ... fixed mercator bug
c-----------------------------------------------------------------------
c
      integer igtype,ierror
      real    grid(6),dx,dy,gridmv(6)
c
      ierror=0
c
      do i=1,6
	gridmv(i)=grid(i)
      end do
c
      if(dx.eq.0. .and. dy.eq.0.) then
	continue
      elseif(igtype.eq.1 .or. igtype.eq.4) then
	gridmv(1)=gridmv(1)-dx
	gridmv(2)=gridmv(2)-dy
      elseif(igtype.eq.2 .or. igtype.eq.3) then
	gridmv(1)=gridmv(1)+gridmv(3)*dx
	gridmv(2)=gridmv(2)+gridmv(4)*dy
      elseif(igtype.eq.5) then
	zpir18 = 2.0*asin(1.0)/180.
	xwa = gridmv(1)*zpir18
	ysa = gridmv(2)*zpir18
	dxa = gridmv(3)*1000.
	dya = gridmv(4)*1000.
	yca = gridmv(5)*zpir18
	x   = 1.+dx
	y   = 1.+dy
	call mer2sph(+1,1,x,y,xwa,ysa,dxa,dya,yca,ierror)
	gridmv(1)=x/zpir18
	gridmv(2)=y/zpir18
      else
	ierror=1
      end if
c
      return
      end

      subroutine grv2pos (icall,igtypa,ga,
     +                    ur,vr,npos,vturn,geopos,udef,ierror)
c
c****************************************************************
c
c     grv2pos - turn velocity components for positions to true
c		east/west and north/south components
c		(i.e. components as in a unrotated spherical grid)
c
c  purpose:
c
c     turn velocity components from different grids; spherical,
c     spherical rotated, polar stereographic or mercator.
c
c  input/output parameters:
c
c     icall  - switch for preparations/turning
c              = 1: preparations
c              = 2: preparations and turning
c              = 3: turning
c     igtypa - input grid type
c     ga     - input grid parameters
c     ur     - input  velocity component in grid x-direction,
c              output velocity component in east/west-direction
c     vr     - input  velocity component in grid y-direction,
c              output velocity component in north/south-direction
c     npos   - no. of positions
c     vturn  - preparation matrix, vturn(4,npos)
c     geopos - matrix as returned from grd2pos or pos2pos
c                   geopos(1,n)=latitude
c                   geopos(2,n)=longitude
c     udef   - value in 'undefined' points
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
c     no turning is done.
c
c  externals:
c
c     uvconvert - here used to compute constants for the turning
c
c  history:
c
c     j.e. haugen/dnmi      nov -94 ... grv2grv
c     a.   foss   dnmi   02.02.1995 ... no size limits
c     a.   foss   dnmi   06.03.1995 ... grv2pos
c     a.   foss   dnmi   12.06.1996 ... new call, using uvconvert,...
c
c****************************************************************
c
      implicit none
c
      integer icall, igtypa, npos, ierror
      real vturn(4,npos), ur(npos), vr(npos), geopos(2,npos),
     +     ga(6), udef
c
      integer nmax
c
      parameter ( nmax = 1000 )
c
      real x1(nmax), y1(nmax), u1(nmax), v1(nmax)
c
      real ua, va, gg(6)
      integer igtypg, j, j0, j1, j2, nxyr
c
      data igtypg/2/
      data gg/1.,1.,1.,1.,0.,0./
c
      if (igtypa.eq.0) return
c
c  preparations:
c
      if (icall.le.2) then
c
	do j1 = 1,npos,nmax
c
	j0 = j1 - 1
	j2 = min(j0+nmax,npos)
	nxyr = j2 - j0
c
	do j = 1,nxyr
           x1(j) = geopos(2,j0+j)
           y1(j) = geopos(1,j0+j)
	enddo
c
	do j = 1,nxyr
           u1(j) = 1.
           v1(j) = 0.
	enddo
	call uvconvert(nxyr,x1(1),y1(1),u1(1),v1(1),
     +                 igtypa,ga,igtypg,gg,udef,ierror)
	if(ierror.ne.0) then
	  write(6,*) '**GRV2POS** UVCONVERT ERROR: ',ierror
	  return
	end if
	do j = 1,nxyr
           vturn(1,j0+j) = u1(j)
           vturn(3,j0+j) = v1(j)
	enddo
c
	do j = 1,nxyr
           u1(j) = 0.
           v1(j) = 1.
	enddo
	call uvconvert(nxyr,x1(1),y1(1),u1(1),v1(1),
     +                 igtypa,ga,igtypg,gg,udef,ierror)
	if(ierror.ne.0) then
	  write(6,*) '**GRV2POS** UVCONVERT ERROR: ',ierror
	  return
	end if
	do j = 1,nxyr
           vturn(2,j0+j) = u1(j)
           vturn(4,j0+j) = v1(j)
	enddo
c
	enddo
c
      end if
c
c  turning
c
      if (icall.ge.2) then
c
      do j = 1,npos
         if (ur(j).ne.udef .and. vr(j).ne.udef) then
	    ua = ur(j)
	    va = vr(j)
            ur(j) = vturn(1,j)*ua + vturn(2,j)*va
            vr(j) = vturn(3,j)*ua + vturn(4,j)*va
         else
            ur(j) = udef
            vr(j) = udef
         endif
      enddo
c
      endif
c
      return
      end

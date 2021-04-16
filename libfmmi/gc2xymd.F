      subroutine gc2xymd(icall,ndim,npos,xydata,igtype,grid,ierror)
c
c  NAME:
c     gc2xymd
c
c  PURPOSE:
c     Convertion between geographic coordinates and misc. grid coordinates.
c     Positions in multidimensioned array.
c     Calling xyconvert.
c
c  SYNOPSIS:
c     subroutine gc2xymd(icall,ndim,npos,xydata,igtype,grid,ierror)
c     integer icall,ndim,npos,igtype,ierror
c     real    xydata(ndim,npos),grid(6)
c
c  INPUT:
c     icall       - +1 : from geographic coordinates to grid coordinates
c		    -1 : from grid coordinates to geographic coordinates
c     ndim        - first diemnsion of xydata
c     npos        - no. of positions
c     igtype      - grid type: 1 = polarstereographic (true at 60N)
c		  	       2 = geographic
c			       3 = spherical (rotated)
c			       4 = polarstereographic
c			       5 = mercator (unrotated)
c     grid        - grid specifications
c
c  INPUT/OUTPUT:
c     xydata(1,n) - geographic latitude  (decimal, +/- = north/south)
c     		    or x position
c     xydata(2,n) - geographic longitude (decimal, +/- =  east/west)
c     		    or y position
c
c  OUTPUT:
c     ierror      - exit status, 0=no error
c
c  NOTES:
c     For spherical (rotated) grid (igtype=2,3):
c       grid(1) - western boundary (degrees)
c       grid(2) - southern boundary (degrees)
c       grid(3) - longitude increment (degrees)
c       grid(4) - latitude increment (degrees)
c       grid(5) - xcen: longitude position of rotated equator (degrees)
c       grid(6) - ycen: latitude  position of rotated equator (degrees)
c                 (lamda,theta)=(xcen,ycen) at (lamda',theta')=(0,0),
c                 where (lamda,theta) are usual spherical coord. and
c                 (lamda',theta') are rotated spherical coord.
c                 xcen = ycen = 0 for usual spherical coord.
c     For polar stereographic grid (igtype=1,4):
c       grid(1) - x-position of north pole
c       grid(2) - y-position of north pole
c       grid(3) - number of grid units between pole and equator
c       grid(4) - rotation angle of the grid (degrees)
c       grid(5) - projection latitude (degrees)
c	          (usually 60 degrees north)
c       grid(6) - 0. (not used)
c     For polar mercator (unrotated) grid (igtype=5):
c       grid(1) - western boundary (longitude for x=1) (degrees)
c       grid(2) - southern boundary (latitude for y=1) (degrees)
c       grid(3) - x (longitude) increment (km)
c       grid(4) - y (latitude)  increment (km)
c       grid(5) - reference (construction) latitude (degrees)
c       grid(6) - 0.  (not used)
c
c----------------------------------------------------------------------
c  DNMI/FoU  15.05.1995  Anstein Foss ... calling pos2pos
c  DNMI/FoU  12.06.1996  Anstein Foss ... new call, xyconvert, mercator
c----------------------------------------------------------------------
c
      implicit none
c
      integer icall,ndim,npos,igtype,ierror
      real    xydata(ndim,npos),grid(6)
c
      integer nmax
      parameter ( nmax = 1000 )
      real    x(nmax),y(nmax)
      integer i,n0,n1,n2,nn,igeogrid
      real    geogrid(6)
c
      data igeogrid/2/
      data  geogrid/1.,1.,1.,1.,0.,0./
c
      ierror=0
c
      if(npos.lt.1) return
c
      if(ndim.lt.2) then
	ierror=1
	return
      end if
c
      if(icall.eq.+1) then
c
        do n1=1,npos,nmax
c
	  n0=n1-1
	  n2=min(n0+nmax,npos)
	  nn=n2-n0
c
	  do i=1,nn
	    x(i)=xydata(2,n0+i)
	    y(i)=xydata(1,n0+i)
	  end do
c
	  call xyconvert(nn,x,y,igeogrid,geogrid,igtype,grid,ierror)
	  if (ierror.ne.0) return
c
	  do i=1,nn
	    xydata(1,n0+i)=x(i)
	    xydata(2,n0+i)=y(i)
	  end do
c
	end do
c
      elseif(icall.eq.-1) then
c
        do n1=1,npos,nmax
c
	  n0=n1-1
	  n2=min(n0+nmax,npos)
	  nn=n2-n0
c
	  do i=1,nn
	    x(i)=xydata(1,n0+i)
	    y(i)=xydata(2,n0+i)
	  end do
c
	  call xyconvert(nn,x,y,igtype,grid,igeogrid,geogrid,ierror)
	  if (ierror.ne.0) return
c
	  do i=1,nn
c
c..check xyconvert output
	    if (x(i).le.-180.) x(i) = x(i) + 360.
	    if (x(i).gt.+180.) x(i) = x(i) - 360.
c
	    xydata(1,n0+i)=y(i)
	    xydata(2,n0+i)=x(i)
	  end do
c
	end do
c
      else
c
	ierror=1
c
      end if
c
      return
      end

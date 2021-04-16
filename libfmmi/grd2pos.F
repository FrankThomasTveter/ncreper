      subroutine grd2pos (icall,igtypa,ga,mvx,mvy,
     +                    fa,nxa,nya,npos,pos,geopos,gxypos,fr,
     +                    jinter,rinter,ndimri,
     +                    udef,inter,intopt,iundef,ierror)
c
c****************************************************************
c
c     grd2pos - interpolation in grids to selected positions
c
c  purpose:
c
c     interpolate in different grids; spherical,
c     spherical rotated, polar stereographic or mercator.
c
c     an undefined value is set in the grid points of the output
c     fields situated outside the input area.
c
c  input/output parameters:
c
c     icall  - switch for preparations/interpolation
c              = 1: preparations
c              = 2: preparations and interpolation
c              = 3: interpolation
c     igtypa - input grid type
c     ga     - input grid parameters
c     mvx    - no. of gridpoints extended at the lower x boundary
c	       of the inpu grid
c     mvy    - no. of gridpoints extended at the lower y boundary
c	       of the input grid
c     fa     - input field
c     nxa    - x-dimension of input field
c     nya    - y-dimension of input field
c     npos   - no. of positions
c     pos    - input positions,
c		    pos(1,n) : position type, 
c				  1 = x,y position in grid
c			          2 = geographic latitude,longitude
c				      (decimal)
c				 (any other value treated as 1)
c		    pos(2,n) : x or latitude  position
c		    pos(3,n) : y or longitude position
c     geopos - output positions,
c		    geopos(1,n) : geographic latitude
c		    geopos(2,n) : geographic longitude
c     gxypos - output positions,
c		    gxypos(1,n) : x position in grid
c		    gxypos(2,n) : y position in grid
c     fr     - output data array
c     udef   - value in 'undefined' points
c     jinter - preparation matrix, jinter(npos,3)
c     rinter - preparation matrix, rinter(ndimri,npos)
c     inter  - interpolation type
c              = 1 - bilinear ............ requires ndimri >= 4
c              = 2 - biquadratic ......... requires ndimri >= 9
c              = 3 - bicubic ............. requires ndimri >= 16
c              = 4 - nearest gridpoint ... requires ndimri >= 0
c              = 5 - bilinear/nearest/ ... requires ndimri >= 6
c                    nearest existing
c                    ('extrapolation' in coast areas, for WAM etc.)
c              =44 - mean of the surrounding 4x4 gridpoints, ndimri >= 0
c     intopt - for interpolation type 3:
c                 check of values after interpolation,
c                 if the 4 nearest gridpoints all are zero then the
c                 the interpolated value will be reset to zero,
c                 intopt = 0 - check disabled
c                 intopt = 1 - check enabled
c              for interpolation type 5:
c                 intopt = 0 - bilinear interpolation / nearest existing
c                 intopt = 1 - nearest gridpoint / nearest existing
c                 (1 typically used for wave peak directions etc.)
c              not used for other interpolation types
c     iundef - switch for checking undefined values
c	       = 0 - not check undefined values
c	       = 1 - check undefined values when icall=2   (first field)
c	       = 2 - check undefined values when icall=2,3 (each  field)
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
c  DNMI/FoU  xx.11.1994  J.E. Haugen .... grd2grd
c  DNMI/FoU  02.02.1995  Anstein Foss ... no size limits, interp speedup
c  DNMI/FoU  06.02.1995  Anstein Foss ... grd2pos
c  DNMI/FoU  12.06.1996  Anstein Foss ... new call, using xyconvert,...
c  DNMI/FoU  03.06.1997  Anstein Foss ... inter=44 (4x4 mean)
c  DNMI/FoU  22.09.1998  Anstein Foss ... inter=5 for WAM etc.
c****************************************************************
c
      implicit none
c
      integer icall, igtypa, mvx, mvy, nxa, nya, npos,
     +        ndimri, inter, intopt, iundef, ierror
      integer jinter(npos,3)
      real fa(nxa*nya), pos(3,npos), geopos(2,npos), gxypos(2,npos),
     +     fr(npos), rinter(ndimri,npos),
     +     ga(6), udef
c
      integer nmax
c
      parameter ( nmax = 1000 )
c
      real x1(nmax), y1(nmax)
c
      integer igtypg, ipos, j, n, iudef, j0, j1, j2, ja, nxyr, ix, iy,
     +        imin, imax, jmin, jmax,
     +        im2m2, im1m2, im0m2, ip1m2, im2m1, im1m1, im0m1, ip1m1,
     +        im2m0, im1m0, ip1m0, im2p1, im1p1, im0p1, ip1p1
      real znx, zny, z100, zr4, zr12, zr36, za1, za2, za3, za4, zbb,
     +     zxmin, zxmax, zymin, zymax, zimin, zimax, zjmin, zjmax,
     +     zmvx, zmvy, zx, zy
      real zeps, zfix(5)
      real gg(6)
c
      data zeps /1.e-2/
      data zfix /100.0,100.5,100.0,100.5,100.0/
c
      data igtypg/2/
      data gg/1.,1.,1.,1.,0.,0./
c
      if ((inter.lt.1 .or. inter.gt.5) .and. inter.ne.44) then
         write(*,*) '**GRD2POS** invalid interpolation. inter = ',inter
         ierror=1
	 return
      endif
c
      if ((inter.eq.1 .and. ndimri.lt. 4) .or.
     +    (inter.eq.2 .and. ndimri.lt. 9) .or.
     +    (inter.eq.3 .and. ndimri.lt.16) .or.
     +    (inter.eq.5 .and. ndimri.lt. 6)) then
         write(*,*) '**GRD2POS** invalid ndimri. inter,ndimri: ',
     +						 inter,ndimri
         ierror=1
	 return
      endif
c
c  preparations:
c
      if (icall.le.2) then
c
        zmvx = float(mvx)
        zmvy = float(mvy)
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
	    write(6,*) '**GRD2POS** XYCONVERT ERROR: ',ierror
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
	    write(6,*) '**GRD2POS** XYCONVERT ERROR: ',ierror
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
      if (inter.eq.1 .or. inter.eq.4 .or. inter.eq.5) then
	 imin = 1
	 imax = nxa
	 jmin = 1
	 jmax = nya
      else
	 imin = 2
	 imax = nxa-1
	 jmin = 2
	 jmax = nya-1
      end if
c
c  avoid problems due to inaccuracy at boundaries
      zimin = float(imin) - zeps
      zimax = float(imax) + zeps
      zjmin = float(jmin) - zeps
      zjmax = float(jmax) + zeps
c
c  avoid addressing problems at boundaries in interpolation
      zxmin = float(imin) + zeps
      zxmax = float(imax) - zeps
      zymin = float(jmin) + zeps
      zymax = float(jmax) - zeps
c
      znx  = float(nxa)
      zny  = float(nya)
      if (inter.eq.44) then
         z100 = zfix(3)
      else
         z100 = zfix(inter)
      end if
c
      if (inter.eq.5) then
         do j = 1,npos
            rinter(5,j)=gxypos(1,j) + zmvx
            rinter(6,j)=gxypos(2,j) + zmvy
         enddo
      end if
c
      do j1 = 1,npos,nmax
c
      j0 = j1 - 1
      j2 = min(j0+nmax,npos)
      nxyr = j2 - j0
c
      do j = 1,nxyr
	 x1(j) = gxypos(1,j0+j) + zmvx
	 y1(j) = gxypos(2,j0+j) + zmvy
         jinter(j0+j,2) = 1
         if (x1(j).lt.zimin .or.
     +       x1(j).gt.zimax .or.
     +       y1(j).lt.zjmin .or.
     +       y1(j).gt.zjmax) jinter(j0+j,2) = 0
         x1(j) = max(x1(j),zxmin)
         x1(j) = min(x1(j),zxmax)
         y1(j) = max(y1(j),zymin)
         y1(j) = min(y1(j),zymax)
      enddo
c
c  preparation for interpolation.
c
c  inter = 1 : bilinear
c  inter = 2 : biquadratic
c  inter = 3 : bicubic
c  inter = 4 : nearest gridpoint
c  inter = 5 : bilinear / nearest gridpoint / nearest existing gridpoint
c  inter =44 : mean of nearest 4x4 gridpoints
c
      do j = 1,nxyr
         x1(j) = znx - x1(j)
         y1(j) = zny - y1(j)
         ix    = int (x1(j) + z100) - 100
         iy    = int (y1(j) + z100) - 100
         x1(j) = x1(j) - float(ix)
         y1(j) = y1(j) - float(iy)
         ix    = nxa - ix
         iy    = nya - iy
         jinter(j0+j,1) = (iy - 1)*nxa + ix
      enddo
c
      if (inter.eq.1 .or. inter.eq.5) then
c
c  constants for  bilinear interpolation
c
      do j = 1,nxyr
         za1            = 1. - x1(j)
         zbb            = 1. - y1(j)
         rinter(1,j0+j) = zbb * za1
         rinter(2,j0+j) = zbb * x1(j)
         rinter(3,j0+j) = y1(j) * za1
         rinter(4,j0+j) = y1(j) * x1(j)
      enddo
c
      elseif (inter.eq.2) then
c
c  constants for biquadratic interpolation
c
      do j = 1,nxyr
         za1            = x1(j) * (1. + x1(j))
         za2            = 1. - x1(j) * x1(j)
         za3            = x1(j) * (1. - x1(j))
         zbb            = y1(j) * (1. + y1(j))
         rinter(1,j0+j) =   zbb * 0.25 * za1
         rinter(2,j0+j) =   zbb * 0.5  * za2
         rinter(3,j0+j) = - zbb * 0.25 * za3
         zbb            = (1. - y1(j) * y1(j))
         rinter(4,j0+j) =   zbb * 0.5  * za1
         rinter(5,j0+j) =   zbb *        za2
         rinter(6,j0+j) = - zbb * 0.5  * za3
         zbb            = y1(j) * (1. - y1(j))
         rinter(7,j0+j) = - zbb * 0.25 * za1
         rinter(8,j0+j) = - zbb * 0.5  * za2
         rinter(9,j0+j) =   zbb * 0.25 * za3
      enddo
c
      elseif (inter.eq.3) then
c
c  constants for bicubic interpolation
c
      zr4   = 1./4.
      zr12  = 1./12.
      zr36  = 1./36.
c
      do j = 1,nxyr
         za1             = x1(j) * (1. - x1(j) * x1(j))
         za2             = x1(j) * (1. + x1(j)) * (2. - x1(j))
         za3             = (1. - x1(j) * x1(j)) * (2. - x1(j))
         za4             = x1(j) * (1. - x1(j)) * (2. - x1(j))
         zbb             = y1(j) * (1. - y1(j) * y1(j))
         rinter( 1,j0+j) =   zbb * zr36 * za1
         rinter( 2,j0+j) = - zbb * zr12 * za2
         rinter( 3,j0+j) = - zbb * zr12 * za3
         rinter( 4,j0+j) =   zbb * zr36 * za4
         zbb             = y1(j) * (1. + y1(j)) * (2. - y1(j))
         rinter( 5,j0+j) = - zbb * zr12 * za1
         rinter( 6,j0+j) =   zbb * zr4  * za2
         rinter( 7,j0+j) =   zbb * zr4  * za3
         rinter( 8,j0+j) = - zbb * zr12 * za4
         zbb             = (1. - y1(j) * y1(j)) * (2. - y1(j))
         rinter( 9,j0+j) = - zbb * zr12 * za1
         rinter(10,j0+j) =   zbb * zr4  * za2
         rinter(11,j0+j) =   zbb * zr4  * za3
         rinter(12,j0+j) = - zbb * zr12 * za4
         zbb             = y1(j) * (1. - y1(j)) * (2. - y1(j))
         rinter(13,j0+j) =   zbb * zr36 * za1
         rinter(14,j0+j) = - zbb * zr12 * za2
         rinter(15,j0+j) = - zbb * zr12 * za3
         rinter(16,j0+j) =   zbb * zr36 * za4
      enddo
c
      end if
c
      enddo
c
      do j = 1,npos
	 jinter(j,3) = jinter(j,2)
      enddo
c
      endif
c
c  check undefined values
c
      iudef = 0
      if (icall.eq.2 .and. iundef.eq.1) iudef = 1
      if (icall.eq.2 .and. iundef.eq.2) iudef = 1
      if (icall.eq.3 .and. iundef.eq.2) iudef = 1
      if (inter.eq.4) iudef = 0
      if (icall.eq.3 .and. inter.eq.5 .and. intopt.eq.1) iudef = 0
c  
      if (iudef.eq.1) then
c
      iudef = 0
      do j = 1,nxa*nya
         if (fa(j).eq.udef) iudef = 1
      enddo
c
      if (iudef.eq.1) then
c
c  set jinter(i,3)=0 in point with udef value in some argument points.
c
      do j = 1,npos
	 jinter(j,3) = jinter(j,2)
      enddo
c
      if (inter.eq.1 .or. inter.eq.5) then
c
c  bilinear interpolation
c
      im0m1 = - nxa
      im1m1 = - nxa - 1
c
      do j = 1,npos
         ix  = jinter(j,1)
         if(fa(ix      ).eq.udef .or.
     +      fa(ix-1    ).eq.udef .or.
     +      fa(ix+im0m1).eq.udef .or.
     +      fa(ix+im1m1).eq.udef) jinter(j,3) = 0
      enddo
c
      elseif (inter.eq.2) then
c
c  biquadratic interpolation
c
      im1m1 = - nxa - 1
      im0m1 = - nxa
      ip1m1 = - nxa + 1
      im1m0 = - 1
      ip1m0 = 1
      im1p1 = nxa - 1
      im0p1 = nxa
      ip1p1 = nxa + 1
c
      do j = 1,npos
         ix  = jinter(j,1)
         if (fa(ix+im1m1).eq.udef .or.
     +       fa(ix+im0m1).eq.udef .or.
     +       fa(ix+ip1m1).eq.udef .or.
     +       fa(ix+im1m0).eq.udef .or.
     +       fa(ix      ).eq.udef .or.
     +       fa(ix+ip1m0).eq.udef .or.
     +       fa(ix+im1p1).eq.udef .or.
     +       fa(ix+im0p1).eq.udef .or.
     +       fa(ix+ip1p1).eq.udef) jinter(j,3) = 0
      enddo
c
      elseif (inter.eq.3 .or. inter.eq.44) then
c
c  bicubic interpolation  or  4x4 mean
c
      im2m2 = - nxa - nxa - 2
      im1m2 = - nxa - nxa - 1
      im0m2 = - nxa - nxa
      ip1m2 = - nxa - nxa + 1
      im2m1 = - nxa - 2
      im1m1 = - nxa - 1
      im0m1 = - nxa
      ip1m1 = - nxa + 1
      im2m0 = - 2
      im1m0 = - 1
      ip1m0 = 1
      im2p1 = nxa - 2
      im1p1 = nxa - 1
      im0p1 = nxa
      ip1p1 = nxa + 1
c
      do j = 1,npos
         ix  = jinter(j,1)
         if (fa(ix+im2m2).eq.udef .or.
     +       fa(ix+im1m2).eq.udef .or.
     +       fa(ix+im0m2).eq.udef .or.
     +       fa(ix+ip1m2).eq.udef .or.
     +       fa(ix+im2m1).eq.udef .or.
     +       fa(ix+im1m1).eq.udef .or.
     +       fa(ix+im0m1).eq.udef .or.
     +       fa(ix+ip1m1).eq.udef .or.
     +       fa(ix+im2m0).eq.udef .or.
     +       fa(ix+im1m0).eq.udef .or.
     +       fa(ix      ).eq.udef .or.
     +       fa(ix+ip1m0).eq.udef .or.
     +       fa(ix+im2p1).eq.udef .or.
     +       fa(ix+im1p1).eq.udef .or.
     +       fa(ix+im0p1).eq.udef .or.
     +       fa(ix+ip1p1).eq.udef) jinter(j,3) = 0
      enddo
c
      end if
c
      elseif (icall.eq.3) then
	 do j = 1,npos
	    jinter(j,3) = jinter(j,2)
	 enddo
      endif
c
      end if
c
      if (icall.ge.2) then
c
c  interpolation
c
c  interpolation from values in fa to fr.
c  preparation for this is done above.
c
c  inter = 1 : bilinear
c  inter = 2 : biquadratic
c  inter = 3 : bicubic
c  inter = 4 : nearest gridpoint
c  inter = 5 : bilinear / nearest gridpoint / nearest existing gridpoint
c  inter =44 : 4x4 mean
c  rinter    : array prepared above
c  jinter    : array prepared above
c  intopt    : for interpolation type 3:
c                 check of values after interpolation,
c                 if the 4 nearest gridpoints all are zero then the
c                 the interpolated value will be reset to zero,
c                 intopt = 0 - check disabled
c                 intopt = 1 - check enabled
c              for interpolation type 5:
c                 intopt = 0 - bilinear interpolation / nearest existing
c                 intopt = 1 - nearest gridpoint / nearest existing
c                 (1 typically used for wave peak directions etc.)
c              not used for other interpolation types
c
      if (inter.eq.1 .or. (inter.eq.5 .and. intopt.eq.0)) then
c
c  bilinear interpolation
c
      im0m1 = - nxa
      im1m1 = - nxa - 1
c
      do j = 1,npos
         ix = jinter(j,1)
         fr(j) =   rinter(1,j) * fa(ix      )
     +           + rinter(2,j) * fa(ix-1    )
     +           + rinter(3,j) * fa(ix+im0m1)
     +           + rinter(4,j) * fa(ix+im1m1)
      enddo
c
      elseif (inter.eq.2) then
c
c  biquadratic interpolation
c
      im1m1 = - nxa - 1
      im0m1 = - nxa
      ip1m1 = - nxa + 1
      im1m0 = - 1
      ip1m0 = 1
      im1p1 = nxa - 1
      im0p1 = nxa
      ip1p1 = nxa + 1
c
      do j = 1,npos
         ix = jinter(j,1)
         fr(j) =   rinter(1,j) * fa(ix+im1m1)
     +           + rinter(2,j) * fa(ix+im0m1)
     +           + rinter(3,j) * fa(ix+ip1m1)
     +           + rinter(4,j) * fa(ix+im1m0)
     +           + rinter(5,j) * fa(ix      )
     +           + rinter(6,j) * fa(ix+ip1m0)
     +           + rinter(7,j) * fa(ix+im1p1)
     +           + rinter(8,j) * fa(ix+im0p1)
     +           + rinter(9,j) * fa(ix+ip1p1)
      enddo
c
      elseif (inter.eq.3) then
c
c  bicubic interpolation
c
      im2m2 = - nxa - nxa - 2
      im1m2 = - nxa - nxa - 1
      im0m2 = - nxa - nxa
      ip1m2 = - nxa - nxa + 1
      im2m1 = - nxa - 2
      im1m1 = - nxa - 1
      im0m1 = - nxa
      ip1m1 = - nxa + 1
      im2m0 = - 2
      im1m0 = - 1
      ip1m0 = 1
      im2p1 = nxa - 2
      im1p1 = nxa - 1
      im0p1 = nxa
      ip1p1 = nxa + 1
c
      do j = 1,npos
         ix = jinter(j,1)
         fr(j) =   rinter( 1,j) * fa(ix+im2m2)
     +           + rinter( 2,j) * fa(ix+im1m2)
     +           + rinter( 3,j) * fa(ix+im0m2)
     +           + rinter( 4,j) * fa(ix+ip1m2)
     +           + rinter( 5,j) * fa(ix+im2m1)
     +           + rinter( 6,j) * fa(ix+im1m1)
     +           + rinter( 7,j) * fa(ix+im0m1)
     +           + rinter( 8,j) * fa(ix+ip1m1)
     +           + rinter( 9,j) * fa(ix+im2m0)
     +           + rinter(10,j) * fa(ix+im1m0)
     +           + rinter(11,j) * fa(ix      )
     +           + rinter(12,j) * fa(ix+ip1m0)
     +           + rinter(13,j) * fa(ix+im2p1)
     +           + rinter(14,j) * fa(ix+im1p1)
     +           + rinter(15,j) * fa(ix+im0p1)
     +           + rinter(16,j) * fa(ix+ip1p1)
      enddo
c
      elseif (inter.eq.4) then
c
c  nearest gridpoint
c
      do j = 1,npos
         ix = jinter(j,1)
         fr(j) = fa(ix)
      enddo
c
      elseif (inter.eq.5 .and. intopt.ne.0) then
c
c  nearest gridpoint
c
      do j = 1,npos
         ix = nint(rinter(5,j))
         iy = nint(rinter(6,j))
         if (ix.ge.1 .and. ix.le.nxa .and.
     +       iy.ge.1 .and. iy.le.nya) then
            ix = ix + (iy-1)*nxa
            fr(j) = fa(ix)
         else
            fr(j) = udef
         end if
      enddo
c
      elseif (inter.eq.44) then
c
c  4x4 mean
c
      im2m2 = - nxa - nxa - 2
      im1m2 = - nxa - nxa - 1
      im0m2 = - nxa - nxa
      ip1m2 = - nxa - nxa + 1
      im2m1 = - nxa - 2
      im1m1 = - nxa - 1
      im0m1 = - nxa
      ip1m1 = - nxa + 1
      im2m0 = - 2
      im1m0 = - 1
      ip1m0 = 1
      im2p1 = nxa - 2
      im1p1 = nxa - 1
      im0p1 = nxa
      ip1p1 = nxa + 1
c
      do j = 1,npos
         ix = jinter(j,1)
         fr(j) = (  fa(ix+im2m2) + fa(ix+im1m2)
     +            + fa(ix+im0m2) + fa(ix+ip1m2)
     +            + fa(ix+im2m1) + fa(ix+im1m1)
     +            + fa(ix+im0m1) + fa(ix+ip1m1)
     +            + fa(ix+im2m0) + fa(ix+im1m0)
     +            + fa(ix      ) + fa(ix+ip1m0)
     +            + fa(ix+im2p1) + fa(ix+im1p1)
     +            + fa(ix+im0p1) + fa(ix+ip1p1)) / 16.
      enddo
c
      end if
c
      if (inter.eq.3 .and. intopt.eq.1) then
c
         im1m1 = - nxa - 1
         im0m1 = - nxa
         im1m0 = - 1
c
         do j = 1,npos
            ix = jinter(j,1)
	    if (fa(ix+im1m1).eq.0. .and.
     +	        fa(ix+im0m1).eq.0. .and.
     +	        fa(ix+im1m0).eq.0. .and.
     +	        fa(ix      ).eq.0.) fr(j) = 0.
         enddo
c
      end if
c
      if (inter.ne.5 .or. intopt.eq.0) then
c
c  insert undefined value udef in points outside data area
c  (also when input field contains undefined values)
c
         do j = 1,npos
            if (jinter(j,3).eq.0) fr(j) = udef
         enddo
c
      endif
c
      if (inter.eq.5) then
c
c  if value is undefined after interpolation above,
c  but a 'real' value is requested, the nearest existing point
c  is used.  Only the nearest 5x5 input grid points are searched!
c
         do j = 1,npos
            if (fr(j).eq.udef) then
               zx = rinter(5,j)
               zy = rinter(6,j)
               ix = nint(zx)
               iy = nint(zy)
               imin= max(ix-2,1)
               imax= min(ix+2,nxa)
               jmin= max(iy-2,1)
               jmax= min(iy+2,nya)
               za1= udef
               zbb= udef
               do iy=jmin,jmax
                  do ix=imin,imax
                     ja=(iy-1)*nxa+ix
                     if(fa(ja).ne.udef) then
                        za2=(zx-float(ix))**2 + (zy-float(iy))**2
                        if(za2.lt.za1) then
                           za1 = za2
                           zbb = fa(ja)
                        end if
                     end if
                  end do
               end do
               fr(j) = zbb
            end if
         enddo
c
      endif
c
      endif
c
      return
      end

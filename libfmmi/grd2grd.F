      subroutine grd2grd (icall,igtypa,ga,igtypr,gr,
     +                    fa,nxa,nya,fr,nxr,nyr,
     +                    jinter,rinter,ndimri,
     +                    udef,inter,intopt,iundef,ierror)
c
c****************************************************************
c
c     grd2grd - interpolation between grids
c
c  purpose:
c
c     interpolate fields between different grids; spherical,
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
c     igtypr - output grid type
c     gr     - output grid parameters
c     fa     - input field
c     nxa    - x-dimension of input field
c     nya    - y-dimension of input field
c     fr     - output field
c     nxr    - x-dimension of output field
c     nyr    - y-dimension of output field
c     udef   - value in 'undefined' points
c     jinter - preparation matrix, jinter(nxr*nyr,3) for inter=1,2,3,4
c                                  jinter(nxr,nyr,4) for inter=5
c              for inter=5 jinter(i,4) must be set before call:
c                 jinter(i,4)=0 : undefined output point
c                 jinter(i,4)=1 : point to be interpolated,
c                                 possibly nearest existing point
c     rinter - preparation matrix, rinter(ndimri,nxr*nyr)
c     inter  - interpolation type
c              = 1 - bilinear ............ requires ndimri >= 4
c              = 2 - biquadratic ......... requires ndimri >= 9
c              = 3 - bicubic ............. requires ndimri >= 16
c              = 4 - nearest gridpoint ... requires ndimri >= 0
c              = 5 - bilinear/nearest/ ... requires ndimri >= 6
c                    nearest existing
c                    ('extrapolation' in coast areas, for WAM etc.)
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
c  description of g = ga and gr (for igtype = igtypa and igtypr):
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
c     g(3) - number of grid distances between pole and equator
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
c  externals:
c
c     xyconvert - performs all coordinate conversions
c
c  history:
c
c  DNMI/FoU  xx.11.1994  J.E. Haugen
c  DNMI/FoU  02.02.1995  Anstein Foss ... no size limits, interp speedup
c  DNMI/FoU  04.06.1996  Anstein Foss ... new call, using xyconvert,...
c  DNMI/FoU  21.09.1998  Anstein Foss ... inter=5 for WAM etc.
c  DNMI/FoU  23.02.2000  Anstein Foss ... geographic longitude wrapping
c****************************************************************
c
      implicit none
c
      integer icall, igtypa, igtypr, nxa, nya, nxr, nyr,
     +        ndimri, inter, intopt, iundef, ierror
      integer jinter(nxr*nyr,4)
      real fa(nxa*nya), fr(nxr*nyr), rinter(ndimri,nxr*nyr),
     +     ga(6), gr(6), udef
c
      integer nmax
c
      parameter ( nmax = 1000 )
c
      real x1(nmax), y1(nmax)
c
      integer j, n, iudef, j0, j1, j2, ja, nxyr, ix, iy,
     +        imin, imax, jmin, jmax,
     +        im2m2, im1m2, im0m2, ip1m2, im2m1, im1m1, im0m1, ip1m1,
     +        im2m0, im1m0, ip1m0, im2p1, im1p1, im0p1, ip1p1
      real znx, zny, z100, zr4, zr12, zr36, za1, za2, za3, za4, zbb,
     +     zxmin, zxmax, zymin, zymax, zimin, zimax, zjmin, zjmax,
     +     zx, zy
      real zwrap, zeps, zfix(5)
c
      data zeps /1.e-2/
      data zfix /100.0,100.5,100.0,100.5,100.0/
c
c  Implementation notes:
c  - for inter=5 the rinter and jinter arrays are prepared for
c    bilinear interpolation (intopt=0), in addition rinter(5:6,n)
c    is used for nearest point and nearest existing point.
c
      if (inter.lt.1 .or. inter.gt.5) then
         write(*,*) '**GRD2GRD** invalid interpolation. inter = ',inter
         ierror=1
	 return
      endif
c
      if ((inter.eq.1 .and. ndimri.lt. 4) .or.
     +    (inter.eq.2 .and. ndimri.lt. 9) .or.
     +    (inter.eq.3 .and. ndimri.lt.16) .or.
     +    (inter.eq.5 .and. ndimri.lt. 6)) then
         write(*,*) '**GRD2GRD** invalid ndimri. inter,ndimri: ',
     +						 inter,ndimri
         ierror=1
	 return
      endif
c
c  preparations:
c
      if (icall.le.2) then
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
        z100 = zfix(inter)
c
	do j1 = 1,nxr*nyr,nmax
c
	j0 = j1 - 1
	j2 = min(j0+nmax,nxr*nyr)
	nxyr = j2 - j0
c
	do j = 1,nxyr
	   iy = (j0+j+nxr-1)/nxr
	   ix = j0+j - (iy-1)*nxr
           x1(j) = float(ix)
           y1(j) = float(iy)
	enddo
c
	call xyconvert(nxyr,x1(1),y1(1),igtypr,gr,igtypa,ga,ierror)
	if(ierror.ne.0) then
	  write(6,*) '**GRD2GRD** XYCONVERT ERROR: ',ierror
	  return
	end if
c
c..xyconvert does not know about the fields dimensions,
c..avoid some of the "longitude" problems for input geographic grid...
	if((igtypa.eq.2 .or. igtypa.eq.3) .and. ga(5).eq.0.
     +				          .and. ga(6).eq.0.) then
	   zwrap=360./ga(3)
	   do j = 1,nxyr
              if(x1(j).lt.zxmin .and. x1(j)+zwrap.lt.zxmax) then
                x1(j)=x1(j)+zwrap
              elseif(x1(j).gt.zxmax .and. x1(j)-zwrap.gt.zxmin) then
                x1(j)=x1(j)-zwrap
	      end if
	   enddo
	end if
c
        if (inter.eq.5) then
           do j = 1,nxyr
              rinter(5,j0+j)=x1(j)
              rinter(6,j0+j)=y1(j)
           enddo
        end if
c
      do j = 1,nxyr
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
      do j = 1,nxr*nyr
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
      n = nxr*nyr
c
      do j = 1,n
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
      do j = 1,n
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
      do j = 1,n
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
      do j = 1,n
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
	 do j = 1,nxr*nyr
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
      n = nxr*nyr
c
      if (inter.eq.1 .or. (inter.eq.5 .and. intopt.eq.0)) then
c
c  bilinear interpolation
c
      im0m1 = - nxa
      im1m1 = - nxa - 1
c
      do j = 1,n
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
      do j = 1,n
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
      do j = 1,n
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
      do j = 1,n
         ix = jinter(j,1)
         fr(j) = fa(ix)
      enddo
c
      elseif (inter.eq.5 .and. intopt.ne.0) then
c
c  nearest gridpoint
c
      do j = 1,n
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
      end if
c
      if (inter.eq.3 .and. intopt.eq.1) then
c
         im1m1 = - nxa - 1
         im0m1 = - nxa
         im1m0 = - 1
c
         do j = 1,n
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
c  or when input field contains undefined values
c
         do j = 1,nxr*nyr
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
         do j = 1,nxr*nyr
            if (jinter(j,4).ne.0 .and. fr(j).eq.udef) then
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
            elseif (jinter(j,4).eq.0 .and. fr(j).ne.udef) then
               fr(j) = udef
            end if
         enddo
c
      endif
c
      endif
c
      return
      end

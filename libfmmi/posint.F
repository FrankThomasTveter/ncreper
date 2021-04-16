      subroutine posint (arg,nx,ny,res,gxypos,npos,inter,ipzero)
c
c  interpolation from values in arg to res.
c  simplified version of interpolation in grd2pos,
c  without preparation arrays.
c  all positions must be inside the grid area.
c
c  inter = 1 : bilinear
c  inter = 2 : biquadratic
c  inter = 3 : bicubic (bessel)
c  inter = 4 : nearest gridpoint
c  ipzero    : check of values after interpolation type 3,
c	       if the 4 nearest gridpoints all are zero then the
c	       the interpolated value will be reset to zero,
c	       = 0 check disabled
c	       = 1 check enabled
c
c  history:
c
c     j.e. haugen/dnmi   nov -94 ...... calcpq,chkpos,interp
c     a.   foss   dnmi   25.04.1995 ... posint
c     a.   foss   dnmi   18.06.1996
c
      integer nx, ny, npos, inter, ipzero
      real    arg(nx*ny), res(npos), gxypos(2,npos)
c
      real xfix(4)
      data xfix /100.0,100.5,100.0,100.5/
c
      if (inter.lt.1 .or. inter.gt.4) then
         write(*,*) 'posint: invalid interpolation , inter = ',inter
         stop
      endif
c
      znx  = float(nx)
      zny  = float(ny)
      z100 = xfix(inter)
c
      if (inter.eq.1) then
c
c  bilinear interpolation
c
      im0m1 = - nx
      im1m1 = - nx - 1
c
      do j = 1,npos
         alfh = znx - gxypos(1,j)
         beth = zny - gxypos(2,j)
         kx   = ifix (alfh + z100) - 100
         ky   = ifix (beth + z100) - 100
         alfh = alfh - float(kx)
         beth = beth - float(ky)
         kx   = nx - kx
         ky   = ny - ky
         ix   = (ky - 1)*nx + kx
         za1  = 1. - alfh
         zbb  = 1. - beth
         zw01 = zbb * za1
         zw02 = zbb * alfh
         zw03 = beth * za1
         zw04 = beth * alfh
c
         res(j) =   zw01 * arg(ix      )
     +            + zw02 * arg(ix-1    )
     +            + zw03 * arg(ix+im0m1)
     +            + zw04 * arg(ix+im1m1)
      enddo
c
      elseif (inter.eq.2) then
c
c  biquadratic interpolation
c
      im1m1 = - nx - 1
      im0m1 = - nx
      ip1m1 = - nx + 1
      im1m0 = - 1
      ip1m0 = 1
      im1p1 = nx - 1
      im0p1 = nx
      ip1p1 = nx + 1
c
      do j = 1,npos
         alfh = znx - gxypos(1,j)
         beth = zny - gxypos(2,j)
         kx   = ifix (alfh + z100) - 100
         ky   = ifix (beth + z100) - 100
         alfh = alfh - float(kx)
         beth = beth - float(ky)
         kx   = nx - kx
         ky   = ny - ky
         ix   = (ky - 1)*nx + kx
         za1  = alfh * (1. + alfh)
         za2  = 1. - alfh * alfh
         za3  = alfh * (1. - alfh)
         zbb  = beth * (1. + beth)
         zw01 =   zbb * 0.25 * za1
         zw02 =   zbb * 0.5  * za2
         zw03 = - zbb * 0.25 * za3
         zbb  = (1. - beth * beth)
         zw04 =   zbb * 0.5  * za1
         zw05 =   zbb *        za2
         zw06 = - zbb * 0.5  * za3
         zbb  = beth * (1. - beth)
         zw07 = - zbb * 0.25 * za1
         zw08 = - zbb * 0.5  * za2
         zw09 =   zbb * 0.25 * za3
c
         res(j) =   zw01 * arg(ix+im1m1)
     +            + zw02 * arg(ix+im0m1)
     +            + zw03 * arg(ix+ip1m1)
     +            + zw04 * arg(ix+im1m0)
     +            + zw05 * arg(ix      )
     +            + zw06 * arg(ix+ip1m0)
     +            + zw07 * arg(ix+im1p1)
     +            + zw08 * arg(ix+im0p1)
     +            + zw09 * arg(ix+ip1p1)
      enddo
c
      elseif (inter.eq.3) then
c
c  bicubic interpolation
c
      zr4   = 1./4.
      zr12  = 1./12.
      zr36  = 1./36.
c
      im2m2 = - nx - nx - 2
      im1m2 = - nx - nx - 1
      im0m2 = - nx - nx
      ip1m2 = - nx - nx + 1
      im2m1 = - nx - 2
      im1m1 = - nx - 1
      im0m1 = - nx
      ip1m1 = - nx + 1
      im2m0 = - 2
      im1m0 = - 1
      ip1m0 = 1
      im2p1 = nx - 2
      im1p1 = nx - 1
      im0p1 = nx
      ip1p1 = nx + 1
c
      do j = 1,npos
         alfh = znx - gxypos(1,j)
         beth = zny - gxypos(2,j)
         kx   = ifix (alfh + z100) - 100
         ky   = ifix (beth + z100) - 100
         alfh = alfh - float(kx)
         beth = beth - float(ky)
         kx   = nx - kx
         ky   = ny - ky
         ix   = (ky - 1)*nx + kx
         za1  = alfh * (1. - alfh * alfh)
         za2  = alfh * (1. + alfh) * (2. - alfh)
         za3  = (1. - alfh * alfh) * (2. - alfh)
         za4  = alfh * (1. - alfh) * (2. - alfh)
         zbb  = beth * (1. - beth * beth)
         zw01 =   zbb * zr36 * za1
         zw02 = - zbb * zr12 * za2
         zw03 = - zbb * zr12 * za3
         zw04 =   zbb * zr36 * za4
         zbb  = beth * (1. + beth) * (2. - beth)
         zw05 = - zbb * zr12 * za1
         zw06 =   zbb * zr4  * za2
         zw07 =   zbb * zr4  * za3
         zw08 = - zbb * zr12 * za4
         zbb  = (1. - beth * beth) * (2. - beth)
         zw09 = - zbb * zr12 * za1
         zw10 =   zbb * zr4  * za2
         zw11 =   zbb * zr4  * za3
         zw12 = - zbb * zr12 * za4
         zbb  = beth * (1. - beth) * (2. - beth)
         zw13 =   zbb * zr36 * za1
         zw14 = - zbb * zr12 * za2
         zw15 = - zbb * zr12 * za3
         zw16 =   zbb * zr36 * za4
c
         res(j) =   zw01 * arg(ix+im2m2)
     +            + zw02 * arg(ix+im1m2)
     +            + zw03 * arg(ix+im0m2)
     +            + zw04 * arg(ix+ip1m2)
     +            + zw05 * arg(ix+im2m1)
     +            + zw06 * arg(ix+im1m1)
     +            + zw07 * arg(ix+im0m1)
     +            + zw08 * arg(ix+ip1m1)
     +            + zw09 * arg(ix+im2m0)
     +            + zw10 * arg(ix+im1m0)
     +            + zw11 * arg(ix      )
     +            + zw12 * arg(ix+ip1m0)
     +            + zw13 * arg(ix+im2p1)
     +            + zw14 * arg(ix+im1p1)
     +            + zw15 * arg(ix+im0p1)
     +            + zw16 * arg(ix+ip1p1)
      enddo
c
      elseif (inter.eq.4) then
c
c  nearest gridpoint
c
      do j = 1,npos
         alfh = znx - gxypos(1,j)
         beth = zny - gxypos(2,j)
         kx   = nx - ifix (alfh + z100) - 100
         ky   = ny - ifix (beth + z100) - 100
         ix = (ky - 1)*nx + kx
c
         res(j) = arg(ix)
      enddo
c
      end if
c
      if (inter.eq.3 .and. ipzero.eq.1) then
c
         im1m1 = - nx - 1
         im0m1 = - nx
         im1m0 = - 1
c
         do j = 1,npos
            alfh = znx - gxypos(1,j)
            beth = zny - gxypos(2,j)
            kx   = nx - ifix (alfh + z100) - 100
            ky   = ny - ifix (beth + z100) - 100
            ix   = (ky - 1)*nx + kx
c
	    if (arg(ix+im1m1).eq.0. .and.
     +	        arg(ix+im0m1).eq.0. .and.
     +	        arg(ix+im1m0).eq.0. .and.
     +	        arg(ix      ).eq.0.) res(j) = 0.
         enddo
c
      end if
c
      return
      end

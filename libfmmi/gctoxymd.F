      subroutine gctoxymd(ndim,npos,xydata,grid)
c
c  NAME:
c     gctoxymd
c
c  PURPOSE:
c     Convert geographic coordinates to polarstereographic coordinates.
c     Positions in multidimensioned array.
c
c  SYNOPSIS:
c     subroutine gctoxymd(ndim,npos,xydata,grid)
c     integer ndim,npos
c     real    xydata(ndim,npos),grid(4)
c
c  INPUT:
c     npos        - no. of positions
c     xydata(1,n) - geographic latitude  (decimal, +/- = north/south)
c     xydata(2,n) - geographic longitude (decimal, +/- =  east/west)
c     grid        - polarstereographic grid specifications
c                      grid(1) : xp, north pole x coordinate
c                      grid(2) : yp, north pole y coordinate
c                      grid(3) : an, no. of grid units between
c                                    north pole and equator
c                      grid(4) : fi, angle of rotation between
c                                    y axis and 0 longitude
c                                    (+/- = east/west)
c
c  OUTPUT:
c     xydata(1,n) - polarstereographic x position
c     xydata(2,n) - polarstereographic y position
c
c----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.1990  Anstein Foss
c  DNMI/FoU  13.03.1994  Anstein Foss
c-------------------------------------------------------------------
c
      integer ndim,npos
      real    xydata(ndim,npos),grid(4)
c
      if(ndim.lt.2) return
c
      xp=grid(1)
      yp=grid(2)
      an=grid(3)
      fi=grid(4)
c
      rad=3.1415927/180.
      vxr=(90.+fi)*rad
      alfa=sin(vxr)
      beta=cos(vxr)
c
      do n=1,npos
        glat=xydata(1,n)*rad
        glon=xydata(2,n)*rad
        rr=an*cos(glat)/(1.+sin(glat))
        xr=+rr*sin(glon)
        yr=-rr*cos(glon)
        xydata(1,n)=xr*alfa-yr*beta+xp
        xydata(2,n)=yr*alfa+xr*beta+yp
      end do
c
      return
      end

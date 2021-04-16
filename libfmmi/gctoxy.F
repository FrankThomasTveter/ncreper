      subroutine gctoxy(npos,x,y,grid)
c
c  NAME:
c     gctoxy
c
c  PURPOSE:
c     Convert geographic coordinates to polarstereographic coordinates.
c
c  SYNOPSIS:
c     subroutine gctoxy(npos,x,y,grid)
c     integer npos
c     real    x(npos),y(npos),grid(4)
c
c  INPUT:
c     npos   - no. of positions
c     x      - geographic latitude  (decimal, +/- = north/south)
c     y      - geographic longitude (decimal, +/- =  east/west)
c     grid   - polarstereographic grid specifications
c                    grid(1) : xp, north pole x coordinate
c                    grid(2) : yp, north pole y coordinate
c                    grid(3) : an, no. of grid units between
c                                  north pole and equator
c                    grid(4) : fi, angle of rotation between
c                                  y axis and 0 longitude
c                                  (+/- = east/west)
c
c  OUTPUT:
c     x      - polarstereographic x position
c     y      - polarstereographic y position
c
c----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.1990  Anstein Foss
c  DNMI/FoU  13.03.1994  Anstein Foss
c-------------------------------------------------------------------
c
      integer npos
      real    x(npos),y(npos),grid(4)
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
        glat=x(n)*rad
        glon=y(n)*rad
        rr=an*cos(glat)/(1.+sin(glat))
        xr=+rr*sin(glon)
        yr=-rr*cos(glon)
        x(n)=xr*alfa-yr*beta+xp
        y(n)=yr*alfa+xr*beta+yp
      end do
c
      return
      end

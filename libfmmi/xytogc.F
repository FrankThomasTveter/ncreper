      subroutine xytogc(npos,x,y,grid)
c
c  NAME:
c     xytogc
c
c  PURPOSE:
c     Convert polarstereographic coordinates to geographic coordinates.
c
c  SYNOPSIS:
c     subroutine xytogc(npos,x,y,grid)
c     integer npos
c     real    x(npos),y(npos),grid(4)
c
c  INPUT:
c     npos   - no. of positions
c     x      - polarstereographic x position
c     y      - polarstereographic y position
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
c     x      - geographic latitude  (decimal, +/- = north/south)
c     y      - geographic longitude (decimal, +/- =  east/west)
c
c----------------------------------------------------------------------
c  DNMI/FoU  30.11.1992  Anstein Foss
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
      deg=180./3.1415927
      deg2=deg*2.
c
      do n=1,npos
        dx=x(n)-xp
        dy=yp-y(n)
        rr=sqrt(dx*dx+dy*dy)
        glat=90.-deg2*atan(rr/an)
        glon=0.
        if(rr.gt.1.e-10)  glon=fi+deg*atan2(dx,dy)
        if(glon.le.-180.) glon=glon+360.
        if(glon.gt.+180.) glon=glon-360.
        x(n)=glat
        y(n)=glon
      end do
c
      return
      end

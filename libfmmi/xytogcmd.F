      subroutine xytogcmd(ndim,npos,xydata,grid)
c
c  NAME:
c     xytogcmd
c
c  PURPOSE:
c     Convert polarstereographic coordinates to geographic coordinates.
c     Positions in multidimensioned array.
c
c  SYNOPSIS:
c     subroutine xytogcmd(ndim,npos,xydata,grid)
c     integer ndim,npos
c     real    xydata(ndim,npos),grid(4)
c
c  INPUT:
c     npos        - no. of positions
c     xydata(1,n) - polarstereographic x position
c     xydata(2,n) - polarstereographic y position
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
c     xydata(1,n) - geographic latitude  (decimal, +/- = north/south)
c     xydata(2,n) - geographic longitude (decimal, +/- =  east/west)
c
c----------------------------------------------------------------------
c  DNMI/FoU  30.11.1992  Anstein Foss
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
      deg=180./3.1415927
      deg2=deg*2.
c
      do n=1,npos
        dx=xydata(1,n)-xp
        dy=yp-xydata(2,n)
        rr=sqrt(dx*dx+dy*dy)
        glat=90.-deg2*atan(rr/an)
        glon=0.
        if(rr.gt.1.e-10)  glon=fi+deg*atan2(dx,dy)
        if(glon.le.-180.) glon=glon+360.
        if(glon.gt.+180.) glon=glon-360.
        xydata(1,n)=glat
        xydata(2,n)=glon
      end do
c
      return
      end

      subroutine mer2sph(icall,n,x,y,xw,ys,yc,ierror)
c
c  conversion between mercator (unrotated) and spherical
c  coordinates. mercator coordinates are given in real
c  numbers, where (xmer,ymer)=(1.0,1.0) in lower left corner of
c  the mercator grid. spherical coordinates are given in radians.
c  xw,ys,yc (x.west,y.south,y.construction.latitude) given in radians.
c
c---------------------------------------------------------------------
c DNMI/FoU  15.05.1996  Anstein Foss
c---------------------------------------------------------------------
c
      integer icall, n, ierror
      real    x(n), y(n), xw, ys, yc
c
ccc   data rearth/6371220.0/
      data rearth/-1.0/
c
      ierror = 0
c
      if(rearth.le.0.) call earthr(rearth)
c
      zpih   = asin(1.0)
      zrcos  = rearth*cos(yc)
      zxmerc = zrcos*xw
      zymerc = zrcos*log((1.+sin(ys))/cos(ys))
c
c  test if possible to avoid heavy computations
c  dependant on y or latitude
      if(n.lt.4) then
	ieq=0
      elseif(y(1).eq.y(2) .or. y(2).eq.y(3) .or. y(3).eq.y(4)) then
	ieq=1
      else
	ieq=0
      end if
c
      if (icall.eq.1) then
c
c  compute spherical coordinates as function of
c  mercator (unrotated) coordinates
c
      if (ieq.eq.0) then
         do j = 1,n
	    xmer = zxmerc + x(j)
	    ymer = zymerc + y(j)
	    xsph = xmer/zrcos
	    ysph = 2. * atan(exp(ymer/zrcos))-zpih
            x(j) = xsph
            y(j) = ysph
         enddo
      else
	 ypos=-1.e+35
	 if(ypos.eq.y(1)) ypos=0.
         do j = 1,n
	    xmer = zxmerc + x(j)
	    xsph = xmer/zrcos
            x(j) = xsph
	    if (y(j).ne.ypos) then
	       ypos = y(j)
	       ymer = zymerc + y(j)
	       ysph = 2. * atan(exp(ymer/zrcos))-zpih
	    endif
            y(j) = ysph
         enddo
      endif
c
      elseif (icall.eq.-1) then
c
c  compute mercator (unrotated) coordinates as function of
c  spherical coordinates
c
      if (ieq.eq.0) then
         do j = 1,n
            xsph = x(j)
            ysph = y(j)
	    xmer = zrcos * xsph
	    ymer = zrcos * log((1.+sin(ysph))/cos(ysph))
	    x(j) = (xmer - zxmerc)
	    y(j) = (ymer - zymerc)
         enddo
      else
	 ysph=-1.e+35
	 if(ysph.eq.y(1)) ysph=0.
         do j = 1,n
            xsph = x(j)
	    xmer = zrcos * xsph
	    x(j) = (xmer - zxmerc)
	    if (y(j).ne.ysph) then
               ysph = y(j)
	       ymer = zrcos * log((1.+sin(ysph))/cos(ysph))
	       ypos = (ymer - zymerc)
	    endif
	    y(j) = ypos
         enddo
      endif
c
      else
c
      ierror = 1
c
      endif
c
      return
      end

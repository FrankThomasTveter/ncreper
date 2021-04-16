      subroutine xyconst(igtyp1,grid1,igtyp2,grid2,cxy1,cxy2,ierror)
c
c        hvis mulig beregnes konstanter for omregning av
c        kooordinater mellom 2 kart-projeksjoner.
c
c        grid typer:  1=polarstereografisk
c                     2=geografisk
c                     3=sfaerisk (rotert)
c                     4=polarstereografisk
c                     5=merkator (urotert) 
c
c        x2 = cxy1(1) + cxy1(2)*x1 + cxy1(3)*y1
c        y2 = cxy1(4) + cxy1(5)*x1 + cxy1(6)*y1
c
c        x1 = cxy2(1) + cxy2(2)*x2 + cxy2(3)*y2
c        y1 = cxy2(4) + cxy2(5)*x2 + cxy2(6)*y2
c
c        input:  igtyp1,grid1(6)
c                igtyp2,grid2(6)
c        output: cxy1(6)
c                cxy2(6)
c                ierror : 0 = ok, cxy1,cxy2 er beregnet
c                         1 = umulig konvertering
c                         2 = ukjent grid type
c
c-------------------------------------------------------------------
c  DNMI/FoU  xx.xx.1990  Anstein Foss ... polarstereografisk
c  DNMI/FoU  06.06.1996  Anstein Foss ... gridtype 1,2,3,4,5
c-------------------------------------------------------------------
c
      integer igtyp1,igtyp2,ierror
      real    grid1(6),grid2(6),cxy1(6),cxy2(6)
c
      real    alfa,calfa,salfa,s,zpir18
c
      ierror=2
c
      if((igtyp1.eq.1 .or. igtyp1.eq.4) .and.
     +   (igtyp2.eq.1 .or. igtyp2.eq.4)) then
c
	ierror=1
c
	if((grid1(5).gt.0. .and. grid2(5).gt.0.) .or.
     +     (grid1(5).lt.0. .and. grid2(5).lt.0.)) then
c
c..polarstereographic (rotation allowed)
c
	  zpir18 = 2.0*asin(1.0)/180.
c
          alfa=(grid1(4)-grid2(4))*zpir18
	  if(grid1(5).lt.0.) alfa=-alfa
          calfa=cos(alfa)
          salfa=sin(alfa)
          s=grid2(3)/grid1(3)
          cxy1(2)=+s*calfa
          cxy1(3)=-s*salfa
          cxy1(5)=+s*salfa
          cxy1(6)=+s*calfa
          cxy1(1)=grid2(1)-cxy1(2)*grid1(1)-cxy1(3)*grid1(2)
          cxy1(4)=grid2(2)-cxy1(5)*grid1(1)-cxy1(6)*grid1(2)
c
          alfa=(grid2(4)-grid1(4))*zpir18
	  if(grid1(5).lt.0.) alfa=-alfa
          calfa=cos(alfa)
          salfa=sin(alfa)
          s=grid1(3)/grid2(3)
          cxy2(2)=+s*calfa
          cxy2(3)=-s*salfa
          cxy2(5)=+s*salfa
          cxy2(6)=+s*calfa
          cxy2(1)=grid1(1)-cxy2(2)*grid2(1)-cxy2(3)*grid2(2)
          cxy2(4)=grid1(2)-cxy2(5)*grid2(1)-cxy2(6)*grid2(2)
c
	  ierror=0
c
	end if
c
      elseif((igtyp1.eq.2 .or. igtyp1.eq.3) .and.
     +       (igtyp2.eq.2 .or. igtyp2.eq.3)) then
c
c..geographic, spherical rotated (no rotation)
c
	ierror=1
c
	if(grid1(5).eq.grid2(5) .and. grid1(6).eq.grid2(6)) then
c
          cxy1(1)=(grid1(1)-grid1(3)-grid2(1)+grid2(3))/grid2(3)
          cxy1(4)=(grid1(2)-grid1(4)-grid2(2)+grid2(4))/grid2(4)
          cxy1(2)=grid1(3)/grid2(3)
          cxy1(6)=grid1(4)/grid2(4)
          cxy1(3)=0.
          cxy1(5)=0.
c
          cxy2(1)=(grid2(1)-grid2(3)-grid1(1)+grid1(3))/grid1(3)
          cxy2(4)=(grid2(2)-grid2(4)-grid1(2)+grid1(4))/grid1(4)
          cxy2(2)=grid2(3)/grid1(3)
          cxy2(6)=grid2(4)/grid1(4)
          cxy2(3)=0.
          cxy2(5)=0.
c
	  ierror=0
c
	end if
c
      elseif(igtyp1.eq.5 .and. igtyp2.eq.5) then
c
c..mercator (unrotated)
c
	zpir18 = 2.0*asin(1.0)/180.
        xw1 = grid1(1)*zpir18
        ys1 = grid1(2)*zpir18
        dx1 = grid1(3)*1000.
        dy1 = grid1(4)*1000.
        yc1 = grid1(5)*zpir18
        xw2 = grid2(1)*zpir18
        ys2 = grid2(2)*zpir18
        dx2 = grid2(3)*1000.
        dy2 = grid2(4)*1000.
        yc2 = grid2(5)*zpir18
	if(xw1.ne.xw2 .or. ys1.ne.ys2) then
	  x1 = xw1
	  y1 = ys1
	  call mer2sph(-1,1,x1,y1,xw2,ys2,dx2,dy2,yc2,ierror)
	  x2 = xw2
	  y2 = ys2
	  call mer2sph(-1,1,x2,y2,xw1,ys1,dx1,dy1,yc1,ierror)
	else
	  x1 = 1.
	  y1 = 1.
	  x2 = 1.
	  y2 = 1.
	end if
	cyc1=cos(yc1)
	cyc2=cos(yc2)
c
	cxy1(2)=(cyc2*dx1)/(cyc1*dx2)
	cxy1(6)=(cyc2*dy1)/(cyc1*dy2)
	cxy1(3)=0.
	cxy1(5)=0.
        cxy1(1)=x1-cxy1(2)-cxy1(3)
        cxy1(4)=y1-cxy1(5)-cxy1(6)
c
	cxy2(2)=(cyc1*dx2)/(cyc2*dx1)
	cxy2(6)=(cyc1*dy2)/(cyc2*dy1)
	cxy2(3)=0.
	cxy2(5)=0.
        cxy2(1)=x2-cxy2(2)-cxy2(3)
        cxy2(4)=y2-cxy2(5)-cxy2(6)
c
	ierror=0
c
      elseif(igtyp1.ge.1 .and. igtyp1.le.5 .and.
     +       igtyp2.ge.1 .and. igtyp2.le.5) then
c
	ierror=1
c
      end if
c
      return
      end

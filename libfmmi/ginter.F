      subroutine ginter(interp,iundef,undef,
     +                  nxi,nyi,fldi,gridin,
     +                  nxo,nyo,fldo,gridot,nundef)
c
c  interpolation in one polarstereografphic grid to another.
c  undefined values in input grid allowed.
c  undefined values in output grid due to area outside input
c  grid or undefined values in input grid.
c  for wind components:
c  use ginter for each component and then gintuv.
c
c  input:
c        interp: 1 = bilinear interpolation ( 4 input gridpoints used).
c                2 = bessel   interpolation (16 input gridpoints used),
c                    bilinear interpolation used near boundaries
c                    and undefined values.
c                3 = use value in nearest input gridpoint,
c                    no interpolation
c        iundef: 0 = no undefined values in input field
c                1 = undefined values in input field
c        undef:      undefined ('high') value (input and output)
c        gridin(4):  xp,yp,an,fi for input  grid
c        gridot(4):  xp,yp,an,fi for output grid
c        nxi,nyi:    dimensions of input field
c        fldi(nxi,nyi): input field
c        nxo,nyo:    dimensions of output field
c
c  output:
c        fldo(nxo,nyo): output field
c        nundef:        no. of undefined values in output field
c
C-------------------------------------------------------------------
c  DNMI/FoU   18.11.1992   Anstein Foss
c  DNMI/FoU   19.08.1993   Anstein Foss
C-------------------------------------------------------------------
c
      integer   interp,iundef,nxi,nyi,nxo,nyo
      real      undef
      real      fldi(nxi,nyi),gridin(4)
      real      fldo(nxo,nyo),gridot(4)
c
      parameter (ninter=512)
c
      integer   ipos(ninter),jpos(ninter),intp(ninter),nintp(0:3)
      real      xpos(ninter),ypos(ninter)
      real      alfa,calfa,salfa,scale
c
c  xin = cx1 + cx2*xout + cx3*yout
c  yin = cy1 + cy2*xout + cy3*yout
c
      alfa=(gridot(4)-gridin(4))*3.141592654/180.
      calfa=cos(alfa)
      salfa=sin(alfa)
      scale=gridin(3)/gridot(3)
      cx2=+scale*calfa
      cx3=-scale*salfa
      cy2=+scale*salfa
      cy3=+scale*calfa
      cx1=gridin(1)-cx2*gridot(1)-cx3*gridot(2)
      cy1=gridin(2)-cy2*gridot(1)-cy3*gridot(2)
c
c  check if interpolation is necessary
      icx1=nint(cx1)
      icx2=nint(cx2)
      icx3=nint(cx3)
      icy1=nint(cy1)
      icy2=nint(cy2)
      icy3=nint(cy3)
      ip=0
      epsilon=1.e-5
      if(abs(float(icx1)-cx1).gt.epsilon) ip=1
      if(abs(float(icx2)-cx2).gt.epsilon) ip=1
      if(abs(float(icx3)-cx3).gt.epsilon) ip=1
      if(abs(float(icy1)-cy1).gt.epsilon) ip=1
      if(abs(float(icy2)-cy2).gt.epsilon) ip=1
      if(abs(float(icy3)-cy3).gt.epsilon) ip=1
c
      if(ip.eq.0) goto 800
c
      nxim1=nxi-1
      nyim1=nyi-1
      nxim2=nxi-2
      nyim2=nyi-2
      xm=float(nxi)
      ym=float(nyi)
c
c  +/- small number due to computational accuracy
      epsilon=1.e-3
      xmin=1. - epsilon
      xmax=xm + epsilon
      ymin=1. - epsilon
      ymax=ym + epsilon
c
      inear =1
      ibilin=1
      ibessl=1
c  check corners in output grid
      do jo=1,nyo,nyo-1
        do io=1,nxo,nxo-1
          x=cx1+cx2*float(io)+cx3*float(jo)
          y=cy1+cy2*float(io)+cy3*float(jo)
          i=x
          j=y
          ni=nint(x)
          nj=nint(y)
          if(ni.lt.1 .or. ni.gt.nxi .or.
     +       nj.lt.1 .or. nj.gt.nyi) inear=0
          if(i.lt.1 .or. i.gt.nxim1 .or.
     +       j.lt.1 .or. j.gt.nyim1) ibilin=0
          if(i.lt.2 .or. i.gt.nxim2 .or.
     +       j.lt.2 .or. j.gt.nyim2) ibessl=0
        end do
      end do
c
      nundef=0
c
c  use fastest method when possible (warning: duplicated code)
      if(interp.eq.1 .and. ibilin.eq.1 .and. iundef.eq.0) goto 200
      if(interp.eq.2 .and. ibessl.eq.1 .and. iundef.eq.0) goto 300
      if(interp.eq.2 .and. ibilin.eq.1 .and. iundef.eq.0) goto 400
      if(interp.eq.3 .and. inear .eq.1 .and. iundef.eq.0) goto 500
c
c  following code covers cases with output gridpoints outside
c  the input grid and undefined values in the input grid
c
      ip=interp
      nud=0
      ud =undef
      udl=undef*0.9
c
      do jo=1,nyo
c
        do io1=1,nxo,ninter
          io2=min0(io1+ninter-1,nxo)
          io=io1-1
          no=io2-io
c
          do n=1,no
            xo=float(io+n)
            yo=float(jo)
            xpos(n)=cx1+cx2*xo+cx3*yo
            ypos(n)=cy1+cy2*xo+cy3*yo
            intp(n)=ip
          end do
          if(ip.eq.1 .or. ip.eq.2) then
            do n=1,no
              ipos(n)=xpos(n)
              jpos(n)=ypos(n)
            end do
          else
            do n=1,no
              ipos(n)=nint(xpos(n))
              jpos(n)=nint(ypos(n))
            end do
          end if
c
          if(ip.eq.2 .and. ibessl.eq.0) then
            do n=1,no
              if(ipos(n).lt.2 .or. ipos(n).gt.nxim2 .or.
     +           jpos(n).lt.2 .or. jpos(n).gt.nyim2) intp(n)=1
            end do
          end if
c
          if((ip.eq.1 .or. ip.eq.2) .and. ibilin.eq.0) then
            do n=1,no
              if(ipos(n).lt.1 .or. ipos(n).gt.nxim1 .or.
     +           jpos(n).lt.1 .or. jpos(n).gt.nyim1) then
                if(xpos(n).gt.xmin .and. xpos(n).lt.xmax .and.
     +             ypos(n).gt.ymin .and. ypos(n).lt.ymax) then
c  position outside input grid, probably due to computational accuracy
                  if(ipos(n).lt.1) then
                    ipos(n)=1
                    xpos(n)=1.
                  end if
                  if(ipos(n).gt.nxim1) then
                    ipos(n)=nxim1
                    xpos(n)=xm
                  end if
                  if(jpos(n).lt.1) then
                    jpos(n)=1
                    ypos(n)=1.
                  end if
                  if(jpos(n).gt.nyim1) then
                    jpos(n)=nyim1
                    ypos(n)=ym
                  end if
                else
                  intp(n)=0
                end if
              end if
            end do
          elseif(ip.eq.3 .and. inear.eq.0) then
            do n=1,no
              if(ipos(n).lt.1 .or. ipos(n).gt.nxi .or.
     +           jpos(n).lt.1 .or. jpos(n).gt.nyi) intp(n)=0
            end do
          end if
c
          if(iundef.ne.0 .and. ip.eq.2) then
            do n=1,no
              if(intp(n).eq.2) then
                i=ipos(n)
                j=jpos(n)
                if(amax1(fldi(i-1,j-1),fldi(i  ,j-1),
     +                   fldi(i+1,j-1),fldi(i+2,j-1),
     +                   fldi(i-1,j  ),fldi(i  ,j  ),
     +                   fldi(i+1,j  ),fldi(i+2,j  ),
     +                   fldi(i-1,j+1),fldi(i  ,j+1),
     +                   fldi(i+1,j+1),fldi(i+2,j+1),
     +                   fldi(i-1,j+2),fldi(i  ,j+2),
     +                   fldi(i+1,j+2),fldi(i+2,j+2)).gt.udl) intp(n)=1
              end if
            end do
          end if
c
          if(iundef.ne.0 .and. (ip.eq.1 .or. ip.eq.2)) then
            do n=1,no
              if(intp(n).eq.1) then
                i=ipos(n)
                j=jpos(n)
                if(amax1(fldi(i,j  ),fldi(i+1,j  ),
     +                   fldi(i,j+1),fldi(i+1,j+1)).gt.udl) intp(n)=0
              end if
            end do
          elseif(iundef.ne.0 .and. ip.eq.3) then
            do n=1,no
              if(intp(n).eq.3) then
                i=ipos(n)
                j=jpos(n)
                if(fldi(i,j).gt.udl) intp(n)=0
              end if
            end do
          end if
c
	  nintp(0)=0
	  nintp(1)=0
	  nintp(2)=0
	  nintp(3)=0
	  do n=1,no
	    nintp(intp(n))=nintp(intp(n))+1
	  end do
c
	  if(nintp(0).gt.0) then
c  not able to interpolate, set undefined
	    do n=1,no
	      if(intp(n).eq.0) fldo(io+n,jo)=ud
	    end do
            nud=nud+nintp(0)
	  end if
c
	  if(nintp(1).gt.0) then
c  bilinear interpolation
            do n=1,no
              if(intp(n).eq.1) then
                i=ipos(n)
                j=jpos(n)
                x=xpos(n)
                y=ypos(n)
                x1=x-float(i)
                y1=y-float(j)
                fldo(io+n,jo)= fldi(i,j)
     +                        +y1*(fldi(i,j+1)-fldi(i,j))
     +                        +x1*(fldi(i+1,j)-fldi(i,j))
     +                        +x1*y1*( fldi(i+1,j+1)-fldi(i+1,j)
     +                                -fldi(i,j+1)  +fldi(i,j))
	      end if  
	    end do
	  end if
c
	  if(nintp(2).gt.0) then
c  bessel interpolation
            do n=1,no
              if(intp(n).eq.2) then
                i=ipos(n)
                j=jpos(n)
                x=xpos(n)
                y=ypos(n)
                x1=x-float(i)
                y1=y-float(j)
                x2=1.-x1
                y2=1.-y1
                y3=-0.25*y1*y2
                y4=-0.1666667*y1*y2*(y1-0.5)
                c1=y3-y4
                c2=y2-y3+3.*y4
                c3=y1-y3-3.*y4
                c4=y3+y4
                t1= c1*fldi(i-1,j-1)+c2*fldi(i-1,j  )
     +             +c3*fldi(i-1,j+1)+c4*fldi(i-1,j+2)
                t2= c1*fldi(i  ,j-1)+c2*fldi(i  ,j  )
     +             +c3*fldi(i  ,j+1)+c4*fldi(i  ,j+2)
                t3= c1*fldi(i+1,j-1)+c2*fldi(i+1,j  )
     +             +c3*fldi(i+1,j+1)+c4*fldi(i+1,j+2)
                t4= c1*fldi(i+2,j-1)+c2*fldi(i+2,j  )
     +             +c3*fldi(i+2,j+1)+c4*fldi(i+2,j+2)
                fldo(io+n,jo)= x2*t2+x1*t3
     +                        -0.25*x2*x1*(t4-t3-t2+t1)
     +                        -0.1666667*x2*x1*(x1-0.5)
     +                                  *(t4-3.*t3+3.*t2-t1)
              end if
            end do
	  end if
c
	  if(nintp(3).gt.0) then
c  nearest gridpoint
            do n=1,no
              if(intp(n).eq.3) then
                i=ipos(n)
                j=jpos(n)
                fldo(io+n,jo)= fldi(i,j)
	      end if  
	    end do
	  end if
c
        end do
c
      end do
c
      nundef=nud
c
      return
c
c  bilinear interpolation for all output grid points
c
  200 do jo=1,nyo
        do io=1,nxo
          x=cx1+cx2*float(io)+cx3*float(jo)
          y=cy1+cy2*float(io)+cy3*float(jo)
          i=x
          j=y
          x1=x-float(i)
          y1=y-float(j)
          fldo(io,jo)=fldi(i,j)+y1*(fldi(i,j+1)-fldi(i,j))
     +                         +x1*(fldi(i+1,j)-fldi(i,j))
     +                         +x1*y1*( fldi(i+1,j+1)-fldi(i+1,j)
     +                                 -fldi(i,j+1)  +fldi(i,j))
        end do
      end do
c
      return
c
c  bessel interpolation for all output gridpoints
c
  300 do jo=1,nyo
        do io=1,nxo
          x=cx1+cx2*float(io)+cx3*float(jo)
          y=cy1+cy2*float(io)+cy3*float(jo)
          i=x
          j=y
          x1=x-float(i)
          y1=y-float(j)
          x2=1.-x1
          y2=1.-y1
          y3=-0.25*y1*y2
          y4=-0.1666667*y1*y2*(y1-0.5)
          c1=y3-y4
          c2=y2-y3+3.*y4
          c3=y1-y3-3.*y4
          c4=y3+y4
          t1= c1*fldi(i-1,j-1)+c2*fldi(i-1,j  )
     +       +c3*fldi(i-1,j+1)+c4*fldi(i-1,j+2)
          t2= c1*fldi(i  ,j-1)+c2*fldi(i  ,j  )
     +       +c3*fldi(i  ,j+1)+c4*fldi(i  ,j+2)
          t3= c1*fldi(i+1,j-1)+c2*fldi(i+1,j  )
     +       +c3*fldi(i+1,j+1)+c4*fldi(i+1,j+2)
          t4= c1*fldi(i+2,j-1)+c2*fldi(i+2,j  )
     +       +c3*fldi(i+2,j+1)+c4*fldi(i+2,j+2)
          fldo(io,jo)= x2*t2+x1*t3
     +                -0.25*x2*x1*(t4-t3-t2+t1)
     +                -0.1666667*x2*x1*(x1-0.5)
     +                          *(t4-3.*t3+3.*t2-t1)
        end do
      end do
c
      return
c
c  bessel interpolation with bilinear interpolation near boundaries
c
  400 do jo=1,nyo
        do io=1,nxo
          x=cx1+cx2*float(io)+cx3*float(jo)
          y=cy1+cy2*float(io)+cy3*float(jo)
          i=x
          j=y
          x1=x-float(i)
          y1=y-float(j)
          if(i.gt.1 .and. i.lt.nxim1 .and.
     +       j.gt.1 .and. j.lt.nyim1) then
c  bessel interpolation
            x2=1.-x1
            y2=1.-y1
            y3=-0.25*y1*y2
            y4=-0.1666667*y1*y2*(y1-0.5)
            c1=y3-y4
            c2=y2-y3+3.*y4
            c3=y1-y3-3.*y4
            c4=y3+y4
            t1= c1*fldi(i-1,j-1)+c2*fldi(i-1,j  )
     +         +c3*fldi(i-1,j+1)+c4*fldi(i-1,j+2)
            t2= c1*fldi(i  ,j-1)+c2*fldi(i  ,j  )
     +         +c3*fldi(i  ,j+1)+c4*fldi(i  ,j+2)
            t3= c1*fldi(i+1,j-1)+c2*fldi(i+1,j  )
     +         +c3*fldi(i+1,j+1)+c4*fldi(i+1,j+2)
            t4= c1*fldi(i+2,j-1)+c2*fldi(i+2,j  )
     +         +c3*fldi(i+2,j+1)+c4*fldi(i+2,j+2)
            fldo(io,jo)= x2*t2+x1*t3
     +                  -0.25*x2*x1*(t4-t3-t2+t1)
     +                  -0.1666667*x2*x1*(x1-0.5)
     +                            *(t4-3.*t3+3.*t2-t1)
          else
c  bilinear interpolation near boundaries
            fldo(io,jo)=fldi(i,j)+y1*(fldi(i,j+1)-fldi(i,j))
     +                           +x1*(fldi(i+1,j)-fldi(i,j))
     +                           +x1*y1*( fldi(i+1,j+1)-fldi(i+1,j)
     +                                   -fldi(i,j+1)  +fldi(i,j))
          end if
        end do
      end do
c
      return
c
c  nearest gridpoint for all output grid points
c
  500 do jo=1,nyo
        do io=1,nxo
          x=cx1+cx2*float(io)+cx3*float(jo)
          y=cy1+cy2*float(io)+cy3*float(jo)
          i=nint(x)
          j=nint(y)
          fldo(io,jo)=fldi(i,j)
        end do
      end do
c
      return
c
c  interpolation not needed
c
c  check corners in output grid
  800 ip=0
      do jo=1,nyo,nyo-1
        do io=1,nxo,nxo-1
	  i=icx1+icx2*io+icx3*jo
	  j=icy1+icy2*io+icy3*jo
          if(i.lt.1 .or. i.gt.nxi .or.
     +       j.lt.1 .or. j.gt.nyi) ip=1
        end do
      end do
c
      nud=0
      ud =undef
      udl=undef*0.9
c
      if(ip.eq.0 .and. icx2.eq.1 .and. icx3.eq.0 
     +           .and. icy2.eq.0 .and. icy3.eq.1) then
        do jo=1,nyo
	  do io=1,nxo
	    fldo(io,jo)=fldi(icx1+io,icy1+jo)
	  end do
	end do
      elseif(ip.eq.0) then
        do jo=1,nyo
	  do io=1,nxo
	    i=icx1+icx2*io+icx3*jo
	    j=icy1+icy2*io+icy3*jo
	    fldo(io,jo)=fldi(i,j)
	  end do
	end do
      else
        do jo=1,nyo
	  do io=1,nxo
	    i=icx1+icx2*io+icx3*jo
	    j=icy1+icy2*io+icy3*jo
	    if(i.gt.0 .and. i.le.nxi .and. 
     +         j.gt.0 .and. j.le.nyi) then
	      fldo(io,jo)=fldi(i,j)
	    else
	      fldo(io,jo)=ud
	    end if
	  end do
	end do
      end if
c
      if(ip.ne.0 .or. iundef.ne.0) then
        do jo=1,nyo
	  do io=1,nxo
	    if(fldo(io,jo).gt.udl) nud=nud+1
	  end do
	end do
      end if
c
      nundef=nud
c
      return
      end

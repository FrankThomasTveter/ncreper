      subroutine int2pos(igtype,grid,nx,ny,nfelt,fdata,kpar,parlim,
     +                   newpos,npos,
     +                   inter,iundef,pos,geopos,gxypos,dat,
     +                   lfelt,felt,
     +                   jinter,ndimri,rinter,vturn,dwork,ierror)
c
c        interpolasjon i (input) felt til posisjoner.
c        polarstereografiske, geografiske, rotert sfaeriske
c        og andre grid kan innleses.
c
c-----------------------------------------------------------------------
c  input:   ( <i>=integer  <i2>=integer*2  <r>=real )
c  ------
c  <i>  igtype           - grid type: 1 = polarstereografisk grid (60N)
c                                     2 = geografisk grid
c                                     3 = sfaerisk rotert grid
c                                     4 = polarstereografisk grid
c                                     5 = merkator (urotert) grid
c                                     0 = ukjent
c  <r>  grid(6)          - polarstereografisk grid: xp,yp,an,fi,fp,0.
c                          geografisk grid:         l1,b1,dl,db,0.,0.
c                          sfaerisk rotert grid:    l1,b1,dl,db,cl,cb
c                          merkator grid:           l1,b1,dx,dy,bc,0.
c                          ukjent grid:             0.,0.,0.,0.,0.,0.
c  <i>  nx               - x-dimensjon av felt
c  <i>  ny               - y-dimensjon av felt
c  <i>  nfelt            - antall felt
c  <r>  fdata(nx,ny,nfelt) - feltene
c  <i>  newpos           - 1=nye posisjoner i pos-arrayen
c                          0=samme posisjoner som sist i pos-arrayen
c                            (hvis samme grid og interpolasjonstype som
c                             sist, beregnes ikke konstanter m.m. for
c                             interpolasjon paa nytt, dvs. arrayene
c                             jinter,rinter,vturn,geopos og gxypos)
c  <i>  npos             - antall posisjoner
c  <i>  inter            - interpolasjons-type:
c                               1 = bilinear
c                               2 = biquadratic
c                               3 = bicubic (bessel)
c                               4 = nearest gridpoint
c  <i>  iundef           - check undefined in input fields,
c                               0 = never
c                               1 = in first field
c                               2 = in each  field
c  <r>  pos(3,npos)      - pos(1,n): posisjons-type;
c                                    1=x,y i grid
c                                    2=geografisk bredde,lengde (desimalt)
c                                    (alle andre verdier behandlet som 1)
c                          pos(2,n): 1. koordinat (x/lat)
c                          pos(3,n): 2. koordinat (y/lon)
c  <i>  kpar(3,nfelt)    - kpar(1,n): parameter-kode:
c                                        0=skalar
c                                        1=x.vektor (u) ) roteres til
c                                        2=y.vektor (v) ) e/w og n/s
c                                        3=retning relativt nord (dd)
c                                        4=retning relativt nord (dd),
c                                          retningen roteres 180 grader
c                          kpar(2,n): 0 = not test minimum value after
c                                         interpolation
c                                     1 = test minimum value after interp.,
c                                         minimum = parlim(1,n)
c                          kpar(3,n): 0 = not test maximum value after
c                                         interpolation
c                                     1 = test maximum value after interp.,
c                                         maximum = parlim(2,n)
c  <r>  parlim(2,nfelt)  - parlim(1,n) = minimum value if kpar(2,n)=1
c                          parlim(2,n) = maximum value if kpar(3,n)=1
c  <i>  lfelt            - lengde av 'felt'
c  <r>  felt(lfelt,2)    - work array
c  <i>  jinter(npos,3)      -
c  <r>  rinter(ndimri,npos) -
c  <r>  vturn(4,npos)       -
c  <r>  dwork(npos)      - work array
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  output:   ( <i>=integer  <i2>=integer*2  <r>=real )
c  ------
c  <r>  geopos(2,npos)   - geopos(1,n) : geografisk bredde
c                          geopos(2,n) : geografisk lengde
c  <r>  gxypos(2,npos)   - gxypos(1,n) : x-posisjon i grid/felt
c                          gxypos(2,n) : y-posisjon i grid/felt
c  <r>  dat(npos,nfelt)  - interpolerte verdier
c  <i>  ierror           - feil-kode, 0=ingen feil
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c       loekker for innlesning og interpolasjon:
c
c       i---- felt ............................ np=1,nfelt
c       i      <innlesning av felt, 'rfelt'>
c       i i-- posisjon ........................ n=1,npos
c       i i    <interpolasjon, 'interp'> ...... dat(n,np)
c       i-i-------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c  Other DNMI library subroutines: extendf
c                                  grd2pos (calling misc. other routines)
c                                  grv2pos (calling misc. other routines)
c
c-----------------------------------------------------------------------
c  DNMI/FoU  11.10.1989  Anstein Foss .. fltint
c  DNMI/FoU  08.03.1990  Anstein Foss
c  DNMI/FoU  27.04.1990  Anstein Foss
c  DNMI/FoU  12.03.1994  Anstein Foss
c  DNMI/FoU  21.03.1995  Anstein Foss .. flt2pos (also for Hirlam grids)
c  DNMI/FoU  25.04.1995  Anstein Foss .. int2pos
c  DNMI/FoU  14.06.1996  Anstein Foss .. new call, new grd2pos,grv2pos
c-----------------------------------------------------------------------
c
      integer   igtype, nx, ny, nfelt
      integer   newpos, npos
      integer   inter, lfelt, ndimri, ierror
      integer   jinter(npos,3)
      integer   kpar(3,nfelt)
      real      fdata(nx*ny,nfelt)
      real      pos(3,npos), geopos(2,npos), gxypos(2,npos)
      real      parlim(2,nfelt), dat(npos,nfelt), grid(6)
      real      felt(lfelt,2)
      real      rinter(ndimri,npos), vturn(4,npos)
      real      dwork(npos)
c
      integer   lastnx,lastny,lastgtype,lastinter,lastnpos
      integer   icall1,icall2,mvx,mvy,nxa,nya
      real      grida(6),gridlast(6)
      logical   lexten
c
      data lastnx,lastny,lastgtype,lastinter,lastnpos/5*-32767/
      data gridlast/6*1.e+35/
      data grida/6*1.e+35/
c
      data icall1,icall2/2,2/
      data mvx,mvy,nxa,nya/0,0,0,0/
      data lexten/.false./
c
      ierror=0
c
      udef=+1.e+35
c
c..sjekk parameter-koder
      k=1
      do np=1,nfelt
        if(kpar(1,np).lt.0 .or. kpar(1,np).gt.4) then
          write(6,*)
          write(6,*) '**int2pos** ukjent parameter-kode: ',kpar(1,np)
          k=0
        elseif(kpar(1,np).eq.1) then
          if(np.eq.nfelt) then
            write(6,*) '**int2pos** feil rekkefoelge u/v beregning'
            k=0
          elseif(kpar(1,np+1).ne.2) then
            write(6,*) '**int2pos** feil rekkefoelge u/v beregning'
            k=0
          end if
        elseif(kpar(1,np).eq.2) then
          if(np.eq.1) then
            write(6,*) '**int2pos** feil rekkefoelge u/v beregning'
            k=0
          elseif(kpar(1,np-1).ne.1) then
            write(6,*) '**int2pos** feil rekkefoelge u/v beregning'
            k=0
          end if
        end if
      end do
      if(k.eq.0) then
        write(6,*)
        write(6,*) '**int2pos** ingen innlesning/interpolasjon.'
        write(6,*)
        return
      end if
c
      newgrd=0
      do j=1,6
        if(grid(j).ne.gridlast(j)) newgrd=1
      end do
c
      if(newpos.eq.1         .or.
     +   inter .ne.lastinter .or.
     +   npos  .ne.lastnpos  .or.
     +   igtype.ne.lastgtype .or.
     +   nx    .ne.lastnx    .or.
     +   ny    .ne.lastny    .or.
     +   newgrd.eq.1)        then
        lastnx = nx
        lastny = ny
        lastinter=inter
        lastnpos =npos
	lastgtype=igtype
        do j=1,6
          grida(j)=grid(j)
          gridlast(j)=grida(j)
        end do
        icall1=2
        icall2=2
        lexten=.false.
        mvx=0
        mvy=0
        newgrd=1
      end if
c
c..field loop
      do 100 np=1,nfelt
c
      do j = 1,nx*ny
        felt(j,1) = fdata(j,np)
      end do
c
      if(newgrd.eq.1) then
c
        icall1=2
        icall2=2
        lexten=.false.
        mvx=0
        mvy=0
c
        nxa = nx
        nya = ny
c
c check if extention of input field for hemispheric data
        if (igtype.eq.2 .or. igtype.eq.3) then
           iwest = nint((grida(1)-grida(3)+360.)*100.)
           ieast = nint((grida(1)+float(nx-1)*grida(3))*100.)
           ilat1 = nint(abs(grida(2)*100.))
           ilat2 = nint(abs((grida(2)+float(ny-1)*grida(4))*100.))
           if (iwest.eq.ieast .and.
     +        (ilat1.eq.9000 .or. ilat2.eq.9000)) then
              if ((nx+6)*(ny+6).le.lfelt) then
                 lexten = .true.
                 mvx=+3
                 mvy=+3
                 nxa = nx+6
                 nya = ny+6
                 grida(1) = grida(1) - 3.*grida(3)
                 grida(2) = grida(2) - 3.*grida(4)
              else
                 write(*,*) 'int2pos: WARNING: not able to extend field'
                 write(*,*) '   extended nx*ny = ',(nx+6)*(ny+6)
                 write(*,*) '      (max) lfelt = ',lfelt
              end if
           end if
        end if
c
        newgrd = 0
c
      end if
c
c extend input field for hemispheric data
      if (lexten) then
         do i = 1,(nxa-6)*(nya-6)
            felt(i,2) = felt(i,1)
         end do
         call extendf(felt(1,2),felt(1,1),nxa-6,nya-6)
      end if
c
c true north direction parameter (dd)
c interpolated as vector/velocity components
      if (kpar(1,np).eq.3 .or. kpar(1,np).eq.4) then
         zpir18 = 2.0*asin(1.0)/180.
         do j = 1,nxa*nya
            if (felt(j,1).ne.udef) then
               dd = felt(j,1)*zpir18
               felt(j,1) = sin(dd)
               felt(j,2) = cos(dd)
            else
               felt(j,2) = udef
            end if
         end do
      end if
c
c possibly check output min and max values
c possibly enable check in interpolation routine
      ipzero=0
      if (kpar(2,np).eq.1 .and. parlim(1,np).eq.0.0) ipzero=1
      if (kpar(3,np).eq.1 .and. parlim(2,np).eq.0.0) ipzero=1
c
c interpolate to positions
      call grd2pos(icall1,igtype,grida,mvx,mvy,
     +             felt(1,1),nxa,nya,
     +             npos,pos(1,1),geopos(1,1),gxypos(1,1),
     +             dat(1,np),jinter(1,1),rinter(1,1),ndimri,
     +             udef,inter,ipzero,iundef,ierror)
      if(ierror.ne.0) write(6,*) 'INT2POS: GRD2POS ERROR: ',ierror
      if(ierror.ne.0) return
      icall1=3
c
      if (kpar(1,np).eq.2 .and. np.gt.1) then
c turn velocity components
         call grv2pos(icall2,igtype,grida,
     +                dat(1,np-1),dat(1,np),npos,
     +                vturn(1,1),geopos(1,1),udef,ierror)
         if(ierror.ne.0) write(6,*) 'INT2POS: GRV2POS ERROR: ',ierror
         if(ierror.ne.0) return
         icall2=3
      end if
c
c true north direction parameter (dd)
      if (kpar(1,np).eq.3 .or. kpar(1,np).eq.4) then
c interpolate v component
         do n = 1,npos
            dwork(n)=dat(n,np)
         end do
         call grd2pos(icall1,igtype,grida,mvx,mvy,
     +                felt(1,2),nxa,nya,
     +                npos,pos(1,1),geopos(1,1),gxypos(1,1),
     +                dat(1,np),jinter(1,1),rinter(1,1),ndimri,
     +                udef,inter,ipzero,iundef,ierror)
         if(ierror.ne.0) write(6,*) 'INT2POS: GRD2POS ERROR: ',ierror
         if(ierror.ne.0) return
         zdeg18 = 180./(2.0*asin(1.0))
ccc      iscale=ifelt(20)
ccc      ddmin = 10.**iscale
         ddturn = 0.
         if (kpar(1,np).eq.4) ddturn=180.
         do n = 1,npos
            if (dwork(n).ne.udef) then
               dd = zdeg18*atan2(dat(n,np),dwork(n)) + ddturn
               if (dd.lt.  0.) dd=dd+360.
               if (dd.gt.360.) dd=dd-360.
c..assure nonzero value when (possibly later) scaling to integer*2
ccc            dd = max(dd,ddmin)
               dat(n,np) = dd
            else
               dat(n,np) = udef
            end if
         end do
      end if
c
      if (kpar(2,np).eq.1) then
c check min value
         do n = 1,npos
            if (dat(n,np).ne.udef)
     +          dat(n,np)=max(dat(n,np),parlim(1,np))
         end do
      end if
      if (kpar(3,np).eq.1) then
c check max value
         do n = 1,npos
            if (dat(n,np).ne.udef)
     +          dat(n,np)=min(dat(n,np),parlim(2,np))
         end do
      end if
c
  100 continue
c
      return
      end

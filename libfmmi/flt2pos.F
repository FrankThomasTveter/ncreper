      subroutine flt2pos(messag,istop1,istop2,
     +                   newpos,npos,npar,nlev,ntim,
     +                   inter,iundef,pos,geopos,gxypos,dat,
     +                   ipar,ilev,itim,iprod,igrid,
     +                   parlim,lfelt,ifelt,felt,itime,
     +                   idsave,nx,ny,igtype,grid,
     +                   jinter,ndimri,rinter,vturn,dwork,ltim,ierror)
c
c        innlesning av felt fra 'felt-file'
c        og interpolasjon til posisjoner.
c        polarstereografiske, geografiske, rotert sfaeriske
c        og andre grid kan innleses.
c
c-----------------------------------------------------------------------
c  input:   ( <i>=integer  <i2>=integer*2  <r>=real )
c  ------
c  <i>  messag           - 0=ingen utskrift (hvis alt o.k.)
c                          1=skriver ut info. om hver tid som behandles
c                          2=skriver ut info. om hvert felt (for test)
c                         -1=ingen utskrift (selv om felt mangler)
c  <i>  istop1           - 1=stopp hvis posisjon er utenfor gridet
c                          0=ikke stopp (jinter(n,3)=0 => pos. utenfor)
c  <i>  istop2           - 0=ikke stopp hvis felt ikke finnes
c                          1=stopp hvis felt ikke finnes
c                          2=hvis felt mangler: hopp til neste tid
c                          3=hvis felt mangler: hopp til neste nivaa
c  <i>  newpos           - 1=nye posisjoner i pos-arrayen
c                          0=samme posisjoner som sist i pos-arrayen
c                            (hvis samme grid og interpolasjonstype som
c                             sist, beregnes ikke konstanter m.m. for
c                             interpolasjon paa nytt, dvs. arrayene
c                             jinter,rinter,vturn,geopos og gxypos)
c  <i>  npos             - antall posisjoner
c  <i>  npar             - antall parametre
c  <i>  nlev             - antall nivaa
c  <i>  ntim             - antall tidspunkt/filer
c  <i>  inter            - interpolasjons-type:
c                               1 = bilinear
c                               2 = biquadratic
c                               3 = bicubic (bessel)
c                               4 = nearest gridpoint
c                               5 = bilinear / nearest / nearest existing
c				    ("interp.ocean", for WAM etc.)
c                              44 = mean of the surrounding 4x4 gridpoints
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
c  <i2> ipar(5,npar)     - ipar(1,n): vertikal-koordinat
c                          ipar(2,n): parameter-nummer
c                          ipar(3,n): nivaa_1 ) brukes hvis nlev=1 !
c                          ipar(4,n): nivaa_2 )
c                          ipar(5,n): parameter-kode:
c                                        0=skalar
c                                        1=x.vektor (u) ) roteres til
c                                        2=y.vektor (v) ) e/w og n/s
c                                        3=retning relativt nord (dd)
c                                        4=retning relativt nord (dd),
c                                          retningen roteres 180 grader
c                                        5=sett alle verdier=0 hvis
c                                          itim(3,n)<=0 (prognose<=+0)
c                                        6=sett alle verdier=0
c                                          (dummy parameter)
c                          ipar(6,n): 0 = not test minimum value after
c                                         interpolation
c                                     1 = test minimum value after interp.,
c                                         minimum = parlim(1,n)
c                          ipar(7,n): 0 = not test maximum value after
c                                         interpolation
c                                     1 = test maximum value after interp.,
c                                         maximum = parlim(2,n)
c  <r>  parlim(2,npar)   - parlim(1,n) = minimum value if ipar(6,n)=1
c                          parlim(2,n) = maximum value if ipar(7,n)=1
c  <i2> ilev(2,nlev)     - ilev(1,n): nivaa_1 ) brukes ikke hvis nlev=1 !
c                          ilev(2,n): nivaa_2 ) nivaa_2=-32767 er tillatt
c  <i2> itim(3,ntim)     - itim(1,n): file unit (felt-file)
c                          itim(2,n): datatype
c                          itim(3,n): tidsparameter (prognose-tid)
c  <i>  iprod            - produsent-nummer
c  <i>  igrid            - grid-nummer
c  <i>  lfelt            - lengde av 'ifelt'
c  <i2> ifelt(20+lfelt)  - array for innlesning av felt
c                          inkl. 20 ord identifikasjon, og evt. ogsaa
c                          ekstra grid identifikasjon bak feltet.
c  <r>  felt(lfelt,2)       - work array
c  <i>  jinter(npos,3)      - work array
c  <r>  rinter(ndimri,npos) - work array
c  <r>  vturn(4,npos)       - work array
c  <r>  dwork(npos)         - work array
c
c       nb! hvis ipar(1,n),ilev(1,n) eller itim(1,n) er gitt verdien
c           -32767 gjoeres det ikke innlesning og interpolasjon
c           for disse.
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  output:   ( <i>=integer  <i2>=integer*2  <r>=real )
c  ------
c  <r>  geopos(2,npos)   - geopos(1,n) : geografisk bredde
c                          geopos(2,n) : geografisk lengde
c  <r>  gxypos(2,npos)   - gxypos(1,n) : x-posisjon i grid/felt
c                          gxypos(2,n) : y-posisjon i grid/felt
c  <r>  dat(npos,npar,nlev,ntim) - interpolerte verdier
c  <i>  itime(5,ntim)            - itime(1,n): aar (-32767: ikke funnet)
c                                  itime(2,n): maaned
c                                  itime(3,n): dag
c                                  itime(4,n): utc (timer), 'analysetid'
c                                  itime(5,n): prognose-tid (timer)
c  <i2> idsave(4,npar,nlev,ntim) - status og litt fra felt-identifikasjon
c                                  idsave(1,np,nl,nt) : 1 = felt funnet
c                                                       0 = felt ikke funnet
c                                  idsave(2,np,nl,nt) : ident(20), skalering
c                                  idsave(3,np,nl,nt) : ident(19) (sigma/B.eta)
c                                  idsave(4,np,nl,nt) : ident( 8) (A.eta)
c  <i>  nx                       - x-dimensjon av felt
c  <i>  ny                       - y-dimensjon av felt
c  <i>  igtype                   - grid type:
c                                     1 = polarstereografisk grid (60N)
c                                     2 = geografisk grid
c                                     3 = sfaerisk rotert grid
c                                     4 = polarstereografisk grid
c                                     5 = merkator (urotert) grid
c                                     0 = ukjent
c  <r>  grid(6)                  - polarstereografisk grid: xp,yp,an,fi,fp,0.
c                                  geografisk grid:         l1,b1,dl,db,0.,0.
c                                  sfaerisk rotert grid:    l1,b1,dl,db,cl,cb
c                                  merkator grid:           l1,b1,dx,dy,bc,0.
c                                  nb! hvis forskjellige grid (med samme
c                                  grid-nr.) gjelder dette sist leste grid.
c  <i>  ltim                     - siste innleste tidspunkt (0-ntim).
c                                  hvis 'istop2'=1 blir 'ltim' satt
c                                  til det siste tidspunkt hvor alle
c                                  felt er funnet, 0 hvis f@rste
c                                  tidspunkt ikke er o.k.
c  <i>  ierror                   - feil-kode, 0=ingen feil
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c       loekker for innlesning og interpolasjon:
c
c       i-------- tid ............................. nt=1,ntim
c       i i------ nivaa ........................... nl=1,nlev
c       i i i---- parameter ....................... np=1,npar
c       i i i      <innlesning av felt, 'rfelt'>
c       i i i i-- posisjon ........................ n=1,npos
c       i i i i    <interpolasjon, 'interp'> ...... dat(n,np,nl,nt)
c       i-i-i-i-------------------
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c  Other DNMI library subroutines: rfelt
c                                  rfturbo
c                                  gridpar
c                                  extendf
c                                  grd2pos (calling misc. other routines)
c                                  grv2pos (calling misc. other routines)
c
c-----------------------------------------------------------------------
c  DNMI/FoU  11.10.1989  Anstein Foss .. fltint
c  DNMI/FoU  08.03.1990  Anstein Foss
c  DNMI/FoU  27.04.1990  Anstein Foss
c  DNMI/FoU  12.03.1994  Anstein Foss .. rfelt(nlev=1) rfturbo(nlev>1)
c  DNMI/FoU  21.03.1995  Anstein Foss .. flt2pos (also for Hirlam grids)
c  DNMI/FoU  09.06.1995  Anstein Foss .. gridpar update
c  DNMI/FoU  29.12.1995  Anstein Foss .. ipar(5,n)=5,6 options
c  DNMI/FoU  14.06.1996  Anstein Foss .. new call, new grd2pos,grv2pos
c  DNMI/FoU  03.06.1997  Anstein Foss .. inter=44 (4x4 mean)
c  DNMI/FoU  22.09.1998  Anstein Foss .. inter=5 (for WAM)
c  DNMI/FoU  22.11.2000  Anstein Foss .. datatype, itim(2,n)=-32767
c-----------------------------------------------------------------------
c
      integer   messag, istop1, istop2
      integer   newpos, npos, npar, nlev, ntim
      integer   inter, iprod, igrid, lfelt, nx, ny, igtype, ndimri
      integer   ierror
      integer   itime(5,ntim)
      integer   jinter(npos,3)
      integer*2 ipar(7,npar), ilev(2,nlev), itim(3,ntim)
      integer*2 ifelt(20+lfelt), idsave(4,npar,nlev,ntim)
      real      pos(3,npos), geopos(2,npos), gxypos(2,npos)
      real      parlim(2,npar), dat(npos,npar,nlev,ntim), grid(6)
      real      felt(lfelt,2)
      real      rinter(ndimri,npos), vturn(4,npos)
      real      dwork(npos)
c
      parameter (maxhlp=50)
c
      integer   kontrl,ingrid
      integer   lastprod,lastgrid,lastinter,lastnpos
      integer   icall1,icall2,mvx,mvy,nxa,nya
      integer   ierr(3), ihelp(6,maxhlp)
      integer*2 in(16), idgrid(6)
      real      grida(6)
      logical   lexten
c
c.rfturbo................................................................
      integer   ihelpr(6)
      integer*2 idrec1r(1024),idrec2r(1024),innhrf(16,64)
c.rfturbo................................................................
c
      data lastprod,lastgrid,lastinter,lastnpos/4*-32767/
      data idgrid/6*-32767/
c
      data kontrl,ingrid/0,0/
      data icall1,icall2/2,2/
      data mvx,mvy,nxa,nya/0,0,0,0/
      data lexten/.false./
c
      ierror=0
c
      udef=+1.e+35
c
c..setter "manglende"
      do nt=1,ntim
        do nl=1,nlev
          do np=1,npar
            idsave(1,np,nl,nt)=0
            idsave(2,np,nl,nt)=0
            idsave(3,np,nl,nt)=0
            idsave(4,np,nl,nt)=0
          end do
        end do
      end do
c
      do nt=1,ntim
        do i=1,5
          itime(i,nt)=-32767
        end do
      end do
c
c..sjekk parameter-koder
      k=1
      do np=1,npar
        if(ipar(1,np).eq.-32767) then
          continue
        elseif(ipar(5,np).lt.0 .or. ipar(5,np).gt.6) then
          write(6,*)
          write(6,*) '**flt2pos** ukjent parameter-kode: ',ipar(5,np)
          k=0
        elseif(ipar(5,np).eq.1) then
          if(np.eq.npar) then
            write(6,*) '**flt2pos** feil rekkefoelge u/v beregning'
            k=0
          elseif(ipar(5,np+1).ne.2) then
            write(6,*) '**flt2pos** feil rekkefoelge u/v beregning'
            k=0
          end if
        elseif(ipar(5,np).eq.2) then
          if(np.eq.1) then
            write(6,*) '**flt2pos** feil rekkefoelge u/v beregning'
            k=0
          elseif(ipar(5,np-1).ne.1) then
            write(6,*) '**flt2pos** feil rekkefoelge u/v beregning'
            k=0
          end if
        end if
      end do
      if(k.eq.0) then
        write(6,*)
        write(6,*) '**flt2pos** ingen innlesning/interpolasjon.'
        write(6,*)
        ltim=0
        return
      end if
c
      ltt=0
      ltl=0
      ltp=0
c
      in(1)=iprod
      in(2)=igrid
c
      if(newpos.eq.1        .or.
     +   iprod.ne.lastprod  .or.
     +   igrid.ne.lastgrid  .or.
     +   inter.ne.lastinter .or.
     +   npos .ne.lastnpos) then
        lastprod =iprod
        lastgrid =igrid
        lastinter=inter
        lastnpos =npos
        kontrl=0
        ingrid=0
        icall1=2
        icall2=2
        lexten=.false.
        mvx=0
        mvy=0
      end if
c
      ifil=-32767
c
c..tids-loekke
      do 100 nt=1,ntim
      if(itim(1,nt).eq.-32767) goto 100
c
      if(ifil.ne.itim(1,nt)) then
        if(nlev.eq.1) then
          mhlp=min(npar,maxhlp)
          do n=1,mhlp
            ihelp(1,n)=-32767
          end do
        else
c.rfturbo................................................................
          do i=1,6
            ihelpr(i)=0
          end do
c.rfturbo................................................................
        end if
      end if
c
      ifil=itim(1,nt)
      in( 9)=itim(2,nt)
      in(10)=itim(3,nt)
      if(messag.eq.1) write(6,1010) ifil,in(1),in(2),in(9),in(10)
 1010 format(' file,produsent,grid,datatype,prog.tid: ',5i6)
c
c..nivaa-loekke
      do 110 nl=1,nlev
      if(ilev(1,nl).eq.-32767 .and. nlev.gt.1) goto 110
      in(13)=ilev(1,nl)
      in(14)=ilev(2,nl)
c
c..parameter-loekke (setter ogs$ inn f@rste niv$)
      do 120 np=1,npar
      if(ipar(1,np).eq.-32767) goto 120
      in( 9)=itim(2,nt)
      in(10)=itim(3,nt)
      in(11)=ipar(1,np)
      in(12)=ipar(2,np)
      if(nlev.eq.1) then
        in(13)=ipar(3,np)
        in(14)=ipar(4,np)
      end if
c
      if(messag.eq.2) write(6,1012) ifil,in(1),in(2),(in(i),i=9,14)
 1012 format(' file,in(1,2,9-14): ',i3,8i6)
c
      if((ipar(5,np).eq.5 .and. itim(3,nt).le.0)
     +                      .or. ipar(5,np).eq.6) then
c..all values = zero
        idsave(1,np,nl,nt)=1
        idsave(2,np,nl,nt)=0
        idsave(3,np,nl,nt)=0
        idsave(4,np,nl,nt)=0
        do n=1,npos
          dat(n,np,nl,nt)=0.
        end do
        ltp=nt
        goto 120
      end if
c
c..innlesning av felt
      if(nlev.eq.1) then
        nhlp=min(np,maxhlp)
        ih1=ihelp(1,nhlp)
        call rfelt(ifil,ip,in,ifelt,20+lfelt,ierr,ihelp(1,nhlp))
      else
        call rfturbo(ifil,ip,in,ifelt,20+lfelt,ierr,ihelpr,
     +               idrec1r,idrec2r,innhrf)
      end if
c
      if(ip.ne.1) then
        if(messag.ne.-1) then
          write(6,1020) ifil,ip,ierr,in(1),in(2),(in(i),i=9,14)
 1020     format(' felt ikke funnet.  file:',i3,'  ip=',i4,'  ierr=',
     +           3z10,/,'   innh.f.(1,2,9-14): ',8i7)
        end if
        if(ierr(1).eq.-1 .and. ierr(2).eq.-1 .and. ierr(3).eq.-1) then
          write(6,*) '**flt2pos** har ikke plass til aa lese feltet'
          write(6,*) '**       "lfelt" for liten.  "lfelt" = ',lfelt
        end if
        if(istop2.eq.1) goto 900
        if(istop2.eq.2) goto 100
        if(istop2.eq.3) goto 110
        goto 120
      end if
c
      if(nlev.eq.1) then
        if(ih1.eq.-32767) then
          mhlp=min(npar,maxhlp)
          do j=1,mhlp
            if(ihelp(1,j).eq.-32767) then
              do i=1,6
                ihelp(i,j)=ihelp(i,nhlp)
              end do
            end if
          end do
        end if
      end if
c
      idsave(1,np,nl,nt)=1
      idsave(2,np,nl,nt)=ifelt(20)
      idsave(3,np,nl,nt)=ifelt(19)
      idsave(4,np,nl,nt)=ifelt( 8)
      iskal=ifelt(20)
      skal=10.**iskal
      nx=ifelt(10)
      ny=ifelt(11)
      do i=1,nx*ny
        if(ifelt(20+i).ne.-32767) then
          felt(i,1)=ifelt(20+i)*skal
        else
          felt(i,1)=udef
        end if
      end do
c
      itime(1,nt)=ifelt(12)
      itime(2,nt)=ifelt(13)/100
      itime(3,nt)=ifelt(13)-ifelt(13)/100*100
      itime(4,nt)=ifelt(14)/100
      itime(5,nt)=ifelt(4)
c
      if(ingrid.eq.0) then
        newgrd=1
        ingrid=1
      elseif(idgrid(1).eq.ifelt(10) .and. idgrid(2).eq.ifelt(11) .and.
     +       idgrid(3).eq.ifelt(15) .and. idgrid(4).eq.ifelt(16) .and.
     +       idgrid(5).eq.ifelt(17) .and. idgrid(6).eq.ifelt(18)) then
        newgrd=0
      else
        newgrd=1
        write(6,*) '**flt2pos**  advarsel: skifter "grid"  ***'
        write(6,*) '**           produsent,grid: ',in(1),in(2)
        write(6,*) '**           ident(10,11,15-18):'
        write(6,*) '**      fra: ',(idgrid(i),i=1,6)
        write(6,*) '**      til: ',(ifelt(i),i=10,11),
     +                             (ifelt(i),i=15,18)
      end if
c
      if(newgrd.eq.1) then
c
        idgrid(1)=ifelt(10)
        idgrid(2)=ifelt(11)
        idgrid(3)=ifelt(15)
        idgrid(4)=ifelt(16)
        idgrid(5)=ifelt(17)
        idgrid(6)=ifelt(18)
        kontrl=0
        icall1=2
        icall2=2
        lexten=.false.
        mvx=0
        mvy=0
c
c..get grid parameters from field identification
        call gridpar(+1,20+lfelt,ifelt,igtype,nx,ny,grida,ierror)
ccc     if(ierror.ne.0) igtype=0
        if(ierror.ne.0) write(6,*) 'FLT2POS: GRIDPAR ERROR: ',ierror
        if(ierror.ne.0) return
c
c print definition of grids
        if(messag.gt.0) then
          write(*,*) 'input grid = ',ifelt(2)
          write(*,*) 'input grid type = ',igtype
          write(*,'(1x,''nx,ny = '',2i10)') nx,ny
          write(*,'(1x,''grid = '',6f10.2)') (grida(i),i=1,6)
        end if
c
        do j=1,6
          grid(j)=grida(j)
        end do
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
                 write(*,*) 'flt2pos: WARNING: not able to extend field'
                 write(*,*) '   extended nx*ny = ',(nx+6)*(ny+6)
                 write(*,*) '      (max) lfelt = ',lfelt
              end if
           end if
        end if
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
      npardir=0
      if (ipar(5,np).eq.3 .or. ipar(5,np).eq.4) then
	 if (inter.eq.4 .or. inter.eq.5) then
	    npardir=-ipar(5,np)
	 else
	    npardir=ipar(5,np)
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
      end if
c
c possibly check output min and max values
c possibly enable check in interpolation routine
      intopt=0
      if (inter.eq.3) then
         if (ipar(6,np).eq.1 .and. parlim(1,np).eq.0.0) intopt=1
         if (ipar(7,np).eq.1 .and. parlim(2,np).eq.0.0) intopt=1
      elseif (npardir.lt.0) then
c only 'nearest gridpoint' interpolation (no u,v rotation needed)
	 if (ipar(5,np).eq.4) then
            iscale=ifelt(20)
            ddmin = 10.**iscale
	    ddturn= 180.
	    do j = 1,nxa*nya
	       if (felt(j,1).ne.udef) then
	          dd = felt(j,1) + ddturn
	          if (dd.le.  0.) dd=dd+360.
	          if (dd.gt.360.) dd=dd-360.
c..assure nonzero value when scaling to integer*2
	          dd = max(dd,ddmin)
	          felt(j,1) = dd
	       end if
	    end do
	 end if
         npardir=0
         intopt=1
      end if
c
c interpolate to positions
      call grd2pos(icall1,igtype,grida,mvx,mvy,
     +             felt(1,1),nxa,nya,
     +             npos,pos(1,1),geopos(1,1),gxypos(1,1),
     +             dat(1,np,nl,nt),jinter(1,1),rinter(1,1),ndimri,
     +             udef,inter,intopt,iundef,ierror)
      if(ierror.ne.0) write(6,*) 'FLT2POS: GRD2POS ERROR: ',ierror
      if(ierror.ne.0) return
      icall1=3
c
      if (ipar(5,np).eq.2 .and. np.gt.1) then
c turn velocity components
         call grv2pos(icall2,igtype,grida,
     +                dat(1,np-1,nl,nt),dat(1,np,nl,nt),npos,
     +                vturn(1,1),geopos(1,1),udef,ierror)
         if(ierror.ne.0) write(6,*) 'FLT2POS: GRV2POS ERROR: ',ierror
         if(ierror.ne.0) return
         icall2=3
      end if
c
c true north direction parameter (dd)
      if (npardir.gt.0) then
c interpolate v component
         do n = 1,npos
            dwork(n)=dat(n,np,nl,nt)
         end do
         call grd2pos(icall1,igtype,grida,mvx,mvy,
     +                felt(1,2),nxa,nya,
     +                npos,pos(1,1),geopos(1,1),gxypos(1,1),
     +                dat(1,np,nl,nt),jinter(1,1),rinter(1,1),ndimri,
     +                udef,inter,intopt,iundef,ierror)
         if(ierror.ne.0) write(6,*) 'FLT2POS: GRD2POS ERROR: ',ierror
         if(ierror.ne.0) return
         zdeg18 = 180./(2.0*asin(1.0))
         iscale=ifelt(20)
         ddmin = 10.**iscale
         ddturn = 0.
         if (ipar(5,np).eq.4) ddturn=180.
         do n = 1,npos
            if (dwork(n).ne.udef) then
               dd = zdeg18*atan2(dwork(n),dat(n,np,nl,nt)) + ddturn
               if (dd.lt.  0.) dd=dd+360.
               if (dd.gt.360.) dd=dd-360.
c..assure nonzero value when (possibly later) scaling to integer*2
               dd = max(dd,ddmin)
               dat(n,np,nl,nt) = dd
            else
               dat(n,np,nl,nt) = udef
            end if
         end do
      end if
c
      if (ipar(6,np).eq.1) then
c check min value
         do n = 1,npos
            if (dat(n,np,nl,nt).ne.udef)
     +          dat(n,np,nl,nt)=max(dat(n,np,nl,nt),parlim(1,np))
         end do
      end if
      if (ipar(7,np).eq.1) then
c check max value
         do n = 1,npos
            if (dat(n,np,nl,nt).ne.udef)
     +          dat(n,np,nl,nt)=min(dat(n,np,nl,nt),parlim(2,np))
         end do
      end if
c
      if(istop1.eq.1 .and. kontrl.eq.0) then
        k=0
        do n=1,npos
          if(jinter(n,3).eq.0.) then
              write(6,*) '**flt2pos** posisjon nr. ',n,' er ikke o.k.'
            end if
          k=k+1
        end do
        if(k.gt.0) then
          write(6,*) '**flt2pos** posisjoner ikke o.k.'
          write(6,*) '**          antall posisjoner totalt: ',npos
          write(6,*) '**          antall ikke o.k:          ',k
          goto 900
        end if
        kontrl=1
      end if
c
      ltp=nt
  120 continue
      ltl=nt
  110 continue
      ltt=nt
  100 continue
c
  900 ltim=ltp
      if(istop2.eq.1) ltim=ltt
      if(istop2.eq.2) ltim=ltt
      if(istop2.eq.3) ltim=ltl
c
      return
      end

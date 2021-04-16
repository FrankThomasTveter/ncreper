      subroutine asr4i2(mode,ndat,dat,idat,iundef,undef,iscale,nundef)
c
c  NAME:
c     asr4i2
c
c  PURPOSE:
c     Automatic scaling from real*4 to integer*2
c     and unscaling from integer*2 to real*4.
c     Find and return the best scaling exponent to be used.
c     Possibly scale or unscale the data. Undefined values are handled.
c
c  SYNOPSIS:
c     subroutine asr4i2(mode,ndat,dat,idat,iundef,undef,iscale,nundef)
c     integer   mode,ndat,iundef,iscale,nundef
c     integer*2 idat(ndat)
c     real*4    dat(ndat),undef
c
c  INPUT:
c     mode   - 0 = return scaling exponent
c	       1 = return scaling exponent and scaled data in array dat
c		   (undefined values set to -32767.0)
c	       2 = return scaling exponent and scaled data in array idat
c		   (undefined values set to -32767)
c	      -1 = use scaling exponent and unscale data in array dat
c	      -2 = use scaling exponent and unscale data in array idat
c		   with output in array dat.
c     ndat   - length of data arrays
c     iundef - 0=no undefined values  1=undefined values may exist
c     undef  - value in an undefined data array element
c
c  INPUT/OUTPUT:
c     dat    - real*4    data array, input  in mode 0,1,2,-1
c				     output in mode 1,-1
c     idat   - integer*2 data array, input  in mode -2
c				     output in mode 2
c     iscale - scaling exponent, input  in mode -1,-2
c				 output in mode 0,1,2
c
c  OUTPUT:
c     nundef - no. of undefine values found
c
c  NOTES:
c     When scaling from real*4 to integer*2:
c	 scale = 10.**(-iscale)
c	 idat(n) = nint(dat(n)*scale)
c     When unscaling from integer*2 to real*4:
c	 scale = 10.**iscale
c	 dat(n) = float(idat(n))*scale
c     If all values are undefined or zero, then this routine
c     returns iscale = 0 .
c
c-----------------------------------------------------------------------
c  DNMI/FoU  27.04.1995  Anstein Foss
c  DNMI/FoU  24.07.1995  Anstein Foss
c-----------------------------------------------------------------------
c
      implicit none
c
      integer   mode,ndat,iundef,iscale,nundef
      integer*2 idat(ndat)
      real*4    dat(ndat),undef,scale
c
      integer   ndef,nudef,i,ifmax
      integer*2 iudef
      real*4    rudef,tudef,fmax
c
      iudef=-32767
      rudef=-32767.0
      tudef=-32766.5
      nudef=0
c
      if(mode.eq.0 .or. mode.eq.1 .or. mode.eq.2) then
c
        fmax=0.
        ndef=0
        if(iundef.eq.0) then
          do i=1,ndat
            fmax=max(fmax,abs(dat(i)))
          end do
	  ndef=ndat
        else
          do i=1,ndat
            if(dat(i).ne.undef) then
	      fmax=max(fmax,abs(dat(i)))
	      ndef=ndef+1
	    end if
          end do
	  nudef=ndat-ndef
        end if
        if(ndef.gt.0 .and. fmax.gt.0.) then
	  iscale=log10(fmax)-4.
	  ifmax=nint(fmax*10.**(-iscale))
	  if(ifmax.lt.3278) then
	    iscale=iscale-1
	    ifmax=nint(fmax*10.**(-iscale))
	  end if
	  if(ifmax.gt.32766) iscale=iscale+1
	  iscale=max(iscale,-30)
        else
	  iscale=0
        end if
c
      end if
c
      if(mode.eq.1) then
c
c..scale: real*4 -> real*4
	scale = 10.**(-iscale)
	if(nudef.eq.0) then
	  do i=1,ndat
	    dat(i)=dat(i)*scale
	  end do
	else
	  do i=1,ndat
	    if(dat(i).ne.undef) then
	      dat(i)=dat(i)*scale
	    else
	      dat(i)=rudef
	    end if
	  end do
	end if
c
      elseif(mode.eq.2) then
c
c..scale: real*4 -> integer*2
	scale = 10.**(-iscale)
	if(nudef.eq.0) then
	  do i=1,ndat
	    idat(i)=nint(dat(i)*scale)
	  end do
	else
	  do i=1,ndat
	    if(dat(i).ne.undef) then
	      idat(i)=nint(dat(i)*scale)
	    else
	      idat(i)=iudef
	    end if
	  end do
	end if
c
      elseif(mode.eq.-1) then
c
c..unscale: real*4 -> real*4
	scale = 10.**iscale
	if(iundef.eq.0) then
	  do i=1,ndat
	    dat(i)=dat(i)*scale
	  end do
	  nudef=0
	else
	  do i=1,ndat
	    if(dat(i).gt.tudef) then
	      dat(i)=dat(i)*scale
	    else
	      dat(i)=undef
	      nudef=nudef+1
	    end if
	  end do
	end if
c
      elseif(mode.eq.-2) then
c
c..unscale: integer*2 -> real*4
	scale = 10.**iscale
	if(iundef.eq.0) then
	  do i=1,ndat
	    dat(i)=float(idat(i))*scale
	  end do
	else
	  do i=1,ndat
	    if(idat(i).eq.iudef) then
	      dat(i)=float(idat(i))*scale
	    else
	      dat(i)=undef
	      nudef=nudef+1
	    end if
	  end do
	end if
c
      end if
c
      nundef=nudef
c
      return
      end

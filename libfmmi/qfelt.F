      subroutine qfelt(iunit,ireq,iexist,nin,in,ifound,nfound,
     +                 iend,ierror,ioerr)
c
c  NAME:
c     qfelt
c
c  PURPOSE:
c     Inquire felt file contents
c
c  SYNOPSIS:
c     subroutine qfelt(iunit,ireq,iexist,nin,in,ifound,nfound,
c    +                 iend,ierror,ioerr)
c     integer    iunit,ireq,iexist,nin,nfound,iend,ierror,ioerr
c     integer    ifound(nin)
c     integer*2  in(16,nin)
c
c  INPUT:
c     iunit         - felt file unit,
c                     if qfelt is called with several files, using
c                     the same file unit no., reset qfelt by a call
c                     with iunit=0.
c     ireq          - request type:
c                       1 = return all occurrences of in(1:16,1),
c                           in(i,1)=-32767 means any value
c                       2 = return different occurrences of one of the
c                           elements in in(1:16,1), the selected
c                           element is the first with value equal
c                           -32766 , while -32767 means any value
c                       3 = return first occurrence of in(1:16,1:nin),
c                           where in(i,n)=-32767 means any value
c                       0 = reset when stopping requests of type
c                           1 or 2 before end is reached and the
c                           next request is identical, but should
c                           start from top of the file index
c                       11,12,13: as 1,2,3 above, but see output ifound
c     iexist        -  control of data existence:
c                      0 = return all innh.fort. found
c                      1 = return only innh.fort. with existing data
c                          (for request types 1 and 2)
c                      2 = return only innh.fort. without existing data
c                          (for request types 1 and 2)
c     nin           - max. no. of innh.fort. to be returned
c     in(16,nin)    - see requirements for each request type (ireq)
c                     and the notes below
c
c  OUTPUT:
c     in(16,...)    - the innh.fort. found
c     ifound(...)   - length of the fields if they exist,
c                     for ireq=1,2,3 the length of the integer*2 buffer
c                     for ireq=11,12,13 the field size (nx*ny)
c                      0 if an innh.fort. is found without data,
c                     -1 if innh.fort. is not found
c     nfound        - no. of innh.fort. found, for request type 1 and 2
c                     (for request type 3 nfound=nin unless some
c     iend          - completion status:
c                      0 = more calls to qfelt is necessary to
c                          return all requested information,
c                          for request types 1 and 2
c                      1 = all requested information is returned
c                          (if same request is repeated in next call,
c                           qfelt starts from top of the file index)
c     ierror        -  0 = no error
c                      1 = read error (file error)
c                      2 = unknown ireq or iexist specified
c                      3 = error in input in array
c     ioerr         -  read error (iostat, 0 = no error)
c
c  NOTES:
c     - qfelt always returns iend=1 and nfound=0 if ierror =/= 0.
c     - a continued call (after iend=0) is identified by unchanged
c       file (iunit), request type (ireq) and in(1:16,1) .
c     - for request type 2 does a continued call not mean that qfelt
c       returns innh.fort. elements different from the previous call(s),
c       but all element values will be returned (possibly more than
c       once)
c     - iend=1 doesn't mean there are more information, rather that
c       the 'in' array is full before the complete file index has neen
c       searched
c     - remeber to reset in(1:16,1) before continued calls.
c     - input date and time, in(3:5,n), are used, set these elements
c       to -32767 if they are unimportant
c     - input in(6:8,n) and in(16,n) are not used
c     - in(1:16,n) means in(i,n) for i= 1 to 16  (fortran90 style)
c     - for performance (i/o) reasons, nin should not be a very small
c       number for request type 1 and 2, rather N*64 (if a small number
c	means many calls to qfelt)
c
c  EXAMPLES:
c	integer    maxnin
c	parameter (maxnin=256)
c	integer    iunit,ireq,iexist,nin,nfound,iend,ierror,ioerr
c	integer    ifound(maxnin)
c	integer*2  in(16,maxnin)
c	integer    i,igrid
c
c          find first (and perhaps only) grid for producer
c
c	ireq=1
c	iexist=0
c	do i=1,16
c	  in(i,1)=-32767
c	end do
c	in(1,1)=88
c	nin=1
c	call qfelt(iunit,ireq,iexist,nin,in,ifound,nfound,
c    +             iend,ierror,ioerr)
c	igrid=in(2,1)
c
c          find all forcast lengths with data for
c	   producer 88, grid 1234 and parameter 5
c
c	ireq=2
c	iexist=1
c	ierror=0
c	iend=0
c	do while (iend.eq.0 .and. ierror.eq.0)
c	  do i=1,16
c	    in(i,1)=-32767
c	  end do
c	  in( 1,1)=    88
c	  in( 2,1)=  1234
c	  in(12,1)=     5
c	  in(10,1)=-32766
c	  nin=maxnin
c	  call qfelt(iunit,ireq,iexist,nin,in,ifound,nfound,
c    +               iend,ierror,ioerr)
c	  do n=1,nfound
c	    <record in(10,n) in a table>
c	  end do
c	end do
c
c          get all 'innh.fort' with existing field data:
c
c	ireq=1
c	iexist=1
c	ierror=0
c	iend=0
c	do while (iend.eq.0 .and. ierror.eq.0)
c	  do i=1,16
c	    in(i,1)=-32767
c	  end do
c	  nin=maxnin
c	  call qfelt(iunit,ireq,iexist,nin,in,ifound,nfound,
c    +               iend,ierror,ioerr)
c	  do n=1,nfound
c	    <use in(i,n) and/or ifound(n)>
c	  end do
c	end do
c
c-----------------------------------------------------------------------
c  DNMI/FoU  30.09.1993  Anstein Foss
c  DNMI/FoU  29.05.1996  Anstein Foss ... iexist=2
c  DNMI/FoU  10.11.2000  Anstein Foss ... automatic byte swap
c  DNMI/FoU  01.02.2001  Anstein Foss ... ireq= 11,12,13 (ifound()=nx*ny)
c-----------------------------------------------------------------------
c
      parameter (lrec=1024,nkeys=16)
      parameter (ninh=lrec/nkeys)
c
c..input/output
      integer   iunit,ireq,nin,nfound,iend,ierror
      integer   ifound(nin)
      integer*2 in(16,nin)
c
c..local
      integer   initi
      integer   igtype, lgeom, irequest
      integer   lunit, lstreq, lstext, lstelm, lstrec, lstinh
      integer*2 minadr,maxadr
      integer*2 inreq(16),idrec1(32),idrec2(256)
      integer*2 inh(16,ninh)
      logical   swapfile,swap
c
      data initi/0/
c
      save lunit, lstreq, lstext, lstelm, lstrec, lstinh
      save inreq, idrec1, idrec2
c
      if(initi.eq.0) then
c..initialization of control variables (for continued calls)
        lunit =0
        lstreq=0
        lstext=0
        lstelm=0
        lstrec=0
        lstinh=0
        do i=1,nkeys
          inreq(i)=-32767
        end do
        initi =1
      end if
c
      nfound=0
      nf=0
      ierror=0
      ioerr=0
c
      irequest=ireq
      if (ireq.ge.11 .and. ireq.le.13) irequest=ireq-10
c
      if(iunit.le.0) then
        lunit =0
	iend=1
        goto 950
      end if
c
      if(irequest.eq.0) then
	iend=1
        goto 950
      end if
c
      if(irequest.lt.1 .or. irequest.gt.3 .or.
     +   iexist.lt.0   .or. iexist.gt.2) then
        ierror=2
        goto 950
      end if
c
c..check if this a continued call,
c..if any elements are specified,
c..find the selected element for request type 2
c..and save the request array to be used in search (type 1 and 2)
c
      icont=1
      if(iunit .ne.lunit ) icont=0
      if(iexist.ne.lstext) icont=0
      if(ireq  .ne.lstreq) icont=0
      ieldef=0
      ielem=0
      if(irequest.ne.3) then
        do i=1,5
          if(in(i,1).ne.inreq(i)) icont=0
          if(in(i,1).ne.-32767) ieldef=ieldef+1
          if(in(i,1).eq.-32766 .and. ielem.eq.0) ielem=i
          inreq(i)=in(i,1)
        end do
        do i=9,15
          if(in(i,1).ne.inreq(i)) icont=0
          if(in(i,1).ne.-32767) ieldef=ieldef+1
          if(in(i,1).eq.-32766 .and. ielem.eq.0) ielem=i
          inreq(i)=in(i,1)
        end do
      end if
      if(irequest.eq.2 .and. ielem.ne.lstelm) icont=0
c
      if(irequest.eq.2 .and. ielem.eq.0) then
        ierror=3
        goto 950
      end if
c
      if(irequest.eq.2) then
	inreq(ielem)=-32767
	ieldef=ieldef-1
      end if
c
      irqtyp=irequest
c..for request type 1 and 2 when no elements are specified
c..(to avoid too much unnecessary testing, and execute faster)
      if(irequest.eq.1 .and. ieldef.eq.0) irqtyp=-1
      if(irequest.eq.2 .and. ieldef.eq.0) irqtyp=-2
c
      ni=1
      if(irequest.eq.3) ni=nin
      minprod=+32767
      maxprod=-32767
      iany=1
      if(irequest.ne.2) then
        iany=0
        do n=1,ni
          if(in(1,n).ne.-32767) then
            if(minprod.gt.in(1,n)) minprod=in(1,n)
            if(maxprod.lt.in(1,n)) maxprod=in(1,n)
          else
            iany=1
          end if
        end do
      end if
      if(iany.eq.0) then
        if(minprod.lt.1 .or. maxprod.gt.99) then
          ierror=3
          goto 950
        end if
      else
        minprod= 1
        maxprod=99
      end if
c
      swap= swapfile(iunit)
c
      if(iunit.ne.lunit) then
        lunit=0
        lstreq=0
        lstext=0
        lstelm=0
        lstrec=0
        lstinh=0
        read(iunit,rec=1,iostat=ios,err=910) idrec1
        read(iunit,rec=2,iostat=ios,err=910) idrec2
	if (swap) call bswap2( 32,idrec1)
	if (swap) call bswap2(256,idrec2)
        lunit=iunit
      end if
c
      if(icont.eq.1) then
        irec1=lstrec
        ninh1=lstinh+1
        if(ninh1.gt.ninh) then
          ninh1=1
          irec1=irec1+1
        end if
      elseif(minprod.eq.1) then
        irec1=3
        ninh1=1
      else
        n=0
        do iprod=minprod,maxprod
          if(n.eq.0 .and. idrec2(iprod*2-1).gt.0) n=iprod
        end do
        if(n.eq.0) then
          iend=1
          goto 950
        end if
        irec1=idrec2(n*2-1)
        ninh1=(idrec2(n*2)+15/16)
      end if
c
      n=0
      if(maxprod.lt.99) then
        do iprod=maxprod+1,99
          if(n.eq.0 .and. idrec2(iprod*2-1).gt.0) n=iprod
        end do
      end if
      if(n.eq.0) then
        irec2=idrec1(9)
        ninh2=ninh
      else
        irec2=idrec2(n*2-1)
        ninh2=(idrec2(n*2)+15/16)-1
        if(ninh2.lt.1) then
          irec2=irec2-1
          ninh2=ninh
        end if
      end if
c
      if(irequest.eq.3) then
	do n=1,nin
	  in( 6,n)=-32767
	  in( 7,n)=-32767
	  in( 8,n)=-32767
	  in(16,n)=-32767
	  ifound(n)=-1
	end do
      end if
c
      minadr=-32767
      maxadr=+32767
      if(iexist.eq.1) minadr=1
      if(iexist.eq.2) maxadr=0
c
      irec=irec1-1
c
      do while (irec.lt.irec2 .and. nf.lt.nin)
c
        irec=irec+1
        read(iunit,rec=irec,iostat=ios,err=910) inh
	if (swap) call bswap2(lrec,inh)
        ni =0
        ni2=ninh
        if(irec.eq.irec1) ni =ninh1-1
        if(irec.eq.irec2) then
	  if(ninh2.eq.ninh) then
	    ninh2=0
	    do n=1,ninh
	      if(inh(1,n).ge.1 .and. inh(1,n).le.99) ninh2=n
	    end do
	  end if
	  ni2=ninh2
	end if
c
        if(irqtyp.eq.-1) then
c..request type 1 with no elements specified
	  do while (ni.lt.ni2 .and. nf.lt.nin)
	    ni=ni+1
            if(inh(6,ni).ge.minadr .and. inh(6,ni).le.maxadr) then
              nf=nf+1
              do i=1,nkeys
                in(i,nf)=inh(i,ni)
              end do
              ifound(nf)=0
	    end if
	  end do
c
        elseif(irqtyp.eq.1) then
c..request type 1 (with one or more elements specified)
	  do while (ni.lt.ni2 .and. nf.lt.nin)
	    ni=ni+1
            if(inh(6,ni).ge.minadr .and. inh(6,ni).le.maxadr) then
	      ihit=1
	      do i=1,nkeys
		if(inreq(i).ne.-32767 .and.
     +		   inreq(i).ne.inh(i,ni)) ihit=0
	      end do
	      if(ihit.eq.1) then
                nf=nf+1
                do i=1,nkeys
                  in(i,nf)=inh(i,ni)
                end do
                ifound(nf)=0
	      end if
	    end if
	  end do
c
        elseif(irqtyp.eq.-2) then
c..request type 2 with no elements specified
	  do while (ni.lt.ni2 .and. nf.lt.nin)
	    ni=ni+1
            if(inh(6,ni).ge.minadr .and. inh(6,ni).le.maxadr) then
              iequal=0
	      n=0
              do while (iequal.eq.0 .and. n.lt.nf)
	        n=n+1
                if(inh(ielem,ni).eq.in(ielem,n)) iequal=1
              end do
              if(iequal.eq.0) then
                nf=nf+1
                do i=1,nkeys
                  in(i,nf)=inh(i,ni)
                end do
                ifound(nf)=0
              end if
	    end if
	  end do
c
        elseif(irqtyp.eq.2) then
c..request type 2 (with one or more elements specified)
	  do while (ni.lt.ni2 .and. nf.lt.nin)
	    ni=ni+1
            if(inh(6,ni).ge.minadr .and. inh(6,ni).le.maxadr) then
	      ihit=1
	      do i=1,nkeys
		if(inreq(i).ne.-32767 .and.
     +		   inreq(i).ne.inh(i,ni)) ihit=0
	      end do
	      if(ihit.eq.1) then
                iequal=0
	        n=0
                do while (iequal.eq.0 .and. n.lt.nf)
	          n=n+1
                  if(inh(ielem,ni).eq.in(ielem,n)) iequal=1
                end do
                if(iequal.eq.0) then
                  nf=nf+1
                  do i=1,nkeys
                    in(i,nf)=inh(i,ni)
                  end do
                  ifound(nf)=0
                end if
	      end if
	    end if
	  end do
c
        elseif(irqtyp.eq.3) then
c..request type 3
	  do while (ni.lt.ni2)
	    ni=ni+1
	    ihit=0
	    n=0
	    do while (n.lt.nin .and. ihit.eq.0)
	      n=n+1
	      if(ifound(n).eq.-1) then
	        ihit=1
	        do i=1,nkeys
		  if(in(i,n).ne.-32767 .and.
     +		     in(i,n).ne.inh(i,ni)) ihit=0
	        end do
	      end if
	    end do
	    if(ihit.eq.1) then
	      do i=1,nkeys
		in(i,n)=inh(i,ni)
	      end do
	      ifound(n)=0
	      nf=nf+1
	    end if
	  end do
	end if
c
      end do
c
      lstrec=irec
      lstinh=ni
      iend=0
      if(irec.eq.irec2 .and. ni.eq.ninh2) iend=1
c
      goto 950
c
  910 ierror=1
      ioerr=ios
c
  950 if(ierror.ne.0) then
	nf=0
	iend=1
      elseif(irequest.eq.3) then
	nf=nin
	iend=1
      end if
      nfound=nf
c
c..length of fields
      if (ireq.eq.irequest) then
c..return length of integer*2 buffer needed to read the field
        do n=1,nf
	  if(in(6,n).gt.0) then
            l1=in( 8,n)
            l2=in(16,n)/100
            ifound(n)=l1+l2*32767
          end if
        end do
      else
c..return length of the field (nx*ny)
        do n=1,nf
	  if(in(6,n).gt.0) then
            l1=in( 8,n)
            l2=in(16,n)/100
	    lgeom=0
	    igtype=in(15,n)
	    if (igtype.ge.1000) lgeom=igtype-(igtype/1000)*1000
	    ifound(n)=l1+l2*32767 - 20 - lgeom
          end if
        end do
      end if
c
      if(iend.eq.0) then
        lstreq=ireq
        lstext=iexist
        lstelm=ielem
      else
        lstreq=0
        lstext=0
        lstelm=0
        lstrec=0
        lstinh=0
      end if
c
      return
      end

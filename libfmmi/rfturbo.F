      subroutine rfturbo(iunit,ip,in,ifelt,limit,ierr,ihelp,
     +                   idrec1,idrec2,innh)
c
c  NAME:
c     rfturbo
c
c  WARNING:
c     Not recommended for 'user' programs.
c
c  PURPOSE:
c     Read field (identification and field data) from Felt file.
c     Turbo version of the standard routine rfelt.
c     For standard utility programs.
c     Faster than rfelt due to less i/o of the index part of the file.
c     Possibly unsafe due to 'transfers' between calls (transfer
c     buffers used to allow reading from more than one file).
c
c  SYNOPSIS:
c     subroutine rfturbo(iunit,ip,in,ifelt,limit,ierr,ihelp,
c    +                   idrec1,idrec2,innh)
c     integer   iunit,ip,limit
c     integer   ierr(3),ihelp(6)
c     integer*2 in(16),ifelt(limit)
c     integer*2 idrec1(1024),idrec2(1024),innh(16,64)
c
c  INPUT:
c     iunit  - Felt file unit no.
c     in     - parameters describing the field ('innholdsfortegnelse'),
c              the following is always required:
c                in( 1) = producer no.
c                in( 2) = grid no.
c                in( 9) = data type
c                in(10) = forecast length (hours)
c                in(11) = vertical coordinate
c                in(12) = parameter no.
c                in(13) = level 1
c                in(14) = level 2
c              if ierr(2)=999 (se below) then also:
c                in( 3) = year
c                in( 4) = month*100+day
c                in( 5) = hour*100+minute
c	       the value -32767 is the same as an unspecified value,
c	       if one or more -32767 then the first existing field with
c	       matching non -32767 elements (if any) will be returned
c	       in(1)=-32767 is not allowed.
c     limit  - length of the ifelt array (maximum length of
c              identification and field data, some gridtypes also
c	       with extra identification after the field data)
c     ierr   - option.
c                check date/time:
c                  ierr(2) =/= 999 : not done
c                  ierr(2)  =  999 : done
c     ihelp  - help array for fast data access (if fields are written
c              in the same sequence as the 'innholdsfortegnelse').
c              initialize ihelp(1:6) to 0 before first call, and set
c              ihelp(1)=0 when a new file is opened.
c     idrec1 - record 1 of the felt file (transfer buffer between calls)
c     idrec2 - record 2 of the felt file (transfer buffer between calls)
c     innh   - one index record          (transfer buffer between calls)
c
c  OUTPUT:
c     ip     - exit status:
c                ip= 1   no error, field read.
c                ip= 0   grid data not read, no data exists
c                        ('innholdsfortegnelse' was found)
c                ip=-1   grid data not read, read error
c                        (i/o error returned in ierr(1:3))
c                ip=-2   'innholdsfortegnelse' not found, therefore
c                        field not read.
c                ip=-3   grid data not read because limit is too small
c                ip=-4   grid data not read, the producer no. is
c                        not defined in the file
c                ip=-5   bad grid data read, corrupt pointer from
c                        index to data
c     in     - the 'innholdsfortegnelse' read (with data pointer etc.)
c     ifelt  - identification (first 20 elements) and field data if ip=1.
c     ierr   - i/o error codes ('iostat'):
c                ierr(1) - errorcode on reading record 1
c                ierr(2) - errorcode on reading 'innholdsfortegnelse'
c                ierr(3) - errorcode on reading field
c     ihelp  - help array, store this until the next call
c     idrec1 - record 1 of the felt file (transfer buffer between calls)
c     idrec2 - record 2 of the felt file (transfer buffer between calls)
c     innh   - one index record          (transfer buffer between calls)
c
c  NOTES:
c     +  Remember to set ierr(2)=999 when reading fields from
c        Archive or Cyclic Archive Felt files
c     +  The ihelp array: ihelp(1) - file unit no.
c                         ihelp(2) - producer no.
c                         ihelp(3) - last read record in 'innh.fort.'
c                         ihelp(4) - first record in 'innh.fort.' for
c                                    producer ihelp(2)
c                         ihelp(5) - pointer to the first 'innh.fort.'
c                                    in the first record
c                         ihelp(6) - last record in 'innh.fort.'
c
c-----------------------------------------------------------------------
c  D N M I
c      written by Rebecca Rudsar    12.05.1982 ... rfelt
c      updated by Rebecca Rudsar     9.08.1982
c                 Anstein Foss      10.05.1984
c                 Anstein Foss      12.03.1986
c                 Anstein Foss      15.10.1992
c                 Anstein Foss      02.12.1993 ... large field & file
c                 Anstein Foss      11.03.1994 ... rfturbo
c                 Anstein Foss      10.05.1994
c                 Anstein Foss      30.04.1995 ... in(i)=-32767
c                 Anstein Foss      09.06.1995 ... ifelt(ldata+1)=-32767
c                 Anstein Foss      10.11.2000 ... automatic byte swap
c                 Anstein Foss      03.03.2003 ... check and ip=-5 return
c-----------------------------------------------------------------------
c
c..lrec   - record length in integer*2 (16 bit) words
c..nin    - no. of innholdsfortegnelser per record
c
      parameter (lrec=1024)
      parameter (nin=lrec/16)
c
c..input/output
      integer   iunit,ip,limit
      integer   ierr(3),ihelp(6)
      integer*2 in(16),ifelt(limit)
      integer*2 innh(16,nin),idrec1(lrec),idrec2(lrec)
c
c..local
      integer   nrb(2,2)
      integer*2 ibuffr(lrec)
      logical   swapfile,swap
c
c..if ierr(2) = 999 then test on date in 'innh'.
      ntest1=2
      if(ierr(2).eq.999) ntest1=5
c
c..unspecified elements
      notex=0
      nuspec=0
      do i=2,ntest1
        if(in(i).eq.-32767) nuspec=nuspec+1
      end do
      do i=9,14
        if(in(i).eq.-32767) nuspec=nuspec+1
      end do
c
      ierr(1)=0
      ierr(2)=0
      ierr(3)=0
c
      in( 6)=-32767
      in( 7)=-32767
      in( 8)=-32767
      in(16)=-32767
c
      if(in(1).lt.1 .or. in(1).gt.99) then
        ip=-4
        goto 1010
      end if
c
      swap= swapfile(iunit)
c
      lastreci=ihelp(3)
c
      newhelp=0
c
   10 continue
c
      if(ihelp(1).ne.iunit .or. ihelp(2).ne.in(1)
     *                     .or.  newhelp.eq.1) then
c
        if(ihelp(1).ne.iunit .or. newhelp.eq.1) then
          read(iunit,rec=1,err=1001,iostat=ios) idrec1
          read(iunit,rec=2,err=1001,iostat=ios) idrec2
	  if (swap) call bswap2(lrec,idrec1)
	  if (swap) call bswap2(lrec,idrec2)
        end if
c..record no. of last 'innh' record.
        ireci2=idrec1(9)
        iprod=in(1)
        ireci1=idrec2(2*iprod-1)
        if(ireci1.le.0) then
          ip=-4
          goto 1010
        end if
        nord=idrec2(2*iprod)
        ji=nord/16+1
c
c..in case we got here after search failure
        if(newhelp.eq.1) then
          if(ihelp(4).eq.ireci1 .and. ihelp(5).eq.ji
     *                          .and. ihelp(6).eq.ireci2) then
c..we have scanned the correct part of the 'innholdsfortegnelse'
	    ihelp(3)=lastreci
            ip=-2
            return
          end if
        end if
c
        ihelp(1)=iunit
        ihelp(2)=in(1)
        ihelp(3)=ireci1
        ihelp(4)=ireci1
        ihelp(5)=ji
        ihelp(6)=ireci2
        lastreci=0
c
        newhelp=1
c
      end if
c
      nloop=2
      if(ihelp(3).lt.ihelp(4) .or.
     *   ihelp(3).gt.ihelp(6)) ihelp(3)=ihelp(4)
      if(ihelp(3).eq.ihelp(4)) nloop=1
      nrb(1,1)=ihelp(3)
      nrb(2,1)=ihelp(6)
      nrb(1,2)=ihelp(4)
      nrb(2,2)=ihelp(3)-1
c
      do 20 iloop=1,nloop
c
        ireci1=nrb(1,iloop)
        ireci2=nrb(2,iloop)
        ji=1
        if(ireci1.eq.ihelp(4)) ji=ihelp(5)
c
        do ireci=ireci1,ireci2
c
          if(ireci.ne.lastreci) then
            read(iunit,rec=ireci,err=1002,iostat=ios) innh
	    if (swap) call bswap2(lrec,innh)
            lastreci=ireci
          end if
c
c..search for correct 'innholdsfortegnelse'
c
	  if(nuspec.eq.0) then
            do j=ji,nin
              k=1
              do i=1,ntest1
                if(innh(i,j).ne.in(i)) k=0
              end do
              do i=9,14
                if(innh(i,j).ne.in(i)) k=0
              end do
              if(k.eq.1) goto 50
            end do
	  else
            do j=ji,nin
              k=1
              do i=1,ntest1
                if(innh(i,j).ne.in(i) .and. in(i).ne.-32767) k=0
              end do
              do i=9,14
                if(innh(i,j).ne.in(i) .and. in(i).ne.-32767) k=0
              end do
              if(k.eq.1 .and. innh(6,j).gt.0
     +			.and. innh(8,j).gt.0) goto 50
              if(k.eq.1) notex=1
            end do
	  end if
c
          if(innh(1,nin).ne.in(1)) goto 20
          ji=1
c
        end do
c
   20 continue
c
c..safety first (ihelp(1) not reset and another file with same unit)
      if(newhelp.eq.0) then
        newhelp=1
        goto 10
      end if
c
      ip=-2
      if(notex.eq.1) ip=0
      ihelp(3)=lastreci
      return
c
c..here we have found the correct grid data.
c
   50 ihelp(3)=ireci
c
      do i=1,16
        in(i)=innh(i,j)
      end do
c
      if(in(6).lt.1 .or. in(8).lt.1) then
        ip=0
        return
      end if
c
c..handle 'old' files well.........................................
      if(in(16).lt.0) in(16)=0
      if(in( 7).lt.1 .or. in(7).gt.lrec) in(7)=1
c..................................................................
c
      irecd1=in(6)
      iword =in(7)
      ldata1=in(8)
      ldata2=in(16)/100
      irecd2=in(16)-ldata2*100
      irec1=irecd1+irecd2*32767
      ldata=ldata1+ldata2*32767
c
      if(ldata.gt.limit) goto 1004
c
      nrec=(iword-1+ldata+lrec-1)/lrec
      irec2=irec1+nrec-1
      l2=0
      if(iword.gt.1) then
        lb=iword-1
        l1=l2+1
        l2=l2+lrec-lb
        if(l2.gt.ldata) l2=ldata
        read(iunit,rec=irec1,err=1003,iostat=ios)
     *               (ibuffr(i),i=1,lb),(ifelt(i),i=l1,l2)
        irec1=irec1+1
      end if
      do irec=irec1,irec2
        l1=l2+1
        l2=l2+lrec
        if(l2.gt.ldata) l2=ldata
        read(iunit,rec=irec,err=1003,iostat=ios) (ifelt(i),i=l1,l2)
      end do
c
      if (swap) call bswap2(ldata,ifelt)
c
c..in cases (i.e. for gridtypes) that should have field identification
c..at the end, this will tell if it is missing or not
      if(ldata.lt.limit) ifelt(ldata+1)=-32767
c
      if(in( 1).ne.-32767 .and. ifelt(1).ne.in( 1)) goto 1005
      if(in( 2).ne.-32767 .and. ifelt(2).ne.in( 2)) goto 1005
      if(in( 9).ne.-32767 .and. ifelt(3).ne.in( 9)) goto 1005
      if(in(10).ne.-32767 .and. ifelt(4).ne.in(10)) goto 1005
      if(in(11).ne.-32767 .and. ifelt(5).ne.in(11)) goto 1005
      if(in(12).ne.-32767 .and. ifelt(6).ne.in(12)) goto 1005
      if(in(13).ne.-32767 .and. ifelt(7).ne.in(13)) goto 1005
      if(in(14).ne.-32767 .and. ifelt(8).ne.in(14)) goto 1005
c
      ip=1
      return
c
 1001 ip=-1
      ierr(1)=ios
      return
c
 1002 ip=-1
      ierr(2)=ios
      return
c
 1003 ip=-1
      ierr(3)=ios
      return
c
 1004 ip=-3
      return
c
 1005 ip=-5
      return
c
 1010 continue
c
      return
      end

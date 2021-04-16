      subroutine wfelt(iunit,ip,in,idrec1,idrec2,lfelt,ifelt,ierr,ihelp)
c
c  NAME:
c     wfelt
c
c  PURPOSE:
c     Write field (identification and field data) to Felt file
c
c  SYNOPSIS:
c     subroutine wfelt(iunit,ip,in,idrec1,idrec2,lfelt,ifelt,ierr,ihelp)
c     integer   iunit,ip,lfelt
c     integer   ierr(3),ihelp(6)
c     integer*2 in(16),idrec1(1024),idrec2(1024),ifelt(lfelt)
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
c              if ierr(2)=999 (se below) or the Felt file is an Archive
c              or Cyclic Archive Felt file then also:
c                in( 3) = year
c                in( 4) = month*100+day
c                in( 5) = hour*100+minute
c     idrec1 - record 1 of the felt file
c     idrec2 - record 2 of the felt file
c     lfelt  - the length (used) of the ifelt array
c              (20 + size of field)
c     ifelt  - array holding field identification and field data.
c              the 20 first word is the field identification,
c              some of these elements are copied from the in array,
c              and the following has to be specified.
c              polarstereographic grid (standard):
c                ifelt( 9) = 1, grid type
c                ifelt(10) = no. of points in x-direction
c                ifelt(11) = no. of points in y-direction
c                ifelt(15) = x-coordinate of north pole * 100
c                ifelt(16) = y-coordinate of north pole * 100
c                ifelt(17) = grid distance unit, km * 10
c                ifelt(18) = grid rotation angle (longitude parallel
c                            y-axis) in degrees (+ = east, - = west)
c              polarstereographic grid (when large distance to the pole):
c                ifelt( 9) = 1, grid type
c                ifelt(10) = no. of points in x-direction
c                ifelt(11) = no. of points in y-direction
c                ifelt(15) = x-coordinate of north pole
c                ifelt(16) = y-coordinate of north pole
c                ifelt(17) = grid distance unit, km * -10 (negative)
c                ifelt(18) = grid rotation angle (longitude parallel
c                            y-axis) in degrees (+ = east, - = west)
c              geographic grid:
c                ifelt( 9) = 2, grid type
c                ifelt(10) = no. of points in x-direction, longitude
c                ifelt(11) = no. of points in y-direction, latitude
c                ifelt(15) = origo (x=1,y=1) latitude,    degrees * 100
c                ifelt(16) = origo (x=1,y=1) longitude,   degrees * 100
c                ifelt(17) = grid distance, latitude (y), degrees * 100
c                ifelt(18) = grid distance, longitude(x), degrees * 100
c              other grid:
c                ifelt( 9) = 0, grid type
c                ifelt(10) = no. of points in x-direction
c                ifelt(11) = no. of points in y-direction
c              and for all grid types:
c                ifelt(19) = 0 (or something else if it is used elsewhere)
c                ifelt(20) = field scaling exponent
c                            (real_value = ifelt(..) * 10.**ifelt(20) )
c     ierr   - options.
c                overwrite existing data:
c                  ierr(1)  =  0 : not allowed
c                  ierr(1) =/= 0 : allowed if the sizes of the existing
c                                  and the new field are equal
c                check date/time:
c                  ierr(2) =/= 999 : not done (unless Archive Felt file)
c                  ierr(2)  =  999 : done
c     ihelp  - help array for fast data access (if fields are written
c              in the same sequence as the 'innholdsfortegnelse').
c              initialize ihelp(1:6) to 0 before first call, and set
c              ihelp(1)=0 when a new file is opened.
c
c  OUTPUT:
c     ip     - exit status:
c                ip= 1   no error, field written.
c                ip= 0   'innholdsfortegnelse' not found, therefore
c                        field not written.
c                ip=-1   grid data not written, read or write error
c                        (i/o error returned in ierr(1:3))
c                ip=-2   grid data not written because of existing data.
c                ip=-3   grid data not overwritten because data length
c                        do not agree.
c                ip=-4   grid data not written, the producer no. is
c                        not defined in the file
c                ip=-5   grid data not written because input data length
c                        less than 21 (the minimum)
c                ip=-6   too big field (lfelt too big)
c                ip=-7   maximum file size exceeded
c     in     - the 'innholdsfortegnelse' written (with data pointer etc.)
c     ierr   - i/o error codes ('iostat'):
c                ierr(1) - errorcode on reading 'innholdsfortegnelse'
c                ierr(2) - errorcode on writing 'innholdsfortegnelse'
c                ierr(3) - errorcode on writing field
c                          (or reading existing data when idrec1(14)=1)
c     ihelp  - help array, store this until the next call
c
c  NOTES:
c     +  If idrec1(14)=1 then fields (identification+field) are
c        written to the field without gaps, i.e. the data will
c        usually not start in the first word of a record, and
c        not occupy a whole number of records).
c        Set this switch with options to programs that create
c        felt files (do not set idrec1(14)=1 in user programs).
c     +  The ihelp array: ihelp(1) - file unit no.
c                         ihelp(2) - producer no.
c                         ihelp(3) - last read record in 'innh.fort.'
c                         ihelp(4) - first record in 'innh.fort.' for
c                                    producer ihelp(2)
c                         ihelp(5) - pointer to the first 'innh.fort.'
c                                    in the first record
c                         ihelp(6) - last record in 'innh.fort.'
c
c---------------------------------------------------------------------
c  D N M I
c      written by Rebecca Rudsar    20.09.1982
c      updated by Anstein Foss      10.05.1984
c                 Anstein Foss      12.03.1986
c                 Anstein Foss      20.03.1991
c                 Anstein Foss      15.10.1992
c                 Anstein Foss      02.12.1993
c                 Anstein Foss      10.05.1994
c                 Anstein Foss      09.06.1995 ... not testing in(15)
c                 Anstein Foss      24.07.1995 ... innh(i,j)=-32767
c                 Anstein Foss      22.03.2003 ... automatic byte swap
c                 Anstein Foss      03.04.2003 ... 0-padding of records
c---------------------------------------------------------------------
c
c
c..lrec   - record length in integer*2 (16 bit) words
c..nin    - no. of innholdsfortegnelser per record
c..maxlen - maximum field length (words)
c..maxrec - maximum record no. where a field can start
c
      parameter (lrec=1024)
      parameter (nin=lrec/16)
      parameter (maxlen=327*32767)
      parameter (maxrec=100*32767)
c
c..input/output
      integer   iunit,ip,lfelt
      integer   ierr(3),ihelp(6)
      integer*2 in(16),idrec1(lrec),idrec2(lrec),ifelt(lfelt)
c
c..local
      integer   nrb(2,2)
      integer*2 innh(16,nin),ibuffr(lrec)
      integer*2 izero(lrec)
      logical   swapfile,swap
      logical   setzero
c
      data      setzero/.true./
c
      save izero
c
      swap= swapfile(iunit)
c
c..overwrite existing field: 0=no 1=yes
      ioverw=ierr(1)
c
c..if ierr(2) = 999 then test on date in 'innh'.
      ntest1=2
      if(ierr(2).eq.999) ntest1=5
      if(idrec1(1).eq.997 .or. idrec1(1).eq.998) ntest1=5
c
c..take care of 'reserved space' (after termin reset) on Cyclic Archive
      iarchive=0
      if(idrec1(1).eq.997) iarchive=1
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
      if(lfelt.lt.21) then
        ip=-5
        goto 1004
      end if
      if(lfelt.gt.maxlen) then
        ip=-6
        goto 1004
      end if
c
      if(in(1).lt.1 .or. in(1).gt.99) then
        ip=-4
        goto 1004
      end if
c
c..safety first (ihelp(1) not reset and another file with same unit)
      if(ihelp(1).eq.iunit .and. ihelp(2).eq.in(1)) then
        iprod=in(1)
        if(idrec2(2*iprod-1).lt.1 .or.
     *     idrec2(2*iprod).lt.1 .or. idrec2(2*iprod).gt.1024) then
          ip=-4
          goto 1004
        end if
        ihelp4=idrec2(2*iprod-1)
        ihelp5=idrec2(2*iprod)/16+1
        ihelp6=idrec1(9)
        if(ihelp4.ne.ihelp(4) .or. ihelp5.ne.ihelp(5)
     *                        .or. ihelp6.ne.ihelp(6)) ihelp(1)=0
      end if
c
      if(ihelp(1).ne.iunit .or. ihelp(2).ne.in(1)) then
        iprod=in(1)
        if(idrec2(2*iprod-1).lt.1 .or.
     *     idrec2(2*iprod).lt.1 .or. idrec2(2*iprod).gt.1024) then
          ip=-4
          goto 1004
        end if
        ihelp(1)=iunit
        ihelp(2)=iprod
        ihelp(3)=idrec2(2*iprod-1)
        ihelp(4)=ihelp(3)
        ihelp(5)=idrec2(2*iprod)/16+1
        ihelp(6)=idrec1(9)
      end if
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
          read(iunit,rec=ireci,iostat=ios,err=1001) innh
c
	  if(swap) call bswap2(16*nin,innh)
c
c..search for correct 'innholdsfortegnelse'
          do j=ji,nin
            k=1
            do i=1,ntest1
              if(innh(i,j).ne.in(i)) k=0
            end do
            do i=9,14
              if(innh(i,j).ne.-32767 .and. innh(i,j).ne.in(i)) k=0
            end do
            if(k.eq.1) goto 50
          end do
c
          if(innh(1,nin).ne.in(1)) goto 20
          ji=1
c
        end do
c
 20   continue
c
c..'innh' not found.
      ip=0
C
C     write innh
C
      write(*,*) '***************** Found in file:'
      do j=ji,nin
         if (innh(1,j).gt.0) then
            write(*,'(100(I10,X))') 
     &           (innh(i,j),i=1,ntest1),
     &           (innh(i,j),i=9,14)
         end if
      end do
      write(*,*) '***************** Searched for:'
      write(*,'(100(I10,X))') (in(i),i=1,ntest1),
     &     (in(i),i=9,14)
      write(*,*) '***************** NO MATCH\!'

      return
c
c..here we have found the correct 'innholdsfortegnelse'
c
   50 ihelp(3)=ireci
c
c..handle 'old' files well.........................................
      if(innh(6,j).gt.0 .and. innh(8,j).gt.0) then
        if(innh( 7,j).lt.1 .or. innh(7,j).gt.lrec) innh(7,j)=1
        if(innh(16,j).lt.0) innh(16,j)=0
      end if
      if(iarchive.eq.1) then
c..cyclic archive felt file .. pointer to reserved space
        if(innh(6,j).lt.0 .and. innh(6,j).ne.-32767 .and.
     *     innh(8,j).lt.0 .and. innh(8,j).ne.-32767) then
          if(innh( 7,j).ge.0 .or. innh(7,j).lt.-lrec) innh(7,j)=-1
          if(innh(16,j).gt.0) innh(16,j)= 0
        end if
      end if
c..................................................................
c
      iupdih=0
      iexist=1
      if(innh(6,j).lt.1 .or. innh( 7,j).lt.1 .or.
     *   innh(8,j).lt.1 .or. innh(16,j).lt.0) iexist=0
c
      if(iarchive.eq.1 .and. iexist.eq.0) then
c..cyclic archive felt file .. test if pointer to reserved space
        if(innh( 6,j).lt.0 .and.
     *     innh( 7,j).lt.0 .and. innh( 7,j).ge.-lrec .and.
     *     innh( 8,j).lt.0 .and.
     *     innh(16,j).le.0 .and. innh(16,j).ge.-32699) then
          innh( 6,j)=-innh( 6,j)
          innh( 7,j)=-innh( 7,j)
          innh( 8,j)=-innh( 8,j)
          innh(16,j)=-innh(16,j)
          iexist=1
          iupdih=1
          ioverw=1
        end if
      end if
c
c..handle 'old' files well.........................................
      if(idrec1(13).lt.1 .or. idrec1(13).gt.lrec) idrec1(13)=lrec
      if(idrec1(12).lt.0) idrec1(12)=0
c..................................................................
      nrecd1=idrec1( 8)
      nrecd2=idrec1(12)
      nrecd=nrecd1+nrecd2*32767
      lword=idrec1(13)
      nogap=idrec1(14)
c
      if(iexist.eq.0) then
        if(nogap.ne.1) then
          irecd=nrecd+1
          iword=1
        else
          irecd=nrecd
          iword=lword+1
          if(iword.gt.lrec) then
            irecd=irecd+1
            iword=1
          end if
        end if
        if(irecd.gt.maxrec) then
          ip=-7
          goto 1004
        end if
        irecd2=(irecd-1)/32767
        irecd1=irecd-irecd2*32767
        lfelt2=(lfelt-1)/32767
        lfelt1=lfelt-lfelt2*32767
        innh( 6,j)=irecd1
        innh( 7,j)=iword
        innh( 8,j)=lfelt1
        innh(16,j)=lfelt2*100+irecd2
        iupdih=1
      elseif(ioverw.eq.0) then
        ip=-2
        return
      else
        irecd1=innh( 6,j)
        iword =innh( 7,j)
        lfelt1=innh( 8,j)
        lfelt2=innh(16,j)/100
        irecd2=innh(16,j)-lfelt2*100
        if(lfelt1+lfelt2*32767.ne.lfelt) then
          ip=-3
          return
        end if
        irecd=irecd1+irecd2*32767
      end if
c
c..possibly update undefined elements in 'innholdsfortegnelse'
      do i=9,14
	if(innh(i,j).eq.-32767) then
	  innh(i,j)=in(i)
	  iupdih=1
	end if
      end do
c
c..possibly correct grid type in 'innholdsfortegnelse'
      if(innh(15,j).ne.ifelt(9)) then
	innh(15,j)=ifelt(9)
	iupdih=1
      end if
c
c..fill in the first 20 words of the data part.
      ifelt(1)=in( 1)
      ifelt(2)=in( 2)
      ifelt(3)=in( 9)
      ifelt(4)=in(10)
      ifelt(5)=in(11)
      ifelt(6)=in(12)
      ifelt(7)=in(13)
      ifelt(8)=in(14)
c.old ifelt(9)=in(15)
c
      ifelt(12)=innh(3,j)
      ifelt(13)=innh(4,j)
      ifelt(14)=innh(5,j)
c
      if(setzero) then
        do i=1,lrec
	  izero(i)=0
	end do
	setzero=.false.
      end if
c
      irecdl=irecd+(iword-1+lfelt-1)/lrec
      l2=0
c
      if(swap) call bswap2(lfelt,ifelt)
c
      do irec=irecd,irecdl
c
        l1=l2+1
        l2=min0(l2+lrec-iword+1,lfelt)
c
        lb1=0
        lb2=0
        lbz=0
        lb =0
        if(irec.eq.irecd .and. iword.gt.1) then
          lb1=iword-1
          lb =lb1
        end if
        if(irec.eq.irecdl .and. nogap.eq.1 .and. iexist.eq.1) then
          lb2=iword+l2-l1+1
	  if(lb2.le.lrec) then
            lb =lrec
	  else
	    lb2=0
	  end if
        elseif(irec.eq.irecdl) then
          lbz=iword+l2-l1+1
	  if(lbz.gt.lrec) lbz=0
        end if
c
        if(lb.eq.0) then
	  if(lbz.eq.0) then
            write(iunit,rec=irec,iostat=ios,err=1003) (ifelt(i),i=l1,l2)
	  else
            write(iunit,rec=irec,iostat=ios,err=1003)
     *			(ifelt(i),i=l1,l2),(izero(i),i=lbz,lrec)
	  end if
          lword=l2-l1+1
        else
          read(iunit,rec=irec,iostat=ios,err=1003) (ibuffr(i),i=1,lb)
          if(lb1.gt.0 .and. lb2.eq.0) then
            write(iunit,rec=irec,iostat=ios,err=1003)
     *              (ibuffr(i),i=1,lb1),(ifelt(i),i=l1,l2)
            lword=lb1+l2-l1+1
          elseif(lb1.eq.0 .and. lb2.gt.0) then
            write(iunit,rec=irec,iostat=ios,err=1003)
     *              (ifelt(i),i=l1,l2),(ibuffr(i),i=lb2,lrec)
            lword=lrec
          else
            write(iunit,rec=irec,iostat=ios,err=1003)
     *              (ibuffr(i),i=1,lb1),(ifelt(i),i=l1,l2),
     *              (ibuffr(i),i=lb2,lrec)
            lword=lrec
          end if
        end if
c
        iword=1
c
      end do
c
      if(swap) call bswap2(lfelt,ifelt)
c
      if(iupdih.eq.1) then
        if(swap) call bswap2(16*nin,innh)
        write(iunit,rec=ireci,iostat=ios,err=1002) innh
        if(swap) call bswap2(16*nin,innh)
      end if
c
      if(iexist.eq.0) then
        nrecd=irecdl
        nrecd2=(nrecd-1)/32767
        nrecd1=nrecd-nrecd2*32767
        idrec1( 8)=nrecd1
        idrec1(10)=idrec1(10)+1
        idrec1(12)=nrecd2
        if(nogap.ne.1) then
          idrec1(13)=lrec
        else
          idrec1(13)=lword
        end if
      end if
c
      do i=1,16
        in(i)=innh(i,j)
      end do
c
      ip=1
      return
c
 1001 ip=-1
      ierr(1)=ios
      ihelp(3)=ireci+1
      return
c
 1002 ip=-1
      ierr(2)=ios
      ihelp(3)=ireci+1
      if(swap) call bswap2(16*nin,innh)
      return
c
 1003 ip=-1
      ierr(3)=ios
      if(swap) call bswap2(lfelt,ifelt)
      return
c
 1004 ihelp(1)=0
c
      return
      end

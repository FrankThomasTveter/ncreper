      subroutine wcfelt(filnam,iunit,iretur)
c
c  NAME:
c     wcfelt
c
c  PURPOSE:
c     Check data pointers in a Felt file.
c     If update_flag set in the first record (word 15 = 1),
c     the 'index' will be scanned to check if any data pointer
c     (+datalength) is grater than the last_data pointer in the
c     first record (word 8,12 and 13).
c     Also checking the number of existing fields (word 10).
c     Programs writing to a felt file should set the update_flag
c     in record one, write the record, write fields and at last
c     reset and rewrite the first record.
c     (Note that this update_flag is not used as a write lock.
c      It is used to prevent wrong pointers after a program has
c      been aborted before updating the first record.)
c     The file is opened and closed.
c
c  SYNOPSIS:
c     subroutine wcfelt(filnam,iunit,iretur)
c     integer       iunit,iretur
c     character*(*) filnam
c
c  INPUT:
c     filnam  - file name (if filename(1:1)='*' the file is opened
c			    without 'file=filnam', only with unit no.)
c     iunit   - file unit no.
c
c  OUTPUT:
c     iretur  - exit status
c                 0 = update_flag was not set, file is not checked
c                 1 = file checked, no correction
c                 2 = file is corrected (record 1)
c                -1 = open/read error
c                -2 = write error (updating record 1)
c
c-----------------------------------------------------------------------
c  DNMI/FoU  28.07.1993  Anstein Foss
c  DNMI/FoU  02.12.1993  Anstein Foss
c  DNMI/FoU  22.03.2003  Anstein Foss ... automatic byte swap
c-----------------------------------------------------------------------
c
c
c..lrec   - record length in integer*2 (16 bit) words
c..ninh   - no. of innholdsfortegnelser per record
c
      parameter (lrec=1024)
      parameter (ninh=lrec/16)
c
c..input/output
      integer       iunit,iretur
      character*(*) filnam
c
c..local
      integer*2 idrec1(lrec),innh(16,ninh)
      logical   swapfile,swap
c
      call rlunit(lrunit)
c
      if(filnam(1:1).ne.'*') then
        open(iunit,file=filnam,
     *             access='direct',form='unformatted',
     *             recl=2048/lrunit,
     *             status='old',iostat=ios)
      else
        open(iunit,access='direct',form='unformatted',
     *             recl=2048/lrunit,
     *             status='old',iostat=ios)
      end if
c
      if(ios.eq.0) swap= swapfile(-iunit)
c
      if(ios.eq.0) read(iunit,rec=1,iostat=ios) idrec1
c
      if(ios.ne.0) then
        write(6,*) '**wcfelt** open/read error.  iostat= ',ios
        call pioerr(6,ios)
        write(6,*) '** ',filnam
        write(6,*) '** file unit: ',iunit
        iretur=-1
        return
      end if
c
      if(idrec1(15).eq.0) then
        close(iunit)
        iretur=0
        return
      end if
c
      if(swap) call bswap2(lrec,idrec1)
c
      write(6,*) '**wcfelt** PREVIOUS WRITE FIELD JOB HAS NOT COMPLETED'
      write(6,*) '**wcfelt** Now checking data pointers in file:'
      write(6,*) '** ',filnam
      write(6,*) '** file unit: ',iunit
c
c..Cyclic Archive Felt file ... check 'reserved field space'
      iarchive=0
      if(idrec1(1).eq.997) iarchive=1
c
c..last_record pointer
      if(idrec1(12).lt.0) idrec1(12)=0
      nrecd1=idrec1( 8)
      nrecd2=idrec1(12)
      nrecdf=nrecd1+nrecd2*32767
      nrecd=nrecdf
      lword =idrec1(13)
      if(lword.lt.1 .or. lword.gt.lrec) lword=lrec
      lwordf=lword
      nfield=idrec1(10)
      nfound=0
      nogap=idrec1(14)
c
      ireci1=3
      ireci2=idrec1(9)
c
      do ireci=ireci1,ireci2
c
        read(iunit,rec=1,iostat=ios) innh
c
        if(ios.ne.0) then
          write(6,*) '**wcfelt** read error.  iostat= ',ios
          call pioerr(6,ios)
          write(6,*) '** ',filnam
          write(6,*) '** file unit: ',iunit
          close(iunit)
          iretur=-1
          return
        end if
c
        if(swap) call bswap2(16*ninh,innh)
c
        if(iarchive.eq.1) then
c..Cyclic Archive Felt file ... check pointer to 'reserved field space'
          do n=1,ninh
            if(innh( 1,n).gt.0) then
              if(innh(6,n).lt.0 .and. innh(6,n).ne.-32767 .and.
     *           innh(8,n).lt.0 .and. innh(8,n).ne.-32767) then
                if(innh( 7,n).ge.0 .or.
     *             innh( 7,n).lt.-lrec) innh(7,n)=-1
                if(innh(16,n).gt.0) innh(16,n)= 0
              end if
              if(innh( 6,n).lt.0 .and.
     *           innh( 7,n).lt.0 .and. innh( 7,n).ge.-lrec .and.
     *           innh( 8,n).lt.0 .and.
     *           innh(16,n).le.0 .and. innh(16,n).ge.-32699) then
                innh( 6,n)=-innh( 6,n)
                innh( 7,n)=-innh( 7,n)
                innh( 8,n)=-innh( 8,n)
                innh(16,n)=-innh(16,n)
              end if
            end if
          end do
        end if
c
        do n=1,ninh
          if(innh(1,n).gt.0 .and.
     *       innh(6,n).gt.0 .and. innh(8,n).gt.0) then
            if(innh( 7,n).lt.1 .or. innh(7,n).gt.lrec) innh(7,n)=1
            if(innh(16,n).lt.0) innh(16,n)=0
            irecd1=innh( 6,n)
            iword =innh( 7,n)
            lfelt1=innh( 8,n)
            lfelt2=innh(16,n)/100
            irecd2=innh(16,n)-lfelt2*100
            irecd=irecd1+irecd2*32767
            lfelt=lfelt1+lfelt2*32767
            nrec=(iword-1+lfelt+lrec-1)/lrec
            iword=iword-1+lfelt-(nrec-1)*lrec
            if(nrecd.lt.irecd+nrec-1) then
              nrecd=irecd+nrec-1
              lword=iword
            elseif(nrecd.eq.irecd+nrec-1 .and. lword.lt.iword) then
              lword=iword
            end if
            nfound=nfound+1
          end if
        end do
c
      end do
c
c..if not file without 'data gaps'
      if(nogap.ne.1) lword=lrec
c
      iretur=0
c
      if(nrecd.eq.nrecdf .and. lword.eq.lwordf) then
        write(6,*) '**wcfelt** data pointers ok'
        iretur=1
      else
        write(6,*) '**wcfelt** last_data pointer corrected'
        write(6,*) '**         record from,to: ',nrecdf,nrecd
        write(6,*) '**         word   from,to: ',lwordf,lword
        nrecd2=(nrecd-1)/32767
        nrecd1=nrecd-nrecd2*32767
        idrec1( 8)=nrecd1
        idrec1(12)=nrecd2
        idrec1(13)=lword
        iretur=2
      end if
c
      if(nfound.eq.nfield) then
        write(6,*) '**wcfelt** no. of fields ok'
      else
        write(6,*) '**wcfelt** no. of fields corrected'
        write(6,*) '**         from,to: ',nfield,nfound
        idrec1(10)=nfound
      end if
c
c..reset
      idrec1(15)=0
c
      if(swap) call bswap2(lrec,idrec1)
c
      write(iunit,rec=1,iostat=ios) idrec1
c
      if(ios.ne.0) then
        write(6,*) '**wcfelt** write error (record 1).  iostat= ',ios
        call pioerr(6,ios)
        write(6,*) '** ',filnam
        write(6,*) '** file unit: ',iunit
        iretur=-2
      end if
c
      close(iunit)
c
      return
      end

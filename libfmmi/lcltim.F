      subroutine lcltim(itime,iunit,ierror)
c
c        'itime' is updated to give 'valid' local (national) date,time
c                                           (with prog.time = 0)
c        Input 'itime' is UTC time.
c        information about local (national) time is read from the first
c        readable file in this order:
c          1) ./local.time
c          2) /usr/local/default/local.time
c
c        Input:  itime(5) - itime(1): year
c                           itime(2): month (1-12)
c                           itime(3): day (1-28/29/30/31)
c                           itime(4): time in hours (00-23)
c                           itime(5): time in hours of prognosis
c                                     (negative, zero or positive)
c                iunit    - unit no. used for 'local.time'
c                           (file is opened and closed)
c        Output: itime(5) - as above, with itime(5)=0
c                ierror   - exit status:
c                             0 = o.k.
c                             1 = not o.k. input date/time
c                                 ('itime' not changed)
c                             2 = 'local.time' file not found
c                                 ('itime' returned as valid UTC)
c
c       NOTE. The 'local.time' file is only read the first time this
c             routine is called. If the file contains information
c             about more than 32 (summer) 'corrections' only the last
c             32 will be used.
c
c-----------------------------------------------------------------------
c The local.time file:
c --------------------
c ********* local.time ********** Norsk normal- og sommer-tid **********
c +1                                <- standard forskjell fra UTC, timer
c +2, 1993,03,28,01, 1993,09,26,02  <- UTC, sommer-tid 1993
c +2, 1994,03,28,01, 1994,09,26,02  <- UTC, sommer-tid 1994 ????????
c +2, 1995,03,28,01, 1995,09,26,02  <- UTC, sommer-tid 1995 ????????
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c  DNMI/FoU  01.06.1993  Anstein Foss
c  DNMI/FoU  28.08.2003  Anstein Foss ... added /metno/local/default/...
c-----------------------------------------------------------------------
c
c..input/output
      integer itime(5),iunit,ierror
c
c..local
      parameter (nfiles=3)
      parameter (maxtim=50)
c
      integer numtim,itimex(4,2,maxtim),ihdiff(maxtim),itimey(4,2)
c
      character*40 files(nfiles)
c
      data files/'./local.time',
     *           '/metno/local/default/local.time',
     *           '/usr/local/default/local.time'/
c
      save itimex
      save ihdiff
c
      data numtim/-1/
      data nhdiff/0/
c
      call vtime(itime,ierror)
      if(ierror.ne.0) then
        ierror=1
        return
      end if
c
      if(numtim.eq.-1) then
c..read the 'local.time' file
c
        ios=-1
        n=0
        do while (ios.ne.0 .and. n.lt.nfiles)
          n=n+1
          open(iunit,file=files(n),
     *               access='sequential',form='formatted',
     *               status='old',iostat=ios)
          if(ios.eq.0) read(iunit,fmt='(a)',iostat=ios) string
          if(ios.eq.0) read(iunit,*,iostat=ios) nhdiff
          if(ios.ne.0) close(iunit)
        end do
c
        if(ios.ne.0) then
          write(6,*) '**LCLTIM** No o.k.  local.time  file was found.'
          do n=1,nfiles
            write(6,*) '**LCLTIM** ',files(n)
          end do
          ierror=2
          return
        end if
c
        numtim=0
        do while (ios.eq.0)
          read(iunit,*,iostat=ios) ihd,itimey
          if(ios.eq.0) then
            if(numtim.eq.maxtim) then
              numtim=numtim-1
              do n=1,numtim
                ihdiff(n)=ihdiff(n+1)
                do j=1,2
                  do i=1,4
                    itimex(i,j,n)=itimex(i,j,n+1)
                  end do
                end do
              end do
            end if
            numtim=numtim+1
            ihdiff(numtim)=ihd
            do j=1,2
              do i=1,4
                itimex(i,j,numtim)=itimey(i,j)
              end do
            end do
          end if
        end do
c
        close(iunit)
c
      end if
c
      ntim=0
      n=0
      do while (ntim.eq.0 .and. n.lt.numtim)
        n=n+1
        igt1=0
        igt2=0
        do i=1,4
          if(igt1.eq.0 .and. itime(i).lt.itimex(i,1,n)) igt1=-1
          if(igt1.eq.0 .and. itime(i).gt.itimex(i,1,n)) igt1=+1
          if(igt2.eq.0 .and. itime(i).lt.itimex(i,2,n)) igt2=-1
          if(igt2.eq.0 .and. itime(i).gt.itimex(i,2,n)) igt2=+1
        end do
        if(igt1.ge.0 .and. igt2.le.0) ntim=n
      end do
c
      ierror=0
c
      nhours=nhdiff
      if(ntim.gt.0) nhours=ihdiff(ntim)
      if(itime(4)+nhours.ge.00 .and. itime(4)+nhours.le.23) then
        itime(4)=itime(4)+nhours
      else
        itime(5)=nhours
        call vtime(itime,ierror)
      end if
c
      return
      end

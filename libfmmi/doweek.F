      subroutine doweek(itime,numday,ierror)
c
c        Day Of WEEK.
c        Returns day no. 1-7, 1 is Monday and 7 is Sunday.
c
c        Input:  itime(5) - itime(1): year
c                           itime(2): month (1-12)
c                           itime(3): day (1-28/29/30/31)
c                           itime(4): time in hours (00-23)
c                           itime(5): time in hours of prognosis
c                                     (negative, zero or positive)
c        Output: numday   - day no. 1-7 meaning Monday-Sunday
c                           (0 if ierror =/= 0)
c                ierror   - exit status:
c                             0 = o.k.
c                             1 = not o.k. input date/time
c
c-----------------------------------------------------------------------
c  DNMI/FoU  01.06.1993  Anstein Foss
c-----------------------------------------------------------------------
c
c..input/output
      integer itime(5),numday,ierror
c
c..local
      integer idayr,itimer(5)
c
c..reference day: Friday 01.01.1993
      data idayr/5/
      data itimer/1993,01,01,00,+0/
c
      call hrdiff(0,0,itimer,itime,nhours,ierr1,ierr2)
      if(ierr1.ne.0 .or. ierr2.ne.0) then
        numday=0
        ierror=1
        return
      end if
c
      if(nhours.ge.0) then
        nday=idayr+nhours/24
      else
        nday=idayr+(nhours-23)/24
      end if
c
      if(nday.gt.0) then
        numday=nday-((nday-1)/7)*7
      else
        numday=nday-(nday/7)*7+7
      end if
c
      ierror=0
c
      return
      end

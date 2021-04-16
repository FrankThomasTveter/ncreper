      subroutine daytim(year,month,day,hour,minute,second)
c
c  PURPOSE: Return machine date and time.
c
c  OUTPUT:  year   - integer -  year  (4 digits)
c           month  - integer -  month
c           day    - integer -  day
c           hour   - integer -  hour
c           minute - integer -  minute
c           second - integer -  second
c
c  WARNING: System routine 'IDATE' only returns year with 2 digits.
c           Fixup to 4 digits here.
c
c-----------------------------------------------------------------------
c  DNMI/FoU  24.10.1992  Anstein Foss
c-----------------------------------------------------------------------
c
c..output
      integer year,month,day,hour,minute,second
c..local
      integer time(3)
c
      integer date(3)
      call idate(date)
      day=date(1)
      month=date(2)
      year=date(3)
c
      if (year.lt.100) then
         if(year.ge.70) then
            year=year+1900
         else
            year=year+2000
         end if
      end if
c
      call itime(time)
c
      hour  =time(1)
      minute=time(2)
      second=time(3)
c
      return
      end

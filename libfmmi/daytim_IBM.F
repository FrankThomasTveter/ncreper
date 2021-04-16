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
c  DNMI/FoU  30.11.1995  Anstein Foss ... IBM/RS 6000 xl fortran
c-----------------------------------------------------------------------
c
c..output
      integer year,month,day,hour,minute,second
c
c..standard unix version.............................................
c..local
cc    integer time(3)
c
cc    call idate(month,day,year)
c
cc    if(year.ge.70) then
cc      year=year+1900
cc    else
cc      year=year+2000
cc    end if
c
cc    call itime(time)
c
cc    hour  =time(1)
cc    minute=time(2)
cc    second=time(3)
c..standard unix version.............................................
c
c..IBM RS/6000 xl fortran version....................................
c..local
      type iardate
	sequence
	integer*4 iday
	integer*4 imonth
	integer*4 iyear
      end type
      type iartime
	sequence
	integer*4 ihour
	integer*4 imin
	integer*4 isec
      end type
      type (iardate) idate_struct
      type (iartime) itime_struct
c
      call idate_(idate_struct)
      call itime_(itime_struct)
c
      year   = iyear
      month  = imonth
      day    = iday
      hour   = ihour
      minute = imin
      second = isec
c
      if(year.lt.100) year=year+1900
c..IBM RS/6000 xl fortran version....................................
c
      return
      end

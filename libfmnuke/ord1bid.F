       integer function ord1bid(saf250,l1bscnum,jdate)
C+
CNOM/NAME  ord1bid.F version 1.1 on 2/12/98
C
CBUT/PURPOSE
C       delivers the order satellite number for a given
C       level 1b id number (PACS number) and a julian date.
C
CDESCRIPTION
C     FORTRAN 77 witth UKMO/METEO-FRANCE accepted extensions
C       The function searches in the array avh_sc_num the element which
C       has the same value as the input argument. A date is given to
C       avoid the conflicts between satellites
C
C       the function returns 0 if the satellite is not found or
C       if an error occur
C
C
CPROPRIETAIRE/OWNER
C       M. Derrien
C
CHISTORIQUE/HISTORY (MODIFICATIONS/CHANGES)
C(version   date      comment           <author>)
C 1.0     13/02/95   initial            M. Derrien
C
C
CFUNCTIONS AND SUBROUTINES CALLS  (RETURN CODES FOR FUNCTIONS)
C
C
CCOPYRIGHT (c) Meteo-France (EVOLUNIX - ATOVS - EUMETSAT project <1994>)
C
      implicit none
      save
C
CARGUMENTS
       integer l1bscnum !  satellite id read from a 1B file
       integer jdate     ! julian date of the data set
C
CDECLARATIONS
C
C identification variable for "what"
      character*62 sccsid
C
      data sccsid
     &/'@(#) ord1bid.F version 1.1 on 2/12/98'/
C
C Common declarations:
       integer maxsatid
       parameter (maxsatid=15)
       integer noa_sc_num(maxsatid)       ! noaa number
       integer l1b_sc_num(maxsatid)       ! level 1b number (PACS)
       integer avh_sc_num(maxsatid)       ! avhrr telemetry number
       integer tip_sc_num(2,maxsatid)     ! tip telemetry number
       integer sc_date_range(2,maxsatid)  ! operational date range
       character*8 launch_name(maxsatid)  ! pre-launch name
       character*9 sc_int_des(maxsatid)   ! s/c international designation
       character*8 noaa_name(maxsatid)    ! noaa name (noaaxx)
       character*8 alias1(maxsatid)       ! first alias of noaa name (noaax)
       character*8 alias2(maxsatid)       ! second alias of noaa name (noaa-xx)
       character*8 alias3(maxsatid)       ! third alias of noaa name (noaa-x)
       character*250 saf250
C
       common/c_satid/noa_sc_num, l1b_sc_num, avh_sc_num,
     &                 tip_sc_num, sc_date_range, launch_name,
     &                 sc_int_des, noaa_name, alias1, alias2, alias3
C
C Local variables:
       integer ierr         ! rdnoaaid return code =0 OK  else I/O error
       integer isat         ! loop satellite order number
       logical first
       data first /.true./
C
C-
C
C       checks if it is the first calling sequence and then
C       call the file reading routine
       ierr = 0
       if(first) then
         call rdnoaaid(saf250,ierr)
         first = .false.
       endif
C
C       checks the return code
       if(ierr.eq.0) then
C
C         search the argument value in the common array for the
C         right avhrr telem number in the good range of dates.
         isat=1
          do while ( (l1b_sc_num(isat).ne.l1bscnum  .or.
     &                jdate.lt.sc_date_range(1,isat) .or.
     &                jdate.gt.sc_date_range(2,isat) ).and.
     &                           isat.le.maxsatid )
           isat=isat+1
         end do
C
C         the satellite id is found
         if(isat.le.maxsatid) then
           ord1bid = isat
         else
C           the satellite id is not present in the array
           ord1bid = 0
         endif
C
       else
C         an error occured during the file reading
         ord1bid = 0
       endif
C
       return
       end

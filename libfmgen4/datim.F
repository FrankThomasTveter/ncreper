      SUBROUTINE DATIM(Systim)

      IMPLICIT NONE

! ----------------------------------------------------------------------
! --- Program      : DATIM subroutine.
! ---
! --- Description  : This subroutine gets the ststem date and time and
! ---              : passes it back to the calling program.
! ---
! --- Arguments    : Systim (o/p)   : INTEGER*4 Systim(8).
! ---              :
! ---              :                : Systim(1)  - not set.
! ---              :                : Systim(2)  - not set.
! ---              :                : Systim(3)  - second  (0-59) 
! ---              :                : Systim(4)  - minute  (0-59)
! ---              :                : Systim(5)  - hour    (0-23)
! ---              :                : Systim(6)  - day     (0-31)
! ---              :                : Systim(7)  - month   (1-12)
! ---              :                : Systim(8)  - year    (yyyy)
! ---
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! --- Declare variables.
! ----------------------------------------------------------------------

       INTEGER           Systim(8)          !- system date time returned
      integer           values(8)

      call date_and_time(VALUES=values)

! ----------------------------------------------------------------------
! --- Internal read to put time into Systim array.
! ----------------------------------------------------------------------
      
      Systim(3)=values(7)       !- second. 0-59
      Systim(4)=values(6)       !- minute. 0-59
      Systim(5)=values(5)       !- hour    (0-23)
      Systim(6)=values(3)       !- day     (0-31)
      Systim(7)=values(2)       !- month   (1-12)
      Systim(8)=values(1)       !- year    (yyyy)
 !
      RETURN
      END

SUBROUTINE MNCREPER(UNITI,IRC)
  ! 
  ! ***********************************************************************
  ! +                                                                     *
  ! +  UNITI = UNIT NUMBER FOR INPUT FILE                                 *
  ! +  IRC=ERROR RETURN CODE (0=OK)                                       *
  ! +                                                                     *
  ! +                                                                     *
  ! VERSION                      : 01/01/00                               *
  ! +                                                                     *
  ! WRITTEN/MODIFIED BY:                                                  *
  ! --------------------------------------------------------------------- *
  ! |    NAME      |   DATE   |                 REASON                  | *
  ! --------------------------------------------------------------------- *
  ! | F. TVETER    | 01/08/14 | NEW                                     | *
  ! |              |          |                                         | *
  ! --------------------------------------------------------------------- *
  ! ***********************************************************************
  ! 
  use ncf
  use sort
  IMPLICIT NONE
  SAVE
  ! 
  ! INTERFACE VARIABLES
  ! 
  INTEGER  UNITI
  INTEGER IRC
  ! 
  CHARACTER*14 MYNAME
  DATA MYNAME /'MNCREPER'/
  logical :: bdeb=.false.
  ! 
  type varptr
     type(variable), pointer :: ptr=>null()
  end type varptr
  !
  integer, external :: length
  !
  INTEGER, PARAMETER :: NRHDR=100
  INTEGER, PARAMETER :: MAXDIFF=100
  INTEGER, PARAMETER :: MAXFILTER=100
  INTEGER, PARAMETER :: MAXCOR=10
  INTEGER, PARAMETER :: MAXIMP=10
  INTEGER, PARAMETER :: MAXSPD=10, MAXRP=10
  INTEGER, PARAMETER :: MAXPST=10
  INTEGER, PARAMETER :: MAXGRP = 100
  INTEGER, PARAMETER :: MAXVAR = 100
  INTEGER, PARAMETER :: MAXRET = 100
  INTEGER, PARAMETER :: MAXAUX = 5
  INTEGER, PARAMETER :: MAXAUD = 5
  REAL, PARAMETER :: eps=1.0D-10
  LOGICAL, PARAMETER :: uniq=.false.
  !
  type init
     LOGICAL :: LFLDAT(NRHDR)
     character*250 :: inp250, out250, aux250
     character*250 :: dim250, gus250, guv250, gud250
     character*250 :: rpd250
     integer :: leni,lena
     !
     integer :: nraux = 0
     character*250 :: avar250(maxaux)
     integer       :: lenavar(maxaux)
     integer       :: nraud(maxaux)
     character*250 :: adim250(maxaux,maxaud)
     !
     integer :: nrdiff
     character*250 :: niff250(maxdiff),diff250(maxdiff)
     real :: hdiff(maxdiff), vmin(maxdiff),vmax(maxdiff)
     logical :: vset(maxdiff)
     integer :: lenniff(maxdiff),lendiff(maxdiff)
     !
     integer :: nriflt
     character*250 :: iflt250(maxfilter)
     real :: miniflt(maxfilter), maxiflt(maxfilter), minf,maxf
     integer :: leniflt(maxfilter),iiflt
     logical :: biflt
     integer ::nrcor
     character*250 :: cor250(maxcor)
     integer :: nrimp
     character*250 :: imp250(maximp)
     !
     integer :: nrspd, nrrp, islice
     character*250 :: spd250(maxspd),wx250(maxspd),wy250(maxspd),rp250(maxrp)
     integer :: gul, gudEntry
     !
     integer :: nrpst
     real :: trg, zero, dimpst(maxpst), outpst(maxpst), vvmax(maxspd)
     integer :: ii, jj, ss, tt, uu, ik, jk
     logical :: lavg=.false.
     logical :: lpst=.false.
     !
     integer :: nrgrp = 0
     character*250 :: grp250(maxgrp),ogr250(maxgrp)
     integer :: lengrp(maxgrp),lenogr(maxgrp)
     !
     integer :: nrvar = 0
     character*250 :: var250(maxvar),new250(maxvar)
     integer :: lenv(maxvar),lenn(maxvar)
     !
     integer :: nrret = 0
     character*250 :: net250(maxret),ret250(maxret),per250(maxret)
     integer :: lennet(maxret),lenret(maxret),lenper(maxret)
     !
  end type init
  !
  type run
     logical :: initWrite
     logical brp,bok,bop,aop
     !
     type(dimension), pointer :: dslice=>null() ! slice-dimension
     type(dimension), pointer :: aslice=>null() ! auxiliary slice-dimension
     integer :: cll,nll,nslice,dnslice
     !
     integer :: itime
     type(variable), pointer  :: iFrtVariable=>null()
     type(variable), pointer  :: iTimeVariable=>null()
     type(dimension), pointer  :: iTimeDim=>null()
     type(dimensionOrder), pointer  ::  ixydo=>null()
     type(dimensionOrder), pointer  ::  axydo=>null()
     !
     integer :: tnrrp
     real :: t2000
     !
     type(inventory), pointer  :: ret => null()   ! return period...
     type(dimensionOrder), pointer  ::  odo=>null(),ido=>null(),ado=>null()
     type(dimension), pointer :: rz=>null() ! return level dimension
     type(dimension), pointer :: ro=>null() ! return level dimension in output file
     !
  end type run
  !
  type(inventory), pointer  :: inp => null()   ! input...
  type(inventory), pointer  :: aux => null()   ! auxiliary...
  type(inventory), pointer  :: out => null()   ! output...
  !
  type(init) :: i
  type(run)  :: r
  integer :: ll
  ! 
  ! debug pointers
  !type(variable), pointer :: debugvar=>null()
  !
  !#     include "netcdf.inc"
  !
  IRC=0
  ! 
  if(bdeb)write(*,*) myname,' Debug: Routine starts.',irc
  !
  call readInit(uniti,i,irc)
  if(irc.ne.0) then
     write(*,*)myname,' Unable to allocate readInit.',irc
     return
  end if
  !
  ! allocate file inventories
  !
  allocate(inp,stat=irc)
  if(irc.ne.0) then
     write(*,*)myname,' Unable to allocate inventory.',irc
     return
  end if
  call ncf_setLabel(inp,"input")
  !
  ! open input file and read inventory
  !
  r%bop=.true.
  call readInput(inp,i,r,irc)
  if(irc.ne.0) then
     write(*,*)myname,' Error return from readInput.',irc
     return
  end if
  !
  ! process inventory
  !
  call processInput(inp,i,r,irc)
  if(irc.ne.0) then
     write(*,*)myname,' Error return from processInput.',irc
     return
  end if
  !
  ! make output file inventory
  !
  out => ncf_copyFile(inp,irc)
  if (irc.ne.0) then
     write(*,*)myname,"Error return from ncf_copyFile.",irc
     return
  end if
  call ncf_setLabel(out,"output")
  !
  ! prepare return periods...
  !
  call prepareRP(inp,out,i,r,irc)
  if (irc.ne.0) then
     write(*,*)myname,"Error return from prepareRP.",irc
     return
  end if
  !
  ! allocate auxiliary file
  !
  if (i%lfldat(15)) then
     allocate(aux,stat=irc)
     if(irc.ne.0) then
        write(*,*)myname,' Unable to allocate inventory.',irc
        return
     end if
     call ncf_setLabel(aux,"aux")
     r%aop=.true.
     call readAuxiliary(aux,i,r,irc)
     if(irc.ne.0) then
        write(*,*)myname,' Error return from readAuxiliary.',irc
        return
     end if
     call makeAuxiliary(aux,out,i,r,irc)
     if(irc.ne.0) then
        write(*,*)myname,' Error return from readAuxiliary.',irc
        return
     end if
  end if
  !
  ! prepare slicing
  !
  call prepareSlice(inp,out,aux,i,r,irc)
  if (irc.ne.0) then
     write(*,*)myname,"Error return from prepareSlice.",irc
     return
  end if
  !
  SLICE: do ll=1,r%nslice,r%dnslice ! r%nslice
     !
     call initSlice(inp,out,aux,i,r,ll,irc)
     if (irc.ne.0) then
        write(*,*)myname,"Error return from initSlice.",irc
        return
     end if
     !
     ! call ncf_checkinventory(out,irc)
     ! if (irc.ne.0) then
     !   write(*,*)myname,' Error return from ncf_checkinventory.',irc
     !   return
     ! end if
     !
     ! process variables
     !
     if (r%bop) then
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! process wind x,y -> spd, if requested
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        if (i%lfldat(40)) then
           call addWindSpd(inp,out,i,r,irc)
           if (irc.ne.0) then
              write(*,*)myname,"Error return from processWindSpd.",irc
              return
           end if
        end if
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! process wind gust, if requested
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        if (i%lfldat(45)) then
           call addWindGust(inp,out,i,r,irc)
           if (irc.ne.0) then
              write(*,*)myname,"Error return from processWindGust.",irc
              return
           end if
        end if
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! process differentiated variables
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        if (i%lfldat(60))then
           call addDifferentiated(inp,out,i,r,irc)
           if (irc.ne.0) then
              write(*,*)myname,"Error return from addDifferentiated.",irc
              return
           end if
        end if
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Add return period
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        if (i%lfldat(70).and.i%lfldat(30)) then
           call addReturnPeriod(inp,out,i,r,irc)
           if (irc.ne.0) then
              write(*,*)myname,' Error return from addReturnPeriod.',irc
              return
           end if
        end if
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ADD AVG/PST VARIABLE...
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        if (i%lfldat(88).and.(i%lavg.or.i%lpst)) then
           call addExtra(inp,out,i,r,irc)
           if (irc.ne.0) then
              write(*,*)myname,' Error return from addExtra.',irc
              return
           end if
        end if
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ADD AUXILIARY OUTPUT
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        if (i%lfldat(15).and.r%aop) then
           call addAuxiliary(aux,out,i,r,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from addAuxiliary.',irc
              return
           end if
        end if
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! IGNORE UNWANTED OR LOAD WANTED VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        if (i%lfldat(90)) then
           call readyVariables(inp,out,aux,i,irc)
           if (irc.ne.0) then
              write(*,*)myname,' Error return from readyVariables.',irc
              return
           end if
        end if
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! filter input variable
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !
        if (i%lfldat(50)) then
           !call filterInput(v,irc)
           !if (irc.ne.0) then
           !   write(*,*)myname,'Error return from filter.',irc
           !   return
           !end if
        end if
        !
        !
        !debugvar => ncf_getVariable(out,"x",bok,irc)
        !call ncf_printVariable(debugvar)
        !
     end if
     if(bdeb)write(*,*)myname,'Wrapping up.'
     !call ncf_checkInventory(inp,irc)
     !
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! DUMP SLICE TO OUTPUT FILE
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     if (r%bop) then
        call writeSlice(inp,out,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from writeSlice.',irc
           return
        end if
     else
        write(*,*)myname,' No data available for output file.'
     end if
     !
  end do SLICE
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WRAP UP AND WRITE OUTPUT FILE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  call wrapUp(inp,out,aux,i,r,irc)
  if (irc.ne.0) then
     write(*,*)myname,'Error return from wrapUp.',irc
     return
  end if
  if (associated(out)) deallocate(out)
  if (associated(inp)) deallocate(inp)
  if (associated(aux)) deallocate(aux)
  !
contains
  !
  subroutine readInit(uniti,i,irc)
    implicit none
    integer :: uniti
    type(init) :: i
    integer :: irc
    INTEGER  KODE, LINE
    INTEGER  II
    CHARACTER*250 HDR250(NRHDR),DAT250, buff250
    CHARACTER*250 NUKEHEAD
    EXTERNAL NUKEHEAD
    INTEGER  INTOUT
    LOGICAL  ENDOFF
    CHARACTER*14 MYNAME
    DATA MYNAME /'readInit'/
    integer :: lenb,lend, lenh
    ! 
    !
    i%nraux=0
    i%islice=1
    i%nrdiff=0
    i%nriflt=0
    i%nrcor=0
    i%nrspd=0
    i%nrpst=0
    i%nrrp=0
    i%nrimp=0
    i%nrvar=0
    i%zero=0.0D0
    !
    DO II=1,NRHDR
       HDR250(II) = ''
       I%LFLDAT(II)=.FALSE.
    ENDDO
    !
    HDR250(1)   = 'NCREPER V1.0 [0]VFLR'
    HDR250(10)  = 'INPUT FILE (NETCDF) [1]VFLR &'
    HDR250(15)  = 'AUXILIARY FILE (NETCDF) [1]VFLR %3&'
    HDR250(16)  = 'USE AUXILIARY VAR, REMOVE DIMENSIONS [*]VFLR %3%'
    HDR250(20)  = 'OUTPUT FILE (NETCDF) [1]VFLR &'
    HDR250(25)  = 'NUMBER OF DATA SLICES [1]VFLR'
    HDR250(30)  = 'RETPER FILE (NETCDF) [*]VFLR%1&'
    HDR250(31)  = 'ZERO RETPER LEVEL [1]VFLR%1&$'
    HDR250(35)  = 'PUT RETPER VARS IN OUTPUT FILE [0]VFLR%1&$'
    HDR250(39)  = 'IMPORT VARIABLES FROM RETPER FILE [*]VFLR%1&$'
    HDR250(40)  = 'SPEED FROM VARX, VARY, MAX [*]VFLR'
    HDR250(45)  = 'WIND GUST FROM VAR, DIM, DIMVALUE [1]VFLR'
    HDR250(50)  = 'FILTER INPUT VAR, MIN, MAX [*]VFLRM'
    HDR250(55)  = 'OUTER COORDINATE DIMENSIONS [*]VFLR'
    HDR250(60)  = 'NEW VAR, OLD VAR, DIFFERENTIATE OVER HOURS, MIN, MAX [*]VFLRM'
    HDR250(70)  = 'NEW VAR, OLD VAR, RETPER VAR [*]VFLR%1&'
    HDR250(80)  = 'RETPER DIM[1]VFLR%1&'
    HDR250(85)  = 'ENSEMBLE DIMENSION [1]VFLR%2&1'
    HDR250(86)  = 'MAKE ENSEMBLE AVERAGE [0]VFLR%2%1'
    HDR250(87)  = 'ENSEMBLE PERCENTILE, ID [*]VFLR%2%1'
    HDR250(88)  = 'NEW VAR, OLD ENSEMBLE VAR [*]VFLR%2&'
    HDR250(90)  = 'KEEP VARS [*]VFLR'
    ! 
    WRITE(*,*) MYNAME,' ----------------------------------------'
    ! 
    ! READ DATA FROM INPUT FILE..............................
    !
    KODE=-1
    CALL NUKEM(KODE,UNITI,HDR250,INTOUT,DAT250,ENDOFF,IRC)
    IF (IRC.NE.0) THEN
       WRITE(*,*) MYNAME,' -1 Error return from NUKEM.',IRC
       RETURN
    END IF
    LINE=INTOUT
    ! 
    KODE=0
    lend=1
    DO WHILE (.NOT.ENDOFF)
       ! 
       ! READ NEXT DATA LINE
       ! 
       CALL NUKEM(KODE,UNITI,HDR250,INTOUT,DAT250,ENDOFF,IRC)
       IF (IRC.NE.0) THEN
          WRITE(*,*) MYNAME,' 0 Error return from NUKEM.',IRC
          RETURN
       END IF
       LINE=INTOUT
       lend=length(dat250,250,lend)
       !     WRITE(*,*) MYNAME,' line=',line
       ! 
       IF (BDEB) WRITE(*,*) MYNAME,' Debug: Read header:',LINE
       ! 
       ! CHECK WHAT LINE WE JUST READ
       ! 
       IF (LINE.EQ.1) THEN
       ELSEIF (LINE.EQ.10) THEN ! input file
          i%inp250=dat250
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.15) THEN ! auxiliary file
          i%aux250=dat250
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.16) THEN ! aux variable and dimensions...
          i%nraux=min(maxaux,i%nraux+1)
          i%avar250(i%nraux)=nukehead(dat250,250)
          i%lenavar(i%nraux)=length(i%avar250(i%nraux),250,10)
          i%nraud(i%nraux)=0
          lend=length(dat250,250,10)
          do while (lend.ne.0)
             i%nraud(i%nraux)=min(maxaud,i%nraud(i%nraux)+1)
             i%adim250(i%nraux,i%nraud(i%nraux))=nukehead(dat250,250)
             lend=length(dat250,250,10)
          end do
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.20) THEN ! output file
          i%out250=dat250
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.25) THEN ! number of data slices
          lend=length(dat250,250,10)
          read(dat250(1:lend),*,err=201) i%islice
          i%islice=max(1,min(100,i%islice))
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.30) THEN ! retper file
          i%nrrp=min(maxrp,i%nrrp+1)
          i%rp250(i%nrrp)=dat250
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.31) THEN ! reper zero level
          read(dat250(1:lend),*,err=201) i%zero
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.35) THEN ! put retper in output file
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.39) THEN ! import form retper
          i%nrimp=min(maximp,i%nrimp+1)
          i%imp250(i%nrimp)=dat250
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.40) THEN ! wind speed from x, y
          i%nrspd=min(maxspd,i%nrspd+1)
          i%spd250(i%nrspd)=nukehead(dat250,250)
          i%wx250(i%nrspd)=nukehead(dat250,250)
          i%wy250(i%nrspd)=nukehead(dat250,250)
          lend=length(dat250,250,10)
          read(dat250(1:lend),*,err=201) i%vvmax(i%nrspd)
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.45) THEN ! wind gust
          i%gus250=nukehead(dat250,250) ! new name
          i%guv250=nukehead(dat250,250) ! variable
          i%gud250=nukehead(dat250,250) ! dimension
          lend=length(dat250,250,10)
          read(dat250(1:lend),*,err=201) i%gul ! dimension level
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.50) THEN ! time filter var, hours
          i%nriflt=min(i%nriflt+1,maxfilter)
          i%iflt250(i%nriflt)=nukehead(dat250,250)
          i%leniflt(i%nriflt)=length(i%iflt250(i%nriflt),250,10)
          lend=length(dat250,250,10)
          read(dat250(1:lend),*,err=201) i%miniflt(i%nriflt),&
               & i%maxiflt(i%nriflt)
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.55) THEN ! ignore essential errors
          i%nrcor=min(i%nrcor+1,maxcor)
          i%cor250(i%nrcor)=dat250
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.60) THEN ! time diff var, hours
          i%nrdiff=min(i%nrdiff+1,maxdiff)
          i%niff250(i%nrdiff)=nukehead(dat250,250) ! new variable
          i%diff250(i%nrdiff)=nukehead(dat250,250) ! old variable
          i%lenniff(i%nrdiff)=length(i%niff250(i%nrdiff),250,10)
          i%lendiff(i%nrdiff)=length(i%diff250(i%nrdiff),250,10)
          lend=length(dat250,250,10)
          buff250=nukehead(dat250,250) ! old variable
          lenb=length(buff250,250,10)
          read(buff250(1:lenb),*,err=201) i%hdiff(i%nrdiff)
          lend=length(dat250,250,10)
          i%vset(i%nrdiff)=(lend.ne.0)
          if (i%vset(i%nrdiff)) then
             read(dat250(1:lend),*,err=201) i%vmin(i%nrdiff),i%vmax(i%nrdiff)
          end if
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.70) THEN ! new var, old var, retper var
          i%nrret=min(i%nrret+1,maxret)
          i%net250(i%nrret)=nukehead(dat250,250) ! new variable
          i%ret250(i%nrret)=nukehead(dat250,250) ! old variable
          i%per250(i%nrret)=nukehead(dat250,250) ! retper variable
          i%lennet(i%nrret)=length(i%net250(i%nrret),250,10)
          i%lenret(i%nrret)=length(i%ret250(i%nrret),250,10)
          i%lenper(i%nrret)=length(i%per250(i%nrret),250,10)
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.80) THEN ! return period dim
          i%rpd250=dat250
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.85) THEN ! ensemble dimension
          i%dim250=dat250
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.86) THEN ! make average 
          i%lavg=.true.
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.87) THEN ! make percentile
          i%lpst=.true.
          i%nrpst=min(i%nrpst+1,maxpst)
          buff250=nukehead(dat250,250) ! new variable
          lenb=length(buff250,250,10)
          read(buff250(1:lenb),*,err=201) i%dimpst(i%nrpst)
          i%dimpst(i%nrpst)=max(0.0D0,min(100.0D0,i%dimpst(i%nrpst)))
          lend=length(dat250,250,10)
          if (lend.eq.0) then
             i%outpst(i%nrpst)=i%dimpst(i%nrpst)
          else
             read(dat250(1:lend),*,err=201) i%outpst(i%nrpst)
          end if
          i%outpst(i%nrpst)=max(0.0D0,min(100.0D0,i%outpst(i%nrpst)))
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.88) THEN ! new group var, old var
          i%nrgrp=min(i%nrgrp+1,maxgrp)
          i%grp250(i%nrgrp)=nukehead(dat250,250) ! new variable
          i%ogr250(i%nrgrp)=nukehead(dat250,250) ! old variable
          i%lengrp(i%nrgrp)=length(i%grp250(i%nrgrp),250,10)
          i%lenogr(i%nrgrp)=length(i%ogr250(i%nrgrp),250,10)
          I%LFLDAT(LINE) = .TRUE.
       ELSEIF (LINE.EQ.90) THEN ! keep variables
          i%nrvar=min(i%nrvar+1,maxvar)
          i%var250(i%nrvar)=nukehead(dat250,250) ! new variable
          i%new250(i%nrvar)=dat250
          i%lenv(i%nrvar)=length(i%var250(i%nrvar),250,10)
          i%lenn(i%nrvar)=length(i%new250(i%nrvar),250,10)
          I%LFLDAT(LINE) = .TRUE.
       ELSE IF (LINE.NE.0) THEN ! LINE.EQ.0 IMPLIES SOMETIMES EOF
          WRITE(*,*) MYNAME,' System error, line not implemented:',LINE
          IRC=999
          RETURN
       ENDIF
    ENDDO
    ! 
    KODE=1
    CALL NUKEM(KODE,UNITI,HDR250,INTOUT,DAT250,ENDOFF,IRC)
    IF (IRC.NE.0) THEN
       WRITE(*,*) MYNAME,' Error return from NUKEM.',IRC
       RETURN
    END IF
    ! 
    GOTO 202
201 CONTINUE
    LEND=LENGTH(DAT250,250,1)
    CALL CHOP0(HDR250(LINE),250)
    LENH=LENGTH(HDR250(LINE),250,1)
    WRITE(*,*) MYNAME,' unable to read: ',DAT250(1:LEND),' (HDR='//HDR250(LINE)(1:LENH)//')'
    IRC=522
    RETURN
202 CONTINUE
  ! 
    WRITE(*,*) MYNAME,' ----------------------------------------'
    return
  end subroutine readInit
  !
  ! read input file
  !
  subroutine readInput(inp,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'readInput'/
    !
    type(dimension), pointer  :: xyDim=>null()
    integer :: ii,jj, lenc
    character*250 :: dum250
    ! 
    i%leni=length(i%inp250,250,10)
    write(*,*)myname,' Scanning: ',i%inp250(1:i%leni)
    call ncf_openFile(i%inp250,inp,r%bop,irc)
    if (irc.ne.0) then
       write(*,*)myname,' Error return from ncf_openFile.',irc
       return
    end if
    if (r%bop) then
       call ncf_readInventory(inp,r%bop,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readInventory.',irc
          return
       end if
       dum250=""
       call chop0(dum250,250)
       call ncf_checkContents(inp,dum250,r%bop,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_checkContents.',irc
       end if
       if (i%lfldat(55)) then
          r%bop=.true.
       end if
       if (.not.r%bop) then
          write(*,*) myname,'missing essential information...'
          irc=845
          return
       end if
       !
       if (i%lfldat(55)) then
          ! make coordinate dimorder
          call ncf_clearDimOrder(r%ixydo) ! remove all dimensions in zdo
          r%ixydo => ncf_newDimOrder(inp,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_newDimOrder.',irc
             return
          end if
          do ii=1,i%nrcor
             lenc=length(i%cor250(ii),250,10)
             jj=ncf_getDimEntry(inp,i%cor250(ii)(1:lenc))
             xyDim => ncf_getInventoryDimension(inp,jj,irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_getInvDim.',irc
                return
             end if
             call ncf_addDimOrderDim(r%ixydo,xyDim,irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_addDimOrdDim.',irc
                return
             end if
          end do
       else
          r%ixydo => ncf_getPrimaryDimOrder(inp,irc)
          if(irc.ne.0) then
             write(*,*)myname,' Error return from ncf_getPrimaryDimOrder.',irc
             return
          end if
       end if
    end if
  end subroutine readInput
  !
  !
  ! read input file
  !
  subroutine readAuxiliary(aux,i,r,irc)
    implicit none
    type(inventory), pointer :: aux
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'readAuxiliary'/
    !
    character*250 :: dum250
    ! 
    i%lena=length(i%aux250,250,10)
    write(*,*)myname,' Scanning: ',i%aux250(1:i%lena)
    call ncf_openFile(i%aux250,aux,r%bop,irc)
    if (irc.ne.0) then
       write(*,*)myname,' Error return from ncf_openFile.',irc
       return
    end if
    if (r%aop) then
       call ncf_readInventory(aux,r%aop,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readInventory.',irc
          return
       end if
       dum250=""
       call chop0(dum250,250)
       call ncf_checkContents(aux,dum250,r%aop,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_checkContents.',irc
          return
       end if
       if (.not.r%aop) then
          write(*,*) myname,'missing essential information...'
          irc=845
          return
       end if
       r%axydo => ncf_getPrimaryDimOrder(aux,irc)
       if(irc.ne.0) then
          write(*,*)myname,' Error return from ncf_getPrimaryDimOrder.',irc
          return
       end if
    end if
  end subroutine readAuxiliary
  !
  subroutine processInput(inp,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'processInput'/
    type(attribute), pointer :: a=>null()
    real :: tanfact
    integer :: lenb
    character*250 :: buff250
    !
    if (r%bop) then
       ! get analysis time...
       r%ifrtVariable => ncf_getVariable(inp,"forecast_reference_time",r%bop,irc)
       if(irc.ne.0) then
          write(*,*)myname,' Error return from ncf_getVariable forecast_reference_time.',irc
          return
       end if
    end if
    if (r%bop) then
       call ncf_readRealData(r%ifrtVariable,r%bop,irc)
       if (.not.r%bop) irc=999
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readRealData A.',irc
          return
       end if
       !write(*,*)myname,'Forecast Reference Time:',r%ifrtVariable%fd(1)
    end if
    if (r%bop) then ! get analysis time
       r%t2000=0.0D0
       a => ncf_getAttribute(r%ifrtVariable,"units",irc)
       if (associated(a)) then
          buff250=ncf_getAttributeText(a)
          call chop0(buff250,250)
          lenb=length(buff250,250,10)
          !write(*,*)myname,"Analysis is in "//buff250(1:lenb), r%t2000
          if (buff250(1:4).eq."days") then ! convert days to seconds
             tanfact=86400.0D0
          else
             tanfact=1.0D0
          end if
          r%t2000=j2000(tanfact*ncf_valuePosition(r%ifrtVariable,irc))
          if (irc.ne.0) then
             write(*,*) myname,'Error return from valuePosition (tanid).',irc
             return
          end if
       end if
    end if
    if (r%bop) then
       r%itime=ncf_getDimEntry(inp,"time")
       if (r%itime.eq.0) r%bop=.false.
    end if
    if (r%bop) then
       r%itimeVariable => ncf_getVariable(inp,"time",r%bop,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_getVariable time.',irc
          return
       end if
    end if
    if (r%bop) then
       call ncf_readRealData(r%itimeVariable,r%bop,irc)
       if (.not.r%bop) irc=999
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readRealData B.',irc
          return
       end if
    end if
    if (r%bop) then
       !if (i%lfldat(70).and.i%lfldat(30)) then
       !write(*,*) myname,'Lonid:',associated(inp%lonid)
       call ncf_readRealData(inp%lonid,r%bop,irc)
       if (.not.r%bop) irc=999
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readRealData C.',irc
          return
       end if
       call ncf_readRealData(inp%latid,r%bop,irc)
       if (.not.r%bop) irc=999
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readRealData D.',irc
          return
       end if
       !end if
    end if
  end subroutine processInput
  !
  subroutine prepareRP(inp,out,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'prepareRP'/
    ! 
    integer :: lenp,II,lenfrom, lento
    integer YY,MM,DD,HH,MI,MM1,DD1,HH1,MI1,SS1,MM2,DD2,HH2,MI2,SS2
    real sec
    real :: ts2000,te2000
    character*250 :: cto250, cfrom250
    type(attribute), pointer :: ato=>null(), afrom=>null()
    
    if (i%lfldat(70).and.i%lfldat(30)) then
       r%tnrrp=1 ! target return-period file...
       do ii=1,i%nrrp
          lenp=length(i%rp250(r%tnrrp),250,10)
          write(*,*)myname,' Scanning: ',i%rp250(ii)(1:lenp)
          call ncf_openFile(i%rp250(ii),r%ret,r%brp,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_openFile.',irc
             return
          end if
          if (r%brp) then
             call ncf_readInventory(r%ret,r%brp,irc)
             if (r%brp) then ! check if this is the correct file...
                ato => ncf_getGlobalAttribute(r%ret,"valid_to",irc)
                afrom => ncf_getGlobalAttribute(r%ret,"valid_from",irc)
                cto250=ncf_getAttributeText(ato)
                call chop0(cto250,250)
                lento=length(cto250,250,10)
                cfrom250=ncf_getAttributeText(afrom)
                call chop0(cfrom250,250)
                lenfrom=length(cfrom250,250,10)
                ! get current analysis date
                call dj2000(r%t2000,yy,mm,dd,hh,mi,sec)
                ! get start+end rp-date
                write(*,*)myname,'From:',cfrom250(1:lenfrom),&
                     & ' to:',cto250(1:lento),' file:',i%rp250(ii)(1:lenp)
                read(cfrom250(1:lenfrom),'(I2,X,I2,X,I2,X,I2,X,I2)',iostat=irc)dd1,mm1,hh1,mi1,ss1
                sec=ss1
                call jd2000(ts2000,yy,mm1,dd1,hh1,mi1,sec)
                read(cto250(1:lento),'(I2,X,I2,X,I2,X,I2,X,I2)',iostat=irc)dd2,mm2,hh2,mi2,ss2
                sec=ss2+1.0D0
                call jd2000(te2000,yy,mm2,dd2,hh2,mi2,sec)
                ! check if analysis is in range
                if (ts2000 .lt. te2000) then
                   if (r%t2000 .ge. ts2000 .and. r%t2000 .le. te2000) then
                      write(*,'(2(X,A),X,I2.2,"/",I2.2,X,"<",I2.2,"/",I2.2,",",I2.2,"/",I2.2,">")')&
                           & myname,'Found valid s-range:',dd,mm,dd1,mm1,dd2,mm2
                      r%tnrrp=ii
                   end if
                else
                   if (.not.(r%t2000 .lt. ts2000 .and. r%t2000 .gt. te2000)) then
                      write(*,'(2(X,A),X,I2.2,"/",I2.2,X,"<",I2.2,"/",I2.2,",",I2.2,"/",I2.2,">")')&
                           & myname,'Found valid w-range:',dd,mm,dd1,mm1,dd2,mm2
                      r%tnrrp=ii
                   end if
                end if
             end if
          end if
          call ncf_closeFile(r%ret,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_closeFile.',irc
             return
          end if
          call ncf_clearInventory(r%ret,irc) ! delete internally stored data...
          if (irc.ne.0) then
             write(*,*) myname,"Error return from ncf_clearInventory.",irc
             return
          end if
       end do
       r%brp=.true.
       lenp=length(i%rp250(r%tnrrp),250,10)
       write(*,*)myname,' Scanning: ',i%rp250(r%tnrrp)(1:lenp)
       call ncf_openFile(i%rp250(r%tnrrp),r%ret,r%brp,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_openFile.',irc
          return
       end if
       if (r%brp) then
          call ncf_readInventory(r%ret,r%brp,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_readInventory.',irc
             return
          end if
       else
          irc=845
          write(*,*)myname,' Corrupt file: ',i%rp250(r%tnrrp)(1:lenp)
          return
       end if
    else
       r%brp=.false. ! no return period information...
    end if

  end subroutine prepareRP
  !
  subroutine prepareSlice(inp,out,aux,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(inventory), pointer :: aux
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'prepareSlice'/
    integer :: ii
    !
    !call ncf_printinventory(out)
    !
    ! make slices (in latid)
    !
    r%dslice=>ncf_getDimensionOrderDimension(r%ixydo,1)
    if (r%bop.and.associated(r%dslice)) then
       r%nslice=r%dslice%lim
       r%dnslice=max(1,1+int(float(r%nslice)/float(i%islice)))
    else
       r%nslice=1
       r%dnslice=1
    end if
    r%nll=ceiling(float(r%nslice)/float(r%dnslice))
    r%cll=0
    r%initWrite=.true.
    if (i%lfldat(15).and.r%aop) then
       ii=ncf_getDimEntry(aux,inp%dim250(r%dslice%ind)(1:inp%lend(r%dslice%ind)))       
       r%aslice=>ncf_makeDimension(aux,ii)
    end if
    return
  end subroutine prepareSlice
  !
  subroutine initSlice(inp,out,aux,i,r,ll,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(inventory), pointer :: aux
    type(init) :: i
    type(run) :: r
    integer :: irc
    integer :: ll
    CHARACTER*14 MYNAME
    DATA MYNAME /'initSlice'/
    integer :: lmax
    r%cll=r%cll+1
    lmax=min(r%nslice,ll+r%dnslice-1)
    write(*,'(X,A,4(A,I0))') myname,&
         & '>>>>>>>>>>>>>>>>>>> SLICE <<<<<<<<<<<<<<<<< [',&
         & ll,",",lmax,"] ",r%cll,"|",r%nll
    call ncf_setSlice(inp,r%dslice,ll,lmax)
    call ncf_setSlice(out,r%dslice,ll,lmax)
    if (i%lfldat(15).and.r%aop) then
       call ncf_setSlice(aux,r%aslice,ll,lmax)
    end if
    return
  end subroutine initSlice
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! process wind x,y -> spd
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine addWindSpd(inp,out,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'addWindSpd'/
    logical :: btmp
    integer :: ss,lenw,lens, lx, ly
    real :: vx, vy, vv
    type(variable), pointer  :: iwxVariable=>null(), iwyVariable=>null(),o=>null()
    type(dimensionOrder), pointer  ::  idx=>null()
    do ss=1,i%nrspd
       call chop0(i%wx250(ss),250)
       lenw=length(i%wx250(ss),250,10)
       !
       if(bdeb)write(*,*)myname,' Processing wind x-spd: ',i%wx250(ss)(1:lenw)
       !
       iwxVariable => ncf_getVariable(inp,i%wx250(ss)(1:lenw),r%bok,irc)
       if (.not.r%bok) irc=702
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_getVariable wx.',irc
          return
       end if
       call ncf_clearDimOrder(idx)
       idx => ncf_makeDimOrder(iwxVariable,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_makeDimOrder.',irc
          return
       end if
       r%bok=.true.
       call ncf_readRealData(iwxVariable,r%bok,irc)
       if (.not.r%bok) irc=999
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readRealData E.',irc
          return
       end if
       !
       call chop0(i%wy250(ss),250)
       lenw=length(i%wy250(ss),250,10)
       if(bdeb)write(*,*)myname,' Processing wind y-spd: ',i%wy250(ss)(1:lenw)
       !
       iwyVariable => ncf_getVariable(inp,i%wy250(ss)(1:lenw),r%bok,irc)
       if (.not.r%bok) irc=702
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_getVariable wy.',irc
          return
       end if
       r%bok=.true.
       call ncf_readRealData(iwyVariable,r%bok,irc)
       if (.not.r%bok) irc=999
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readRealData F.',irc
          return
       end if
       !
       lens=length(i%spd250(ss),250,10)
       o => ncf_getVariable(out,i%spd250(ss)(1:lens),btmp,irc)
       if (.not.btmp) then
          o => ncf_copyVariable(iwxVariable,irc) ! the new ouput variable
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_copyVariable.',irc
             return
          end if
          write(*,*)myname,' Creating: ',i%spd250(ss)(1:lens)
          call ncf_setVariableName(o,i%spd250(ss))
          call ncf_setTextAttribute(o,"standard_name",i%spd250(ss)(1:lens))
          if(bdeb)write(*,*)myname,"### Prepending A:'"//o%var250(1:o%lenv)//"'",o%lend
          call ncf_prependVariable(out,o)
       end if
       call ncf_initField(o,o%filld,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_initField.',irc
          return
       end if
       ! make sure positions are the same...
       call ncf_importVariable(inp,o,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       !
       call ncf_resetPos(idx,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_resetPos.',irc
          return
       end if
       WGRID: do while (ncf_increment(inp,idx,irc))
          lx = ncf_getLocation(iwxVariable)
          ly = ncf_getLocation(iwyVariable)
          vx=iwxVariable%fd(lx)
          vy=iwyVariable%fd(ly)
          if (vx.ne.iwxVariable%filld .and. vy.ne.iwyVariable%filld) then
             vv=dsqrt(vx*vx+vy*vy)
             if (vv.le.i%vvmax(ss)) then
                o%fd(lx)=vv
             end if
          end if
       end do WGRID
       !
       call ncf_importVariable(out,o,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       if (bdeb) then
          call ncf_checkInventory(inp,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_checkInventory A.',irc
             return
          end if
          if(bdeb)write(*,*)myname,' Size:',o%var250(1:o%lenv),o%lend
       end if
       !
    end do
    nullify(o)
  end subroutine addWindSpd
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! process wind gust
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine addWindGust(inp,out,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'addWindGust'/
    type(variable), pointer :: g=>null(),o=>null()
    integer :: leng, lend,gusEntry, li, lo
    type(dimensionOrder), pointer  ::  odo=>null()
    type(dimension), pointer  ::  gusDim=>null()
    !
    call chop0(i%guv250,250)
    leng=length(i%guv250,250,10)
    !
    if(bdeb)write(*,*)myname,' Processing wind gust: ',i%guv250(1:leng)
    !
    g => ncf_getVariable(inp,i%guv250(1:leng),r%bok,irc)
    if (.not.r%bok) irc=702
    if (irc.ne.0) then
       write(*,*)myname,' Error return from ncf_getVariable wx.',irc
       return
    end if
    r%bok=.true.
    call ncf_readRealData(g,r%bok,irc)
    if (.not.r%bok) irc=999
    if (irc.ne.0) then
       write(*,*)myname,' Error return from ncf_readRealData G.',irc
       return
    end if
    gusEntry=ncf_getDimEntry(inp,i%gud250)
    call ncf_clearDimOrder(odo)
    odo => ncf_makeDimOrder(g,irc)
    gusDim => ncf_removeDimOrderEntry(odo,gusEntry)
    !
    if (associated(gusDim)) then ! average necessary
       o => ncf_copyVariable(g,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_copyVariable.',irc
          return
       end if
       call ncf_setVariableDimOrder(o,odo,irc)
       if (irc.ne.0) return
       leng=length(i%gus250,250,10)
       write(*,*)myname,' Creating: ',i%gus250(1:leng)
       call ncf_setVariableName(o,i%gus250)
       call ncf_setTextAttribute(o,"standard_name",i%gus250(1:leng))
       call ncf_initField(o,o%filld,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_initField.',irc
          return
       end if
       call ncf_setDimensionValue(inp,gusDim,i%gul)       ! forecast time in output-file
       call ncf_resetPos(odo,irc)
       do while (ncf_increment(inp,odo,irc))
          li = ncf_getLocation(g) ! location in output array
          lo = ncf_getLocation(o) ! location in output array
          o%fd(lo)=g%fd(li)
       end do
    else
       lend=length(i%gud250,250,10)
       write(*,*)myname,' Missing dimension: ',i%gud250(1:lend), &
            & ' in Gust-variable: '//i%guv250(1:leng)
    end if
    !
    if (bdeb) then
       call ncf_checkInventory(inp,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_checkInventory A.',irc
          return
       end if
       write(*,*)myname,' Size:',o%var250(1:o%lenv),o%lend
    end if
    !
    if(bdeb)write(*,*)myname,"### Prepending B:'"//o%var250(1:o%lenv)//"'",o%lend
    call ncf_prependVariable(out,o)
    nullify(o)
    if (bdeb) then
       call ncf_checkInventory(out,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_checkInventory C.',irc
          return
       end if
    end if
    return
  end subroutine addWindGust
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! process differentiated variables
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine addDifferentiated(inp,out,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*18 MYNAME
    DATA MYNAME /'addDifferentiated'/
    type(dimensionOrder), pointer  ::  ido=>null()
    integer ihrs, ik,jk, li,lj,lk,lo, lenc
    real :: fact,hrs,trg,val
    character*10 chrs10
    type(variable), pointer :: v=>null(),n=>null(),o=>null()
    if(bdeb)write(*,*)myname,' Processing differentiated variables.',i%nrdiff
    LHRS: do ihrs=1,i%nrdiff
       n => inp%firstVariable%next
       LDIFF: do while (.not.associated(n,target=inp%lastVariable))
          v => n
          if(bdeb)write(*,*) myname,' Processing diff: ',v%lend,v%var250(1:v%lenv)
          n => n%next
          if (.not.i%diff250(ihrs)(1:i%lendiff(ihrs)).eq.v%var250(1:v%lenv))cycle LDIFF ! not target variable
          if (.not.ncf_variableContainsDim(v,r%itime)) cycle LDIFF ! no time dimension
          !
          if (.not.v%readyfield) then
             if(bdeb)write(*,*) myname,' Reading diff: ',v%var250(1:v%lenv),v%lend,associated(v%fd),r%bok
             r%bok=.true.
             call ncf_readRealData(v,r%bok,irc)
             if (.not.r%bok) irc=999
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_readRealData H.',irc,v%lend
                return
             end if
             !!write(*,*) myname,' Read: ',v%var250(1:v%lenv),v%lend,associated(v%fd)
          end if
          !
          hrs=max(1.0D-10,abs(i%hdiff(ihrs)))
          write(chrs10,'(F10.2)')hrs
          call chop0(chrs10,10)
          lenc=length(chrs10,10,10)
          write(*,*)myname,' Differentiating: ',v%var250(1:v%lenv),v%lend
          !
          ido => ncf_makeDimOrder(v,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_makeDimOrder.',irc
             return
          end if
          r%iTimeDim => ncf_removeDimOrderEntry(ido,r%itime)
          !
          o => ncf_getVariable(out,i%niff250(ihrs)(1:i%lenniff(ihrs)),r%bok,irc)
          if (r%bok) then
             write(*,*)myname," Found '"//i%niff250(ihrs)(1:i%lenniff(ihrs))//"'"
          else
             write(*,*)myname," Creating '"//i%niff250(ihrs)(1:i%lenniff(ihrs))//"'"
          end if
          !
          ! make variable
          !
          if (.not.r%bok) then ! make import variable
             ! create differentiated output variable
             o => ncf_copyVariable(v,irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_copyVariable.',irc
                return
             end if
             ! copy coordinate attributes
             call ncf_copyVariableAttribute(o,v,"grid_mapping",irc)
             call ncf_copyVariableAttribute(o,v,"coordinates",irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_copyVariableAttribute.',irc
                return
             end if
             !
             call ncf_setTextAttribute(o,"accumulation_hours",chrs10(1:lenc))
             ! initialise field to "undefined"
             ! set variable name
             call ncf_setVariableName(o,i%niff250(ihrs))
             call ncf_setTextAttribute(o,"standard_name",i%niff250(ihrs)(1:i%lenniff(ihrs)))
             !
             ! add output variable
             if(bdeb)write(*,*)myname,"### Prepending C:'"//o%var250(1:o%lenv)//"'",o%lend
             call ncf_prependVariable(out,o)
          end if
          !
          call ncf_importVariable(inp,o,irc)
          if (irc.ne.0) then
             write(*,*)myname,&
                  & ' Error return from ncf_importvariable.',irc
             return
          end if
          !
          call ncf_initField(o,o%filld,irc)
          if (irc.ne.0) then
             write(*,*)myname,&
                  & ' Error return from ncf_initField.',irc
             return
          end if
          ! loop over times
          LITIMELOOP: do ik=2,r%iTimeVariable%lend ! last i-time index
             !if (bdeb) write(*,*)myname,' I-time-loop start',ik
             ! search for start-time...
             jk=1
             trg=r%iTimeVariable%fd(ik)-hrs*3600.0D0
             fact=0.0D0
             LJTIMELOOP: do while (jk.lt.ik)! find j-time index (hrs earlier than i-time)
                if (r%iTimeVariable%fd(jk).le.trg.and.r%iTimeVariable%fd(jk+1).gt.trg) then
                   fact=(trg-r%iTimeVariable%fd(jk))/(r%iTimeVariable%fd(jk+1)-r%iTimeVariable%fd(jk))
                   exit LJTIMELOOP
                end if
                jk=jk+1
             end do LJTIMELOOP
             !if (bdeb) write(*,*)myname,' J-time',ik,r%iTimeVariable%fd(ik),trg,jk
             if (jk.lt.ik) then ! did we find a valid start-time?
                if(bdeb)write(*,*) myname,' Found diff-time-match for:',&
                     & v%var250(1:v%lenv),ik,jk,r%iTimeVariable%lend,associated(v%fd)
                !
                call ncf_resetPos(ido,irc)
                !
                do while (ncf_increment(inp,ido,irc))
                   !
                   call ncf_setDimensionValue(inp,r%iTimeDim,ik)       ! forecast time in input-file
                   li = ncf_getLocation(v) ! location in input array
                   lo = ncf_getLocation(o) ! location in input array
                   !
                   call ncf_setDimensionValue(inp,r%iTimeDim,jk)       ! forecast time in output-file
                   lj = ncf_getLocation(v) ! location in output array
                   !
                   call ncf_setDimensionValue(inp,r%iTimeDim,jk+1)       ! forecast time in output-file
                   lk = ncf_getLocation(v) ! location in output array
                   !
                   !if(bdeb)write(*,*)myname,' Inner loop:',v%var250(1:v%lenv),li,lj,v%lend
                   if ( v%fd(li).ne.v%filld .and. &
                        & v%fd(lj).ne.v%filld .and. &
                        & v%fd(lk).ne.v%filld ) then ! we have a real value
                      val = (v%fd(li)-(v%fd(lk)*fact+(1.0D0-fact)*v%fd(lj)))! /hrs: difference, not derivative...
                      !if (val.lt.0.0D0) then
                      !   write(*,*)myname,'Strange:',ik,v%fd(li),jk,v%fd(lj),v%fd(lk),fact,val
                      !endif
                      if (i%vset(ihrs)) then
                         val=max(min(val,&
                              & max(i%vmin(ihrs),i%vmax(ihrs))),&
                              & min(i%vmax(ihrs),i%vmin(ihrs)))
                      end if
                   else
                      val = o%filld
                   end if
                   !write(*,*)myname,' Inner loop:',v%var250(1:v%lenv),ik,jk,li,lj,v%fd(li),v%fd(lj),val
                   o%fd(lo) = val
                   !if (o%fd(li).lt.0.0D0) write(*,*)myname,' Sanity test:',li,o%fd(li),val
                end do
                !if (bdeb) write(*,*)myname,' I-time-loop end'
             else
                write(*,*)myname,' No diff-match for: ',v%var250(1:v%lenv),ik,r%iTimeVariable%lend,&
                     & r%iTimeVariable%fd(ik),hrs
             end if
          end do LITIMELOOP
          ! ! sanity check
          ! if (bdeb) then
          ! write(*,*)myname,"Starting Sanit check '"//o%var250(1:o%lenv)//"'"
          ! do ii=1,o%lend
          !    if (o%fd(ii).lt.0.0D0) then
          !       write(*,*)myname,'Insane value:',ii,o%fd(ii)
          !    end if
          ! end do
          ! write(*,*)myname,'Sanit check done.'
          ! end if
          call ncf_importVariable(inp,o,irc)
          if (irc.ne.0) then
             write(*,*)myname,&
                  & ' Error return from ncf_importvariable.',irc
             return
          end if
          !
          nullify(o)
          call ncf_clearDimOrder(ido)
          if(irc.ne.0)return
          ! exit LDIFF
       end do LDIFF
    end do LHRS
  end subroutine addDifferentiated
  !
  subroutine addReturnPeriod(inp,out,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*16 MYNAME
    DATA MYNAME /'addReturnPeriod'/
    !
    integer :: nretper=0
    real, allocatable :: aretper(:)
    integer, allocatable :: aretind(:),aretmax(:),aretmin(:)
    integer :: dbg, ii,jj,left,right
    type(weight), pointer :: wgt => null()   ! weight of surrounding grid points
    !
    integer :: ivar,lz,lr,li,lo
    logical :: brok, biok, binp, maxevent,increasing !
    type(variable), pointer  :: o=>null(),v=>null(),zi=>null(),rp=>null(),ri=>null()
    type(varptr), pointer  :: impVariable(:)=>null()
    type(dimensionOrder), pointer  ::  ido=>null()
    integer :: nn=0,newnn
    integer :: ook(2)
    integer :: orm(2)
    real :: drf,pst(2)
    logical :: first,valid_first
    real :: valid_min,valid_max,lat,lon
    logical :: bbrp
    character*250 :: unit250
    integer :: lenu=0,lenr, iz,lenii, lenp
    type(dimension), pointer :: rx=>null() ! return level x-dimension
    type(dimension), pointer :: ry=>null() ! return level y-dimension
    type(dimensionOrder), pointer  ::  rdo=>null(),zdo=>null(),xydo=>null()
    logical :: boz,bozz,found
    !
    dbg=0
    !
    ! open return period file and read inventory
    !
    if (r%brp) then
       bbrp=r%brp
       lenr=length(i%rpd250,250,10)
       call ncf_checkParid(r%ret,i%rpd250,bbrp,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_checkContents.',irc
          return
       end if
       !
       ! read latitude and longitude into memory
       !
       call ncf_readRealData(r%ret%lonid,r%bop,irc)
       if (.not.r%bop) irc=999
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readRealData J.',irc
          return
       end if
       call ncf_readRealData(r%ret%latid,r%bop,irc)
       if (.not.r%bop) irc=999
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readRealData K.',irc
          return
       end if
       ! read return-period variable into memory
       call ncf_readRealData(r%ret%parid,bbrp,irc)
       if (.not.bbrp) irc=999
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readRealData I.',irc,v%lend,&
               & " par='"//i%rpd250(1:lenr)//"'"
          return
       end if
       write(*,'(X,A,A,100(X,F0.3))')myname,'Return periods ('//i%rpd250(1:lenr)//'):',r%ret%parid%fd
       unit250 = ncf_getTextAttribute(r%ret%parid,"units",irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_getAttribute.',irc
          return
       end if
       call chop0(unit250,250)
       lenu=length(unit250,250,10)
       if(bdeb)write(*,*)myname,' Found return-parameter unit:',unit250(1:lenu)
       ! read import variables into memory...
       if (.not. associated(impVariable)) then
          allocate(impVariable(maximp),stat=irc)
       end if
       ! read import variables into memory...
       do ii=1,i%nrimp
          r%bok=.true.
          lenii=length(i%imp250(ii),250,10)
          impVariable(ii)%ptr=>ncf_getVariable(r%ret,i%imp250(ii)(1:lenii),r%bok,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_getVariable.',irc
             return
          end if
          if (associated(impVariable(ii)%ptr)) then
             if (.not. impVariable(ii)%ptr%readyfield) then
                call ncf_readRealData(impVariable(ii)%ptr,r%bok,irc)
                if (.not.r%bok) irc=999
                if (irc.ne.0) then
                   write(*,*)myname,' Error return from ncf_readRealData KI.',irc
                   return
                end if
                write(*,*)myname," Importing: "//i%imp250(ii)(1:lenii),impVariable(ii)%ptr%lend
             else
                write(*,*)myname," Already imported: "//i%imp250(ii)(1:lenii)
             end if
          else
             write(*,*)myname," Import not found: "//i%imp250(ii)(1:lenii)
             irc=845
             return
          end if
       end do
    end if
    !
    call ncf_clearDimOrder(xydo)
    if (bbrp) then
       ! x,y dimensions
       xydo => ncf_getPrimaryDimOrder(r%ret,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_getPrimaryDimensions L.',irc
          return
       end if
       rx=>ncf_getDimensionOrderDimension(xydo,1)
       ry=>ncf_getDimensionOrderDimension(xydo,2)
       if (bdeb) then
          write(*,*)myname,'Return period xy-dimorder:'
          call ncf_printDimOrder(xydo)
       end if
       if (bdeb) then
          call ncf_printDimension(r%ret,rx)
          call ncf_printDimension(r%ret,ry)
       end if
    else
       irc=341
       WRITE(*,*) MYNAME," Missing XY-navigation in RetPer."
       RETURN
    end if
    !
    ! import variables....
    !
    if (r%brp) then
       ! make import variables in output file...
       do ii=1,i%nrimp
          ! check if variable already exists
          r%bok=.false. ! we dont expect to find this variable (silent)
          lenii=length(i%imp250(ii),250,10)
          zi => ncf_getVariable(out,i%imp250(ii)(1:lenii),r%bok,irc)
          if (r%bok) then
             write(*,*)myname," Variable exists '"//i%imp250(ii)(1:lenii)//"'"
          else
             write(*,*)myname," Interpolating '"//i%imp250(ii)(1:lenii)//"'"
          end if
          !
          ! make import-variable in output file...
          !
          if (.not.r%bok) then ! make import variable
             !write(*,*)myname," Interpolating '"//i%imp250(ii)(1:lenii)//"'"
             ri => impVariable(ii)%ptr;
             zi => ncf_copyVariable(ri,irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_copyVariable.',irc
                return
             end if
             call ncf_clearDimOrder(zdo)
             zdo=>ncf_newDimOrder(inp,irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_newDimOrder.',irc
                return
             end if
             call ncf_addDimOrder(zdo,r%ixydo,irc) ! only position dims...
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_addDimOrder.',irc
                return
             end if
             call ncf_importDimOrder(out,zdo,irc)
             if (irc.ne.0) then
                write(*,*)myname,&
                     & ' Error return from ncf_importdimorder.',irc
                return
             end if
             call ncf_setVariableDimOrder(zi,zdo,irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_setVariableDimOrder.',irc
                return
             end if
             if(bdeb)write(*,*)myname,"### Prepending D:'"//zi%var250(1:zi%lenv)//"'",zi%lend
             call ncf_prependVariable(out,zi)
          end if
          !
          ! make sure positions are the same...
          call ncf_importVariable(inp,zi,irc)
          if (irc.ne.0) then
             write(*,*)myname,&
                  & ' Error return from ncf_importvariable.',irc
             return
          end if
          !
          ! initialise field to "undefined"
          call ncf_initField(zi,zi%filld,irc)
          if (irc.ne.0) then
             write(*,*)myname,&
                  & ' Error return from ncf_initField.',irc
             return
          end if
          !
          ! get lat/lon dimensions
          !
          ! copy import variable to output file...
          !
          wgt => ncf_makeWeight4(r%ret%nrdim,irc)
          ! loop over output-grid
          call ncf_resetPos(r%ixydo,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_resetPos.',irc
             return
          end if
          LLIMP: do while (ncf_increment(inp,r%ixydo,irc))
             ! interpolate...
             biok=.true.
             call ncf_interpolate2D(inp,r%ret,rx,ry,wgt,biok,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from interpolate2d.',irc
                return
             end if
             ! assign...
             if (biok) then
                lz = ncf_getLocation(zi) ! location in input array
                lr = ncf_getLocation(ri) ! location in input array
                zi%fd(lz)=ncf_valueWeighted(ri,wgt,irc)
                !if (0 .eq. mod(lz,1000)) write(*,*)myname,'Lz:',lz,zi%fd(lz),&
                !     & lr,r%fd(lr)
             end if
          end do LLIMP
          !call ncf_printVariable(zi)
          call ncf_importVariable(out,zi,irc)
          if (irc.ne.0) then
             write(*,*)myname,&
                  & ' Error return from ncf_importvariable.',irc
             return
          end if
          call ncf_clearWeight(wgt)
          call ncf_clearDimOrder(zdo)
          nullify(zi)
          !end if ! import variable
       end do ! nrimp
    end if
    !
    ! make return period variables...
    !
    if (r%brp) then
       lenp=length(i%rp250(r%tnrrp),250,10)
       write(*,*)myname,' processing: ',i%rp250(r%tnrrp)(1:lenp),r%brp
       ook(1)=0
       ook(2)=0
       orm(1)=0
       orm(2)=0
       wgt => ncf_makeWeight4(r%ret%nrdim,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_makeWeight4.',irc
          return
       end if
       !
       ! make dimorder for retper
       rdo=>ncf_makeDimOrder(r%ret%parid,irc)
       nretper=ncf_getDimOrderLength(rdo) ! number of return periods
       nn=nretper                         ! actual number of return periods
       increasing=(r%ret%parid%fd(1).le.r%ret%parid%fd(nretper))
       !we assume that the return periods are sorted...
       if (allocated(aretper)) deallocate(aretper)
       if (allocated(aretind)) deallocate(aretind)
       if (allocated(aretmax)) deallocate(aretmax)
       if (allocated(aretmin)) deallocate(aretmin)
       allocate(aretper(nretper),aretind(nretper),aretmax(nretper),&
            & aretmin(nretper),stat=irc)
       if(irc.ne.0) then
          write(*,*)myname,' Unable to allocate ARETPER.',irc
          return
       end if
       ! The return_levels can be for either "minimum events" or "maximum events"
       ! If the return_level increases with increasing return_period, it is a "maximum event"...
       do ii=1,nn
          aretind(ii)=ii ! return_period index
       end do
       call sort_heapsort1r(nretper,r%ret%parid%fd,eps,newnn,nn,aretind,uniq)
       do ii=1,newnn                        ! index to valid return periods
          aretmax(ii)=aretind(ii)           ! maxevent index
          aretmin(ii)=aretind(newnn-ii+1)   ! minevent index
       end do
       !
       ! add return-period variable to output file
       if (i%lfldat(35)) then ! create return-level variable
          ! check if return-level variable already exists
          r%bok=.false. ! we dont expect to find this variable (silent)
          zi => ncf_getVariable(inp,i%rpd250(1:lenr),r%bok,irc)
          if (r%bok) then
             write(*,*)myname," Variable exists '"//i%rpd250(1:lenr)//"'"
          else
             write(*,*)myname," Interpolating '"//i%rpd250(1:lenr)//"'"
          end if
          bozz=(.not.r%bok) ! should we make return-level variable?
       else
          bozz=.false.
       end if
       if (bozz) then ! make return-level variable
          !write(*,*)myname," Interpolating '"//i%rpd250(1:lenr)//"'"
          zi =>ncf_copyVariable(r%ret%parid,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_copyVariable.',irc
             return
          end if
          call ncf_copyVariableField(zi,r%ret%parid,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_copyVariableField.',irc
             return
          end if
          !
          ! make return-period dimension
          r%ro=>ncf_createDimension(out,i%rpd250(1:lenr),nretper)
          r%rz=>ncf_createDimension(inp,i%rpd250(1:lenr),nretper)
          iz=ncf_getDimEntry(inp,i%rpd250(1:lenr))
          write(*,*)myname,"Created dimension '"//i%rpd250(1:lenr)//"'",r%rz%ind
          !
          call ncf_clearDimOrder(zdo) ! remove all dimensions in zdo
          ! add return-period dimension (r%rz)
          zdo=>ncf_newDimOrder(inp,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_newDimOrder.',irc
             return
          end if
          call ncf_addDimOrderEntry(zdo,iz,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_DimOrderDim.',irc
             return
          end if
          call ncf_setVariableDimOrder(zi,zdo,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_DimOrderDim.',irc
             return
          end if
          call ncf_importVariable(out,zi,irc)
          if (irc.ne.0) then
             write(*,*)myname,&
                  & ' Error return from ncf_importvariable.',irc
             return
          end if
          !
          if(bdeb)write(*,*)myname,"### Prepending E:'"//zi%var250(1:zi%lenv)//"'",zi%lend
          call ncf_prependVariable(inp,zi)
          nullify(zi)
          call ncf_clearDimOrder(zdo)
          ! at this point the return-level dimension RZ is defined...
       end if
       if (bdeb) then
          call ncf_checkInventory(inp,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_checkInventory B.',irc
             return
          end if
       end if
       ! loop over output return variables
       if(bdeb)write(*,*)myname,' Return periods: ',i%nrret
       LREPER: do ivar=1,i%nrret
          first=.true.
          if(bdeb)WRITE(*,*) MYNAME," Creating return-period: '"//i%net250(ivar)(1:i%lennet(ivar))//"'"
          ! initialise return period...
          brok=.true.
          rp => ncf_getVariable(r%ret,i%per250(ivar)(1:i%lenper(ivar)),brok,irc)
          IF (IRC.NE.0) THEN
             WRITE(*,*) MYNAME,' Error return from ncf_getVariable.',IRC
             RETURN
          END IF
          if (.not.brok) then
             irc=345
             WRITE(*,*) MYNAME," Missing return-level variable '"//i%per250(ivar)(1:i%lenper(ivar))//"'"
             RETURN
          end if
          if (.not.rp%readyfield) then ! read data into memory if not already there
             call ncf_readRealData(rp,brok,irc)
             if (.not.brok) irc=999
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_readRealData L.',irc,v%lend
                return
             end if
          end if
          call ncf_clearDimOrder(xydo)
          if (brok) then
             ! x,y dimensions
             xydo => ncf_getPrimaryDimOrder(r%ret,irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_getPrimaryDimensions L.',irc
                return
             end if
             rx=>ncf_getDimensionOrderDimension(xydo,1)
             ry=>ncf_getDimensionOrderDimension(xydo,2)
             if (bdeb) then
                write(*,*)myname,'Return period xy-dimorder:'
                call ncf_printDimOrder(xydo)
             end if
             if (bdeb) then
                call ncf_printDimension(r%ret,rx)
                call ncf_printDimension(r%ret,ry)
             end if
          else
             irc=345
             WRITE(*,*) MYNAME," Unable to read variable: '"//i%per250(ivar)(1:i%lenper(ivar))//"'"
             RETURN
          end if
          !
          found=.false.
          v => ncf_getVariable(inp,i%ret250(ivar)(1:i%lenret(ivar)),binp,irc)
          if (irc.ne.0) then
             binp=.false.
          end if
          if (.not.binp) then
             v => ncf_getVariable(out,i%ret250(ivar)(1:i%lenret(ivar)),r%bok,irc)
             ! make sure positions are the same...
             call ncf_importVariable(inp,v,irc)
             if (irc.ne.0) then
                write(*,*)myname,&
                     & ' Error return from ncf_importvariable.',irc
                return
             end if
             if (.not.r%bok) irc=701
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_getVariable...',&
                     & irc,i%ret250(ivar)(1:i%lenret(ivar))
                return
             end if
          else
             r%bok=.true.
          end if
          if (associated(v)) then
             !n => inp%firstVariable%next
             !LRET: do while (.not.associated(n,target=inp%lastVariable))
             !   v => n
             !   n => n%next
             !   if (.not.i%ret250(ivar)(1:i%lenret(ivar)).eq.v%var250(1:v%lenv)) cycle LRET ! not target variable
             !   if (.not.ncf_variableContainsDim(v,itime)) cycle LRET ! no time dimension
             !
             found=.true.
             if(bdeb)write(*,*) myname,' Return-variable: ',v%var250(1:v%lenv),v%lend
             if (.not.v%readyfield) then ! read data into memory if not already there
                brok=.true. ! always true...
                call ncf_readRealData(v,brok,irc)
                if (.not.brok) irc=999
                if (irc.ne.0) then
                   write(*,*)myname,' Error return from ncf_readRealData M.',irc,v%lend
                   return
                end if
             end if
             !
             call ncf_clearDimOrder(ido)
             ido=>ncf_makeDimOrder(v,irc)
             call ncf_setInventoryDimOrder(inp,ido,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from setInventoryDimOrder.',irc
                return
             end if
             ! remove coordinate dimensions from variable dim-order
             call ncf_removeDimOrder(ido,r%ixydo)
             !
             ! get any old variable that has the same name as output variable...
             !
             r%bok=.false. ! we dont expect to find this variable (silent)
             !
             ! create differentiated output variable
             !
             o => ncf_getVariable(out,i%net250(ivar)(1:i%lennet(ivar)),r%bok,irc)
             if (.not.r%bok) then
                o => ncf_copyVariable(v,irc)
                if (irc.ne.0) then
                   write(*,*)myname,' Error return from ncf_copyVariable.',irc
                   return
                end if
                ! set correct unit
                !write(*,*)myname," Setting unit '"//unit250(1:lenu)//"'"
                call ncf_setTextAttribute(o,"units",unit250(1:lenu))
                ! set variable name
                call ncf_setVariableName(o,i%net250(ivar))
                call ncf_setTextAttribute(o,"standard_name",i%net250(ivar)(1:i%lennet(ivar)))
                if(bdeb)write(*,*)myname,"### Prepending F:'"//o%var250(1:o%lenv)//"'",o%lend
                call ncf_prependVariable(out,o)
             end if
             ! initialise field to "undefined"
             call ncf_initField(o,o%filld,irc)
             if (irc.ne.0) then
                write(*,*)myname,&
                     & ' Error return from ncf_initField.',irc
                return
             end if
             ! make sure positions are the same...
             call ncf_importVariable(inp,o,irc)
             if (irc.ne.0) then
                write(*,*)myname,&
                     & ' Error return from ncf_importvariable.',irc
                return
             end if
             !
             ! create return-level variable
             !
             if (bozz) then ! create return-level variable
                ! check if return-level variable already exists
                r%bok=.false. ! we dont expect to find this variable (silent)
                zi => ncf_getVariable(out,i%per250(ivar)(1:i%lenper(ivar)),r%bok,irc)
                boz=(.not.r%bok) ! should we make return-level variable?
             else
                boz=.false.
             end if
             if (boz) then ! make return-level variable
                if(bdeb)write(*,*)myname," Interpolating '"//i%per250(ivar)(1:i%lenper(ivar))//"'"
                zi =>ncf_copyVariable(rp,irc)
                if (irc.ne.0) then
                   write(*,*)myname,' Error return from ncf_copyVariable.',irc
                   return
                end if
                call ncf_copyVariableAttribute(zi,v,"grid_mapping",irc)
                call ncf_copyVariableAttribute(zi,v,"coordinates",irc)
                if (irc.ne.0) then
                   write(*,*)myname,' Error return from ncf_copyVariableAttribute.',irc
                   return
                end if
                !
                call ncf_clearDimOrder(zdo)
                zdo=>ncf_copyDimOrder(r%ixydo,irc)
                if (irc.ne.0) then
                   write(*,*)myname,' Error return from ncf_copyOrderDim.',irc
                   return
                end if
                ! add return-period dimension (r%rz)
                call ncf_addDimOrderEntry(zdo,iz,irc)
                if (irc.ne.0) then
                   write(*,*)myname,' Error return from ncf_DimOrderDim.',irc
                   return
                end if
                call ncf_setVariableDimOrder(zi,zdo,irc)
                if (irc.ne.0) then
                   write(*,*)myname,' Error return from ncf_DimOrderDim.',irc
                   return
                end if
                if(bdeb)write(*,*)myname,"### Prepending G:'"//&
                     & zi%var250(1:zi%lenv)//"'",zi%lend
                call ncf_prependVariable(out,zi)
                ! initialise field to "undefined"
                call ncf_initField(zi,zi%filld,irc)
                if (irc.ne.0) then
                   write(*,*)myname,&
                        & ' Error return from ncf_initField.',irc
                   return
                end if
                ! call ncf_printVariable(zi)
                ! make sure we use inp-positions
                call ncf_importVariable(inp,zi,irc)
                if (irc.ne.0) then
                   write(*,*)myname,&
                        & ' Error return from ncf_importvariable.',irc
                   return
                end if
             end if
             !
             valid_first=.true.
             valid_max=o%filld
             valid_min=o%filld
             !write(*,*)myname,'Ret pos: '
             if (bdeb) then
                call ncf_firstDimension(r%ret,rx)
                call ncf_firstDimension(r%ret,ry)
                call ncf_printPos(r%ret%pos)
                ! print corners...
                li = ncf_getLocation(r%ret%lonid)
                lon=r%ret%lonid%fd(li)
                li = ncf_getLocation(r%ret%latid)
                lat=r%ret%latid%fd(li)
                write(*,*)myname,'Corner (1,1):',lon,lat
                call ncf_lastDimension(r%ret,rx)
                call ncf_lastDimension(r%ret,ry)
                call ncf_printPos(r%ret%pos)
                ! print corners...
                li = ncf_getLocation(r%ret%lonid)
                lon=r%ret%lonid%fd(li)
                li = ncf_getLocation(r%ret%latid)
                lat=r%ret%latid%fd(li)
                write(*,*)myname,'Corner (+,+):',lon,lat
                call ncf_printDimOrder(r%ixydo)
                ! loop over coordinates dimension-order
                call ncf_firstDimension(r%ret,rx)
                call ncf_firstDimension(r%ret,ry)
             end if
             call ncf_resetPos(r%ixydo,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_resetPos.',irc
                return
             end if
             LLGRID: do while (ncf_increment(inp,r%ixydo,irc))
                call ncf_resetPos(ido,irc)
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_resetPos.',irc
                   return
                end if
                ! interpolate positions...
                dbg=dbg+1
                biok=.true.
                !write(*,*)myname," Interpolating."
                !dbt=949*381+481
                ncf_bdeb=.false. ! (dbg.eq.dbt)
                call ncf_interpolate2D(inp,r%ret,rx,ry,wgt,biok,irc)
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from interpolate2d.',irc
                   return
                end if
                if (ncf_bdeb) then
                   write(*,*)myname," ========Input pos."
                   call ncf_printPos(inp%pos)
                   write(*,*)myname," ========Ret pos."
                   call ncf_printPos(r%ret%pos)
                   write(*,*)myname," ========Get Retper.",biok
                   !if (dbg.eq.dbt) return
                end if

                call ncf_resetPos(rdo,irc)
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_resetPos.',irc
                   return
                end if
                if (biok) then
                   jj=0
                   LRETPER: do while (ncf_increment(r%ret,rdo,irc))
                      jj=jj+1
                      !call ncf_printPos(rp%f%pos)
                      !call ncf_printWeight(wgt)
                      aretper(jj)=ncf_valueWeighted(rp,wgt,irc)
                      if (aretper(jj).eq.rp%filld) then
                         biok=.false. ! undefined return_levels detected
                         exit LRETPER
                      end if
                      if (boz) then! set return-period variable
                         call ncf_setDimensionValue(inp,r%rz,jj)       ! forecast time in input-file
                         lz = ncf_getLocation(zi) ! location in input array
                         zi%fd(lz)=aretper(jj)
                         if (ncf_bdeb) then
                            write(*,*)'Retper:',jj,lz,zi%fd(lz)
                         end if
                      end if
                   end do LRETPER
                end if
                !write(*,*)myname," Sub loop.",biok
                if (biok) then
                   ! loop over other dimensions...
                   LRGRID: do while (ncf_increment(inp,ido,irc))
                      if (bdeb)then
                         write(*,*)myname,'Grid loop...',biok
                         call ncf_printPos(inp%pos)
                         call ncf_printDimOrder(ido)
                      end if
                      !
                      ! loop over return-period dimension
                      !
                      left=-1
                      right=-1
                      li = ncf_getLocation(v)
                      lo = ncf_getLocation(o)
                      biok=(v%fd(li).ne.v%filld)
                      if (biok) then
                         if (jj.ne.nretper) then
                            irc=956
                            write(*,*)myname,'Invalid NRETPER:',jj,nretper
                            return
                         end if
                         if (jj.gt.1) then
                            maxevent=(aretper(1).le.aretper(nretper))
                         else if (jj.eq.1) then ! assume maxevent
                            maxevent=.true.
                         else ! no information..
                            irc=845
                            write(*,*)myname,'No return-information available.',irc
                            return
                         end if
                         ! store in output grid
                         first=.false.
                         !write(*,*)myname," Store.",maxevent,li,lo
                         if (maxevent) then        ! maximum event
                            if (v%fd(li).lt.aretper(aretmax(1))) then          ! below lowest return_level
                               if (i%lfldat(31)) then ! we have zero level
                                  if (aretper(aretmax(1)).ge.0.0D0.and.v%fd(li).ge.i%zero) then
                                     right=aretmax(1)
                                     drf=max(1.0D-10, aretper(right))
                                     o%fd(lo)=((v%fd(li)-i%zero)/drf) * r%ret%parid%fd(right)
                                     if (v%fd(li).gt.i%zero.and.o%fd(lo).lt.1.0D-10) then
                                        write(*,*) myname,'Zero:',v%fd(li),i%zero,drf," rp:",r%ret%parid%fd," ps:",aretper
                                     end if
                                  else
                                     !write(*,*)myname,'Out of range:',aretper(aretmax(1)),&
                                     !     &li,v%fd(li),lo,associated(o%fd)
                                     o%fd(lo)=o%filld
                                  end if
                               else
                                  o%fd(lo)=o%filld
                               end if
                            else if (v%fd(li).gt.aretper(aretmax(newnn))) then ! above
                               left=1
                               right=newnn
                               left=aretmax(left) 
                               right=aretmax(right)
                               o%fd(lo)=extrapolate(aretper(left),r%ret%parid%fd(left),&
                                    & aretper(right),r%ret%parid%fd(right), &
                                    & v%fd(li))
                               ! o%fd(lo)=r%ret%parid%fd(aretmax(newnn))
                            else                                               ! between
                               call sort_heapsearch1r(nretper,aretper,eps,newnn,aretmax,v%fd(li),left,right)
                               if (right.eq.0) then
                                  write(*,*)myname," Heapsearch1r system error:'",v%var250(1:v%lenv)//"'",&
                                       & v%fd(li),newnn,left,right,eps
                                  write(*,*)myname," aretper:",aretper
                                  write(*,*)myname," aretper:",aretmax
                                  call ncf_printPos(v%f%pos)
                                  irc=911
                                  return
                               end if
                               left=aretmax(left) 
                               right=aretmax(right)
                               drf=max(1.0D-10, aretper(right)-aretper(left) )
                               o%fd(lo)=((v%fd(li)-aretper(left))/drf)&
                                    & * (r%ret%parid%fd(right)-r%ret%parid%fd(left)) &
                                    & + r%ret%parid%fd(left)
                            end if
                         else                        ! minimum event
                            if (v%fd(li).lt.aretper(aretmin(1))) then ! below
                               left=newnn
                               right=1
                               left=aretmin(left)
                               right=aretmin(right)
                               o%fd(lo)=extrapolate(aretper(left),r%ret%parid%fd(left),&
                                    & aretper(right),r%ret%parid%fd(right), &
                                    & v%fd(li))
                               ! o%fd(lo)=r%ret%parid%fd(aretmin(1))
                            else if (v%fd(li).gt.aretper(aretmin(newnn))) then ! above highest return_level
                               if (i%lfldat(31)) then
                                  if (aretper(aretmin(newnn)).le.0.0D0.and.v%fd(li).le.i%zero) then
                                     left=aretmin(newnn)
                                     drf=min(-1.0D-10, aretper(left))
                                     o%fd(lo)=((v%fd(li)-i%zero)/drf) * r%ret%parid%fd(left)
                                  else
                                     o%fd(lo)=o%filld
                                  end if
                               else
                                  o%fd(lo)=o%filld
                               end if
                            else
                               call sort_heapsearch1r(nretper,aretper,eps,newnn,aretmin,v%fd(li),left,right)
                               if (right.eq.0) then
                                  write(*,*)myname," Heapsearch1r system error:'",v%var250(1:v%lenv)//"'",&
                                       & v%fd(li),newnn,left,right,eps
                                  write(*,*)myname," aretper:",aretper
                                  call ncf_printPos(v%f%pos)
                                  irc=912
                                  return
                               end if
                               left=aretmin(left)
                               right=aretmin(right)
                               drf=max(1.0D-10, aretper(right)-aretper(left) )
                               o%fd(lo)=((v%fd(li)-aretper(left))/drf)&
                                    & * (r%ret%parid%fd(right)-r%ret%parid%fd(left)) &
                                    & + r%ret%parid%fd(left)
                            end if
                         end if
                         if (o%fd(lo).ne.o%filld) then
                            if (valid_first) then
                               valid_first=.false.
                               valid_max=o%fd(lo)
                               valid_min=o%fd(lo)
                            else
                               valid_max=max(valid_max,o%fd(lo))
                               valid_min=min(valid_min,o%fd(lo))
                            end if
                         end if
                         if (ncf_bdeb) then
                            write(*,*)'Output:',li,v%fd(li),lo,o%fd(lo),maxevent,left,right
                         end if
                         ook(2)=ook(2)+1
                      else
                         orm(2)=orm(2)+1
                      end if
                   end do LRGRID
                   ook(1)=ook(1)+1
                else
                   orm(1)=orm(1)+1
                end if
             end do LLGRID
             if (boz) then ! make return-level variable
                ! make sure positions are the same in input/output...
                call ncf_importVariable(out,zi,irc)
                if (irc.ne.0) then
                   write(*,*)myname,&
                        & ' Error return from ncf_importvariable.',irc
                   return
                end if
             end if
             call ncf_importVariable(out,o,irc) ! new variable...
             if (irc.ne.0) then
                write(*,*)myname,&
                     & ' Error return from ncf_importvariable.',irc
                return
             end if
             !
             ! set max/min attributes...
             if (.not.valid_first) then
                call ncf_setAttribute(o,"valid_min",valid_min,o%type)
                call ncf_setAttribute(o,"valid_max",valid_max,o%type)
             end if
             ! clean up
             call ncf_clearDimOrder(ido)
             ! add return level variable
             if (boz) then
                nullify(zi)
                call ncf_clearDimOrder(zdo)
             end if
             if(irc.ne.0)return
             ! exit LRET
          end if
          if (.not. binp) then ! variable is in out, not inp...
             call ncf_importVariable(out,v,irc)
             if (irc.ne.0) then
                write(*,*)myname,&
                     & ' Error return from ncf_importvariable.',irc
                return
             end if
          end if
          !end do LRET
          if (.not.found) then
             write(*,*)myname,"Missing variable '"//i%ret250(ivar)(1:i%lenret(ivar))//"'"
             irc=845
             return
          end if
          if (first)WRITE(*,*) MYNAME," No data for: '"//i%net250(ivar)(1:i%lennet(ivar))//"'"
       end do LREPER
       if (allocated(aretper)) deallocate(aretper)
       if (allocated(aretind)) deallocate(aretind)
       if (allocated(aretmax)) deallocate(aretmax)
       if (allocated(aretmin)) deallocate(aretmin)
       call ncf_clearWeight(wgt)
       pst(1)=100.0D0*real(ook(1))/real(max(ook(1)+orm(1),1))
       pst(2)=100.0D0*real(ook(2))/real(max(ook(2)+orm(2),1))
       write(*,'(X,A,A,I0,"(",F0.1,"%)")') myname,'Valid interpolation: ',ook(1),pst(1)
       write(*,'(X,A,A,I0,"(",F0.1,"%)")') myname,'Accepted grid points:',ook(2),pst(2)
    end if
    if (bdeb) then
       call ncf_checkInventory(inp,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_checkInventory B.',irc
          return
       end if
    end if
    write(*,*)myname,' Closing: ',i%rp250(r%tnrrp)(1:lenp)
  end subroutine addReturnPeriod
  !
  ! extrapolate return period
  !
  real function extrapolate(val0, yrp0, val1, yrp1, val)
    real :: val0, yrp0
    real :: val1, yrp1
    real :: val, yrp
    real dval, dyrp, sigma, factor
    real,parameter :: gamma=10.0D0
    real,parameter :: eps=0.01D0
    dval=val0-val1
    dyrp=Log(yrp0)-Log(yrp1)
    if (abs(dval).lt.1.0D-10) then
       yrp=yrp1 + eps
    else
       !! Generalized Pareto Distribution
       !sigma=dyrp/dval
       !factor=(1.0D0+gamma*(val-val1)/sigma)**(1.0D0/gamma)
       !yrp=yrp1*factor
       !!write(*,*)'Extrapolate:',val0,yrp0,val1,yrp1,val,yrp,sigma,factor
       !!yrp=yrp1*Exp(dyrp*(val-val1)/dval)
       yrp=yrp1 + eps
    end if
    extrapolate=min(yrp,yrp1*10.0D0)
    return
  end function extrapolate
  !
  ! Add PST or AVG variables
  !
  subroutine addExtra(inp,out,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'addExtra'/
    integer :: ivar
    type(variable), pointer :: v=>null(),oa=>null()
    logical :: btmp
    do ivar=1,i%nrgrp
       v => ncf_getvariable(inp,i%ogr250(ivar)(1:i%lenogr(ivar)),btmp,irc)
       if (.not.btmp) then
          v => ncf_getvariable(out,i%ogr250(ivar)(1:i%lenogr(ivar)),btmp,irc)
       end if
       ! make sure variable is loaded (it could be new)
       if (ncf_hasNoData(v)) then
          call ncf_readRealData(v,r%bop,irc)
          if (.not.r%bop) irc=999
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_readRealData X.',irc
             return
          end if
       end if
       if (btmp) then
          if(bdeb)write(*,*)myname,'Generating "'//i%grp250(ivar)(1:i%lengrp(ivar))//'"',&
               & ivar,i%lengrp(ivar)
          if (i%lavg) then
             call makeAvgVariable(inp,out,v,oa,&
                  & i%grp250(ivar)(1:i%lengrp(ivar)),i%dim250,irc)
             if(irc.ne.0)return
          else if (i%lpst) then
             call makePstVariable(inp,out,v,oa,i%nrpst,i%dimpst,i%outpst,&
                  & i%grp250(ivar)(1:i%lengrp(ivar)),i%dim250,irc)
             if(irc.ne.0)return
          end if
       end if
    end do
  end subroutine addExtra
  !
  ! Add auxiliary variables
  !
  subroutine makeAuxiliary(aux,out,i,r,irc)
    implicit none
    type(inventory), pointer :: aux
    type(inventory), pointer :: out
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'makeAuxiliary'/
    integer :: ivar, ii, jj, lend
    type(dimensionorder), pointer :: ado => null(),xdo => null(),odo => null()
    type(dimension), pointer :: dx => null()
    type(variable), pointer :: v=>null(),ov=>null()
    logical :: btmp, otmp
    do ivar =1,i%nraux
       btmp=.false.
       otmp=.false.
       v => ncf_getvariable(aux,i%avar250(ivar)(1:i%lenavar(ivar)),btmp,irc)
       if (btmp) then
          ov => ncf_getvariable(out,i%avar250(ivar)(1:i%lenavar(ivar)),otmp,irc)
          if (otmp) then
             write(*,*)myname,' Auxiliary variable exists:',i%avar250(ivar)(1:i%lenavar(ivar))
             irc=112
             return
          end if
       else
          write(*,*)myname,' Auxiliary variable is missing:',i%avar250(ivar)(1:i%lenavar(ivar))
          irc=112
          return
       end if
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_getVariable.',irc
          return
       end if
       ! create output variable
       if (btmp.and..not.otmp) then
          call ncf_clearDimOrder(ado)
          ado => ncf_makeDimOrder(v,irc)
          ! remove extra dimensions
          call ncf_clearDimOrder(xdo)
          xdo => ncf_createDimOrder(irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from createDimOrder.',irc
             return
          end if
          call ncf_setDimOrderInventory(xdo,aux)
          do ii=1,i%nraud(ivar)
             lend=length(i%adim250(ivar,ii),250,10)
             !write(*,*)myname,"Checking dimension=",i%adim250(ivar,ii)(1:lend),lend
             jj=ncf_getDimEntry(aux,i%adim250(ivar,ii)(1:lend))
             if (jj.eq.0) then
                write(*,*)myname,"Missing auxiliary dim=",i%adim250(ivar,ii)(1:lend)
             else
                dx => ncf_makeDimension(aux,jj)
                ! sanity check
                if (dx%siz .ne. 1) then ! can not ignore non-unit dimensions
                   write(*,*)myname,"Auxiliary '",i%avar250(ivar)(1:i%lenavar(ivar)),&
                        & "' dim=",i%adim250(ivar,ii)(1:lend)," has size=",dx%siz
                   irc=845
                   return
                else
                   call ncf_setDimensionValue(aux,dx,1)
                   call ncf_addDimOrderDim(xdo,dx,irc)
                   if (irc.ne.0) then
                      write(*,*)myname,' Error return from ncf_addDimOrdDim.',irc
                      return
                   end if
                end if
             end if
          end do
          call ncf_removeDimOrder(ado,xdo) ! remove ignored dimensions from auxiliary
          call ncf_clearDimOrder(odo)
          odo => ncf_mapDimOrder(out,ado,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncd_mapDimOrder.',irc
             return
          end if
          ov => ncf_copyVariable(v,irc) ! the new ouput variable
          call ncf_setVariableDimOrder(ov,odo,irc)
          call ncf_appendVariable(out,ov)
       end if
    end do
    return
  end subroutine makeAuxiliary
  !
  ! Add auxiliary variables
  !
  subroutine addAuxiliary(aux,out,i,r,irc)
    implicit none
    type(inventory), pointer :: aux
    type(inventory), pointer :: out
    type(init) :: i
    type(run) :: r
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'addAuxiliary'/
    integer :: ivar
    type(variable), pointer :: av=>null(),ov=>null()
    logical :: btmp
    do ivar=1,i%nraux
       av => ncf_getvariable(aux,i%avar250(ivar)(1:i%lenavar(ivar)),btmp,irc)
       if (btmp) then
          ov => ncf_getvariable(out,i%avar250(ivar)(1:i%lenavar(ivar)),btmp,irc)
       end if
       ! make sure variable is loaded (it could be new)
       if (btmp.and.ncf_hasNoData(av)) then
          call ncf_readRealData(av,btmp,irc)
          if (.not.btmp) irc=999
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_readRealData A.',irc
             return
          end if
       end if
       if (btmp) then ! copy to output variable...
          call ncf_copyVariableField(ov,av,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_copyVariableField.',irc
             return
          end if
       else
          write(*,*)myname,'Unable to add:',i%avar250(ivar)(1:i%lenavar(ivar))
       end if
    end do
  end subroutine addAuxiliary
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! IGNORE UNWANTED OR LOAD WANTED VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine readyVariables(inp,out,aux,i,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(inventory), pointer :: aux
    type(init) :: i
    type(run) :: r
    integer :: irc
    type(variable), pointer  :: v=>null(),n=>null(),iv=>null()
    integer :: ifound,ii
    logical :: found
    CHARACTER*16 MYNAME
    DATA MYNAME /'readyVariables'/
    !
    found=.false.
    v => out%firstVariable%next
    DEL: do while (.not.associated(v,target=out%lastVariable))
       n => v%next
       ifound = 0
       do ii=1,i%nrvar
          if (i%var250(ii)(1:i%lenv(ii)).eq.v%var250(1:v%lenv)) then
             ifound=ii
             found=.true.
          end if
       end do
       if (ifound.eq.0) then
          if(bdeb)write(*,*) myname,'Deleting: ',v%lend,v%var250(1:v%lenv)
          call ncf_ignoreVariable(v)
          call ncf_clearVariable(v,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_clearVariable.',irc
             return
          end if
       else
          if (.not. v%readyfield) then
             if(bdeb)write(*,*) myname,'Loading:  ',v%lend,v%var250(1:v%lenv)
             iv => ncf_getVariable(inp,v%var250(1:v%lenv),r%bok,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_getVariable X.',irc,iv%lend
                return
             end if
             if (associated(iv)) then
                if(bdeb)write(*,*) myname,'Reading:  ',iv%lend,iv%var250(1:iv%lenv)
                r%bok=.true.
                call ncf_readData(iv,r%bok,irc)
                if (.not.r%bok) irc=999
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_readRealData N.',irc,v%lend
                   return
                end if

                ! this wont work for some reason.... ??????????????????

                !write(*,*)myname,'Copyfield:',iv%initialised,iv%readyfield,&
                 !    & v%initialised,v%readyfield,size(iv%fd),associated(iv%fd),&
                     ! & size(v%fd),associated(v%fd),r%bok
                call ncf_copyVariableField(v,iv,irc)
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_copyVariableField',irc,v%lend
                   return
                end if
             else
                write(*,*) myname,'Missing:  ',v%lend,v%var250(1:v%lenv)
             end if
          else
             if(bdeb)write(*,*) myname,'Keeping:  ',v%lend,v%var250(1:v%lenv),v%readyField
          end if
          if (i%lenn(ifound).ne.0) then  ! rename variable
             if(bdeb)write(*,*) myname,'Renaming '//v%var250(1:v%lenv)//' => '//i%new250(ifound)(1:i%lenn(ifound))
             call ncf_setOutVariableName(v,i%new250(ifound)(1:i%lenn(ifound)))
          endif
       end if
       v=>n
       ! exit DEL
    end do DEL
    if (.not.found) then
       write(*,*)myname,'Keep cnt:',i%nrvar
       do ii=1,i%nrvar
          write(*,*)myname,'No variable match:',ii,i%var250(ii)(1:i%lenv(ii))
       end do
    end if
  end subroutine readyVariables
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DUMP SLICE TO OUTPUT FILE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine writeSlice(inp,out,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    integer :: irc
    CHARACTER*14 MYNAME
    DATA MYNAME /'writeSlice'/
    type(variable), pointer :: v=>null()
    !call  ncf_printInventory(out)
    if(bdeb)write(*,*) myname,' Compressing variables.'
    call ncf_compressVariables(out,irc) ! compress variables
    if (irc.ne.0) then
       write(*,*)myname,' Error return from ncf_compressVariables.',irc
       return
    end if
    !if(bdeb)write(*,*) myname,' Compressing dimensions.'
    !call ncf_compressDimensions(out) ! compress dimensions
    if(bdeb)write(*,*) myname,' Writing output.'
    if (r%initWrite) then
       r%initWrite=.false.
       out%fn250=i%out250
       out%lenf=length(i%out250,250,10)
       write(*,*)myname,' Writing: ',out%fn250(1:out%lenf)
       !call ncf_checkInventory(out,irc)
       call ncf_writeNcOpen(out,i%out250,irc) ! open, write, close file
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_writeNcOut.',irc
          return
       end if
       ! do not re-write variables that do not have slice dimension...
       v => out%firstVariable%next
       REWRITE: do while (.not.associated(v,target=out%lastVariable))
          v%writefield=(ncf_hasDimension(v,r%dslice))
          v => v%next
       end do REWRITE
    else
       call ncf_writeNcData(out,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_writeNcData.',irc
          return
       end if
    end if
    if(bdeb)write(*,*)myname,'Clear variables...',irc
    !call ncf_clearVariables(out,irc)
    !if (irc.ne.0) then
    !   write(*,*)myname,' Error return from ncf_clearVariables.',irc
    !   return
    !end if
    if (bdeb) then
       call ncf_checkInventory(out,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from checkInventory.',irc
          return
       end if
    end if
  end subroutine writeSlice
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WRAP UP AND WRITE OUTPUT FILE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine wrapUp(inp,out,aux,i,r,irc)
    implicit none
    type(inventory), pointer :: inp
    type(inventory), pointer :: out
    type(inventory), pointer :: aux
    type(init) :: i
    type(run) :: r
    integer :: irc
    if (r%brp) then
       if(bdeb)write(*,*)myname,'Cleaning up...',irc
       ! close return-period file
       call ncf_closeFile(r%ret,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_closeFile.',irc
          return
       end if
       !
       call ncf_clearInventory(r%ret,irc) ! delete internally stored data...
       if (irc.ne.0) then
          write(*,*) myname,"Error return from ncf_clearInventory.",irc
          return
       end if
       if (associated(r%ret)) deallocate(r%ret)
    end if
    !
    if (.not. r%initWrite) then
       call ncf_closeNcFile(out,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_closeNcFile.',irc
          return
       end if
    end if
    !
    ! wrap up
    if (r%bop) then
       write(*,*)myname,' Closing: ',i%inp250(1:i%leni),irc
       call ncf_closeFile(inp,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_closeFile.',irc
          return
       end if
    end if
    !
    ! wrap up
    if (r%aop) then
       write(*,*)myname,' Closing: ',i%aux250(1:i%lena),irc
       call ncf_closeFile(aux,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_closeFile.',irc
          return
       end if
    end if
    !
    if(bdeb)write(*,*) myname,' Cleaning memory.',irc
    call ncf_clearInventory(out,irc) ! delete internally stored data...
    if (irc.ne.0) then
       write(*,*) myname,"Error return from ncf_clearInventory.",irc
       return
    end if
    if(bdeb)write(*,*) myname,' Clearing memory.'
    call ncf_clearInventory(inp,irc) ! delete internally stored data...
    if (irc.ne.0) then
       write(*,*) myname,"Error return from ncf_clearInventory.",irc
       return
    end if
    if(bdeb)write(*,*) myname,' Clearing memory.'
    call ncf_clearInventory(aux,irc) ! delete internally stored data...
    if (irc.ne.0) then
       write(*,*) myname,"Error return from ncf_clearInventory.",irc
       return
    end if
  end subroutine wrapUp
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FILTER VARIABLE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine filterInput(i,r,v,irc)
    implicit none
    type(init) :: i
    type(run) :: r
    type(variable), pointer :: v
    integer :: irc
    i%iiflt=filter(v,i%nriflt,i%iflt250,i%leniflt,irc)
    if (irc.ne.0) then
       write(*,*)myname,' Error return from filter.',irc
       return
    end if
    if (i%iiflt.ne.0) then
       i%biflt=.true.
       i%minf=i%miniflt(i%iiflt)
       i%maxf=i%maxiflt(i%iiflt)
    else
       i%biflt=.false.
    end if
    !
    ! if (i%biflt.and.v%fd(li).ne.v%filld) then
    !    if (v%fd(li).lt.i%minf.or.v%fd(li).gt.i%maxf) then
    !       v%fd(li)=v%filld
    !    end if
    ! end if
  end subroutine filterInput
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! OTHER SUBROUTINES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  integer function differentiate(v,nrdiff,diff250,lendiff,irc)
    implicit none
    type(variable), pointer :: v
    integer :: nrdiff
    character*250 :: diff250(nrdiff)
    integer :: lendiff(nrdiff)
    integer :: irc
    integer :: ii
    differentiate=0
    do ii=1,nrdiff
       if (diff250(ii)(1:lendiff(ii)).eq.v%var250(1:v%lenv)) then
          differentiate=ii
          return
       end if
    end do
    return
  end function differentiate
  !
  integer function filter(v,nrflt,flt250,leniflt,irc)
    implicit none
    type(variable), pointer :: v
    integer :: nrflt
    character*250 :: flt250(nrflt)
    integer :: leniflt(nrflt)
    integer :: irc
    integer :: ii
    filter=0
    do ii=1,nrflt
       if (flt250(ii)(1:leniflt(ii)).eq.v%var250(1:v%lenv)) then
          filter=ii
          return
       end if
    end do
    return
  end function filter
  !
  
  subroutine makeAvgVariable(inp,out,v,o,net,dim250,irc)
    implicit none
    type(inventory), pointer  :: inp
    type(inventory), pointer  :: out
    type(variable), pointer  :: v
    type(variable), pointer  :: o ! output variable
    character*(*) :: net
    character*250 :: dim250
    integer :: irc
    character*14 :: myname="makeAvgVariable"
    integer, external :: length
    integer :: lend
    type(dimensionOrder), pointer  ::  odo=>null(),ido=>null(),ado=>null()
    type(inventory), pointer :: ginv=>null()
    character*250 :: net250
    real,allocatable :: avgval(:)
    integer,allocatable :: avgind(:)
    type(dimension), pointer  :: avgdim=>null()
    integer :: avgentry
    logical :: btmp
    integer :: li,lo, tcnt, cnt
    !
    ! REMEMBER THAT WE IMPORT THE OUTPUT FIELD TO THE INPUT FILE,
    ! AND LOOP OVER THE INPUT DIMENSIONS...
    !
    ginv => v%f
    avgEntry=ncf_getDimEntry(ginv,dim250)
    !write(*,*)myname,'Entering...',avgentry
    !
    lend=length(dim250,250,10)

    !write(*,*)myname," Processing dim '"//dim250(1:lend)//"'",associated(v)
    if (.not.associated(v)) then
       write(*,*)myname,'Variable not associated.'
       irc=812
       return
    end if
    if(bdeb)write(*,*) myname,' Processing: ',v%var250(1:v%lenv),v%lend
    !
    ido=>ncf_makeDimOrder(v,irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
       return
    end if
    !
    odo => ncf_copyDimOrder(ido,irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
       return
    end if
    !
    avgdim => ncf_removeDimOrderEntry(odo,avgEntry)

    if (associated(avgdim)) then ! average necessary

       ! make averaging dimension order
       ado => ncf_newDimOrder(out,irc)
       if(irc.ne.0)return
       !
       call ncf_addDimOrderDim(ado,avgdim,irc)
       if(irc.ne.0)then
          write(*,*)myname,'Error return from addDimOrderDim.',irc
          return
       end if
       !
       tcnt=0
       ! insert values into field
       if(bdeb)write(*,*)myname,'Looking for "'//net//'"'
       o => ncf_getvariable(out,net,btmp,irc)
       if (.not. btmp) then ! create output variable
          o => ncf_copyVariable(v,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyVariable.',irc
             return
          end if
          net250=net
          call ncf_setVariableName(o,net250)
          call ncf_setTextAttribute(o,"standard_name",net)
          !
          call ncf_setVariableDimOrder(o,odo,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_setVariableDimOrder.',irc
             return
          end if
          call ncf_addAttribute(o,"Averaged_over_dim",dim250(1:lend))
          if(bdeb)write(*,*)myname,"### Prepending H:'"//&
               & o%var250(1:o%lenv)//"'",o%lend
          call ncf_prependVariable(out,o)
       end if
       ! make sure positions are the same...
       call ncf_importVariable(inp,o,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       !
       call ncf_initField(o,o%filld,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_initField.',irc
          return
       end if
       if(bdeb)write(*,*)myname,' Averaging: ',v%var250(1:v%lenv),v%lend
       ! make sure positions are the same...
       call ncf_importDimOrder(inp,odo,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       !
       ! make sure positions are the same...
       call ncf_importDimOrder(inp,ado,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       !
       call ncf_importVariable(inp,v,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       !call ncf_printdimorder(odo)
       if(bdeb)write(*,*)myname,'Starting outer loop...'
       !call ncf_printdimorder(ado)
       call ncf_resetPos(odo,irc)
       if(irc.ne.0)return
       OTHR:do while (ncf_increment(inp,odo,irc))
          ! loop over input data, updating output arrays
          call ncf_resetPos(ado,irc)
          if(irc.ne.0)return
          lo = ncf_getLocation(o)
          o%fd(lo)=0.0D0
          cnt=0
          AVG:do while (ncf_increment(inp,ado,irc))
             li = ncf_getLocation(v)
             if (v%fd(li).ne.v%filld) then
                !if(bdeb)write(*,*)myname,' Location:',lo,o%lend,associated(o%fd)
                o%fd(lo)=o%fd(lo) + v%fd(li)
                cnt=cnt+1
                !call ncf_printPos(o%f%pos)
                !write(*,*)myname,' Location:',v%var250(1:v%lenv),loc,lic
                ! if (cnt .lt. 50) then
                ! end if
             end if
             tcnt=tcnt+1
             ! if (mod(tcnt,100000).eq.0) then
             !    write(*,*)myname,' Location:',v%var250(1:v%lenv),lo,li,v%fd(li),o%fd(lo),cnt
             ! end if
          end do AVG
          if (cnt.ne.0) then
             o%fd(lo)=o%fd(lo)/max(1.0D0,real(cnt))
          else
             o%fd(lo)=o%filld
          end if
          if(irc.ne.0)return
       end do OTHR
       call ncf_importVariable(ginv,v,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       call ncf_importVariable(out,o,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       if (allocated(avgval)) deallocate(avgval)
       if (allocated(avgind)) deallocate(avgind)
       !call ncf_checkInventory(out,irc)
    else
       if(bdeb)write(*,*)myname,' Not processed: ',v%var250(1:v%lenv),v%lend
    end if
    !call ncf_checkInventory(out,irc)
    !write(*,*)myname,'Done...',avgentry
    call ncf_clearDimOrder(ido)
    call ncf_clearDimOrder(odo)
  end subroutine makeAvgVariable
  !
  subroutine makePstVariable(inp,out,v,o,nrpst,pst,name,net,dim250,irc)
    implicit none
    type(inventory), pointer  :: inp
    type(inventory), pointer  :: out
    type(variable), pointer  :: v
    type(variable), pointer  :: o ! output variable
    character*(*) :: net
    integer :: nrpst
    real :: pst(nrpst), name(nrpst)
    character*250 :: dim250
    integer :: irc
    character*14 :: myname="makePstVariable"
    integer, external :: length
    integer :: lend
    type(dimensionOrder), pointer  ::  odo=>null(),ido=>null(),pstdo=>null(),ldo=>null()
    integer :: pstlen, opstlen=0
    type(inventory), pointer :: ginv=>null()
    real :: rind
    integer :: nind
    character*250 :: net250
    type(variable), pointer  :: pstvar=>null() ! output variable
    real,allocatable :: pstval(:)
    integer,allocatable :: pstind(:)
    integer :: nn=0,newnn
    type(dimension), pointer  :: ensdim=>null(), pstdim=>null()
    integer :: ensentry,pstentry
    character*(*), parameter :: pstname="percentile"
    logical :: btmp
    integer :: li,lo, tcnt,ii
    !
    ! REMEMBER THAT WE IMPORT THE OUTPUT FIELD TO THE INPUT FILE,
    ! AND LOOP OVER THE INPUT DIMENSIONS...
    !
    ginv => v%f
    ensentry=ncf_getDimEntry(ginv,dim250)
    !write(*,*)myname,'Entering...',ensentry
    !
    lend=length(dim250,250,10)

    !write(*,*)myname," Processing dim '"//dim250(1:lend)//"'",associated(v)
    if (.not.associated(v)) then
       write(*,*)myname,'Variable not associated.'
       irc=812
       return
    end if
    if(bdeb)write(*,*) myname,' Processing: ',v%var250(1:v%lenv),v%lend
    !
    odo => ncf_makeDimOrder(v,irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
       return
    end if
    !
    ! look for percentile dimension
    ensdim => ncf_removeDimOrderEntry(odo,ensentry)
    if (nrpst.ge.1) then
       pstentry=ncf_getDimEntry(out,pstname)
       if (pstentry.eq.0) then ! create percentile dimension...
          ! create percentile dimension in input file...
          ! ...so we can import output variables
          pstdim => ncf_createDimension(inp,pstname,nrpst)
          ! create percentile dimension in output file..
          pstdim => ncf_createDimension(out,pstname,nrpst)
          pstentry=pstdim%ind
          ! create variable also...
          pstvar => ncf_getvariable(out,pstname,btmp,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_getVariable.',irc
             return
          end if
          if (.not.btmp) then ! create percentile-variable
             pstdo => ncf_newDimOrder(out,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_newDimOrder.',irc
                return
             end if
             call ncf_addDimOrderEntry(pstdo,pstentry,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_addDimOrder.',irc
                return
             end if
             pstvar => ncf_createVariable(out,pstname,nf_int,pstdo,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_createVariable.',irc
                return
             end if
             call ncf_setVariableName(pstvar,pstname)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_setVariableName.',irc
                return
             end if
             call ncf_setTextAttribute(pstvar,"standard_name",pstname)
             call ncf_initField(pstvar,pstvar%filld,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_initField.',irc
                return
             end if
             do ii=1,pstvar%lend
                pstvar%fd(ii)=name(ii)
             end do
             if(bdeb)write(*,*)myname,"### Prepending I:'"//&
                  & pstvar%var250(1:pstvar%lenv)//"'",pstvar%lend
             call ncf_prependVariable(out,pstvar)
          end if
       else
          pstdim => ncf_makeDimension(out,pstentry)
       end if
    else
       pstentry=0
    end if
    ! make variable...
    if (associated(ensdim)) then ! average necessary
       ldo => ncf_copyDimOrder(odo,irc) ! outer loop
       if (pstentry.ne.0) then ! we have a percentile dimension
          ido=>ncf_newDimOrder(inp,irc) ! inner loop
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
             return
          end if
          ! add percentile dimension to odo
          call ncf_addDimOrderEntry(odo,pstentry,irc)
          call ncf_addDimOrderEntry(ido,pstentry,irc)
       end if
       ! make averaging dimension order
       pstdo => ncf_newDimOrder(out,irc)
       if(irc.ne.0)return
       !
       call ncf_addDimOrderDim(pstdo,ensdim,irc)
       if(irc.ne.0)then
          write(*,*)myname,'Error return from addDimOrderDim.',irc
          return
       end if
       !
       tcnt=0
       ! insert values into field
       if(bdeb)write(*,*)myname,'Looking for "'//net//'"'
       o => ncf_getvariable(out,net,btmp,irc)
       if (.not. btmp) then ! create output variable
          o => ncf_copyVariable(v,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyVariable.',irc
             return
          end if
          net250=net
          call ncf_setVariableName(o,net250)
          call ncf_setTextAttribute(o,"standard_name",net)
          !
          call ncf_setVariableDimOrder(o,odo,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_setVariableDimOrder.',irc
             return
          end if
          lend=length(dim250,250,10)
          call ncf_addAttribute(o,"Percentile_dim",dim250(1:lend))
          if(bdeb)write(*,*)myname,"### Prepending J:'"//&
               & o%var250(1:o%lenv)//"'",o%lend
          call ncf_prependVariable(out,o)
       end if
       ! make sure positions are the same...
       call ncf_importVariable(inp,o,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       !
       call ncf_initField(o,o%filld,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_initField.',irc
          return
       end if
       if(bdeb)write(*,*)myname,' Percentile: ',v%var250(1:v%lenv),v%lend
       pstlen=ncf_getDimOrderLength(pstdo)
       if (pstlen.ne.opstlen) then
          if (allocated(pstval)) deallocate(pstval)
          if (allocated(pstind)) deallocate(pstind)
          opstlen=pstlen
          allocate(pstval(pstlen),pstind(pstlen),stat=irc)
          if(irc.ne.0) then
             write(*,*)myname,'Unable to allocate pstind, pstval.',pstlen
             return
          end if
       end if
       ! make sure positions are the same...
       call ncf_importDimOrder(inp,odo,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       !
       ! make sure positions are the same...
       call ncf_importDimOrder(inp,pstdo,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       !
       call ncf_importVariable(inp,v,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       !call ncf_printVariable(v)
       !call ncf_printdimorder(ldo)
       if(bdeb)write(*,*)myname,'Starting outer loop...'
       !call ncf_printdimorder(pstdo)
       call ncf_resetPos(ldo,irc)
       if(irc.ne.0)return
       OTHR:do while (ncf_increment(inp,ldo,irc))
          ! loop over input data, updating output arrays
          call ncf_resetPos(pstdo,irc)
          if(irc.ne.0)return
          nn=0
          tcnt=tcnt+1
          PER:do while (ncf_increment(inp,pstdo,irc))
             li = ncf_getLocation(v)
             if (i%biflt.and.v%fd(li).ne.v%filld) then
                if (v%fd(li).lt.i%minf.or.v%fd(li).gt.i%maxf) then
                   v%fd(li)=v%filld
                end if
             end if
             if (v%fd(li).ne.v%filld) then
                nn=nn+1
                pstind(nn)=nn
                pstval(nn)=v%fd(li)
             end if
             ! if (mod(tcnt,100000).eq.0) then
             !    write(*,*)myname,' Location:',v%var250(1:v%lenv),&
             !         & li,v%fd(li),(v%fd(li).ne.v%filld),nn
             ! end if
          end do PER
          if (nn.ne.0) then
             call sort_heapsort1r(pstlen,pstval,eps,newnn,nn,pstind,uniq)
             if (pstentry.ne.0) then
                call ncf_resetPos(ido,irc)
                if(irc.ne.0)return
                INNER:do while (ncf_increment(inp,ido,irc))
                   lo = ncf_getLocation(o)
                   ! loop over pstdim
                   ii=ncf_getDimensionValue(inp,pstdim)
                   rind=((real(newnn-1)*i%dimpst(ii))/100.0D0)+1.0D0
                   nind=int(rind)
                   rind=rind-real(nind)
                   if (rind.gt.0.0D0.and.nind.lt.newnn) then
                      o%fd(lo)=pstval(pstind(nind))*(1.0D0-rind)+ &
                           & pstval(pstind(nind+1))*rind
                   else
                      o%fd(lo)=pstval(pstind(nind))
                   end if
                   !!!write(*,*)myname,'PST:',ii,lo,rind,o%fd(lo),tcnt
                end do INNER
                !!!if (tcnt>2) stop("Debug");
             else
                lo = ncf_getLocation(o)
                rind=((real(newnn-1)*i%dimpst(1))/100.0D0)+1.0D0
                nind=int(rind)
                rind=rind-real(nind)
                if (rind.gt.0.0D0.and.nind.lt.newnn) then
                   o%fd(lo)=pstval(pstind(nind))*(1.0D0-rind)+ &
                        & pstval(pstind(nind+1))*rind
                else
                   o%fd(lo)=pstval(pstind(nind))
                end if
             end if
          end if
          if(irc.ne.0)return
       end do OTHR
       call ncf_importVariable(ginv,v,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       call ncf_importVariable(out,o,irc)
       if (irc.ne.0) then
          write(*,*)myname,&
               & ' Error return from ncf_importvariable.',irc
          return
       end if
       if (allocated(pstval)) deallocate(pstval)
       if (allocated(pstind)) deallocate(pstind)
       opstlen=0
       !call ncf_printVariable(o)
       !call ncf_checkInventory(out,irc)
    else
       if(bdeb)write(*,*)myname,' Not processed: ',v%var250(1:v%lenv),v%lend
    end if
    !call ncf_checkInventory(out,irc)
    !write(*,*)myname,'Done...',ensentry
    call ncf_clearDimOrder(ido)
    call ncf_clearDimOrder(ldo)
    call ncf_clearDimOrder(odo)
  end subroutine makePstVariable
  real function j2000(sec)
    implicit none
    real sec
    j2000=sec/86400.0D0 - 10957.0D0  ! convert from seconds to days
    return
  end function j2000
end subroutine mncreper

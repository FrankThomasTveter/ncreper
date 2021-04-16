! ecmwf netcdf filer...   ls /fou/nwparc2/ec/2013/.nc
SUBROUTINE MNCVERIFY(UNITI,IRC)
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
  IMPLICIT NONE
  SAVE
  ! 
  ! INTERFACE VARIABLES
  ! 
  INTEGER  UNITI
  INTEGER IRC
  ! 
  integer nrpred, maxpred
  REAL, allocatable ::   PRED(:)
  INTEGER  PCAT
  ! 
  CHARACTER*12 MYNAME
  DATA MYNAME /'MNCVERIFY'/
  ! 
  INTEGER  IKODE,OKODE
  LOGICAL  OK,LFLDAT(100),BDEB,ACTIVE,found
  DATA BDEB   /.false./
  DATA ACTIVE /.FALSE./
  ! 
  INTEGER  KODE, length, lena, lend, lenb, lenh, lene, lenr
  external length
  ! 
  INTEGER  LINE
  INTEGER  NRHDR
  PARAMETER (NRHDR=100)
  CHARACTER*250 HDR250(NRHDR),DAT250, INF250,buff250,buffx250
  CHARACTER*250 FILENM,NUKEHEAD
  EXTERNAL NUKEHEAD
  INTEGER  INTOUT
  LOGICAL  ENDOFF
  ! 
  type inventory
     logical       :: initialised = .false.
     integer       :: ncid
     character*250 :: fn250
     integer       :: lenf
     logical       :: opened = .false.
     integer       :: nrdim
     integer       :: nrvar
     integer       :: nrgatt
     integer       :: unlimdimid
     character*250, pointer    :: dim250(:) => null()
     integer, pointer          :: lend(:) => null()
     integer, pointer          :: dimid(:) => null()
     type(position), pointer   :: pos
     type(variable), pointer   :: firstVariable
     type(variable), pointer   :: lastVariable
     type(variable), pointer   :: latid => null()
     type(variable), pointer   :: lonid => null()
     type(variable), pointer   :: timid => null()
     type(variable), pointer   :: terid => null()
     type(variable), pointer   :: parid => null()
     type(variable), pointer   :: xcoid => null()
     type(variable), pointer   :: ycoid => null()
     type(variable), pointer   :: gridid => null()
     type(attribute), pointer  :: firstAttribute
     type(attribute), pointer  :: lastAttribute
  end type inventory

  type variable
     logical       :: initialised = .false.
     type(inventory), pointer  :: f => null()
     type(variable), pointer  :: prev => null()
     type(variable), pointer  :: next => null()
     integer       :: ncid
     integer       :: varid
     character*250 :: var250
     integer       :: lenv
     integer       :: type
     integer       :: nrdim
     integer       :: nratt
     real          :: scale
     integer,    pointer :: ind(:) => null()
     integer,    pointer :: sta(:) => null()
     integer,    pointer :: lim(:) => null()
     character*1,pointer::  fc(:) => null()
     integer*1,  pointer::  f1(:) => null()
     integer*2,  pointer::  f2(:) => null()
     integer*4,  pointer::  f4(:) => null()
     real*4,     pointer::  fr(:) => null()
     real*8,     pointer::  fd(:) => null()
     integer :: len=0
     integer :: lenc=0
     integer :: len1=0
     integer :: len2=0
     integer :: len4=0
     integer :: lenr=0
     integer :: lend=0
     type(attribute),pointer  :: firstAttribute
     type(attribute),pointer  :: lastAttribute
     type(attribute),pointer  :: fillAttribute => null()
  end type variable

  type attribute
     logical       :: initialised = .false.
     type(variable), pointer  :: v => null()
     type(attribute), pointer  :: prev => null()
     type(attribute), pointer  :: next => null()
     integer       :: ncid
     integer       :: varid
     integer       :: attid
     character*250  :: att250
     integer :: lena
     integer :: type
     integer :: len
     character*1, pointer ::  ac(:) => null()
     integer*1,  pointer::  a1(:) => null()
     integer*2,  pointer::  a2(:) => null()
     integer*4,  pointer::  a4(:) => null()
     real*4,     pointer::  ar(:) => null()
     real*8,     pointer::  ad(:) => null()
  end type attribute

  type position
     integer :: maxdim=0
     integer :: nrdim
     integer, pointer :: pos(:)
     integer, pointer :: sta(:)
     integer, pointer :: lim(:)
     integer        :: loc =0
  end type position

  type dimension
     logical          :: initialised = .false.
     type(dimension), pointer  :: prev => null()
     type(dimension), pointer  :: next => null()
     integer          :: ind
     integer          :: lim
     integer          :: sta
     logical          :: increase=.true.
  end type dimension

  type index
     integer nrdim
     type(dimension), pointer  :: firstDimension
     type(dimension), pointer  :: lastDimension
     logical :: resetIndex=.true.
  end type index
  !
  type weight
     integer :: nweight = 0
     type(position), pointer :: pos(:)
     real, pointer           :: w(:)
  end type weight
  ! 
  type(inventory),pointer  :: exp, e(:) => null()
  type(inventory),pointer  :: ref, r(:) => null()
  type(inventory),pointer  :: g => null()
  type(inventory) :: out
  type(variable), pointer  :: v
  type(index),pointer  ::  expLat,expLon,expTim,expPar,expAll,expLatlon,refFoot,expFoot,expUnion
  type(index),pointer  ::  refLat,refLon,refTim,refPar,refAll,refLatlon,refOut,refUnion
  type(index),pointer  ::  outLatlon,outLat,outLon,outExpGroup,outOffset,outLatlonFinal
  type(dimension), pointer :: refI, refJ
  type(dimension), pointer :: expI, expJ
  type(weight), pointer :: refwgt, expwgt
  type(variable), pointer :: sr1,se1,sr2,se2,ser,sae,sum
  type(variable), pointer :: stdv,bias,rms,sesr,stdve,meane,sese,stdvr,meanr,srsr
  type(variable), pointer :: gStart, gEnd, oStart, oEnd
  type(Position),pointer :: rpos,epos
  ! 
  character*250 var250,epar250,rpar250,outnc250,outr250,grid250,gmap,dum250
  integer nexp, nref, ee, rr, tt, ii,jj, leng
  real latc,lonc
  integer maxfl, group
  parameter (maxfl=10000)
  character*250 exp250(maxfl), ref250(maxfl),dat250x
  integer nf,no, nn, maxnf
  parameter (maxnf=100)
  real diam, dtmin(0:maxnf), dtmax(0:maxnf),offmin(maxnf),offmax(maxnf)
  real e2000, ea2000, r2000, ra2000, t2000, dt, b2000, x2000
  real estart,estop,rstart,rstop
  integer yy,mm,dd,hh,mi
  real sec
  character*24 estart24,estop24,rstart24,rstop24, b24, x24
  ! 
#     include "netcdf.inc"
  ! 
  INTEGER ret,CHUNKSIZEHINT
  logical proceed, proceedb, bok, bbok, first, igrid, ifield
  integer cntr,cnte,cntm, maxcnt
  parameter (maxcnt=maxfl*10)
  integer fe(maxcnt)
  integer te(maxcnt)
  real je(maxcnt)
  real ae(maxcnt)
  integer ge(maxcnt)
  integer, allocatable :: ind(:)
  integer ncnt,left,right, ntime
  real evalue,rvalue,xvalue,rcnt,ecnt
  integer cnt
  integer fexp
  real eps
  logical uniq
  ! 
  character*1::  fillc
  integer*1  ::  fill1
  integer*2  ::  fill2
  integer*4  ::  fill4
  real*4     ::  fillr
  real*8     ::  filld
  real val
!
  real rscale, escale,fact
!

  integer eok(10), erm(10),rok(10), rrm(10)
  real pst(10)
  ! 
  IRC=0
  ! 
  ! DEBUG SYSTEM
  ! 
  ! IF (.NOT.ACTIVE) CALL DEBUG(MYNAME,BDEB,ACTIVE)
  ! bdeb=.true.
  ! bdeb=.FALSE.
  ! 
  nexp=0
  nref=0
  rscale=1.0D0
  escale=1.0D0
  ! 
  fillc=char(nf_fill_char)
  fill1=nf_fill_int1
  fill2=nf_fill_int2
  fill4=nf_fill_int
  fillr=nf_fill_real
  filld=nf_fill_double
  ! 
  BDEB = .false.
  IF (BDEB) WRITE(*,*) MYNAME,'Debug: Routine starts.',IRC
  nf=0
  no=0
  dtmin(:)=0.0D0
  dtmax(:)=0.0D0
  offmin(:)=0.0D0
  offmax(:)=0.0D0
  eok(:)=0
  erm(:)=0
  rok(:)=0
  rrm(:)=0
  ! 
  DO II=1,NRHDR
     HDR250(II) = ''
     LFLDAT(II)=.FALSE.
  ENDDO

  HDR250(1)   = 'NCVERIFY V1.0 [0]VFLR'
  HDR250(10)  = 'EXPERIMENT FILES (NETCDF)[*]VFLR &'
  HDR250(15)  = 'REFERENCE FILES (NETCDF) [*]VFLR &'
  HDR250(20)  = 'EXPERIMENT PARAMETER     [1]VFLR'
  HDR250(25)  = 'REFERENCE PARAMETER      [1]VFLR'
  HDR250(30)  = 'FOOTPRINT DIAMETER (KM)  [1]VFLR'
  HDR250(34)  = 'MIN, MAX TIME OFFSET, EXP-REF (HOURS) [*]VFLR'
  HDR250(35)  = 'MIN, MAX REFERENCE FORECAST (HOURS) [1]VFLR'
  HDR250(36)  = 'MIN, MAX EXPERIMENT FORECAST (HOURS) [*]VFLR'
  HDR250(40)  = 'OUTPUT FILE (NETCDF)     [1]VFLR'
  HDR250(50)  = 'OUTPUT FILE (S-PLUS)     [1]VFLR'
  HDR250(60)  = 'USE GRID FROM FILE (NETCDF) [1]VFLR'
  ! 
  ! NF_FILL_CHAR, NF_FILL_BYTE, NF_FILL_SHORT, NF_FILL_INT, NF_FILL_FLOAT, and NF_FILL_DOUBLE
  ! 
  WRITE(*,*) MYNAME,'----------------------------------------'
  ! 
  ! READ DATA FROM INPUT FILE..............................
  ! 
  KODE=-1
  CALL NUKEM(KODE,UNITI,HDR250,INTOUT,DAT250,ENDOFF,IRC)
  IF (IRC.NE.0) THEN
     WRITE(*,*) MYNAME,'-1 Error return from NUKEM.',IRC
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
     !     WRITE(*,*) MYNAME,'line=',line
     ! 
     IF (BDEB) WRITE(*,*) MYNAME,'Debug: Read header:',LINE
     ! 
     ! CHECK WHAT LINE WE JUST READ
     ! 
     IF (LINE.EQ.1) THEN
     ELSEIF (LINE.EQ.10) THEN ! experiment files
        nexp = nexp+1
        if (nexp .gt. maxfl) then
           write(*,*)myname,'Too many experiment files.',nexp
           irc=347
           return
        end if
        exp250(nexp)=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.15) THEN ! reference files
        nref = nref+1
        if (nref .gt. maxfl) then
           write(*,*)myname,'Too many reference files.',nref
           irc=347
           return
        end if
        ref250(nref)=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.20) THEN ! experiment parameter
        dat250x=dat250
        buff250=nukehead(dat250,250)
        lend=length(dat250,250,lend)
        lenb=length(buff250,250,lend)
        read(buff250(1:lenb),*,iostat=irc)escale
        if (irc.ne.0.or.lend.eq.0) then
           irc=0
           epar250=dat250x
           escale=1.0D0
        else
           epar250=dat250
        end if
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.25) THEN ! reference parameter
        dat250x=dat250
        buff250=nukehead(dat250,250)
        lend=length(dat250,250,lend)
        lenb=length(buff250,250,lend)
        read(buff250(1:lenb),*,iostat=irc)rscale
        if (irc.ne.0.or.lend.eq.0) then
           irc=0
           rpar250=dat250x
           rscale=1.0D0
        else
           rpar250=dat250
        end if
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.30) THEN ! footprint diameter (km)
        dat250x=dat250
        buff250=nukehead(dat250,250)
        lenb=length(buff250,250,lend)
        read(buff250(1:lenb),*,iostat=irc)diam
        if (irc.ne.0) then
           write(*,*) myname,'Unable to read diameter:',dat250x(1:lend)
           return
        end if
        diam=diam*360.0D0/6378.137D0     ! convert to degrees
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.34) THEN ! min, max time offset (hours)
        no=no+1
        if (no.gt.maxnf) then
           write(*,*)myname,'Too many experiment forecast intervals.',no
           irc=232
           return
        end if
        dat250x=dat250
        buff250=nukehead(dat250,250)
        lenb=length(buff250,250,lend)
        read(buff250(1:lenb),*,iostat=irc)offmin(no)
        if (irc.ne.0) then
           write(*,*) myname,'Unable to read min time offset:',dat250x(1:lend)
           return
        end if
        buff250=nukehead(dat250,250)
        lenb=length(buff250,250,lend)
        read(buff250(1:lenb),*,iostat=irc)offmax(no)
        if (irc.ne.0) then
           offmax(no)=offmin(no)
!           write(*,*) myname,'Unable to read max time offset:',dat250x(1:lend)
!           return
        end if
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.35) THEN ! min, max ref forecast (hours)
        dat250x=dat250
        buff250=nukehead(dat250,250)
        lenb=length(buff250,250,lend)
        read(buff250(1:lenb),*,iostat=irc)dtmin(0)
        if (irc.ne.0) then
           write(*,*) myname,'Unable to read min time offset:',dat250x(1:lend)
           return
        end if
        buff250=nukehead(dat250,250)
        lenb=length(buff250,250,lend)
        read(buff250(1:lenb),*,iostat=irc)dtmax(0)
        if (irc.ne.0) then
           dtmax(0)=dtmin(0)
           irc=0
!           write(*,*) myname,'Unable to read max time offset:',dat250x(1:lend)
!           return
        end if
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.36) THEN ! min, max exp forecast (hours)
        nf=nf+1
        if (nf.gt.maxnf) then
           write(*,*)myname,'Too many experiment forecast intervals.',nf
           irc=237
           return
        end if
        dat250x=dat250
        buff250=nukehead(dat250,250)
        lenb=length(buff250,250,lend)
        read(buff250(1:lenb),*,iostat=irc)dtmin(nf)
        if (irc.ne.0) then
           write(*,*) myname,'Unable to read min time offset:',dat250x(1:lend)
           return
        end if
        buff250=nukehead(dat250,250)
        lenb=length(buff250,250,lend)
        read(buff250(1:lenb),*,iostat=irc)dtmax(nf)
        if (irc.ne.0) then
           dtmax(nf)=dtmin(nf)
           irc=0
           !           write(*,*) myname,'Unable to read max time offset:',dat250x(1:lend)
           !           return
        end if
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.40) THEN ! output file netcdf
        outnc250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.50) THEN ! output file s-plus
        outr250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.60) THEN ! use grid from file
        grid250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSE IF (LINE.NE.0) THEN ! LINE.EQ.0 IMPLIES SOMETIMES EOF
        WRITE(*,*) MYNAME,'System error, line not implemented:',LINE
        IRC=999
        RETURN
     ENDIF
  ENDDO
  ! 
  KODE=1
  CALL NUKEM(KODE,UNITI,HDR250,INTOUT,DAT250,ENDOFF,IRC)
  IF (IRC.NE.0) THEN
     WRITE(*,*) MYNAME,'Error return from NUKEM.',IRC
     RETURN
  END IF
  ! 
  GOTO 202
201 CONTINUE
  LEND=LENGTH(DAT250,250,1)
  CALL CHOP0(HDR250(LINE),250)
  LENH=LENGTH(HDR250(LINE),250,1)
  WRITE(*,*) MYNAME,'unable to read: ',DAT250(1:LEND),' (HDR='//HDR250(LINE)(1:LENH)//')'
  IRC=522
  RETURN
202 CONTINUE
  ! 
  WRITE(*,*) MYNAME,'----------------------------------------'
  !
  if (.not.lfldat(34)) then ! we must always have an offset
     no=1
     offmin(no) = 0.0D0
     offmax(no) = 0.0D0
  end if
  !
  if (.not.lfldat(36)) then ! if no forecast is available, we still need a forecast group...
     nf=0
     dtmin(nf) = 0.0D0
     dtmax(nf) = 0.0D0
     nf=nf+1
     dtmin(nf) = 0.0D0
     dtmax(nf) = 0.0D0
  end if
  !
  ifield=.true.
  igrid=.true.
  if (lfldat(60)) then ! read output grid into memory
     leng=length(grid250,250,10)
     write(*,*)myname,'Reading grid: ',grid250(1:leng)
     bok=.true.
     allocate(g,stat=irc)
     if (irc.ne.0) then
        write(*,*)myname,'Unable to allocate GRID.',irc
        return
     end if
     call openFile(grid250,g,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from openFile.',irc
        return
     end if
     if (bok) then
        call readInventory(g,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from readInventory.',irc
           return
        end if
        dum250=""
        call chop0(dum250,250)
        call checkContents(g,dum250,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from checkContents.',irc
           return
        end if
     end if
     if (bok) then
        ! read reference latitude field into memory
        v=>g%latid
        bbok=.true.
        if (v%lend.le.0) then ! not in memory
           call readRealData(v,bbok,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from readData.',irc
              return
           end if
        end if
        ! read reference longitude field into memory
        v=>g%lonid
        if (v%lend.le.0) then ! not in memory
           call readRealData(v,bbok,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from readData.',irc
              return
           end if
        end if
        call initGrid(out,g,irc)
        if (irc.ne.0) then
           write(*,*) myname,'Error return from INITGRID.',irc
        end if
        call closeFile(g,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from closeFile.',irc
           return
        end if
        igrid=.false.
     end if
     call clearInventory(g,irc)
     deallocate(g,stat=irc)
  end if
  !
  allocate(e(nexp),r(nref),stat=irc)
  if (irc.ne.0) then
     write(*,*)myname,'Unable to allocate INVENTORY',nexp, nref
     return
  end if

  cnte=0                     ! number of relevant experiment fields found
  cntr=0                     ! number of relevant reference fields found
  cntm=0                     ! number of matches

  ! loop over experiment files and make a register with sorted index
  do ee=1,nexp
     !
     ! open file and read inventory
     !
     lene=length(exp250(ee),250,10)
     write(*,*)myname,'Opening EXP: ',exp250(ee)(1:lene)
     call openFile(exp250(ee),e(ee),bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from openFile.',irc
        return
     end if
     if (bok) then
        eok(1)=eok(1)+1
     else
        erm(1)=erm(1)+1 ! unable to read file
     end if
     if (bok) then
        call readInventory(e(ee),bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from readInventory.',irc
           return
        end if
     end if
     !
     ! check that we have lat/lon/analysis-time/times/parameter
     !
     if (bok) then
        call checkContents(e(ee),epar250,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from checkContents.',irc
           return
        end if
        if (bok) then
           eok(2)=eok(2)+1
        else
           erm(2)=erm(2)+1 ! invalid contents
        end if
     end if
     !
     ! read analysis-time (must be present at this point)
     !
     if (bok) then
        v=>e(ee)%terid
        call readRealData(v,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from readData.',irc
           return
        end if
        ea2000=v%fd(1)
        call clearIndex(expLat)
        call clearIndex(expLat)
        call clearIndex(expLon)
        call clearIndex(expTim)
        call clearIndex(expPar)
        call clearIndex(expAll)
        expLat=>getDimIndex(e(ee)%latid,irc)
        expLon=>getDimIndex(e(ee)%lonid,irc)
        expTim=>getDimIndex(e(ee)%timid,irc)
        expPar=>getDimIndex(e(ee)%parid,irc)
        expAll=>getAllIndex(e(ee),irc)
        call clearIndex(expUnion)
        if (expTim%nrdim.eq.1) then
           expUnion => unionIndex(expTim,expPar,irc)
           if (expUnion%nrdim.eq.1) then
              ! determine number of times
              ntime=0
              if (associated(e(ee)%timid)) then
                 v=>e(ee)%timid
                 ntime=v%len
                 call readRealData(v,bok,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from readData.',irc
                    return
                 end if
              end if
           else !  parameter does not have any time dimension
              ntime=0
           end if
        else !  time variable does not have any dimension
           ntime=0
        end if
        ! loop over times
        bok=.false.
        do tt=min(1,ntime),ntime ! analysis if tt==0
           proceedb=.true.
           ! determine time of each parameter - field, use analysis if flag is set
           if (tt.eq.0) then ! use analysis time
              e2000=ea2000
           else
              v=>e(ee)%timid
              e2000=v%fd(tt)
           end if
           ! apply any experiment filter (from/to date etc.)
           if (proceedb) then
              dt=(e2000-ea2000)/3600.0D0
              group=0
              do jj=1,nf
                 if (dt.ge.dtmin(jj).and.dt.le.dtmax(jj)) then
                    group=jj
                 end if
              end do
              if (group.eq.0) proceedb=.false.
           end if
           !
           ! store file, time and dimension values in field-register
           !
           if (proceedb) then
              cnte=cnte+1
              if (cnte.gt.maxcnt) then
                 write(*,*) myname,'Too many fields:',cnte
                 irc=734
                 return
              end if
              fe(cnte)=ee  ! file number
              te(cnte)=tt  ! time number, 0==analysis
              je(cnte)=e2000 ! time value
              ae(cnte)=ea2000 ! forecast time
              ge(cnte)=group
              if (cnte.eq.1) then
                 estart=e2000
                 estop=e2000
              else
                 estart=min(estart,e2000)
                 estop=max(estop,e2000)
              end if
              bok=.true.
           end if
        end do              ! end loop over time
        if (bok) then
           eok(3)=eok(3)+1
        else
           erm(3)=erm(3)+1 ! no forecasts in valid forecast group
        end if
     end if                 ! proceed
     call closeFile(e(ee),irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from closeFile.',irc
        return
     end if
  end do                    ! end loop over experiment files
  ! write experiment statistics
  do ii=1,10
     pst(ii)=dfloat(erm(ii))/max(1.0d0,dfloat(erm(ii)+eok(ii)))*100
  end do
  WRITE(*,*)
  WRITE(*,998) MYNAME,                 'Available EXP files:   ', eok(1)+erm(1)
  IF (ERM(1).NE.0) WRITE(*,999) MYNAME,'Unable to read file:   ', -ERM(1),PST(1)
  IF (ERM(2).NE.0) WRITE(*,999) MYNAME,'Invalid contents:      ', -ERM(2),PST(2)
  IF (ERM(3).NE.0) WRITE(*,999) MYNAME,'No valid forecasts:    ', -ERM(3),PST(3)
  WRITE(*,*)   MYNAME,     '-------------------------------------------'
  WRITE(*,998) MYNAME,                 'Accepted EXP files:    ', eok(3)
  WRITE(*,998) MYNAME,                 'Accepted EXP fields:   ', cnte
  !
  ! sort experiment-field-register chronologically...
  !
  allocate(ind(cnte),stat=irc)
  do ee=1,cnte
     ind(ee)=ee
  end do
  eps=1.0D-10
  uniq=.false.
  call heapsort1(cnte,je,ncnt,cnte,ind,eps,uniq)
  if (ncnt.le.0) then
     write(*,*)myname,'No relevant experiment fields found.'
     irc=145
     return
  end if
  !
  ! loop over reference files
  !
  first=.true.
  do rr=1,nref
     ! open file and read inventory
     lenr=length(ref250(rr),250,10)
     write(*,*)myname,'Opening REF: ',ref250(rr)(1:lenr)
     call openFile(ref250(rr),r(rr),bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from READINVENTORY.',irc
        return
     end if
     ref=>r(rr)
     if (bok) then
        rok(1)=rok(1)+1
     else
        rrm(1)=rrm(1)+1 ! unable to open file
     end if
     if (bok) then
        call readInventory(ref,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from readInventory.',irc
           return
        end if
        if (bok) then
           rok(2)=rok(2)+1
        else
           rrm(2)=rrm(2)+1 ! unable to read inventory
        end if
     end if
     if (bok) then
        call checkContents(ref,rpar250,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from checkContents.',irc
           return
        end if
        if (bok) then
           rok(3)=rok(3)+1
        else
           rrm(3)=rrm(3)+1 ! contents error
        end if
     end if
     ! read analysis time (must be present)
     if (bok) then
        v=>ref%terid
        call readRealData(v,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from readData.',irc
           return
        end if
        if (bok) then
           ra2000=v%fd(1)
           rok(4)=rok(4)+1
        else
           rrm(4)=rrm(4)+1 ! unable to read analysis
        end if
     end if
     if (bok) then
        call clearIndex(refLat)
        call clearIndex(refLat)
        call clearIndex(refLon)
        call clearIndex(refTim)
        call clearIndex(refPar)
        call clearIndex(refAll)
        refLat=>getDimIndex(ref%latid,irc)
        refLon=>getDimIndex(ref%lonid,irc)
        refTim=>getDimIndex(ref%timid,irc)
        refPar=>getDimIndex(ref%parid,irc)
        refAll=>getAllIndex(ref,irc)
        if (refTim%nrdim.eq.1) then !  time does not have any (time) dimension
           call clearIndex(refUnion)
           refUnion => unionIndex(refTim,refPar,irc)
           if (refUnion%nrdim.eq.1) then !  parameter does not have any time dimension
              ! determine number of times in reference file
              ntime=0
              if (associated(ref%timid)) then
                 v=>ref%timid
                 ntime=v%len
                 call readRealData(v,bok,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from readData.',irc
                    return
                 end if
              end if
           else ! parameter does not contain time dimension
              ntime=0
           end if
        else ! no time dimension in file
           ntime=0
        end if
        ! loop over times
        bok=.false.
        do tt=min(1,ntime),ntime ! analysis if tt==0
           ! determine time of each parameter - field, use analysis if flag is set
           if (tt.eq.0) then ! use analysis time
              r2000=ra2000
           else
              v=>ref%timid
              r2000=v%fd(tt)
           end if
           do nn=1,no            ! offset number
              proceedb=.true.
              !
              ! apply any reference filter (from/to date etc.)
              !
              if (proceedb) then
                 ! find relevant experiment fiels...
                 dt=(r2000-ra2000)/3600.0D0
                 if (dt.lt.dtmin(0) .or. dt.gt.dtmax(0)) then
                    proceedb=.false.
                 else if (nn.eq.1) then
                    b2000=r2000/86400.0D0 - 10957.0D0  ! convert from seconds to days
                    call dj2000(b2000,yy,mm,dd,hh,mi,sec)
                    write (b24,'(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",F7.4)') yy,mm,dd,hh,mi,sec
                    write(*,*)myname,'Processing REF time: '//b24
                 end if
              end if
              !
              ! update start/stop times and count
              !
              if (proceedb) then
                 cntr=cntr+1
                 if (cntr.eq.1) then
                    rstart=r2000
                    rstop=r2000
                 else
                    rstart=min(rstart,r2000)
                    rstop=max(rstop,r2000)
                 end if
              end if
              !
              ! find releavnt experiment fields
              !
              if (proceedb) then
                 t2000=r2000+(offmax(nn)+offmin(nn))*1800.0D0
                 eps=abs(offmax(nn)-offmin(nn))*1800.0D0 + 1.0D-10
                 call binsearch(maxcnt,ncnt,je,ind,t2000,eps,left,right)
                 proceedb=(left.le.right) ! we now know that there will be some calculations...
              end if
              !
              ! read fields into memory
              !
              if (proceedb) then
                 first=.false.
                 ! read reference latitude field into memory
                 v=>ref%latid
                 bbok=.true.
                 if (v%lend.le.0) then ! not in memory
                    call readRealData(v,bbok,irc)
                    if (irc.ne.0) then
                       write(*,*)myname,'Error return from readData.',irc
                       return
                    end if
                 end if
                 ! read reference longitude field into memory
                 v=>ref%lonid
                 if (v%lend.le.0) then ! not in memory
                    call readRealData(v,bbok,irc)
                    if (irc.ne.0) then
                       write(*,*)myname,'Error return from readData.',irc
                       return
                    end if
                 end if
                 ! read parameter field into memory
                 v=>ref%parid
                 if (v%lend.le.0) then ! not in memory
                    call readRealData(v,bbok,irc)
                    if (irc.ne.0) then
                       write(*,*)myname,'Error return from readData.',irc
                       return
                    end if
                 end if
                 proceedb=bbok
              end if
              if (proceedb) then
                 b2000=t2000/86400.0D0 - 10957.0D0   ! convert from seconds to days
                 call dj2000(b2000,yy,mm,dd,hh,mi,sec)
                 write (b24,'(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",F7.4)') yy,mm,dd,hh,mi,sec
                 write(*,*)myname,'Processing GROUP time: '//b24
                 cntm=cntm+1
                 do jj=left,right
                    !
                    b2000=ae(ind(jj))/86400.0D0 - 10957.0D0   ! convert from seconds to days
                    call dj2000(b2000,yy,mm,dd,hh,mi,sec)
                    write (b24,'(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",F7.4)') yy,mm,dd,hh,mi,sec
                    x2000=je(ind(jj))/86400.0D0 - 10957.0D0   ! convert from seconds to days
                    call dj2000(x2000,yy,mm,dd,hh,mi,sec)
                    write (x24,'(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",F7.4)') yy,mm,dd,hh,mi,sec
                    write(*,'(X,A,"Processing EXP time:",A," + ",F5.2,"h = ",A," FcGroup=",I2.2," EmRtime=",I2.2)') &
                         & myname,b24,(je(ind(jj))-ae(ind(jj)))/3600.0D0,x24,ge(ind(jj)),nn
                    !(je(ind(jj))-ae(ind(jj)))/3600.0D0,x24,ge(ind(jj))
                    fexp=fe(ind(jj))  ! file inventory index
                    group=ge(ind(jj))
                    ! clear old experiment data not used any more
                    bbok=.true.
                    call reopenFile(e(fexp),bbok,irc)
                    if (irc.ne.0) then
                       write(*,*)myname,'Unable to reopen.'
                       irc=845
                       return
                    end if
                    exp => e(fexp)    ! experiment file inventory
                    ! read experiment latitude field into memory
                    if (bbok) then
                       v=>exp%latid
                       if (v%lend.le.0) then ! not in memory
                          call readRealData(v,bbok,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from readData.',irc
                             return
                          end if
                       end if
                       ! read experiment longitude field into memory
                       v=>exp%lonid
                       if (v%lend.le.0) then ! not in memory
                          call readRealData(v,bbok,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from readData.',irc
                             return
                          end if
                       end if
                       ! read parameter field into memory
                       v=>exp%parid
                       if (v%lend.le.0) then ! not in memory
                          call readRealData(v,bbok,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from readData.',irc
                             return
                          end if
                       end if
                    end if
                    if (bbok) then
                       !
                       ! get latlon dimensions in experiment
                       !
                       call clearIndex(expLat)
                       call clearIndex(expLon)
                       call clearIndex(expTim)
                       call clearIndex(expPar)
                       call clearIndex(expAll)
                       expLat=>getDimIndex(exp%latid,irc)
                       expLon=>getDimIndex(exp%lonid,irc)
                       expTim=>getDimIndex(exp%timid,irc)
                       expPar=>getDimIndex(exp%parid,irc)
                       expAll=>getAllIndex(exp,irc)
                       !
                       ! set time position
                       !
                       call setPos(exp,expTim,te(ind(jj))) ! not used if parameter is independent of time
                       call setPos(ref,refTim,max(1,tt))   ! not used if parameter is independent of time
                       !
                       ! make index of outer dimensions in reference
                       !
                       refOut => copyIndex(refPar)
                       call delIndex(refOut,refLat)
                       call delIndex(refOut,refLon)
                       call delIndex(refOut,refTim)
                       !
                       ! make index of inner dimensions (lat/lon dimensions)
                       !
                       refLatlon => copyIndex(refLon)
                       call addIndex(refLatlon,refLat)
                       !
                       expLatlon => copyIndex(expLon)
                       call addIndex(expLatlon,expLat)
                       !
                       ! check that lat and lon are 2-D
                       !
                       if (refLatlon%nrdim.ne.2) then
                          write(*,*)myname,'Invalid number of REF lat/lon dimensions:',refLatlon%nrdim
                          irc=645
                          return
                       end if
                       refI=>getdim(refLatlon,1)
                       refJ=>getdim(refLatlon,2)
                       !
                       if (expLatlon%nrdim.ne.2) then
                          write(*,*)myname,'Invalid number of EXP lat/lon dimensions:',expLatlon%nrdim
                          irc=646
                          return
                       end if
                       expI=>getdim(expLatlon,1)
                       expJ=>getdim(expLatlon,2)
                       !
                       ! store output grid
                       !
                       if (ifield) then
                          ifield=.false.
                          if (igrid) then
                             call initGrid(out,ref,irc)
                             if (irc.ne.0) then
                                write(*,*)myname,'Error return from initGrid.',irc
                                return
                             end if
                             igrid=.false.
                          end if
                          !
                          call clearIndex(outLat)
                          call clearIndex(outLon)
                          outLat=>getDimIndex(out%latid,irc)
                          outLon=>getDimIndex(out%lonid,irc)
                          outLatlon => copyIndex(outLon)
                          call addIndex(outLatlon,outLat)
                          
                          outLatlonFinal=>copyIndex(outLatlon)
                          if (lfldat(36)) then ! forecast group specified
                             outExpGroup=>createDimension(out,"experimentForecastGroup",nf)
                             call addIndex(outLatlonFinal,outExpGroup)
                          end if
                          outOffset=>createDimension(out,"timeOffsetGroup",no)
                          call addIndex(outLatlonFinal,outOffset)

                          
                          if (outLatlon%nrdim.ne.2) then
                             write(*,*)myname,'Invalid number of OUT lat/lon dimensions:',outLatlon%nrdim
                             irc=645
                             return
                          end if
                          !
                          ! insert statistics fields
                          !
                          sr1=>insertField(out,"sr1",outLatlonFinal,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from insertField.',irc
                             return
                          end if
                          se1=>insertField(out,"se1",outLatlonFinal,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from insertField.',irc
                             return
                          end if
                          sr2=>insertField(out,"sr2",outLatlonFinal,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from insertField.',irc
                             return
                          end if
                          se2=>insertField(out,"se2",outLatlonFinal,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from insertField.',irc
                             return
                          end if
                          ser=>insertField(out,"ser",outLatlonFinal,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from insertField.',irc
                             return
                          end if
                          sae=>insertField(out,"mae",outLatlonFinal,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from insertField.',irc
                             return
                          end if


                          if (lfldat(36)) then ! forecast group specified
                             gStart=>insertField(out,"forecastStart",outExpGroup,irc)
                             if (irc.ne.0) then
                                write(*,*)myname,'Error return from insertField (gStart).',irc
                                return
                             end if
                             gEnd=>insertField(out,"forecastEnd",outExpGroup,irc)
                             if (irc.ne.0) then
                                write(*,*)myname,'Error return from insertField (gEnd).',irc
                                return
                             end if
                             call addAttribute(gStart,"long_name", "start of forecasts")
                             call addAttribute(gStart,"units",     "hour")
                             call addAttribute(gEnd,  "long_name", "end of forecasts")
                             call addAttribute(gEnd,  "units",     "hour")
                          end if
                          oStart=>insertField(out,"timeOffsetStart",outOffset,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from insertField (oStart).',irc
                             return
                          end if
                          oEnd=>insertField(out,"timeOffsetEnd",outOffset,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from insertField (oEnd).',irc
                             return
                          end if
                          call addAttribute(oStart,"long_name", "start of offset")
                          call addAttribute(oStart,"units",     "hour")
                          call addAttribute(oEnd,  "long_name", "end of offset")
                          call addAttribute(oEnd,  "units",     "hour")
                          sum=>insertField(out,"cnt",outLatlonFinal,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from insertField.',irc
                             return
                          end if

                          call addAttribute(sr1,"long_name","mean reference value")
                          call addAttribute(se1,"long_name","mean experiment value")
                          call addAttribute(sr2,"long_name","mean squared reference value")
                          call addAttribute(se2,"long_name","mean squared experiment value")
                          call addAttribute(ser,"long_name","mean experiment*reference value")
                          call addAttribute(sum,"long_name","number of matches")
                          call addAttribute(sae,"long_name","mean absolute error (reference-experiment)")

                          if (associated(out%gridid)) then
                             gmap=out%gridid%var250
                             leng=out%gridid%lenv
                             if (leng.gt.0) then
                                call addAttribute(sr1,"grid_mapping",gmap(1:leng))
                                call addAttribute(se1,"grid_mapping",gmap(1:leng))
                                call addAttribute(sr2,"grid_mapping",gmap(1:leng))
                                call addAttribute(se2,"grid_mapping",gmap(1:leng))
                                call addAttribute(ser,"grid_mapping",gmap(1:leng))
                                call addAttribute(sum,"grid_mapping",gmap(1:leng))
                                call addAttribute(sae,"grid_mapping",gmap(1:leng))
                             end if
                          end if

                          refwgt => makeWeight4(ref%nrdim,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from makeWeight4.',irc
                             return
                          end if
                          expwgt => makeWeight4(exp%nrdim,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from makeWeight4.',irc
                             return
                          end if
                          allocate(rpos,stat=irc)
                          call allocatePos(ref%nrdim,rpos,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from allocatePos.',irc
                             return
                          end if
                          allocate(epos,stat=irc)
                          call allocatePos(exp%nrdim,epos,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from allocatePos.',irc
                             return
                          end if
                          ! make group start/end values...
                          if (lfldat(36)) then ! forecast group specified
                             do ii=1,nf
                                call setPos(out,outExpGroup,ii)
                                call setValue(gStart, dtmin(ii), irc)
                                call setValue(gEnd,   dtmax(ii), irc)
                             end do
                          end if
                          ! make offset start/end values...
                          do ii=1,no
                             call setPos(out,outOffset,ii)
                             call setValue(oStart, offmin(ii), irc)
                             call setValue(oEnd,   offmax(ii), irc)
                          end do
                       end if
                       !
                       ! Check that reference grid equals output grid
                       !
                       call checkGrid(out,ref,irc)
                       if (irc.ne.0) then
                          write(*,*)myname,'Error return from checkGrid.',irc
                          return
                       end if
                       !
                       ! set group position
                       !
                       if (lfldat(36)) then ! forecast group specified
                          call setPos(out,outExpGroup,ge(ind(jj)))
                       end if
                       !
                       ! set offset position
                       !
                       call setPos(out,outOffset,nn)
                       !
                       ! make sub index for footprint
                       !
                       expFoot=>copyIndex(expLatlon)
                       refFoot=>copyIndex(refLatlon)
                       !
                       ! loop over output grid and match grid points
                       !
                       cnt=0
                       ! write(*,*) myname,'Initialising outer REF'
                       call resetIndex(refOut,irc)
                       do while (increment(ref,refOut,irc))
                          !
                           write(*,*) myname,'Initialising inner REF'
                          call resetIndex(outLatlon,irc)

                          do while (increment(out,outLatlon,irc))
                             !
                             !                          if (mod(cnt,10000).eq.0) then
                             !                             write(*,*) myname,'Inside loop...'
                             !                             call printpos(out%pos)
                             !                          end if
                             cnt=cnt+1
                             rcnt=0.0D0
                             ecnt=0.0D0
                             if (out%nrdim .ne. out%pos%nrdim) then
                                write(*,*) myname,'Pos mismatch.',out%nrdim,out%pos%nrdim
                                irc=999
                                return
                             end if
                             bbok=.true.         ! data looks fine so far...
                             if (bbok) then
                                ! interpolate reference field
                                call interpolate2D(out,ref,refI,refJ,refwgt,bbok,irc)
                                if (irc.ne.0) then
                                   write(*,*)myname,'Error return from interpolate2d.',irc
                                   return
                                end if
                                ! call printpos(ref%pos)
                             end if
                             if (bbok) then
                                ! interpolate experiment field
                                call interpolate2D(out,exp,expI,expJ,expwgt,bbok,irc)
                                if (irc.ne.0) then
                                   write(*,*)myname,'Error return from interpolate2d.',irc
                                   return
                                end if
                                ! call printpos(exp%pos)
                             end if
                             ! handle footprint
                             if (bbok) then
                                !
                                ! reference value
                                ! ...we calculate the radius from lower left corner..
                                !
                                call copyPos(rpos,ref%pos,irc) ! make copy of center position
                                if (getFootIndex(ref,refLatlon,refFoot,latc,lonc,diam,irc)) then
                                   rvalue=0.0D0
                                   call resetIndex(refFoot,irc)
                                   do while (incrementFoot(ref,refFoot,latc,lonc,diam,irc))
                                      val=valueWeighted(ref%parid,refwgt,irc)
                                      if (val.eq.nf_fill_double) then
                                         rcnt=-1
                                      else if (rcnt.ge.0.0D0) then
                                         rvalue=rvalue+val
                                         rcnt=rcnt+1.0D0
                                      end if
                                   end do
                                   if (rcnt.gt.0.0D0) then
                                      rvalue=rscale*rvalue/max(1.0D0,rcnt)
                                   else
                                      rvalue=nf_fill_double
                                   end if
                                else
                                   rvalue=nf_fill_double
                                end if
                                call copyPos(ref%pos,rpos,irc) ! recover position
                                !
                                ! experiment value
                                ! ...we calculate the radius from lower left corner..
                                !
                                call copyPos(epos,exp%pos,irc) ! make copy of center position
                                if (getFootIndex(exp,expLatlon,expFoot,latc,lonc,diam,irc)) then
                                   evalue=0.0D0
                                   call resetIndex(expFoot,irc)
                                   do while (incrementFoot(exp,expFoot,latc,lonc,diam,irc))
                                      val=valueWeighted(exp%parid,expwgt,irc)
                                      if (val.eq.nf_fill_double) then
                                         ecnt=-1
                                      else if (ecnt.ge.0.0D0) then
                                         evalue=evalue+val
                                         ecnt=ecnt+1.0D0
                                      end if
                                   end do
                                   if (ecnt.gt.0.0D0) then
                                      evalue=escale*evalue/max(1.0D0,ecnt)
                                   else
                                      evalue=nf_fill_double
                                   end if
                                else
                                   evalue=nf_fill_double
                                end if
                                call copyPos(exp%pos,epos,irc) ! recover position
                             end if
                             if (bbok) then
                                bbok=(ecnt.gt.0 .and. rcnt.gt.0) ! we are in business...
                             end if
                             if (bbok) then
                                
                                if (abs(evalue).gt.1000.0.or.abs(rvalue).gt.1000.0) write(*,'(X,A,A,2(X,F17.6,X,I3))') &
                                     & myname,'Matched:',rvalue,nint(rcnt),evalue,nint(ecnt)
                                
                                !
                                ! update fields
                                !
                                call addValue(sr1, rvalue,       irc)
                                if (irc.ne.0) then
                                   write(*,*) myname,'Error return from addValue.',irc
                                   return
                                end if
                                call addValue(se1, evalue,       irc)
                                if (irc.ne.0) then
                                   write(*,*) myname,'Error return from addValue.',irc
                                   return
                                end if
                                call addValue(sr2, rvalue*rvalue,irc)
                                if (irc.ne.0) then
                                   write(*,*) myname,'Error return from addValue.',irc
                                   return
                                end if
                                call addValue(se2, evalue*evalue,irc)
                                if (irc.ne.0) then
                                   write(*,*) myname,'Error return from addValue.',irc
                                   return
                                end if
                                call addValue(ser, evalue*rvalue,irc)
                                if (irc.ne.0) then
                                   write(*,*) myname,'Error return from addValue.',irc
                                   return
                                end if
!                                if (ge(ind(jj)).eq.1.and.abs(rvalue-evalue).gt.1.0D-10) then
!                                   write(*,*) 'Invalid abs error:',rvalue,evalue
!                                end if
                                call addValue(sae, abs(rvalue-evalue),irc)
                                if (irc.ne.0) then
                                   write(*,*) myname,'Error return from addValue.',irc
                                   return
                                end if
                                call addValue(sum, 1.0D0,        irc)
                                if (irc.ne.0) then
                                   write(*,*) myname,'Error return from addValue.',irc
                                   return
                                end if
                                bok=.true.
                             end if
                          end do
                       end do    ! outer loop
                       call closeFile(exp,irc)
                       if (irc.ne.0) then
                          write(*,*)myname,'Error return from closeFile.',irc
                          return
                       end if
                    end if     ! bbok
                 end do        ! left,right
              end if    ! proceedb
           end do       ! nn
        end do    ! tt        
        if (bok) then
           rok(5)=rok(5)+1
        else
           rrm(5)=rrm(5)+1 ! no match found
        end if
     end if       ! bok
     call clearData(ref,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from cleanMemory.',irc
        return
     end if
     call closeFile(ref,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from closeFile.',irc
        return
     end if
  end do             ! nref
  !
  ! write output files...
  !
  if (lfldat(40)) then ! write netcdf file
     if (out%initialised .and. .not.ifield) then
        !
        ! normalise
        !
        call divideField(sr1,sum)
        call divideField(sr2,sum)
        call divideField(se1,sum)
        call divideField(se2,sum)
        call divideField(ser,sum)
        call divideField(sae,sum)
        !
        ! get experiment stdve and mean
        !
        stdve=>insertField(out,"stdv_exp",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        meane=>insertField(out,"mean_exp",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        sese=>insertField(out,"mean2_exp",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        call addAttribute(stdve,"long_name","standard deviation (experiment)")
        call addAttribute(meane,"long_name","mean (experiment)")
        !
        call copyfield(stdve,se2)
        call copyField(meane,se1)
        call copyField(sese,se1)
        call multiplyField(sese,sese)
        call subtractField(stdve,sese)
        call sqrtField(stdve)
        !
        ! get reference stdve and mean
        !
        stdvr=>insertField(out,"stdv_ref",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        meanr=>insertField(out,"mean_ref",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        srsr=>insertField(out,"mean2_ref",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        call addAttribute(stdvr,"long_name","standard deviation (reference)")
        call addAttribute(meanr,"long_name","mean (reference)")
        !
        call copyfield(stdvr,sr2)
        call copyField(meanr,sr1)
        call copyField(srsr,sr1)
        call multiplyField(srsr,srsr)
        call subtractField(stdvr,srsr)
        call sqrtField(stdvr)
        !
        ! get stdv and bias between experiment and reference
        !
        stdv=>insertField(out,"stdv_E-R",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        rms=>insertField(out,"rms_E-R",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        bias=>insertField(out,"bias_E-R",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        sesr=>insertField(out,"sesr",outLatlonFinal,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from insertField.',irc
           return
        end if
        call addAttribute(stdv,"long_name","standard deviation (reference-experiment)")
        call addAttribute(rms,"long_name","root mean squared (reference-experiment)")
        call addAttribute(bias,"long_name","bias (reference-experiment)")
        !
        call copyField(bias,sr1)
        call subtractField(bias,se1)
        !
        call copyField(sesr,se1)
        call multiplyField(sesr,sr1)
        !
        call copyfield(stdv,ser)
        call subtractField(stdv,sesr)
        call scalefield(stdv,-2.0D0)
        call addField(stdv,sr2)
        call addField(stdv,se2)
        call subtractField(stdv,srsr)
        call subtractField(stdv,sese)
        call sqrtField(stdv)
        call roundField(stdv,1.0D-5)
        !
        call copyfield(rms,ser)
        call scalefield(rms,-2.0D0)
        call addField(rms,sr2)
        call addField(rms,se2)
        call sqrtField(rms)
        call roundField(rms,1.0D-5)
        !
        call removeField(sese)
        call removeField(srsr)
        call removeField(sesr)
        !
        ! recover non-real data
        !
        call unmakeRealData(out%latid,irc)
        call unmakeRealData(out%lonid,irc)
        !
        ! remove unused dimensions and write to file
        !
        if (associated(out%gridid)) then
           gmap=out%gridid%var250
           leng=out%gridid%lenv
           if (leng.gt.0) then
              call addAttribute(stdve,"grid_mapping",gmap(1:leng))
              call addAttribute(meane,"grid_mapping",gmap(1:leng))
              call addAttribute(stdvr,"grid_mapping",gmap(1:leng))
              call addAttribute(meanr,"grid_mapping",gmap(1:leng))
              call addAttribute(stdv,"grid_mapping",gmap(1:leng))
              call addAttribute(bias,"grid_mapping",gmap(1:leng))
           end if
        end if
        
        call compressDimensions(out)
        call writeNcOut(outnc250,out,irc)
        if (irc.ne.0) then
           write(*,*) myname,'Error return from writeNcOut.',irc
           return
        end if
     else
        write(*,*) myname,'No output data available.'
     end if
  end if
  if (lfldat(50)) then ! write splus file
     !
  end if
  !
  ! write summary   
  !
  !convert to days since 2000/01/01
  estart=estart/86400.0D0 - 10957.0D0   ! convert from seconds to days
  estop=estop/86400.0D0 - 10957.0D0   ! convert from seconds to days
  rstart=rstart/86400.0D0 - 10957.0D0   ! convert from seconds to days
  rstop=rstop/86400.0D0 - 10957.0D0   ! convert from seconds to days

  call dj2000(estart,yy,mm,dd,hh,mi,sec)
  write (estart24,'(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",F7.4)') yy,mm,dd,hh,mi,sec

  call dj2000(estop,yy,mm,dd,hh,mi,sec)
  write (estop24,'(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",F7.4)') yy,mm,dd,hh,mi,sec

  call dj2000(rstart,yy,mm,dd,hh,mi,sec)
  write (rstart24,'(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",F7.4)') yy,mm,dd,hh,mi,sec

  call dj2000(rstop,yy,mm,dd,hh,mi,sec)
  write (rstop24,'(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",F7.4)') yy,mm,dd,hh,mi,sec

  write(*,'(X,A,A,I6)')myname, 'Found EXP data from  '//estart24//'  to  '//estop24//',  fields=',cnte
  write(*,'(X,A,A,I6)')myname, 'Found REF data from  '//rstart24//'  to  '//rstop24//',  fields=',cntr
  !
  do ii=1,10
     pst(ii)=dfloat(rrm(ii))/max(1.0d0,dfloat(rrm(ii)+rok(ii)))*100
  end do
  WRITE(*,*)
  WRITE(*,998) MYNAME,                  'Available REF files:   ', rok(1)+rrm(1)
  IF (RRM(1).NE.0)  WRITE(*,999) MYNAME,'Unable to read file:   ', -RRM(1),PST(1)
  IF (RRM(2).NE.0)  WRITE(*,999) MYNAME,'Invalid inventory:     ', -RRM(2),PST(2)
  IF (RRM(3).NE.0)  WRITE(*,999) MYNAME,'Invalid contents:      ', -RRM(3),PST(3)
  IF (RRM(4).NE.0)  WRITE(*,999) MYNAME,'No valid analysis:     ', -RRM(4),PST(4)
  IF (RRM(5).NE.0)  WRITE(*,999) MYNAME,'No matching EXP:       ', -RRM(5),PST(5)
  WRITE(*,*)   MYNAME,     '-------------------------------------------'
  WRITE(*,998) MYNAME,                  'Processed REF files:   ', rok(5)
  WRITE(*,998) MYNAME,                  'Accepted REF fields:   ', cntr
  WRITE(*,*)   MYNAME,     '-------------------------------------------'
  WRITE(*,998) MYNAME,                  'EXP-REF-field matches: ', cntm
  !
999 FORMAT(X,A12,X,A23,I10,' (',F5.1,'%)')
998 FORMAT(X,A12,X,A23,I10)
997 FORMAT(X,A12,X,A,10(X,I10))
  !
  ! exit with error if no matches were found...
  !
  if (cntm.eq.0) irc=143
  return

contains
  ! 
  ! read netcdf inventory
  ! 
  subroutine openFile(fn250,f,biok,irc)
    implicit none
    ! 
    character*250 fn250
    type (inventory) f
    logical biok
    integer irc
    ! 
    character*250 buff250
    integer lenb, length
    external length
    !
    f%fn250=fn250
    f%lenf=length(f%fn250,250,10)
    ! open nceriment file
    proceed=.true.
    ! -> open NETCDF file
    chunksizehint= 1024*1024*1024
    ! chunksizehint= 1024*1024*1024
    RET = NF__OPEN(f%fn250(1:f%lenf),nf_nowrite,chunksizehint,f%ncid)
    ! RET = NF_OPEN(IN250(1:NCLENI),nf_nowrite,f%NCID)
    if (ret .ne. NF_NOERR) then
       write(buff250,"(A)") " Unable to open: "//f%fn250(1:f%lenf)
       call chop0(buff250,250)
       lenb=length(buff250,250,10)
       write(*,*)'openFile '//buff250(1:lenb)
       biok=.false.
       return
    else
       write(buff250,"(A)") "  Opened: "//f%fn250(1:f%lenf)
       call chop0(buff250,250)
       lenb=length(buff250,250,10)
       ! write(*,*)'openFile '//buff250(1:lenb)
       biok=.true.
       f%opened=.true.
    endif
  end subroutine openFile

  subroutine reopenFile(f,biok,irc)
    implicit none
    ! 
    character*250 fn250
    type (inventory) f
    logical biok
    integer irc
    ! 
    ! open file
    proceed=.true.
    ! -> open NETCDF file
    chunksizehint= 1024*1024*1024
    ! chunksizehint= 1024*1024*1024
    RET = NF__OPEN(f%fn250(1:f%lenf),nf_nowrite,chunksizehint,f%ncid)
    ! RET = NF_OPEN(IN250(1:NCLENI),nf_nowrite,f%NCID)
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_OPEN:",nf_strerror(ret),f%ncid
       write(buff250,"(A)") " Unable to open: "//f%fn250(1:f%lenf)
       call chop0(buff250,250)
       lenb=length(buff250,250,10)
       write(*,*)myname,buff250(1:lenb)
       biok=.false.
       return
    else
       write(buff250,"(A)") "  Re-opening: "//f%fn250(1:f%lenf)
       call chop0(buff250,250)
       lenb=length(buff250,250,10)
       write(*,*)myname,buff250(1:lenb)
       biok=.true.
       f%opened=.true.
    endif
    v=>f%firstVariable%next
    do while(.not.associated(v,target=f%lastVariable))
       v%ncid=f%ncid
       v=>v%next
    end do
  end subroutine reopenFile


  subroutine readInventory(f,biok,irc)
    implicit none
    ! 
    type (inventory) :: f
    logical biok
    integer irc
    ! 
    type(variable), pointer :: v
    type(attribute), pointer :: a
    character*250 buff250
    integer lenb, length
    integer kk,ll,mm
    external length
    !
    ! get number of dimension, variables, attributes, unlimited dimension id
    RET = NF_INQ(f%ncid, f%nrdim, f%nrvar, f%nrgatt, f%unlimdimid)
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_INQ:",nf_strerror(ret),f%ncid
       irc=790
       return
    end if
    ! allocate dimensions in file
    allocate(f%dim250(f%nrdim),f%lend(f%nrdim),stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,'Unable to allocate F-dims.',irc
       return
    end if
    allocate(f%pos,stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,'Unable to allocate F-pos.',irc
       return
    end if
    f%pos%nrdim=f%nrdim
    allocate(f%pos%sta(max(1,f%pos%nrdim)),f%pos%lim(max(1,f%pos%nrdim)),f%pos%pos(max(1,f%pos%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,'Unable to allocate F-pos-pos.',irc
       return
    end if
    f%pos%pos(1)=1
    f%pos%lim(1)=1
    f%pos%sta(1)=1
    do ll=1,f%pos%nrdim
       f%pos%pos(ll)=1
       f%pos%lim(ll)=1
       f%pos%sta(ll)=1
    end do
    ! 
    ! -> store dimension names
    do kk=1,f%nrdim
       call readDimension(f%ncid,kk,f%dim250(kk),f%lend(kk),f%pos%lim(kk),irc)
       if (irc.ne.0) then
          write(*,*) 'READINVENTORY Error return from READDIMENSION.',irc
          return
       end if
    end do
    !
    allocate(f%firstVariable,f%lastVariable, stat=irc)
    f%firstVariable%next=>f%lastVariable
    f%lastVariable%prev=>f%firstVariable
    ! -> store variables
    do kk=1,f%nrvar
       allocate(v,stat=irc)
       call readVariableAttributes(f%ncid,kk,v,irc)
       if (irc.ne.0) then
          write(*,*) 'READINVENTORY Error return from READVARIABLEATTRIBUTES.',irc
          return
       end if
       call setVariableDimension(f,v,irc)
       if (irc.ne.0) then
          write(*,*) 'READINVENTORY Error return from SetVariableDimension.',irc
          return
       end if
       v%prev=>f%lastVariable%prev
       v%next=>f%lastVariable
       f%lastVariable%prev%next=>v
       f%lastVariable%prev=>v
       nullify(v)
    end do
    !
    allocate(f%firstAttribute,f%lastAttribute, stat=irc)
    f%firstAttribute%next=>f%lastAttribute
    f%lastAttribute%prev=>f%firstAttribute
    !
    ! -> store global attributes
    do kk=1,f%nrgatt
       allocate(a,stat=irc)
       call readAttribute(f%ncid,nf_global,kk,a,irc)
       if (irc.ne.0) then
          write(*,*) 'READINVENTORY Error return from READATTRIBUTE.',irc
          return
       end if
       a%prev=>f%lastAttribute%prev
       a%next=>f%lastAttribute
       f%lastAttribute%prev%next=>a
       f%lastAttribute%prev=>a
       nullify(a)
    end do
    ! reset indexes...
    nullify(f%latid)
    nullify(f%lonid)
    nullify(f%timid)
    nullify(f%terid)
    nullify(f%parid)
    nullify(f%xcoid)
    nullify(f%ycoid)
    nullify(f%gridid)
    f%initialised=.true.
  end subroutine readInventory

  subroutine clearInventory(f,irc)
    implicit none
    type (inventory) f
    integer irc
    ! 
    type(variable),pointer :: v,vx
    type(attribute),pointer :: a,ax
    !
    if (f%initialised) then
       v=>f%firstVariable%next
       do while(.not.associated(v,target=f%lastVariable))
          call clearVariable(v,irc)
          vx=>v
          v=>v%next
          deallocate(vx,stat=irc)
       end do
       a=>f%firstAttribute%next
       do while(.not.associated(a,target=f%lastAttribute))
          call clearAttribute(a,irc)
          ax=>a
          a=>a%next
          deallocate(ax,stat=irc)
       end do
       if (associated(f%dim250)) deallocate(f%dim250,stat=irc)
       if (associated(f%lend)) deallocate(f%lend,stat=irc)
       if (associated(f%pos%lim)) deallocate(f%pos%lim,stat=irc)
       if (associated(f%pos%pos)) deallocate(f%pos%pos,stat=irc)
       if (associated(f%pos)) deallocate(f%pos,stat=irc)
       f%initialised=.false.
       f%opened=.false.
    end if
  end subroutine clearInventory


  subroutine clearVariable(v,irc)
    implicit none
    type (Variable) v
    integer irc
    type(attribute), pointer :: a
    if (v%initialised) then
       a=>v%firstAttribute%next
       do while (.not.associated(a,target=v%lastAttribute))
          call clearAttribute(a,irc)
          a=>a%next
       end do
       if (associated(v%ind)) deallocate(v%ind,stat=irc)
       if (associated(v%lim)) deallocate(v%lim,stat=irc)
       if (associated(v%sta)) deallocate(v%sta,stat=irc)
       if (associated(v%fc)) deallocate(v%f1,stat=irc)
       if (associated(v%f1)) deallocate(v%f2,stat=irc)
       if (associated(v%f2)) deallocate(v%f2,stat=irc)
       if (associated(v%f4)) deallocate(v%f4,stat=irc)
       if (associated(v%fr)) deallocate(v%fr,stat=irc)
       if (associated(v%fd)) deallocate(v%fd,stat=irc)
       v%initialised=.false.
    end if
  end subroutine clearVariable


  subroutine clearAttribute(a,irc)
    implicit none
    type (attribute) a
    integer irc
    if (a%initialised) then
       if (associated(a%ac)) deallocate(a%a1,stat=irc)
       if (associated(a%a1)) deallocate(a%a2,stat=irc)
       if (associated(a%a2)) deallocate(a%a2,stat=irc)
       if (associated(a%a4)) deallocate(a%a4,stat=irc)
       if (associated(a%ar)) deallocate(a%ar,stat=irc)
       if (associated(a%ad)) deallocate(a%ad,stat=irc)
       a%initialised=.false.
    end if
  end subroutine clearAttribute


  subroutine readDimension(ncid,dimid,dim250,lend,lim,irc)
    implicit none
    integer ncid
    integer dimid
    character*250 dim250
    integer lend
    integer lim
    integer irc
    integer ret
    character*250 buff250
    RET = NF_INQ_DIM(ncid,dimid,dim250,lim)
    if (ret .ne. NF_NOERR) then
       write(*,*) "READDIMENSION ERROR from NF_INQ_DIM:",nf_strerror(ret)
       irc=801
       return
    end if
    call chop0(dim250,250)
    lend=length(dim250,250,10)
    write(buff250,'(X,A20,5X,X,I5)')dim250(1:lend),lim
    call chop0(buff250,250)
    lenb=length(buff250,250,10)
    ! write(*,'(10X,A)') buff250(1:lenb)
    return
  end subroutine readDimension

  subroutine readVariableAttributes(ncid,varid,v,irc)
    implicit none
    integer ncid
    integer varid
    type(variable), target :: v
    integer irc
    type(attribute), pointer :: a
    integer jj
    v%ncid=ncid
    v%varid=varid
    ret = NF_INQ_VARNDIMS (v%ncid, v%varid, v%nrdim);
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_INQ_VARNDIMS:",nf_strerror(ret)
       irc=802
       return
    end if
    allocate(v%ind(v%nrdim), stat=irc)
    if (irc.ne.0) then
       write(*,*) 'READVARIABLE Unable to allocate DIMENSIONS.',irc
       return
    end if
    RET = NF_INQ_VAR(v%ncid,v%varid,v%var250,v%type,v%nrdim,v%ind,v%nratt)
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_INQ_VAR:",nf_strerror(ret)
       irc=802
       return
    end if
    call chop(v%var250,250)
    v%lenv=length(v%var250,250,10)
    allocate(v%firstAttribute,v%lastAttribute,v%fillAttribute, stat=irc)
    v%firstAttribute%next=>v%lastAttribute
    v%lastAttribute%prev=>v%firstAttribute
    ! read attributes
    do jj=1,v%nratt
       allocate(a,stat=irc)
       if (irc.ne.0) then
          write(*,*) 'READVARIABLEATTRIBUTES Unable to allocate ATTRIBUTE.',irc
          return
       end if
       call readAttribute(v%ncid,v%varid,jj,a,irc)
       if (irc.ne.0) then
          write(*,*) 'READVARIABLEATTRIBUTES Error return from READATTRIBUTE.',irc
          return
       end if
       a%v => v
       a%prev=>v%lastAttribute%prev
       a%next=>v%lastAttribute
       v%lastAttribute%prev%next=>a
       v%lastAttribute%prev=>a
       nullify(a)
    end do
    v%initialised=.true.
  end subroutine readVariableAttributes


  subroutine setVariableDimension(f,v,irc)
    implicit none
    type(inventory), target :: f
    type(variable) :: v
    integer irc
    integer jj
    v%f=>f
    if (associated(v%lim)) deallocate(v%lim,stat=irc)
    allocate(v%lim(max(1,v%nrdim)),stat=irc)
    if (associated(v%sta)) deallocate(v%sta,stat=irc)
    allocate(v%sta(max(1,v%nrdim)),stat=irc)
    v%len=1
    v%lim(1)=1
    v%sta(1)=1
    do jj=1,v%nrdim
       v%lim(jj)=f%pos%lim(v%ind(jj))
       v%sta(jj)=1
       v%len=v%len*v%lim(jj)
    end do
  end subroutine setVariableDimension

  subroutine readAttribute(ncid,varid,attid,a,irc)
    implicit none
    integer ncid
    integer varid
    integer attid
    type(attribute) a
    integer irc
    integer length
    external length
    a%ncid=ncid
    a%varid=varid
    a%attid=attid
    RET=NF_INQ_ATTNAME(ncid,varid,attid,a%att250)
    call chop(a%att250,250)
    a%lena=length(a%att250,250,10)
    if (ret .ne. NF_NOERR) then
       write(*,*) 'READATTRIBUTE Error return from NF_INQ_ATTNAME.',ret,nf_strerror(ret)
       irc=834
       return
    end if
    RET=NF_INQ_ATT(ncid,varid,a%att250,a%type,a%len)
    if (ret .ne. NF_NOERR) then
       write(*,*) 'READATTRIBUTE Error return from NF_INQ_ATT.',ret,nf_strerror(ret)
       irc=835
       return
    end if
    if (a%type.eq.nf_char)then
       allocate(a%ac(a%len),stat=irc)
       ret = nf_get_att_text (ncid,varid,a%att250,a%ac)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_TEXT.',ret
          irc=836
          return
       end if
    else if (a%type   .eq.nf_int1)then
       allocate(a%a1(a%len),stat=irc)
       ret = nf_get_att_int1 (ncid,varid,a%att250,a%a1)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_INT1.',ret
          irc=837
          return
       end if
    elseif (a%type   .eq.nf_int2)then
       allocate(a%a2(a%len),stat=irc)
       ret = nf_get_att_int2 (ncid,varid,a%att250,a%a2)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_INT2.',ret
          irc=838
          return
       end if
    elseif (a%type   .eq.nf_int)then
       allocate(a%a4(a%len),stat=irc)
       ret = nf_get_att_int (ncid,varid,a%att250,a%a4)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_INT.',ret
          irc=839
          return
       end if
    elseif (a%type   .eq.nf_real)then
       allocate(a%ar(a%len),stat=irc)
       ret = nf_get_att_real (ncid,varid,a%att250,a%ar)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_REAL.',ret
          irc=840
          return
       end if
    elseif (a%type   .eq.nf_double)then
       allocate(a%ad(a%len),stat=irc)
       ret = nf_get_att_double (ncid,varid,a%att250,a%ad)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_DOUBLE.',ret
          irc=841
          return
       end if
    else
    end if
    call chop0(a%att250,250)
    lena=length(a%att250,250,10)
    if (a%att250(1:lena).eq."history") then
    end if
    a%initialised=.true.
  end subroutine readAttribute

  function getAttribute(v,att,irc)
    implicit none
    type(attribute), pointer :: getAttribute
    type(variable) v
    character*(*) att
    integer lena
    type(attribute),pointer :: a
    integer irc
    nullify(getAttribute)
    a=>v%firstAttribute%next
    do while (.not.associated(a,target=v%lastAttribute))
       if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
          getAttribute => a
          a => v%lastAttribute
       else
          a => a%next
       end if
    end do
  end function getAttribute

  logical function hasAttribute(v,att,val,irc)
    implicit none
    type(variable) v
    character*(*) att, val
    integer lena
    type(attribute),pointer :: a
    integer irc
    integer lenb,length
    external length
    character*250 buff250
    hasAttribute=.false.
    a=>v%firstAttribute%next
    do while (.not.associated(a,target=v%lastAttribute))
       if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
          buff250=getAttributeText(a)
          lenb=length(buff250,250,1)
          if (val(1:len(val)).eq.buff250(1:lenb)) then
             hasAttribute=.true.
          end if
          a => v%lastAttribute
       else
          a => a%next
       end if
    end do
  end function hasAttribute

  character*250 function getAttributeText(a)
    implicit none
    type(attribute), pointer :: a
    integer ii, lens, lenb, length
    external length
    character*250 buff250,sum250
    lens=0
    sum250=""
    do ii=1,a%len
       if (a%type.eq.nf_int1) then
          write(buff250,'(I20)')a%a1(ii) 
       else if (a%type.eq.nf_int2) then
          write(buff250,'(I20)')a%a2(ii) 
       else if (a%type.eq.nf_int) then
          write(buff250,'(I20)')a%a4(ii) 
       else if (a%type.eq.nf_real) then
          write(buff250,'(F12.5)')a%ar(ii) 
       else if (a%type.eq.nf_double) then
          write(buff250,'(F17.10)')a%ad(ii) 
       else if (a%type.eq.nf_char) then
          buff250=a%ac(ii) 
       end if
       call chop(buff250,25)
       lenb=length(buff250,25,10)
       sum250=sum250(1:lens)//buff250(1:lenb)
       lens=min(250,lens+lenb)
    end do
    getAttributeText=sum250
  end function getAttributeText
  
  subroutine checkContents(f,par250,biok,irc)
    implicit none
    type(inventory) :: f
    character*250 par250,att250
    logical biok
    integer irc
!
    integer ii
    type(variable), pointer :: v => null()
    type(variable), pointer :: x => null()
    type(attribute), pointer :: a => null()
    integer lenp, length
    external length 
    integer gridatt, leng
    character*250 grid250
    biok=.true.
    ! identify lat/lon/time/parameter variables
    ! grid definition attribute
    grid250="grid_mapping"
    call chop(grid250,250)
    leng=length(grid250,250,10)
    lenp=length(par250,250,10)
    nullify(a)
    v=> f%firstVariable%next
    do while (.not.associated(v,target=f%lastVariable))
       if (.not.associated(a)) then
          a=>getAttribute(v,grid250(1:leng),irc)
          ! find grid variable
          if (associated(a)) then
             att250 = getAttributeText(a)
             call chop(att250,250)
             lena=length(att250,250,10)
             x=>f%firstVariable%next
             do while (.not.associated(x,target=f%lastVariable))
                if (att250(1:lena).eq.x%var250(1:x%lenv)) then
                   f%gridid=>x
                   x=>f%lastVariable
                else
                   x=>x%next
                end if
             end do
             if (.not.associated(f%gridid)) then
                write(*,*)myname,'Found INVALID grid attribute:',att250(1:lena)
                nullify(a) ! error in nc-file, search again...
             end if
          end if
       end if
       if ( hasAttribute(v,"standard_name","projection_x_coordinate",irc) .and. .not.associated(f%xcoid) ) then
          f%xcoid=>v
          v%fillAttribute => getAttribute(v,"_FillValue",irc)
       else if ( hasAttribute(v,"standard_name","projection_y_coordinate",irc) .and. .not.associated(f%ycoid) ) then
          f%ycoid=>v
          v%fillAttribute => getAttribute(v,"_FillValue",irc)
       else if ( (hasAttribute(v,"standard_name","latitude",irc).or. &
            &     hasAttribute(v,"long_name","latitude",irc)) .and. .not.associated(f%latid) ) then
          f%latid=>v
          v%fillAttribute => getAttribute(v,"_FillValue",irc)
       else if ( (hasAttribute(v,"standard_name","longitude",irc).or. &
            &     hasAttribute(v,"long_name","longitude",irc)) .and.  .not.associated(f%lonid)) then
          f%lonid=>v
          v%fillAttribute => getAttribute(v,"_FillValue",irc)
       else if (hasAttribute(v,"standard_name","time",irc) .and. .not.associated(f%timid)) then
          f%timid=>v
          v%fillAttribute => getAttribute(v,"_FillValue",irc)
       else if (hasAttribute(v,"standard_name","forecast_reference_time",irc) .and. .not.associated(f%terid)) then
          f%terid=>v
          if (v%len.ne.1) then
             write(*,*) 'checkContents Invalid reference_time length.',v%len
          end if
          v%fillAttribute => getAttribute(v,"_FillValue",irc)
       else if (v%var250(1:v%lenv) .eq. par250(1:lenp)   .and. .not.associated(f%parid)) then
          f%parid=>v
          v%fillAttribute => getAttribute(v,"_FillValue",irc)
       end if
       v=>v%next
    end do                 ! variable loop
    ! check if we have lat, lon and parameter...
    if (.not.associated(f%latid)) then
       write(*,*)'checkContents Missing LATITUDE in file: ', f%fn250(1:f%lenf)
       biok=.false.
    end if
    if (.not.associated(f%lonid)) then
       write(*,*)'checkContents Missing LONGITUDE in file: ', f%fn250(1:f%lenf)
       biok=.false.
    end if
    if (lenp.ne.0.and..not.associated(f%terid)) then
       write(*,*)'checkContents Missing REFERENCE_TIME in file: ', f%fn250(1:f%lenf)
       biok=.false.
    end if
    if (lenp.ne.0.and..not.associated(f%timid)) then
       write(*,*)'checkContents Missing TIME in file: ', f%fn250(1:f%lenf)
       ! biok=.false. ! we dont need time as long as we have reference time
    end if
    if (lenp.ne.0.and..not.associated(f%parid)) then
       write(*,*)'checkContents Missing PARAMETER "', par250(1:lenp),'" in file: ',f%fn250(1:f%lenf)
       biok=.false.
    end if
    if (leng.ne.0.and..not.associated(f%gridid)) then
       write(*,*)'checkContents Missing GRID "',grid250(1:leng),'" in file: ',f%fn250(1:f%lenf)
       biok=.false.
    end if
    ! check something else...
    ! if (biok) then
    !    ! loop over dimension, check if they equal lat or lon id.
    !    do dd=1,f%pardim
    !       do jj=1,f%latdim
    !          if (f%latdims(jj).eq.f%pardims(dd)) then
    !          end if
    !       end do
    !       do jj=1,f%londim
    !          if (f%londims(jj).eq.f%pardims(dd)) then
    !          end if
    !       end do
    !       do jj=1,f%timdim
    !          if (f%timdims(jj).eq.f%pardims(dd)) then
    !          end if
    !       end do
    !    end do
    ! end if
    return
  end subroutine checkContents
  
  subroutine readRealData(v,biok,irc)
    implicit none
    type(variable),pointer :: v
    logical biok
    integer irc
    if (associated(v)) then
       call readData(v,biok,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from readData.',irc
          return
       end if
       call makeRealData(v,biok,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from makeRealData.',irc
          return
       end if
       ! write(*,*)'readRealData sample:',v%fd(1:min(v%lend,5))
    else
       write(*,*)'System error: attempt to read undefined variable.'
       biok=.false.
    end if
  end subroutine readRealData


  subroutine readData(v,biok,irc)
    implicit none
    type(variable),pointer :: v
    logical biok
    integer irc
    character*1  fillc
    integer*1    fill1
    integer*2    fill2
    integer*4    fill4
    real*4       fillr
    real*8       filld
    ! init
    fillc=char(nf_fill_char)
    fill1=nf_fill_int1
    fill2=nf_fill_int2
    fill4=nf_fill_int
    fillr=nf_fill_real
    filld=nf_fill_double
    ! allocate real data array
    ! write(*,*)'readData ',gettype(v%type),"  Reading data: ",v%var250(1:v%lenv)
    if (.not.biok) return
    if (v%type.eq.nf_int1) then
       if (v%len1.ne.v%len) then
          if (associated(v%f1)) deallocate(v%f1,stat=irc)
          v%len1=v%len
          allocate(v%f1(v%len1),stat=irc)
          if(irc.ne.0) then
             write(*,*)'readData Unable to allocate work array:',v%len
             irc=940
             return
          end if
       end if
       ret = nf_get_vara_INT1(v%ncid,v%varid,v%sta,v%lim,v%f1)
       if (ret .ne. NF_NOERR) then
          write(*,*)'readData Unable to load array from file.'
          irc=941
          return
       end if
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%a1)) then
             fill1=v%fillAttribute%a1(1)
             do ii=1,v%len1
                if (v%f1(ii).eq.fill1) v%f1(ii)=nf_fill_int1
             end do
          end if
       end if
    else if (v%type.eq.nf_int2) then
       if (v%len2.ne.v%len) then
          v%len2=v%len
          if (associated(v%f2)) deallocate(v%f2,stat=irc)
          allocate(v%f2(v%len2),stat=irc)
          if(irc.ne.0) then
             write(*,*)'readData Unable to allocate work array:',v%len
             irc=942
             return
          end if
       end if
       ret = nf_get_vara_INT2(v%ncid,v%varid,v%sta,v%lim,v%f2)
       if (ret .ne. NF_NOERR) then
          write(*,*)'readData Unable to load array from file.'
          irc=943
          return
       end if
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%a2)) then
             fill2=v%fillAttribute%a2(1)
             do ii=1,v%len2
                if (v%f2(ii).eq.fill2) v%f2(ii)=nf_fill_int2
             end do
          end if
       end if
    else if (v%type.eq.nf_int) then
       if (v%len4.ne.v%len) then
          v%len4=v%len
          if (associated(v%f4)) deallocate(v%f4,stat=irc)
          allocate(v%f4(v%len4),stat=irc)
          if(irc.ne.0) then
             write(*,*)'readData Unable to allocate work array:',v%len
             irc=944
             return
          end if
       end if
       !       write(*,*)'READFIELD: ',v%nrdim,v%len4,v%len
       !       if (v%nrdim.eq.0) then
       !          ret = nf_get_var_INT(v%ncid,v%varid,v%f4)
       !       else
       ret = nf_get_vara_INT(v%ncid,v%varid,v%sta,v%lim,v%f4)
       !       end if
       if (ret .ne. NF_NOERR) then
          write(*,*)'readData Unable to load array from file.'
          irc=945
          return
       end if
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%a4)) then
             fill4=v%fillAttribute%a4(1)
             do ii=1,v%len4
                if (v%f4(ii).eq.fill4) v%f4(ii)=nf_fill_int
             end do
          end if
       end if
    else if (v%type.eq.nf_real) then
       if (v%lenr.ne.v%len) then
          v%lenr=v%len
          if (associated(v%fr)) deallocate(v%fr,stat=irc)
          allocate(v%fr(v%lenr),stat=irc)
          if(irc.ne.0) then
             write(*,*)'readData Unable to allocate work array:',v%len
             irc=946
             return
          end if
       end if
       ret = nf_get_vara_real(v%ncid,v%varid,v%sta,v%lim,v%fr)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'Variable=',v%var250(1:v%lenv),'   ',v%f%fn250(1:v%f%lenf),v%lenr
          write(*,*) 'Start, lim=',v%sta,v%lim
          write(*,*)'readData Unable to load array from file. ',nf_strerror(ret)
          irc=947
          return
       end if
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%ar)) then
             fillr=v%fillAttribute%ar(1)
             do ii=1,v%lenr
                if (v%fr(ii).eq.fillr) v%fr(ii)=nf_fill_real
             end do
          end if
       end if
    else if (v%type.eq.nf_double) then
       if (v%lend.ne.v%len) then
          v%lend=v%len
          if (associated(v%fd)) deallocate(v%fd,stat=irc)
          allocate(v%fd(v%lend),stat=irc)
          if(irc.ne.0) then
             write(*,*)'readData Unable to allocate work array:',v%len
             irc=948
             return
          end if
       end if
       ret = nf_get_vara_double(v%ncid,v%varid,v%sta,v%lim,v%fd)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'Variable=',v%var250(1:v%lenv),'   ',v%f%fn250(1:v%f%lenf),v%lend
          write(*,*) 'Start, lim=',v%sta,v%lim
          write(*,*)'readData Unable to load array from file.',nf_strerror(ret)
          irc=949
          return
       end if
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%ad)) then
             filld=v%fillAttribute%ad(1)
             do ii=1,v%lend
                if (v%fd(ii).eq.filld) v%fd(ii)=nf_fill_double
             end do
          end if
       end if
    end if
  end subroutine readData

  subroutine makeRealData(v,biok,irc)
    implicit none
    type(variable) :: v
    logical biok
    integer irc
    integer ii
    ! allocate real data array
    ! write(*,*)'makeRealData ',gettype(v%type),"  making real data: ",v%var250(1:v%lenv)
    if (.not.biok) return
    if (v%type.ne.nf_double) then
       if (v%lend.ne.v%len) then
          v%lend=v%len
          if (associated(v%fd)) deallocate(v%fd,stat=irc)
          allocate(v%fd(v%lend),stat=irc)
          if(irc.ne.0) then
             write(*,*)'makeRealData Unable to allocate work array:',v%len
             irc=945
             return
          end if
       end if
    end if
    if (v%type.eq.nf_int1) then
       do ii=1,v%len1
          if (v%f1(ii) .eq. fill1) then
             v%fd(ii)=filld
          else
             v%fd(ii)=v%f1(ii)
          end if
       end do
       deallocate(v%f1,stat=irc)
       if(irc.ne.0) then
          write(*,*)'makeRealData Unable to deallocate work array.'
          irc=945
          return
       end if
    else if (v%type.eq.nf_int2) then
       do ii=1,v%len2
          if (v%f2(ii) .eq. fill2) then
             v%fd(ii)=filld
          else
             v%fd(ii)=v%f2(ii)
          end if
       end do
       deallocate(v%f2,stat=irc)
       if(irc.ne.0) then
          write(*,*)'makeRealData Unable to deallocate work array.'
          irc=945
          return
       end if
    else if (v%type.eq.nf_int) then
       do ii=1,v%len4
          if (v%f4(ii) .eq. fill4) then
             v%fd(ii)=filld
          else
             v%fd(ii)=v%f4(ii)
          end if
       end do
       deallocate(v%f4,stat=irc)
       if(irc.ne.0) then
          write(*,*)'makeRealData Unable to deallocate work array.'
          irc=945
          return
       end if
    else if (v%type.eq.nf_real) then
       do ii=1,v%lenr
          if (v%fr(ii) .eq. fillr) then
             v%fd(ii)=filld
          else
             v%fd(ii)=v%fr(ii)
          end if
       end do
       deallocate(v%fr,stat=irc)
       if(irc.ne.0) then
          write(*,*)'makeRealData Unable to deallocate work array.'
          irc=945
          return
       end if
    else if (v%type.eq.nf_double) then
    end if
  end subroutine makeRealData

  subroutine unmakeRealData(v,irc)
    implicit none
    type(variable) :: v
    integer irc
    integer ii
 ! allocate real data array
    if (associated(v%fd)) then
       if (v%type.eq.nf_int1) then
          do ii=1,v%len1
             allocate(v%f1(v%len1),stat=irc)
             if(irc.ne.0) then
                write(*,*)'unmakeRealData Unable to allocate work array:',v%len1
                irc=945
                return
             end if
             if (v%fd(ii) .eq. filld) then
                v%f1(ii)=fill1
             else
                v%f1(ii)=v%fd(ii)
             end if
          end do
          deallocate(v%fd,stat=irc)
          if(irc.ne.0) then
             write(*,*)'unmakeRealData Unable to deallocate work array.'
             irc=945
             return
          end if
          v%lend=0
       else if (v%type.eq.nf_int2) then
          allocate(v%f2(v%len2),stat=irc)
          if(irc.ne.0) then
             write(*,*)'unmakeRealData Unable to allocate work array:',v%len2
             irc=945
             return
          end if
          do ii=1,v%len2
             if (v%fd(ii) .eq. filld) then
                v%f2(ii)=fill2
             else
                v%f2(ii)=v%fd(ii)
             end if
          end do
          deallocate(v%fd,stat=irc)
          if(irc.ne.0) then
             write(*,*)'unmakeRealData Unable to deallocate work array.'
             irc=945
             return
          end if
          v%lend=0
       else if (v%type.eq.nf_int) then
          allocate(v%f4(v%len4),stat=irc)
          if(irc.ne.0) then
             write(*,*)'unmakeRealData Unable to allocate work array:',v%len4
             irc=945
             return
          end if
          do ii=1,v%len4
             if (v%fd(ii) .eq. filld) then
                v%f4(ii)=fill4
             else
                v%f4(ii)=v%fd(ii)
             end if
          end do
          deallocate(v%fd,stat=irc)
          if(irc.ne.0) then
             write(*,*)'unmakeRealData Unable to deallocate work array.'
             irc=945
             return
          end if
          v%lend=0
       else if (v%type.eq.nf_real) then
          allocate(v%fr(v%lenr),stat=irc)
          if(irc.ne.0) then
             write(*,*)'unmakeRealData Unable to allocate work array:',v%lenr
             irc=945
             return
          end if
          do ii=1,v%lenr
             if (v%fd(ii) .eq. filld) then
                v%fr(ii)=fillr
             else
                v%fr(ii)=v%fd(ii)
             end if
          end do
          deallocate(v%fd,stat=irc)
          if(irc.ne.0) then
             write(*,*)'unmakeRealData Unable to deallocate work array.'
             irc=945
             return
          end if
          v%lend=0
       else if (v%type.eq.nf_double) then
       end if
    else
       write(*,*)'UNMAKEREALDATA No real data was found.'
       irc=944
       return
    end if
    
  end subroutine unmakeRealData
  
  subroutine clearData(f,irc)
    implicit none
    type(inventory) :: f
    integer irc
    type(variable), pointer  :: v
    integer ii
    ! loop over variables, remove data (keep attributes etc)
    if (f%initialised) then
       v=>f%firstVariable%next
       do while(.not.associated(v,target=f%lastVariable))
          call clearField(v,irc)
          v=>v%next
       end do
    end if
  end subroutine clearData

  subroutine clearField(v,irc)
    implicit none
    type(variable) :: v
    integer irc
    if (associated(v%f1)) deallocate(v%f1,stat=irc)
    if (associated(v%f2)) deallocate(v%f2,stat=irc)
    if (associated(v%f4)) deallocate(v%f4,stat=irc)
    if (associated(v%fr)) deallocate(v%fr,stat=irc)
    if (associated(v%fd)) deallocate(v%fd,stat=irc)
    v%len1=0
    v%len2=0
    v%len4=0
    v%lenr=0
    v%lend=0
    return
  end subroutine clearField

  subroutine removeField(v)
    implicit none
    type(variable),pointer :: v
    v%prev%next => v%next
    v%next%prev => v%prev
    return
  end subroutine removeField


  subroutine closeFile(f,irc)
    implicit none
    type(inventory) :: f
    integer irc
    integer ret
    integer lenb, length
    character*250 buff250
    external length
    if (f%opened) then
       f%opened=.false.
       ! delete all field data...
       call clearData(f,irc)
!       write(buff250,"(A)") "  Closing: "//f%fn250(1:f%lenf)
!       call chop0(buff250,250)
!       lenb=length(buff250,250,10)
!       write(*,*)'closeFile '//buff250(1:lenb)
       ret=NF_CLOSE(f%ncid)        ! end definitions: leave define mode
       if (ret .ne. NF_NOERR) then
          write(*,*) "closeFile ERROR return from NF_CLOSE:",nf_strerror(ret),f%ncid,nf_strerror(ret)
          irc=790
          return
       end if
    end if
    return
  end subroutine closeFile

  type(inventory) function copyInventory(f,irc)
    implicit none
    type(inventory) f
    integer irc
    integer ii
    type(variable),pointer :: v,fv
    type(attribute),pointer :: a,fa
    copyInventory%ncid=f%ncid
    copyInventory%fn250=f%fn250
    copyInventory%lenf=f%lenf
    copyInventory%nrdim=f%nrdim
    copyInventory%nrvar=f%nrvar
    copyInventory%nrgatt=f%nrgatt
    allocate(copyInventory%dim250(copyInventory%nrdim),stat=irc)
    allocate(copyInventory%lend(copyInventory%nrdim),stat=irc)
    allocate(copyInventory%pos,stat=irc)
    allocate(copyInventory%pos%pos(copyInventory%nrdim),stat=irc)
    allocate(copyInventory%pos%sta(copyInventory%nrdim),stat=irc)
    allocate(copyInventory%pos%lim(copyInventory%nrdim),stat=irc)
    do ii=1,copyInventory%nrdim
       copyInventory%dim250(ii)=f%dim250(ii)
       copyInventory%lend(ii)=f%lend(ii)
    end do
    copyInventory%pos%nrdim=f%pos%nrdim
    do ii=1,copyInventory%pos%nrdim
       copyInventory%pos%pos(ii)=f%pos%pos(ii)
       copyInventory%pos%sta(ii)=f%pos%sta(ii)
       copyInventory%pos%lim(ii)=f%pos%lim(ii)
    end do
    copyInventory%pos%loc=f%pos%loc
    copyInventory%unlimdimid=f%unlimdimid
    fv=f%firstVariable%next
    do while (.not.associated(fv,target=f%lastVariable))
       v=>copyVariable(fv,irc)
       if (associated(fv,target=f%latid)) then
          copyInventory%latid=>v
       else if (associated(fv,target=f%lonid)) then
          copyInventory%lonid=>v
       else if (associated(fv,target=f%timid)) then
          copyInventory%timid=>v
       else if (associated(fv,target=f%terid)) then
          copyInventory%terid=>v
       else if (associated(fv,target=f%parid)) then
          copyInventory%parid=>v
       else if (associated(fv,target=f%xcoid)) then
          copyInventory%xcoid=>v
       else if (associated(fv,target=f%ycoid)) then
          copyInventory%ycoid=>v
       else if (associated(fv,target=f%gridid)) then
          copyInventory%gridid=>v
       end if
       v%prev=>copyInventory%lastVariable%prev
       v%next=>copyInventory%lastVariable
       copyInventory%lastVariable%prev%next=>v
       copyInventory%lastVariable%prev=>v
       nullify(v)
       fv=>fv%next
    end do
    fa=f%firstAttribute%next
    do while (.not.associated(fa,target=f%lastAttribute))
       a=copyAttribute(fa,irc)
       a%prev=>f%lastAttribute%prev
       a%next=>f%lastAttribute
       f%lastAttribute%prev%next=>a
       f%lastAttribute%prev=>a
       nullify(a)
       fa=>fa%next
    end do
  end function copyInventory

  subroutine copyGrid(o,f,irc)
    type(inventory), target :: o
    type(inventory) f
    integer irc
    integer ii,jj,nrdim
    type(variable), pointer :: sv,v
    type(attribute), pointer :: sa,a
    integer, allocatable :: indx(:)
    logical used
    o%ncid=f%ncid
    o%fn250=f%fn250
    o%lenf=f%lenf
    !
    ! copy dimensions (can be compressed later)
    o%nrdim=f%nrdim
    allocate(o%dim250(o%nrdim),stat=irc)
    allocate(o%lend(o%nrdim),stat=irc)
    do ii=1,o%nrdim
       o%dim250(ii)=f%dim250(ii)
       o%lend(ii)=f%lend(ii)
    end do
    allocate(o%pos,stat=irc)
    o%pos%nrdim=f%pos%nrdim
    allocate(o%pos%pos(o%pos%nrdim),stat=irc)
    allocate(o%pos%sta(o%pos%nrdim),stat=irc)
    allocate(o%pos%lim(o%pos%nrdim),stat=irc)
    do ii=1,o%pos%nrdim
       o%pos%pos(ii)=f%pos%pos(ii)
       o%pos%sta(ii)=f%pos%sta(ii)
       o%pos%lim(ii)=f%pos%lim(ii)
    end do
    o%unlimdimid=f%unlimdimid 
    o%pos%loc=f%pos%loc
    !
    allocate(o%firstVariable,o%lastVariable, stat=irc)
    o%firstVariable%next=>o%lastVariable
    o%lastVariable%prev=>o%firstVariable
    !
    o%nrvar=0
    sv=> f%firstVariable%next
    do while (.not.associated(sv,target=f%lastVariable))
       used=.true.
       if (associated(sv,target=f%latid)) then
          o%nrvar=o%nrvar+1
          v=>copyVariable(sv,irc)
          o%latid => v
       else if (associated(sv,target=f%lonid)) then
          o%nrvar=o%nrvar+1
          v=>copyVariable(sv,irc)
          o%lonid => v
       else if (associated(sv,target=f%xcoid)) then
          o%nrvar=o%nrvar+1
          v=>copyVariable(sv,irc)
          o%xcoid => v
       else if (associated(sv,target=f%ycoid)) then
          o%nrvar=o%nrvar+1
          v=>copyVariable(sv,irc)
          o%ycoid => v
       else if (associated(sv,target=f%gridid)) then
          o%nrvar=o%nrvar+1
          v=>copyVariable(sv,irc)
          o%gridid => v
       else
          used=.false.
       end if
       if (used) then
          ! write(*,*) myname,'Copying variable Field:',v%var250(1:v%lenv)
          call copyVariableField(v,sv,irc)
          v%prev=>o%lastVariable%prev
          v%next=>o%lastVariable
          o%lastVariable%prev%next=>v
          o%lastVariable%prev=>v
          v%f=>o
          nullify(v)
       end if
       sv=>sv%next
    end do
    o%nrgatt=f%nrgatt
    !
    allocate(o%firstAttribute,o%lastAttribute, stat=irc)
    o%firstAttribute%next=>o%lastAttribute
    o%lastAttribute%prev=>o%firstAttribute
    !
    sa=>f%firstAttribute%next
    do while (.not.associated(sa,f%lastAttribute))
       a=>copyAttribute(sa,irc)
       a%prev=>o%lastAttribute%prev
       a%next=>o%lastAttribute
       o%lastAttribute%prev%next=>a
       o%lastAttribute%prev=>a
       nullify(a)
       sa=>sa%next
    end do
    o%initialised=f%initialised
  end subroutine copyGrid

  subroutine compressDimensions(f)
    ! remove dimensions that are not used...
    implicit none
    type(inventory) :: f
    type(variable), pointer :: v
    integer ii, nrdim, unlimdimid, indx(f%nrdim), cnt(f%nrdim)
    do jj=1,f%nrdim
       cnt(jj)=0
       indx(jj)=0
    end do
    v=>f%firstVariable%next
    do while (.not.associated(v,f%lastVariable))
       do jj=1,v%nrdim
          cnt(v%ind(jj))=cnt(v%ind(jj))+1
       end do
       v=>v%next
    end do
    unlimdimid=0
    nrdim=0
    do jj=1,f%nrdim
       if (cnt(jj).ne.0) then
          nrdim=nrdim+1
          indx(jj)=nrdim ! indx(jj) <= jj
          if (jj.eq.f%unlimdimid) then
             unlimdimid=indx(jj)
          end if
       end if
    end do
    v=>f%firstVariable%next
    do while (.not.associated(v,f%lastVariable))
       do jj=1,v%nrdim
          v%ind(jj)=indx(v%ind(jj))
       end do
       v=>v%next
    end do
    f%unlimdimid=unlimdimid
    do jj=1,f%nrdim
       if (indx(jj).ne.0) then
          f%dim250(indx(jj))=f%dim250(jj) ! indx(jj) <= jj
          f%lend(indx(jj))=f%lend(jj) ! indx(jj) <= jj
          f%pos%pos(indx(jj))=f%pos%pos(jj) ! indx(jj) <= jj
          f%pos%sta(indx(jj))=f%pos%sta(jj) ! indx(jj) <= jj
          f%pos%lim(indx(jj))=f%pos%lim(jj) ! indx(jj) <= jj
       end if
    end do
    f%pos%loc=0
    f%pos%nrdim=nrdim
    f%nrdim=nrdim
  end subroutine compressDimensions

  subroutine copyVariableField(v,s,irc)
    implicit none
    type(variable), pointer :: s,v
    integer irc
    integer jj
    ! copy field
    if (s%initialised) then
       v%len=s%len
       v%len1=s%len1
       v%len2=s%len2
       v%len4=s%len4
       v%lenr=s%lenr
       v%lend=s%lend
       v%lenc=s%lenc
       v%type=s%type
       if (v%len1.gt.0) then
          if (associated(v%f1)) deallocate(v%f1)
          allocate(v%f1(v%len1),stat=irc)
          if (associated(s%f1)) then
             do jj=1,v%len1
                v%f1(jj)=s%f1(jj)
             end do
          end if
       end if
       if (v%len2.gt.0) then
          if (associated(v%f2)) deallocate(v%f2)
          allocate(v%f2(v%len2),stat=irc)
          if (associated(s%f2)) then
             do jj=1,v%len2
                v%f2(jj)=s%f2(jj)
             end do
          end if
       end if
       if (v%len4.gt.0) then
          if (associated(v%f4)) deallocate(v%f4)
          allocate(v%f4(v%len4),stat=irc)
          if (associated(s%f4)) then
             do jj=1,v%len4
                v%f4(jj)=s%f4(jj)
             end do
          end if
       end if
       if (v%lenr.gt.0) then
          if (associated(v%fr)) deallocate(v%fr)
          allocate(v%fr(v%lenr),stat=irc)
          if (associated(s%fr)) then
             do jj=1,v%lenr
                v%fr(jj)=s%fr(jj)
             end do
          end if
       end if
       if (v%lend.gt.0) then
          if (associated(v%fd)) deallocate(v%fd)
          allocate(v%fd(v%lend),stat=irc)
          if (associated(s%fd)) then
             do jj=1,v%lend
                v%fd(jj)=s%fd(jj)
             end do
          end if
       end if
       v%initialised=.true.
    end if
  end subroutine copyVariableField

  
  function copyVariable(v,irc)
    type(variable), pointer :: copyVariable
    type(variable) v
    integer irc
    integer ii
    type(attribute),pointer :: a,va
    allocate(copyVariable,stat=irc)
    copyVariable%ncid=v%ncid
    copyVariable%varid=v%varid
    copyVariable%var250=v%var250
    copyVariable%lenv=v%lenv
    copyVariable%type=v%type
    copyVariable%nrdim=v%nrdim
    copyVariable%nratt=v%nratt
    copyVariable%scale=v%scale
    allocate(copyVariable%ind(copyVariable%nrdim),stat=irc)
    allocate(copyVariable%sta(copyVariable%nrdim),stat=irc)
    allocate(copyVariable%lim(copyVariable%nrdim),stat=irc)
    do ii=1,v%nrdim
       copyVariable%ind(ii)=v%ind(ii)
       copyVariable%sta(ii)=v%sta(ii)
       copyVariable%lim(ii)=v%lim(ii)
    end do
    copyVariable%len=0
    copyVariable%lenc=0
    copyVariable%len1=0
    copyVariable%len2=0
    copyVariable%len4=0
    copyVariable%lenr=0
    copyVariable%lend=0
    allocate(copyVariable%firstAttribute,copyVariable%lastAttribute,copyVariable%fillAttribute, stat=irc)
    copyVariable%firstAttribute%next=>copyVariable%lastAttribute
    copyVariable%lastAttribute%prev=>copyVariable%firstAttribute
    va=>v%firstAttribute%next
    do while (.not.associated(va,v%lastAttribute))
       a=>copyAttribute(va,irc)
       a%prev=>copyVariable%lastAttribute%prev
       a%next=>copyVariable%lastAttribute
       copyVariable%lastAttribute%prev%next=>a
       copyVariable%lastAttribute%prev=>a
       a%v => copyVariable
       nullify(a)
       va=>va%next
    end do
  end function copyVariable

  function copyAttribute(a,irc)
    type(attribute), pointer :: copyAttribute
    type(attribute) a
    integer irc
    integer ii
    allocate(copyAttribute,stat=irc)
    copyAttribute%ncid=a%ncid
    copyAttribute%varid=a%varid
    copyAttribute%attid=a%attid
    copyAttribute%att250=a%att250
    copyAttribute%lena=a%lena
    copyAttribute%type=a%type
    copyAttribute%len=a%len
    if (a%len.gt.0) then
       if (copyAttribute%type.eq.nf_char) then
          allocate(copyAttribute%ac(a%len),stat=irc)
          do ii=1,a%len
             copyAttribute%ac(ii)=a%ac(ii)
          end do
       else if (copyAttribute%type.eq.nf_int1) then
          allocate(copyAttribute%a1(a%len),stat=irc)
          do ii=1,a%len
             copyAttribute%a1(ii)=a%a1(ii)
          end do
       else if (copyAttribute%type.eq.nf_int2) then
          allocate(copyAttribute%a2(a%len),stat=irc)
          do ii=1,a%len
             copyAttribute%a2(ii)=a%a2(ii)
          end do
       else if (copyAttribute%type.eq.nf_int) then
          allocate(copyAttribute%a4(a%len),stat=irc)
          do ii=1,a%len
             copyAttribute%a4(ii)=a%a4(ii)
          end do
       else if (copyAttribute%type.eq.nf_real) then
          allocate(copyAttribute%ar(a%len),stat=irc)
          do ii=1,a%len
             copyAttribute%ar(ii)=a%ar(ii)
          end do
       else if (copyAttribute%type.eq.nf_double) then
          allocate(copyAttribute%ad(a%len),stat=irc)
          do ii=1,a%len
             copyAttribute%ad(ii)=a%ad(ii)
          end do
       end if
    end if
  end function copyAttribute

  function insertField(f,nam,i,irc)
    implicit none
    type(variable), pointer ::  insertField
    type(inventory),target :: f
    character*(*) nam
    character*250 nam250
    type(index) :: i
    integer irc
    !
    type(variable), pointer :: v
    type(attribute),pointer :: a,va
    integer length, lenn, kk
    type(dimension), pointer :: d
    external length
    !
    nam250=nam(1:len(nam))
    call chop(nam250,250)
    lenn=length(nam250,250,10)
    allocate(insertField,stat=irc)
    insertField%f=>f
    insertField%ncid=f%ncid
    insertField%varid=0
    insertField%var250=nam250
    insertField%lenv=lenn
    insertField%type=nf_double
    insertField%nrdim=i%nrdim
    insertField%nratt=0
    insertField%scale=1.0D0
    allocate(insertField%ind(insertField%nrdim),stat=irc)
    if (irc.ne.0) then
       write(*,*)'insertField Unable to allocate:',insertField%nrdim
       return
    end if
    allocate(insertField%lim(insertField%nrdim),stat=irc)
    if (irc.ne.0) then
       write(*,*)'insertField Unable to allocate:',insertField%nrdim
       return
    end if
    allocate(insertField%sta(insertField%nrdim),stat=irc)
    if (irc.ne.0) then
       write(*,*)'insertField Unable to allocate:',insertField%nrdim
       return
    end if
    insertField%len=1
    kk=0
    d=>i%firstDimension%next
    do while (.not.associated(d,target=i%lastDimension))
       if (d%ind.lt.1.or.d%ind.gt.f%nrdim) then
          write(*,*)'insertField Invalid dimension index:',kk,d%ind,f%nrdim
          irc=834
          return
       end if
       kk=kk+1
       insertField%ind(kk)=d%ind
       insertField%sta(kk)=1
       insertField%lim(kk)=d%lim
       insertField%len=insertField%len*insertField%lim(kk)
       d=>d%next
    end do
    if (kk.ne.insertField%nrdim) then
       write(*,*)'insertField Dimension mismatch:',kk,insertField%nrdim
       d=>i%firstDimension%next
       do while (.not.associated(d,target=i%lastDimension))
          write(*,*) 'insertField Dimension:',d%ind, d%lim
          d=>d%next
       end do
       irc=945
       return
    end if
    insertField%lenc=0
    insertField%len1=0
    insertField%len2=0
    insertField%len4=0
    insertField%lenr=0
    insertField%lend=insertField%len
    allocate(insertField%fd(insertField%lend), stat=irc)
    if (irc.ne.0) then
       write(*,*)'insertField Unable to allocate:',insertField%lend
       return
    end if
    do kk=1,insertField%lend
       insertField%fd(kk)=0.0D0
    end do
    allocate(insertField%firstAttribute,insertField%lastAttribute,insertField%fillAttribute, stat=irc)
    insertField%firstAttribute%next=>insertField%lastAttribute
    insertField%lastAttribute%prev=>insertField%firstAttribute
    ! insert into inventory
    insertField%prev=>f%lastVariable%prev
    insertField%next=>f%lastVariable
    f%lastVariable%prev%next=>insertField
    f%lastVariable%prev=>insertField
  end function insertField

  function getDimIndex(v,irc)
    type(index), pointer :: getDimIndex
    type(variable),pointer :: v
    integer irc
    integer kk
    type(dimension), pointer :: d
    allocate(getDimIndex,stat=irc)
    allocate(getDimIndex%firstDimension,getDimIndex%lastDimension, stat=irc)
    getDimIndex%firstDimension%next => getDimIndex%lastDimension
    getDimIndex%lastDimension%prev  => getDimIndex%firstDimension
    if (associated(v)) then
       getDimIndex%nrdim=v%nrdim
       do kk=1,v%nrdim
          allocate(d,stat=irc)
          d%ind=v%ind(kk)
          d%sta=v%sta(kk)
          d%lim=v%lim(kk)
          d%prev=>getDimIndex%lastDimension%prev
          d%next=>getDimIndex%lastDimension
          getDimIndex%lastDimension%prev%next=>d
          getDimIndex%lastDimension%prev=>d
          nullify(d)
       end do
    else
       getDimIndex%nrdim=0
    end if
  end function getDimIndex

  function getAllIndex(f,irc)
    type(index), pointer :: getAllIndex
    type(inventory) :: f
    integer irc
    integer kk
    type(dimension), pointer :: d
    allocate(getAllIndex,stat=irc)
    getAllIndex%nrdim=f%nrdim
    allocate(getAllIndex%firstDimension,getAllIndex%lastDimension, stat=irc)
    getAllIndex%firstDimension%next => getAllIndex%lastDimension
    getAllIndex%lastDimension%prev  => getAllIndex%firstDimension
    do kk=1,f%nrdim
       allocate(d,stat=irc)
       d%ind=kk
       d%sta=1
       d%lim=f%pos%lim(kk)
       d%prev=>getAllIndex%lastDimension%prev
       d%next=>getAllIndex%lastDimension
       getAllIndex%lastDimension%prev%next=>d
       getAllIndex%lastDimension%prev=>d
       nullify(d)
    end do
  end function getAllIndex

  subroutine setValue(v,val,irc)
    implicit none
    type(variable), pointer :: v
    real val
    integer irc
    integer loc
    type(index),pointer :: i
    loc=getLocation(v,irc)
    if (irc.ne.0) then
       write(*,*) 'addValue Error return from getLocation.',irc
       return
    end if
    ! write(*,*)'addValue location:',loc
    if (loc.lt.1.or.loc.gt.v%lend) then
       write(*,*)myname,'Invalid position.',loc,v%lend
       i=>getdimindex(v,irc)
       call printindex(i)
       call printpos(v%f%pos)
    end if
    v%fd(loc)=val
  end subroutine setValue
  
  subroutine addValue(v,val,irc)
    implicit none
    type(variable), pointer :: v
    real val
    integer irc
    integer loc
    type(index),pointer :: i
    loc=getLocation(v,irc)
    if (irc.ne.0) then
       write(*,*) 'addValue Error return from getLocation.',irc
       return
    end if
    ! write(*,*)'addValue location:',loc
    if (loc.lt.1.or.loc.gt.v%lend) then
       write(*,*)myname,'Invalid position.',loc,v%lend
       i=>getdimindex(v,irc)
       call printindex(i)
       call printpos(v%f%pos)
    end if
    v%fd(loc)=v%fd(loc)+val
  end subroutine addValue
  
  real function valuePosition(v,irc)
    implicit none
    type(variable) :: v
    type(index) :: i
    integer irc
    integer loc
    loc=getLocation(v,irc)
    if (irc.ne.0) then
       write(*,*) 'valuePosition Error return from getLocation.',irc
       return
    end if
    !write(*,*) 'ValuePosition loc:',loc,v%fd(loc)
    valuePosition=v%fd(loc)
  end function valuePosition
  
  real function valueWeighted(v,w,irc)
    implicit none
    type(variable) :: v
    type(weight) :: w
    integer irc
    integer kk
    integer loc
    logical fill
    real val
    fill=.false.
    valueWeighted=0.0D0
    do kk=1,w%nweight
       call addPos(v%f%pos,w%pos(kk),irc)
       loc=getLocation(v,irc)
       if (irc.ne.0) then
          write(*,*) 'valueWeighted Error return from getLocation.',irc
          return
       end if
       ! call printPos(v%f%pos)
       ! write(*,*)'ValueWeighted LOC:',loc, v%len, v%lend, v%lenr, v%type, nf_real, nf_double
       val=v%fd(loc)
       if (val.eq.nf_fill_double) then
          fill=.true.
       else
          if (abs(val).gt.1000.0D0) then
             write(*,*) myname,'Invalid exp:',val
          end if
          valueWeighted=valueWeighted + w%w(kk) * val
       end if
    end do
    if (fill) then
       valueWeighted=nf_fill_double
    end if
  end function valueWeighted

  integer function getLocation(v,irc)
    type(variable) :: v
    integer irc
    integer kk
    if (.not.associated(v%f)) then
       write(*,*)'getLocation Invalid Inventory!'
       irc=991
       return
    end if
!    if (v%f%pos%loc.ne.0) then
!       getLocation=v%f%pos%loc
!    else
       if (v%nrdim .gt. 0) then
          getLocation=v%f%pos%pos(v%ind(v%nrdim))-1
          ! write(*,*) 'getLocation A:',getLocation
          ! write(*,*) 'getLocation B:',v%nrdim
          ! write(*,*) 'getLocation C:',v%ind(v%nrdim)
          ! call printPos(v%f%pos)
          do kk=v%nrdim-1,1,-1
             getLocation=getLocation*v%lim(kk)+v%f%pos%pos(v%ind(kk))-1
          end do
          getLocation=getLocation+1
       else
          getLocation=1
       end if
       v%f%pos%loc=getLocation
!    end if
    ! write(*,*) 'getLocation X:',getLocation
    return
  end function getLocation
  
  subroutine setPos(f,i,value)
    type(inventory) :: f
    type(index)  :: i
    integer :: value
    type(dimension), pointer :: d
    d=>i%firstDimension%next
    do while (.not.associated(d,target=i%lastDimension))
       if (d%ind.lt.1.or.d%ind.gt.f%nrdim) then
          write(*,*) 'setPos Ivalid index:',jj,d%ind,f%nrdim
       else
          f%pos%pos(d%ind)=value
       end if
       d=>d%next
    end do
    f%pos%loc=0
  end subroutine setPos

  subroutine addPos(p,t,irc)
    type(position) :: p
    type(position) :: t
    integer irc
    integer kk
    do kk=1,p%nrdim
       p%pos(kk)=p%pos(kk)+t%pos(kk)
    end do
    p%loc=t%loc
  end subroutine addPos
  
  subroutine copyPos(p,t,irc)
    type(position), pointer :: p
    type(position) :: t
    integer irc
    integer kk
    type(position), pointer :: b
    if (t%nrdim .gt. p%nrdim) then
       allocate(b,stat=irc)
       call allocatePos(t%nrdim,b,irc)
       p=>b
    end if
    p%nrdim=t%nrdim
    do kk=1,p%nrdim
       p%pos(kk)=t%pos(kk)
       p%sta(kk)=t%sta(kk)
       p%lim(kk)=t%lim(kk)
    end do
    p%loc=t%loc
  end subroutine copyPos
    
  function copyIndex(i)
    ! copies index "i" to "t"
    type(index),pointer :: copyIndex
    type(index) :: i
    type(dimension),pointer  :: d, cd
    integer irc
    allocate(copyIndex,stat=irc)
    allocate(copyIndex%firstDimension,copyIndex%lastDimension, stat=irc)
    copyIndex%firstDimension%next => copyIndex%lastDimension
    copyIndex%lastDimension%prev  => copyIndex%firstDimension
    copyIndex%nrdim =0
    d=>i%firstDimension%next
    do while (.not.associated(d,i%lastDimension))
       allocate(cd,stat=irc)
       cd%ind=d%ind
       cd%lim=d%lim
       cd%sta=d%sta
       cd%prev=>copyIndex%lastDimension%prev
       cd%next=>copyIndex%lastDimension
       copyIndex%lastDimension%prev%next=>cd
       copyIndex%lastDimension%prev=>cd
       copyIndex%nrdim = copyIndex%nrdim + 1
       nullify(cd)
       d=>d%next
    end do
  end function copyIndex

  function createDimension(f,dim,lim)
    ! creates new dimension and puts in into index
    type(index),pointer :: createDimension
    type(inventory) :: f
    character*(*) dim
    integer lim
    type(dimension),pointer  :: d
    integer irc
    character*250 xdim250(f%nrdim+1)
    integer xlend(f%nrdim+1), xpos(f%nrdim+1), xsta(f%nrdim+1), xlim(f%nrdim+1)
    integer ii
    do ii=1,f%nrdim
       xdim250(ii)=f%dim250(ii)
       xlend(ii)=f%lend(ii)
       xpos(ii)=f%pos%pos(ii)
       xsta(ii)=f%pos%sta(ii)
       xlim(ii)=f%pos%lim(ii)
    end do
    if (associated(f%dim250)) deallocate(f%dim250,stat=irc)
    if (associated(f%lend)) deallocate(f%lend,stat=irc)
    if (associated(f%pos%pos)) deallocate(f%pos%pos,stat=irc)
    if (associated(f%pos%sta)) deallocate(f%pos%sta,stat=irc)
    if (associated(f%pos%lim)) deallocate(f%pos%lim,stat=irc)
    f%nrdim=f%nrdim+1
    xdim250(f%nrdim)=dim
    xlend(f%nrdim)=len(dim)
    xpos(f%nrdim)=1
    xsta(f%nrdim)=1
    xlim(f%nrdim)=lim
    f%pos%nrdim=f%nrdim
    f%pos%maxdim=f%nrdim
    allocate(f%dim250(f%nrdim),stat=irc)
    allocate(f%lend(f%nrdim),stat=irc)
    allocate(f%pos%pos(f%pos%maxdim),stat=irc)
    allocate(f%pos%sta(f%pos%maxdim),stat=irc)
    allocate(f%pos%lim(f%pos%maxdim),stat=irc)
    do ii=1,f%nrdim
       f%dim250(ii)  = xdim250(ii)
       f%lend(ii)    = xlend(ii)
       f%pos%pos(ii) = xpos(ii)
       f%pos%sta(ii) = xsta(ii)
       f%pos%lim(ii) = xlim(ii)
    end do
    allocate(createDimension,stat=irc)
    createDimension%nrdim=1
    allocate(createDimension%firstDimension,createDimension%lastDimension, stat=irc)
    createDimension%firstDimension%next => createDimension%lastDimension
    createDimension%lastDimension%prev  => createDimension%firstDimension
    allocate(d,stat=irc)
    d%ind=f%nrdim
    d%sta=f%pos%sta(f%nrdim)
    d%lim=f%pos%lim(f%nrdim)
    d%prev=>createDimension%lastDimension%prev
    d%next=>createDimension%lastDimension
    createDimension%lastDimension%prev%next=>d
    createDimension%lastDimension%prev=>d
    nullify(d)
  end function createDimension

  subroutine delIndex(t,i)
    ! deletes index "i" from "t"
    type(index) :: t
    type(index) :: i
    type(dimension),pointer  :: td, id, x
    integer irc
    id=>i%firstDimension%next
    do while (.not.associated(id,target=i%lastDimension))
       td=>t%firstDimension%next
       do while (.not.associated(td,target=t%lastDimension))
          if (id%ind .eq. td%ind) then
             x => td
             td=>td%next
             t%nrdim = t%nrdim - 1
             x%prev%next=> x%next
             x%next%prev=> x%prev 
             nullify(x%prev,x%next)
             deallocate(x,stat=irc)
             nullify(x)
          else
             td=>td%next
          end if
       end do
       id=>id%next
    end do
  end subroutine delIndex

  subroutine printIndex(i)
    type(index) :: i
    type(dimension),pointer  :: id
    id=>i%firstDimension%next
    do while (.not.associated(id,i%lastDimension))
       write(*,*) 'printIndex dim:',id%ind,id%sta,id%lim
       id=>id%next
    end do
  end subroutine printIndex

  subroutine addIndex(t,i)
    ! adds index "i" to "t"
    type(index) :: t
    type(index) :: i
    type(dimension),pointer  :: td, id, x
    integer irc
    id=>i%firstDimension%next
    do while (.not.associated(id,i%lastDimension))
       found=.false.
       td=>t%firstDimension%next
       do while (.not.associated(td,t%lastDimension))
          if (id%ind .eq. td%ind) then
             found=.true.
             td=>t%lastDimension
          else
             td=>td%next
          end if
       end do
       if (.not.found) then ! add dimension to the end
          allocate(x,stat=irc)
          x%ind=id%ind
          x%lim=id%lim
          x%sta=id%sta
          x%prev=>t%lastDimension%prev
          x%next=>t%lastDimension
          t%lastDimension%prev%next=>x
          t%lastDimension%prev=>x
          t%nrdim = t%nrdim + 1
          nullify(x)
       end if
       id=>id%next
    end do
  end subroutine addIndex

  function unionIndex(t,i,irc)
    type(index),pointer :: unionIndex
    ! find union dimensions in t and i
    type(index),pointer  :: t
    type(index),pointer  :: i
    type(dimension),pointer  :: td, id, x
    integer irc
    allocate(unionIndex,stat=irc)
    allocate(unionIndex%firstDimension,unionIndex%lastDimension, stat=irc)
    unionIndex%firstDimension%next => unionIndex%lastDimension
    unionIndex%lastDimension%prev  => unionIndex%firstDimension
    unionIndex%nrdim =0
    if (associated(t).and.associated(i)) then
       id=>i%firstDimension%next
       do while (.not.associated(id,i%lastDimension))
          found=.false.
          td=>t%firstDimension%next
          do while (.not.associated(td,t%lastDimension))
             if (id%ind .eq. td%ind) then
                found=.true.
                td=>t%lastDimension
             else
                td=>td%next
             end if
          end do
          if (found) then ! add dimension to the end
             allocate(x,stat=irc)
             x%ind=id%ind
             x%lim=id%lim
             x%sta=id%sta
             x%prev=>unionIndex%lastDimension%prev
             x%next=>unionIndex%lastDimension
             unionIndex%lastDimension%prev%next=>x
             unionIndex%lastDimension%prev=>x
             unionIndex%nrdim = unionIndex%nrdim + 1
             nullify(x)
          end if
          id=>id%next
       end do
    end if
  end function unionIndex

  function makeWeight4(nrdim,irc)
    type(weight), pointer :: makeWeight4
    integer nrdim
    integer irc
    integer kk
    allocate(makeWeight4,stat=irc)
    makeWeight4%nweight=4
    allocate(makeWeight4%pos(makeWeight4%nweight),stat=irc)
    do kk=1,makeWeight4%nweight
       call allocatePos(nrdim,makeWeight4%pos(kk),irc)
    end do
    allocate(makeWeight4%w(makeWeight4%nweight),stat=irc)
  end function makeWeight4

  subroutine allocatePos(nrdim,pos,irc)
    integer nrdim
    type(position) :: pos
    integer irc
    pos%maxdim=nrdim
    pos%nrdim=nrdim
    allocate(pos%pos(pos%nrdim),stat=irc)
    allocate(pos%sta(pos%nrdim),stat=irc)
    allocate(pos%lim(pos%nrdim),stat=irc)
    return
  end subroutine allocatePos

  subroutine checkGrid(f,t,irc)
    implicit none
    type(inventory) :: f
    type(inventory) :: t
    integer irc
    integer kk
    if (.not.associated(f%gridid)) then
       write(*,*)'checkGrid No grid defined in F.'
       irc=876
       return
    end if
    if (.not.associated(t%gridid)) then
       write(*,*)'checkGrid No grid defined in T.'
       irc=877
       return
    end if
    if (f%gridid%nrdim.ne.t%gridid%nrdim) then
       write(*,*) 'checkGrid Grid dimension mismatch:',f%gridid%nrdim,t%gridid%nrdim
       irc=878
       return
    end if
    do kk=1,f%gridid%nrdim
       if (f%gridid%lim(kk).ne.f%gridid%lim(kk)) then
          write(*,*) 'checkGrid Dimension length mismatch:',kk,f%gridid%lim(kk),t%gridid%lim(kk)
          irc=879
          return
       end if
    end do
  end subroutine checkGrid

  logical function incrementDim(f,d,biok,inc,buff,irc)
    implicit none
    type(inventory) :: f
    type(dimension) :: d
    logical biok
    integer inc
    integer buff
    integer irc
    logical dok
    if (biok) then
       f%pos%pos(d%ind)=f%pos%pos(d%ind)+inc
       if (f%pos%pos(d%ind).gt.d%lim-buff) then
          f%pos%pos(d%ind)=d%lim-buff
          biok=.false.
       else
          biok=.true.
       end if
       f%pos%loc=0
    end if
    incrementDim=biok
  end function incrementDim

  logical function decrementDim(f,d,biok,inc,buff,irc)
    implicit none
    type(inventory) :: f
    type(dimension) :: d
    logical biok
    integer inc
    integer buff
    integer irc
    if (biok) then
       f%pos%pos(d%ind)=f%pos%pos(d%ind)-inc
       if (f%pos%pos(d%ind).lt.d%sta+buff) then
          f%pos%pos(d%ind)=d%sta+buff
          biok=.false.
       else
          biok=.true.
       end if
       f%pos%loc=0
    end if
    decrementDim=biok
  end function decrementDim

  logical function increment(f,i,irc)
    implicit none
    type(inventory) :: f
    type(index) :: i
    integer irc
    type(dimension),pointer :: d
    increment=.false.
    if (i%resetIndex) then
       d=>i%firstDimension%next
       do while (.not.associated(d,target=i%lastDimension))
          f%pos%pos(d%ind)=d%sta
          d%increase=.true.
          d=>d%next
       end do
       increment=.true. ! there is always one valid return
       i%resetIndex=.false.
    else
       d=>i%firstDimension%next
       do while (.not.associated(d,target=i%lastDimension))
          if (d%ind.le.f%nrdim.and.d%ind.ge.1) then
             if (d%increase) then
                f%pos%pos(d%ind)=f%pos%pos(d%ind)+1
                if (f%pos%pos(d%ind).gt.d%lim) then
                   d%increase=.false.
                   f%pos%pos(d%ind)=d%lim
                   d=>d%next
                else
                   increment=.true.
                   d=>i%lastDimension
                end if
             else
                f%pos%pos(d%ind)=f%pos%pos(d%ind)-1
                if (f%pos%pos(d%ind).lt.d%sta) then
                   d%increase=.true.
                   f%pos%pos(d%ind)=d%sta
                   d=>d%next
                else
                   increment=.true.
                   d=>i%lastDimension
                end if
             end if
          else
             irc=945
             return
          end if
       end do
    end if
    f%pos%loc=0
  end function increment

  subroutine resetIndex(i,irc)
    implicit none
    type(index) :: i
    integer irc
    i%resetIndex=.true.
  end subroutine resetIndex
!     
!     find closest i and j
!     
  real function xycross(x1,y1,x2,y2,x3,y3)
    real :: x1,y1,x2,y2,x3,y3
    real :: dx1,dy1,dx2,dy2
    dx1=x3-x1
    dy1=y3-y1
    dx2=x2-x1
    dy2=y2-y1
    xycross = dx1 * dy2 - dx2 * dy1
  end function xycross
  
  real function degtor(x)
    implicit none
    real x
    real pi
    parameter (pi=3.14159265359)
    degtor=x*pi/180.
  end function degtor
  
  real function rtodeg(x)
    implicit none
    real x
    real pi
    parameter (pi=3.14159265359)
    rtodeg=x*180./pi
  end function rtodeg
  
  real function sindeg(x)
    implicit none
    real x,degtor
    sindeg=sin(degtor(x))
  end function sindeg
  
  real function cosdeg(x)
    implicit none
    real x,degtor
    cosdeg=cos(degtor(x))
  end function cosdeg
  
  real function tandeg(x)
    implicit none
    real x,degtor
    tandeg=tan(degtor(x))
  end function tandeg
  
  real function asindeg(x)
    implicit none
    real x,rtodeg
    asindeg=rtodeg(asin(x))
  end function asindeg
  
  real function acosdeg(x)
    implicit none
    real x,rtodeg
    acosdeg=rtodeg(acos(x))
  end function acosdeg
  
  real function atandeg(x)
    implicit none
    real x,rtodeg
    atandeg=rtodeg(atan(x))
  end function atandeg
  
  real function atan2deg(y,x)
    implicit none
    real y,x,rtodeg
    atan2deg=rtodeg(atan2(y,x))
  end function atan2deg
  
  function getDim(i,num)
    type(dimension), pointer :: getDim
    type(index) :: i
    integer num
    integer kk
    type(dimension), pointer :: d
    nullify(getDim)
    kk=0
    d=>i%firstDimension%next
    do while (.not.associated(d,i%lastDimension))
       kk=kk+1
       if (kk.eq.num) then
          getDim=>d
          d=>i%lastDimension
       else
          d=>d%next
       end if
    end do
  end function getDim

  subroutine interpolate2D(r,e,ix,iy,wgt,biok,irc)
    implicit none
    type(inventory) r
    type(inventory) e
    type(dimension) ix
    type(dimension) iy
    type(weight) wgt
    integer irc
    real :: xf,yf
    integer :: nxs,nys,nzs
    real :: xs(2,2), ys(2,2)
    logical :: biok, xok,yok
    logical :: bdone,bdeb, changed, inside
    real :: bc,rc,tc,lc,fc,tbc,lrc
    real :: brc,rtc,tlc,lbc
    integer :: jbc,jrc,jtc,jlc,jbrc,jrtc,jtlc,jlbc,ii,jj
    integer :: iterations = 0
    !e%pos%pos(ix%ind)=max(ix%sta,min(ix%lim-1,e%pos%pos(ix%ind)))
    !e%pos%pos(iy%ind)=max(iy%sta,min(iy%lim-1,e%pos%pos(iy%ind)))
    ! search for match where position is within gridcell
    xf=valuePosition(r%lonid,irc)
    if (irc.ne.0) then
       write(*,*)'interpolate2D error return from valuePosition.',irc
       return
    end if
    yf=valuePosition(r%latid,irc)
    if (irc.ne.0) then
       write(*,*)'interpolate2D error return from valuePosition.',irc
       return
    end if
    bdone=.false.
    bdeb=.false.
    ! iterations=0
    do while (.not. bdone)
       iterations=iterations+1
       xok=.true.
       yok=.true.
       xok=incrementDim(e,ix,xok,1,0,irc) ! -> bottom right
       xs(2,1)=valuePosition(e%lonid,irc)
       ys(2,1)=valuePosition(e%latid,irc)

       yok=incrementDim(e,iy,yok,1,0,irc) ! -> top right
       xs(2,2)=valuePosition(e%lonid,irc)
       ys(2,2)=valuePosition(e%latid,irc)

       xok=decrementDim(e,ix,xok,1,0,irc) ! -> top left
       xs(1,2)=valuePosition(e%lonid,irc)
       ys(1,2)=valuePosition(e%latid,irc)

       yok=decrementDim(e,iy,yok,1,0,irc) ! -> bottom left
       xs(1,1)=valuePosition(e%lonid,irc)
       ys(1,1)=valuePosition(e%latid,irc)

       if (.not. xok.or. .not.yok) then
          biok=.false.
          xok=decrementDim(e,ix,xok,1,0,irc) ! -> top left
          yok=decrementDim(e,iy,yok,1,0,irc) ! -> bottom left
          call printPos(e%pos)
          write(*,*) 'interpolate2D Invalid starting position.'
          irc=945
          return
       else
          !write(*,*)'XS:',xf,xs
          !write(*,*)'YS:',yf,ys
          !call exit(0)
          
          ! walk around grid cell border and calculate cross product
          bc=xycross(xs(1,1),ys(1,1), xf,yf, xs(2,1),ys(2,1)) ! bottom
          rc=xycross(xs(2,1),ys(2,1), xf,yf, xs(2,2),ys(2,2)) ! right
          tc=xycross(xs(2,2),ys(2,2), xf,yf, xs(1,2),ys(1,2)) ! top
          lc=xycross(xs(1,2),ys(1,2), xf,yf, xs(1,1),ys(1,1)) ! left
          
          if (abs(bc).lt.1.0D-10) then
             jbc=0
          else
             jbc=nint(sign(1.0D0,bc))
          end if
          if (abs(rc).lt.1.0D-10) then
             jrc=0
          else
             jrc=nint(sign(1.0D0,rc))
          end if
          if (abs(tc).lt.1.0D-10) then
             jtc=0
          else
             jtc=nint(sign(1.0D0,tc))
          end if
          if (abs(lc).lt.1.0D-10) then
             jlc=0
          else
             jlc=nint(sign(1.0D0,lc))
          end if
          
          brc=xycross(xs(1,1),ys(1,1), xs(2,2),ys(2,2), xs(2,1),ys(2,1)) ! bottom-right
          rtc=xycross(xs(2,1),ys(2,1), xs(1,2),ys(1,2), xs(2,2),ys(2,2)) ! right-top
          tlc=xycross(xs(2,2),ys(2,2), xs(1,1),ys(1,1), xs(1,2),ys(1,2)) ! top-left
          lbc=xycross(xs(1,2),ys(1,2), xs(2,1),ys(2,1), xs(1,1),ys(1,1)) ! left-bottom
          
          jbrc=nint(sign(1.0D0,brc))
          jrtc=nint(sign(1.0D0,rtc))
          jtlc=nint(sign(1.0D0,tlc))
          jlbc=nint(sign(1.0D0,lbc))
          
          if (bdeb) then
             write(*,'("XYSEARCH ",2(3X,"(",4(X,F13.5),")"),2(3X,"(",4(X,I2),")"))') &
                  &bc,rc,tc,lc, brc,rtc,tlc,lbc,jbc,jrc,jtc,jlc, jbrc,jrtc,jtlc,jlbc
             write(*,'("XYSEARCH ",5(" (",F10.5,",",F10.5,")"))') &
                  & xf,yf, &
                  & xs(1,1),ys(1,1), &
                  & xs(2,1),ys(2,1), &
                  & xs(2,2),ys(2,2), &
                  & xs(1,2),ys(1,2)
          end if
          !     
          !     if cross product has same sign as axis-sign => inside, else outside
          !     (make sure we do not cross border while searching...)
          !     
          inside=.true.
          changed=.false.
          if (jbc.eq.-jbrc) then  ! decrease y
             if (decrementDim(e,iy,yok,1,0, irc)) then
                changed=.true.

                !write(*,*)'interpolate2d Y-',e%pos%pos(iy%ind),iterations,bc,brc

             end if
             inside=.false.
          else if (jtc.eq.-jtlc) then ! increase y
             if (incrementDim(e,iy,yok,1,1, irc)) then
                changed=.true.


                !write(*,*)'interpolate2d Y+',e%pos%pos(iy%ind),iterations,tc,tlc

             end if
             inside=.false.
          end if
          if (jlc.eq.-jlbc) then  ! decrease x
             if (decrementDim(e,ix,xok,1,0, irc)) then
                changed=.true.

                !write(*,*)'interpolate2d X-',e%pos%pos(ix%ind),iterations,lc,lbc

             end if
             inside=.false.
          else if (jrc.eq.-jrtc) then ! increase x
             if (incrementDim(e,ix,xok,1,1, irc)) then
                changed=.true.

                !write(*,*)'interpolate2d X+',e%pos%pos(ix%ind),iterations,rc,rtc

             end if
             inside=.false.
          end if
          if (inside) then       ! we are inside the cell
             biok=.true.
             bdone=.true.
             tbc=tc+bc
             lrc=lc+rc
             fc=(tbc*lrc)
             fc=max(fc,1.0D-10)
             if (abs(fc).lt.1.0D-10) then
                call printPos(e%pos)
                write(*,*) myname,'Invalid grid:',rc,tc,lc,bc,fc
                irc=945
                return
             end if
             wgt%w(1)=lc*tc/fc ! bottom right
             wgt%w(2)=lc*bc/fc ! top right
             wgt%w(3)=rc*bc/fc ! top left
             wgt%w(4)=rc*tc/fc ! bottom left
             ! store position vector
             do ii=1,wgt%nweight
                if (wgt%pos(ii)%maxdim .lt. e%pos%nrdim) then ! sufficient dimensions?
                   call allocatePos(e%pos%nrdim,wgt%pos(ii),irc)
                else
                   wgt%pos(ii)%nrdim=e%pos%nrdim
                end if
                do jj=1,wgt%pos(ii)%nrdim
                   wgt%pos(ii)%pos(jj)=0
                end do
             end do
             wgt%pos(1)%pos(ix%ind)=+1 ! bottom right
             wgt%pos(2)%pos(iy%ind)=+1 ! top right
             wgt%pos(3)%pos(ix%ind)=-1 ! top left
             wgt%pos(4)%pos(iy%ind)=-1 ! bottom left (origo)
             ! write(*,'(X,A,6(F10.3,X))') 'Search done:',xs(1,1),xs(2,2),xf,ys(1,1),ys(2,2),yf
             if (bdeb) write(*,*)"XYSEARCH Done:",wgt%w
          else if (.not.changed) then ! not inside + no valid step
             if (bdeb) then
                write(*,'(X,A,6(F10.3,X))') 'Search fail:',xs(1,1),xs(2,2),xf,ys(1,1),ys(2,2),yf
             end if
             bdone=.true.
             biok=.false.
             return
          end if
       end if
    end do    

!    if (mod(iterations,10000).eq.0) write(*,*) 'Interpolate2d iterations:',iterations
!    if (mod(iterations,100000).eq.0) call printpos(e%pos)

    return
  end subroutine interpolate2D
  
  subroutine printWeight(wgt)
    implicit none
    type(weight) :: wgt
    integer ii
    do jj=1,wgt%nweight
       write(*,*) 'printWeight Weight:',jj,wgt%w(jj)
       do ii=1,wgt%pos(jj)%nrdim
          write(*,*) 'printWeight      Pos:',ii,wgt%pos(jj)%pos(ii)
       end do
    end do
  end subroutine printWeight

  subroutine printPos(pos)
    implicit none
    type(position) :: pos
    integer ii
    if (pos%nrdim .gt.0) then
       do ii=1,pos%nrdim
          write(*,'(X,A," (",I3,")=",I8,"  (",I4,I4," )")') 'printPos',ii,pos%pos(ii),pos%sta(ii),pos%lim(ii)
       end do
    else
       write(*,*)'printPos: Position has no dimensions!'
    end if
    return
  end subroutine printPos


  character*6 function gettype(type)
    integer type
    if (type.eq.nf_char) then
       gettype="char  "
    elseif (type.eq.nf_int1) then
       gettype="int1  "
    elseif (type.eq.nf_int2) then
       gettype="int2  "
    elseif (type.eq.nf_int) then
       gettype="int   "
    elseif (type.eq.nf_real) then
       gettype="real  "
    elseif (type.eq.nf_double) then
       gettype="double"
    else
       gettype="any   "
    end if
    return
  end function gettype

  real function getDist(LATA,LONA,LATB,LONB)
    IMPLICIT NONE
    SAVE
    REAL   LONA,LATA,LONB,LATB
    REAL   CDIFF
    real sindeg, cosdeg,acosdeg
    external sindeg,cosdeg,acosdeg
    CHARACTER*8 MYNAME 
    DATA MYNAME /'DIST'/
    CDIFF=SINDEG(LATA)*SINDEG(LATB)+COSDEG(LATA)*COSDEG(LATB)*COSDEG(LONA-LONB)
    CDIFF=MAX(-1.0D0,MIN(1.0D0,CDIFF)) ! handle truncation errors
    getDist=ABS(ACOSDEG(CDIFF))
    RETURN
  END function getDist
  

  logical function getFootIndex(f,i,s,latc,lonc,diam,irc)
    implicit none
    type(inventory), pointer    :: f
    type(index), pointer        :: i
    type(index), pointer        :: s
    real latc
    real lonc
    real diam ! diameter in degrees
    integer irc
    integer dn(2),ii
    logical biok
    real latx,lonx,dist
    type(dimension), pointer :: d,e
! find increment distance in each dimension direction
    biok=.true.
    latc=valuePosition(f%latid,irc)
    lonc=valuePosition(f%lonid,irc)
    do ii=1,i%nrdim ! must be 2 (checked earlier)
       d=>getdim(i,ii)
       e=>getdim(s,ii)
       if (d%ind .ne. e%ind) then
          write(*,*)'System error:',d%ind,e%ind
          irc=845
          return
       end if
       e%sta=f%pos%pos(e%ind)
       if (biok) biok=incrementDim(f,d,biok,1,0,irc)
       if (biok) then
          latx=valuePosition(f%latid,irc)
          lonx=valuePosition(f%lonid,irc)
          biok=decrementDim(f,d,biok,1,0,irc)
          dist=2.0D0*getDist(latc,lonc,latx,lonx)
          dn(ii)=max(1,1+int(diam/max(dist,1.0D-10)))
       end if
    end do
! determine lower left corner
    if (biok) then
       do ii=1,i%nrdim ! must be 2 (checked earlier)
          d=>getdim(i,ii)
          biok=decrementDim(f,d,biok,dn(ii),0,irc)
          e=>getdim(s,ii)
          e%sta=e%sta-dn(ii)
          e%lim=e%sta+2*dn(ii)
          if (e%sta.lt.1.or.e%lim.gt.f%pos%lim(e%ind)) biok=.false.
       end do
    end if
    !    call printIndex(s)
    getFootIndex=biok
  end function getFootIndex

  logical function incrementFoot(f,i,latc,lonc,diam,irc)
    implicit none
    type(inventory), pointer    :: f
    type(index), pointer        :: i
    real latc
    real lonc
    real diam
    logical reset
    integer irc
    real latx,lonx,dist
    logical bdone
    bdone=.false.
    do while (.not. bdone)
       if (increment(f,i,irc)) then
          ! call printIndex(i)
          ! call printPos(f%pos)
          latx=valuePosition(f%latid,irc)
          lonx=valuePosition(f%lonid,irc)
          dist=2.0D0*getDist(latc,lonc,latx,lonx)
          if (dist.lt.diam) then ! grid wihin footprint...
             bdone=.true.
          end if
          incrementFoot=.true.
       else !     no more grid points to check
          incrementFoot=.false.
          bdone=.true.
       end if
    end do
  end function incrementFoot

  subroutine clearIndex(i)
    type(index), pointer  :: i
    type(dimension), pointer :: d
    if (associated(i)) then
       d=>i%firstDimension%next
       do while (.not.associated(d,i%lastDimension))
          deallocate(d%prev)
          d=>d%next
       end do
       deallocate(d)
       deallocate(i)
    end if
  end subroutine clearIndex

  subroutine writeNcOut(outnc250,out,irc)
    implicit none
    character*250  outnc250
    type(inventory) :: out
    integer irc
    integer length, leno
    external length
    character*14 :: myname = "writeNcOut"
    call createNcFile(out,outnc250,irc)
    call writeNcDims(out,irc)
    call writeNcVars(out,irc)
    call writeNcAtts(out,out%firstAttribute,out%lastAttribute,irc) ! global attributes
    call endNcDef(out,irc)
    call writeNcData(out,irc)
    call closeNcFile(out,irc)
    return
  end subroutine writeNcOut

  subroutine createNcFile(out,outnc250,irc)
    implicit none
    character*250  outnc250
    type(inventory) :: out
    integer irc
    integer length,leno, ret
    external length
    character*14 :: myname = "createNcFile"
    call chop(outnc250,250)
    leno=length(outnc250,250,10)
    out%fn250=outnc250
    out%lenf=leno
    ! or(NF_NOCLOBBER,NF_SHARE) - do not overwrite existing file but give error message
    ret=NF_CREATE(out%fn250(1:out%lenf),NF_SHARE,out%ncid) ! create netCDF dataset: enter define mode
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_CREATE:",nf_strerror(ret),out%ncid
       irc=790
       return
    end if
  end subroutine createNcFile

  subroutine writeNcDims(out,irc)
    implicit none
    type(inventory) :: out
    integer irc
    integer ii,ret,leng,length
    external length
    character*14 :: myname = "writeNcDims"
    allocate(out%dimid(out%nrdim),stat=irc)
    do ii=1,out%nrdim
       leng=length(out%dim250(ii),250,10)
       if (out%unlimdimid.eq.ii) then
          ret=NF_DEF_DIM(out%ncid,out%dim250(ii)(1:leng),nf_unlimited,out%dimid(ii))
       else
          ret=NF_DEF_DIM(out%ncid,out%dim250(ii)(1:leng),out%pos%lim(ii),out%dimid(ii))
       end if
       if (ret .ne. NF_NOERR) then
          write(*,*) myname,"ERROR from NF_DEF_DIM:",nf_strerror(ret),out%ncid
          irc=790
          return
       end if
    end do
  end subroutine WRITENCDIMS

  subroutine writeNcVars(out,irc)
    implicit none
    type(inventory) :: out
    integer irc
    integer ret
    type(variable), pointer :: v
    integer gdims(out%nrdim)
    character*14 :: myname = "writeNcVars"
    integer ii,jj
    v=>out%firstVariable%next
    do while(.not.associated(v,target=out%lastVariable))
       do ii=1,v%nrdim
          gdims(ii)=out%dimid(v%ind(ii))
       end do
       ret=NF_DEF_VAR(out%NCID,v%var250(1:v%lenv),v%type,v%nrdim,gdims,v%varid)
       if (ret .ne. NF_NOERR) then
          write(*,*) myname,"ERROR from NF_DEF_VAR:",nf_strerror(ret),out%ncid
          irc=790
          return
       end if
       call writeNcAtts(out,v%firstAttribute,v%lastAttribute,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from writeNcAtts.',irc
          return
       end if
       v=>v%next
    end do
  end subroutine writeNcVars

  subroutine writeNcAtts(out,first,last,irc)
    implicit none
    type(inventory) :: out
    type(attribute), pointer :: first,last
    integer irc
    type(attribute), pointer :: a
    integer varid,ret
    character*14 :: myname = "writeNcAtts"
    a => first%next
    do while (.not. associated(a,target=last))
       if (associated(a%v)) then ! variable attribute
          varid=a%v%varid
       else
          varid=NF_GLOBAL
       end if
       if (a%type.eq.nf_char) then
          ret=NF_PUT_ATT_TEXT  (OUT%NCID, varid,a%att250(1:a%lena),a%LEN, a%ac)
       else if (a%type.eq.nf_int1) then
          ret=NF_PUT_ATT_INT1  (OUT%NCID, varid,a%att250(1:a%lena),a%type,a%LEN, a%a1)
       else if (a%type.eq.nf_int2) then
          ret=NF_PUT_ATT_INT2  (OUT%NCID, varid,a%att250(1:a%lena),a%type,a%LEN, a%a2)
       else if (a%type.eq.nf_int) then
          ret=NF_PUT_ATT_INT  (OUT%NCID, varid,a%att250(1:a%lena),a%type,a%LEN, a%a4)
       else if (a%type.eq.nf_real) then
          ret=NF_PUT_ATT_REAL  (OUT%NCID, varid,a%att250(1:a%lena),a%type,a%LEN, a%ar)
       else if (a%type.eq.nf_DOUBLE) then
          ret=NF_PUT_ATT_DOUBLE  (OUT%NCID, varid,a%att250(1:a%lena),a%type,a%LEN, a%ad)
       end if
      if (ret .ne. NF_NOERR) then
         write(*,*) myname,"ERROR-2 from NF_PUT_ATT_*:",nf_strerror(ret),out%ncid
         irc=790
         return
      end if
      a=>a%next
   end do
  end subroutine writeNcAtts

  subroutine endNcDef(out,irc)
    implicit none
    type(inventory) :: out
    integer irc
    integer ret
    character*14 :: myname = "endNcDef"
    ret=NF_ENDDEF(out%ncid)       ! end definitions: leave define mode
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_ENDDEF:",nf_strerror(ret),out%ncid
       irc=790
       return
    end if
  end subroutine endNcDef

  subroutine writeNcData(out,irc)
    implicit none
    type(inventory) :: out
    integer irc
    integer ret
    type(variable), pointer :: v
    integer gdims(out%nrdim)
    character*14 :: myname = "writeNcData"
    v=>out%firstVariable%next
    do while(.not.associated(v,target=out%lastVariable))
       if (v%type.eq.nf_char.and.associated(v%fc)) then
          ret= NF_PUT_VARA_TEXT(out%ncid,v%varid,v%sta,v%lim,v%fc)
       else if (v%type.eq.nf_int1.and.associated(v%f1)) then
          ret= NF_PUT_VARA_INT1(out%ncid,v%varid,v%sta,v%lim,v%f1)
       else if (v%type.eq.nf_int2.and.associated(v%f2)) then
          ret= NF_PUT_VARA_INT2(out%ncid,v%varid,v%sta,v%lim,v%f2)
       else if (v%type.eq.nf_int.and.associated(v%f4)) then
          ret= NF_PUT_VARA_INT(out%ncid,v%varid,v%sta,v%lim,v%f4)
       else if (v%type.eq.nf_real.and.associated(v%fr)) then
          ret= NF_PUT_VARA_REAL(out%ncid,v%varid,v%sta,v%lim,v%fr)
       else if (v%type.eq.nf_double.and.associated(v%fd)) then
          ret= NF_PUT_VARA_DOUBLE(out%ncid,v%varid,v%sta,v%lim,v%fd)
       else
          write(*,*)myname,'Undefined variable data: ',v%var250(1:v%lenv),v%type, &
               &nf_double,nf_real,nf_int,nf_int2,nf_int1, &
               &associated(v%fd),associated(v%fr),associated(v%f4),associated(v%f2),associated(v%f1)
          ret=nf_noerr
       end if
       if (ret .ne. NF_NOERR) then
          write(*,*) myname,"ERROR from NF_PUT_VARA_*:",nf_strerror(ret)
          irc=812
          return
       end if
       v => v%next
    end do
  end subroutine writeNcData

  subroutine closeNcFile(out,irc)
    implicit none
    type(inventory) :: out
    integer irc
    character*14 :: myname = "closeNcFile"
    integer ret
    ret=NF_CLOSE(out%ncid)        ! end definitions: leave define mode
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_CLOSE:",nf_strerror(ret),out%ncid
       irc=790
       return
    end if
  end subroutine closeNcFile

  subroutine copyField(v,c)
    implicit none
    type(variable), pointer :: v
    type(variable), pointer :: c
    integer ll
    do ll=1,v%lend
       v%fd(ll)=c%fd(ll)
    end do
  end subroutine copyField
  
  subroutine addField(v,c)
    implicit none
    type(variable) :: v
    type(variable) :: c
    integer ll
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=v%fd(ll)+c%fd(ll)
       end if
    end do
  end subroutine addField
  
  subroutine subtractField(v,c)
    implicit none
    type(variable), pointer :: v
    type(variable), pointer :: c
    integer ll
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=v%fd(ll)-c%fd(ll)
       end if
    end do
  end subroutine subtractField
  
  subroutine multiplyField(v,c)
    implicit none
    type(variable), pointer :: v
    type(variable), pointer :: c
    integer ll
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=v%fd(ll)*c%fd(ll)
       end if
    end do
  end subroutine multiplyField
  
  subroutine scaleField(v,scale)
    implicit none
    type(variable), pointer :: v
    real scale
    integer ll
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=v%fd(ll)*scale
       end if
    end do
  end subroutine scaleField
  
  subroutine divideField(v,c)
    implicit none
    type(variable), pointer :: v
    type(variable), pointer :: c
    integer ll
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double .and. abs(c%fd(ll)).gt.1.0D-10) then
          v%fd(ll)=v%fd(ll)/c%fd(ll)
       end if
    end do
  end subroutine divideField
  
  subroutine sqrtField(v)
    implicit none
    type(variable), pointer :: v
    integer ll
    do ll=1,v%lend
       v%fd(ll)=sqrt(max(0.0D0,v%fd(ll)))
    end do
  end subroutine sqrtField
  
  subroutine roundField(v,prec)
    implicit none
    type(variable),pointer  :: v
    real prec
    integer ll
    real xprec
    xprec=max(1.0D-10,prec)
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=DBLE(nint(v%fd(ll)/xprec))*xprec
       end if
    end do
  end subroutine roundField
  
  subroutine addAttribute(v,nam,val)
    type(variable), target :: v
    character*(*) :: nam
    character*(*) :: val
    type(attribute), pointer :: a
    integer irc
    integer ll
    allocate(a,stat=irc)
    if (irc.ne.0) then
       write(*,*) 'addAttribute Unable to allocate ATTRIBUTE.',irc
       return
    end if
    a%att250=nam
    a%lena=len(nam)
    a%type=nf_char
    a%len = len(val)
    allocate(a%ac(a%len),stat=irc)
    do ll=1,len(val)
       a%ac(ll)=val(ll:ll)
    end do
    a%v => v
    a%prev=>v%lastAttribute%prev
    a%next=>v%lastAttribute
    v%lastAttribute%prev%next=>a
    v%lastAttribute%prev=>a
    nullify(a)
  end subroutine addAttribute

  subroutine addGlobalAttribute(f,nam,val)
    type(inventory), target :: f
    character*(*) :: nam
    character*(*) :: val
    type(attribute), pointer :: a
    integer irc
    integer ll
    allocate(a,stat=irc)
    if (irc.ne.0) then
       write(*,*) 'addAttribute Unable to allocate ATTRIBUTE.',irc
       return
    end if
    a%att250=nam
    a%lena=len(nam)
    a%type=nf_char
    a%len = len(val)
    allocate(a%ac(a%len),stat=irc)
    do ll=1,len(val)
       a%ac(ll)=val(ll:ll)
    end do
    a%v => v
    a%prev=>f%lastAttribute%prev
    a%next=>f%lastAttribute
    f%lastAttribute%prev%next=>a
    f%lastAttribute%prev=>a
    nullify(a)
  end subroutine addGlobalAttribute

  subroutine initGrid(out,ref,irc)
    implicit none
    type(inventory) :: out
    type(inventory),pointer :: ref
    integer irc
    !
    type(variable), pointer :: v
    logical bbok
    !
    ! make sure reference grid data is read into memory
    !
    v=>ref%xcoid
    bbok=associated(v)
    call readData(v,bbok,irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from readData.',irc
       return
    end if
    v=>ref%ycoid
    bbok=associated(v)
    call readData(v,bbok,irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from readData.',irc
       return
    end if
    v=>ref%gridid
    bbok=associated(v)
    call readData(v,bbok,irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from readData.',irc
       return
    end if
    !
    ! make a copy of the reference grid for storing output
    !
    call copyGrid(out,ref,irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from copyGrid.',irc
       return
    end if
    return
  end subroutine initGrid

end subroutine mncverify

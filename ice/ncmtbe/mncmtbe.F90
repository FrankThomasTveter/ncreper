! ecmwf netcdf filer...   ls /fou/nwparc2/ec/2013/.nc
SUBROUTINE MNCMTBE(UNITI,IRC)
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
  CHARACTER*14 MYNAME
  DATA MYNAME /'MNCMTBE'/
  ! 
  LOGICAL  LFLDAT(100)
  logical :: bdeb=.false.
  ! 
  character*250 :: nam250
  INTEGER  lend, lenh, leni, lene, lenr, leno, lens, lenn
  integer, external :: length
  ! 
  INTEGER  KODE, LINE
  INTEGER  NRHDR
  PARAMETER (NRHDR=100)
  CHARACTER*250 HDR250(NRHDR),DAT250, INF250,buff250,buffx250
  CHARACTER*250 FILENM,NUKEHEAD
  EXTERNAL NUKEHEAD
  INTEGER  INTOUT
  LOGICAL  ENDOFF
  !
  type(inventory),pointer  :: iEnsInv => null()
  type(inventory),pointer  :: oMtbeInv => null()
  type(inventory),pointer  :: iStaInv => null()
  type(variable), pointer  :: oMtbeFrtVar,iStaFrtVar,iEnsFrtVar
  type(variable), pointer  :: oMtbeTimeVar,iStaTimeVar,iEnsTimeVar
  type(variable), pointer  :: iStaThrVar
  type(variable), pointer :: oMtbeVar,iStaVar,iEnsVar,v,vn
  type(dimension), pointer  :: iStaStatDim,iStaResDim,iStaThrDim,iStaAvgDim
  type(dimension), pointer  :: oMtbeTimeDim,iStaTimeDim,iEnsTimeDim
  type(dimension), pointer  :: iEnsDim,iEnsAvgDim
  type(dimensionOrder),pointer  ::  iEnsDo,iEnsTimeDo,iEnsAvgDo,oMtbedo,iStado
  integer :: iStaStaEntry,iStaResEntry,iStaThrEntry,oMtbeTimeEntry,iEnsAvgEntry,iStaAvgEntry
  integer :: iStaTimeEntry,iEnsTimeEntry,ix,iy,rx,ry
  real :: oMtbeFrtValue,iStaFrtValue                          ! latest analysis time
  logical :: bok,bbok,bdone,firstOut, firstStat, firstivar, firstevar
  type(attribute), pointer :: a => null()
  !
  integer :: resolution
  real :: range
  real :: resfact
  real :: rpos, prob, p
  real :: gamma = 1.0D-8
  character*4, parameter :: sthr="_thr"
  !
  integer, parameter :: nt=6
  real :: thres(nt),utcnt(0:nt),ltcnt(0:nt)
  data thres / 0.0D0, 0.1D0, 1.0D0, 10.0D0, 100.0D0, 1000.0D0 /
  ! 
  integer fok(10), frm(10)
  integer ook(20), orm(20)
  real pst(20)
  !
  !#     include "netcdf.inc"
  ! 
  INTEGER ret,CHUNKSIZEHINT
  ! 
  character*1::  fillc
  integer*1  ::  fill1
  integer*2  ::  fill2
  integer*4  ::  fill4
  real*4     ::  fillr
  real*8     ::  filld
  !
  character*250 :: ens250, out250, stat250, dim250
  integer :: lf, lc, lo, li
  integer :: lcnt,lmean,lstdv,lres,lfrt,laa,lbb,lcc,ldd
  integer :: ii, jj, ss, tt, uu, ik, jk, ecnt
  integer, allocatable :: io(:), is(:)
  real, external :: getu
  !
  real :: val,xval,valu,valy,vald,valm,vals,cnts,cmean,cstdv,cnt,minv,maxv,lgg
  real :: scnt,smean,sstdv,xcnt,xmean,xstdv,ucnt,umean,ustdv,caa,cbb,ccc,cdd
  !
  IRC=0
  ! 
  fillc=char(nf_fill_char)
  fill1=nf_fill_int1
  fill2=nf_fill_int2
  fill4=nf_fill_int
  fillr=nf_fill_real
  filld=nf_fill_double
  ! 
  if(bdeb)write(*,*) myname,'Debug: Routine starts.',irc
  !
  fok(:)=0
  frm(:)=0
  ook(:)=0
  orm(:)=0
  ! 
  DO II=1,NRHDR
     HDR250(II) = ''
     LFLDAT(II)=.FALSE.
  ENDDO

  HDR250(1)   = 'NCMTBE V1.0 [0]VFLR'
  HDR250(10)  = 'INPUT ENSEMBLE FILE (NETCDF) [1]VFLR &'
  HDR250(30)  = 'OUTPUT FILE (NETCDF) [1]VFLR &'
  HDR250(40)  = 'STATISTICS FILE (NETCDF) [1]VFLR &'
  HDR250(90)  = 'DUMP DATA [0]'
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
     ELSEIF (LINE.EQ.10) THEN ! input ensemble file
        ens250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.30) THEN ! output file netcdf
        out250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.40) THEN ! statistics file netcdf
        stat250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.90) THEN ! DUMP DATA
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
  firstout=.true.
  !
  ! allocate file inventories
  !
  allocate(oMtbeInv,iStaInv,stat=irc)
  if(irc.ne.0) then
     write(*,*)myname,'Unable to allocate inventory.',irc
     return
  end if
  !
  !
  ! Open input ENS file and read inventory ==========================
  !
  bok=.true.
  if (bok) then
     lene=length(ens250,250,10)
     write(*,*)myname,'Scanning: ',ens250(1:lene)
     call ncf_openFile(ens250,iEnsInv,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_openFile.',irc
        return
     end if
     if (bok) then
        fok(1)=fok(1)+1
     else
        frm(1)=frm(1)+1 ! no file
     end if
  end if
  if (bok) then
     call ncf_readInventory(iEnsInv,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_readInventory.',irc
        return
     end if
     if (bok) then
        fok(2)=fok(2)+1
     else
        frm(2)=frm(2)+1 ! no inventory
     end if
  end if
  if (bok) then
     ! get analysis time...
     iEnsFrtVar => ncf_getVariable(iEnsInv,"forecast_reference_time",bok,irc)
     if(irc.ne.0) then
        write(*,*)myname,'Error return from ncf_getVariable forecast_reference_time.',irc
        return
     end if
     if (bok) then
        fok(3)=fok(3)+1
     else
        frm(3)=frm(3)+1 ! no frt variable
     end if
  end if
  if (bok) then
     call ncf_readRealData(iEnsFrtVar,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_readData.',irc
        return
     end if
     if (bok) then
        fok(4)=fok(4)+1
     else
        frm(4)=frm(4)+1 ! no frt data
     end if
  end if
  if (bok) then
     iEnsTimeEntry=ncf_getDimEntry(iEnsInv,"time")
     if (iEnsTimeEntry.eq.0) bok=.false.
     if (bok) then
        fok(5)=fok(5)+1
     else
        frm(5)=frm(5)+1 ! no time dimension
     end if
  end if
  if (bok) then
     iEnsTimeVar => ncf_getVariable(iEnsInv,"time",bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_getVariable time.',irc
        return
     end if
     if (bok) then
        fok(6)=fok(6)+1
     else
        frm(6)=frm(6)+1 ! no time variable
     end if
  end if
  if (bok) then
     call ncf_readRealData(iEnsTimeVar,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_readData.',irc
        return
     end if
     if (bok) then
        fok(7)=fok(7)+1
     else
        frm(7)=frm(7)+1 ! no time data
     end if
  end if
  !
  ! read statistics file =========================
  !
  if (bok) then
     leno=length(out250,250,10)
     !
     ! open statistics file and read inventory
     !
     lens=length(stat250,250,10)
     write(*,*)myname,'Opening: ',stat250(1:lens)
     call ncf_openFile(stat250,iStaInv,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_openFile.',irc
        return
     end if
     !
     call ncf_readInventory(iStaInv,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_readInventory.',irc
        return
     end if
     !
     ! read all variable fields into memory (as real data)
     v => iStaInv%firstvariable%next
     do while (.not.associated(v,target=iStaInv%lastVariable))
        call ncf_readRealData(v,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_readData.',irc
           return
        end if
        if (.not.bok)then
           write(*,*)myname,"Unable to read forecast_reference_time."
           irc=924
           return
        end if
        v => v%next
     end do
     !
     ! identify analysis time...
     iStaFrtVar => ncf_getVariable(iStaInv,"forecast_reference_time",bok,irc)
     if(irc.ne.0) then
        write(*,*)myname,'Error return from getVariable.',irc
        return
     end if
     if (.not.bok)then
        write(*,*)myname,"Unable to find forecast_reference_time."
        irc=923
        return
     end if
     iStaFrtValue=iStaFrtVar%fd(1) ! latest analysis time
     !
     iStaStaEntry = ncf_getDimEntry(iStaInv,"statistics")
     if (iStaStaEntry.eq.0) then
        write(*,*)myname,'Unable to find statistics dimension.'
        irc=458
        return
     end if
     iStaStatDim => ncf_makeDimension(iStaInv,iStaStaEntry)
     !
     ! identify the resolution dimension...
     iStaResEntry = ncf_getDimEntry(iStaInv,"resolution")
     if (iStaResEntry.eq.0) then
        call ncf_printInventory(iStaInv)
        write(*,*)myname,'Unable to find resolution dimension.'
        irc=459
        return
     end if
     iStaResDim => ncf_makeDimension(iStaInv,iStaResEntry)
     resolution=iStaResDim%lim
     !
     ! identify the threshold dimension...
     iStaThrEntry = ncf_getDimEntry(iStaInv,"threshold")
     if (iStaThrEntry.eq.0) then
        call ncf_printInventory(iStaInv)
        write(*,*)myname,'Unable to find threshold dimension.'
        irc=460
        return
     end if
     iStaThrDim => ncf_makeDimension(iStaInv,iStaThrEntry)
      !
     iStaTimeEntry = ncf_getDimEntry(iStaInv,"time")
     iStaTimeVar => ncf_getVariable(iStaInv,"time",bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_getVariable.',irc
        return
     end if
     if (bok) then
        fok(8)=fok(8)+1
     else
        frm(8)=frm(8)+1 ! unable to read statistics file
     end if
  end if
  !
  ! loop over input file data and make output file ======================
  !
  if (bok) then
     ! check if we should make output file
     write(*,*)myname,'Creating output from: ',ens250(1:lene)
     
     if(bdeb)write(*,*)myname,'Firstivar:',iEnsInv%firstvariable%next%var250(1:iEnsInv%firstvariable%next%lenv)
     
     call ncf_openFile(ens250,oMtbeInv,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_openFile.',irc
        return
     end if

     call ncf_readInventory(oMtbeInv,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_copyInventory.',irc
        return
     end if

     
     ! loop over variables, reset variables with "time"
     oMtbeVar => oMtbeInv%firstvariable%next
     do while (.not.associated(oMtbeVar,oMtbeInv%lastVariable))
        vn => oMtbeVar%next
        if (ncf_variableContainsDim(oMtbeVar,iEnsTimeEntry).and..not.oMtbeVar%var250(1:oMtbeVar%lenv).eq."time") then
           ! remove ensemble dimension
           oMtbeDO => ncf_makeDimOrder(oMtbeVar,irc) ! get dimension order
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
              return
           end if
           iStaVar => ncf_getVariable(iStaInv, oMtbeVar%var250(1:oMtbeVar%lenv), bok, irc)
        if (.not.bok) irc=702
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_getVariable.',irc
              return
           end if
           a => ncf_getAttribute(iStaVar,"averaged_over",irc)
           if (associated(a)) then
              dim250=ncf_getAttributeText(a)
              iEnsAvgEntry=ncf_getDimEntry(oMtbeInv,dim250)
              if (iEnsAvgEntry.eq.0) then
                 lend=length(dim250,250,10)
                 write(*,*)myname,'Unable to find Ensemble dimension "'//dim250(1:lend)&
                      & //'" in Ens-file.'
                 irc=934
                 bok=.false.
                 return
              end if
              ! remove ensemble dimension
              iEnsDim => ncf_removeDimOrderEntry(oMtbeDO,iEnsAvgEntry)
           end if
           ! set new dimension order
           call ncf_setVariableDimOrder(oMtbeVar,oMtbeDO,irc) ! insert dimension order
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_setVariableDimOrder.',irc
              return
           end if
           call ncf_clearDimOrder(oMtbeDO)
           call ncf_initField(oMtbeVar,nf_fill_double,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_initField.',irc
              return
           end if
        else ! just read data from file
           call ncf_readRealData(oMtbeVar,bok,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_readData.',irc
              return
           end if
        end if
        oMtbeVar => vn
     end do

     oMtbeFrtVar => ncf_getVariable(oMtbeInv,"forecast_reference_time",bok,irc)
        if (.not.bok) irc=702
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_getVariable.',irc
        return
     end if
     oMtbeFrtValue=oMtbeFrtVar%fd(1)
     oMtbeTimeEntry = ncf_getDimEntry(oMtbeInv,"time")
     oMtbeTimeVar => ncf_getVariable(oMtbeInv,"time",bok,irc)
        if (.not.bok) irc=702
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_getVariable.',irc
        return
     end if
     oMtbeTimeDim => ncf_getVariableDimension( oMtbeTimeVar, oMtbeTimeEntry )
     !
     call ncf_closeFile(oMtbeInv,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_closeFile.',irc
        return
     end if
     !
     firstOut=.false.
     if (bok) then
        fok(9)=fok(9)+1
     else
        frm(9)=frm(9)+1 ! unable to make output file
     end if
  end if
     !
  if (bok) then ! process variable
     write(*,*)myname,"Processing:",iEnsInv%fn250(1:iEnsInv%lenf)
     firstivar=.true.
     iEnsVar => iEnsInv%firstvariable%next
     IVARLOOP: do while (.not.associated(iEnsVar,iEnsInv%lastVariable))
        ! loop over variables in ensemble-file
        ecnt=0
        bok=.true. ! process variable
        if (bdeb) write(*,*)myname,'I-var-loop start: "'//iEnsVar%var250(1:iEnsVar%lenv)//'"'
        !
        iEnsDo => ncf_makeDimOrder(iEnsVar,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
           return
        end if
        call ncf_clearDimOrder(iEnsTimeDO)
        iEnsTimeDO => ncf_newDimOrder(iEnsInv,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_newDimOrder.',irc
           return
        end if
        iEnsTimeDim => ncf_removeDimOrderEntry(iEnsDo,iEnsTimeEntry)
        if (associated(iEnsTimeDim)) then
           call ncf_addDimOrderDim(iEnsTimeDO,iEnsTimeDim,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_addDimOrderDim.',irc
              return
           end if
           ook(1)=ook(1)+1
        else
           orm(1)=orm(1)+1 ! no time dimension
           bok=.false. ! no time dimension...
        end if
        !
        if (bdeb) write(*,*)myname,'Checking for averaging.',bok
        if (bok) then
           iStaVar => ncf_getVariable(iStaInv, iEnsVar%var250(1:iEnsVar%lenv), bok, irc)
           if (.not.bok) irc=702
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_getVariable.',irc
              return
           end if
           iStado => ncf_makeDimOrder(iStaVar,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
              return
           end if
           !call ncf_printDimOrder(istado)
           if (bdeb) write(*,*)myname,'...remove time.',associated(istado),istatimeentry
           iStaTimeDim => ncf_removeDimOrderEntry(iStaDo,iStaTimeEntry)
           if (bok) then
              ook(2)=ook(2)+1
           else
              orm(2)=orm(2)+1 ! no statistics variable
           end if
        end if
        if (bok) then
           iEnsAvgDO => ncf_newDimOrder(iEnsInv,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_newDimOrder.',irc
              return
           end if
           a => ncf_getAttribute(iStaVar,"averaged_over",irc)
           if (associated(a)) then
              dim250=ncf_getAttributeText(a)
              lend=length(dim250,250,10)
              iEnsAvgEntry=ncf_getDimEntry(iEnsInv,dim250(1:lend))
              if (iEnsAvgEntry.eq.0) then
                 write(*,*)myname,'Unable to find Ensemble dimension "'//dim250(1:lend)&
                      & //'" in Ens-file.'
                 irc=934
                 bok=.false.
                 return
              end if
              iEnsAvgDim => ncf_removeDimOrderEntry(iEnsDo,iEnsAvgEntry)
              ! the average dimension must exist in the statistics file...
              ! ...it will automatically be deleted if not used.
              iStaAvgEntry=ncf_getDimEntry(iStaInv,dim250(1:lend))
              if (iStaAvgEntry.eq.0) then
                 iStaAvgDim => ncf_createDimension(iStaInv,dim250(1:lend),1)
              end if
              if (associated(iEnsAvgDim)) then
                 call ncf_addDimOrderDim(iEnsAvgDO,iEnsAvgDim,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_addDimOrderDim.',irc
                    return
                 end if
              else
                 bok=.false. ! no average dimension in ensemble variable
              end if
           end if
           if (bok) then
              ook(3)=ook(3)+1
           else
              orm(3)=orm(3)+1 ! no average dimension
           end if
        end if
        !
        if (bdeb) write(*,*)myname,'Find output variable.',bok
        if (bok) then
           ! find corresponding variable in output-file
           oMtbeVar => ncf_getVariable(oMtbeInv, iEnsVar%var250(1:iEnsVar%lenv),bok,irc)
           if (.not.bok) irc=702
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_getVariable.',irc
              return
           end if
           if (bdeb) write(*,*)myname,'...make dimension order.'
           oMtbedo => ncf_makeDimOrder(oMtbeVar,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
              return
           end if
           if (bdeb) write(*,*)myname,'...remove time.'
           oMtbeTimeDim => ncf_removeDimOrderEntry(oMtbedo,oMtbeTimeEntry)
           !
           if (bok) then
              ook(4)=ook(4)+1
           else
              orm(4)=orm(4)+1 ! no output variable
           end if
        end if
        if (bok) then
           if (ncf_variableContainsDim(iEnsVar,iEnsTimeEntry).and.&
                & .not.iEnsVar%var250(1:iEnsVar%lenv).eq."time") then
              call ncf_readRealData(iEnsVar,bok,irc)
              if (irc.ne.0) then
                 write(*,*)myname,'Error return from ncf_readData.',irc
                 return
              end if
              if (bok) then
                 ook(6)=ook(6)+1
              else
                 orm(6)=orm(6)+1 ! unable to read time data
              end if
              if (bok) then
                 ! find corresponding distribution variable in statistics-file
                 iStaThrVar => ncf_getVariable(iStaInv, iEnsVar%var250(1:iEnsVar%lenv)//sthr, bok, irc)
                 if (.not.bok) irc=702
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_getVariable.',irc
                    return
                 end if
                 if (bok) then
                    ook(7)=ook(7)+1
                 else
                    orm(7)=orm(7)+1 ! no statistics threshold variable
                 end if
              end if
              if (bok) then
                 range=ncf_getRealAttribute(iStaThrVar,"range")
                 if (bdeb) write(*,*)myname,' Range:',range,resolution,&
                      & iStaThrVar%var250(1:iStaThrVar%lenv)
                 resfact = real(resolution)/range
                 !
                 ! change output variable name
                 !
                 oMtbeVar%var250=iEnsVar%var250(1:iEnsVar%lenv)//"_mtbe"
                 call chop0(oMtbeVar%var250,250)
                 oMtbeVar%lenv=length(oMtbeVar%var250,250,10)
                 !
                 ! reset "units" attribute...
                 !
                 a => ncf_getAttribute(oMtbeVar,"units",irc)
                 if (associated(a)) then
                    call ncf_setAttributeText(a,"1")
                 end if
                 !
                 ! change standard name
                 !
                 a => ncf_getAttribute(oMtbeVar,"standard_name",irc)
                 if (associated(a)) then
                    nam250=ncf_getAttributeText(a)
                    call chop0(nam250,250)
                    lenn=length(nam250,250,10)
                    call ncf_setAttributeText(a,nam250(1:lenn)//"_mtbe")
                    call ncf_setTextAttribute(oMtbeVar,"long_name",nam250(1:lenn)//"_mtb")
                    call ncf_setTextAttribute(oMtbeVar,"description","Mean time between grid events.")
                    call ncf_setTextAttribute(oMtbeVar,"units","years")
                 end if
                 !
                 scnt=0.0D0
                 smean=0.0D0
                 sstdv=0.0D0
                 !
                 xcnt=0.0D0
                 xmean=0.0D0
                 xstdv=0.0D0
                 !
                 ucnt=0.0D0
                 umean=0.0D0
                 ustdv=0.0D0
                 !
                 do ii=0,nt
                    utcnt(ii)=0.0D0
                    ltcnt(ii)=0.0D0
                 end do
                 !
                 call ncf_makeIndex(io,iEnsInv,oMtbeInv,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_makeIndex (io).',irc
                    return
                 end if
                 call ncf_makeIndex(is,iEnsInv,iStaInv,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_makeIndex (is).',irc
                    return
                 end if
                 ITIMELOOP: do ik=1,iEnsTimeVar%lend ! loop over time
                    if (bdeb) write(*,*)myname,'Time loop : ',ik,iEnsTimeVar%lend
                    !
                    jk=iEnsTimeVar%lend
                    JTIMELOOP: do while (jk.ge.1)! find j-time index
                       if (mod(iEnsTimeVar%fd(jk)-iEnsTimeVar%fd(ik),86400.0D0).eq.0) then
                          exit JTIMELOOP
                       end if
                       jk=jk-1
                    end do JTIMELOOP
                    bbok=(jk.ge.1) ! did we find a matching time?
                    if (bbok) then
                       ook(8)=ook(8)+1
                    else
                       orm(8)=orm(8)+1 ! no climatology found
                    end if
                    if (bbok) then
                       minv=0.0D0
                       maxv=0.0D0
                       !
                       call ncf_resetPos(iEnsDo,irc)
                       !
                       call ncf_setDimensionValue(iEnsInv,iEnsTimeDim,ik)
                       call ncf_setDimensionValue(oMtbeInv,oMtbeTimeDim,ik)
                       !
                       LOC: do while (ncf_increment(iEnsInv,iEnsDo,irc)) ! loop over locations
                          call ncf_copyIndexpos(io,iEnsInv,oMtbeInv,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from ncf_copyIndexPos (io).',irc
                             return
                          end if
                          call ncf_copyIndexpos(is,iEnsInv,iStaInv,irc)
                          if (irc.ne.0) then
                             write(*,*)myname,'Error return from ncf_copyIndexPos (is).',irc
                             return
                          end if
                          !
                          lo = ncf_getLocation(oMtbeVar) ! location in output array
                          !
                          call ncf_setDimensionValue(iStaInv,iStaTimeDim,ik)  ! actual time
                          call ncf_setDimensionValue(iStaInv,iStaStatDim,1) ! 1: cnt
                          lcnt = ncf_getLocation(iStaVar)
                          call ncf_setDimensionValue(iStaInv,iStaStatDim,2) ! 2: mean
                          lmean = ncf_getLocation(iStaVar)
                          call ncf_setDimensionValue(iStaInv,iStaStatDim,3) ! 3: stdv
                          lstdv = ncf_getLocation(iStaVar)
                          call ncf_setDimensionValue(iStaInv,iStaStatDim,4) ! 4: stdv residual
                          lres = ncf_getLocation(iStaVar)
                          ! locations in statistics array
                          call ncf_setDimensionValue(iStaInv,iStaStatDim,5) ! 5: time
                          lfrt = ncf_getLocation(iStaVar)
                          call ncf_setDimensionValue(iStaInv,iStaStatDim,6) ! 6: aa
                          laa = ncf_getLocation(iStaVar)
                          call ncf_setDimensionValue(iStaInv,iStaStatDim,7) ! 7: bb
                          lbb = ncf_getLocation(iStaVar)
                          call ncf_setDimensionValue(iStaInv,iStaStatDim,8) ! 8: cc
                          lcc = ncf_getLocation(iStaVar)
                          call ncf_setDimensionValue(iStaInv,iStaStatDim,9) ! 9: dd
                          ldd = ncf_getLocation(iStaVar)
                          !
                          if (iStaVar%fd(lcnt).ne.nf_fill_double.and.&
                               & iStaVar%fd(lstdv).ne.nf_fill_double.and.&
                               & iStaVar%fd(lres).ne.nf_fill_double.and.&
                               & iStaVar%fd(laa).ne.nf_fill_double.and.&
                               & iStaVar%fd(lbb).ne.nf_fill_double.and.&
                               & iStaVar%fd(lcc).ne.nf_fill_double.and.&
                               & iStaVar%fd(ldd).ne.nf_fill_double) then
                             cmean = iStaVar%fd(lmean)
                             cstdv = iStaVar%fd(lstdv)
                             if (cstdv.ne.nf_fill_double.and.&
                                  & cstdv.ne.nf_fill_double) then
                                !call ncf_printVarPos(iExpVar)
                                scnt=scnt+1.0D0
                                smean=smean+cmean
                                sstdv=sstdv+cstdv
                             end if
                             !
                             vals=0.0D0
                             cnts=0.0D0
                             !
                             caa=iStaVar%fd(laa)
                             cbb=iStaVar%fd(lbb)
                             ccc=iStaVar%fd(lcc)
                             cdd=iStaVar%fd(ldd)
                             !
                             bbok=.true.
                             call ncf_resetPos(iEnsAvgDo,irc)
                             ENS: do while (ncf_increment(iEnsInv,iEnsAvgDo,irc)) ! loop over ensembles
                                li = ncf_getLocation(iEnsVar) ! location in env array
                                val = iEnsVar%fd(li) ! input value
                                !
                                if (val.ne.nf_fill_double) then ! we have a real value
                                   !
                                   xval = (val-cmean)/max(1.0D-10,cstdv)
                                   valu= getu(caa,cbb,ccc,cdd,xval,bok,irc)
                                   if (irc.ne.0) then
                                      write(*,*)myname,'Error return from getu.',irc
                                      return
                                   end if
                                   if (bok) then
                                      xcnt=xcnt+1.0D0
                                      xmean=xmean+xval
                                      xstdv=xstdv+xval*xval
                                      ucnt=ucnt+1.0D0
                                      umean=umean+valu
                                      ustdv=ustdv+valu*valu
                                      rpos=max(2.0,min(real(resolution),real(resolution)/2.0D0 + &
                                           & resfact*valu   -  0.5D0)) ! accumulated values are shifted by 0.5
                                      call ncf_setDimensionValue(iStaInv,iStaTimeDim,jk)  ! climatology time
                                      call ncf_setDimensionValue(iStaInv,iStaThrDim,2)  ! reference
                                      call ncf_setDimensionValue(iStaInv,iStaResDim,floor(rpos))     ! probability floor
                                      lf = ncf_getLocation(iStaThrVar)
                                      !call ncf_printVarPos(iStaThrVar)
                                      call ncf_setDimensionValue(iStaInv,iStaResDim,ceiling(rpos))     ! probability ceiling
                                      lc = ncf_getLocation(iStaThrVar)
                                      if (iStaThrVar%fd(lf).ne.nf_fill_double.and.&
                                           & iStaThrVar%fd(lc).ne.nf_fill_double) then
                                         prob=iStaThrVar%fd(lf) + (rpos-floor(rpos))*(iStaThrVar%fd(lc)-iStaThrVar%fd(lf))
                                         if (prob.gt.0.5D0) then
                                            vald=(0.5D0+gamma)/(1.0D0-prob+gamma)-1.0D0
                                         else
                                            vald=-(0.5D0+gamma)/(prob+gamma)+1.0D0
                                         end if
                                         valy = max(min(199.0D3,vald / 365.0D0),-199.0D3) ! days -> years
                                         minv=min(minv,valy)
                                         maxv=max(maxv,valy)
                                         ! if ((abs(valy-minv).lt.1.0D-10.or.abs(valy-maxv).lt.1.0D-10)&
                                         !      & .and.abs(val).gt.10.0D0.and. ecnt.le.2) then
                                         !    ecnt=ecnt+1
                                         !    write(*,*)myname,'Extreme "'//iEnsVar%var250(1:iEnsVar%lenv)//'" ',&
                                         !         & valy,val
                                         ! end if
                                         !
                                         utcnt(0)=utcnt(0)+1.0D0
                                         ltcnt(0)=ltcnt(0)+1.0D0
                                         do ii=1,nt
                                            if (valy.ge.thres(ii)) then
                                               utcnt(ii)=utcnt(ii)+1.0D0
                                            else if (valy.le.-thres(ii)) then
                                               ltcnt(ii)=ltcnt(ii)+1.0D0
                                            end if
                                         end do
                                         ! calculate ensemble mean MTBE
                                         vals=vals+valy
                                         !vals=vals+xval
                                         !vals=vals+max(min(500.0D0,valy/365.0D0),-500.0D0)
                                         cnts=cnts+1.0D0

                                         ! if (abs(valy).gt.1.0D3) then
                                         !    write(*,*)myname,'Extreme value:',valy,vals/cnts,vald,prob,valu
                                         ! end if
                                         
                                         ! if (ook(11).lt.2) then
                                         ! write(*,*)myname,'Trans:',caa,cbb,ccc,cdd,cmean,cstdv
                                         ! !& ,xval,caa+cbb*valu+ccc*valu*valu+cdd*valu*valu*valu
                                         ! write(*,*)myname,'Val:',val,xval,valu,valy
                                         ! write(*,*)myname,'Prob:',iStaThrVar%fd(lf),iStaThrVar%fd(lc),prob,vald
                                         ! write(*,*)myname,'Stat:',cnts,vals
                                         ! else
                                         !    stop ("Debug")
                                         ! end if

                                      end if
                                      ook(12)=ook(12)+1
                                   else
                                      orm(12)=orm(12)+1 ! unable to get transformed variable
                                      bbok=.false.
                                   end if
                                   ook(11)=ook(11)+1
                                else
                                   orm(11)=orm(11)+1 ! undefined ensemble value
                                end if
                             end do ENS
                             if (cnts.gt.0.0D0) then
                                valm = vals/max(1.0D0,cnts)
                             else
                                valm=nf_fill_double
                             end if
                             
                             ! if (ook(10).lt.20) then
                             !    write(*,*)myname,'=======Val:',valm,cnts
                             ! else
                             !    write(*,*)myname,'=======Val:',valm,cnts
                             !    stop ("Debug")
                             ! end if

                             if (bbok) then
                                ook(10)=ook(10)+1
                             else
                                orm(10)=orm(10)+1 ! no valid ensemble data
                             end if
                             ook(9)=ook(9)+1
                          else
                             orm(9)=orm(9)+1 ! no transformation available for location
                             valm=nf_fill_double
                          end if
                          !write(*,*)myname,'ENS:',valm,vals,cnts

                          !write(*,*)myname,'Value:', iEnsVar%var250(1:iEnsVar%lenv),lo,val,cstdv,cnt
                          !if (abs(valm).lt.1.0D-5) then
                          !   write(*,*)myname,'Valm:',valm,vals,cnts
                          !end if
                          oMtbeVar%fd(lo) = valm
                       end do LOC
                       if (lfldat(90)) then ! dump data
                          if (jk.eq.ik) then
                             !write(*,*) iEnsVar%varid,iEnsVar%var250(1:iEnsVar%lenv),ik
                             call ncf_setDimensionValue(iStaInv,iStaTimeDim,jk)  ! climatology time
                             call ncf_setDimensionValue(iStaInv,iStaThrDim,2)  ! reference
                             do ii=1,resolution
                                call ncf_setDimensionValue(iStaInv,iStaResDim,ii)     ! stdv
                                lf = ncf_getLocation(iStaThrVar)
                                prob=iStaThrVar%fd(lf)
                                if (prob.gt.0.5D0) then
                                   val=1.0D0/max(1.0D-8,1.0D0-prob)-2.0D0
                                   lgg=log10(max(val+1.0,1.0D-8))
                                else
                                   val=-1.0D0/max(1.0D-8,prob)+2.0D0
                                   lgg=-log10(max(-val+1.0,1.0D-8))
                                end if
                                val = max(min(199.0D3,val / 365.0D0),-199.0D3) ! days -> years
                                
                                if (abs(val).lt.1.0D5) then
                                   write(15,'(2(I3,X),I5,X,F10.3,2(X,F15.8))') &
                                        & iEnsVar%varid,ik,ii,val,prob,lgg
                                end if
                             end do
                          end if
                       end if
                    end if
                 end do ITIMELOOP
                 if (bdeb) write(*,*)myname,'Time loop stop',scnt
                 smean=smean/max(scnt,1.0D0)
                 !write(*,*) myname,'Stdv:',sstdv,scnt,cstdv
                 sstdv=sstdv/max(scnt,1.0D0)
                 xmean=xmean/max(xcnt,1.0D0)
                 xstdv=xstdv/max(xcnt,1.0D0)
                 xstdv=dsqrt(max(0.0D0,xstdv-xmean*xmean))
                 umean=umean/max(ucnt,1.0D0)
                 ustdv=ustdv/max(ucnt,1.0D0)
                 ustdv=dsqrt(max(0.0D0,ustdv-umean*umean))
                 write(*,'(X,A,2(A,F10.2),A,I0,A)') iEnsVar%var250(1:iEnsVar%lenv),', mean:',&
                      & smean,', stdv:',sstdv,', cnt:',nint(scnt),' (raw)'
                 write(*,'(X,A,2(A,F10.5),A,I0,A)') iEnsVar%var250(1:iEnsVar%lenv),', mean:',&
                      & xmean,', stdv:',xstdv,', cnt:',nint(xcnt),' (corrected)'
                 write(*,'(X,A,2(A,F10.5),A,I0,A)') iEnsVar%var250(1:iEnsVar%lenv),', mean:',&
                      & umean,', stdv:',ustdv,', cnt:',nint(ucnt),' (transformed)'
                 !
                 write(*,'(" +++++ ",10(X,I12," (",F7.3,"%)"))') &
                      & (nint(utcnt(ii)), utcnt(ii)*100.0D0/max(1.0D0,utcnt(0)),ii=1,nt)
                 write(*,'(" ----- ",10(X,I12," (",F7.3,"%)"))') &
                      & (nint(ltcnt(ii)), ltcnt(ii)*100.0D0/max(1.0D0,ltcnt(0)),ii=1,nt)
                 !
                 if (allocated(io)) deallocate(io)
                 if (allocated(is)) deallocate(is)
                 call ncf_clearVariable(iEnsVar,irc)
                 if(irc.ne.0)return
              end if
              ook(5)=ook(5)+1
           else
              orm(5)=orm(5)+1 ! variable is the time variable
           end if
           call ncf_clearDimOrder(iEnsDo)
           call ncf_clearDimOrder(iEnsTimeDo)
           call ncf_clearDimOrder(iEnsAvgDO)
           call ncf_clearDimOrder(oMtbeDO)
           call ncf_clearDimOrder(iStaDO)
           if (bok) then
              ook(6)=ook(6)+1
           else
              orm(6)=orm(6)+1
           end if
           if (bdeb) write(*,*)myname,'I-var-loop end'
        end if
        iEnsVar => iEnsVar%next
     end do IVARLOOP
     ! loop over output variables and prepare writing to file...
     if(bdeb)write(*,*) myname,'End of I-VAR.'
     if (bok) then
        fok(10)=fok(10)+1
     else
        frm(10)=frm(10)+1 ! unable to process
     end if
  end if
  !
  if (bok) then
     call ncf_closeFile(iEnsInv,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_closeFile.',irc
        return
     end if
     !
     !     call ncf_checkInventory(iEnsInv,irc)
     !
     lene=length(ens250,250,10)
     write(*,*)myname,"Clearing:",ens250(1:lene)
     call ncf_clearInventory(iEnsInv,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_clearInventory.',irc
        return
     end if
     if(bdeb)write(*,*) myname,'End of I-FILE.'
  end if
  !
  if(bdeb)write(*,*) myname,'Done processing.'
  !
  if (bok) then
     if (.not.firstout) then
        if(bdeb)write(*,*) myname,'Compressing variables.'
        call ncf_compressVariables(oMtbeInv,irc) ! compress variables
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_compressVariables.',irc
           return
        end if
        if(bdeb)write(*,*) myname,'Compressing dimensions.'
        call ncf_compressDimensions(oMtbeInv) ! compress dimensions
        !
        if(bdeb)write(*,*) myname,'Writing output.'
        oMtbeInv%fn250=out250
        oMtbeInv%lenf=length(out250,250,10)
        !     call ncf_checkInventory(oMtbeInv,irc)
        call ncf_writeNcOut(oMtbeInv,out250,irc) ! open, write, close file
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_writeNcOut.',irc
           return
        end if
        call ncf_clearInventory(oMtbeInv,irc) ! delete internally stored data...
        if (irc.ne.0) then
           write(*,*) myname,"Error return from ncf_clearInventory.",irc
           return
        end if
     else
        write(*,*)myname,'No data available for output file.'
     end if
  end if
  if (associated(iEnsInv)) deallocate(iEnsInv)
  if (associated(oMtbeInv)) deallocate(oMtbeInv)
  if (associated(iStaInv)) deallocate(iStaInv)
  !
  ! write statistics
  !
  if(bdeb)write(*,*) myname,'Writing statistics.'
  do ii=1,10
     pst(ii)=dfloat(frm(ii))/max(1.0d0,dfloat(frm(ii)+fok(ii)))*100
  end do
  WRITE(*,*)
  WRITE(*,998) MYNAME,                 'Available files:           ', fok(1)+frm(1)
  IF (FRM(1).NE.0) WRITE(*,999) MYNAME,'Unable to read file:       ', -FRM(1),PST(1)
  IF (FRM(2).NE.0) WRITE(*,999) MYNAME,'Invalid inventory:         ', -FRM(2),PST(2)
  IF (FRM(3).NE.0) WRITE(*,999) MYNAME,'Missing FRT-variable:      ', -FRM(3),PST(3)
  IF (FRM(4).NE.0) WRITE(*,999) MYNAME,'Missing FRT-value:         ', -FRM(4),PST(4)
  IF (FRM(5).NE.0) WRITE(*,999) MYNAME,'Missing time-dimension:    ', -FRM(5),PST(5)
  IF (FRM(6).NE.0) WRITE(*,999) MYNAME,'Missing time-variable:     ', -FRM(6),PST(6)
  IF (FRM(7).NE.0) WRITE(*,999) MYNAME,'Missing time-value:        ', -FRM(7),PST(7)
  IF (FRM(8).NE.0) WRITE(*,999) MYNAME,'Invalid statistics file:   ', -FRM(8),PST(8)
  IF (FRM(9).NE.0) WRITE(*,999) MYNAME,'Unable to make output file:', -FRM(9),PST(9)
  IF (FRM(10).NE.0) WRITE(*,999) MYNAME,'Unable to process data:    ', -FRM(10),PST(10)
  WRITE(*,*)   MYNAME,     '-----------------------------------------------'
  WRITE(*,998) MYNAME,                 'Accepted files:            ', fok(10)
  write(*,*)
  !
  do ii=1,20
     pst(ii)=dfloat(orm(ii))/max(1.0d0,dfloat(orm(ii)+ook(ii)))*100
  end do
  WRITE(*,998) MYNAME,                 'Available variables        ', ook(1)+orm(1)
  IF (ORM(1).NE.0) WRITE(*,999) MYNAME,'No time dimension:         ', -ORM(1),PST(1)
  IF (ORM(2).NE.0) WRITE(*,999) MYNAME,'No statistics variable:    ', -ORM(2),PST(2)
  IF (ORM(3).NE.0) WRITE(*,999) MYNAME,'No average dimension:      ', -ORM(3),PST(3)
  IF (ORM(4).NE.0) WRITE(*,999) MYNAME,'No output variable         ', -ORM(4),PST(4)
  IF (ORM(5).NE.0) WRITE(*,999) MYNAME,'System variables:          ', -ORM(5),PST(5)
  IF (ORM(6).NE.0) WRITE(*,999) MYNAME,'Unable to read data:       ', -ORM(6),PST(6)
  IF (ORM(7).NE.0) WRITE(*,999) MYNAME,'No threshold statistics:   ', -ORM(7),PST(7)
  WRITE(*,*)   MYNAME,     '-----------------------------------------------'
  WRITE(*,999) MYNAME,                 'Accepted variables:        ', ook(7),&
       & dfloat(ook(7))/max(1.0d0,dfloat(orm(1)+ook(1)))*100.0D0
  write(*,*)
  IF (ORM(8).NE.0) then
     WRITE(*,999) MYNAME,'No climatology time found: ', -ORM(8),PST(8)
     write(*,*)
  end if
  WRITE(*,998) MYNAME,                 'Available data             ', ook(9)+orm(9)
  IF (ORM(9).NE.0) WRITE(*,999) MYNAME,'No valid var-trans:        ', -ORM(9),PST(9)
  IF (ORM(10).NE.0)WRITE(*,999) MYNAME,'No valid ensemble-trans:   ', -ORM(10),PST(10)
  WRITE(*,*)   MYNAME,     '-----------------------------------------------'
  WRITE(*,999) MYNAME,                 'Accepted data:             ', ook(10),&
       & dfloat(ook(10))/max(1.0d0,dfloat(orm(9)+ook(9)))*100.0D0
  write(*,*)
  WRITE(*,998) MYNAME,                 'Available ensembles        ', ook(11)+orm(11)
  IF (ORM(11).NE.0)WRITE(*,999) MYNAME,'No ensemble value:         ', -ORM(11),PST(11)
  IF (ORM(12).NE.0)WRITE(*,999) MYNAME,'No ensemble var-trans:     ', -ORM(12),PST(12)
  WRITE(*,*)   MYNAME,     '-----------------------------------------------'
  WRITE(*,999) MYNAME,                 'Accepted ensembles:        ', ook(12),&
       & dfloat(ook(12))/max(1.0d0,dfloat(orm(11)+ook(11)))*100.0D0
  !
999 FORMAT(X,A12,X,A27,I10,' (',F5.1,'%)')
998 FORMAT(X,A12,X,A27,I10)
997 FORMAT(X,A12,X,A,10(X,I10))
  !
contains
  real function phi(x) ! cumulative normal distribution
    real :: x, t, y, s
    ! constants
    real :: a1 =  0.254829592
    real :: a2 = -0.284496736
    real :: a3 =  1.421413741
    real :: a4 = -1.453152027
    real :: a5 =  1.061405429
    real :: p  =  0.3275911
    ! Save the sign of x
    s = 1.0D0
    if (x < 0) s = -1.0D0
    x = abs(x)/sqrt(2.0)
    ! A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
    phi =  0.5*(1.0 + s*y)
    return
  end function phi
end subroutine mncmtbe

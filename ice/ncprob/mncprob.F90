! ecmwf netcdf filer... 
SUBROUTINE MNCPROB(UNITI,IRC)
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
  DATA MYNAME /'MNCPROB'/
  ! 
  LOGICAL  LFLDAT(100)
  logical :: bdeb=.false.
  ! 
  INTEGER  lend, lenh, leni, leno, lens, lenr
  integer, external :: length
  ! 
  INTEGER  KODE, LINE
  INTEGER  NRHDR
  PARAMETER (NRHDR=100)
  CHARACTER*250 HDR250(NRHDR),DAT250, INF250
  CHARACTER*250 FILENM,NUKEHEAD
  EXTERNAL NUKEHEAD
  INTEGER  INTOUT
  LOGICAL  ENDOFF
  !
  integer, parameter :: resolution = 1000
  integer :: nerr = 0
  integer :: tcnt,icnt
  real, parameter :: range = 200.0D0
  real :: rescnt(0:resolution),resacc
  real, parameter :: resfact = real(resolution)/range
  character*4, parameter :: sthr="_thr"
  character*4, parameter :: smom="_mom"
  REAL :: MINV, MAXV
  logical :: firstv
  !
  type(inventory),pointer  :: iEnsInv => null()
  type(inventory),pointer  :: oStaInv => null()
  type(variable), pointer  :: oStaVar,oStaFrtVar,oStaTimeVar,oStaThrVar,oStaMomVar
  type(variable), pointer  :: iEnsVar,iEnsFrtVar,iEnsTimeVar
  type(dimension), pointer :: oStaTimeDim,oStaResDim,oStaThrDim,oStaAvgDim,oStaDim,oStaMomDim
  type(dimension), pointer :: iEnsTimeDim, iEnsDim
  type(dimensionOrder),pointer  ::  oStaDO, oStaThrDO,oStaMomDO
  type(dimensionOrder), pointer :: iEnsDO, iEnsDimDO
  integer :: oStaResEntry,oStaThrEntry,oStaAvgEntry,oStaMomEntry,oStaTimeEntry
  integer :: iEnsTimeEntry,iEnsEntry,oStaEntry
  integer :: cvar,ivar,ix,iy
  real :: oStaThrFrtValue,oStaFrtValue
  logical :: bok,climOk,transOk,bdone,firstoSta,firstiens
  !
  ! 
  integer ook(10), orm(10)
  real pst(10)
  !
  INTEGER ret
  ! 
  character*250 :: iens250, ista250, osta250, dim250
  integer, parameter :: ne=4
  integer :: ii, jj, ik, jk, mm, mx, my, pos, li
  integer :: lcnt,lmean,lstdv,lres,lfrt,laa,lbb,lcc,ldd,lwgt,ll(ne)
  integer :: lxy(0:ne,0:ne),la,lb
  integer :: clcnt,clmean,clstdv,clres,clfrt,claa,clbb,clcc,cldd,clwgt
  integer :: clxy(0:ne,0:ne),cla,clb
  real :: cmean,cstdv
  real :: scnt,sxy(0:ne,0:ne), sx(0:ne), sy(0:ne)
  logical :: uniq
  integer, allocatable :: io(:)
  real :: dtrelax, dtdist
  real :: cnt
  real :: ym,cres
  real :: val,xval,valu,xx,dx,dt,dtx
  real :: xtx(2,2),xty(2)
  real :: totalWgt
  real :: oldWgt,wgt0,wgt1,twgt
  integer :: ipos,lt1,lt2,clt1,clt2
  real :: rpos, det, detinv, B(2,2),A(2)
  real, parameter :: dpdt=1.0D0/(86400.0D0)
  integer, parameter :: nabc=4
  real :: ee(ne),ff(ne),xabc(nabc),mu,gm
  real, external :: getu

  !
  IRC=0
  ! 
  if(bdeb)write(*,*) myname,'Debug: Routine starts.',irc
  !
  ook(:)=0
  orm(:)=0
  ! 
  dtrelax=10.0D0
  dtdist=20.0D0

  DO II=1,NRHDR
     HDR250(II) = ''
     LFLDAT(II)=.FALSE.
  ENDDO

  HDR250(10)  = 'NCPROB V1.0 [0]VFLR'
  HDR250(20)  = 'INPUT ENSEMBLE FILE (NETCDF) [1]VFLR &'
  HDR250(30)  = 'AVERAGE OVER DIMENSION [1]VFLR &'
  HDR250(50)  = 'INPUT STATISTICS FILE (NETCDF) [1]VFLR'
  HDR250(60)  = 'OUTPUT STATISTICS FILE (NETCDF) [1]VFLR &'
  HDR250(80)  = 'STATISTICAL RELAXATION TIMES (DAYS) [1]VFLR'
  HDR250(90)  = 'FORCE DATA INTO STATISTICS FILE[0]'
  HDR250(91)  = 'INCLUDE MOMENTS IN STATISTICS FILE[0]'
  ! 
  ! NF_FILL_CHAR, NF_FILL_BYTE, NF_FILL_SHORT, NF_FILL_INT, NF_FILL_FLOAT, and NF_FILL_DOUBLE
  ! 
  WRITE(*,*) MYNAME,'----------------------------------------'
  ! 
  ! READ DATA FROM INPUT FILES..............................
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
     IF (LINE.EQ.10) THEN
     ELSEIF (LINE.EQ.20) THEN ! input ensemble file
        iens250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.30) THEN ! average over dimension
        dim250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.50) THEN ! statistics file netcdf
        ista250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.60) THEN ! statistics file netcdf
        osta250=dat250
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.80) THEN ! statistics file netcdf
        read(dat250(1:lend),*, err=201) dtrelax, dtdist
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.90) THEN ! force data
        LFLDAT(LINE) = .TRUE.
     ELSEIF (LINE.EQ.91) THEN ! include moments
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
  ! open statistics file and read inventory
  !
  if (lfldat(50)) then
     lens=length(ista250,250,10)
     write(*,*)myname,'Opening: ',ista250(1:lens)
     call ncf_openFile(ista250,oStaInv,bok,irc)
     if (irc.ne.0) then
        irc=0
        firstoSta=.true.
     else
        firstoSta=.false.
        bok=.true.
        !
        call ncf_readInventory(oStaInv,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_readInventory.',irc
           return
        end if
        !
        ! read all variable fields into memory (as real data)
        oStaVar => oStaInv%firstvariable%next
        do while (.not.associated(oStaVar,target=oStaInv%lastVariable))
           call ncf_readRealData(oStaVar,bok,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_readData.',irc
              return
           end if
           if (.not.bok)then
              write(*,*)myname,"Unable to read forecast_reference_time."
              irc=924
              return
           end if
           oStaVar => oStaVar%next
        end do
        !
        ! identify the statistics dimension...
        oStaEntry = ncf_getDimEntry(oStaInv,"statistics") ! a2, ab, b2, a, b, cnt
        if (oStaEntry.eq.0) then
           call ncf_printInventory(oStaInv)
           write(*,*)myname,'Unable to find statistics dimension.'
           irc=458
           return
        end if
        oStaDim => ncf_makeDimension(oStaInv,oStaEntry) ! a2, ab, b2, a, b, cnt
        !
        ! identify the resolution dimension...
        oStaResEntry = ncf_getDimEntry(oStaInv,"resolution") ! a2, ab, b2, a, b, cnt
        if (oStaResEntry.eq.0) then
           call ncf_printInventory(oStaInv)
           write(*,*)myname,'Unable to find resolution dimension.'
           irc=459
           return
        end if
        oStaResDim => ncf_makeDimension(oStaInv,oStaResEntry) ! a2, ab, b2, a, b, cnt
        !
        ! identify the threshold dimension...
        oStaThrEntry = ncf_getDimEntry(oStaInv,"threshold") ! a2, ab, b2, a, b, cnt
        if (oStaThrEntry.eq.0) then
           call ncf_printInventory(oStaInv)
           write(*,*)myname,'Unable to find threshold dimension.'
           irc=460
           return
        end if
        oStaThrDim => ncf_makeDimension(oStaInv,oStaThrEntry) ! a2, ab, b2, a, b, cnt
        !
        ! identify the time dimension and variable...
        oStaTimeEntry = ncf_getDimEntry(oStaInv,"time")
        oStaTimeVar => ncf_getVariable(oStaInv,"time",bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_getVariable.',irc
           return
        end if
        if (.not.bok) then
           write(*,*) myname,'Unable to file time variable.'
           irc=347
           return
        end if
        !
        ! identify forecast reference time variable...
        oStaFrtVar => ncf_getVariable(oStaInv,"forecast_reference_time",bok,irc)
        if(irc.ne.0) then
           write(*,*)myname,'Error return from getVariable.',irc
           return
        end if
        if (.not.bok)then
           write(*,*)myname,"Unable to find forecast_reference_time."
           irc=923
           return
        end if
        oStaFrtValue=oStaFrtVar%fd(1) ! latest analysis time
     end if
  end if
  !
  ! open input ensemble file and read inventory
  !
  bok=.true.
  allocate(iEnsInv,stat=irc)
  if(irc.ne.0) then
     write(*,*)myname,'Unable to allocate inventory.',irc
     return
  end if
  leni=length(iens250,250,10)
  write(*,*)myname,' Scanning: ',iens250(1:leni)
  call ncf_openFile(iens250,iEnsInv,bok,irc)
  if (irc.ne.0) then
     write(*,*)myname,' Error return from ncf_openFile.',irc
     return
  end if
  if (bok) then
     call ncf_readInventory(iEnsInv,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,' Error return from ncf_readInventory.',irc
        return
     end if
     if (bok) then
        ! get analysis time...
        iEnsFrtVar => ncf_getVariable(iEnsInv,"forecast_reference_time",bok,irc)
        if(irc.ne.0) then
           write(*,*)myname,' Error return from ncf_getVariable forecast_reference_time.',irc
           return
        end if
     end if
     if (bok) then
        call ncf_readRealData(iEnsFrtVar,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,' Error return from ncf_readRealData A.',irc
           return
        end if
     end if
     if (bok) then
        iEnsEntry=ncf_getDimEntry(iEnsInv,dim250)
        if (iEnsEntry.eq.0) then
           lend=length(dim250,250,10)
           write(*,*)myname,'Unable to find Ensemble dimension "'//dim250(1:lend)&
                & //'" in Ens-file.'
           irc=934
           bok=.false.
           return
        end if
     end if
     if (bok) then
        iEnsTimeEntry=ncf_getDimEntry(iEnsInv,"time")
        if (iEnsTimeEntry.eq.0) bok=.false.
     end if
     if (bok) then
        iEnsTimeVar => ncf_getVariable(iEnsInv,"time",bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,' Error return from ncf_getVariable time.',irc
           return
        end if
     end if
     if (bok) then
        call ncf_readRealData(iEnsTimeVar,bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,' Error return from ncf_readRealData B.',irc
           return
        end if
     end if
     if (bok) then
        call ncf_clearDimOrder(iEnsDO)
        iEnsDO => ncf_makeDimOrder(iEnsTimeVar,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
           return
        end if
        ! ensemble time/ensemble dimension
        iEnsTimeDim => ncf_removeDimOrderEntry(iEnsDO,iEnsTimeEntry)
        if (.not.associated(iEnsTimeDim)) then
           lend=length(dim250,250,10)
           write(*,*)myname,'Unable to find Time dimension in ens-file.'
           bok=.false.
           irc=959
           return
        end if
     end if
     !
     ! check if we should make statistics file
     if (firstoSta) then ! create statistics file based on input ensemble file, add statistics dimension
        leni=length(iens250,250,10)
        write(*,*)myname,'Creating statistics from: ',iens250(1:leni)

        if(bdeb)write(*,*)myname,'First:',iEnsInv%firstvariable%next%var250(1:iEnsInv%firstvariable%next%lenv)

        oStaInv => ncf_copyInventory(iEnsInv,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_openFile.',irc
           return
        end if
        oStaTimeEntry = ncf_getDimEntry(oStaInv,"time")

        ! add "statistics" dimension
        oStaDim => ncf_createDimension(oStaInv,"statistics",12 + (1+ne)*(2+ne)/2) !
        oStaEntry = ncf_getDimEntry(oStaInv,"statistics")

        ! add "resolution" dimension
        oStaResDim => ncf_createDimension(oStaInv,"resolution",resolution) !
        oStaResEntry = ncf_getDimEntry(oStaInv,"resolution")

        ! add "threshold" dimension
        oStaThrDim => ncf_createDimension(oStaInv,"threshold",2) !
        oStaThrEntry = ncf_getDimEntry(oStaInv,"threshold")

        ! create threshold dimension order
        oStaThrDO => ncf_newDimOrder(oStaInv,irc) ! make new dimension order
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
           return
        end if

        if (oStaTimeEntry.eq.0.or.oStaEntry.eq.0.or.oStaResEntry.eq.0.or.oStaThrEntry.eq.0) then
           write(*,*)myname,'Dimension error.',oStaTimeEntry,oStaEntry,oStaResEntry,oStaThrEntry
           irc=347
           return
        end if

        !!!! call ncf_setDimOrderInventory(oStaThrDO,oStaInv)
        ! add resolution dimension
        call ncf_addDimOrderEntry(oStaThrDO,oStaResEntry,irc) ! add resolution dimension
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_addDimOrderEntry.',irc
           return
        end if
        ! add threshold dimension
        call ncf_addDimOrderEntry(oStaThrDO,oStaThrEntry,irc) ! add threshold dimension
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_addDimOrderEntry.',irc
           return
        end if
        ! add time dimension
        !call ncf_addDimOrderEntry(oStaThrDO,oStaTimeEntry,irc) ! add time dimension
        !if (irc.ne.0) then
        !   write(*,*)myname,'Error return from ncf_addDimOrderEntry.',irc
        !   return
        !end if
        
        ! loop over variables and add "statistics" dimensions to variables with "time" dimension
        iEnsVar=>iEnsInv%firstvariable%next
        do while (.not.associated(iEnsVar,iEnsInv%lastVariable))
           oStaVar => ncf_getVariable(oStaInv, iEnsVar%var250(1:iEnsVar%lenv), bok, irc)  ! get the corresponding output variable
           if (.not.bok) irc=702
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_getVariable.',irc
              irc=845
              return
           end if
           if (ncf_variableContainsDim(oStaVar,oStaTimeEntry)) then ! variable contains "time" dimension
              if (oStaVar%var250(1:oStaVar%lenv).eq."time") then ! dont add statistics to time itself
                 call ncf_initField(oStaVar,nf_fill_double,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_initField.',irc
                    return
                 end if
              else
                 oStaDO => ncf_makeDimOrder(oStaVar,irc) ! get dimension order
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
                    return
                 end if
                 ! remove ensemble dimension
                 iEnsDim => ncf_removeDimOrderEntry(oStaDO,iEnsEntry)
                 ! add statistics dimension
                 call ncf_addDimOrderEntry(oStaDO,oStaEntry,irc) ! add statistics dimension
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_addDimOrderEntry.',irc
                    return
                 end if
                 call ncf_setVariableDimOrder(oStaVar,oStaDO,irc) ! insert dimension order
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_setVariableDimOrder.',irc
                    return
                 end if
                 if (associated(iEnsDim)) then
                    lend=length(dim250,250,10)
                    call ncf_addAttribute(oStaVar,"averaged_over",dim250(1:lend))
                 end if
                 ! make sure variable is of type double
                 call ncf_setVariableTypeDouble(oStaVar,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_setVariableTypeDouble.',irc
                    return
                 end if
                 call ncf_initField(oStaVar,nf_fill_double,irc) ! initialise field
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_initField.',irc
                    return
                 end if
                 ! create accumulated fields
                 write(*,*)myname,'Creating: "',iEnsVar%var250(1:iEnsVar%lenv)//sthr,'"',&
                      & nf_double,associated(oStaInv),associated(oStaThrDO),irc
                 oStaThrVar => ncf_createVariable(oStaInv,iEnsVar%var250(1:iEnsVar%lenv)//sthr,&
                      & nf_double,oStaThrDO,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_createVariable.',irc
                    return
                 end if
                 ! make sure variable is of type double
                 call ncf_setVariableTypeDouble(oStaThrVar,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_setVariableTypeDouble.',irc
                    return
                 end if
                 call ncf_addAttribute(oStaThrVar,"long_name","threshold_statistics")
                 call ncf_addAttribute(oStaThrVar,"units","1")
                 call ncf_initField(oStaThrVar,nf_fill_double,irc) ! initialise field
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_initField.',irc
                    return
                 end if
                 call   ncf_appendVariable(oStaInv,oStaThrVar)         
              end if
           else ! initialise data... 
              if (iEnsVar%lend.eq.0) then
                 call ncf_readRealData(iEnsVar,bok,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_readData.',irc
                    return
                 end if
              end if
              call ncf_initField(oStaVar,nf_fill_double,irc) ! initialise field
              if (irc.ne.0) then
                 write(*,*)myname,'Error return from ncf_initField.',irc
                 return
              end if
              call ncf_copyField(oStaVar,iEnsVar)
           end if
           iEnsVar => iEnsVar%next
        end do
        !
        oStaFrtVar => ncf_getVariable(oStaInv,"forecast_reference_time",bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_getVariable.',irc
           return
        end if
        if (.not.bok) then
           write(*,*) myname,'Unable to get forecast_reference_time.'
           irc=347
           return
        end if
        call ncf_copyField(oStaFrtVar,iEnsFrtVar)
        oStaFrtValue=oStaFrtVar%fd(1)-1000.0D0
        oStaTimeEntry = ncf_getDimEntry(oStaInv,"time")
        oStaTimeVar => ncf_getVariable(oStaInv,"time",bok,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_getVariable.',irc
           return
        end if
        if (.not.bok) then
           write(*,*) myname,'Unable to get time.'
           irc=347
           return
        end if
        call ncf_copyField(oStaTimeVar,iEnsTimeVar)
        !
        firstoSta=.false.
     else ! check if new file has larger time dimensions... (free dimension) 
        oStaTimeDim => ncf_getVariableDimension( oStaTimeVar, oStaTimeEntry )
        if (oStaTimeDim%lim.lt.iEnsTimeDim%lim) then   ! increase time dimension
           call ncf_increaseDimension(oStaInv,oStaTimeDim,iEnsTimeDim%lim,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_increaseDimension.',irc
              return
           end if
        end if
        deallocate(oStaTimeDim,stat=irc)
     end if
     if (.not.firstoSta.and.iEnsFrtVar%fd(1).le.oStaFrtValue) then ! up to date...
        write(*,*)myname,'Statistics file is up to date.',iEnsFrtVar%fd(1),oStaFrtValue
        bok=.false.
     end if
     call ncf_closeFile(iEnsInv,irc)
     if (irc.ne.0) then
        write(*,*)myname,' Error return from ncf_closeFile.',irc
        return
     end if
  end if
  !
  if (bok) then
     call ncf_clearInventoryAttributes(iEnsInv,irc) ! delete attributes
     if (irc.ne.0) then
        write(*,*) myname,"Error return from ncf_clearInventory.",irc
        return
     end if
  end if
  !
  ! process data...
  !
  if (.not.bok) then
     write(*,*)myname,'Unable to process input ensemble files.'
  else
     if (lfldat(91)) then ! check if we have moment-dimension
        oStaMomEntry = ncf_getDimEntry(oStaInv,"moments")
        if (oStaMomEntry.eq.0) then ! create dimension
           oStaMomDim => ncf_createDimension(oStaInv,"moments",ne) !
           oStaMomEntry = ncf_getDimEntry(oStaInv,"moments")
        end if
        oStaMomDim => ncf_makeDimension(oStaInv,oStaMomEntry) ! moments
     end if
     call ncf_reopenFile(iEnsInv,bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_reopenFile.',irc
        return
     end if
     if (.not.bok) then
        write(*,*) myname,'Unable to reopen file.'
        irc=346
        return
     end if
     write(*,*)myname,"Processing:",iEnsInv%fn250(1:iEnsInv%lenf)
     !
     ! get time dimension index
     !
     iEnsTimeEntry = ncf_getDimEntry(iEnsInv,"time")
     !
     ! get forecast variable
     iEnsFrtVar => ncf_getVariable(iEnsInv,"forecast_reference_time",bok,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_getVariable.',irc
        return
     end if
     if (.not.bok) then
        write(*,*) myname,'Unable to get forecast_reference_time.'
        irc=347
        return
     end if

     ! loop over ENSEMBLE-variables
     iEnsVar=>iEnsInv%firstvariable%next
     do while (.not.associated(iEnsVar,iEnsInv%lastVariable))
        bok=.true.
        if (bdeb) write(*,*)myname,'I-var-loop start'
        if (bok) then
           bok=(ncf_variableContainsDim(iEnsVar,iEnsTimeEntry))! skip if variable does not contain "time" dimension
        end if
        if (bok) then
           bok=(iEnsVar%var250(1:iEnsVar%lenv).ne."time")
        end if
        !if (bok) then
        !   bok=(iEnsVar%var250(1:4).eq."air_")
        !end if
        if (bok) then ! set up variables
           oStaVar => ncf_getVariable(oStaInv, iEnsVar%var250(1:iEnsVar%lenv), bok, irc)
           if (.not.bok) irc=702
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_getVariable.',irc
              irc=845
              return
           end if
           ! make sure variable is of type double
           call ncf_setVariableTypeDouble(oStaVar,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_setVariableTypeDouble.',irc
              return
           end if
           ! threshold variable
           oStaThrVar => ncf_getVariable(oStaInv, iEnsVar%var250(1:iEnsVar%lenv)//sthr, bok, irc)
           if (.not.bok) irc=702
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_getVariable.',irc
              irc=845
              return
           end if
           ! make sure variable is of type double
           call ncf_setVariableTypeDouble(oStaThrVar,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_setVariableTypeDouble.',irc
              return
           end if
           ! forecast reference time variable
           oStaThrFrtValue=ncf_getRealAttribute(oStaThrVar,"frt")
           ! dimension order in output variable
           oStaDO => ncf_makeDimOrder(oStaVar,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
              return
           end if
           oStaDim => ncf_removeDimOrderEntry(oStaDO,oStaEntry)
           if (lfldat(91)) then
              ! moments variable
              oStaMomVar => ncf_getVariable(oStaInv, iEnsVar%var250(1:iEnsVar%lenv)//smom, bok, irc)
              if (irc.ne.0) then
                 write(*,*)myname,'Error return from ncf_getVariable.',irc
                 irc=845
                 return
              end if
              if (.not.bok) then
                 ! create threshold dimension order
                 call ncf_clearDimOrder(oStaMomDO)
                 oStaMomDO => ncf_copyDimOrder(oStaDO,irc) ! make new dimension order
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
                    return
                 end if
                 call ncf_addDimOrderEntry(oStaMomDO,oStaMomEntry,irc) ! add "moments" to dimension order
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_addDimOrderEntry.',irc
                    return
                 end if
                 oStaMomVar => ncf_createVariable(oStaInv,iEnsVar%var250(1:iEnsVar%lenv)//smom,&
                      & nf_double,oStaMomDO,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_createVariable.',irc
                    return
                 end if
                 call ncf_addAttribute(oStaMomVar,"long_name","statistical_moments")
                 call ncf_initField(oStaMomVar,nf_fill_double,irc) ! initialise field
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_initField.',irc
                    return
                 end if
                 call ncf_appendVariable(oStaInv,oStaMomVar)         
              end if
           end if
           ! output time/statistics dimension
           oStaTimeDim => ncf_removeDimOrderEntry(oStaDO,oStaTimeEntry)
           ! ensemble dimension order
           call ncf_clearDimOrder(iEnsDO)
           iEnsDO => ncf_makeDimOrder(iEnsVar,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
              return
           end if
           ! time/ensemble dimension
           iEnsTimeDim => ncf_removeDimOrderEntry(iEnsDO,iEnsTimeEntry)
           iEnsDim => ncf_removeDimOrderEntry(iEnsDO,iEnsEntry)
           if (.not.associated(iEnsDim)) then
              lend=length(dim250,250,10)
              write(*,*)myname,'Unable to find Ensemble dimension "'//dim250(1:lend)&
                   & //'" in Ens-file.'
              irc=958
              return
           end if
           ! define dimension order to only include ensemble dimension
           call ncf_clearDimOrder(iEnsDimDO)
           iEnsDimDO => ncf_newDimOrder(iEnsInv,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_newDimOrder.',irc
              return
           end if
           call ncf_addDimOrderDim(iEnsDimDO,iEnsDim,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_addDimOrderDim.',irc
              return
           end if
           !
           lend=length(dim250,250,10)
           oStaAvgEntry=ncf_getDimEntry(oStaInv,dim250(1:lend))
           if (oStaAvgEntry.eq.0) then
              oStaAvgDim => ncf_createDimension(oStaInv,dim250(1:lend),1)
           end if
           !
           ! make index from ensemble to statistics file
           call ncf_makeIndex(io,iEnsInv,oStaInv,irc)
           if (irc.ne.0) then
              write(*,*)myname,'Error return from ncf_makeIndex.',irc
              return
           end if
           !
           firstv=.true.
           tcnt=0
           icnt=0
           firstiEns=.true.
           do ik=iEnsTimeVar%lend,1,-1 ! loop over all times using "new" climatology
              ! get climatological time...
              jk=iEnsTimeVar%lend
              JTIMELOOP: do while (jk.ge.1)! find j-time index
                 if (mod(iEnsTimeVar%fd(jk)-iEnsTimeVar%fd(ik),86400.0D0).eq.0) then ! same time of day (seconds)
                    exit JTIMELOOP
                 end if
                 jk=jk-1
              end do JTIMELOOP
              bok=(jk.ge.1) ! did we find a matching time?
              !
              do ii=0,resolution
                 rescnt(ii)=0.0D0
              end do
              !
              call ncf_resetPos(iEnsDO,irc)
              call ncf_resetPos(oStaDO,irc)
              LDATA: do while (ncf_increment(iEnsInv,iEnsDO,irc)) ! loop other dimension than time and ensemble
                 call ncf_setDimensionValue(iEnsInv,iEnsTimeDim,ik)
                 call ncf_copyIndexPos(io,iEnsInv,oStaInv,irc)
                 if (irc.ne.0) then
                    write(*,*)myname,'Error return from ncf_copyIndexPos.',irc
                    return
                 end if
                 if (firstiEns)then
                    call ncf_readRealData(iEnsVar,bok,irc)
                    if (irc.ne.0) then
                       write(*,*)myname,'Error return from ncf_readData.',irc
                       return
                    end if
                    if (.not.bok) then
                       write(*,*) myname,'Unable to read ensemble variable. ',iEnsVar%var250(1:iEnsVar%lenv)
                       irc=347
                       return
                    end if
                    firstiEns=.false.
                 end if
                 !
                 ! locations
                 call ncf_setDimensionValue(oStaInv,oStaTimeDim,ik)     ! forecast time in j-file
                 ! locations in statistics array
                 pos=0
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 1: cnt
                 lcnt = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 2: mean
                 lmean = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 3: stdv
                 lstdv = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 4: stdv residual
                 lres = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 5: time
                 lfrt = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 6: aa
                 laa = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 7: bb
                 lbb = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 8: cc
                 lcc = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 9: dd
                 ldd = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 10: weight
                 lwgt = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 11: aa: y=aa+bb*x
                 la = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 12: bb: y=aa+bb*x
                 lb = ncf_getLocation(oStaVar)
                 do mx=0,ne
                    do my=0, ne - mx 
                       pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 11: x**mx * y**my
                       lxy(mx,my) = ncf_getLocation(oStaVar)
                    end do
                 end do
                 ! moments
                 if (lfldat(91)) then
                    do ii=1,ne
                       call ncf_setDimensionValue(oStaInv,oStaMomDim,ii) ! moments
                       ll(ii) = ncf_getLocation(oStaMomVar)
                    end do
                 end if
                 ! climatological locations
                 call ncf_setDimensionValue(oStaInv,oStaTimeDim,jk)     ! climatology time
                 pos=0
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 1: cnt
                 clcnt = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 2: mean
                 clmean = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 3: stdv
                 clstdv = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 4: stdv residual
                 clres = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 5: time
                 clfrt = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 6: aa
                 claa = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 7: bb
                 clbb = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 8: cc
                 clcc = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 9: dd
                 cldd = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 10: weight
                 clwgt = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 11: aa: y=aa+bb*x
                 cla = ncf_getLocation(oStaVar)
                 pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 12: bb: y=aa+bb*x
                 clb = ncf_getLocation(oStaVar)
                 do mx=0,ne
                    do my=0, ne - mx 
                       pos=pos+1; call ncf_setDimensionValue(oStaInv,oStaDim,pos) ! 11: x**mx * y**my
                       clxy(mx,my) = ncf_getLocation(oStaVar)
                    end do
                 end do
                 !
                 xx = (iEnsTimeVar%fd(ik) - iEnsFrtVar%fd(1))*dpdt ! get lead time in days
                 if (xx.lt.0.or.xx.gt.30.0D0) then
                    write(*,*)myname,'Invalid lead time...',xx,ik,iEnsTimeVar%fd(ik) - iEnsFrtVar%fd(1)
                    irc=456
                    return
                 end if
                 sx(0)=1.0D0
                 do mm=1,ne
                    sx(mm) = sx(mm-1) * xx
                 end do
                 !
                 ! offset previous statistics to new analysis time (x=0)...
                 if (oStaVar%fd(lfrt).ne.nf_fill_double) then
                    dtx=(iEnsFrtVar%fd(1)-oStaVar%fd(lfrt))*dpdt
                    dx=-dtx
                    if (oStaVar%fd(lwgt).ne.nf_fill_double) then
                       do mx=0,ne
                          do my=0, ne - mx 
                             sxy(mx,my)=oStaVar%fd(lxy(mx,my))
                          end do
                       end do
                       !
                       do mx=0,ne
                          do my=0, ne - mx 
                             if (mx.eq.0) then
                             else if (mx.eq.1) then
                                sxy(mx,my)=sxy(mx,my)+sxy(mx-1,my)*dx
                             else if (mx.eq.2) then
                                sxy(mx,my)=sxy(mx,my)&
                                     & + 2.0D0*sxy(mx-1,my)*dx &
                                     & + sxy(mx-2,my)*dx*dx
                             else if (mx.eq.3) then
                                sxy(mx,my)=sxy(mx,my) &
                                     & + 2.0D0*sxy(mx-1,my)*dx &
                                     & + sxy(mx-2,my)*dx*dx
                             else if (mx.eq.4) then
                                sxy(mx,my)=sxy(mx,my) &
                                     & + 4.0D0*sxy(mx-1,my)*dx &
                                     & + 6.0D0*sxy(mx-2,my)*dx*dx &
                                     & + 4.0D0*sxy(mx-3,my)*dx*dx*dx &
                                     & + sxy(mx-4,my)*dx*dx*dx*dx
                             end if
                          end do
                       end do
                       ! sanity check
 !                      if (sxy(0,0).gt.1.0D0) then
 !                         write(*,*)myname,'Insane Ex:',sxy(1,0),sxy(0,0),dx,&
 !                              & oStaVar%fd(lxy(1,0)),oStaVar%fd(lxy(0,0))
 !                         irc=342
 !                         return
 !                      end if
                       do mx=0,ne
                          do my=0, ne - mx 
                             oStaVar%fd(lxy(mx,my))=sxy(mx,my)
                          end do
                       end do
!                    else if (oStaVar%fd(lxy(1,0)).ne.nf_fill_double) then
!                       write(*,*)MYNAME,'No weight defined, how strange.',oStaVar%fd(lxy(1,0))
                    end if
                 else
                    dtx=0.0D0
                 end if
                 oStaVar%fd(lfrt)   = iEnsFrtVar%fd(1)
                 !
                 ! get climatological moments
                 if (oStaVar%fd(cla).ne.nf_fill_double.and. &
                      & oStaVar%fd(clb).ne.nf_fill_double) then ! we have a valid climatology
                    dx = xx + (oStaVar%fd(lfrt) - oStaVar%fd(clfrt))*dpdt ! dt==xx since climatology is up to date...
                    mu = oStaVar%fd(cla) + oStaVar%fd(clb) * dx
                    do mx=0,ne
                       do my=0, ne - mx 
                          sxy(mx,my)=oStaVar%fd(clxy(mx,my))
                       end do
                    end do
                    if (sxy(0,0)-1.0D0.gt.1.0D-10) then
                       write(*,*) myname,'Non-normalised moments:',sxy
                       irc=955
                       return
                    end if
                    call smoment(ne,sxy,ff,oStaVar%fd(cla),oStaVar%fd(clb),gm,climOk)
                    if (lfldat(91)) then
                       do ii=1,ne
                          oStaMomVar%fd(ll(ii)) = ff(ii)
                       end do
                    end if
                 else
                    mu=nf_fill_double
                    gm=nf_fill_double
                    if (lfldat(91)) then
                       do ii=1,ne
                          oStaMomVar%fd(ll(ii)) = nf_fill_double
                       end do
                    end if
                    climOk=.false.
                 end if
                 !
                 ! get a,b,c in transormation u = a + b*x + c*x**3
                 if (climOk) then ! we have climatology
                    bok=.true.
                    xabc(1)=oStaVar%fd(laa)
                    xabc(2)=oStaVar%fd(lbb)
                    xabc(3)=oStaVar%fd(lcc)
                    xabc(4)=oStaVar%fd(ldd)
                    call getABC(NABC,XABC,ne,ff,bok,IRC)
                    if (irc.ne.0) then
                       write(*,*) 'ne=',ne
                       do ii=0,ne
                          do jj=0,ne-ii
                             write(*,'(X,A,I0,A,I0,A,F20.10)') &
                                  &'sxy(',ii,',',jj,')=',oStaVar%fd(clxy(ii,jj))
                          end do
                       end do
                       write(*,*)'a=',oStaVar%fd(cla)
                       write(*,*)'b=',oStaVar%fd(clb)
                       do ii=1,ne
                          write(*,'(X,A,I0,A,F20.10)') &
                               & 'E(',ii,')=',ff(ii)
                       end do
                       write(*,*)'mu=',mu
                       write(*,*)'gm=',gm
                       write(*,*)myname,'Error return from getABC.',nabc,xabc,ne,ff,mu,gm
                       if (.not.bok) irc=457
                       return
                    else if (.not.bok) then
                       oStaVar%fd(laa)=nf_fill_double
                       oStaVar%fd(lbb)=nf_fill_double
                       oStaVar%fd(lcc)=nf_fill_double
                       oStaVar%fd(ldd)=nf_fill_double
                       transOk=.false.
                    else
                       oStaVar%fd(laa)=xabc(1)
                       oStaVar%fd(lbb)=xabc(2)
                       oStaVar%fd(lcc)=xabc(3)
                       oStaVar%fd(ldd)=xabc(4)
                       transOk=.true.
                    end if
                 else
                    transOk=.false.
                 end if
                 !
                 ! update count
                 !
                 if (oStaVar%fd(lcnt).eq.nf_fill_double) then ! first time
                    cnt=1.0D0
                 else
                    cnt=1.0D0+oStaVar%fd(lcnt)
                 end if
                 oStaVar%fd(lcnt)   = min(9999.0D0,cnt)
                 !
                 ! update statistics
                 !
                 scnt=0.0D0
                 do mx=0,ne
                    do my=0, ne - mx 
                       sxy(mx,my)=0.0D0
                    end do
                 end do
                 cres=0.0D0
                 call ncf_resetPos(iEnsDimDO,irc)
                 ENS: do while (ncf_increment(iEnsInv,iEnsDimDO,irc)) ! loop over ensemble members
                    li = ncf_getLocation(iEnsVar)
                    ym = iEnsVar%fd(li)                               ! ensemble member
                    if (ym.ne.nf_fill_double) then ! we have valid ensemble value
                       ! update statistics
                       scnt=scnt+1.0D0
                       sy(0)=1.0D0
                       do mm=1,ne
                          sy(mm) = sy(mm-1) * ym
                       end do
                       do mx=0,ne
                          do my=0, ne - mx
                             sxy(mx,my) = sxy(mx,my) + sy(my) * sx(mx)
                          end do
                       end do
                       if (climOk) then
                          if (transOk) then
                             val = (ym - mu)*gm
                             ! calculate transformed variable here... XXXXXXXXXXXXXXXXXXXXXXx
                             valu= getu(xabc(1),xabc(2),xabc(3),xabc(4),val,bok,irc)
                             if (irc.ne.0) then
                                write(*,*)myname,'Error return from getu.',irc
                                return
                             end if
                             if (bok) then
                                ipos=nint(real(resolution)/2.0D0 + resfact*valu)
                                if ((ipos.gt.resolution .or.ipos.lt.1)) then
                                   nerr=nerr+1
                                   if (nerr.lt.100) then
                                      write(*,*)myname,'Trans: ',iEnsVar%var250(1:iEnsVar%lenv),&
                                           & xabc(1),xabc(2),xabc(3),xabc(4),val,valu
                                      write(*,*)myname,'Out of S-range: ',iEnsVar%var250(1:iEnsVar%lenv),&
                                           &ipos,resolution,iEnsVar%fd(li),ym,mu,gm
                                      if (nerr.lt.3) then
                                         write(*,*) 'ne=',ne
                                         do ii=0,ne
                                            do jj=0,ne-ii
                                               write(*,*) 'sxy(',ii,',',jj,')=',sxy(ii,jj)/scnt
                                            end do
                                         end do
                                         write(*,*)'a=',oStaVar%fd(cla)
                                         write(*,*)'b=',oStaVar%fd(clb)
                                         do ii=1,ne
                                            write(*,*) 'ff(',ii,')=',ff(ii)
                                         end do
                                         write(*,*)'mu=',mu
                                         write(*,*)'gm=',gm
                                      end if
                                   end if
                                   orm(5)=orm(5)+1
                                else
                                   ipos=max(2,min(resolution,ipos))
                                   rescnt(ipos)=rescnt(ipos)+1.0D0
                                   rescnt(0)=rescnt(0)+1.0D0
                                   tcnt=tcnt+1
                                   ook(5)=ook(5)+1
                                end if
                                ook(4)=ook(4)+1
                             else
                                orm(4)=orm(4)+1 ! invalid transformation
                             end if
                             ook(3)=ook(3)+1
                          else
                             orm(3)=orm(3)+1
                          end if ! trans
                          ook(2)=ook(2)+1
                       else
                          orm(2)=orm(2)+1 ! missing climatology/trans
                       end if ! clim
                       ook(1)=ook(1)+1
                    else
                       orm(1)=orm(1)+1 ! missing value
                    end if
                    icnt=icnt+1
                 end do ENS
                 ! store results
                 if (scnt.gt.0.0D0) then
                    scnt=max(1.0D0,scnt)
                    do mx=0,ne
                       do my=0, ne - mx
                          sxy(mx,my) = sxy(mx,my)/scnt
                       end do
                    end do
                    !
                    ! old weight
                    if (oStaVar%fd(lwgt).eq.nf_fill_double) then
                       oldWgt=0.0D0
                    else
                       ! forecast reference time (frt) is in seconds
                       oldWgt=oStaVar%fd(lwgt)*exp(-dtx/dtrelax)
                    end if
                    ! get weights
                    totalWgt=1.0D0+oldWgt    ! total weight
                    wgt0=oldWgt/totalWgt     ! old weight
                    wgt1=1.0D0/totalWgt      ! new weight
                    !
                    ! store statistics
                    oStaVar%fd(lwgt)  = totalWgt
                    do mx=0,ne
                       do my=0, ne - mx 
                          oStaVar%fd(lxy(mx,my)) = oStaVar%fd(lxy(mx,my))* wgt0 + sxy(mx,my) * wgt1
                       end do
                    end do
                    !
                    ! make climatological mean, stdv and stdv residual
                    if (climOk) then
                       cmean=mu
                       cstdv=1.0D0/max(1.0D-10,gm)
                       oStaVar%fd(lmean)  = cmean
                       oStaVar%fd(lstdv)  = cstdv
                       ! residual
                       cres=cres + sxy(0,2) - 2.0D0*sxy(0,1)*cmean + cmean*cmean
                       if (oStaVar%fd(lres).eq.nf_fill_double) then
                          oStaVar%fd(lres)=cres ! residual
                       else
                          oStaVar%fd(lres) = (oStaVar%fd(lres) * wgt0) + (cres * wgt1)
                       end if
                    end if
                    !
                    ! make current a,b
                    if (cnt.gt.dtrelax.and.&     ! linear regression
                         & oStaVar%fd(lxy(1,1)).ne.nf_fill_double.and.&
                         & oStaVar%fd(lxy(0,1)).ne.nf_fill_double) then
                       ! solve for a and b in a + b*dx = y**e using linear regression
                       xtx(1,1) = 1
                       xtx(1,2) = oStaVar%fd(lxy(1,0))
                       xtx(2,1) = oStaVar%fd(lxy(1,0))
                       xtx(2,2) = oStaVar%fd(lxy(2,0))
                       xty(1)   = oStaVar%fd(lxy(0,1))
                       xty(2)   = oStaVar%fd(lxy(1,1))
                       ! Calculate the inverse determinant of the matrix
                       det = (XTX(1,1)*XTX(2,2) - XTX(1,2)*XTX(2,1))
                       if (abs(det).gt.1.0D-10) then
                          detinv=1.0D0/det
                          
                          ! Calculate the inverse of the matrix
                          B(1,1) = +detinv * XTX(2,2)
                          B(2,1) = -detinv * XTX(2,1)
                          B(1,2) = -detinv * XTX(1,2)
                          B(2,2) = +detinv * XTX(1,1)
                          A(1)=B(1,1)*XTY(1) + B(1,2)*XTY(2)
                          A(2)=B(2,1)*XTY(1) + B(2,2)*XTY(2)
                       else
                          A(1)=XTY(1)
                          A(2)=0.0D0
                       end if
                       
                       !write(*,*)myname,'AB:',iEnsVar%var250(1:iEnsVar%lenv),ik,mm,a(1),a(2)
                       
                    else
                       A(1)=nf_fill_double        
                       A(2)=nf_fill_double        
                    end if
                    !
                    oStaVar%fd(la)    = A(1)
                    oStaVar%fd(lb)    = A(2)
                 end if
              end do LDATA
              ! update probability distribution...
              if (rescnt(0).gt.0.0D0) then ! we have some data
                 call ncf_setDimensionValue(oStaInv,oStaTimeDim,ik)    ! accumulate in current
                 call ncf_setDimensionValue(oStaInv,oStaResDim,1)
                 call ncf_setDimensionValue(oStaInv,oStaThrDim,1)       ! distribution
                 lwgt = ncf_getLocation(oStaThrVar)
                 tWgt=oStaThrVar%fd(lwgt)
                 if (tWgt.eq.nf_fill_double) then ! first time
                    oldWgt=0.0d0
                 else
                    oldWgt=tWgt*exp(-(iEnsFrtVar%fd(1)-oStaThrFrtValue)/&
                         & (dtdist*86400.0D0))
                 END IF
                 totalWgt=1.0D0+oldWgt     ! total weight
                 oStaThrVar%fd(lwgt)=totalWgt
                 wgt0=oldWgt/totalWgt      ! old weight
                 wgt1=1.0D0/totalWgt       ! new weight
                 resacc=0.0D0
                 do ipos=2,resolution
                    call ncf_setDimensionValue(oStaInv,oStaResDim,ipos)
                    call ncf_setDimensionValue(oStaInv,oStaTimeDim,ik) ! current time
                    call ncf_setDimensionValue(oStaInv,oStaThrDim,1)   ! pdf
                    lt1 = ncf_getLocation(oStaThrVar)
                    call ncf_setDimensionValue(oStaInv,oStaThrDim,2)   ! acc
                    lt2 = ncf_getLocation(oStaThrVar)
                    call ncf_setDimensionValue(oStaInv,oStaTimeDim,jk) ! climatology
                    call ncf_setDimensionValue(oStaInv,oStaThrDim,1)   ! pdf
                    clt1 = ncf_getLocation(oStaThrVar)
                    call ncf_setDimensionValue(oStaInv,oStaThrDim,2)   ! acc
                    clt2 = ncf_getLocation(oStaThrVar)
                    rescnt(ipos)=rescnt(ipos)/max(1.0D0,rescnt(0))
                    if (oStaThrVar%fd(lt1).eq.nf_fill_double) oStaThrVar%fd(lt1)=0.0D0
                    rescnt(ipos)=oStaThrVar%fd(lt1) + wgt1 * (rescnt(ipos) - oStaThrVar%fd(lt1))
                    resacc=resacc+rescnt(ipos)
                    oStaThrVar%fd(lt1)=rescnt(ipos)     ! update current
                    !oStaThrVar%fd(clt1)=rescnt(ipos)   ! update climatology
                    oStaThrVar%fd(lt2)=resacc           ! update current
                    !oStaThrVar%fd(clt2)=resacc         ! update climatology
                 end do ! resolution
                 ! write(*,*)myname,tcnt, oStaThrVar%var250(1:oStaThrVar%lenv),ik
                 if (abs(resacc-1.0).gt.1.0D-10) then
                    write(*,*)myname,'Non-normalised ACC:',resacc,totalWgt,&
                         & oStaThrFrtValue
                 end if
              end if ! do we have probability distribution data?
              oStaThrFrtValue=iEnsFrtVar%fd(1)
           end do ! ik
           !
           write(*,'(X,A,A,X,I0,X,A,X,I0,X,A)')myname,'Data count:',icnt,'->',tcnt, iEnsVar%var250(1:iEnsVar%lenv)
           !
           call ncf_setDoubleAttribute(oStaThrVar,"frt",iEnsFrtVar%fd(1))
           call ncf_setDoubleAttribute(oStaThrVar,"range",range)
           if (allocated(io)) deallocate(io)
           call ncf_clearDimOrder(iEnsDO)
           call ncf_clearDimOrder(oStaDO)
           call ncf_clearDimOrder(oStaThrDO)
           if (.not.firstiEns) then
              call ncf_clearVariable(iEnsVar,irc)
              if(irc.ne.0)return
           end if
        end if
        iEnsVar => iEnsVar%next
        if (bdeb) write(*,*)myname,'I-var-loop end',bok
     end do ! iEnsVar
     ! loop over oStaInv variables and prepare writing to file...
     if(bdeb)write(*,*) myname,'End of I-VAR.'
     ! redefine forecast_reference_time
     oStaFrtValue=max(oStaFrtValue,iEnsFrtVar%fd(1))
     !
     call ncf_closeFile(iEnsInv,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_closeFile.',irc
        return
     end if
     !
     call ncf_closeFile(iEnsInv,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_closeFile.',irc
        return
     end if
     !
     ! call ncf_checkInventory(iEnsInv,irc)
     !
     leni=length(iens250,250,10)
     write(*,*)myname,"Clearing:",iens250(1:leni)
     call ncf_clearInventory(iEnsInv,irc)
     if (irc.ne.0) then
        write(*,*)myname,'Error return from ncf_clearInventory.',irc
        return
     end if
     !
     if(bdeb)write(*,*) myname,'Done processing.'
     !
     if (.not. firstoSta) then
        if(bdeb)write(*,*) myname,'Compressing variables.'
        call ncf_compressVariables(oStaInv,irc) ! compress variables
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_compressVariables.',irc
           return
        end if
        if(bdeb)write(*,*) myname,'Compressing dimensions.'
        call ncf_compressDimensions(oStaInv) ! compress dimensions
        !
        oStaInv%fn250=osta250
        oStaInv%lenf=length(osta250,250,10)
        if(bdeb)write(*,*) myname,'Writing statistics:',oStaInv%fn250(1:oStaInv%lenf)
        !     call ncf_checkInventory(oStaInv,irc)
        call ncf_writeNcOut(oStaInv,osta250,irc) ! open, write, close file
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_writeNcOut.',irc
           return
        end if
        leni=length(osta250,250,10)
        write(*,*)myname,"Clearing:",osta250(1:leni)
        call ncf_clearInventory(oStaInv,irc)
        if (irc.ne.0) then
           write(*,*)myname,'Error return from ncf_clearInventory.',irc
           return
        end if
     else
        write(*,*)myname,'No data available for oStaInvput file.'
     end if
  end if
  !
  ! write statistics
  !
  if(bdeb)write(*,*) myname,'Writing statistics.'
  do jj=1,10
     pst(jj)=dfloat(orm(jj))/max(1.0d0,dfloat(orm(jj)+ook(jj)))*100
  end do
  WRITE(*,*)
  WRITE(*,998) MYNAME,                 'Available data:            ', ook(1)+orm(1)
  IF (ORM(1).NE.0) WRITE(*,999) MYNAME,'Missing value:             ', -ORM(1),PST(1)
  IF (ORM(2).NE.0) WRITE(*,999) MYNAME,'Missing climatology:       ', -ORM(2),PST(2)
  IF (ORM(3).NE.0) WRITE(*,999) MYNAME,'Missing trans:             ', -ORM(3),PST(3)
  IF (ORM(4).NE.0) WRITE(*,999) MYNAME,'Invalid trans:             ', -ORM(4),PST(4)
  IF (ORM(4).NE.0) WRITE(*,999) MYNAME,'Invalid range:             ', -ORM(5),PST(5)
  WRITE(*,*)   MYNAME,     '-----------------------------------------------'
  WRITE(*,998) MYNAME,                 'Accepted data:        ', ook(5)
  !
999 FORMAT(X,A12,X,A27,I10,' (',F6.2,'%)')
998 FORMAT(X,A12,X,A27,I10)
997 FORMAT(X,A12,X,A,10(X,I10))
  !
end subroutine MNCPROB


program ttt
  use ncf
  use sort
  implicit none
  character*250 :: inp250="/lustre/storeB/immutable/archive/"//&
       & "projects/metproduction/meps/2017/03/14/meps_mbr0_extracted_2_5km_20170314T09Z.nc"
  INTEGER  :: leni,lenv
  integer, external :: length
  type(inventory),pointer  :: inp => null()   ! equals output...
  type(variable),pointer :: ivar
  character*250 :: var250
  logical :: bop
  integer irc,ii
  character*10 :: myname = "ttt"
  ncf_bdeb=.true.
  call chop0(inp250,250)
  leni=length(inp250,250,10)
  write(*,*)myname,'Opening:',inp250(1:leni)
  call ncf_openFile(inp250,inp,bop,irc)
  if (irc.ne.0.or..not.bop) then
     write(*,*)myname,' Error return from ncf_openFile.',irc
     stop
  end if
  call ncf_readInventory(inp,bop,irc)
  if (irc.ne.0.or..not.bop) then
     write(*,*)myname,' Error return from ncf_readInventory.',irc
     stop
  end if
  do ii=1,6
     write(*,*)myname,'Looping.',ii
     select case(ii)
     case (1)
        var250="time";
     case (2)
        var250="longitude";
     case (3)
        var250="latitude";
     case (4)
        var250="forecast_reference_time";
     case (5)
        var250="air_temperature_2m";
     case (6)
        var250="relative_humidity_2m";
     end select
     call chop0(var250,250)
     lenv=length(var250,250,10)
     !write(*,*)myname,'Checking.',var250(1:lenv)
     !call ncf_checkContents(inp,var250,bop,irc)
     !if (irc.ne.0.or..not.bop) then
     !   write(*,*)myname,' Error return from ncf_checkContents.',irc
     !   stop
     !end if
     write(*,*)myname,'Reading.',var250(1:lenv)
     ivar => ncf_getVariable(inp,var250(1:lenv),bop,irc)
     call ncf_readRealData(ivar,bop,irc)
     if (irc.ne.0.or..not.bop) then
        write(*,*)myname,' Error return from ncf_readRealData A.',irc
        stop
     end if
     write(*,*) ivar%fd(1:10)
  end do
  write(*,*)myname,'Closing:',inp250(1:leni)
  call ncf_closeFile(inp,irc)
  if (irc.ne.0.or..not.bop) then
     write(*,*)myname,' Error return from ncf_closeFile.',irc
     stop
  end if
end program ttt

Module ncf
# include "netcdf.inc"
  !
  integer*8, parameter :: nf_fill_int64 = -9223372036854775806_8
  integer :: notimp=0
  !
  logical :: ncf_bdeb =.false.
  type inventory
     character*10  :: lab10=""
     integer       :: lenl=0
     logical       :: initialised = .false.
     integer       :: ncid
     character*250 :: fn250
     integer       :: lenf = 0
     logical       :: opened = .false.
     integer       :: nrdim
     integer       :: nrvar
     integer       :: nrgatt
     integer       :: unlimdimid
     character*250, pointer    :: dim250(:)=>null()
     integer, pointer          :: lend(:)=>null()
     integer, pointer          :: dimid(:)=>null()
     type(position), pointer   :: pos => null()
     type(variable), pointer   :: firstVariable
     type(variable), pointer   :: lastVariable
     type(variable), pointer   :: latid=>null()
     type(variable), pointer   :: lonid=>null()
     type(variable), pointer   :: timid=>null()
     type(variable), pointer   :: terid=>null()
     type(variable), pointer   :: parid=>null()
     type(variable), pointer   :: xcoid=>null()
     type(variable), pointer   :: ycoid=>null()
     type(variable), pointer   :: gridid=>null()
     type(attribute), pointer  :: firstAttribute=>null()
     type(attribute), pointer  :: lastAttribute=>null()
  end type inventory

  type variable
     logical       :: initialised = .false.
     logical       :: readyfield  = .false.
     logical       :: ignored     = .false.
     logical       :: writefield  = .true.
     type(variable), pointer  :: prev=>null()
     type(variable), pointer  :: next=>null()
     type(inventory), pointer  :: f=>null()
     integer       :: ncid
     integer       :: varid
     character*250 :: var250
     integer       :: lenv = 0
     logical       :: lout=.false.
     character*250 :: out250
     integer       :: leno = 0
     integer       :: type = -1
     integer       :: nrdim = 0
     integer       :: nratt
     real          :: scale
     integer,    pointer :: ind(:)=>null()
     integer,    pointer :: sta(:)=>null()
     integer,    pointer :: lim(:)=>null()
     integer,    pointer :: siz(:)=>null()
     character*1,pointer::  fc(:)=>null()
     integer*1,  pointer::  f1(:)=>null()
     integer*2,  pointer::  f2(:)=>null()
     integer*4,  pointer::  f4(:)=>null()
     integer*8,  pointer::  f8(:)=>null()
     real*4,     pointer::  fr(:)=>null()
     real*8,     pointer::  fd(:)=>null()
     character*1::  fillc
     integer*1  ::  fill1 = nf_fill_int1
     integer*2  ::  fill2 = nf_fill_int2
     integer*4  ::  fill4 = nf_fill_int
     integer*8  ::  fill8 = nf_fill_int64
     real*4     ::  fillr = nf_fill_real
     real*8     ::  filld = nf_fill_double
     integer :: lens=0
     integer :: len=0
     integer :: lenc=0
     integer :: len1=0
     integer :: len2=0
     integer :: len4=0
     integer :: len8=0
     integer :: lenr=0
     integer :: lend=0
     type(attribute),pointer  :: firstAttribute=>null()
     type(attribute),pointer  :: lastAttribute=>null()
     type(attribute),pointer  :: fillAttribute=>null()
  end type variable

  type attribute
     logical       :: initialised = .false.
     type(variable), pointer  :: v=>null()
     type(attribute), pointer  :: prev=>null()
     type(attribute), pointer  :: next=>null()
     integer       :: ncid
     integer       :: varid
     integer       :: attid
     character*250  :: att250
     integer :: lena
     integer :: type
     integer :: len
     character*1, pointer ::  ac(:)=>null()
     integer*1,  pointer::  a1(:)=>null()
     integer*2,  pointer::  a2(:)=>null()
     integer*4,  pointer::  a4(:)=>null()
     integer*8,  pointer::  a8(:)=>null()
     real*4,     pointer::  ar(:)=>null()
     real*8,     pointer::  ad(:)=>null()
  end type attribute

  type position
     integer :: maxdim=0
     integer :: nrdim =0
     integer, pointer :: pos(:)
     integer, pointer :: sta(:)
     integer, pointer :: lim(:)
     integer, pointer :: siz(:)
     integer        :: loc =0
  end type position

  type dimension
     logical          :: initialised = .false.
     type(dimension), pointer  :: prev=>null()
     type(dimension), pointer  :: next=>null()
     integer          :: ind =0
     integer          :: sta =0
     integer          :: lim =0
     integer          :: siz = 0
     logical          :: increase=.true.
  end type dimension

  type dimensionOrder
     integer :: nrdim = 0
     type(inventory), pointer  :: f=>null()
     type(dimension), pointer  :: firstDimension
     type(dimension), pointer  :: lastDimension
     logical :: resetPos=.true.
  end type dimensionOrder

  type weight
     integer :: nweight = 0
     type(position), pointer :: pos(:)=>null()
     real, pointer           :: w(:)=>null()
  end type weight

contains

  ! MATHEMATICS ******************************************************

  real function ncf_xycross(x1,y1,x2,y2,x3,y3)
    implicit none
    real :: x1,y1,x2,y2,x3,y3
    real :: dx1,dy1,dx2,dy2
    dx1=x3-x1
    dy1=y3-y1
    dx2=x2-x1
    dy2=y2-y1
    ncf_xycross = dx1 * dy2 - dx2 * dy1
  end function ncf_xycross

  real function ncf_degtor(x)
    implicit none
    real x
    real pi
    parameter (pi=3.14159265359)
    ncf_degtor=x*pi/180.
  end function ncf_degtor

  real function ncf_rtodeg(x)
    implicit none
    real x
    real pi
    parameter (pi=3.14159265359)
    ncf_rtodeg=x*180./pi
  end function ncf_rtodeg

  real function ncf_sindeg(x)
    implicit none
    real x
    ncf_sindeg=sin(ncf_degtor(x))
  end function ncf_sindeg

  real function ncf_cosdeg(x)
    implicit none
    real x
    ncf_cosdeg=cos(ncf_degtor(x))
  end function ncf_cosdeg

  real function ncf_tandeg(x)
    implicit none
    real x
    ncf_tandeg=tan(ncf_degtor(x))
  end function ncf_tandeg

  real function ncf_asindeg(x)
    implicit none
    real x
    ncf_asindeg=ncf_rtodeg(asin(x))
  end function ncf_asindeg

  real function ncf_acosdeg(x)
    implicit none
    real x
    ncf_acosdeg=ncf_rtodeg(acos(x))
  end function ncf_acosdeg

  real function ncf_atandeg(x)
    implicit none
    real x
    ncf_atandeg=ncf_rtodeg(atan(x))
  end function ncf_atandeg

  real function ncf_atan2deg(y,x)
    implicit none
    real y,x
    ncf_atan2deg=ncf_rtodeg(atan2(y,x))
  end function ncf_atan2deg
  !
  ! interpolates current reference position to experiment position
  !
  subroutine ncf_interpolate2D(r,e,ix,iy,wgt,biok,irc)
    implicit none
    type(inventory) r  ! reference
    type(inventory) e  ! experiment
    type(dimension) ix ! experiment x-dimension
    type(dimension) iy ! experiment y-dimension
    type(weight) wgt   ! weight of surrounding grid points
    logical :: biok    ! is interpolation successful?
    integer irc        ! error return code
    !
    real :: xf,yf      ! current reference position
    integer :: nxs,nys,nzs
    real :: xs(2,2), ys(2,2)  ! interpolated position
    logical :: xok,yok
    logical :: bdone,changed, inside
    real :: bc,rc,tc,lc,fc,tbc,lrc
    real :: brc,rtc,tlc,lbc
    integer :: jbc,jrc,jtc,jlc,jbrc,jrtc,jtlc,jlbc,ii,jj
    integer :: iterations = 0
    character*20 :: myname = "ncf_interpolate2D"
    !
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ! search for match where position is within gridcell
    xf=ncf_valuePosition(r%lonid,irc)
    if (irc.ne.0) then
       write(*,*)myname,'error return from ncf_valuePosition.',irc
       return
    end if
    yf=ncf_valuePosition(r%latid,irc)
    if (irc.ne.0) then
       write(*,*)myname,'error return from ncf_valuePosition.',irc
       return
    end if
    bdone=.false.
    ! iterations=0
    do while (.not. bdone)
       iterations=iterations+1
       xok=.true.
       yok=.true.
       xok=ncf_incrementDimension(e,ix,xok,1,0,irc) ! -> bottom right
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_incrementDimension.',irc
          return
       end if
       xs(2,1)=ncf_valuePosition(e%lonid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_valuePosition.',irc
          return
       end if
       ys(2,1)=ncf_valuePosition(e%latid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_valuePosition.',irc
          return
       end if
       !
       yok=ncf_incrementDimension(e,iy,yok,1,0,irc) ! -> top right
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_incrementDimension.',irc
          return
       end if
       xs(2,2)=ncf_valuePosition(e%lonid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_valuePosition.',irc
          return
       end if
       ys(2,2)=ncf_valuePosition(e%latid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_valuePosition.',irc
          return
       end if
       !
       xok=ncf_decrementDimension(e,ix,xok,1,0,irc) ! -> top left
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_decrementDimension.',irc
          return
       end if
       xs(1,2)=ncf_valuePosition(e%lonid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_valuePosition.',irc
          return
       end if
       ys(1,2)=ncf_valuePosition(e%latid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_valuePosition.',irc
          return
       end if
       !
       yok=ncf_decrementDimension(e,iy,yok,1,0,irc) ! -> bottom left
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_decrementDimension.',irc
          return
       end if
       xs(1,1)=ncf_valuePosition(e%lonid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_valuePosition.',irc
          return
       end if
       ys(1,1)=ncf_valuePosition(e%latid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'error return from ncf_valuePosition.',irc
          return
       end if
       if (ncf_bdeb) then
          write(*,*)myname,"Experiment position..."
          call ncf_printPos(e%pos)
       end if
       !
       if (.not. xok.or. .not.yok) then
          biok=.false.
          xok=ncf_decrementDimension(e,ix,xok,1,0,irc) ! -> top left
          if (irc.ne.0) then
             write(*,*)myname,'error return from ncf_decrementDimension.',irc
             return
          end if
          yok=ncf_decrementDimension(e,iy,yok,1,0,irc) ! -> bottom left
          if (irc.ne.0) then
             write(*,*)myname,'error return from ncf_decrementDimension.',irc
             return
          end if
          call ncf_printPos(e%pos)
          write(*,*) 'ncf_interpolate2D Invalid starting position.'
          irc=700
          return
       else
          !write(*,*)'XS:',xf,xs
          !write(*,*)'YS:',yf,ys
          !call exit(0)

          ! walk around grid cell border and calculate cross product
          bc=ncf_xycross(xs(1,1),ys(1,1), xf,yf, xs(2,1),ys(2,1)) ! bottom
          rc=ncf_xycross(xs(2,1),ys(2,1), xf,yf, xs(2,2),ys(2,2)) ! right
          tc=ncf_xycross(xs(2,2),ys(2,2), xf,yf, xs(1,2),ys(1,2)) ! top
          lc=ncf_xycross(xs(1,2),ys(1,2), xf,yf, xs(1,1),ys(1,1)) ! left

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

          brc=ncf_xycross(xs(1,1),ys(1,1), xs(2,2),ys(2,2), xs(2,1),ys(2,1)) ! bottom-right
          rtc=ncf_xycross(xs(2,1),ys(2,1), xs(1,2),ys(1,2), xs(2,2),ys(2,2)) ! right-top
          tlc=ncf_xycross(xs(2,2),ys(2,2), xs(1,1),ys(1,1), xs(1,2),ys(1,2)) ! top-left
          lbc=ncf_xycross(xs(1,2),ys(1,2), xs(2,1),ys(2,1), xs(1,1),ys(1,1)) ! left-bottom

          jbrc=nint(sign(1.0D0,brc))
          jrtc=nint(sign(1.0D0,rtc))
          jtlc=nint(sign(1.0D0,tlc))
          jlbc=nint(sign(1.0D0,lbc))

          if (ncf_bdeb) then
             write(*,'("XYSEARCH Bxt:",F10.5,"  Rxt:",F10.5,"  Txt:",F10.5,"  Lxt:",F10.5,X,4(X,I2))') &
                  & bc,rc,tc,lc,       jbc,jrc,jtc,jlc
             write(*,'("XYSEARCH BxR:",F10.5,"  RxT:",F10.5,"  TxL:",F10.5,"  LxB:",F10.5,X,4(X,I2))') &
                  & brc,rtc,tlc,lbc,   jbrc,jrtc,jtlc,jlbc
             write(*,'("XYSEARCH ",A," (",F0.5,",",F0.5,")")') &
                  & "targ =",xf,yf
             write(*,'("XYSEARCH ",A," (",F0.5,",",F0.5,")")') &
                  & "[1,1]=",xs(1,1),ys(1,1)
             write(*,'("XYSEARCH ",A," (",F0.5,",",F0.5,")")') &
                  & "[2,1]=",xs(2,1),ys(2,1)
             write(*,'("XYSEARCH ",A," (",F0.5,",",F0.5,")")') &
                  & "[2,2]=",xs(2,2),ys(2,2)
             write(*,'("XYSEARCH ",A," (",F0.5,",",F0.5,")")') &
                  & "[1,2]=",xs(1,2),ys(1,2)
          end if
          !     
          !     if cross product has same sign as axis-sign=>inside, else outside
          !     (make sure we do not cross border while searching...)
          !     
          inside=.true.
          changed=.false.
          if (jbc.eq.-jbrc) then  ! decrease y
             if (ncf_decrementDimension(e,iy,yok,1,0, irc)) then
                changed=.true.

                !write(*,*)'interpolate2d Y-',e%pos%pos(iy%ind),iterations,bc,brc

             end if
             if (irc.ne.0) then
                write(*,*)myname,'error return from ncf_decrementDimension.',irc
                return
             end if
             if (ncf_bdeb)write(*,*)myname,'Decreasing Y',changed
             inside=.false.
          else if (jtc.eq.-jtlc) then ! increase y
             if (ncf_incrementDimension(e,iy,yok,1,1, irc)) then
                changed=.true.


                !write(*,*)'interpolate2d Y+',e%pos%pos(iy%ind),iterations,tc,tlc

             end if
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_incrementDimension.',irc
                return
             end if
             if (ncf_bdeb)write(*,*)myname,'Increasing Y',changed
             inside=.false.
          end if
          if (jlc.eq.-jlbc) then  ! decrease x
             if (ncf_decrementDimension(e,ix,xok,1,0, irc)) then
                changed=.true.

                !write(*,*)'interpolate2d X-',e%pos%pos(ix%ind),iterations,lc,lbc

             end if
             if (irc.ne.0) then
                write(*,*)myname,'error return from ncf_decrementDimension.',irc
                return
             end if
             if (ncf_bdeb)write(*,*)myname,'Decreasing X',changed
             inside=.false.
          else if (jrc.eq.-jrtc) then ! increase x
             if (ncf_incrementDimension(e,ix,xok,1,1, irc)) then
                changed=.true.

                !write(*,*)'interpolate2d X+',e%pos%pos(ix%ind),iterations,rc,rtc

             end if
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_incrementDimension.',irc
                return
             end if
             if (ncf_bdeb)write(*,*)myname,'Increasing X',changed
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
                call ncf_printPos(e%pos)
                write(*,*) myname,' Invalid grid:',rc,tc,lc,bc,fc
                irc=701
                return
             end if
             wgt%w(1)=lc*tc/fc ! bottom right
             wgt%w(2)=lc*bc/fc ! top right
             wgt%w(3)=rc*bc/fc ! top left
             wgt%w(4)=rc*tc/fc ! bottom left
             ! store position vector
             do ii=1,wgt%nweight
                if (wgt%pos(ii)%maxdim .lt. e%pos%nrdim) then ! sufficient dimensions?
                   call ncf_allocatePos(e%pos%nrdim,wgt%pos(ii),irc)
                   if (irc.ne.0) then
                      write(*,*)myname,'error return from ncf_allocatePos.',irc
                      return
                   end if
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
             if (ncf_bdeb) write(*,*)"XYSEARCH Done:",wgt%w
          else if (.not.changed) then ! not inside + no valid step
             if (ncf_bdeb) then
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

    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_interpolate2D

  real function ncf_getDist(LATA,LONA,LATB,LONB)
    IMPLICIT NONE
    SAVE
    REAL   LONA,LATA,LONB,LATB
    REAL   CDIFF
    character*20 :: myname = "ncf_getDist"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    CDIFF=NCF_SINDEG(LATA)*NCF_SINDEG(LATB)+NCF_COSDEG(LATA)*&
         & NCF_COSDEG(LATB)*NCF_COSDEG(LONA-LONB)
    CDIFF=MAX(-1.0D0,MIN(1.0D0,CDIFF)) ! handle truncation errors
    ncf_getDist=ABS(NCF_ACOSDEG(CDIFF))
    if (ncf_bdeb) write(*,*)myname,' Done.'
    RETURN
  END function ncf_getDist

  !@#  OPEN/CLOSE FILE ******************************************************

  subroutine ncf_openFile(fn250,f,biok,irc)
    implicit none
    ! 
    character*250 :: fn250
    type (inventory), pointer :: f
    logical :: biok
    integer :: irc
    ! 
    character*250 :: buff250
    integer :: lenb, lenf
    integer, external :: length
    INTEGER :: ret,CHUNKSIZEHINT
    logical :: proceed
    character*20 :: myname = "ncf_openFile"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    !
    lenf=length(fn250,250,10)
    biok=.false.
    if (.not.associated(f)) allocate(f,stat=irc)
    if (irc.ne.0) then
       write(buff250,"(A)") " Unable to allocate inventory for: "//fn250(1:lenf)
       call chop0(buff250,250)
       lenb=length(buff250,250,10)
       write(*,*)'ncf_openFile '//buff250(1:lenb)
       f%opened=.false.
       irc=833
       return
    end if
    f%fn250=fn250
    f%lenf=lenf
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
       write(*,*)myname,buff250(1:lenb)
       f%opened=.false.
       irc=834
       return
    else
       write(buff250,"(A)") "  Opened: "//f%fn250(1:f%lenf)
       call chop0(buff250,250)
       lenb=length(buff250,250,10)
       ! write(*,*)'ncf_openFile '//buff250(1:lenb)
       biok=.true.
       f%opened=.true.
    endif
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_openFile

  subroutine ncf_reopenFile(f,biok,irc)
    implicit none
    ! 
    character*250 fn250
    type (inventory) f
    logical biok
    integer irc
    INTEGER ret,CHUNKSIZEHINT
    character*250 :: buff250
    integer :: lenb
    integer, external :: length
    type(variable), pointer :: v
    logical :: proceed
    character*20 :: myname = "ncf_reopenFile"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
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
       write(buff250,"(A)") "  Re-opening: ",f%fn250(1:f%lenf),f%lenf
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
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_reopenFile

  subroutine ncf_closeFile(f,irc)
    implicit none
    type(inventory), pointer :: f
    integer irc
    integer ret
    integer lenb
    character*250 buff250
    integer, external :: length
    character*20 :: myname = "ncf_closeFile"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (f%opened) then
       f%opened=.false.
       ret=NF_CLOSE(f%ncid)        ! end definitions: leave define mode
       if (ret .ne. NF_NOERR) then
          write(*,*) "ncf_closeFile ERROR return from NF_CLOSE:",nf_strerror(ret),f%ncid,nf_strerror(ret)
          irc=790
          return
       end if
    end if
    return
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_closeFile

  function ncf_copyFile(f,irc)
    implicit none
    type(inventory), pointer :: ncf_copyFile
    type(inventory), pointer :: f
    integer :: irc
    integer :: ii
    type(variable),pointer :: v,fv
    type(attribute),pointer :: a,fa
    character*20 :: myname = "ncf_copeyInventory"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(ncf_copyFile,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate ncf_copyFile.'
       return
    end if
    ncf_copyFile%ncid=f%ncid
    ncf_copyFile%fn250=f%fn250
    ncf_copyFile%lenf=f%lenf
    ncf_copyFile%nrdim=f%nrdim
    ncf_copyFile%nrvar=f%nrvar
    ncf_copyFile%nrgatt=f%nrgatt
    allocate(ncf_copyFile%dim250(max(1,ncf_copyFile%nrdim)),&
         & ncf_copyFile%lend(max(1,ncf_copyFile%nrdim)),&
         & ncf_copyFile%pos,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,'Unable to allocate.',irc
       return
    end if
    allocate(ncf_copyFile%pos%pos(max(1,ncf_copyFile%nrdim)),&
         & ncf_copyFile%pos%sta(max(1,ncf_copyFile%nrdim)),&
         & ncf_copyFile%pos%lim(max(1,ncf_copyFile%nrdim)),&
         & ncf_copyFile%pos%siz(max(1,ncf_copyFile%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,'Unable to allocate.',irc
       return
    end if
    do ii=1,ncf_copyFile%nrdim
       ncf_copyFile%dim250(ii)=f%dim250(ii)
       ncf_copyFile%lend(ii)=f%lend(ii)
    end do
    ncf_copyFile%pos%nrdim=f%pos%nrdim
    do ii=1,ncf_copyFile%pos%nrdim
       ncf_copyFile%pos%pos(ii)=f%pos%pos(ii)
       ncf_copyFile%pos%sta(ii)=f%pos%sta(ii)
       ncf_copyFile%pos%lim(ii)=f%pos%lim(ii)
       ncf_copyFile%pos%siz(ii)=f%pos%siz(ii)
    end do
    ncf_copyFile%pos%loc=f%pos%loc
    ncf_copyFile%unlimdimid=f%unlimdimid
    !
    allocate(ncf_copyFile%firstVariable,ncf_copyFile%lastVariable, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,'Unable to allocate.',irc
       return
    end if
    ncf_copyFile%firstVariable%next=>ncf_copyFile%lastVariable
    ncf_copyFile%lastVariable%prev=>ncf_copyFile%firstVariable
    !
    fv=>f%firstVariable%next
    do while (.not.associated(fv,target=f%lastVariable))
       v=>ncf_copyVariable(fv,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from copyVariable.',irc
          return
       end if
       call ncf_copyVariableField(v,fv,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from copyVariableField.',irc
          return
       end if
       v%f=>ncf_copyFile
       v%prev=>ncf_copyFile%lastVariable%prev
       v%next=>ncf_copyFile%lastVariable
       ncf_copyFile%lastVariable%prev%next=>v
       ncf_copyFile%lastVariable%prev=>v
       nullify(v)
       fv=>fv%next
    end do
    !
    allocate(ncf_copyFile%firstAttribute,ncf_copyFile%lastAttribute, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,'Unable to allocate.',irc
       return
    end if
    ncf_copyFile%firstAttribute%next=>ncf_copyFile%lastAttribute
    ncf_copyFile%lastAttribute%prev=>ncf_copyFile%firstAttribute
    !
    fa=>f%firstAttribute%next
    do while (.not.associated(fa,target=f%lastAttribute))
       a=>ncf_copyAttribute(fa,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_copyAttribute.',irc
          return
       end if
       a%prev=>ncf_copyFile%lastAttribute%prev
       a%next=>ncf_copyFile%lastAttribute
       ncf_copyFile%lastAttribute%prev%next=>a
       ncf_copyFile%lastAttribute%prev=>a
       nullify(a)
       fa=>fa%next
    end do
    ncf_copyFile%initialised=f%initialised
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end function ncf_copyFile

  !@#  INVENTORY ******************************************************

  subroutine ncf_setLabel(f,lab)
    implicit none
    type (inventory), pointer :: f
    character*(*) :: lab
    integer, external :: length
    if (associated(f)) then
       f%lab10=lab
       call chop0(f%lab10,10)
       f%lenl=length(f%lab10,10,5)
    end if
    return
  end subroutine ncf_setLabel
  !
  character*10 function ncf_getLabel(v)
    implicit none
    type (variable), pointer :: v
    character*10 :: b10
    ncf_getLabel="          ";
    if (associated(v)) then
       if (associated(v%f)) then
          b10=v%f%lab10
          call fill0(b10,10)
          ncf_getLabel=b10
       end if
    end if
    return
  end function ncf_getLabel
  !
  subroutine ncf_readInventory(f,biok,irc)
    implicit none
    ! 
    type (inventory) :: f
    logical biok
    integer irc
    ! 
    type(variable), pointer :: v
    type(attribute), pointer :: a
    character*250 buff250
    integer lenb
    integer kk,ll,mm
    integer, external :: length
    integer :: ret
    character*20 :: myname = "ncf_readInventory"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    !
    ! get number of dimension, variables, attributes, unlimited dimension id
    RET = NF_INQ(f%ncid, f%nrdim, f%nrvar, f%nrgatt, f%unlimdimid)
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_INQ:",nf_strerror(ret),f%ncid
       irc=790
       return
    end if
    ! allocate dimensions in file
    allocate(f%dim250(max(1,f%nrdim)),f%lend(max(1,f%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate F-dims.',irc
       return
    end if
    allocate(f%pos,stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate F-pos.',irc
       return
    end if
    f%pos%nrdim=f%nrdim
    allocate(f%pos%sta(max(1,f%pos%nrdim)),&
         & f%pos%lim(max(1,f%pos%nrdim)),&
         & f%pos%siz(max(1,f%pos%nrdim)),&
         & f%pos%pos(max(1,f%pos%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate F-pos-pos.',irc
       return
    end if
    f%pos%pos(1)=1
    f%pos%lim(1)=1
    f%pos%siz(1)=1
    f%pos%sta(1)=1
    do ll=1,f%pos%nrdim
       f%pos%pos(ll)=1
       f%pos%lim(ll)=1
       f%pos%siz(ll)=1
       f%pos%sta(ll)=1
    end do
    ! 
    ! -> store dimension names
    do kk=1,f%nrdim
       call ncf_readDimension(f%ncid,kk,f%dim250(kk),f%lend(kk),&
            & f%pos%siz(kk),irc)
       if (irc.ne.0) then
          write(*,*) 'READINVENTORY Error return from READDIMENSION.',irc
          return
       end if
       f%pos%sta(kk)=1
       f%pos%lim(kk)=f%pos%siz(kk)
    end do
    !
    allocate(f%firstVariable,f%lastVariable, stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate lastvariable.',irc
       return
    end if
    f%firstVariable%next=>f%lastVariable
    f%lastVariable%prev=>f%firstVariable
    ! -> store variables
    do kk=1,f%nrvar
       allocate(v,stat=irc)
       if (irc.ne.0) then
          write(*,*) myname,' Unable to allocate v.',irc
          return
       end if
       call ncf_readVariableAttributes(f%ncid,kk,v,irc)
       if (irc.ne.0) then
          write(*,*) 'READINVENTORY Error return from READVARIABLEATTRIBUTES.',irc
          return
       end if
       call ncf_setVariableDimInv(f,v,irc)
       if (irc.ne.0) then
          write(*,*) 'READINVENTORY Error return from SetVariableDimension.',irc
          return
       end if
       if (ncf_bdeb)write(*,*)myname,"Var:'"//v%var250(1:v%lenv)//"'",v%nrdim,associated(v%ind),v%len
       v%prev=>f%lastVariable%prev
       v%next=>f%lastVariable
       f%lastVariable%prev%next=>v
       f%lastVariable%prev=>v
       nullify(v)
    end do
    !
    allocate(f%firstAttribute,f%lastAttribute, stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate firstattribute.',irc
       return
    end if
    f%firstAttribute%next=>f%lastAttribute
    f%lastAttribute%prev=>f%firstAttribute
    !
    ! -> store global attributes
    do kk=1,f%nrgatt
       allocate(a,stat=irc)
       if (irc.ne.0) then
          write(*,*) myname,' Unable to allocate a.',irc
          return
       end if
       call ncf_readAttribute(f%ncid,nf_global,kk,a,irc)
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
    f%initialised=.true.
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_readInventory

  subroutine ncf_checkInventory(f,irc)
    implicit none
    ! 
    type (inventory), target :: f
    integer irc
    ! 
    type(variable), pointer :: v
    type(attribute), pointer :: a
    character*250 buff250
    integer lenb
    integer kk,ll,mm
    integer, external :: length
    integer :: ret
    character*20 :: myname = "ncf_checkInventory"
    write(*,*)myname,' Entering. '
    if (f%lenf.ne.0) then
       write(*,*) myname,' File: ',f%fn250(1:f%lenf)
       if (associated(f%firstVariable).and.associated(f%lastVariable)) then
          v=>f%firstVariable%next
          do while (associated(v).and..not.associated(v,target=f%lastVariable))
             if (v%lenv.le.0) then
                write(*,*)myname,'*** Zero variable name found.',v%nrdim
             elseif (v%lenv.gt.250) then
                write(*,*)myname,'*** Infinite variable name found.',v%nrdim,v%lend
             else
                write(*,*)myname,v%lenv,v%type
                write(*,*)myname,v%var250(1:min(250,v%lenv))
                write(*,*)myname,ncf_gettype(v%type)
                write(*,*)myname,v%nrdim
                write(*,*)myname,v%len
                write(*,*)myname,v%lend
                write(*,*)myname,v%readyfield
                write(*,*)myname,associated(v%fd)

                write(*,'(X,A,A,A,A,I0,A,I0,X,I0,2(A,L1))')&
                     & myname,' Checking variable: ',&
                     & v%var250(1:min(250,v%lenv))//" "//ncf_gettype(v%type),&
                     & " ndim=",v%nrdim," len=",v%len,v%lend,&
                     & " rf=",v%readyfield
                write(*,'(X,A,A,L1)') v%var250(1:min(250,v%lenv))//" "//ncf_gettype(v%type),&
                     & " a(fd)=",associated(v%fd)
             end if
             if (v%nrdim.ne.0.and..not.associated(v%ind)) then
                write(*,*)myname,'*** Variable with non-associated dimensions:',&
                     & v%var250(1:min(250,v%lenv)),&
                     & v%nrdim,associated(v%ind),associated(v%sta),associated(v%lim)
             end if
             if (.not.associated(v%f,target=f)) then
                write(*,*)myname,'*** Coocoo child found: ',v%var250(1:min(250,v%lenv))
             end if
             if (.not.associated(v%prev).or..not.associated(v%next)) then
                write(*,*)myname,'*** Unlinked variable chain: ',&
                     & "'"//v%var250(1:min(250,v%lenv))//"'", &
                     & associated(v%prev),associated(v%next)
             end if
             v=>v%next
          end do
       else
          write(*,*)myname,'*** Unlinked Inventory variable chain: ', &
               & associated(f%firstVariable),associated(f%lastVariable)
       end if
    else
       write(*,*)myname,' Undefined inventory.'
    end if
    ! reset indexes...
    nullify(f%latid)
    nullify(f%lonid)
    nullify(f%timid)
    nullify(f%terid)
    nullify(f%parid)
    nullify(f%xcoid)
    nullify(f%ycoid)
    nullify(f%gridid)
    write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_checkInventory

  function ncf_copyInventory(f,irc)
    implicit none
    type(inventory), pointer :: ncf_copyInventory
    type(inventory) :: f
    integer :: irc
    integer :: ii
    type(variable),pointer :: v,fv
    type(attribute),pointer :: a,fa
    character*20 :: myname = "ncf_copeyInventory"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(ncf_copyInventory,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate ncf_copyInventory.'
       return
    end if
    ncf_copyInventory%ncid=f%ncid
    ncf_copyInventory%fn250=f%fn250
    ncf_copyInventory%lenf=f%lenf
    ncf_copyInventory%nrdim=f%nrdim
    ncf_copyInventory%nrvar=f%nrvar
    ncf_copyInventory%nrgatt=f%nrgatt
    allocate(ncf_copyInventory%dim250(max(1,ncf_copyInventory%nrdim)),&
         & ncf_copyInventory%lend(max(1,ncf_copyInventory%nrdim)),&
         & ncf_copyInventory%pos,stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate.',irc
       return
    end if
    allocate(ncf_copyInventory%pos%pos(max(1,ncf_copyInventory%nrdim)),&
         & ncf_copyInventory%pos%sta(max(1,ncf_copyInventory%nrdim)),&
         & ncf_copyInventory%pos%lim(max(1,ncf_copyInventory%nrdim)),&
         & ncf_copyInventory%pos%siz(max(1,ncf_copyInventory%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate.',irc
       return
    end if
    do ii=1,ncf_copyInventory%nrdim
       ncf_copyInventory%dim250(ii)=f%dim250(ii)
       ncf_copyInventory%lend(ii)=f%lend(ii)
    end do
    ncf_copyInventory%pos%nrdim=f%pos%nrdim
    do ii=1,ncf_copyInventory%pos%nrdim
       ncf_copyInventory%pos%pos(ii)=f%pos%pos(ii)
       ncf_copyInventory%pos%sta(ii)=f%pos%sta(ii)
       ncf_copyInventory%pos%lim(ii)=f%pos%lim(ii)
       ncf_copyInventory%pos%siz(ii)=f%pos%siz(ii)
    end do
    ncf_copyInventory%pos%loc=f%pos%loc
    ncf_copyInventory%unlimdimid=f%unlimdimid
    !
    allocate(ncf_copyInventory%firstVariable,ncf_copyInventory%lastVariable, stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate.',irc
       return
    end if
    ncf_copyInventory%firstVariable%next=>ncf_copyInventory%lastVariable
    ncf_copyInventory%lastVariable%prev=>ncf_copyInventory%firstVariable
    !
    fv=>f%firstVariable%next
    do while (.not.associated(fv,target=f%lastVariable))
       v=>ncf_copyVariable(fv,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_copyVariable.',irc
          return
       end if
       v%f=>ncf_copyInventory
       if (associated(fv,target=f%latid)) then
          ncf_copyInventory%latid=>v
       else if (associated(fv,target=f%lonid)) then
          ncf_copyInventory%lonid=>v
       else if (associated(fv,target=f%timid)) then
          ncf_copyInventory%timid=>v
       else if (associated(fv,target=f%terid)) then
          ncf_copyInventory%terid=>v
       else if (associated(fv,target=f%parid)) then
          ncf_copyInventory%parid=>v
       else if (associated(fv,target=f%xcoid)) then
          ncf_copyInventory%xcoid=>v
       else if (associated(fv,target=f%ycoid)) then
          ncf_copyInventory%ycoid=>v
       else if (associated(fv,target=f%gridid)) then
          ncf_copyInventory%gridid=>v
       end if
       v%prev=>ncf_copyInventory%lastVariable%prev
       v%next=>ncf_copyInventory%lastVariable
       ncf_copyInventory%lastVariable%prev%next=>v
       ncf_copyInventory%lastVariable%prev=>v
       nullify(v)
       fv=>fv%next
    end do
    !
    allocate(ncf_copyInventory%firstAttribute,ncf_copyInventory%lastAttribute, stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate.',irc
       return
    end if
    ncf_copyInventory%firstAttribute%next=>ncf_copyInventory%lastAttribute
    ncf_copyInventory%lastAttribute%prev=>ncf_copyInventory%firstAttribute
    !
    fa=>f%firstAttribute%next
    do while (.not.associated(fa,target=f%lastAttribute))
       a=>ncf_copyAttribute(fa,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_copyAttribute.',irc
          return
       end if
       a%prev=>ncf_copyInventory%lastAttribute%prev
       a%next=>ncf_copyInventory%lastAttribute
       ncf_copyInventory%lastAttribute%prev%next=>a
       ncf_copyInventory%lastAttribute%prev=>a
       nullify(a)
       fa=>fa%next
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_copyInventory

  subroutine ncf_copyGrid(o,f,irc)
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
    allocate(o%dim250(max(1,o%nrdim)),&
         & o%lend(max(1,o%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate.',irc
       return
    end if
    do ii=1,o%nrdim
       o%dim250(ii)=f%dim250(ii)
       o%lend(ii)=f%lend(ii)
    end do
    allocate(o%pos,stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate pos.',irc
       return
    end if
    o%pos%nrdim=f%pos%nrdim
    allocate(o%pos%pos(max(1,o%pos%nrdim)),&
         & o%pos%sta(max(1,o%pos%nrdim)),&
         & o%pos%lim(max(1,o%pos%nrdim)),&
         & o%pos%siz(max(1,o%pos%nrdim)),&
         & stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate pos%pos.',irc
       return
    end if
    do ii=1,o%pos%nrdim
       o%pos%pos(ii)=f%pos%pos(ii)
       o%pos%sta(ii)=f%pos%sta(ii)
       o%pos%lim(ii)=f%pos%lim(ii)
       o%pos%siz(ii)=f%pos%siz(ii)
    end do
    o%unlimdimid=f%unlimdimid 
    o%pos%loc=f%pos%loc
    !
    allocate(o%firstVariable,o%lastVariable, stat=irc)
    if (irc.ne.0) then
       write(*,*) myname,' Unable to allocate firstvariable.',irc
       return
    end if
    o%firstVariable%next=>o%lastVariable
    o%lastVariable%prev=>o%firstVariable
    !
    o%nrvar=0
    sv=>f%firstVariable%next
    do while (.not.associated(sv,target=f%lastVariable))
       used=.true.
       if (associated(sv,target=f%latid)) then
          o%nrvar=o%nrvar+1
          v=>ncf_copyVariable(sv,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyVariable.',irc
             return
          end if
          o%latid=>v
       else if (associated(sv,target=f%lonid)) then
          o%nrvar=o%nrvar+1
          v=>ncf_copyVariable(sv,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyVariable.',irc
             return
          end if
          o%lonid=>v
       else if (associated(sv,target=f%xcoid)) then
          o%nrvar=o%nrvar+1
          v=>ncf_copyVariable(sv,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyVariable.',irc
             return
          end if
          o%xcoid=>v
       else if (associated(sv,target=f%ycoid)) then
          o%nrvar=o%nrvar+1
          v=>ncf_copyVariable(sv,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyVariable.',irc
             return
          end if
          o%ycoid=>v
       else if (associated(sv,target=f%gridid)) then
          o%nrvar=o%nrvar+1
          v=>ncf_copyVariable(sv,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyVariable.',irc
             return
          end if
          o%gridid=>v
       else
          used=.false.
       end if
       if (used) then
          ! write(*,*) myname,'Copying variable Field:',v%var250(1:v%lenv)
          call ncf_copyVariableField(v,sv,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyVariableFied.',irc
             return
          end if
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
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate firstattribute.'
       return
    end if
    o%firstAttribute%next=>o%lastAttribute
    o%lastAttribute%prev=>o%firstAttribute
    !
    sa=>f%firstAttribute%next
    do while (.not.associated(sa,f%lastAttribute))
       a=>ncf_copyAttribute(sa,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_copyAttribute.',irc
          return
       end if
       a%prev=>o%lastAttribute%prev
       a%next=>o%lastAttribute
       o%lastAttribute%prev%next=>a
       o%lastAttribute%prev=>a
       nullify(a)
       sa=>sa%next
    end do
    o%initialised=f%initialised
  end subroutine ncf_copyGrid

  ! remove all data
  subroutine ncf_clearVariables(f,irc)
    implicit none
    type (inventory),pointer :: f
    integer irc
    type(variable),pointer :: v,vx
    character*20 :: myname = "ncf_clearVariables"
    if (f%initialised) then
       write(*,*)myname,' Clearing variables.' ! if (ncf_bdeb) 
       v=>f%firstVariable%next
       do while(.not.associated(v,target=f%lastVariable))
          write(*,*)myname,"Clearing:'"//v%var250(1:v%lenv)//"'"
          call ncf_clearVariable(v,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_clearVariable.',irc
             return
          end if
          v=>v%next
       end do
       nullify(v)
    end if
    return
  end subroutine ncf_clearVariables

  ! remove data, attributes and dimension references
  subroutine ncf_cleanVariables(f,irc)
    implicit none
    type (inventory),pointer :: f
    integer irc
    type(variable),pointer :: v,vx
    character*20 :: myname = "ncf_cleanVariables"
    if (f%initialised) then
       if (ncf_bdeb) write(*,*)myname,' Removing variables.'
       v=>f%firstVariable%next
       do while(.not.associated(v,target=f%lastVariable))
          call ncf_cleanVariable(v,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_cleanVariable.',irc
             return
          end if
          v=>v%next
       end do
       nullify(v)
    end if
    return
  end subroutine ncf_cleanVariables

  ! remove variable completely, also remove it from inventory
  subroutine ncf_removeVariables(f,irc)
    implicit none
    type (inventory),pointer :: f
    integer irc
    type(variable),pointer :: v,vx
    character*20 :: myname = "ncf_removeVariables"
    if (f%initialised) then
       if (ncf_bdeb) write(*,*)myname,' Removing variables.'
       v=>f%firstVariable%next
       do while(.not.associated(v,target=f%lastVariable))
          vx=>v%next
          call ncf_removeVariable(v,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_removeVariable.',irc
             return
          end if
          v=>vx
       end do
       nullify(v)
       nullify(vx)
    end if
    return
  end subroutine ncf_removeVariables

  subroutine ncf_clearInventory(f,irc)
    implicit none
    type (inventory),pointer :: f
    integer irc
    ! 
    type(variable),pointer :: v,vx
    type(attribute),pointer :: a,ax
    character*20 :: myname = "ncf_clearInventory"
    if (ncf_bdeb) write(*,*)myname,' Entering.',associated(f)
    if (associated(f)) then
       if (f%initialised) then
          call ncf_removeVariables(f,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_removeVariables.',irc
             return
          end if
          if(associated(f%firstVariable))deallocate(f%firstVariable,stat=irc)
          if(associated(f%lastVariable))deallocate(f%lastVariable,stat=irc)
          if (ncf_bdeb) write(*,*)myname,' Removing attributes.'
          a=>f%firstAttribute%next
          do while(.not.associated(a,target=f%lastAttribute))
             ax=>a%next
             call ncf_removeAttribute(a,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_removeAttribute.',irc
                return
             end if
             a=>ax
          end do
          nullify(a)
          nullify(ax)
          if(associated(f%firstAttribute)) deallocate(f%firstAttribute,stat=irc)
          if(associated(f%lastAttribute)) deallocate(f%lastAttribute,stat=irc)
          if (ncf_bdeb) write(*,*)myname,' Removing rest.'
          if (associated(f%dim250)) deallocate(f%dim250,stat=irc)
          if (associated(f%lend)) deallocate(f%lend,stat=irc)
          if (associated(f%dimid)) deallocate(f%dimid,stat=irc)
          call ncf_deallocatePos(f%pos,irc)
          if (associated(f%pos)) deallocate(f%pos,stat=irc)
          f%initialised=.false.
          f%opened=.false.
       end if
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_clearInventory

  subroutine ncf_removeVariablesWithPrefix(pre,f,irc)
    implicit none
    character*(*):: pre
    type (inventory),pointer :: f
    integer irc
    ! 
    type(variable),pointer :: v,vx
    type(attribute),pointer :: a,ax
    character*20 :: myname = "ncf_removeVariables..."
    if (ncf_bdeb) write(*,*)myname,' Entering.',associated(f)
    if (associated(f)) then
       if (f%initialised) then
          if (ncf_bdeb) write(*,*)myname,' Removing variables.'
          v=>f%firstVariable%next
          do while(.not.associated(v,target=f%lastVariable))
             vx=>v%next
             if (v%var250(1:len(pre)).eq.pre(1:len(pre))) then
                write(*,*)myname,' Removing:',v%var250(1:v%lenv)
                call ncf_removeVariable(v,irc)
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_removeVariable.',irc
                   return
                end if
             end if
             v=>vx
          end do
          nullify(v)
          nullify(vx)
       end if
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_removeVariablesWithPrefix

  subroutine ncf_clearInventoryAttributes(f,irc)
    implicit none
    type (inventory) f
    integer irc
    ! 
    type(variable),pointer :: v,vx
    type(attribute),pointer :: a,ax
    character*20 :: myname = "ncf_clearInventoryAttributes"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    !
    if (f%initialised) then
       v=>f%firstVariable%next
       do while(.not.associated(v,target=f%lastVariable))
          vx=>v%next
          call ncf_clearAttributes(v,irc)
          if (irc.ne.0) then
             write(*,*) myname,"Error return from ncf_clearAttributes.",irc
             return
          end if
          v=>vx
       end do
       a=>f%firstAttribute%next
       do while(.not.associated(a,target=f%lastAttribute))
          ax=>a%next
          call ncf_removeAttribute(a,irc)
          if (irc.ne.0) then
             write(*,*) myname,"Error return from ncf_removeAttribute.",irc
             return
          end if
          a=>ax
       end do
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_clearInventoryAttributes

  integer function ncf_maxVarid(f,irc)
    implicit none
    type(inventory) :: f
    integer :: irc
    type(variable),pointer :: v,fv
    character*20 :: myname = "ncf_maxVarid"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    !
    ncf_maxVarid=0
    fv=>f%firstVariable%next
    do while (.not.associated(fv,target=f%lastVariable))
       ncf_maxVarid=max(ncf_maxVarid,fv%varid)
       fv=>fv%next
    end do
  end function ncf_maxVarid
  !@#  DIMENSIONS ******************************************************

  subroutine ncf_readDimension(ncid,dimid,dim250,lend,siz,irc)
    ! reads dimension data from netcdf file
    implicit none
    integer ncid
    integer dimid
    character*250 dim250
    integer lend,lenb
    integer, external :: length
    integer siz
    integer irc
    integer ret
    character*250 buff250
    character*20 :: myname = "ncf_readDimension"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    RET = NF_INQ_DIM(ncid,dimid,dim250,siz)
    if (ret .ne. NF_NOERR) then
       write(*,*) "READDIMENSION ERROR from NF_INQ_DIM:",nf_strerror(ret)
       irc=801
       return
    end if
    call chop0(dim250,250)
    lend=length(dim250,250,10)
    write(buff250,'(X,A20,5X,X,I5)')dim250(1:lend),siz
    call chop0(buff250,250)
    lenb=length(buff250,250,10)
    ! write(*,'(10X,A)') buff250(1:lenb)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_readDimension

  subroutine ncf_setDimensionValue(i,d,ival)
    ! copies dimension start/stop from invetory to variable
    implicit none
    type(inventory), target :: i
    type(dimension), target :: d
    integer :: ival
    i%pos%pos(d%ind)=ival
  end subroutine ncf_setDimensionValue

  function ncf_getDimensionValue(i,d)
    ! copies dimension start/stop from invetory to variable
    implicit none
    integer :: ncf_getDimensionValue
    type(inventory), target :: i
    type(dimension), target :: d
    integer :: ival
    ncf_getDimensionValue=i%pos%pos(d%ind)
    return
  end function ncf_getDimensionValue

  subroutine ncf_setVariableDimInv(f,v,irc)
    ! copies dimension start/stop from invetory to variable
    implicit none
    type(inventory), target :: f
    type(variable) :: v
    integer irc
    integer jj
    character*20 :: myname = "ncf_setVariableDimInv"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    v%f=>f
    if (associated(v%lim)) deallocate(v%lim,stat=irc)
    allocate(v%lim(max(1,v%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate v%lim.'
       return
    end if
    if (associated(v%sta)) deallocate(v%sta,stat=irc)
    allocate(v%sta(max(1,v%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate v%sta.'
       return
    end if
    if (associated(v%siz)) deallocate(v%siz,stat=irc)
    allocate(v%siz(max(1,v%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate v%siz.'
       return
    end if
    call ncf_setVariableDim(f,v)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc,v%len
  end subroutine ncf_setVariableDimInv

  logical function ncf_hasDimension(v,d)
    implicit none
    type(variable),pointer :: v
    type(dimension),pointer :: d
    integer :: jj
    if (associated(v).and.associated(d)) then
       if (v%nrdim.gt.0) then
          do jj=1,v%nrdim
             if (v%ind(jj).eq.d%ind) then
                ncf_hasDimension=.true.
                return
             end if
          end do
       end if
    end if
    ncf_hasDimension=.false.
    return
  end function ncf_hasDimension

  subroutine ncf_setVariableDim(f,v)
    implicit none
    type(inventory), target :: f
    type(variable) :: v
    integer :: jj
    v%len=1
    v%lens=1
    if (v%nrdim.gt.0) then
       v%sta(1)=1
       v%lim(1)=1
       v%siz(1)=1
       do jj=1,v%nrdim
          v%sta(jj)=f%pos%sta(v%ind(jj))
          v%lim(jj)=f%pos%lim(v%ind(jj))
          v%siz(jj)=f%pos%siz(v%ind(jj))
          v%len=v%len*v%lim(jj)
          v%lens=v%lens*v%siz(jj)
       end do
    end if
    return
  end subroutine ncf_setVariableDim
  
  subroutine ncf_setInventoryDimOrder(f,o,irc)
    ! sets variable dimensions from index
    implicit none
    type(inventory), pointer ::f
    type(dimensionOrder),pointer :: o
    integer :: irc
    character*20 :: myname = "ncf_setInventoryDimOrder"
    o%f => f
    return
  end subroutine ncf_setInventoryDimOrder

  subroutine ncf_setVariableDimOrder(v,o,irc)
    ! sets variable dimensions from index
    implicit none
    type(variable), pointer ::v
    type(dimensionOrder),pointer :: o
    integer :: irc
    type(dimension),pointer :: d
    integer :: ii
    character*20 :: myname = "ncf_setVariableDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    !
    if (associated(v%ind)) deallocate(v%ind)
    if (associated(v%sta)) deallocate(v%sta)
    if (associated(v%lim)) deallocate(v%lim)
    if (associated(v%siz)) deallocate(v%siz)
    if (associated(o%f))v%f=>o%f
    v%nrdim=o%nrdim
    allocate(v%ind(max(1,v%nrdim)),v%sta(max(1,v%nrdim)),&
         & v%lim(max(1,v%nrdim)),v%siz(max(1,v%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate variable dimensinos.',v%nrdim
       return
    end if
    v%len=1
    v%lens=1
    d=>o%firstDimension%next
    do ii=1,o%nrdim
       if (.not.associated(d,target=o%lastDimension)) then
          v%ind(ii)=d%ind
          v%sta(ii)=d%sta
          v%lim(ii)=d%lim
          v%siz(ii)=d%siz
          v%len=v%len*v%lim(ii)
          v%lens=v%lens*v%siz(ii)
          d=>d%next
       else
          d=>o%firstDimension%next
          do while (.not.associated(d,target=o%lastDimension))
             write(*,*)myname,' Dim:',d%ind, d%sta, d%lim
             d=>d%next
          end do
          write(*,*)myname,"system error.",ii,o%nrdim
          irc=845
          return
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_setVariableDimOrder

  function ncf_createDimension(f,dim,lim)
    implicit none
    ! creates new dimension, creates dimIndex, updates inventory
    type(dimension),pointer :: ncf_createDimension
    type(inventory) :: f
    character*(*) dim
    integer lim
    integer irc
    character*250 xdim250(f%nrdim+1)
    integer xlend(f%nrdim+1), xpos(f%nrdim+1), xsta(f%nrdim+1)
    integer xlim(f%nrdim+1), xsiz(f%nrdim+1)
    integer ii
    character*20 :: myname = "ncf_createDimension"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    do ii=1,f%nrdim
       xdim250(ii)=f%dim250(ii)
       xlend(ii)=f%lend(ii)
       xpos(ii)=f%pos%pos(ii)
       xsta(ii)=f%pos%sta(ii)
       xlim(ii)=f%pos%lim(ii)
       xsiz(ii)=f%pos%siz(ii)
    end do
    if (associated(f%dim250)) deallocate(f%dim250,stat=irc)
    if (associated(f%lend)) deallocate(f%lend,stat=irc)
    if (associated(f%pos%pos)) deallocate(f%pos%pos,stat=irc)
    if (associated(f%pos%sta)) deallocate(f%pos%sta,stat=irc)
    if (associated(f%pos%lim)) deallocate(f%pos%lim,stat=irc)
    if (associated(f%pos%siz)) deallocate(f%pos%siz,stat=irc)
    f%nrdim=f%nrdim+1
    xdim250(f%nrdim)=dim
    xlend(f%nrdim)=len(dim)
    xpos(f%nrdim)=1
    xsta(f%nrdim)=1
    xlim(f%nrdim)=lim
    xsiz(f%nrdim)=lim
    f%pos%nrdim=f%nrdim
    f%pos%maxdim=f%nrdim
    allocate(f%dim250(max(1,f%nrdim)),&
         & f%lend(max(1,f%nrdim)),&
         & f%pos%pos(max(1,f%pos%maxdim)),&
         & f%pos%sta(max(1,f%pos%maxdim)),&
         & f%pos%lim(max(1,f%pos%maxdim)),&
         & f%pos%siz(max(1,f%pos%maxdim)),&
         & stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    do ii=1,f%nrdim
       f%dim250(ii)  = xdim250(ii)
       f%lend(ii)    = xlend(ii)
       f%pos%pos(ii) = xpos(ii)
       f%pos%sta(ii) = xsta(ii)
       f%pos%lim(ii) = xlim(ii)
       f%pos%siz(ii) = xsiz(ii)
    end do
    allocate(ncf_createDimension,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_createDimension%ind=f%nrdim
    ncf_createDimension%sta=f%pos%sta(f%nrdim)
    ncf_createDimension%lim=f%pos%lim(f%nrdim)
    ncf_createDimension%siz=f%pos%siz(f%nrdim)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_createDimension

  subroutine ncf_increaseDimension(f,d,lim,irc)
    ! remove dimensions that are not used from inventory...
    implicit none
    type(inventory), pointer :: f
    type(dimension), pointer :: d
    integer :: lim
    integer :: irc 
    integer indDim,indVar,lv,lvx
    type(variable), pointer :: v, vn, vx
    type(dimensionOrder), pointer :: dv,dvx
    character*20 :: myname = "ncf_increaseDimension"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(f)) then
       indDim=d%ind
       if (indDim.ne.0) then
          ! loop over variables
          v=>f%firstVariable%next
          do while (.not.associated(v,f%lastVariable))
             vn=>v%next 
             ! check if variable uses dimension
             indVar=ncf_getVariableDimEntry(v,indDim)
             if (indVar.ne.0) then  ! create new variable with new dimension limit
                vx=>ncf_copyVariable(v,irc) ! create new variable
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_copyVariable.',irc
                   return
                end if
                vx%lim(indVar)=lim ! set new dimension limit
                vx%siz(indVar)=lim ! set new dimension limit
                call ncf_initField(vx,vx%filld,irc)
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_initField.',irc
                   return
                end if
                ! make dimension orders
                dv=>ncf_makeDimOrder(v,irc)
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
                   return
                end if
                dvx=>ncf_makeDimOrder(vx,irc)
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
                   return
                end if
                ! loop over old dimension index
                call ncf_resetPos(dv,irc)
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_resetPos.',irc
                   return
                end if
                do while (ncf_increment(f,dv,irc))
                   lv = ncf_getLocation(v)
                   lvx = ncf_getLocation(vx)
                   vx%fd(lvx)=v%fd(lv) ! copy array
                end do
                if (irc.ne.0) then
                   write(*,*)myname,'Error return from ncf_incrementDimension.',irc
                   return
                end if
                ! replace old variable by new variable
                call ncf_replaceVariable(v,vx)
             end if
             v=>vn
          end do
          ! increase dimension value
          if (ncf_bdeb) write(*,*)myname,' Increasing dimension: ',&
               & f%dim250(indDim)(1:f%lend(inddim)), f%pos%lim(indDim),' -> ',lim
          f%pos%lim(indDim)=lim
          f%pos%siz(indDim)=max(lim,f%pos%siz(indDim))
       end if
    end if
    return
  end subroutine ncf_increaseDimension

  subroutine ncf_compressVariables(f,irc)
    ! remove dimensions that are not used from inventory...
    implicit none
    type(inventory), pointer :: f
    integer :: irc
    type(variable), pointer :: v, vn
    integer ii, jj, nrdim, unlimdimid, indx(f%nrdim), cnt(f%nrdim)
    character*20 :: myname = "ncf_compressVariables"
    if (ncf_bdeb) write(*,*)myname,' Entering.',associated(f)
    if (associated(f)) then
       if (ncf_bdeb) write(*,*)myname,'  Init loop.',associated(f%firstVariable),associated(f%lastVariable),&
            & associated(f%firstVariable%next)
       v=>f%firstVariable%next
       do while (.not.associated(v,f%lastVariable))
          if (ncf_bdeb) write(*,*)myname,'  Staring Loop.',associated(v)
          vn=>v%next
          if (ncf_bdeb) write(*,*)myname,'  Looping.',associated(vn)
          if (.not.v%readyfield) then ! no data available
             if (ncf_bdeb) write(*,'(3(X,A),10(X,I0))')myname,&
                  & ' Removing variable:',v%var250(1:v%lenv),v%lenc,&
                  & v%len1,v%len2,v%len4,v%len8,v%lenr,v%lend,v%len,v%lens
             call ncf_ignoreVariable(v)
             call ncf_clearVariable(v,irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_clearVariable.',irc
                return
             end if
          else if (v%lend.ne.0) then
             call ncf_unmakeRealData(v,irc)
             if (irc.ne.0) then
                write(*,*)myname,' Error return from ncf_unmakeRealData.',irc
                return
             end if
          end if
          v=>vn
       end do
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_compressVariables

  subroutine ncf_compressDimensions(f)
    ! remove dimensions that are not used from inventory...
    implicit none
    type(inventory), pointer :: f
    type(variable), pointer :: v
    integer ii, jj, nrdim, unlimdimid, indx(f%nrdim), cnt(f%nrdim)
    character*20 :: myname = "ncf_compressDimension"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(f)) then
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
             f%pos%siz(indx(jj))=f%pos%siz(jj) ! indx(jj) <= jj
          end if
       end do
       f%pos%loc=0
       f%pos%nrdim=nrdim
       f%nrdim=nrdim
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_compressDimensions

  !@#  DIMENSION INDEXES ******************************************************
  !
  subroutine ncf_importDimOrder(i,o,irc)
    type(inventory), pointer :: i
    type(dimensionOrder), pointer :: o
    integer :: irc
    type(dimension), pointer :: d
    integer :: ind
    character*20 :: myname = "ncf_importDimOrder"
    character*50 :: o50,n50
    d => o%firstdimension%next
    do while (.not.associated(d,target=o%lastdimension))
       write(o50,'(I0,"<",I0,",",I0,"|",I0,"> ",A)') d%ind,d%sta,d%lim,d%siz,&
            & "     "
       ind=ncf_getDimEntry(i,o%f%dim250(d%ind)(1:o%f%lend(d%ind)))
       if (ind.eq.0) then
          write(*,*)myname,'Missing dimension:',o%f%dim250(d%ind)(1:o%f%lend(d%ind))
          irc=477
          return
       end if
       d%ind=ind
       d%sta=i%pos%sta(ind)
       d%lim=i%pos%lim(ind)
       d%siz=i%pos%siz(ind)
       write(n50,'(I0,"<",I0,",",I0,"|",I0,"> ",A)') d%ind,d%sta,d%lim,d%siz,&
            & i%dim250(d%ind)(1:i%lend(d%ind))
       !write(*,*) myname,"Importing: ",trim(o50)//"->"//trim(n50)
       d => d%next
    end do
    o%f=>i
    return
  end subroutine ncf_importDimOrder
  !
  subroutine ncf_importVariable(i,v,irc)
    type(inventory), pointer :: i
    type(variable), pointer :: v
    integer :: irc
    integer :: ii,ind
    character*20 :: myname = "ncf_importVariable"
    character*50 :: o50,n50
    ! import dimensions
    if (associated(v)) then
       do ii=1,v%nrdim
          write(o50,'(I0,"<",I0,",",I0,"|",I0,"> ",A)') v%ind(ii),v%sta(ii),v%lim(ii),&
               & v%siz(ii),"     "
          ind=ncf_getDimEntry(i,v%f%dim250(v%ind(ii))(1:v%f%lend(v%ind(ii))))
          if (ind.eq.0) then
             call ncf_printDimensions(i)
             call ncf_printvariable(v)
             write(*,*)myname,'Missing dimension:',&
                  & v%f%dim250(v%ind(ii))(1:v%f%lend(v%ind(ii)))
             irc=477
             return
          end if
          v%ind(ii)=ind
          v%sta(ii)=i%pos%sta(ind)
          v%lim(ii)=i%pos%lim(ind)
          v%siz(ii)=i%pos%siz(ind)
          write(n50,'(I0,"<",I0,",",I0,"|",I0,"> ",A)') v%ind(ii),v%sta(ii),v%lim(ii),&
               & v%siz(ii),i%dim250(ind)(1:i%lend(ind))
          !write(*,*) myname,"Importing: ",trim(o50)//"->"//trim(n50)
       end do
       v%f=>i
    else
       write(*,*)myname,'Undefined variable...'
    end if
    return
  end subroutine ncf_importVariable
  !
  function ncf_newDimOrder(i,irc)
    implicit none
    type(dimensionOrder), pointer :: ncf_newDimOrder
    type(inventory), pointer :: i
    integer irc
    integer kk
    character*20 :: myname = "ncf_newDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(ncf_newDimOrder,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_newDimOrder%f=>i
    allocate(ncf_newDimOrder%firstDimension,ncf_newDimOrder%lastDimension, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_newDimOrder%firstDimension%next=>ncf_newDimOrder%lastDimension
    ncf_newDimOrder%lastDimension%prev=>ncf_newDimOrder%firstDimension
    ncf_newDimOrder%nrdim=0
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_newDimOrder

  subroutine ncf_addDimOrderDim(o,d,irc)
    implicit none
    type(dimensionOrder), pointer :: o
    type(dimension), pointer :: d
    integer irc
    integer kk
    character*20 :: myname = "ncf_addDimOrderDim"
    if (ncf_bdeb) write(*,*)myname,' Entering.',associated(o),associated(d)
    if (.not.associated(d)) then
       write(*,*)myname,'Invalid dimension.'
       irc=945
       return
    else if (.not.associated(o)) then
       write(*,*)myname,'Dimorder does not exist.'
       irc=946
       return
    else
       d%prev=>o%lastDimension%prev
       d%next=>o%lastDimension
       o%lastDimension%prev%next=>d
       o%lastDimension%prev=>d
       o%nrdim=o%nrdim+1
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_addDimOrderDim
  !
  function ncf_mapDimOrder(f,ido,irc)
    implicit none
    type(inventory), pointer :: f
    type(dimensionOrder), pointer :: ido, ncf_mapDimOrder
    type(variable),pointer :: v
    integer irc
    integer kk
    type(dimension), pointer :: d, dx
    character*20 :: myname = "ncf_mapDimOrder"
    type(dimensionOrder), pointer :: dd
    integer :: ii
    dd => ncf_createDimOrder(irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from createDimOrder.',irc
       return
    end if
    call ncf_setDimOrderInventory(dd,f)
    ! loop over ido
    d=>ido%firstDimension%next
    do while (.not.associated(d, target= ido%lastDimension))
       ! get dimension in f
       ii = ncf_getDimEntry(f,ido%f%dim250(d%ind)(1:ido%f%lend(d%ind)))
       if (ii.eq.0) then ! dimension does not exist in output
          if (d%siz.ne.1) then
             write(*,*)myname,"Missing output dimension.",&
                  & ido%f%dim250(d%ind)(1:ido%f%lend(d%ind))
             call ncf_printDimension(ido%f,d)
          else
             write(*,*)myname,"Ignoring dimension:",&
                  & ido%f%dim250(d%ind)(1:ido%f%lend(d%ind))
          end if
       else ! found output dimension
          dx => ncf_makeDimension(f,ii)
          ! sanity check
          if (dx%siz .ne. d%siz) then
             irc=288
             write(*,*)myname,"Dimension mismatch.",&
                  & ido%f%dim250(d%ind)(1:ido%f%lend(d%ind))
             call ncf_printDimension(ido%f,d)
             call ncf_printDimension(f,dx)
             return
          end if
          ! add dimension to dd
          call ncf_addDimOrderDim(dd,dx,irc)
       end if
       d=>d%next
    end do
    ncf_mapDimOrder => dd
    return
  end function ncf_mapDimOrder
  !
  function ncf_makeDimOrder(v,irc)
    implicit none
    type(dimensionOrder), pointer :: dd, ncf_makeDimOrder
    type(variable),pointer :: v
    integer irc
    integer kk
    type(dimension), pointer :: d
    character*20 :: myname = "ncf_makeDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.',associated(v%f)
    dd => ncf_createDimOrder(irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from createDimOrder.',irc
       return
    end if
    if (associated(v)) then
       if(ncf_bdeb)write(*,*)myname,' Variable:',v%var250(1:v%lenv),v%nrdim,v%ind
       dd%nrdim=v%nrdim
       do kk=1,v%nrdim
          allocate(d,stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate dimension: ',kk,v%ind
             return
          end if
          if (associated(d,target=d%prev) .or. &
               & associated(d,target=d%next).or. &
               & associated(d,target=dd%firstDimension).or. &
               & associated(d,target=dd%lastDimension)) then
             write(*,*)myname,'Memory corruption:',&
                  & associated(d,target=d%prev), &
                  & associated(d,target=d%next), &
                  & associated(d,target=dd%firstDimension), &
                  & associated(d,target=dd%lastDimension)
          end if
          d%ind=v%ind(kk)
          d%sta=v%sta(kk)
          d%lim=v%lim(kk)
          d%siz=v%siz(kk)
          d%prev=>dd%lastDimension%prev
          d%next=>dd%lastDimension
          dd%lastDimension%prev%next=>d
          dd%lastDimension%prev=>d
          call ncf_setDimOrderInventory(dd,v%f)
          nullify(d)
       end do
    else
       dd%nrdim=0
    end if
    ncf_makeDimOrder=>dd
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_makeDimOrder
  !
  subroutine ncf_checkDimOrder(v,do,irc)
    implicit none
    type(variable),pointer :: v
    type(dimensionOrder), pointer :: do
    integer irc
    integer kk
    type(dimension), pointer :: d
    character*20 :: myname = "ncf_checkDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(v)) then
       d=>do%firstDimension%next
       do kk=1,v%nrdim
          if (.not.associated(d,target=do%lastDimension)) then
             if (d%ind .ne. v%ind(kk) .or. &
                  d%sta .ne. v%sta(kk) .or. &
                  d%lim .ne. v%lim(kk)) then
                write(*,'(X,A,A,X,I0,3(X,I0,X,"(",I0,")"))') &
                     & myname,' Dimension error:',kk,&
                     & d%ind,v%ind(kk),d%sta,v%sta(kk),d%lim,v%lim(kk)
                irc=348
                d=>d%next
             end if
          end if
       end do
    else
       write(*,*)myname,' Invalid Variable.'
       irc=349
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_checkDimOrder
  !
  subroutine ncf_setDimOrderInventory(d,f)
    type(dimensionOrder), pointer :: d
    type(inventory), pointer :: f
    d%f=>f
    return
  end subroutine ncf_setDimOrderInventory
  !
  function ncf_createDimOrder(irc)
    implicit none
    integer irc
    character*20 :: myname = "ncf_createDimOrder"
    type(dimensionOrder), pointer :: dd=>null(), ncf_createDimOrder
    allocate(dd,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate dimOrder.'
       return
    end if
    allocate(dd%firstDimension,dd%lastDimension, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate first/last-dimensions.'
       return
    end if
    dd%firstDimension%next=>dd%lastDimension
    dd%lastDimension%prev=>dd%firstDimension
    ncf_createDimOrder => dd
    nullify(dd)
    return
  end function ncf_createDimOrder
  !
  function ncf_createVariable(i,name,type,ido,irc)
    implicit none
    type(variable), pointer :: ncf_createVariable
    type(inventory), pointer :: i
    character*(*) :: name
    integer :: type
    type(dimensionorder), pointer :: ido
    integer :: irc
    type(variable),pointer :: v
    type(dimension),pointer :: d
    integer :: lenr,lenrs,ii
    integer, external :: length
    character*20 :: myname = "ncf_createVariable"
    allocate(ncf_createVariable,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,'Unable to allocate.',irc
       return
    end if
    ncf_createVariable%f=>i
    ncf_createVariable%ncid=i%ncid
    v=>i%lastVariable%prev
    if (.not.associated(v,target=i%firstVariable)) then
       ncf_createVariable%varid=1+ncf_maxVarid(i,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_maxVarid.',irc
          return
       end if
    else
       ncf_createVariable%varid=1
    end if
    ncf_createVariable%var250=name
    call chop0(ncf_createVariable%var250,250)
    ncf_createVariable%lenv=length(ncf_createVariable%var250,250,10)
    ncf_createVariable%type=type
    ncf_createVariable%nrdim=ido%nrdim
    ncf_createVariable%nratt=0.0D0
    ncf_createVariable%scale=1.0D0
    ncf_createVariable%fillc=char(nf_fill_char)
    ncf_createVariable%fill1=nf_fill_int1
    ncf_createVariable%fill2=nf_fill_int2
    ncf_createVariable%fill4=nf_fill_int
    ncf_createVariable%fill8=nf_fill_int64
    ncf_createVariable%fillr=nf_fill_real
    ncf_createVariable%filld=nf_fill_double
    allocate(ncf_createVariable%ind(max(1,ncf_createVariable%nrdim)),&
         & ncf_createVariable%sta(max(1,ncf_createVariable%nrdim)),&
         & ncf_createVariable%lim(max(1,ncf_createVariable%nrdim)),&
         & ncf_createVariable%siz(max(1,ncf_createVariable%nrdim)),&
         & stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    if (ncf_bdeb)write(*,*)myname,"Var:'"//name//"'"
    lenr=1
    lenrs=1
    ii=0
    d=>ido%firstDimension%next
    do while (.not.associated(d,target=ido%lastDimension))
       ii=ii+1
       !write(*,*)myname,' Adding dimension:',d%ind, &
       !     & ido%f%dim250(d%ind)(1:ido%f%lend(d%ind)),d%sta,d%lim
       if (ii.gt.ido%nrdim) then
          call ncf_printDimOrder(ido)
          write(*,*)myname,' Corrupt dimension order.',ii,ido%nrdim
          irc=954
          return
       end if
       ncf_createVariable%ind(ii)=d%ind
       ncf_createVariable%sta(ii)=d%sta
       ncf_createVariable%lim(ii)=d%lim
       ncf_createVariable%siz(ii)=d%siz
       lenr=lenr*max(1,d%lim)
       lenrs=lenrs*max(1,d%siz)
       d=>d%next
    end do
    if (ii.ne.ido%nrdim) then
       write(*,*)myname,' CORRUPT dimension order.',ii,ido%nrdim
       irc=955
       return
    end if
    ncf_createVariable%len=lenr
    ncf_createVariable%lens=lenrs
    ncf_createVariable%lenc=0
    ncf_createVariable%len1=0
    ncf_createVariable%len2=0
    ncf_createVariable%len4=0
    ncf_createVariable%len8=0
    ncf_createVariable%lenr=0
    ncf_createVariable%lend=0
    allocate(ncf_createVariable%firstAttribute,ncf_createVariable%lastAttribute, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_createVariable%firstAttribute%next=>ncf_createVariable%lastAttribute
    ncf_createVariable%lastAttribute%prev=>ncf_createVariable%firstAttribute
    ncf_createVariable%initialised=.true.
    ncf_createVariable%readyfield=.false.
    return
  end function ncf_createVariable


  subroutine ncf_setVariableTypeDouble(v,irc)
    implicit none
    type(variable),pointer :: v
    integer :: type
    real :: fill
    integer :: irc
    v%type=nf_double
    call ncf_setDoubleAttribute(v,"_FillValue",nf_fill_double)
    return
  end subroutine ncf_setVariableTypeDouble

  function ncf_getVariable(i,name,bok,irc)
    implicit none
    type(variable), pointer :: ncf_getVariable
    type(inventory), pointer :: i
    character*(*) :: name
    logical :: bok
    integer :: irc
    character*20 :: myname = "ncf_getVariable"
    bok=.true.
    ncf_getVariable=>i%firstVariable%next
    do while (.not.associated(ncf_getVariable,i%lastVariable))
       if (name(1:len(name)).eq.ncf_getVariable%var250(1:ncf_getVariable%lenv)) then
          bok=.true.
          return
       else
          ncf_getVariable=>ncf_getVariable%next
       end if
    end do
    nullify(ncf_getVariable)
    if (ncf_bdeb) then
       ncf_getVariable=>i%firstVariable%next
       do while (.not.associated(ncf_getVariable,i%lastVariable))
          write(*,*)myname,"No match:'"//name(1:len(name))//"' <> '"//ncf_getVariable%var250(1:ncf_getVariable%lenv)//"'",&
               & (name.eq.ncf_getVariable%var250(1:ncf_getVariable%lenv)),len(name),ncf_getVariable%lenv
          ncf_getVariable=>ncf_getVariable%next
       end do
    end if
    nullify(ncf_getVariable)
    if (bok)write(*,*)myname,"Unable to find '"//name//"' ("//i%lab10(1:i%lenl)//")"
    bok=.false.
    return
  end function ncf_getVariable

  subroutine ncf_setVariableName(v,name)
    implicit none
    type(variable),pointer :: v
    character*(*) :: name
    integer, external :: length
    if (associated(v)) then
       v%var250=name
       call chop0(v%var250,250)
       v%lenv=length(v%var250,250,10)
    end if
    return
  end subroutine ncf_setVariableName
  !
  subroutine ncf_setOutVariableName(v,name)
    implicit none
    type(variable),pointer :: v
    character*(*) :: name
    integer, external :: length
    if (associated(v)) then
       v%out250=name
       call chop0(v%out250,250)
       v%leno=length(v%out250,250,10)
       v%lout= (v%leno.ne.0)
    end if
    return
  end subroutine ncf_setOutVariableName
  !
  integer function ncf_getVariableDimEntry(v,iavg)
    implicit none
    type(variable),pointer :: v
    integer :: iavg
    integer :: ii
    character*20 :: myname = "ncf_getVariableDimEntry"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ncf_getVariableDimEntry=0
    do ii=1,v%nrdim
       if (v%ind(ii) .eq. iavg) then
          ncf_getVariableDimEntry=ii
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end function ncf_getVariableDimEntry

  function ncf_getVariableDimension(v,indDim)
    implicit none
    type(dimension),pointer :: ncf_getVariableDimension
    type(variable),pointer :: v
    integer :: indDim
    integer :: ii
    character*20 :: myname = "ncf_getVariableDimension"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    do ii=1,v%nrdim
       if (v%ind(ii) .eq. indDim) then
          allocate(ncf_getVariableDimension)
          ncf_getVariableDimension%ind=indDim
          ncf_getVariableDimension%sta=v%sta(indDim)
          ncf_getVariableDimension%lim=v%lim(indDim)
          ncf_getVariableDimension%siz=v%siz(indDim)
          return
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end function ncf_getVariableDimension

  function ncf_makeDimOrderInv(f,irc)
    implicit none
    ! returns a dimIndex of all dimensions in inventory
    type(dimensionOrder), pointer :: ncf_makeDimOrderInv
    type(inventory) :: f
    integer irc
    integer kk
    type(dimension), pointer :: d
    character*20 :: myname = "ncf_makeDimOrderInv"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(ncf_makeDimOrderInv,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_makeDimOrderInv%nrdim=f%nrdim
    allocate(ncf_makeDimOrderInv%firstDimension,ncf_makeDimOrderInv%lastDimension, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_makeDimOrderInv%firstDimension%next=>ncf_makeDimOrderInv%lastDimension
    ncf_makeDimOrderInv%lastDimension%prev=>ncf_makeDimOrderInv%firstDimension
    do kk=1,f%nrdim
       allocate(d,stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       d%ind=kk
       d%sta=1
       d%lim=f%pos%lim(kk)
       d%siz=f%pos%siz(kk)
       d%prev=>ncf_makeDimOrderInv%lastDimension%prev
       d%next=>ncf_makeDimOrderInv%lastDimension
       ncf_makeDimOrderInv%lastDimension%prev%next=>d
       ncf_makeDimOrderInv%lastDimension%prev=>d
       nullify(d)
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_makeDimOrderInv

  subroutine ncf_makeIndex(ind,f,i,irc)
    implicit none
    ! returns a dimIndex of all dimensions in inventory
    integer, allocatable :: ind(:)
    type(inventory) :: f, i
    integer irc
    integer kk,jj
    type(dimension), pointer :: d
    character*20 :: myname = "ncf_makeIndex"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (allocated(ind)) deallocate(ind)
    allocate(ind(max(1,f%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    do kk=1,f%nrdim
       ind(kk)=0
       do jj=1,i%nrdim
          if (f%dim250(kk)(1:f%lend(kk)).eq.i%dim250(jj)(1:i%lend(jj))) then
             ind(kk)=jj
          end if
       end do
       if (ind(kk).eq.0) then
          write(*,*)myname,'No match for dim :"'//f%dim250(kk)(1:f%lend(kk))//'"'
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_makeIndex
  !
  subroutine ncf_copyIndexPos(ind,f,i,irc)
    integer, allocatable :: ind(:)
    type(inventory),pointer :: f, i
    integer :: irc
    character*20 :: myname = "ncf_copyIndexPos"
    do kk=1,size(ind)
       if (ind(kk).le.0) then
          write(*,*)myname,'Invalid index:',kk,ind(kk)
          irc=847
          return
       else
          i%pos%pos(ind(kk))=f%pos%pos(kk)
       end if
    end do
    return
  end subroutine ncf_copyIndexPos
  !
  subroutine ncf_clearDimOrder(i)
    implicit none
    ! deletes the contents in the dimIndex
    type(dimensionOrder), pointer  :: i
    type(dimension), pointer :: d, nd
    character*20 :: myname = "ncf_clearDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(i)) then
       if (associated(i%firstDimension).and.associated(i%lastDimension)) then
          d=>i%firstDimension%next
          do while (.not.associated(d,i%lastDimension))
             nd=>d%next
             deallocate(d)
             d=>nd
          end do
          deallocate(i%firstDimension,i%lastDimension)
       end if
       deallocate(i)
       nullify(i)
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_clearDimOrder

  subroutine ncf_addDimOrder(t,i,irc)
    ! adds dimensionOrder "i" to "t"
    implicit none
    type(dimensionOrder) :: t
    type(dimensionOrder) :: i
    integer :: irc
    type(dimension),pointer  :: td, id, x
    logical :: found
    character*20 :: myname = "ncf_addDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
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
          if (irc.ne.0) then
             write(*,*)myname,"Unable to allocate dimension.",irc
             return
          end if
          x%ind=id%ind
          x%lim=id%lim
          x%siz=id%siz
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
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_addDimOrder
  !
  subroutine ncf_removeDimOrder(t,i)
    implicit none
    ! removes contents matching index "i" from index "t"
    type(dimensionOrder),pointer :: t
    type(dimensionOrder),pointer :: i
    type(dimension),pointer  :: td, id, x
    integer irc
    character*20 :: myname = "ncf_removeDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (.not.associated(t) .or. .not.associated(i)) then
       write(*,*)myname,'Missing dimensionOrder.'
    elseif (.not.associated(t%f) .or. .not.associated(i%f)) then
       write(*,*)myname,'Missing parents.'
    elseif (.not.associated(t%f,target=i%f)) then
       write(*,*)myname,'Differing parents.',&
            & t%f%lab10(1:t%f%lenl)//" != "//i%f%lab10(1:i%f%lenl)
    else
       id=>i%firstDimension%next
       do while (.not.associated(id,target=i%lastDimension))
          td=>t%firstDimension%next
          do while (.not.associated(td,target=t%lastDimension))
             if (i%f%dim250(id%ind)(1:i%f%lend(id%ind)) .eq. &
                  & t%f%dim250(td%ind)(1:t%f%lend(td%ind))) then
                x=>td
                td=>td%next
                t%nrdim = t%nrdim - 1
                x%prev%next=>x%next
                x%next%prev=>x%prev 
                nullify(x%prev,x%next)
                deallocate(x,stat=irc)
                nullify(x)
             else
                td=>td%next
             end if
          end do
          id=>id%next
       end do
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_removeDimOrder
  !
  subroutine ncf_setSlice(f,d,sta,sto)
    implicit none
    type(inventory), pointer :: f
    type(dimension), pointer :: d
    integer :: sta
    integer :: sto
    type(variable), pointer :: v
    integer :: irc
    character*20 :: myname = "ncf_setSlice"
    if (associated(d)) then
       !call ncf_printinventory(f)
       d%sta=sta
       d%lim=sto-sta+1
       f%pos%sta(d%ind)=sta
       f%pos%lim(d%ind)=sto-sta+1
       write(*,'(X,A,4(A,I0),A)') myname," Slicing '"//f%lab10(1:f%lenl)// &
            & "': dim='"//f%dim250(d%ind)(1:f%lend(d%ind))// &
            & "' (indx=",d%ind,&
            & ") start=",sta," stop=",sto," (size=",d%siz,")"
       ! set variable-dimension limits if not read into memory
       v => f%firstvariable%next
       do while (.not.associated(v,target=f%lastvariable))
          if (v%len.ne.v%lens.or..not. v%readyfield) then
             !call ncf_printvariable(v)
             if (v%readyfield) call ncf_clearVariable(v,irc)
             call ncf_setVariableDim(f,v)
             if (ncf_bdeb) write(*,*)myname,'Sliced:',v%var250(1:v%lenv),&
                  & v%len,v%lens,v%nrdim,v%readyfield
          end if
          v=>v%next
       end do
    end if
    return
  end subroutine ncf_setSlice
  !
  integer function ncf_getDimEntry(i,dim)
    implicit none
    ! get global dimension index
    type(inventory),pointer :: i
    character*(*) :: dim
    integer, external :: length
    integer :: lend,ii
    character*20 :: myname = "ncf_getDimEntry"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ncf_getDimEntry=0
    lend=length(dim,len(dim),1)
    do ii=1,i%nrdim
       if (i%dim250(ii)(1:i%lend(ii)).eq.dim(1:lend)) then
          ncf_getDimEntry=ii
          if(ncf_bdeb)write(*,*) myname,"'"//i%dim250(ii)(1:i%lend(ii))//"' == '"//dim(1:len(dim))//"'"
       end if
    end do
    if (ncf_bdeb .and. ncf_getDimEntry.eq.0) then
       do ii=1,i%nrdim
          write(*,*)myname," Checking:",ii,"'"//i%dim250(ii)(1:i%lend(ii))//"' '"//dim(1:lend)//"'",&
               & (i%dim250(ii)(1:i%lend(ii)).eq.dim),i%lend(ii),len(dim),lend
       end do
    end if
    if (ncf_bdeb) write(*,*)myname,' Done. '//dim//" index = ",ncf_getDimEntry
    return
  end function ncf_getDimEntry
  !
  subroutine ncf_addDimOrderEntry(o,indDim,irc)
    implicit none
    ! deletes dimension entry from dimIndex
    type(dimensionOrder),pointer :: o
    integer :: indDim
    integer :: irc
    type(dimension), pointer :: d,n
    character*20 :: myname = "ncf_addDimOrderEntry"
    logical :: found
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (.not.ncf_dimOrderContainsDim(o,indDim)) then
       allocate(d,stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       d%ind=indDim
       d%sta=o%f%pos%sta(indDim)
       d%lim=o%f%pos%lim(indDim)
       d%siz=o%f%pos%siz(indDim)
       call ncf_addDimOrderDim(o,d,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_addDimOrderDim.',irc
          return
       end if
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end subroutine ncf_addDimOrderEntry
  !
  function ncf_removeDimOrderEntry(i,indDim)
    implicit none
    type(dimension), pointer :: ncf_removeDimOrderEntry
    ! deletes dimension entry from dimIndex
    type(dimensionOrder),pointer :: i
    integer :: indDim
    type(dimension), pointer :: d,n
    character*20 :: myname = "ncf_removeDimOrderEntry"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    nullify(ncf_removeDimOrderEntry)
    d=>i%firstDimension%next
    do while (.not. associated(d,target=i%lastDimension))
       n=>d%next
       if (d%ind.eq.indDim) then
          d%prev%next=>d%next
          d%next%prev=>d%prev
          ncf_removeDimOrderEntry=>d
          i%nrdim=i%nrdim-1
          d=>i%lastDimension
       else
          d=>n
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end function ncf_removeDimOrderEntry
  !
  ! get dimension from dimIndex
  !
  function ncf_getDimensionOrderDimension(ido,indDim)
    implicit none
    ! returns the indicated dimension
    type(dimension), pointer :: ncf_getDimensionOrderDimension
    type(dimensionOrder) :: ido
    integer :: indDim
    integer :: kk
    type(dimension), pointer :: d
    character*20 :: myname = "ncf_getDimensionOrderDimension"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    nullify(ncf_getDimensionOrderDimension)
    kk=0
    d=>ido%firstDimension%next
    do while (.not.associated(d,ido%lastDimension))
       kk=kk+1
       if (kk.eq.indDim) then
          ncf_getDimensionOrderDimension=>d
          d=>ido%lastDimension
       else
          d=>d%next
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end function ncf_getDimensionOrderDimension
  !
  integer function ncf_getDimOrderCount(ido)
    implicit none
    ! returns the indicated dimension
    type(dimensionOrder), pointer :: ido
    integer :: kk
    type(dimension), pointer :: d
    character*20 :: myname = "ncf_getDimOrderCount"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    kk=0
    if (associated(ido)) then
       d=>ido%firstDimension%next
       do while (.not.associated(d,ido%lastDimension))
          kk=kk+1
          d=>d%next
       end do
    end if
    ncf_getDimOrderCount=kk
    if (ncf_bdeb) write(*,*)myname,' Done.',kk
  end function ncf_getDimOrderCount
  !
  function ncf_getInventoryDimension(i,indDim,irc)
    implicit none
    ! returns the indicated dimension
    type(dimension), pointer :: ncf_getInventoryDimension
    type(inventory) :: i
    integer :: indDim
    integer :: irc
    character*20 :: myname = "ncf_getInventoryDimension"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    nullify(ncf_getInventoryDimension)
    if (indDim.le.i%nrdim) then
       allocate(ncf_getInventoryDimension)
       ncf_getInventoryDimension%ind=indDim
       ncf_getInventoryDimension%sta=i%pos%sta(indDim)
       ncf_getInventoryDimension%lim=i%pos%lim(indDim)
       ncf_getInventoryDimension%siz=i%pos%siz(indDim)
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end function ncf_getInventoryDimension
  !
  function ncf_getPrimaryDimOrder(ret,irc)
    implicit none
    type(dimensionOrder), pointer :: ncf_getPrimaryDimOrder
    type(inventory),pointer :: ret
    integer :: irc
    type(dimensionOrder), pointer :: xdo,ydo
    type(dimension), pointer :: dx,dy
    integer :: nx,ny
    character*20 :: myname = "ncf_getPrimaryDimOrder"
    nullify(ncf_getPrimaryDimOrder)
    ! try x-y oordinates...
    if (associated(ret%xcoid)) then
       xdo=>ncf_makeDimOrder(ret%xcoid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
          return
       end if
       nx=ncf_getDimOrderCount(xdo)
       if (nx.eq.2) then
          ncf_getPrimaryDimOrder=>xdo
          return
       else if (nx.eq.1.and.associated(ret%ycoid)) then
          ydo=>ncf_makeDimOrder(ret%ycoid,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
             return
          end if
          ny=ncf_getDimOrderCount(ydo)
          if (ny.eq.1) then
             call ncf_addDimOrder(xdo,ydo,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_addDimOrder.',irc
                return
             end if
             call ncf_clearDimOrder(ydo)
             ncf_getPrimaryDimOrder=>xdo
             return
          end if
          call ncf_clearDimOrder(ydo)
       end if
       call ncf_clearDimOrder(xdo)
    end if
    ! try lat-lon
    if (associated(ret%lonid)) then
       xdo=>ncf_makeDimOrder(ret%lonid,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
          return
       end if
       nx=ncf_getDimOrderCount(xdo)
       if (nx.eq.2) then
          ncf_getPrimaryDimOrder=>xdo
          return
       else if (nx.eq.1.and.associated(ret%latid)) then
          ydo=>ncf_makeDimOrder(ret%latid,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
             return
          end if
          ny=ncf_getDimOrderCount(ydo)
          if (ny.eq.1) then
             call ncf_addDimOrder(xdo,ydo,irc)
             if (irc.ne.0) then
                write(*,*)myname,'Error return from ncf_addDimOrder.',irc
                return
             end if
             call ncf_clearDimOrder(ydo)
             ncf_getPrimaryDimOrder=>xdo
             return
          end if
          call ncf_clearDimOrder(ydo)
       end if
       call ncf_clearDimOrder(xdo)
    end if
    ! nothing is working...
    irc=345
    WRITE(*,*) MYNAME," Unable to find primary dimensions."
    RETURN
    return
  end function ncf_getPrimaryDimOrder
  !
  function ncf_makeDimension(i,indDim)
    implicit none
    ! returns the indicated dimension
    type(dimension), pointer :: ncf_makeDimension
    type(inventory) :: i
    integer indDim
    character*20 :: myname = "ncf_makeDimension"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(ncf_makeDimension)
    ncf_makeDimension%ind=indDim
    ncf_makeDimension%sta=i%pos%sta(indDim)
    ncf_makeDimension%lim=i%pos%lim(indDim)
    ncf_makeDimension%siz=i%pos%siz(indDim)
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end function ncf_makeDimension
  !
  function ncf_dimOrderContainsDim(o,indDim)
    implicit none
    ! returns the indicated dimension
    logical :: ncf_dimOrderContainsDim
    type(dimensionOrder) :: o
    type(Dimension), pointer :: d
    integer indDim
    integer kk
    character*20 :: myname = "ncf_dimOrderContainsDim"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    d=>o%firstDimension%next
    do while (.not. associated(d, target=o%lastDimension))
       if (d%ind.eq.indDim) then
          ncf_dimOrderContainsDim=.true.
          if (ncf_bdeb) write(*,*)myname,' Done.',ncf_dimOrderContainsDim
          return
       end if
       d=>d%next
    end do
    ncf_dimOrderContainsDim=.false.
    if (ncf_bdeb) write(*,*)myname,' Done.',ncf_dimOrderContainsDim
    return
  end function ncf_dimOrderContainsDim
  !
  ! check if dimIndex contains a specific index
  !
  function ncf_dimIndexContains(i,indDim,irc)
    implicit none
    ! checks in a dimIndex contains ind-th dimension
    logical :: ncf_dimIndexContains
    type(dimensionOrder),pointer :: i
    integer :: indDim
    integer irc
    integer kk
    type(dimension), pointer :: d
    character*20 :: myname = "ncf_dimIndexContains"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(i)) then
       d=>i%firstDimension%next
       do while (.not. associated(d,target=i%lastDimension))
          if (d%ind.eq.indDim) then
             ncf_dimIndexContains=.true.
             return
          end if
          d=>d%next
       end do
    end if
    ncf_dimIndexContains=.false.
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end function ncf_dimIndexContains
  !
  !
  function ncf_variableContainsDim(v,indDim)
    implicit none
    ! returns the indicated dimension
    logical :: ncf_variableContainsDim
    type(variable) :: v
    integer indDim
    integer kk
    character*20 :: myname = "ncf_variableContainsDim"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    do kk=1,v%nrdim
       if (v%ind(kk).eq.indDim) then
          ncf_variableContainsDim=.true.
          if (ncf_bdeb) write(*,*)myname,' Done.',ncf_variableContainsDim
          return
       end if
    end do
    ncf_variableContainsDim=.false.
    if (ncf_bdeb) write(*,*)myname,' Done.',ncf_variableContainsDim
    return
  end function ncf_variableContainsDim
  !
  function ncf_makeUnionDimOrder(t,i,irc)
    implicit none
    ! returns an dimIndex of the union between dimIndex "i" and dimIndex "t"
    type(dimensionOrder),pointer :: ncf_makeUnionDimOrder
    type(dimensionOrder),pointer  :: t
    type(dimensionOrder),pointer  :: i
    type(dimension),pointer  :: td, id, x
    logical :: found
    integer irc
    character*20 :: myname = "ncf_makeUnionDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(ncf_makeUnionDimOrder,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    allocate(ncf_makeUnionDimOrder%firstDimension,ncf_makeUnionDimOrder%lastDimension, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_makeUnionDimOrder%firstDimension%next=>ncf_makeUnionDimOrder%lastDimension
    ncf_makeUnionDimOrder%lastDimension%prev=>ncf_makeUnionDimOrder%firstDimension
    ncf_makeUnionDimOrder%nrdim =0
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
             if (irc.ne.0) then
                write(*,*)myname,' Unable to allocate.'
                return
             end if
             x%ind=id%ind
             x%lim=id%lim
             x%siz=id%siz
             x%sta=id%sta
             x%prev=>ncf_makeUnionDimOrder%lastDimension%prev
             x%next=>ncf_makeUnionDimOrder%lastDimension
             ncf_makeUnionDimOrder%lastDimension%prev%next=>x
             ncf_makeUnionDimOrder%lastDimension%prev=>x
             ncf_makeUnionDimOrder%nrdim = ncf_makeUnionDimOrder%nrdim + 1
             nullify(x)
          end if
          id=>id%next
       end do
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_makeUnionDimOrder

  function ncf_copyDimOrder(i,irc)
    implicit none
    ! copies index "i" to "t"
    type(dimensionOrder),pointer :: ncf_copyDimOrder
    type(dimensionOrder) :: i
    integer :: irc
    type(dimension),pointer  :: d, cd
    character*20 :: myname = "ncf_copyDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(ncf_copyDimOrder,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.',irc
       return
    end if
    allocate(ncf_copyDimOrder%firstDimension,ncf_copyDimOrder%lastDimension, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate dimension-chain.',irc
       return
    end if
    ncf_copyDimOrder%firstDimension%next=>ncf_copyDimOrder%lastDimension
    ncf_copyDimOrder%lastDimension%prev=>ncf_copyDimOrder%firstDimension
    ncf_copyDimOrder%nrdim=0
    d=>i%firstDimension%next
    do while (.not.associated(d,i%lastDimension))
       allocate(cd,stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate dimension.',irc
          return
       end if
       cd%ind=d%ind
       cd%lim=d%lim
       cd%siz=d%siz
       cd%sta=d%sta
       cd%prev=>ncf_copyDimOrder%lastDimension%prev
       cd%next=>ncf_copyDimOrder%lastDimension
       ncf_copyDimOrder%lastDimension%prev%next=>cd
       ncf_copyDimOrder%lastDimension%prev=>cd
       ncf_copyDimOrder%nrdim = ncf_copyDimOrder%nrdim + 1
       nullify(cd)
       d=>d%next
    end do
    ncf_copyDimOrder%f=>i%f
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_copyDimOrder

  integer function ncf_getDimOrderLength(i)
    implicit none
    type(dimensionOrder), pointer :: i
    integer :: lenc
    type(dimension), pointer :: d
    character*20 :: myname = "ncf_getDimOrderLength"
    ncf_getDimOrderLength=1
    d=>i%firstDimension%next
    do while (.not.associated(d,target=i%lastDimension))
       ncf_getDimOrderLength = ncf_getDimOrderLength * d%lim
       d=>d%next
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end function ncf_getDimOrderLength
  !
  !@#  POSITION ******************************************************

  subroutine ncf_allocatePos(nrdim,pos,irc)
    implicit none
    integer nrdim
    type(position) :: pos
    integer irc
    character*20 :: myname = "ncf_allocatePos"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    pos%maxdim=nrdim
    pos%nrdim=nrdim
    allocate(pos%pos(max(1,pos%nrdim)),pos%sta(max(1,pos%nrdim)),&
         & pos%lim(max(1,pos%nrdim)),pos%siz(max(1,pos%nrdim)),&
         & stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate pos/sta/lim.',irc
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_allocatePos

  subroutine ncf_deallocatePos(pos,irc)
    implicit none
    type(position) :: pos
    integer irc
    character*20 :: myname = "ncf_deallocatePos"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(pos%sta)) deallocate(pos%sta,stat=irc)
    if (associated(pos%lim)) deallocate(pos%lim,stat=irc)
    if (associated(pos%siz)) deallocate(pos%siz,stat=irc)
    if (associated(pos%pos)) deallocate(pos%pos,stat=irc)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_deallocatePos

  subroutine ncf_setPos(f,i,value)
    implicit none
    type(inventory) :: f
    type(dimensionOrder)  :: i
    integer :: value
    type(dimension), pointer :: d
    character*20 :: myname = "ncf_setPos"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    d=>i%firstDimension%next
    do while (.not.associated(d,target=i%lastDimension))
       if (d%ind.lt.1.or.d%ind.gt.f%nrdim) then
          write(*,*) 'ncf_setPos Ivalid index:',d%ind,f%nrdim
       else
          f%pos%pos(d%ind)=value
       end if
       d=>d%next
    end do
    f%pos%loc=0
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_setPos

  subroutine ncf_addPos(p,t,irc)
    implicit none
    type(position) :: p
    type(position) :: t
    integer irc
    integer kk
    character*20 :: myname = "ncf_addPos"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    do kk=1,p%nrdim
       p%pos(kk)=p%pos(kk)+t%pos(kk)
    end do
    p%loc=t%loc
    !if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_addPos

  subroutine ncf_copyPos(p,t,irc)
    implicit none
    type(position), pointer :: p
    type(position) :: t
    integer irc
    integer kk
    type(position), pointer :: b
    character*20 :: myname = "ncf_copyPos"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (t%nrdim .gt. p%nrdim) then
       allocate(b,stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       call ncf_allocatePos(t%nrdim,b,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_allocatePos.',irc
          return
       end if
       p=>b
    end if
    p%nrdim=t%nrdim
    do kk=1,p%nrdim
       p%pos(kk)=t%pos(kk)
       p%sta(kk)=t%sta(kk)
       p%lim(kk)=t%lim(kk)
       p%siz(kk)=t%siz(kk)
    end do
    p%loc=t%loc
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_copyPos

  integer function ncf_getRemoteLocation(f,v)
    implicit none
    type(inventory), pointer :: f
    type(variable) :: v
    integer kk
    character*20 :: myname = "ncf_getRemoteLocation"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (.not.associated(f)) then
       ncf_getRemoteLocation=0
    else
       if (v%nrdim .gt. 0) then
          ncf_getRemoteLocation=f%pos%pos(v%ind(v%nrdim))-v%sta(v%nrdim)
          ! write(*,*) 'ncf_getLocation A:',ncf_getLocation
          ! write(*,*) 'ncf_getLocation B:',v%nrdim
          ! write(*,*) 'ncf_getLocation C:',v%ind(v%nrdim)
          ! call ncf_printPos(f%pos)
          do kk=v%nrdim-1,1,-1
             ncf_getRemoteLocation=ncf_getRemoteLocation*v%lim(kk)+f%pos%pos(v%ind(kk))-v%sta(kk)
          end do
          ncf_getRemoteLocation=ncf_getRemoteLocation+1
       else
          ncf_getRemoteLocation=1
       end if
       f%pos%loc=ncf_getRemoteLocation
    end if
    !if (ncf_bdeb) write(*,*)myname,' Done.',ncf_getRemoteLocation
    return
  end function ncf_getRemoteLocation

  integer function ncf_getLocation(v)
    implicit none
    type(variable) :: v
    integer kk
    character*20 :: myname = "ncf_getLocation"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (.not.associated(v%f)) then
       write(*,*)myname,'Zombie variable:',v%var250(1:v%lenv)
       ncf_getLocation=0
    else
       if (v%nrdim .gt. 0) then
          ncf_getLocation=v%f%pos%pos(v%ind(v%nrdim))-v%sta(v%nrdim)
          ! write(*,*) 'ncf_getLocation A:',ncf_getLocation
          ! write(*,*) 'ncf_getLocation B:',v%nrdim
          ! write(*,*) 'ncf_getLocation C:',v%ind(v%nrdim)
          ! call ncf_printPos(v%f%pos)
          do kk=v%nrdim-1,1,-1
             ncf_getLocation=ncf_getLocation*v%lim(kk)+v%f%pos%pos(v%ind(kk))-v%sta(kk)
          end do
          ncf_getLocation=ncf_getLocation+1
       else
          ncf_getLocation=1
       end if
       if (ncf_getlocation.eq.0) then
          write(*,*)myname,'Strange variable:',v%var250(1:v%lenv)
          do kk=v%nrdim-1,1,-1
             write(*,*)myname,'*** ',kk,v%ind(kk),v%sta(kk),v%lim(kk),v%f%pos%pos(v%ind(kk))
          end do
       end if
       v%f%pos%loc=ncf_getLocation
    end if
    !if (ncf_bdeb) write(*,*)myname,' Done.',ncf_getLocation
    return
  end function ncf_getLocation
  !
  subroutine ncf_resetPos(i,irc)
    implicit none
    type(dimensionOrder) :: i
    integer irc
    type(dimension),pointer :: d
    character*20 :: myname = "ncf_resetPos"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    i%resetPos=.true.
    d => i%firstDimension%next
    do while (.not. associated(d,target=i%lastDimension))
       d%sta=i%f%pos%sta(d%ind)
       d%lim=i%f%pos%lim(d%ind)
       d=>d%next
    end do
    !if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_resetPos
  !
  subroutine ncf_lastDimension(f,d)
    implicit none
    type(inventory),pointer :: f
    type(dimension),pointer :: d
    character*20 :: myname = "ncf_lastDimension"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    f%pos%pos(d%ind)=d%lim
    f%pos%loc=0
    !if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_lastDimension
  !
  subroutine ncf_firstDimension(f,d)
    implicit none
    type(inventory),pointer :: f
    type(dimension),pointer :: d
    character*20 :: myname = "ncf_firstDimension"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    f%pos%pos(d%ind)=1
    f%pos%loc=0
    !if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_firstDimension
  !
  logical function ncf_incrementDimension(f,d,biok,inc,buff,irc)
    implicit none
    type(inventory) :: f
    type(dimension) :: d
    logical biok
    integer inc
    integer buff
    integer irc
    logical dok
    character*20 :: myname = "ncf_incrementDimension"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
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
    ncf_incrementDimension=biok
    !if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_incrementDimension
  !
  logical function ncf_decrementDimension(f,d,biok,inc,buff,irc)
    implicit none
    type(inventory) :: f
    type(dimension) :: d
    logical biok
    integer inc
    integer buff
    integer irc
    character*20 :: myname = "ncf_decrementDimension"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
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
    ncf_decrementDimension=biok
    !if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_decrementDimension
  !
  logical function ncf_increment(f,i,irc)
    implicit none
    type(inventory) :: f
    type(dimensionOrder) :: i
    integer irc
    type(dimension),pointer :: d
    character*20 :: myname = "ncf_increment"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    ncf_increment=.false.
    if (i%resetPos) then
       d=>i%firstDimension%next
       do while (.not.associated(d,target=i%lastDimension))
          f%pos%pos(d%ind)=d%sta
          d%increase=.true.
          d=>d%next
       end do
       ncf_increment=.true. ! there is always one valid return
       i%resetPos=.false.
    else
       d=>i%firstDimension%next
       do while (.not.associated(d,target=i%lastDimension))
          if (d%ind.le.f%nrdim.and.d%ind.ge.1) then
             if (d%increase) then
                f%pos%pos(d%ind)=f%pos%pos(d%ind)+1
                if (f%pos%pos(d%ind).gt.d%sta+d%lim-1) then
                   d%increase=.false.
                   f%pos%pos(d%ind)=d%sta+d%lim-1
                   d=>d%next
                else
                   ncf_increment=.true.
                   d=>i%lastDimension
                end if
             else
                f%pos%pos(d%ind)=f%pos%pos(d%ind)-1
                if (f%pos%pos(d%ind).lt.d%sta) then
                   d%increase=.true.
                   f%pos%pos(d%ind)=d%sta
                   d=>d%next
                else
                   ncf_increment=.true.
                   d=>i%lastDimension
                end if
             end if
          else
             irc=703
             return
          end if
       end do
    end if
    f%pos%loc=0
    !if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_increment

  ! !@#  FOOTPRINT ******************************************************

  ! logical function ncf_getFootIndex(f,i,s,latc,lonc,diam,irc)
  !   implicit none
  !   type(inventory), pointer    :: f
  !   type(dimensionOrder), pointer        :: i
  !   type(dimensionOrder), pointer        :: s
  !   real latc
  !   real lonc
  !   real diam ! diameter in degrees
  !   integer irc
  !   integer dn(2),ii
  !   logical biok
  !   real latx,lonx,dist
  !   type(dimension), pointer :: d,e
  !   character*20 :: myname = "ncf_getFootIndex"
  !   if (ncf_bdeb) write(*,*)myname,' Entering.'
  !   ! find increment distance in each dimension direction
  !   biok=.true.
  !   latc=ncf_valuePosition(f%latid,irc)
  !   lonc=ncf_valuePosition(f%lonid,irc)
  !   do ii=1,i%nrdim ! must be 2 (checked earlier)
  !      d=>ncf_getDimensionOrderDimension(i,ii)
  !      e=>ncf_getDimensionOrderDimension(s,ii)
  !      if (d%ind .ne. e%ind) then
  !         write(*,*)'System error:',d%ind,e%ind
  !         irc=845
  !         return
  !      end if
  !      e%sta=f%pos%pos(e%ind)
  !      if (biok) biok=ncf_incrementDimension(f,d,biok,1,0,irc)
  !      if (biok) then
  !         latx=ncf_valuePosition(f%latid,irc)
  !         lonx=ncf_valuePosition(f%lonid,irc)
  !         biok=ncf_decrementDimension(f,d,biok,1,0,irc)
  !         dist=2.0D0*ncf_getDist(latc,lonc,latx,lonx)
  !         dn(ii)=max(1,1+int(diam/max(dist,1.0D-10)))
  !      end if
  !   end do
  !   ! determine lower left corner
  !   if (biok) then
  !      do ii=1,i%nrdim ! must be 2 (checked earlier)
  !         d=>ncf_getDimensionOrderDimension(i,ii)
  !         biok=ncf_decrementDimension(f,d,biok,dn(ii),0,irc)
  !         e=>ncf_getDimensionOrderDimension(s,ii)
  !         e%sta=e%sta-dn(ii)
  !         e%lim=e%sta+2*dn(ii)
  !         if (e%sta.lt.1.or.e%lim.gt.f%pos%lim(e%ind)) biok=.false.
  !      end do
  !   end if
  !   !    call ncf_printDimOrder(s)
  !   ncf_getFootIndex=biok
  !   if (ncf_bdeb) write(*,*)myname,' Done.',irc
  ! end function ncf_getFootIndex

  ! logical function ncf_incrementFoot(f,i,latc,lonc,diam,irc)
  !   implicit none
  !   type(inventory), pointer    :: f
  !   type(dimensionOrder), pointer        :: i
  !   real latc
  !   real lonc
  !   real diam
  !   logical reset
  !   integer irc
  !   real latx,lonx,dist
  !   logical bdone
  !   character*20 :: myname = "ncf_incrementFoot"
  !   if (ncf_bdeb) write(*,*)myname,' Entering.'
  !   bdone=.false.
  !   do while (.not. bdone)
  !      if (ncf_increment(f,i,irc)) then
  !         ! call ncf_printDimOrder(i)
  !         ! call ncf_printPos(f%pos)
  !         latx=ncf_valuePosition(f%latid,irc)
  !         lonx=ncf_valuePosition(f%lonid,irc)
  !         dist=2.0D0*ncf_getDist(latc,lonc,latx,lonx)
  !         if (dist.lt.diam) then ! grid wihin footprint...
  !            bdone=.true.
  !         end if
  !         ncf_incrementFoot=.true.
  !      else !     no more grid points to check
  !         ncf_incrementFoot=.false.
  !         bdone=.true.
  !      end if
  !   end do
  !   if (ncf_bdeb) write(*,*)myname,' Done.',irc
  ! end function ncf_incrementFoot

  !@#  VARIABLES ******************************************************

  subroutine ncf_replaceVariable(i,o) ! replace variable "i" by "o" in inventory "f"...
    type(variable),pointer :: i,o
    character*20 :: myname = "ncf_replaceVariable"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(o%next).and. associated(o%prev))then ! remove from old chain
       o%next%prev=>o%prev
       o%prev%next=>o%next
    end if
    o%next=>i%next
    o%prev=>i%prev
    o%next%prev=>o
    o%prev%next=>o
    nullify(i%prev)
    nullify(i%next)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_replaceVariable

  subroutine ncf_prependVariable(inv,o) ! replace variable "i" by "o" in inventory "f"...
    type(inventory),pointer :: inv
    type(variable),pointer :: o
    character*20 :: myname = "ncf_replaceVariable"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(o,target=inv%firstVariable%next)) return
    if (associated(o%next)) o%next%prev => o%prev
    if (associated(o%prev)) o%prev%next => o%next
    o%next=>inv%firstVariable%next
    o%prev=>inv%firstVariable
    o%next%prev=>o
    o%prev%next=>o
    o%varid=1+ncf_maxVarid(inv,irc)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_prependVariable

  subroutine ncf_appendVariable(inv,o) ! replace variable "i" by "o" in inventory "f"...
    type(inventory),pointer :: inv
    type(variable),pointer :: o
    character*20 :: myname = "ncf_replaceVariable"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(o,target=inv%lastVariable%prev)) return
    if (associated(o%next)) o%next%prev => o%prev
    if (associated(o%prev)) o%prev%next => o%next
    o%next=>inv%lastVariable
    o%prev=>inv%lastVariable%prev
    o%next%prev=>o
    o%prev%next=>o
    o%varid=1+ncf_maxVarid(inv,irc)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_appendVariable

  subroutine ncf_cleanVariable(v,irc)
    implicit none
    type (Variable),pointer :: v
    integer irc
    type(attribute), pointer :: a,ax
    character*20 :: myname = "ncf_cleanVariable"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (.not.associated(v)) return
    if (v%initialised) then
       a=>v%firstAttribute%next
       do while (.not.associated(a,target=v%lastAttribute))
          ax=>a%next
          call ncf_removeAttribute(a,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_removeAttribute.',irc
             return
          end if
          a=>ax
       end do
       call ncf_removeAttribute(v%firstAttribute,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_removeAttribute.',irc
          return
       end if
       call ncf_removeAttribute(v%lastAttribute,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_removeAttribute.',irc
          return
       end if
       nullify(v%fillattribute)
       if (associated(v%ind)) deallocate(v%ind,stat=irc)
       if (associated(v%lim)) deallocate(v%lim,stat=irc)
       if (associated(v%siz)) deallocate(v%siz,stat=irc)
       if (associated(v%sta)) deallocate(v%sta,stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,'Unable to deallocate...',irc
          return
       end if
       call ncf_clearVariable(v,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_clearVariable.',irc
          return
       end if
       v%initialised=.false.
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_cleanVariable

  subroutine ncf_removeVariable(v,irc)
    implicit none
    type (Variable),pointer :: v
    integer irc
    type(attribute), pointer :: a,ax
    character*20 :: myname = "ncf_removeVariable"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (.not.associated(v)) return
    if (associated(v%prev))v%prev%next=>v%next
    if (associated(v%next))v%next%prev=>v%prev
    call ncf_cleanVariable(v,irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from ncf_cleanVariable.',irc
       return
    end if
    deallocate(v,stat=irc)
    nullify(v)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_removeVariable

  function ncf_copyVariable(v,irc)
    implicit none
    type(variable), pointer :: ncf_copyVariable
    type(variable) v
    integer irc
    integer ii
    type(attribute),pointer :: a,va
    character*20 :: myname = "ncf_copyVariable"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(ncf_copyVariable,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_copyVariable%initialised=v%initialised
    ncf_copyVariable%ignored=v%ignored
    ncf_copyVariable%f=>v%f
    ncf_copyVariable%ncid=v%ncid
    ncf_copyVariable%varid=v%varid
    ncf_copyVariable%var250=v%var250
    ncf_copyVariable%lenv=v%lenv
    ncf_copyVariable%type=v%type
    ncf_copyVariable%nrdim=v%nrdim
    ncf_copyVariable%nratt=v%nratt
    ncf_copyVariable%scale=v%scale
    ncf_copyVariable%fillc=v%fillc
    ncf_copyVariable%fill1=v%fill1
    ncf_copyVariable%fill2=v%fill2
    ncf_copyVariable%fill4=v%fill4
    ncf_copyVariable%fillr=v%fillr
    ncf_copyVariable%filld=v%filld
    allocate(ncf_copyVariable%ind(max(1,ncf_copyVariable%nrdim)),&
         & ncf_copyVariable%sta(max(1,ncf_copyVariable%nrdim)),&
         & ncf_copyVariable%lim(max(1,ncf_copyVariable%nrdim)),&
         & ncf_copyVariable%siz(max(1,ncf_copyVariable%nrdim)),&
         & stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    if (ncf_bdeb)write(*,*)myname,"Var:'"//v%var250(1:v%lenv)//"'",v%nrdim,associated(v%ind),v%len
    do ii=1,v%nrdim
       ncf_copyVariable%ind(ii)=v%ind(ii)
       ncf_copyVariable%sta(ii)=v%sta(ii)
       ncf_copyVariable%lim(ii)=v%lim(ii)
       ncf_copyVariable%siz(ii)=v%siz(ii)
    end do
    ncf_copyVariable%len=v%len
    ncf_copyVariable%lens=v%lens
    ncf_copyVariable%lenc=0
    ncf_copyVariable%len1=0
    ncf_copyVariable%len2=0
    ncf_copyVariable%len4=0
    ncf_copyVariable%len8=0
    ncf_copyVariable%lenr=0
    ncf_copyVariable%lend=0
    allocate(ncf_copyVariable%firstAttribute,ncf_copyVariable%lastAttribute, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_copyVariable%firstAttribute%next=>ncf_copyVariable%lastAttribute
    ncf_copyVariable%lastAttribute%prev=>ncf_copyVariable%firstAttribute
    if (associated(v%firstAttribute)) then
       va=>v%firstAttribute%next
       do while (.not.associated(va,v%lastAttribute))
          a=>ncf_copyAttribute(va,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyAttribute.',irc
             return
          end if
          a%prev=>ncf_copyVariable%lastAttribute%prev
          a%next=>ncf_copyVariable%lastAttribute
          ncf_copyVariable%lastAttribute%prev%next=>a
          ncf_copyVariable%lastAttribute%prev=>a
          a%v=>ncf_copyVariable
          nullify(a)
          va=>va%next
       end do
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_copyVariable

  integer function ncf_getVariableLength(v)
    implicit none
    type(variable), pointer :: v
    integer :: ii
    character*20 :: myname = "ncf_getVariableLength"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ncf_getVariableLength=1
    do ii=1,v%nrdim
       ncf_getVariableLength = ncf_getVariableLength * v%lim(ii)
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end function ncf_getVariableLength
  !@#  VALUES/WEIGHT ******************************************************

  subroutine ncf_setValue(v,val,irc)
    implicit none
    type(variable), pointer :: v
    real val
    integer irc
    integer loc
    type(dimensionOrder),pointer :: i
    character*20 :: myname = "ncf_setValue"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    loc=ncf_getLocation(v)
    if (loc.lt.1.or.loc.gt.v%lend) then
       write(*,*)myname,' Invalid position.',loc,v%lend
       i=>ncf_makeDimOrder(v,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
          return
       end if
       call ncf_printDimOrder(i)
       call ncf_printpos(v%f%pos)
    end if
    if(ncf_bdeb)write(*,*)'ncf_addValue location:',loc,val
    v%fd(loc)=val
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_setValue

  subroutine ncf_addValue(v,val,irc)
    implicit none
    type(variable), pointer :: v
    real val
    integer irc
    integer loc
    type(dimensionOrder),pointer :: i
    character*20 :: myname = "ncf_addValue"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    loc=ncf_getLocation(v)
    ! write(*,*)'ncf_addValue location:',loc
    if (loc.lt.1.or.loc.gt.v%lend) then
       write(*,*)myname,' Invalid position.',loc,v%lend
       i=>ncf_makeDimOrder(v,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_makeDimOrder.',irc
          return
       end if
       call ncf_printDimOrder(i)
       call ncf_printpos(v%f%pos)
    end if
    v%fd(loc)=v%fd(loc)+val
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_addValue

  real function ncf_valuePosition(v,irc)
    implicit none
    type(variable) :: v
    type(dimensionOrder) :: i
    integer irc
    integer loc
    character*20 :: myname = "ncf_valuePosition"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    loc=ncf_getLocation(v)
    if (.not.associated(v%fd)) then
       irc=999
       write(*,*)myname," Variable not allocated: '"//v%var250(1:v%lenv)//"'",v%lend
       return
    end if
    ncf_valuePosition=v%fd(loc)
   ! if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_valuePosition

  real function ncf_valueWeighted(v,w,irc)
    implicit none
    type(variable) :: v
    type(weight) :: w
    integer irc
    integer kk
    integer loc
    logical fill
    real val
    character*20 :: myname = "ncf_valueWeighted"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    fill=.false.
    ncf_valueWeighted=0.0D0
    do kk=1,w%nweight
       call ncf_addPos(v%f%pos,w%pos(kk),irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_addPos.',irc
          return
       end if
       loc=ncf_getLocation(v)
       !call ncf_printPos(v%f%pos)
       !write(*,*)'ValueWeighted LOC:',loc, v%len, v%lend, v%lenr, v%type, nf_real, nf_double
       val=v%fd(loc)
       if (val.eq.v%filld) then
          fill=.true.
       else
          if (abs(val).gt.1000.0D0) then
             write(*,*) myname,' Invalid ens:',val
          end if
          ncf_valueWeighted=ncf_valueWeighted + w%w(kk) * val
       end if
    end do
    if (fill) then
       ncf_valueWeighted=v%filld
    end if
    !if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_valueWeighted

  function ncf_makeWeight4(nrdim,irc)
    implicit none
    type(weight), pointer :: ncf_makeWeight4
    integer nrdim
    integer irc
    integer kk
    character*20 :: myname = "ncf_makeWeight4"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(ncf_makeWeight4,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_makeWeight4%nweight=4
    allocate(ncf_makeWeight4%pos(ncf_makeWeight4%nweight),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    do kk=1,ncf_makeWeight4%nweight
       call ncf_allocatePos(nrdim,ncf_makeWeight4%pos(kk),irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_allocatePos.',irc
          return
       end if
    end do
    allocate(ncf_makeWeight4%w(ncf_makeWeight4%nweight),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_makeWeight4

  subroutine ncf_clearWeight(wgt)
    type(weight),pointer :: wgt
    if (associated(wgt)) then
       if (associated(wgt%pos)) deallocate(wgt%pos)
       if (associated(wgt%w)) deallocate(wgt%w)
       nullify(wgt)
    end if
    return
  end subroutine ncf_clearWeight

  !@#  ATTRIBUTES ******************************************************

  subroutine ncf_removeAttribute(a,irc)
    implicit none
    type (attribute),pointer :: a
    integer irc
    character*20 :: myname = "ncf_clearAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(a%ac)) deallocate(a%ac,stat=irc)
    if (associated(a%a1)) deallocate(a%a1,stat=irc)
    if (associated(a%a2)) deallocate(a%a2,stat=irc)
    if (associated(a%a4)) deallocate(a%a4,stat=irc)
    if (associated(a%a8)) deallocate(a%a8,stat=irc)
    if (associated(a%ar)) deallocate(a%ar,stat=irc)
    if (associated(a%ad)) deallocate(a%ad,stat=irc)
    if (associated(a%prev))a%prev%next=>a%next
    if (associated(a%next))a%next%prev=>a%prev
    a%initialised=.false.
    deallocate(a,stat=irc)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_removeAttribute
  !
  subroutine ncf_clearAttributes(v,irc)
    implicit none
    type(variable), pointer :: v
    integer :: irc
    type(attribute), pointer :: a, ax
    character*20 :: myname = "ncf_clearAttributes"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(v%firstAttribute)) then
       a=>v%firstAttribute%next
       do while(.not.associated(a,target=v%lastAttribute))
          ax=>a%next
          call ncf_removeAttribute(a,irc)
          if (irc.ne.0) then
             write(*,*) myname,"Error return from ncf_removeAttributes.",irc
             return
          end if
          a=>ax
       end do
    else
       allocate(v%firstAttribute%next,v%lastAttribute%next, stat=irc)
       if(irc.ne.0) then
          write(*,*) myname,' Unable to allocate ATTRIBUTE.',irc
          return
       end if
       v%firstAttribute%next=>v%lastAttribute
       v%lastAttribute%prev=>v%firstAttribute
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_clearAttributes

  subroutine ncf_readAttribute(ncid,varid,attid,a,irc)
    implicit none
    integer :: ncid
    integer :: varid
    integer :: attid
    type(attribute) :: a
    integer :: irc
    integer :: lena
    integer, external :: length
    integer :: ret
    character*20 :: myname = "ncf_readAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    a%ncid=ncid
    a%varid=varid
    a%attid=attid
    if (associated(a%ac)) deallocate(a%ac,stat=irc)
    if (associated(a%a1)) deallocate(a%a1,stat=irc)
    if (associated(a%a2)) deallocate(a%a2,stat=irc)
    if (associated(a%a4)) deallocate(a%a4,stat=irc)
    if (associated(a%a8)) deallocate(a%a8,stat=irc)
    if (associated(a%ar)) deallocate(a%ar,stat=irc)
    if (associated(a%ad)) deallocate(a%ad,stat=irc)
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
       allocate(a%ac(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       ret = nf_get_att_text (ncid,varid,a%att250,a%ac)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_TEXT.',ret
          irc=836
          return
       end if
    else if (a%type   .eq.nf_int1)then
       allocate(a%a1(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       ret = nf_get_att_int1 (ncid,varid,a%att250,a%a1)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_INT1.',ret
          irc=837
          return
       end if
    elseif (a%type   .eq.nf_int2)then
       allocate(a%a2(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       ret = nf_get_att_int2 (ncid,varid,a%att250,a%a2)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_INT2.',ret
          irc=838
          return
       end if
    elseif (a%type   .eq.nf_int)then
       allocate(a%a4(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       ret = nf_get_att_int (ncid,varid,a%att250,a%a4)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_INT.',ret
          irc=839
          return
       end if
    elseif (a%type   .eq.nf_int64)then
       allocate(a%a8(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       notimp=notimp+1
       if (notimp.eq.1) write(*,*)myname,'Int64 not implemented...',a%att250(1:a%lena)
       !       ret = nf_get_att_int64 (ncid,varid,a%att250,a%a8)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_INT.',ret
          irc=839
          return
       end if
    elseif (a%type   .eq.nf_real)then
       allocate(a%ar(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       ret = nf_get_att_real (ncid,varid,a%att250,a%ar)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'READATTRIBUTE Error return from NF_GET_ATT_REAL.',ret
          irc=840
          return
       end if
    elseif (a%type   .eq.nf_double)then
       allocate(a%ad(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
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
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_readAttribute

  function ncf_copyAttribute(a,irc)
    implicit none
    type(attribute), pointer :: ncf_copyAttribute
    type(attribute) a
    integer irc
    integer ii
    character*20 :: myname = "ncf_copyAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering:',a%att250(1:a%lena)
    allocate(ncf_copyAttribute,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_copyAttribute%ncid=a%ncid
    ncf_copyAttribute%varid=a%varid
    ncf_copyAttribute%attid=a%attid
    ncf_copyAttribute%att250=a%att250
    ncf_copyAttribute%lena=a%lena
    ncf_copyAttribute%type=a%type
    ncf_copyAttribute%len=a%len
    if (a%len.gt.0) then
       if (ncf_copyAttribute%type.eq.nf_char) then
          allocate(ncf_copyAttribute%ac(max(1,a%len)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          do ii=1,a%len
             ncf_copyAttribute%ac(ii)=a%ac(ii)
          end do
       else if (ncf_copyAttribute%type.eq.nf_int1) then
          allocate(ncf_copyAttribute%a1(max(1,a%len)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          do ii=1,a%len
             ncf_copyAttribute%a1(ii)=a%a1(ii)
          end do
       else if (ncf_copyAttribute%type.eq.nf_int2) then
          allocate(ncf_copyAttribute%a2(max(1,a%len)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          do ii=1,a%len
             ncf_copyAttribute%a2(ii)=a%a2(ii)
          end do
       else if (ncf_copyAttribute%type.eq.nf_int) then
          allocate(ncf_copyAttribute%a4(max(1,a%len)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          do ii=1,a%len
             ncf_copyAttribute%a4(ii)=a%a4(ii)
          end do
       else if (ncf_copyAttribute%type.eq.nf_int64) then
          allocate(ncf_copyAttribute%a8(max(1,a%len)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          do ii=1,a%len
             ncf_copyAttribute%a8(ii)=a%a8(ii)
          end do
       else if (ncf_copyAttribute%type.eq.nf_real) then
          allocate(ncf_copyAttribute%ar(max(1,a%len)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          do ii=1,a%len
             ncf_copyAttribute%ar(ii)=a%ar(ii)
          end do
       else if (ncf_copyAttribute%type.eq.nf_double) then
          allocate(ncf_copyAttribute%ad(max(1,a%len)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          do ii=1,a%len
             ncf_copyAttribute%ad(ii)=a%ad(ii)
          end do
       end if
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_copyAttribute

  subroutine ncf_readVariableAttributes(ncid,varid,v,irc)
    implicit none
    integer ncid
    integer varid
    type(variable), target :: v
    integer irc
    type(attribute), pointer :: a
    integer jj,ret
    integer, external :: length
    character*20 :: myname = "ncf_readVariableAttributes"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    v%ncid=ncid
    v%varid=varid
    ret = NF_INQ_VARNDIMS (v%ncid, v%varid, v%nrdim);
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_INQ_VARNDIMS:",nf_strerror(ret)
       irc=802
       return
    end if
    allocate(v%ind(max(1,v%nrdim)), stat=irc)
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
    allocate(v%firstAttribute,v%lastAttribute, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    v%firstAttribute%next=>v%lastAttribute
    v%lastAttribute%prev=>v%firstAttribute
    ! read attributes
    do jj=1,v%nratt
       allocate(a,stat=irc)
       if (irc.ne.0) then
          write(*,*) 'READVARIABLEATTRIBUTES Unable to allocate ATTRIBUTE.',irc
          return
       end if
       call ncf_readAttribute(v%ncid,v%varid,jj,a,irc)
       if (irc.ne.0) then
          write(*,*) 'READVARIABLEATTRIBUTES Error return from READATTRIBUTE.',irc
          return
       end if
       a%v=>v
       a%prev=>v%lastAttribute%prev
       a%next=>v%lastAttribute
       v%lastAttribute%prev%next=>a
       v%lastAttribute%prev=>a
       nullify(a)
    end do
    v%initialised=.true.
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_readVariableAttributes

  subroutine ncf_copyVariableAttribute(v,s,att,irc)
    ! copies attributes from "s" to "v"
    implicit none
    type(variable), pointer :: v ! variable
    type(variable), pointer :: s ! source variable
    character*(*) att
    integer :: irc
    type(attribute), pointer :: a,n
    character*250 :: buff250
    integer,external :: length
    integer :: jj, lenb
    character*20 :: myname = "ncf_copyVariableAttributes"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(v).and.associated(s)) then
       a=>ncf_getAttribute(v,att,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_getAttribute.',irc
          return
       end if
       if (associated(a)) then ! unlink and remove
          a%next%prev => a%prev
          a%prev%next => a%next
          if (associated(a%ac)) deallocate(a%ac)
          deallocate(a)
       end if
       a=>ncf_getAttribute(s,att,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_getAttribute.',irc
          return
       end if
       if (associated(a)) then
          n=>ncf_copyAttribute(a,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_copyAttribute.',irc
             return
          end if
          n%v=>v
          n%ncid = v%ncid
          n%varid = v%varid
          n%prev=>v%lastAttribute%prev
          n%next=>v%lastAttribute
          v%lastAttribute%prev%next=>n
          v%lastAttribute%prev=>n
          nullify(n)
       end if
    else
       write(*,*)myname,'Invalid Variable pointers:',associated(v),associated(s)
       irc=915
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_copyVariableAttribute

  subroutine ncf_copyVariableAttributes(v,s,irc)
    ! copies attributes from "s" to "v"
    implicit none
    type(variable), pointer :: v ! variable
    type(variable), target :: s ! source variable
    integer irc
    type(attribute), pointer :: a,newatt
    integer jj
    character*20 :: myname = "ncf_copyVariableAttributes"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    call ncf_clearAttributes(v,irc)
    if (irc.ne.0) then
       write(*,*)myname,'Error return from ncf_clearAttribute.',irc
       return
    end if
    ! read attributes
    v%nratt=s%nratt
    a=>s%firstAttribute%next
    do while (.not.associated(a,target=s%lastAttribute))
       newatt=>ncf_copyAttribute(a,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_copyAttribute.',irc
          return
       end if
       newatt%v=>v
       newatt%ncid = v%ncid
       newatt%varid = v%varid
       newatt%prev=>v%lastAttribute%prev
       newatt%next=>v%lastAttribute
       v%lastAttribute%prev%next=>newatt
       v%lastAttribute%prev=>newatt
       nullify(newatt)
       a=>a%next
    end do
    v%initialised=.true.
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_copyVariableAttributes

  function ncf_getGlobalAttribute(f,att,irc)
    implicit none
    type(attribute), pointer :: ncf_getGlobalAttribute
    type(inventory),pointer :: f
    character*(*) att
    integer lena
    type(attribute),pointer :: a
    integer irc
    character*20 :: myname = "ncf_getGlobalAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    nullify(ncf_getGlobalAttribute)
    if (associated(f)) then
       a=>f%firstAttribute%next
       do while (.not.associated(a,target=f%lastAttribute))
          if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
             ncf_getGlobalAttribute=>a
             a=>f%lastAttribute
          else
             a=>a%next
          end if
       end do
    else
       write(*,*)myname,'Invalid Inventory pointer:',associated(f)
       irc=934
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_getGlobalAttribute

  function ncf_getAttribute(v,att,irc)
    implicit none
    type(attribute), pointer :: ncf_getAttribute
    type(variable),pointer :: v
    character*(*) att
    integer lena
    type(attribute),pointer :: a
    integer irc
    character*20 :: myname = "ncf_getAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    nullify(ncf_getAttribute)
    if (associated(v)) then
       a=>v%firstAttribute%next
       do while (.not.associated(a,target=v%lastAttribute))
          if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
             ncf_getAttribute=>a
             a=>v%lastAttribute
          else
             a=>a%next
          end if
       end do
    else
       write(*,*)myname,'Invalid Variable pointers:',associated(v)
       irc=914
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_getAttribute

  subroutine ncf_setAttribute(v,att,val,type)
    implicit none
    type(variable), pointer :: v
    character*(*) att
    real :: val
    integer :: type
    integer lena
    type(attribute),pointer :: a
    integer ii,irc
    character*20 :: myname = "ncf_setAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    a=>v%firstAttribute%next
    LOOP: do while (.not.associated(a,target=v%lastAttribute))
       if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
          exit LOOP
       else
          a=>a%next
       end if
    end do LOOP
    ! did not find attribute, add new attribute
    if (associated(a,target=v%lastAttribute)) then
       allocate(a,stat=irc)
       if (irc.ne.0) then
          write(*,*) 'ncf_addAttribute Unable to allocate ATTRIBUTE.',irc
          return
       end if
    end if
    ! assign
    a%att250=att
    a%lena=len(att)
    if (associated(a%ac)) deallocate(a%ac,stat=irc)
    if (associated(a%a1)) deallocate(a%a1,stat=irc)
    if (associated(a%a2)) deallocate(a%a2,stat=irc)
    if (associated(a%a4)) deallocate(a%a4,stat=irc)
    if (associated(a%a8)) deallocate(a%a8,stat=irc)
    if (associated(a%ar)) deallocate(a%ar,stat=irc)
    if (associated(a%ad)) deallocate(a%ad,stat=irc)
    a%type=type
    a%len=1
    if (type.eq.nf_int1) then
       allocate(a%a1(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate attribute.',len(att)
          return
       end if
       a%a1(1)=nint(val)
    elseif (type.eq.nf_int2) then
       allocate(a%a2(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate attribute.',len(att)
          return
       end if
       a%a2(1)=nint(val)
    elseif (type.eq.nf_int) then
       allocate(a%a4(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate attribute.',len(att)
          return
       end if
       a%a4(1)=nint(val)
    elseif (type.eq.nf_int64) then
       allocate(a%a8(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate attribute.',len(att)
          return
       end if
       a%a8(1)=nint(val)
    elseif (type.eq.nf_real) then
       allocate(a%ar(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate attribute.',len(att)
          return
       end if
       a%ar(1)=real(val)
    elseif (type.eq.nf_double) then
       allocate(a%ad(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate attribute.',len(att)
          return
       end if
       a%ad(1)=val
    else
       write(*,*)myname,' Invalid type.',type
       irc=999
       return
    end if
    a%v=>v
    a%prev=>v%lastAttribute%prev
    a%next=>v%lastAttribute
    v%lastAttribute%prev%next=>a
    v%lastAttribute%prev=>a
    if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end subroutine ncf_setAttribute

  real function ncf_getRealAttribute(v,att)
    implicit none
    type(variable) v
    character*(*) att
    integer lena
    type(attribute),pointer :: a
    character*20 :: myname = "ncf_getRealAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ncf_getRealAttribute=nf_fill_double
    a=>v%firstAttribute%next
    do while (.not.associated(a,target=v%lastAttribute))
       if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
          if (a%len .gt. 0) then
             if (associated(a%ad)) then
                ncf_getRealAttribute = a%ad(1)
             elseif (associated(a%ar)) then
                ncf_getRealAttribute = a%ar(1)
             elseif (associated(a%a4)) then
                ncf_getRealAttribute = a%a4(1)
             elseif (associated(a%a8)) then
                ncf_getRealAttribute = a%a8(1)
             elseif (associated(a%a2)) then
                ncf_getRealAttribute = a%a2(1)
             elseif (associated(a%a1)) then
                ncf_getRealAttribute = a%a1(1)
             end if
          end if
          a=>v%lastAttribute
       else
          a=>a%next
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end function ncf_getRealAttribute

  subroutine ncf_setDoubleAttribute(v,att,val)
    implicit none
    type(variable), pointer :: v
    character*(*) att
    real :: val
    integer lena
    type(attribute),pointer :: a
    integer ii,irc
    character*20 :: myname = "ncf_setDoubleAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    a=>v%firstAttribute%next
    do while (.not.associated(a,target=v%lastAttribute))
       if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
          if (associated(a%ac)) deallocate(a%ac,stat=irc)
          if (associated(a%a1)) deallocate(a%a1,stat=irc)
          if (associated(a%a2)) deallocate(a%a2,stat=irc)
          if (associated(a%a4)) deallocate(a%a4,stat=irc)
          if (associated(a%a8)) deallocate(a%a8,stat=irc)
          if (associated(a%ar)) deallocate(a%ar,stat=irc)
          if (associated(a%ad)) deallocate(a%ad,stat=irc)
          a%type=nf_double
          a%len=1
          allocate(a%ad(max(1,a%len)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate attribute.',len(att)
          end if
          a%ad(1)=val
          a=>v%lastAttribute
          return ! we are done
       else
          a=>a%next
       end if
    end do
    ! did not find attribute, add new attribute
    allocate(a,stat=irc)
    if (irc.ne.0) then
       write(*,*) 'ncf_addAttribute Unable to allocate ATTRIBUTE.',irc
       return
    end if
    a%att250=att
    a%lena=len(att)
    a%type=nf_double
    a%len=1
    allocate(a%ad(max(1,a%len)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    a%ad(1)=val
    a%v=>v
    a%prev=>v%lastAttribute%prev
    a%next=>v%lastAttribute
    v%lastAttribute%prev%next=>a
    v%lastAttribute%prev=>a
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_setDoubleAttribute


  character*250 function ncf_getTextAttribute(v,att,irc)
    implicit none
    type(variable) :: v
    character*(*) :: att
    integer :: irc
    integer lena,ii
    type(attribute),pointer :: a
    character*20 :: myname = "ncf_getTextAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ncf_getTextAttribute=""
    a=>v%firstAttribute%next
    do while (.not.associated(a,target=v%lastAttribute))
       if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
          if (a%len .gt. 0) then
             if (associated(a%ac)) then
                do ii=1,a%len
                   ncf_getTextAttribute(ii:ii) = a%ac(ii)
                end do
                return
             end if
          end if
          a=>v%lastAttribute
       else
          a=>a%next
       end if
    end do
    irc=458
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end function ncf_getTextAttribute

  subroutine ncf_setTextAttribute(v,att,text)
    implicit none
    type(variable), pointer :: v
    character*(*) att
    character*(*) text
    integer lena
    type(attribute),pointer :: a
    integer ii,irc
    character*20 :: myname = "ncf_setTextAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    a=>v%firstAttribute%next
    do while (.not.associated(a,target=v%lastAttribute))
       if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
          if (associated(a%ac)) deallocate(a%ac,stat=irc)
          if (associated(a%a1)) deallocate(a%a1,stat=irc)
          if (associated(a%a2)) deallocate(a%a2,stat=irc)
          if (associated(a%a4)) deallocate(a%a4,stat=irc)
          if (associated(a%a8)) deallocate(a%a8,stat=irc)
          if (associated(a%ar)) deallocate(a%ar,stat=irc)
          if (associated(a%ad)) deallocate(a%ad,stat=irc)
          a%type=nf_char
          a%len=len(text)
          allocate(a%ac(max(1,a%len)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate attribute.',len(att)
          end if
          do ii=1,a%len
             a%ac(ii)=text(ii:ii)
          end do
          a=>v%lastAttribute
          return ! we are done
       else
          a=>a%next
       end if
    end do
    ! did not find attribute, add new attribute
    allocate(a,stat=irc)
    if (irc.ne.0) then
       write(*,*) 'ncf_addAttribute Unable to allocate ATTRIBUTE.',irc
       return
    end if
    a%att250=att
    a%lena=len(att)
    a%type=nf_char
    a%len=len(text)
    allocate(a%ac(max(1,a%len)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    do ii=1,a%len
       a%ac(ii)=text(ii:ii)
    end do
    a%v=>v
    a%prev=>v%lastAttribute%prev
    a%next=>v%lastAttribute
    v%lastAttribute%prev%next=>a
    v%lastAttribute%prev=>a
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_setTextAttribute


  logical function ncf_hasAttribute(v,att,val,irc)
    implicit none
    type(variable) v
    character*(*) att, val
    integer lena
    type(attribute),pointer :: a
    integer irc
    integer lenb,length
    external length
    character*250 buff250
    character*20 :: myname = "ncf_hasAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ncf_hasAttribute=.false.
    a=>v%firstAttribute%next
    do while (.not.associated(a,target=v%lastAttribute))
       if (att(1:len(att)) .eq. a%att250(1:a%lena)) then
          buff250=ncf_getAttributeText(a)
          lenb=length(buff250,250,1)
          if (val(1:len(val)).eq.buff250(1:lenb)) then
             ncf_hasAttribute=.true.
          end if
          a=>v%lastAttribute
       else
          a=>a%next
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_hasAttribute

  character*250 function ncf_getAttributeText(a)
    implicit none
    type(attribute), pointer :: a
    integer ii, lens, lenb, length
    external length
    character*250 buff250,sum250
    character*20 :: myname = "ncf_getAttributeText"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    lens=0
    sum250=""
    do ii=1,a%len
       if (a%type.eq.nf_int1) then
          write(buff250,'(I20)')a%a1(ii) 
       else if (a%type.eq.nf_int2) then
          write(buff250,'(I20)')a%a2(ii) 
       else if (a%type.eq.nf_int) then
          write(buff250,'(I20)')a%a4(ii) 
       else if (a%type.eq.nf_int64) then
          write(buff250,'(I20)')a%a8(ii) 
       else if (a%type.eq.nf_real) then
          write(buff250,'(F12.5)')a%ar(ii) 
       else if (a%type.eq.nf_double) then
          write(buff250,'(F17.10)')a%ad(ii) 
       else if (a%type.eq.nf_char) then
          buff250=a%ac(ii) 
       end if
       call chop(buff250,25)
       lenb=length(buff250,25,10)
       if (lenb.eq.0) then
          sum250=sum250(1:lens)//" "
          lens=min(250,lens+1)
       else
          sum250=sum250(1:lens)//buff250(1:lenb)
          lens=min(250,lens+lenb)
       end if
    end do
    ncf_getAttributeText=sum250
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end function ncf_getAttributeText

  subroutine ncf_setAttributeText(a,att)
    implicit none
    type(attribute), pointer :: a
    character*(*) :: att
    integer ii,irc
    character*20 :: myname = "ncf_setAttributeText"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(a%ac)) deallocate(a%ac,stat=irc)
    if (associated(a%a1)) deallocate(a%a1,stat=irc)
    if (associated(a%a2)) deallocate(a%a2,stat=irc)
    if (associated(a%a4)) deallocate(a%a4,stat=irc)
    if (associated(a%a8)) deallocate(a%a8,stat=irc)
    if (associated(a%ar)) deallocate(a%ar,stat=irc)
    if (associated(a%ad)) deallocate(a%ad,stat=irc)
    a%type=nf_char
    a%len=len(att)
    allocate(a%ac(max(1,a%len)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate attribute.',len(att)
    end if
    do ii=1,a%len
       a%ac(ii)=att(ii:ii)
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_setAttributeText

  subroutine ncf_addAttribute(v,nam,val)
    implicit none
    type(variable), target :: v
    character*(*) :: nam
    character*(*) :: val
    type(attribute), pointer :: a, o
    integer irc
    integer ll
    logical :: found
    character*20 :: myname = "ncf_addAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ! find any old attribute with same name
    o=>v%firstAttribute%next
    found=(o%att250.eq.nam)
    do while (.not.found .and. .not.associated(o,target=v%lastAttribute))
       o=>o%next
       found=(o%att250.eq.nam)
    end do
    if (len(val).eq.0) then ! remove attribute
       if ( .not.associated(o,target=v%lastAttribute)) then ! remove o
          o%prev%next=>o%next
          o%next%prev=>o%prev
          if (associated(o%ac)) deallocate(o%ac)
          deallocate(o)
       end if
    else if (len(nam).ne.0) then ! add attribute
       allocate(a,stat=irc)
       if (irc.ne.0) then
          write(*,*) 'ncf_addAttribute Unable to allocate ATTRIBUTE.',irc
          return
       end if
       a%att250=nam
       a%lena=len(nam)
       a%type=nf_char
       a%len=len(val)
       allocate(a%ac(max(1,a%len)),stat=irc)
       if (irc.ne.0) then
          write(*,*)myname,' Unable to allocate.'
          return
       end if
       do ll=1,len(val)
          a%ac(ll)=val(ll:ll)
       end do
       a%v=>v
       if ( .not.associated(o,target=v%lastAttribute)) then ! replace o
          a%prev=>o%prev
          a%next=>o%next
          a%prev%next=>a
          a%next%prev=>a
          if (associated(o%ac)) deallocate(o%ac)
          deallocate(o)
       else ! add to end of attribute chain
          a%prev=>v%lastAttribute%prev
          a%next=>v%lastAttribute
          v%lastAttribute%prev%next=>a
          v%lastAttribute%prev=>a
       end if
       nullify(a)
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_addAttribute

  subroutine ncf_addGlobalAttribute(f,nam,val)
    implicit none
    type(inventory), target :: f
    character*(*) :: nam
    character*(*) :: val
    type(attribute), pointer :: a
    integer irc
    integer ll
    character*20 :: myname = "ncf_addGlobalAttribute"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    allocate(a,stat=irc)
    if (irc.ne.0) then
       write(*,*) 'ncf_addAttribute Unable to allocate ATTRIBUTE.',irc
       return
    end if
    a%att250=nam
    a%lena=len(nam)
    a%type=nf_char
    a%len=len(val)
    allocate(a%ac(max(1,a%len)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    do ll=1,len(val)
       a%ac(ll)=val(ll:ll)
    end do
    a%prev=>f%lastAttribute%prev
    a%next=>f%lastAttribute
    f%lastAttribute%prev%next=>a
    f%lastAttribute%prev=>a
    nullify(a)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_addGlobalAttribute

  !@#  DATA/GRID ******************************************************

  subroutine ncf_checkContents(f,par250,biok,irc)
    implicit none
    type(inventory) :: f
    character*250 par250,att250
    logical biok
    integer irc
    integer ii
    type(variable), pointer :: v=>null()
    type(variable), pointer :: x=>null()
    type(attribute), pointer :: a=>null()
    integer lenp, length
    external length 
    integer gridatt, leng,lena
    character*250 grid250
    character*20 :: myname = "ncf_checkContents"
    biok=.true.
    ! identify lat/lon/time/parameter variables
    ! grid definition attribute
    grid250="grid_mapping"
    call chop(grid250,250)
    leng=length(grid250,250,10)
    lenp=length(par250,250,10)
    nullify(a)
    v=>f%firstVariable%next
    do while (.not.associated(v,target=f%lastVariable))
       if (.not.associated(a)) then
          a=>ncf_getAttribute(v,grid250(1:leng),irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_getAttribute.',irc
             return
          end if
          ! find grid variable
          if (associated(a)) then
             att250=ncf_getAttributeText(a)
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
       v%fillAttribute=>ncf_getAttribute(v,"_FillValue",irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_getAttribute.',irc
          return
       end if
       if ( ncf_hasAttribute(v,"standard_name","projection_x_coordinate",irc) .and. .not.associated(f%xcoid) ) then
          f%xcoid=>v
       else if ( ncf_hasAttribute(v,"standard_name","projection_y_coordinate",irc) .and. .not.associated(f%ycoid) ) then
          f%ycoid=>v
       else if ( (ncf_hasAttribute(v,"standard_name","latitude",irc).or. &
            &     ncf_hasAttribute(v,"long_name","latitude",irc)) .and. .not.associated(f%latid) ) then
          f%latid=>v
       else if ( (ncf_hasAttribute(v,"standard_name","longitude",irc).or. &
            &     ncf_hasAttribute(v,"long_name","longitude",irc)) .and.  .not.associated(f%lonid)) then
          f%lonid=>v
       else if (ncf_hasAttribute(v,"standard_name","time",irc).or. &
            &     ncf_hasAttribute(v,"long_name","time",irc) .and. .not.associated(f%timid)) then
          f%timid=>v
       else if (ncf_hasAttribute(v,"standard_name","forecast_reference_time",irc) .and. .not.associated(f%terid)) then
          f%terid=>v
          if (v%len.ne.1) then
             write(*,*) 'checkContents Invalid reference_time length.',v%len
          end if
       else if (v%var250(1:v%lenv) .eq. par250(1:lenp)   .and. .not.associated(f%parid)) then
          f%parid=>v
       end if
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_hasAttribute.',irc
          return
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
  end subroutine ncf_checkContents

  subroutine ncf_checkParid(f,par250,biok,irc)
    implicit none
    type(inventory) :: f
    character*250 par250,att250
    logical biok
    integer irc
    integer ii
    type(variable), pointer :: v=>null()
    integer lenp, length
    external length 
    integer lena
    character*20 :: myname = "ncf_checkParid"
    biok=.true.
    ! identify lat/lon/time/parameter variables
    ! grid definition attribute
    lenp=length(par250,250,10)
    nullify(f%parid)
    nullify(f%latid)
    nullify(f%lonid)
    biok=.true.
    v=>f%firstVariable%next
    LOOP: do while (.not.associated(v,target=f%lastVariable))
       v%fillAttribute=>ncf_getAttribute(v,"_FillValue",irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_getAttribute.',irc
          return
       end if
       if (v%var250(1:v%lenv) .eq. par250(1:lenp)   .and. .not.associated(f%parid)) then
          f%parid=>v
       else if ( ncf_hasAttribute(v,"standard_name","projection_x_coordinate",irc) .and. .not.associated(f%xcoid) ) then
          f%xcoid=>v
       else if ( ncf_hasAttribute(v,"standard_name","projection_y_coordinate",irc) .and. .not.associated(f%ycoid) ) then
          f%ycoid=>v
       else if ( (v%var250(1:v%lenv).eq."latitude".or.&
            &     ncf_hasAttribute(v,"standard_name","latitude",irc).or. &
            &     ncf_hasAttribute(v,"long_name","latitude",irc).or. &
            &     ncf_hasAttribute(v,"standard_name","Lat",irc).or. &
            &     ncf_hasAttribute(v,"long_name","Lat",irc)) .and. .not.associated(f%latid) ) then
          f%latid=>v
       else if ( (v%var250(1:v%lenv).eq."longitude".or.&
            &     ncf_hasAttribute(v,"standard_name","longitude",irc).or. &
            &     ncf_hasAttribute(v,"long_name","longitude",irc).or. &
            &     ncf_hasAttribute(v,"standard_name","Lon",irc).or. &
            &     ncf_hasAttribute(v,"long_name","Lon",irc)) .and.  .not.associated(f%lonid)) then
          f%lonid=>v
       end if
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_hasAttribute.',irc
          return
       end if
       v=>v%next
    end do LOOP
    if (.not.associated(f%latid)) then
       write(*,*)'checkContents Missing LATITUDE in file: ', f%fn250(1:f%lenf)
       biok=.false.
    end if
    if (.not.associated(f%lonid)) then
       write(*,*)'checkContents Missing LONGITUDE in file: ', f%fn250(1:f%lenf)
       biok=.false.
    end if
    if (.not.associated(f%xcoid)) then
       f%xcoid=>f%lonid
    end if
    if (.not.associated(f%ycoid)) then
       f%ycoid=>f%latid
    end if
    if (lenp.ne.0.and..not.associated(f%parid)) then
       write(*,*)'checkContents Missing PARAMETER "', par250(1:lenp),'" in file: ',f%fn250(1:f%lenf)
       biok=.false.
    end if
    return
  end subroutine ncf_checkParid

  logical function ncf_hasNoData(v)
    implicit none
    type(variable),pointer :: v
    ncf_hasNoData=.true.
    if (.not.associated(v)) return
    if (v%len.ne.0.and.&
         & v%len1.eq.0.and.v%len2.eq.0.and.&
         & v%len4.eq.0.and.v%len8.eq.0.and.&
         & v%lenr.eq.0.and.v%lend.eq.0&
         & ) return
    ncf_hasNoData=.false.
    return
  end function ncf_hasNoData
  
  subroutine ncf_readRealData(v,biok,irc)
    implicit none
    type(variable),pointer :: v
    logical biok
    integer irc
    character*20 :: myname = "ncf_readRealData"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(v)) then
       if (ncf_bdeb) write(*,*)myname,' Variable:',v%var250(1:v%lenv)
       call ncf_readData(v,biok,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_readData.',irc
          return
       end if
       if (ncf_bdeb) write(*,*)myname,' Making.'
       call ncf_makeRealData(v,biok,irc)
       if (irc.ne.0) then
          write(*,*)myname,' Error return from ncf_makeRealData.',irc
          return
       end if
       ! write(*,*)'ncf_readRealData sample:',v%fd(1:min(v%lend,5))
    else
       write(*,*)'System error: attempt to read undefined variable.'
       biok=.false.
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_readRealData

  subroutine ncf_readData(v,biok,irc)
    implicit none
    type(variable),pointer :: v
    logical biok
    integer irc
    integer :: ii,ret
    integer, allocatable :: gsta(:),glim(:)
    character*20 :: myname = "ncf_readData"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (.not.associated(v)) then
       irc=704
       return
    end if
    if (allocated(gsta)) deallocate(gsta)
    if (allocated(glim)) deallocate(glim)
    allocate(gsta(max(1,v%nrdim)),glim(max(1,v%nrdim)),stat=irc)
    gsta(1)=1
    glim(1)=1
    do ii=1,v%nrdim
       gsta(ii)=v%sta(ii)
       glim(ii)=v%lim(ii) ! +v%sta(ii)-1
    end do
    if (ncf_bdeb) write(*,*)myname,'Read: ',v%varid,v%var250(1:v%lenv),v%nrdim," len=",v%len
    !
    v%fillc=char(nf_fill_char)
    v%fill1=nf_fill_int1
    v%fill2=nf_fill_int2
    v%fill4=nf_fill_int
    v%fill8=nf_fill_int64
    v%fillr=nf_fill_real
    v%filld=nf_fill_double
    if (ncf_bdeb)write(*,*)'ncf_readData Here...',v%var250(1:v%lenv),v%type
    ! allocate real data array
    !if (ncf_bdeb)write(*,*)'ncf_readData ',ncf_gettype(v%type),&
    !     & "  Reading data: ",v%var250(1:v%lenv)
    if (ncf_bdeb)write(*,*)'ncf_readData There...'
    if (.not.biok) return
    if (v%type.eq.nf_char) then
       if (v%lenc.ne.v%len) then
          v%lenc=v%len
          if (associated(v%fc)) deallocate(v%fc,stat=irc)
          if (ncf_bdeb)write(*,*)myname,'Allocating '//v%var250(1:v%lenv)//'%fc(',v%len,')'
          allocate(v%fc(max(1,v%lenc)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_readData Unable to allocate work array:',v%len
             irc=940
             return
          end if
       end if
       ret = nf_get_vara_TEXT(v%ncid,v%varid,gsta,glim,v%fc)
       if (ret .ne. NF_NOERR) then
          write(*,*)'ncf_readData Unable to load array from file.'
          irc=941
          return
       end if
       v%readyfield=.true.
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%ac)) then
             v%fillc=v%fillAttribute%ac(1)
             do ii=1,v%lenc
                if (v%fc(ii).eq.v%fillc) v%fc(ii)=" " ! nf_fill_char
             end do
          end if
       end if
    else if (v%type.eq.nf_int1) then
       if (v%len1.ne.v%len) then
          v%len1=v%len
          if (associated(v%f1)) deallocate(v%f1,stat=irc)
          allocate(v%f1(max(1,v%len1)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_readData Unable to allocate work array:',v%len
             irc=940
             return
          end if
       end if
       ret = nf_get_vara_INT1(v%ncid,v%varid,gsta,glim,v%f1)
       if (ret .ne. NF_NOERR) then
          write(*,*)'ncf_readData Unable to load array from file.'
          irc=941
          return
       end if
       v%readyfield=.true.
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%a1)) then
             v%fill1=v%fillAttribute%a1(1)
             do ii=1,v%len1
                if (v%f1(ii).eq.v%fill1) v%f1(ii)=nf_fill_int1
             end do
          end if
       end if
    else if (v%type.eq.nf_int2) then
       if (v%len2.ne.v%len) then
          v%len2=v%len
          if (associated(v%f2)) deallocate(v%f2,stat=irc)
          allocate(v%f2(max(1,v%len2)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_readData Unable to allocate work array:',v%len
             irc=942
             return
          end if
       end if
       ret = nf_get_vara_INT2(v%ncid,v%varid,gsta,glim,v%f2)
       if (ret .ne. NF_NOERR) then
          write(*,*)'ncf_readData Unable to load array from file.'
          irc=943
          return
       end if
       v%readyfield=.true.
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%a2)) then
             v%fill2=v%fillAttribute%a2(1)
             do ii=1,v%len2
                if (v%f2(ii).eq.v%fill2) v%f2(ii)=nf_fill_int2
             end do
          end if
       end if
    else if (v%type.eq.nf_int) then
       if (v%len4.ne.v%len) then
          v%len4=v%len
          if (associated(v%f4)) deallocate(v%f4,stat=irc)
          allocate(v%f4(max(1,v%len4)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_readData Unable to allocate work array:',v%len
             irc=944
             return
          end if
       end if
       !       write(*,*)'READFIELD: ',v%nrdim,v%len4,v%len
       !       if (v%nrdim.eq.0) then
       !          ret = nf_get_var_INT(v%ncid,v%varid,v%f4)
       !       else
       ret = nf_get_vara_INT(v%ncid,v%varid,gsta,glim,v%f4)
       !       end if
       if (ret .ne. NF_NOERR) then
          write(*,*)'ncf_readData Unable to load array from file.'
          irc=705
          return
       end if
       v%readyfield=.true.
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%a4)) then
             v%fill4=v%fillAttribute%a4(1)
             do ii=1,v%len4
                if (v%f4(ii).eq.v%fill4) v%f4(ii)=nf_fill_int
             end do
          end if
       end if
    else if (v%type.eq.nf_int64) then
       if (v%len8.ne.v%len) then
          v%len8=v%len
          if (associated(v%f8)) deallocate(v%f8,stat=irc)
          allocate(v%f8(max(1,v%len8)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_readData Unable to allocate work array:',v%len
             irc=944
             return
          end if
       end if
       !       write(*,*)'READFIELD: ',v%nrdim,v%len8,v%len
       !       if (v%nrdim.eq.0) then
       !          ret = nf_get_var_INT(v%ncid,v%varid,v%f4)
       !       else
       ret = nf_get_vara_INT64(v%ncid,v%varid,gsta,glim,v%f8)
       !       end if
       if (ret .ne. NF_NOERR) then
          write(*,*)'ncf_readData Unable to load array from file.'
          irc=705
          return
       end if
       v%readyfield=.true.
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%a8)) then
             v%fill8=v%fillAttribute%a8(1)
             do ii=1,v%len8
                if (v%f8(ii).eq.v%fill8) v%f8(ii)=nf_fill_int64
             end do
          end if
       end if
    else if (v%type.eq.nf_real) then
       !write(*,*)myname,'Real A:',associated(v),v%len,v%lenr,associated(v%fr)
       if (v%lenr.ne.v%len) then
          v%lenr=v%len
          if (associated(v%fr)) deallocate(v%fr,stat=irc)
          allocate(v%fr(max(1,v%lenr)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_readData Unable to allocate work array:',v%len
             irc=946
             return
          end if
       end if
       if (ncf_bdeb) write(*,*) myname,'nf_get_vara_real:'//v%var250(1:v%lenv)//": gsta=",&
            & gsta," glim=",glim," varid=",v%ncid,v%varid," array=",associated(v%fr),size(v%fr)
       !call ncf_printvariable(v)
       !write(*,*)myname,'Real B:',v%lenr, v%len,associated(v%fr),v%var250(1:v%lenv)
       ret = nf_get_vara_real(v%ncid,v%varid,gsta,glim,v%fr)
       if (ret .ne. NF_NOERR) then
          write(*,*) myname,'Variable=',v%var250(1:v%lenv),'   ',&
               & v%f%fn250(1:v%f%lenf),v%lenr,v%varid
          write(*,*) "Start=",v%sta," lim=", v%lim," siz=",v%siz," act=",size(v%fd)
          write(*,*)'ncf_readData Unable to load real array from file. ',nf_strerror(ret)
          irc=947
          return
       end if
       v%readyfield=.true.
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%ar)) then
             v%fillr=v%fillAttribute%ar(1)
             do ii=1,v%lenr
                if (v%fr(ii).eq.v%fillr) v%fr(ii)=nf_fill_real
             end do
          end if
       end if
    else if (v%type.eq.nf_double) then
       if (v%lend.ne.v%len) then
          v%lend=v%len
          if (associated(v%fd)) then
             !write(*,*)myname,'--Deallocating:'//v%var250(1:v%lenv)//ncf_getLabel(v)
             deallocate(v%fd,stat=irc)
          end if
          !write(*,*)myname,'++Allocating:'//v%var250(1:v%lenv),v%lend
          allocate(v%fd(max(1,v%lend)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_readData Unable to allocate work array:',v%len
             irc=948
             return
          end if
       end if
       if (ncf_bdeb) write(*,*) myname,'nf_get_vara_double:'//v%var250(1:v%lenv)//": gsta=",&
            & gsta," glim=",glim," varid=",v%ncid,v%varid," array=",associated(v%fd),size(v%fd)
       ret = nf_get_vara_double(v%ncid,v%varid,gsta,glim,v%fd)
       if (ret .ne. NF_NOERR) then
          write(*,*) 'Variable=',v%var250(1:v%lenv),'   ',v%f%fn250(1:v%f%lenf),v%lend
          write(*,*) 'Start, lim=',v%sta,v%lim,size(v%fd)
          write(*,*)'ncf_readData Unable to load double array from file.',nf_strerror(ret)
          irc=949
          return
       end if
       v%readyfield=.true.
       if (associated(v%fillAttribute)) then
          if (associated(v%fillAttribute%ad)) then
             v%filld=v%fillAttribute%ad(1)
             do ii=1,v%lend
                if (v%fd(ii).eq.v%filld) v%fd(ii)=nf_fill_double
             end do
          end if
       end if
    else
       v%readyfield=.false.
       write(*,*) 'Unknown data type:',v%type,' Variable=',v%var250(1:v%lenv),&
            & '   ',v%f%fn250(1:v%f%lenf),v%lend
       irc=345
       return
    end if
    deallocate(gsta,glim,stat=irc)
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_readData
  subroutine ncf_makeRealData(v,biok,irc)
    implicit none
    type(variable) :: v
    logical biok
    integer irc
    integer ii
    character*20 :: myname = "ncf_makeRealData"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ! allocate real data array
    ! write(*,*)'ncf_makeRealData ',ncf_gettype(v%type),"  making real data: ",v%var250(1:v%lenv)
    if (.not.biok) return
    if (v%type.ne.nf_double) then
       if (v%lend.ne.v%len) then
          v%lend=v%len
          if (associated(v%fd)) then
             !write(*,*)myname,'--Deallocating:'//v%var250(1:v%lenv)
             deallocate(v%fd,stat=irc)
          end if
          !write(*,*)myname,'++Allocating:'//v%var250(1:v%lenv),v%lend
          allocate(v%fd(max(1,v%lend)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_makeRealData Unable to allocate work array:',v%len
             irc=706
             return
          end if
       end if
    end if
    if (v%type.eq.nf_int1) then
       do ii=1,v%len1
          if (v%f1(ii) .eq. v%fill1) then
             v%fd(ii)=v%filld
          else
             v%fd(ii)=v%f1(ii)
          end if
       end do
       deallocate(v%f1,stat=irc)
       if(irc.ne.0) then
          write(*,*)'ncf_makeRealData Unable to deallocate work array.'
          irc=707
          return
       end if
       nullify(v%f1)
       v%len1=0
    else if (v%type.eq.nf_int2) then
       do ii=1,v%len2
          if (v%f2(ii) .eq. v%fill2) then
             v%fd(ii)=v%filld
          else
             v%fd(ii)=v%f2(ii)
          end if
       end do
       deallocate(v%f2,stat=irc)
       if(irc.ne.0) then
          write(*,*)'ncf_makeRealData Unable to deallocate work array.'
          irc=708
          return
       end if
       nullify(v%f2)
       v%len2=0
    else if (v%type.eq.nf_int) then
       do ii=1,v%len4
          if (v%f4(ii) .eq. v%fill4) then
             v%fd(ii)=v%filld
          else
             v%fd(ii)=v%f4(ii)
          end if
       end do
       deallocate(v%f4,stat=irc)
       if(irc.ne.0) then
          write(*,*)'ncf_makeRealData Unable to deallocate work array.'
          irc=709
          return
       end if
       nullify(v%f4)
       v%len4=0
    else if (v%type.eq.nf_int64) then
       do ii=1,v%len8
          if (v%f8(ii) .eq. v%fill8) then
             v%fd(ii)=v%filld
          else
             v%fd(ii)=v%f8(ii)
          end if
       end do
       deallocate(v%f8,stat=irc)
       if(irc.ne.0) then
          write(*,*)'ncf_makeRealData Unable to deallocate work array.'
          irc=709
          return
       end if
       nullify(v%f8)
       v%len8=0
    else if (v%type.eq.nf_real) then
       do ii=1,v%lenr
          if (v%fr(ii) .eq. v%fillr) then
             v%fd(ii)=v%filld
          else
             v%fd(ii)=v%fr(ii)
          end if
       end do
       deallocate(v%fr,stat=irc)
       if(irc.ne.0) then
          write(*,*)'ncf_makeRealData Unable to deallocate work array.'
          irc=710
          return
       end if
       nullify(v%fr)
       v%lenr=0
    else if (v%type.eq.nf_double) then
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_makeRealData

  subroutine ncf_unmakeRealData(v,irc)
    implicit none
    type(variable) :: v
    integer irc
    integer ii
    character*20 :: myname = "ncf_unmakeRealData"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ! allocate real data array
    if (associated(v%fd)) then
       if (v%type.eq.nf_int1) then
          v%len1=v%lend
          if (associated(v%f1)) deallocate(v%f1,stat=irc)
          allocate(v%f1(max(1,v%len1)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to allocate work array:',v%len1
             irc=711
             return
          end if
          do ii=1,v%len1
             if (v%fd(ii) .eq. v%filld) then
                v%f1(ii)=v%fill1
             else
                v%f1(ii)=nint(v%fd(ii))
             end if
          end do
          !write(*,*)myname,'--deallocating:'//v%var250(1:v%lenv),v%lend
          deallocate(v%fd,stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to deallocate work array.'
             irc=712
             return
          end if
          nullify(v%fd)
          v%lend=0
       else if (v%type.eq.nf_int2) then
          v%len2=v%lend
          if (associated(v%f2)) deallocate(v%f2,stat=irc)
          allocate(v%f2(max(1,v%len2)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to allocate work array:',v%len2
             irc=713
             return
          end if
          do ii=1,v%len2
             if (v%fd(ii) .eq. v%filld) then
                v%f2(ii)=v%fill2
             else
                v%f2(ii)=nint(v%fd(ii))
             end if
          end do
          !write(*,*)myname,'--deallocating:'//v%var250(1:v%lenv),v%lend
          deallocate(v%fd,stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to deallocate work array.'
             irc=714
             return
          end if
          nullify(v%fd)
          v%lend=0
       else if (v%type.eq.nf_int) then
          v%len4=v%lend
          if (associated(v%f4)) deallocate(v%f4,stat=irc)
          allocate(v%f4(max(1,v%len4)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to allocate work array:',v%len4
             irc=715
             return
          end if
          do ii=1,v%len4
             if (v%fd(ii) .eq. v%filld) then
                v%f4(ii)=v%fill4
             else
                v%f4(ii)=nint(v%fd(ii))
             end if
          end do
          !write(*,*)myname,'--deallocating:'//v%var250(1:v%lenv),v%lend
          deallocate(v%fd,stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to deallocate work array.'
             irc=716
             return
          end if
          nullify(v%fd)
          v%lend=0
       else if (v%type.eq.nf_int64) then
          v%len8=v%lend
          if (associated(v%f8)) deallocate(v%f8,stat=irc)
          allocate(v%f8(max(1,v%len8)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to allocate work array:',v%len4
             irc=715
             return
          end if
          do ii=1,v%len8
             if (v%fd(ii) .eq. v%filld) then
                v%f8(ii)=v%fill8
             else
                v%f8(ii)=nint(v%fd(ii))
             end if
          end do
          !write(*,*)myname,'--deallocating:'//v%var250(1:v%lenv),v%lend
          deallocate(v%fd,stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to deallocate work array.'
             irc=716
             return
          end if
          nullify(v%fd)
          v%lend=0
       else if (v%type.eq.nf_real) then
          v%lenr=v%lend
          if (associated(v%fr)) deallocate(v%fr,stat=irc)
          allocate(v%fr(max(1,v%lenr)),stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to allocate work array:',v%lenr
             irc=717
             return
          end if
          do ii=1,v%lenr
             if (v%fd(ii) .eq. v%filld) then
                v%fr(ii)=v%fillr
             else
                v%fr(ii)=v%fd(ii)
             end if
          end do
          !write(*,*)myname,'--deallocating:'//v%var250(1:v%lenv),v%lend
          deallocate(v%fd,stat=irc)
          if(irc.ne.0) then
             write(*,*)'ncf_unmakeRealData Unable to deallocate work array.'
             irc=718
             return
          end if
          nullify(v%fd)
          v%lend=0
       else if (v%type.eq.nf_double) then
       end if
    else
       write(*,*)'UNMAKEREALDATA No real data was found.'
       irc=944
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_unmakeRealData

  subroutine ncf_clearData(f,irc)
    implicit none
    type(inventory) :: f
    integer irc
    type(variable), pointer  :: v
    integer ii
    character*20 :: myname = "ncf_clearData"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ! loop over variables, remove data (keep attributes etc)
    if (f%initialised) then
       v=>f%firstVariable%next
       do while(.not.associated(v,target=f%lastVariable))
          call ncf_clearVariable(v,irc)
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_clearVariable.',irc
             return
          end if
          v=>v%next
       end do
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_clearData
  !
  subroutine ncf_ignoreVariable(v)
    implicit none
    type(variable),pointer :: v
    if (associated(v)) v%ignored=.true.
    return
  end subroutine ncf_ignoreVariable
  !
  subroutine ncf_clearVariable(v,irc)
    implicit none
    type(variable),pointer :: v
    integer irc
    character*20 :: myname = "ncf_clearVariable"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(v)) then
       !if (ncf_bdeb)
       !write(*,*)myname,"Var:'"//v%var250(1:v%lenv)//"'",v%nrdim,associated(v%ind),v%len
       if (associated(v%fc)) then
          deallocate(v%fc,stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,'Unable to deallocate fc.',irc,v%var250(1:v%lenv)
          end if
          nullify(v%fc)
       end if
       if (associated(v%f1)) then
          deallocate(v%f1,stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,'Unable to deallocate f1.',irc,v%var250(1:v%lenv)
          end if
          nullify(v%f1)
       end if
       if (associated(v%f2)) then
          deallocate(v%f2,stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,'Unable to deallocate f2.',irc,v%var250(1:v%lenv)
          end if
          nullify(v%f2)
       end if
       if (associated(v%f4)) then
          deallocate(v%f4,stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,'Unable to deallocate f4.',irc,v%var250(1:v%lenv)
          end if
          nullify(v%f4)
       end if
       if (associated(v%f8)) then
          deallocate(v%f8,stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,'Unable to deallocate f8.',irc,v%var250(1:v%lenv)
          end if
          nullify(v%f8)
       end if
       if (associated(v%fr)) then
          deallocate(v%fr,stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,'Unable to deallocate fr.',irc,v%var250(1:v%lenv)
          end if
          nullify(v%fr)
       end if
       if (associated(v%fd)) then
          !write(*,*)myname,'--Deallocating:'//v%var250(1:v%lenv),&
           !       & associated(v%fd), v%lend
          deallocate(v%fd,stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,'Unable to deallocate fd.',irc,v%var250(1:v%lenv)
          end if
          !write(*,*)myname,'--Deallocated:'//v%var250(1:v%lenv),associated(v%fd)
          nullify(v%fd)
       end if
       v%lenc=0
       v%len1=0
       v%len2=0
       v%len4=0
       v%len8=0
       v%lenr=0
       v%lend=0
       v%readyfield=.false.
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_clearVariable

  !@#  FIELD ******************************************************

  subroutine ncf_copyField(v,c)
    implicit none
    type(variable), pointer :: v
    type(variable), pointer :: c
    integer ll
    character*20 :: myname = "ncf_copyField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(v%fd).and.associated(c%fd)) then
       do ll=1,v%lend
          v%fd(ll)=c%fd(ll)
       end do
    else
       write(*,*)myname,' Fields not allocated.',associated(v%fd),associated(c%fd)
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_copyField

  subroutine ncf_copyInventoryField(f,s,irc)
    implicit none
    type(inventory),pointer :: f
    type(variable), pointer :: s
    integer irc
    character*20 :: myname = "ncf_copyInventoryField"
    type(variable), pointer :: v
    logical :: bok
    v=>ncf_getVariable(f,s%var250(1:s%lenv),bok,irc)
    if (.not.bok) then
       write(*,*)myname,' Variable not found:',s%var250(1:s%lenv)
       irc=394
       return
    else if (irc.ne.0) then
       write(*,*)myname,' Error return from ncf_getVariable.',irc
       return
    end if
    call ncf_copyVariableField(v,s,irc)
    if (irc.ne.0) then
       write(*,*)myname,' Error return from ncf_copyVariableField.',irc
       return
    end if
    return
  end subroutine ncf_copyInventoryField

  subroutine ncf_copyVariableField(v,s,irc)
    implicit none
    type(variable), pointer :: s,v
    integer irc
    integer jj
    character*20 :: myname = "ncf_copyVariableField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (s%initialised) then
       v%len=s%len
       v%lens=s%lens
       v%lenc=s%lenc
       v%len1=s%len1
       v%len2=s%len2
       v%len4=s%len4
       v%len8=s%len8
       v%lenr=s%lenr
       v%lend=s%lend
       v%type=s%type
       !write(*,*)myname,'Copying data:',&
       !     & v%var250(1:v%lenv),v%readyfield
       if (v%lenc.gt.0) then
          if (associated(v%fc)) deallocate(v%fc)
          allocate(v%fc(max(1,v%lenc)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          if (associated(s%fc)) then
             do jj=1,v%lenc
                v%fc(jj)=s%fc(jj)
             end do
          end if
       end if
       if (v%len1.gt.0) then
          if (associated(v%f1)) deallocate(v%f1)
          allocate(v%f1(max(1,v%len1)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          if (associated(s%f1)) then
             do jj=1,v%len1
                v%f1(jj)=s%f1(jj)
             end do
          end if
       end if
       if (v%len2.gt.0) then
          if (associated(v%f2)) deallocate(v%f2)
          allocate(v%f2(max(1,v%len2)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          if (associated(s%f2)) then
             do jj=1,v%len2
                v%f2(jj)=s%f2(jj)
             end do
          end if
       end if
       if (v%len4.gt.0) then
          if (associated(v%f4)) deallocate(v%f4)
          allocate(v%f4(max(1,v%len4)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          if (associated(s%f4)) then
             do jj=1,v%len4
                v%f4(jj)=s%f4(jj)
             end do
          end if
       end if
       if (v%len8.gt.0) then
          if (associated(v%f8)) deallocate(v%f8)
          allocate(v%f8(max(1,v%len8)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          if (associated(s%f8)) then
             do jj=1,v%len8
                v%f8(jj)=s%f8(jj)
             end do
          end if
       end if
       if (v%lenr.gt.0) then
          if (associated(v%fr)) deallocate(v%fr)
          allocate(v%fr(max(1,v%lenr)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          if (associated(s%fr)) then
             do jj=1,v%lenr
                v%fr(jj)=s%fr(jj)
             end do
          end if
       end if
       if (v%lend.gt.0) then
          if (associated(v%fd)) then
             !write(*,*)myname,'--Deallocating:'//v%var250(1:v%lenv)
             deallocate(v%fd)
          end if
          !write(*,*)myname,'++Allocating:'//v%var250(1:v%lenv),v%lend
          allocate(v%fd(max(1,v%lend)),stat=irc)
          if (irc.ne.0) then
             write(*,*)myname,' Unable to allocate.'
             return
          end if
          if (associated(s%fd)) then
             do jj=1,v%lend
                v%fd(jj)=s%fd(jj)
             end do
          end if
       end if
       v%readyfield=s%readyfield
       v%initialised=.true.
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_copyVariableField

  function ncf_insertField(f,nam,i,irc)
    implicit none
    type(variable), pointer ::  ncf_insertField
    type(inventory),target :: f
    character*(*) nam
    character*250 nam250
    type(dimensionOrder) :: i
    integer irc
    !
    type(variable), pointer :: v
    type(attribute),pointer :: a,va
    integer length, lenn, kk
    type(dimension), pointer :: d
    external length
    character*20 :: myname = "ncf_insertField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    !
    nam250=nam(1:len(nam))
    call chop(nam250,250)
    lenn=length(nam250,250,10)
    allocate(ncf_insertField,stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_insertField%f=>f
    ncf_insertField%ncid=f%ncid
    ncf_insertField%varid=0
    ncf_insertField%var250=nam250
    ncf_insertField%lenv=lenn
    ncf_insertField%type=nf_double
    ncf_insertField%nrdim=i%nrdim
    ncf_insertField%nratt=0
    ncf_insertField%scale=1.0D0
    allocate(ncf_insertField%ind(max(1,ncf_insertField%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)'ncf_insertField Unable to allocate:',ncf_insertField%nrdim
       return
    end if
    allocate(ncf_insertField%sta(max(1,ncf_insertField%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)'ncf_insertField Unable to allocate:',ncf_insertField%nrdim
       return
    end if
    allocate(ncf_insertField%lim(max(1,ncf_insertField%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)'ncf_insertField Unable to allocate:',ncf_insertField%nrdim
       return
    end if
    allocate(ncf_insertField%siz(max(1,ncf_insertField%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)'ncf_insertField Unable to allocate:',ncf_insertField%nrdim
       return
    end if
    ncf_insertField%len=1
    ncf_insertField%lens=1
    kk=0
    d=>i%firstDimension%next
    do while (.not.associated(d,target=i%lastDimension))
       if (d%ind.lt.1.or.d%ind.gt.f%nrdim) then
          write(*,*)'ncf_insertField Invalid dimension index:',kk,d%ind,f%nrdim
          irc=834
          return
       end if
       kk=kk+1
       ncf_insertField%ind(kk)=d%ind
       ncf_insertField%sta(kk)=d%sta
       ncf_insertField%lim(kk)=d%lim
       ncf_insertField%siz(kk)=d%siz
       ncf_insertField%len=ncf_insertField%len*ncf_insertField%lim(kk)
       ncf_insertField%lens=ncf_insertField%lens*ncf_insertField%siz(kk)
       d=>d%next
    end do
    if (kk.ne.ncf_insertField%nrdim) then
       write(*,*)'ncf_insertField Dimension mismatch:',kk,ncf_insertField%nrdim
       d=>i%firstDimension%next
       do while (.not.associated(d,target=i%lastDimension))
          write(*,*) 'ncf_insertField Dimension:',d%ind, d%lim
          d=>d%next
       end do
       irc=719
       return
    end if
    ncf_insertField%lenc=0
    ncf_insertField%len1=0
    ncf_insertField%len2=0
    ncf_insertField%len4=0
    ncf_insertField%len8=0
    ncf_insertField%lenr=0
    ncf_insertField%lend=ncf_insertField%len
    allocate(ncf_insertField%fd(max(1,ncf_insertField%lend)), stat=irc)
    if (irc.ne.0) then
       write(*,*)'ncf_insertField Unable to allocate:',ncf_insertField%lend
       return
    end if
    do kk=1,ncf_insertField%lend
       ncf_insertField%fd(kk)=0.0D0
    end do
    allocate(ncf_insertField%firstAttribute,ncf_insertField%lastAttribute, stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    ncf_insertField%firstAttribute%next=>ncf_insertField%lastAttribute
    ncf_insertField%lastAttribute%prev=>ncf_insertField%firstAttribute
    ! insert into inventory
    ncf_insertField%prev=>f%lastVariable%prev
    ncf_insertField%next=>f%lastVariable
    f%lastVariable%prev%next=>ncf_insertField
    f%lastVariable%prev=>ncf_insertField
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end function ncf_insertField

  subroutine ncf_allocateField(v,irc)
    implicit none
    type(variable),pointer :: v
    real :: val
    integer :: irc
    integer :: ii
    character*20 :: myname = "ncf_allocateField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ! get field length
    v%lend=ncf_getVariableLength(v)
    ! allocate field
    if (associated(v%fd)) then
       !write(*,*)myname,'--Deallocating:'//v%var250(1:v%lenv)
       deallocate(v%fd)
    end if
    !write(*,*)myname,'++Allocating:'//v%var250(1:v%lenv),v%lend
    allocate(v%fd(max(1,v%lend)), stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_allocateField

  subroutine ncf_initField(v,val,irc)
    implicit none
    type(variable),pointer :: v
    real :: val
    integer :: irc
    integer :: ii
    character*20 :: myname = "ncf_initField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ! get field length
    v%lend=ncf_getVariableLength(v)
    ! allocate field
    if (associated(v%fd)) then
       !write(*,*)myname,'--Deallocating:'//v%var250(1:v%lenv)
       deallocate(v%fd)
    end if
    !write(*,*)myname,'++Allocating:'//v%var250(1:v%lenv),v%lend
    allocate(v%fd(max(1,v%lend)), stat=irc)
    if (irc.ne.0) then
       write(*,'(X,A,A,I0,A)')myname,&
            & " Unable to allocate '"//v%var250(1:v%lenv)//"(",v%lend,")'"
       return
    end if
    do ii=1,v%lend
       v%fd(ii)=val
    end do
    v%readyfield=.true.
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_initField

  subroutine ncf_addField(v,c)
    implicit none
    type(variable) :: v
    type(variable) :: c
    integer ll
    character*20 :: myname = "ncf_addField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=v%fd(ll)+c%fd(ll)
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_addField

  subroutine ncf_subtractField(v,c)
    implicit none
    type(variable), pointer :: v
    type(variable), pointer :: c
    integer ll
    character*20 :: myname = "ncf_subtractField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=v%fd(ll)-c%fd(ll)
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_subtractField

  subroutine ncf_multiplyField(v,c)
    implicit none
    type(variable), pointer :: v
    type(variable), pointer :: c
    integer ll
    character*20 :: myname = "ncf_multiplyField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=v%fd(ll)*c%fd(ll)
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_multiplyField

  subroutine ncf_scaleField(v,scale)
    implicit none
    type(variable), pointer :: v
    real scale
    integer ll
    character*20 :: myname = "ncf_scaleField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=v%fd(ll)*scale
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_scaleField

  subroutine ncf_divideField(v,c)
    implicit none
    type(variable), pointer :: v
    type(variable), pointer :: c
    integer ll
    character*20 :: myname = "ncf_divideField"
    if (ncf_bdeb) write(*,*)myname,' Entering.',v%lend,c%lend
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double .and. abs(c%fd(ll)).gt.1.0D-10) then
          v%fd(ll)=v%fd(ll)/c%fd(ll)
       else
          v%fd(ll)=nf_fill_double
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_divideField

  subroutine ncf_sqrtField(v)
    implicit none
    type(variable), pointer :: v
    integer ll
    character*20 :: myname = "ncf_sqrtField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    do ll=1,v%lend
       v%fd(ll)=sqrt(max(0.0D0,v%fd(ll)))
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_sqrtField

  subroutine ncf_roundField(v,prec)
    implicit none
    type(variable),pointer  :: v
    real prec
    integer ll
    real xprec
    character*20 :: myname = "ncf_roundField"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    xprec=max(1.0D-10,prec)
    do ll=1,v%lend
       if (v%fd(ll).ne.nf_fill_double) then
          v%fd(ll)=DBLE(nint(v%fd(ll)/xprec))*xprec
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_roundField

  !@#  WRITE TO FILE ******************************************************

  subroutine ncf_writeNcOpen(out,outnc250,irc)
    implicit none
    character*250  outnc250
    type(inventory),pointer :: out
    integer irc
    integer length, leno
    external length
    character*14 :: myname = "ncf_writeNcOpen"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(out)) then
       call ncf_createNcFile(out,outnc250,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_createNcFile.',irc
          return
       end if
       call ncf_writeNcDims(out,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_writeNcDims.',irc
          return
       end if
       call ncf_writeNcVars(out,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_writeNcVars.',irc
          return
       end if
       call ncf_writeNcAtts(out,out%firstAttribute,out%lastAttribute,irc) ! global attributes
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_writeNcAtts.',irc
          return
       end if
       call ncf_endNcDef(out,irc)
       if (irc.ne.0) then
          write(*,*)myname,'Error return from ncf_endNcDef.',irc
          return
       end if
       call ncf_writeNcData(out,irc)
       if (associated(out)) then
          if (irc.ne.0) then
             write(*,*)myname,'Error return from ncf_writeNcData.',irc
             return
          end if
       end if
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_writeNcOpen
  !
  subroutine ncf_writeNcOut(out,outnc250,irc)
    implicit none
    character*250  outnc250
    type(inventory),pointer :: out
    integer irc
    integer length, leno
    external length
    character*14 :: myname = "ncf_writeNcOut"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(out)) then
       call ncf_writeNcOpen(out,outnc250,irc)
       if (irc.ne.0) then
          write(*,*)myname,"Error return from writeNcOpen.",irc
          return
       end if
       call ncf_closeNcFile(out,irc)
       if (irc.ne.0) then
          write(*,*)myname,"Error return from closeNcFile.",irc
          return
       end if
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
    return
  end subroutine ncf_writeNcOut

  subroutine ncf_createNcFile(out,outnc250,irc)
    implicit none
    type(inventory),pointer :: out
    character*250  outnc250
    integer irc
    integer length,leno, ret
    external length
    character*14 :: myname = "ncf_createNcFile"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    call chop(outnc250,250)
    leno=length(outnc250,250,10)
    if (.not.associated(out)) then
       irc=945
       write(*,*)myname,' Invalid inventory: ',outnc250(1:leno)
       return
    end if
    if (leno.eq.0) then
       irc=946
       write(*,*)myname,' Zero-length file name.'
       return
    end if
    out%ncid=0
    out%fn250=outnc250
    out%lenf=leno
    write(*,*)myname,'  Opening: ',out%fn250(1:out%lenf), out%lenf
    ! or(NF_NOCLOBBER,NF_SHARE) - do not overwrite existing file but give error message
    !    ret=NF_CREATE(out%fn250(1:out%lenf),NF_SHARE,out%ncid) ! create netCDF dataset: enter define mode
    ret=NF_CREATE(out%fn250(1:out%lenf),NF_NETCDF4,out%ncid) ! create netCDF dataset: enter define mode
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_CREATE:",nf_strerror(ret),ret,"( "//outnc250(1:leno)//" )"
       irc=790
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_createNcFile

  subroutine ncf_writeNcDims(out,irc)
    implicit none
    type(inventory),pointer :: out
    integer irc
    integer ii,ret,leng,length
    external length
    character*14 :: myname = "ncf_writeNcDims"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (associated(out%dimid)) deallocate(out%dimid,stat=irc)
    allocate(out%dimid(max(1,out%nrdim)),stat=irc)
    if (irc.ne.0) then
       write(*,*)myname,' Unable to allocate.'
       return
    end if
    do ii=1,out%nrdim
       leng=length(out%dim250(ii),250,10)
       if (out%unlimdimid.eq.ii.and.ii.eq.out%nrdim) then ! only "last" dimension can be unlimited
          ret=NF_DEF_DIM(out%ncid,out%dim250(ii)(1:leng),nf_unlimited,out%dimid(ii))
       else
          ret=NF_DEF_DIM(out%ncid,out%dim250(ii)(1:leng),out%pos%siz(ii),out%dimid(ii))
       end if
       if (ret .ne. NF_NOERR) then
          write(*,*) myname,"ERROR from NF_DEF_DIM:",nf_strerror(ret),out%ncid
          irc=790
          return
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine NCF_WRITENCDIMS
  !
  subroutine ncf_writeNcVars(out,irc)
    implicit none
    type(inventory), pointer :: out
    integer irc
    integer ret
    type(variable), pointer :: v
    integer gdims(out%nrdim)
    character*250 :: buff250
    integer ii,jj,len,siz,lenb
    integer, external :: length
    character*16 :: myname = "ncf_writeNcVars"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    write(*,*)myname,'  Defining: ',out%fn250(1:out%lenf), out%lenf
    v=>out%firstVariable%next
    do while(.not.associated(v,target=out%lastVariable))
       if (.not.v%ignored) then
          len=1
          siz=1
          do ii=1,v%nrdim
             if (v%ind(ii).gt.out%nrdim) then
                write(*,*)myname,'Invalid dim:',v%ind(ii)
                call ncf_printinventory(out)
                call ncf_printvariable(v)
                irc=998
                return
             end if
             gdims(ii)=out%dimid(v%ind(ii))
             len=len*out%pos%lim(v%ind(ii))
             siz=siz*out%pos%siz(v%ind(ii))
          end do
          write(buff250,'(X,A,I3,X,A,I0,A,I0,A)')'Variable: ',&
               & v%varid,v%var250(1:v%lenv)//'->'//v%out250(1:v%leno)//&
               & ' (',len,'|',siz,') '//ncf_gettype(v%type)
          call chop0(buff250,250)
          lenb=length(buff250,250,20)
          write(*,*)myname,buff250(1:lenb)
          if (v%lout) then
             ret=NF_DEF_VAR(out%NCID,v%out250(1:v%leno),v%type,v%nrdim,gdims,v%varid)
          else
             ret=NF_DEF_VAR(out%NCID,v%var250(1:v%lenv),v%type,v%nrdim,gdims,v%varid)
          end if
          if (ret .ne. NF_NOERR) then
             write(*,*) myname,"ERROR from NF_DEF_VAR:",trim(nf_strerror(ret)),out%ncid,v%var250(1:v%lenv)
             irc=790
             return
          end if
          call ncf_writeNcAtts(out,v%firstAttribute,v%lastAttribute,irc)
          if (irc.ne.0) then
             write(*,*)myname,' Error return from ncf_writeNcAtts.',irc
             return
          end if
          ret = NF_DEF_VAR_deflate(out%NCID,v%varid, 0, 1, 3) ! 4 works well also
          if (ret .ne. NF_NOERR) then
             write(*,*) myname,"ERROR from NF_DEF_VAR:",nf_strerror(ret),out%ncid
             irc=790
             return
          end if
       end if
       v=>v%next
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_writeNcVars

  subroutine ncf_writeNcAtts(out,first,last,irc)
    implicit none
    type(inventory), pointer :: out
    type(attribute), pointer :: first,last
    integer irc
    type(attribute), pointer :: a
    integer varid,ret
    character*14 :: myname = "ncf_writeNcAtts"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    a=>first%next
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
       else if (a%type.eq.nf_int64) then
          ret=NF_NOERR
          notimp=notimp+1
          if(notimp.eq.1)write(*,*)myname,'Int64 not implemented...',a%att250(1:a%lena)
!          ret=NF_PUT_ATT_INT64(OUT%NCID, varid,a%att250(1:a%lena),a%type,a%LEN, a%a8)
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
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_writeNcAtts

  subroutine ncf_endNcDef(out,irc)
    implicit none
    type(inventory), pointer :: out
    integer irc
    integer ret
    character*14 :: myname = "ncf_endNcDef"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    ret=NF_ENDDEF(out%ncid)       ! end definitions: leave define mode
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_ENDDEF:",nf_strerror(ret),out%ncid
       irc=790
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_endNcDef

  subroutine ncf_writeNcData(out,irc)
    implicit none
    type(inventory), pointer :: out
    integer irc
    integer ret
    type(variable), pointer :: v
    integer, allocatable :: gsta(:),glim(:)
    integer :: ii
    character*16 :: myname = "ncf_writeNcData"
    if (ncf_bdeb) write(*,*)myname,'Entering.'
    write(*,*)myname,'Updating: ',out%fn250(1:out%lenf), out%lenf
    v=>out%firstVariable%next
    do while(.not.associated(v,target=out%lastVariable))
       if (v%readyfield.and. .not.v%ignored.and.v%writefield) then
          if (ncf_bdeb)write(*,'(X,A,A,X,I3,X,A,X,I0,"|",I0,X,L1)')&
               & myname,'Writing:',v%varid,&
               & v%var250(1:v%lenv)//"  ("//ncf_gettype(v%type)//")",&
               & v%len,v%lens,v%readyfield
          allocate(gsta(max(1,v%nrdim)),glim(max(1,v%nrdim)),stat=irc)
          gsta(1)=1
          glim(1)=1
          do ii=1,v%nrdim
             gsta(ii)=v%sta(ii)
             glim(ii)=v%lim(ii) ! +v%sta(ii) -1
          end do
          if (.not. v%readyfield) then
             write(*,*)myname,'No data for: ',v%var250(1:v%lenv),ncf_gettype(v%type)
          else if (v%type.eq.nf_char.and.associated(v%fc)) then
             ret= NF_PUT_VARA_TEXT(out%ncid,v%varid,gsta,glim,v%fc)
          else if (v%type.eq.nf_int1.and.associated(v%f1)) then
             ret= NF_PUT_VARA_INT1(out%ncid,v%varid,gsta,glim,v%f1)
          else if (v%type.eq.nf_int2.and.associated(v%f2)) then
             !write(*,*)myname,'Int2:',out%ncid,size(v%f2),v%f2(size(v%f2)/2)
             ret= NF_PUT_VARA_INT2(out%ncid,v%varid,gsta,glim,v%f2)
          else if (v%type.eq.nf_int.and.associated(v%f4)) then
             ret= NF_PUT_VARA_INT(out%ncid,v%varid,gsta,glim,v%f4)
          else if (v%type.eq.nf_int.and.associated(v%f8)) then
             ret= NF_PUT_VARA_INT64(out%ncid,v%varid,gsta,glim,v%f8)
          else if (v%type.eq.nf_real.and.associated(v%fr)) then
             if (ncf_bdeb)write(*,*)myname,'put_vara_float:',&
                  & out%ncid,v%varid," sta=",gsta,&
                  & " lim=",glim,associated(v%fr),size(v%fr),v%lend
             ret= NF_PUT_VARA_REAL(out%ncid,v%varid,gsta,glim,v%fr)
          else if (v%type.eq.nf_double.and.associated(v%fd)) then
             if (ncf_bdeb)write(*,*)myname,'put_vara_double:',&
                  & out%ncid,v%varid," sta=",gsta,&
                  & " lim=",glim,associated(v%fd),size(v%fd),v%lend
             ret= NF_PUT_VARA_DOUBLE(out%ncid,v%varid,gsta,glim,v%fd)
          else
             write(*,*)myname,'Undefined variable data: ',v%var250(1:v%lenv),v%type, &
                  &nf_double,nf_real,nf_int,nf_int2,nf_int1, &
                  &associated(v%fd),associated(v%fr),associated(v%f4),&
                  & associated(v%f8),associated(v%f2),associated(v%f1)
             ret=nf_noerr
          end if
          if (ret .ne. NF_NOERR) then
             write(*,*)myname,'sta:',v%sta,'lim:',v%lim,&
                  & size(v%fc),size(v%f1),size(v%f2), &
                  & size(v%f4),size(v%f8),size(v%fr),size(v%fd),v%varid
             call ncf_printvariable(v)
             write(*,*) myname,"ERROR from NF_PUT_VARA_*:",nf_strerror(ret)
             irc=812
             return
          end if
          deallocate(gsta,glim,stat=irc)
       end if
       v=>v%next
    end do
    if (ncf_bdeb) write(*,*)myname,'Done.',irc
  end subroutine ncf_writeNcData

  subroutine ncf_closeNcFile(out,irc)
    implicit none
    type(inventory), pointer :: out
    integer irc
    integer ret
    character*14 :: myname = "ncf_closeNcFile"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    write(*,*)myname,'  Closing: ',out%fn250(1:out%lenf), out%lenf
    ret=NF_CLOSE(out%ncid)        ! end definitions: leave define mode
    if (ret .ne. NF_NOERR) then
       write(*,*) myname,"ERROR from NF_CLOSE:",nf_strerror(ret),out%ncid
       irc=790
       return
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.',irc
  end subroutine ncf_closeNcFile

  !@#  DEBUG/PRINT ******************************************************

  subroutine ncf_printInventory(i)
    implicit none
    type(inventory), pointer :: i
    type(variable), pointer :: v
    integer :: ii
    character*20 :: myname = "ncf_printInventory"
    do ii=1,i%nrdim
       write(*,*)myname,ii,i%dim250(ii)(1:i%lend(ii))
    end do
    call ncf_printPos(i%pos)
    v=>i%firstvariable%next
    do while (.not.associated(v,target=i%lastVariable))
       call ncf_printVariable(v)
       v=>v%next
    end do
    return
  end subroutine ncf_printInventory

  subroutine ncf_printDimensions(i)
    implicit none
    type(inventory), pointer :: i
    integer :: ii
    character*20 :: myname = "ncf_printDimensions"
    do ii=1,i%nrdim
       write(*,*)myname," Dim:",ii,"'"//i%dim250(ii)(1:i%lend(ii))//"'"
    end do
  end subroutine ncf_printDimensions

  function ncf_gettype(type)
    implicit none
    character*6 ncf_gettype
    integer type
    character*20 :: myname = "ncf_gettype"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (type.eq.nf_char) then
       ncf_gettype="char  "
    elseif (type.eq.nf_int1) then
       ncf_gettype="int1  "
    elseif (type.eq.nf_int2) then
       ncf_gettype="int2  "
    elseif (type.eq.nf_int) then
       ncf_gettype="int   "
    elseif (type.eq.nf_real) then
       ncf_gettype="real  "
    elseif (type.eq.nf_double) then
       ncf_gettype="double"
    else
       ncf_gettype="any   "
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end function ncf_gettype

  subroutine ncf_printVariable(v)
    implicit none
    type(variable),pointer :: v
    type(dimension),pointer  :: id
    character*10 :: ctype10
    character*50 :: f50
    character*20 :: val20
    integer :: ii
    character*20 :: myname = "ncf_printVariable"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (.not.associated(v)) then
       write(*,*)myname,"***Attempt to print undefined variable***"
    else
       if (associated(v%f)) then
          f50=" "//v%f%lab10(1:v%f%lenl)//" '"//v%f%fn250(1:v%f%lenf)//"'"
       else
          f50=""
       end if
       if (associated(v%fd)) then
          write(val20,'(F20.5)') v%fd(1)
       else if (associated(v%fr)) then
          write(val20,'(F20.5)') v%fr(1)
       else
          val20="  --??--            "
       end if
       ctype10=ncf_gettype(v%type)
       write(*,*)myname,"Variable '"//v%var250(1:v%lenv)//"' ",trim(ctype10),&
            & associated(v%fc),associated(v%f1),associated(v%f2), &
            & associated(v%f4),associated(v%f8),associated(v%fr),&
            & associated(v%fd)," rf=",v%readyfield, &
            & trim(f50)
       write(*,'(3(X,A),10(2X,I0))')myname,"  ",&
            & val20,v%lenc,v%len1,v%len2,v%len4,v%len8,&
            & v%lenr,v%lend,v%len,v%lens,size(v%fd)
       do ii=1,v%nrdim
          if (v%ind(ii).gt.0.and.v%ind(ii).le.v%f%nrdim) then
             write(*,'(X,A," i=",I0,2X,A,2X,I0,2X,I0," <",I0,">")') myname,v%ind(ii), &
                  & "'"//v%f%dim250(v%ind(ii))(1:v%f%lend(v%ind(ii)))//"'",  &
                  & v%sta(ii),v%lim(ii),v%siz(ii)
          else
             write(*,*)myname,'Invalid dimension:',v%ind(ii),v%f%nrdim,&
                  & v%sta(ii),v%lim(ii)
          end if
       end do
    end if
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_printVariable

  subroutine ncf_printDimension(i,id)
    implicit none
    type(inventory),pointer  :: i
    type(dimension),pointer  :: id
    character*20 :: myname = "ncf_printDimension"
    if (associated(i)) then
       write(*,'(X,A," ind=",I3,X,A20,"->[",I0,",",I0,"]")') myname,id%ind, &
            & i%dim250(id%ind)(1:i%lend(id%ind)),id%sta,id%lim
    else
       write(*,'(X,A," ind=",I3,X,A20,"->[",I0,",",I0,"]")') myname,id%ind, &
            & "----",id%sta,id%lim
    end if
    return
  end subroutine ncf_printDimension
  
  subroutine ncf_printDimOrder(i)
    implicit none
    type(dimensionOrder) :: i
    type(dimension),pointer  :: id
    character*20 :: myname = "ncf_printDimOrder"
    integer :: cnt
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    write(*,*)myname,' Dimensions:',i%nrdim,i%resetPos
    cnt=0
    id=>i%firstDimension%next
    do while (.not.associated(id,i%lastDimension))
       call ncf_printDimension(i%f,id)
       id=>id%next
       cnt=cnt+1
       if (cnt.gt.i%nrdim) then
          write(*,*)myname,'Corrupt dimorder...',cnt,i%nrdim,&
               & associated(id,target=id%next)
          id=>i%lastDimension
       end if
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_printDimOrder

  subroutine ncf_printVarPos(v)
    implicit none
    type(variable), pointer  :: v
    integer :: ii
    character*20 :: myname = "ncf_printDimOrder"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    write(*,*)myname,' Position:',v%var250(1:v%lenv),v%nrdim
    do ii=1,v%nrdim
       write(*,'(X,A,X,I3,X,A20,3(X,I8))') myname,v%ind(ii), &
            & v%f%dim250(v%ind(ii))(1:v%f%lend(v%ind(ii)))//" = ",&
            & v%f%pos%pos(v%ind(ii)),v%sta(ii),v%lim(ii)
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
  end subroutine ncf_printVarPos

  subroutine ncf_printWeight(wgt)
    implicit none
    type(weight) :: wgt
    integer ii,jj
    character*20 :: myname = "ncf_printWeight"
    if (ncf_bdeb) write(*,*)myname,' Entering.'
    do jj=1,wgt%nweight
       write(*,*) 'ncf_printWeight Weight:',jj,wgt%w(jj)
       do ii=1,wgt%pos(jj)%nrdim
          write(*,*) 'ncf_printWeight      Pos:',ii,wgt%pos(jj)%pos(ii)
       end do
    end do
    if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end subroutine ncf_printWeight

  subroutine ncf_printPos(pos)
    implicit none
    type(position) :: pos
    integer ii
    character*250 :: buff250
    character*10 :: n10,s10,e10,a10
    integer :: lenb,lenn,lens,lene,lena
    integer, external :: length
    character*20 :: myname = "ncf_printPos"
    !if (ncf_bdeb) write(*,*)myname,' Entering.'
    if (pos%nrdim .gt.0) then
       buff250=myname
       call chop0(buff250,250)
       lenb=length(buff250,250,10)
       do ii=1,pos%nrdim
          write(n10,'(I10)') pos%pos(ii)
          write(s10,'(I10)') pos%sta(ii)
          write(e10,'(I10)') pos%lim(ii)
          write(a10,'(I10)') pos%siz(ii)
          call chop0(n10,10)
          call chop0(s10,10)
          call chop0(e10,10)
          call chop0(a10,10)
          lenn=length(n10,10,10)
          lens=length(s10,10,10)
          lene=length(e10,10,10)
          lena=length(a10,10,10)
          buff250=buff250(1:lenb)//" "//n10(1:lenn)//" ("//s10(1:lens)&
               & //"-"//e10(1:lene)//") <"//a10(1:lena)//">"
          call chop0(buff250,250)
          lenb=length(buff250,250,10)
       end do
       write(*,*)buff250(1:lenb)
    else
       write(*,*)'ncf_printPos: Position has no dimensions!'
    end if
    !if (ncf_bdeb) write(*,*)myname,' Done.'
    return
  end subroutine ncf_printPos

end module ncf

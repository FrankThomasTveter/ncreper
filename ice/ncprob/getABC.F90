SUBROUTINE GETABC (NABC,XABC,ne,e,bok,IRC)
  !
  IMPLICIT NONE
  SAVE
  !
  INTEGER NABC ! number of abcstant elements
  REAL XABC(NABC) ! abcstant elements
  integer :: ne ! number of moments
  real :: e(ne) ! moments
  logical bok   ! all ok?
  INTEGER IRC   ! fatal error if (irc.ne.0)...
  !
  INTEGER INFORM, ITER, N,NBND, &
       &     NCLIN,NCNLN,LIWORK,LWORK,NROWA, &
       &     NROWJ,NROWR,MAXBND 
  DOUBLE PRECISION   OBJF
  PARAMETER (&
       &     N = 4,&
       &     NCLIN  =  0,&
       &     NCNLN  =  1,&
       &     NROWA  =  NCLIN, &
       &     NROWJ  =  1, &
       &     NROWR =   N,&
       &     LIWORK =  3*n+NCLIN+2*NCNLN,&
       &     LWORK  =  2*N*N + N*NCLIN + 2*N*NCNLN + 20*N + 11*NCLIN + 21*NCNLN,&
       &     MAXBND =  N + NROWA + NROWJ)
  INTEGER            ISTATE(MAXBND)
  INTEGER            IWORK(LIWORK)
  DOUBLE PRECISION   A(NROWA,N)
  DOUBLE PRECISION   BL(MAXBND) 
  DOUBLE PRECISION   BU(MAXBND)
  DOUBLE PRECISION   C(NROWJ) 
  DOUBLE PRECISION   CJAC(NROWJ,N) 
  DOUBLE PRECISION   CLAMDA(MAXBND)
  DOUBLE PRECISION   OBJGRD(N) 
  DOUBLE PRECISION   R(NROWR,N) 
  DOUBLE PRECISION   X(N)
  DOUBLE PRECISION   WORK(LWORK)
  EXTERNAL           OBJFN, CONFN
  !
  INTEGER II,JJ,XPOS
  !
  CHARACTER*16 MYNAME
  DATA MYNAME /'GETABC'/
  !
  !     MINIMISE
  !
  CALL NPOPTN( 'NOLIST')
  CALL NPOPTN( 'Major print level              0')
  CALL NPOPTN( 'Minor print level              0')
  !CALL NPOPTN( '   Derivative level               0')
  !CALL NPOPTN( '   Step Level                   1.0')
  CALL NPOPTN( '   Verify                       Yes')
  ! CALL NPOPTN( '   Optimal Tolerance      0.0000001')
  ! CALL NPOPTN( '   Warm Start')
  ! CALL NPOPTN( '   Cold Start')
  ! CALL NPOPTN( '   Major iteration limit        100')
  !
  NBND   = N + NCLIN + NCNLN
  A(:,:)=0.0D0
  !
  BL(1)=-1.0D0
  BL(2)=0.5D0
  BL(3)=-1.0D0
  BL(4)=1.0D-3
  !
  BU(1)=1.0D0
  BU(2)=1.5D0
  BU(3)=1.0D0
  BU(4)=1.0D0
  !
  BU(n+1)=10.0D0 ! nonlinear limits...
  BL(n+1)=0.5D0
  ! first guess is: V = e(1) + U*sqrt(e(2)-e(1)*e(1))
  if (xabc(2).gt.0.1D0.and.xabc(2).lt.10.0D0) then ! use provided first guess
     do ii=1,n
        x(ii)=xabc(ii)
     end do
  else
     x(:)=1.0D-3
     X(2)=1.0D0
  end if
  !
  CALL IOBJFN(e,irc)
  if (irc.ne.0) then
     write(*,*) myname,'Error return from IOBJFN.',irc
     return
  end if
  CALL NPSOL ( N, NCLIN, NCNLN, NROWA, NROWJ, NROWR,&
       &           A, BL, BU,&
       &           CONFN, OBJFN,&
       &           INFORM, ITER, ISTATE,&
       &           C, CJAC, CLAMDA, OBJF, OBJGRD, R, X,&
       &           IWORK, LIWORK, WORK, LWORK )
  ! inform 1: no further improvement possible
  ! inform 2: no feasible point within linear constraints
  ! inform 3: no feasible point within non-linear constraints
  ! inform 4: reached iteration limit (Major iteration limit)
  ! inform 6: too small tolerance?
  ! inform 7: incorrect derivatives
  ! inform 9: invalid input parameter
  ! overflow: increase "Linear Feasibility Tolerance" or "Nonlinear Feasibility Tolerance"
  IF (INFORM .GT. 1) THEN ! ignore iteration limit errors  
     if (inform.ne.4.and.inform.ne.6) then
        !IRC=INFORM
        WRITE(*,*) MYNAME,'------------------------------'
        WRITE(*,*) MYNAME,'Error return from NPSOL.',inform, x
        write(*,*) ' X   =',(X(ii),ii=1,n)
        write(*,*) ' Obj =',objf
        write(*,*) ' ObjG=',(objgrd(ii),ii=1,n)
        CALL NPOPTN( '   Major print level             10')
        CALL NPOPTN( '   Minor print level              5')
        CALL NPOPTN( '   Verify                       Yes')
        !
        x(:)=1.0D-3
        X(2)=1.0D0
        !
        write(*,*) ' nabc =',NABC
        write(*,*) ' ne   =',NE
        do ii=1,ne
           write(*,'(X,A,I0,A,F20.10)') &
                & 'E(',II,')=',E(ii)
        end do
        !
        CALL NPSOL ( N, NCLIN, NCNLN, NROWA, NROWJ, NROWR,&
             &           A, BL, BU,&
             &           CONFN, OBJFN,&
             &           INFORM, ITER, ISTATE,&
             &           C, CJAC, CLAMDA, OBJF, OBJGRD, R, X,&
             &           IWORK, LIWORK, WORK, LWORK )
        WRITE(*,*) MYNAME,'Second return from NPSOL.',inform
        write(*,*) ' X   =',(X(ii),ii=1,n)
        write(*,*) ' Obj =',objf
        write(*,*) ' ObjG=',(objgrd(ii),ii=1,n)
        do ii=1,min(n,nabc)
           xabc(ii)=x(ii)
        end do
     else
     do ii=1,min(n,nabc)
        xabc(ii)=0.0D0
     end do
     end if
     bok=.false.
  else
     do ii=1,min(n,nabc)
        xabc(ii)=x(ii)
     end do
     !     
     bok=.true.
  END IF
  
  RETURN
END SUBROUTINE GETABC


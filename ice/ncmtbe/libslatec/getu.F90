real function getu(a,b,c,d,x,bok,irc)
  ! solves: val = a + b*u  + c*u**2 + d*u**3
  implicit none
  real a,b,c,d,x
  logical :: bok
  integer :: irc
  !
  real aa,bb,cc,dd,val
  common /cmb_u /aa,bb,cc,dd,val
  ! initialise common blocks
  integer, parameter :: n=1, ldfjac=n,LR=(N*(N+1))/2
  INTEGER :: IOPT,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,NJEV
  REAL :: XTOL,EPSFCN,FACTOR
  REAL :: U(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N),&
       & WA1(N),WA2(N),WA3(N),WA4(N)
  character*6, parameter :: myname="getu"
  !
  ! set common block
  aa=a
  bb=b
  cc=c
  dd=d
  val=x
  !  
  IOPT=1        ! we have a Jacobian
  u(1)=x        ! first guess
  XTOL=1.0D-3
  maxfev=25
  ML=1
  MU=1
  EPSFCN=1.0D-4
  MODE=1
  DIAG(1)=1.0D0
  FACTOR=100.0D0
  NPRINT=0
  INFO=0
  NJEV=10
  
  call SNSQ (FCN, JAC, IOPT, N, U, FVEC, FJAC, LDFJAC, XTOL, &
       &   MAXFEV, ML, MU, EPSFCN, DIAG, MODE, FACTOR, NPRINT, INFO, NFEV, &
       &   NJEV, R, LR, QTF, WA1, WA2, WA3, WA4)
  if (info.ne.1) then
     bok=.false.
     irc=info
     write(*,*)myname,'Error return from SNSQ.',info,x,u(1)
     return
  end if
  getu=u(1)
  return
contains 
  !
  SUBROUTINE FCN(N,X,FVEC,IFLAG)
    implicit none
    INTEGER N,IFLAG
    REAL X(N),FVEC(N)
    real aa,bb,cc,dd,val
    common /cmb_u /aa,bb,cc,dd,val
    fvec(1)= aa + X(1)*(bb + X(1)*(cc + X(1)*dd)) - val
    !write(*,*) 'FCN:',x(1),fvec(1)

    RETURN
  END SUBROUTINE FCN
  !
  SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
    implicit none
    INTEGER :: N,LDFJAC,IFLAG
    REAL :: X(N),FVEC(N),FJAC(LDFJAC,N)
    real :: aa,bb,cc,dd,val
    common /cmb_u /aa,bb,cc,dd,val
    fjac(1,1)= bb + 2.0D0*cc*x(1)+ 3.0D0*dd*x(1)*x(1)
    !write(*,*) 'JAC:',x(1),fvec(1),fjac(1,1)
    RETURN
  END SUBROUTINE JAC
  !
end function getu

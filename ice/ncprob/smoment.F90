  subroutine smoment(n,sxy,f,aa,bb,gm,bok)
    ! shift moments by mu/gm, i.e. f(1)=(e(1)+mu)*gm, f(2)=(e(2)+2*mu*e(1)+mu**2)*gm**2
    implicit none
    ! shifts moments from   u -> (u+mu)*gm
    integer :: n
    real :: sxy(0:n,0:n)
    real :: f(n)
    real :: aa
    real :: bb
    real :: gm
    logical :: bok
    real*16 :: gm2,gm3,gm4,gm5,gm6,g(n)
    real*16 :: aa2,aa3,aa4,bb2,bb3,bb4,buff
    integer :: ii,jj
    !
    aa2=aa*aa
    aa3=aa2*aa
    aa4=aa3*aa
    bb2=bb*bb
    bb3=bb2*bb
    bb4=bb3*bb
    g(1) = -aa - bb*sxy(1,0) + sxy(0,1)
    if (n.ge.2) g(2)= -2.0D0*aa*sxy(0,1) &
         & + aa2 &
         & - 2.0D0*bb*sxy(1,1) &
         & + 2.0D0*aa*bb*sxy(1,0) &
         & + bb2 * sxy(2,0) &
         & + sxy(0,2)
    if (n.ge.3) g(3)=-3.0D0*aa*sxy(0,2) &
         & + 3.0D0*aa2*sxy(0,1) &
         & - aa3 &
         & - 3.0D0*bb*sxy(1,2) &
         & + 6.0D0*aa*bb*sxy(1,1) &
         & - 3.0D0*aa2*bb*sxy(1,0) &
         & + 3.0D0*bb2*sxy(2,1) &
         & - 3.0D0*aa*bb2*SXY(2,0) &
         & - bb3*sxy(3,0) &
         & + sxy(0,3)
    if (n.ge.4) g(4)=-4.0D0*aa*sxy(0,3) &
         & + 6.0D0*aa2*sxy(0,2) &
         & - 4.0D0*aa3*sxy(0,1) + aa4 &
         & - 4.0D0*bb*sxy(1,3) &
         & + 12.0D0*aa*bb*sxy(1,2) &
         & - 12.0D0*aa2*bb*sxy(1,1) &
         & + 4.0D0*aa3*bb*sxy(1,0) &
         & + 6.0D0*bb2*sxy(2,2) &
         & - 12.0D0*aa*bb2*sxy(2,1) &
         & + 6.0D0*aa2*bb2*sxy(2,0) &
         & - 4.0D0*bb3*sxy(3,1) &
         & + 4.0D0*aa*bb3*sxy(3,0) &
         & + bb4*sxy(4,0) &
         & + sxy(0,4)
    if (g(2).lt.0.0D0.or.g(4).lt.0.0D0) then !sanity check
       !write(*,*)'AB=',aa,bb
       !do ii=0,n
       !   write(*,'(X,A,I0,A,5(X,F15.3))')'SXY(',ii,',*)=',(sxy(ii,jj),jj=0,n)
       !end do
       !do ii=1,n
       !   write(*,'(X,A,I0,A,5(X,F15.3))')'E(',ii,')=',g(ii)
       !end do
       !write(*,*)'SMOMENT Negative moments (2&4):',g(2),g(4)
       bok=.false.
    end if
    gm=g(2)-g(1)*g(1)
    if (gm.lt.1.0D-10) then
       write(*,*)'SMOMENT: Invalid Stdv, ', gm,f
       do ii=1,n
          write(*,'(100(X,F15.3))')(sxy(ii,jj),jj=1,n)
       end do
       bok=.false.
    else
       bok=.true.
       gm=1.0D0/sqrt(gm)
       gm2=gm*gm
       gm3=gm2*gm
       gm4=gm3*gm
       f(1)=g(1)*gm
       if (n.ge.2) f(2)=g(2)*gm2
       if (n.ge.3) f(3)=g(3)*gm3
       if (n.ge.4) f(4)=g(4)*gm4
    end if
    return
  end subroutine smoment
  !

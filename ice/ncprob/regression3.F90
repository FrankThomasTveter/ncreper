  subroutine regression3(xtx,xty,xm,cmean,order)
    implicit none
    real :: xtx(3,3) ! xtx(1,1)=1,xtx(2,2)=sum(X*X),xtx(3,3)=sum(x*x*x*x),xtx(1,2)=xtx(2,1)=sum(x), xtx(3,1)=xtx(1,3)=sum(x*x)
    real :: xty(3) ! xty(1)=sum(y), xty(2)=sum(x*y), xty(3)=sum(x*x*y)
    real :: xm    ! x-value
    real :: cmean ! mean at x-value
    integer :: order ! required order -> used order
    !
    real :: det, detinv, B(3,3),A(3)
    det = XTX(1,1)*XTX(2,2)*XTX(3,3) - XTX(1,1)*XTX(2,3)*XTX(3,2)&
         - XTX(1,2)*XTX(2,1)*XTX(3,3) + XTX(1,2)*XTX(2,3)*XTX(3,1)&
         + XTX(1,3)*XTX(2,1)*XTX(3,2) - XTX(1,3)*XTX(2,2)*XTX(3,1)
    if (order.ge.3.and.abs(det).gt.1.0D-10) then
       detinv=1.0D0/det
       ! Calculate the inverse of the matrix
       B(1,1) = +detinv * (XTX(2,2)*XTX(3,3) - XTX(2,3)*XTX(3,2))
       B(2,1) = -detinv * (XTX(2,1)*XTX(3,3) - XTX(2,3)*XTX(3,1))
       B(3,1) = +detinv * (XTX(2,1)*XTX(3,2) - XTX(2,2)*XTX(3,1))
       B(1,2) = -detinv * (XTX(1,2)*XTX(3,3) - XTX(1,3)*XTX(3,2))
       B(2,2) = +detinv * (XTX(1,1)*XTX(3,3) - XTX(1,3)*XTX(3,1))
       B(3,2) = -detinv * (XTX(1,1)*XTX(3,2) - XTX(1,2)*XTX(3,1))
       B(1,3) = +detinv * (XTX(1,2)*XTX(2,3) - XTX(1,3)*XTX(2,2))
       B(2,3) = -detinv * (XTX(1,1)*XTX(2,3) - XTX(1,3)*XTX(2,1))
       B(3,3) = +detinv * (XTX(1,1)*XTX(2,2) - XTX(1,2)*XTX(2,1))
       A(1)=B(1,1)*XTY(1) + B(1,2)*XTY(2) + B(1,3)*XTY(3)
       A(2)=B(2,1)*XTY(1) + B(2,2)*XTY(2) + B(2,3)*XTY(3)
       A(3)=B(3,1)*XTY(1) + B(3,2)*XTY(2) + B(3,3)*XTY(3)
       order=3
       cmean=a(1) + a(2)*xm + a(3)*xm*xm
    else
       det = XTX(1,1)*XTX(2,2) - XTX(1,2)*XTX(2,1)
       if (order.ge.2.and.abs(det).gt.1.0D-10) then
          detinv=1.0D0/det
          B(1,1) = +detinv * XTX(2,2)
          B(2,1) = -detinv * XTX(2,1)
          B(1,2) = -detinv * XTX(1,2)
          B(2,2) = +detinv * XTX(1,1)
          A(1)=B(1,1)*XTY(1) + B(1,2)*XTY(2)
          A(2)=B(2,1)*XTY(1) + B(2,2)*XTY(2)
          A(3)=B(3,1)*XTY(1) + B(3,2)*XTY(2)
          order=2
          cmean=a(1) + a(2)*xm
       else
          order=1
          cmean=xty(1)
       end if
    end if
    return
  end subroutine regression3

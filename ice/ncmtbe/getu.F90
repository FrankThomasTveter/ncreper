real function getU(a,b,c,d,x,bok,irc)
  implicit none
  real :: a,b,c,d,x
  logical :: bok
  integer :: irc
  real :: c2,c3,c4,c5
  real :: d2,d3,d4,d5
  complex :: z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,i
  complex :: u1,u2,u3
  i=(0.0D0,1.0D0)
  c2=c*c
  c3=c2*c
  c4=c3*c
  d2=d*d
  d3=d2*d
  d4=d3*d
  if (d.lt.1.0D-10) then
     if (b.gt.1.0D-10) then
        getu=(x-a)/b
        return
     else
        bok=.false.
        irc=348
        write(*,*)'GetU invalid transformation.',a,b,c,d,x,&
             & b-c*c/max(1.0D-10,3.0D0*d)-0.1D0
        return
     end if
  end if
  z1=(-1.0D0-I*Sqrt(3.0D0))/2.0D0
  z2=(-1.0D0+I*Sqrt(3.0D0))/2.0D0
  z3=-9.0D0*b*c/(2.0D0*d2)+c3/d3+27.0D0*(a-x)/(2.0D0*d)
  z4=(-9.0D0*b*c/d2+2.0D0*c3/d3+27.0D0*(a-x)/d)
  z5=(-3.0D0*b/d+c2/d2)
  z6=zSqrt(z4*z4-4.0D0*z5*z5*z5)
  z7=(z3+z6/2.0D0)**(1.0D0/3.0D0)
  if (abs(z7).lt.1.0D-10) then
     write(*,*) 'GetU z1:',z1
     write(*,*) 'GetU z2:',z2
     write(*,*) 'GetU z3:',z3
     write(*,*) 'GetU z4:',z4
     write(*,*) 'GetU z5:',z5
     write(*,*) 'GetU z6:',z6
     write(*,*) 'GetU z7:',z7
     write(*,*) 'GetU Invalid transformation.',a,b,c,d,x,&
             & b-c2/max(1.0D-10,3.0D0*d)-0.1D0
     bok=.false.
     irc=448
     return
  end if
  z8=c/(3.0D0*d)
  z9=z1*z7
  z10=z2*z7
  u1=-z8-z5/(3*z7)-z7/3
  u2=-z8-z9/3-z5/(3*z9)
  u3=-z8-z10/3-z5/(3*z10)
  if (abs(imagpart(u1)).lt.1.0D-10) then
     getU=realpart(u1)
  elseif (abs(imagpart(u2)).lt.1.0D-10) then
     getU=realpart(u2)
  elseif (abs(imagpart(u3)).lt.1.0D-10) then
     getU=realpart(u3)
  else
     write(*,*)'GetU Transformation:',a,b,c,d,x,&
             & b-c*c/max(1.0D-10,3.0D0*d)-0.1D0
     write(*,*)'GetU Products:',c2,c3,c4,d2,d3,d4
     write(*,*) 'GetU z1:',z1
     write(*,*) 'GetU z2:',z2
     write(*,*) 'GetU z3:',z3
     write(*,*) 'GetU z4:',z4
     write(*,*) 'GetU z5:',z5
     write(*,*) 'GetU z6:',z6
     write(*,*) 'GetU z7:',z7
     write(*,*) 'GetU z8:',z8
     write(*,*) 'GetU z9:',z9
     write(*,*) 'GetU z10:',z10
     write(*,*)'GetU No solution found.',u1,u2,u3
     bok=.false.
     irc=548
  end if
  bok=.true.
  return
end function getU

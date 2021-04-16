program ttc
  implicit none
  integer iunit
  integer,parameter :: ne=4
  integer :: ii,jj,kk,mm,irc
  integer, parameter :: nabc=4
  real :: ff(ne),sxy(0:ne,0:ne),xabc(nabc),mu,gm,a,b
  real, external :: getu
  logical :: bdone,bok
  bok=.true.
  sxy(           0 ,           0 )=   1.0000000000000000     
  sxy(           0 ,           1 )=   2.4894873155669000E-002
  sxy(           0 ,           2 )=   1.8572365505036363E-003
  sxy(           0 ,           3 )=   3.3188937072259106E-004
  sxy(           0 ,           4 )=   7.6711578916613410E-005
  sxy(           1 ,           0 )=   9.0000000000000000     
  sxy(           1 ,           1 )=  0.22405385840102099     
  sxy(           1 ,           2 )=   1.6715128954532723E-002
  sxy(           1 ,           3 )=   2.9870043365033197E-003
  sxy(           2 ,           0 )=   81.000000000000000     
  sxy(           2 ,           1 )=   2.0164847256091889     
  sxy(           2 ,           2 )=  0.15043616059079451     
  sxy(           3 ,           0 )=   729.00000000000000     
  sxy(           3 ,           1 )=   18.148362530482700     
  sxy(           4 ,           0 )=   6561.0000000000000     
  a=   4.3348887733604652E-002
  b=  -1.3675576273356384E-003
  ff(           1 )=  -3.4032420042174958E-015
  ff(           2 )=   490.45886150545709     
  ff(           3 )=   31266.380863243805     
  ff(           4 )=   3359649.0311528305     
  mu=   3.1040869087583904E-002
  gm=   490.45886150545709     
  call smoment(ne,sxy,ff,a,b,gm,bok)
  do ii=1,ne
     write(*,*) 'ff(',ii,')=',ff(ii)
  end do
  
#ifdef JNK   
   bok=.true.
   irc=0
   ff(1)=0.0D0
   ff(2)=1.0D0
   ff(3)=0.0D0
   ff(4)=3.0D0
   call getABC(NABC,XABC,ne,ff,bok,IRC)
   write(*,'(X,A,X,L1,X,I0,10(X,F15.8))')'Return with:',bok,irc,&
        & (xabc(mm),mm=1,nabc)
#endif
 end program ttc
 

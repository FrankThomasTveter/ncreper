program ttt
  implicit none
  integer iunit
  integer,parameter :: mpar=3,mtime=41,mpos=1000,ne=4
  real :: aval(mpar,mtime,mpos),aprob(mpar,mtime,mpos),aacc(mpar,mtime,mpos)
  integer :: npar,ii,jj,kk,mm,irc
  integer :: ipar,itime,ipos, opar,otime
  integer :: ctime,cpos
  integer :: ntime(mpar),npos(mpar,mtime)
  integer, parameter :: nabc=4
  real :: ee(ne),ff(ne),gg(ne),val,prob,acc,xval,xabc(nabc),mu,gm,dxdu,sprob
  real, external :: getu
  logical :: bdone,bok
  
  open(50,file='fort.10')
  opar=0
  otime=0
  npar=0
  do ii=1,mpar
     ntime(ii)=0
     do jj=1,mtime
        npos(ii,jj)=0
     end do
  end do
  bdone=.false.
  do while(.not.bdone)
     read(50,*,iostat=irc) ipar,itime,ipos,val,prob,acc
     if (irc.ne.0) then
        bdone=.true.
     else
        if (ipar.ne.opar) then
           npar=npar+1
           ctime=0
           opar=ipar
        end if
        if (itime.ne.otime) then
           ctime=ctime+1
           otime=itime
           cpos=0
        end if
        cpos=min(mpos,cpos+1)
        ntime(npar)=ctime ! ntime
        npos(npar,ntime(npar))=min(mpos,npos(npar,ntime(npar))+1)
        !write(*,*)'Debug:',npar,ntime(npar),npos(npar,ntime(npar))
        aval(npar,ntime(npar),npos(npar,ntime(npar)))=val+100.0D0
        aprob(npar,ntime(npar),npos(npar,ntime(npar)))=prob
        aacc(npar,ntime(npar),npos(npar,ntime(npar)))=acc
     end if
  end do
  close (50)
  ! get moments
  do ii=1,npar
     do jj=1,ntime(ii)
        do mm=1,ne
           ee(mm)=0.0D0
        end do
        sprob=(aval(ii,jj,npos(ii,jj))-aval(ii,jj,1))/real(npos(ii,jj))
        do kk=1,npos(ii,jj)
           val=aval(ii,jj,kk)
           prob=aprob(ii,jj,kk)
           do mm=1,ne
              ee(mm)=ee(mm) + (val**mm)*prob
           end do
        end do
        write(*,*)'=====> SPROB:',sprob
        mu=ee(1)
        gm=1.0D0/Sqrt(max(0.0D0,ee(2)-ee(1)*ee(1)))
        write(*,*)'Calling Trans.',mu,gm,(ee(mm),mm=1,ne)
        call smoment(ne,ee,ff,mu,gm)
        !stop ("Remember to transform so that ff(1)=0..., or correct cost function!!!")
        call getABC(NABC,XABC,ne,ff,bok,IRC)
        write(*,'(X,A,10(X,F15.8))')'Recalc with:',(xabc(mm),mm=1,nabc)
        write(52,'(2(X,I3),100(X,F27.15))') ii,jj,(xabc(mm),mm=1,nabc),(ff(mm),mm=1,ne)
        do mm=1,ne
           gg(mm)=0.0D0
        end do
        do kk=1,npos(ii,jj)
           val=(aval(ii,jj,kk)-mu)*gm
           prob=aprob(ii,jj,kk)

           write(60,'(3(X,I0),3(X,F27.10))',iostat=irc) &
                & ii,jj,kk,val,prob/(gm*sprob),&
                & exp(-val*val*0.5D0)/(Sqrt(2.0D0*3.141592654))

           !xval=mu+(getu(xabc(1),xabc(2),xabc(3),xabc(4),val)/gm)
           xval=getu(xabc(1),xabc(2),xabc(3),xabc(4),val)


           dxdu=xabc(2)+xabc(3)*2.0D0*xval+xabc(4)*3.0D0*xval*xval
           write(80,'(3(X,I0),3(X,F27.10))',iostat=irc) &
                & ii,jj,kk,xval,aprob(ii,jj,kk)*dxdu/(gm*sprob),&
                & exp(-xval*xval*0.5D0)/(Sqrt(2.0D0*3.141592654))

           write(81,'(3(X,I4),4(X,F27.15))')&
                & ii,jj,kk,val,xval,aprob(ii,jj,kk),aacc(ii,jj,kk)
            do mm=1,ne
              gg(mm)=gg(mm) + (xval**mm)*prob
           end do
        end do
        write(*,*)'Recalc done.',mu,gm,(ff(mm),mm=1,ne),(gg(mm),mm=1,ne)

        write(82,'(2(X,I3),100(X,F27.15))') ii,jj,(xabc(mm),mm=1,nabc),(ee(mm),mm=1,ne)

     end do
  end do
  
contains
  !
end program ttt

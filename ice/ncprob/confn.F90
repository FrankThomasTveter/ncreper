      SUBROUTINE CONFN(MODE, NCNLN, N, NROWJ,&
           &     NEEDC, X, C, CJAC, NSTATE )
        !
        ! INPUT
        ! mode: can be ignored if def der level is being used
        ! ncnln : dimension of c, same as in npsol
        ! n : number of control variables
        ! nrowj (>= 1 && >=ncnln), leading dimension of cJac
        ! needc : can be ignored by unsophisticated user
        ! x : control variables
        ! nstate : same as for funobj
        ! OUTPUT
        ! mode : error return code, -1 => npsol terminates
        ! c : array of dimension minimum ncnln contain values of constraints
        ! cJac : array of declared dimension (nrowj,k) k>=n, jacobians evaluated at x
        !
        IMPLICIT           none
        SAVE
        integer :: mode, ncnln, n, nrowj
        INTEGER :: NEEDC(*)
        DOUBLE PRECISION   X(N), C(*), CJAC(NROWJ,*)
        integer :: nstate
        real aa,bb,cc,dd
        !
        CHARACTER*16 MYNAME
        DATA MYNAME /'CONFN'/
        !
        ! We must constraint the solution to a strictly growing transformation...
        !
        mode=0.0D0
        bb=x(2)
        cc=x(3)
        dd=max(1.0D-10,x(4))
        c(1) = bb-cc*cc/(3.0D0*dd)
        CJAC(1,1)=0.0D0
        CJAC(1,2)=1.0D0
        CJAC(1,3)=-2.0D0*cc/(3.0D0*dd)
        CJAC(1,4)=cc*cc/(3.0D0*dd*dd)
        !     
        !     DUMMY ROUTINE
        !
        !WRITE(*,*) MYNAME,'This routine should never be called.'
        !
        RETURN
      END SUBROUTINE CONFN

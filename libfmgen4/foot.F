      SUBROUTINE FOOT (KODE,
     &     AA,VEC,YCEN,XCEN,
     &     YMIN,YMAX,Y,XMIN,XMAX,
     &     IRC)
C     ***********************************************************************
C     WRITTEN/MODIFIED BY:                                                  *
C     --------------------------------------------------------------------- *
C     |    NAME      |   DATE   |                 REASON                  | *
C     --------------------------------------------------------------------- *
C     | F. TVETER    | 20/08/02 | NEW                                     | *
C     |              |          |                                         | *
C     --------------------------------------------------------------------- *
C     ***********************************************************************
C
      IMPLICIT NONE
      SAVE
C
      INTEGER  KODE            ! -1=GET GLOBAL YMIN, YMAX; 0=GET XMIN, XMAX
      REAL   AA(2),VEC(2,2)     ! FOOTPRINT SPECIFICATIONS (INPUT)
      REAL   XCEN,YCEN          ! CENTER POSITION
      REAL   YMIN,YMAX          ! GLOBAL MIN/MAX VALUES FOR Y (OUTPUT)
      REAL   Y                  ! CURRENT Y VALUE (INPUT)
      REAL   XMIN,XMAX          ! RESULTING MIN/MAX X (OUTPUT)
      INTEGER  IRC
C
      CHARACTER*8 MYNAME
      DATA MYNAME /'FOOT'/
      LOGICAL  INITIALISED
      DATA INITIALISED /.FALSE./
      REAL   Z(50),DY,DX
C
      IF (KODE.EQ.-1) THEN
         INITIALISED=.TRUE.
         Z(1)=AA(2)
         Z(2)=Z(1)**2
         Z(3)=VEC(1,1)
         Z(4)=Z(3)**2
         Z(5)=-(Z(2)*Z(4))
         Z(6)=AA(1)
         Z(7)=Z(6)**2
         Z(8)=VEC(1,2)
         Z(9)=Z(8)**2
         Z(10)=-(Z(7)*Z(9))
         Z(11)=Z(5) + Z(10)
         Z(12)=SQRT(-Z(11))
         Z(13)=VEC(2,1)
         Z(14)=Z(13)**2
         Z(15)=-(Z(9)*Z(14))
         Z(16)=VEC(2,2)
         Z(17)=2*Z(3)*Z(8)*Z(13)*Z(16)
         Z(18)=Z(16)**2
         Z(19)=-(Z(4)*Z(18))
         Z(20)=Z(15) + Z(17) + Z(19)
         Z(21)=SQRT(max(0.,-Z(20)))
         Z(22)=1/Z(21)
         Z(23)=Z(2)*Z(4)
         Z(24)=Z(7)*Z(9)
         Z(25)=Z(23) + Z(24)
         Z(26)=1/Z(25)
         Z(29)=Z(1)**4
         Z(30)=Z(4)*Z(7)*Z(29)
         Z(31)=Z(6)**4
         Z(32)=Z(2)*Z(9)*Z(31)
C
C     GET YMIN AND YMAX
C
         YMIN=YCEN-(Z(12)*Z(22))
         YMAX=YCEN+Z(12)*Z(22)
C
      ELSE IF (KODE.EQ.0.AND.INITIALISED) THEN
C
C     GET XMIN AND XMAX
C
C
         DY=Y-YCEN
         Z(27)=-(DY*Z(2)*Z(3)*Z(13)*Z(26))
         Z(28)=-(DY*Z(7)*Z(8)*Z(16)*Z(26))
         Z(33)=DY*DY
         Z(34)=-(Z(2)*Z(7)*Z(9)*Z(14)*Z(33))
         Z(35)=2*Z(2)*Z(3)*Z(7)*Z(8)*Z(13)*Z(16)*Z(33)
         Z(36)=-(Z(2)*Z(4)*Z(7)*Z(18)*Z(33))
         Z(37)=Z(30) + Z(32) + Z(34) + Z(35) + Z(36)
         if (z(37).lt.0.) then
            XMIN=XCEN + Z(27) + Z(28)
            XMAX=XCEN + Z(27) + Z(28)
         else
            Z(38)=SQRT(Z(37))
            XMIN=XCEN + Z(27) + Z(28) - Z(26)*Z(38)
            XMAX=XCEN + Z(27) + Z(28) + Z(26)*Z(38)
         end if
C
      ELSE IF (KODE.EQ.1) THEN
         INITIALISED=.FALSE.
      ELSE
         WRITE(*,*) MYNAME,'Invalid KODE or call sequence.',KODE
         IRC=933
         RETURN
      END IF
C
      RETURN
      END

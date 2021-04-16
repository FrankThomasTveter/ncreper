      SUBROUTINE FOOTROT (XCM,YCM,XAM,YAM,XBM,YBM,AA,VEC,IRC)
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
C     CALCULATES ELIPSE-PARAMETERS IN LINEARLY-TRANSFORMED COORDINATE SYSTEM
C
      IMPLICIT NONE
      SAVE
C
      REAL   XCM,YCM            ! POSITION OF CENTER (IN TRANSFORMED SYSTEM)
      REAL   XAM,YAM            ! POSITION OF ELIPSE-SEMI-MAJOR AXIS A
      REAL   XBM,YBM            ! POSITION OF ELIPSE-SEMI-MAJOR AXIS B
      REAL   AA(2)              ! OUT:SEMI MAJOR AXIS IN TRANSFORMED SYSTEM
      REAL   VEC(2,2)           ! OUT:DIRECTION OF SEMI-MAJOR AXIS' -"-
      INTEGER  IRC
C
      CHARACTER*8 MYNAME
      DATA MYNAME /'FOOTROT'/
      REAL   Z(23)              ! WORKSPACE
      REAL   DXAM,DXBM,DYAM,DYBM,RR(2)
      INTEGER  II
C
      DXAM=XAM-XCM
      DYAM=YAM-YCM
      DXBM=XBM-XCM
      DYBM=YBM-YCM
C
      Z(1)=DSQRT(2.0D0)
      Z(2)=1.0D0/Z(1)
      Z(3)=DXAM*DXAM
      Z(4)=DXBM*DXBM
      Z(5)=DYAM*DYAM
      Z(6)=DYBM*DYBM
      Z(7)=-Z(3)
      Z(8)=-Z(4)
      Z(9)=-Z(5)
      Z(10)=-Z(6)
      Z(11)=Z(7) + Z(8) + Z(9) + Z(10)
      Z(12)=Z(11)*Z(11)
      Z(13)=Z(4)*Z(5)
      Z(14)=-2.0D0*DXAM*DXBM*DYAM*DYBM
      Z(15)=Z(3)*Z(6)
      Z(16)=Z(13) + Z(14) + Z(15)
      Z(17)=-4*Z(16)
      Z(18)=Z(12) + Z(17)
      Z(19)=SQRT(Z(18))
      Z(20)=-Z(19)
      Z(21)=DXAM*DYAM
      Z(22)=DXBM*DYBM
      Z(23)=Z(21) + Z(22)
      AA(1)=Z(2)*SQRT(Z(3) + Z(4) + Z(5) + Z(6) + Z(20))
      AA(2)=Z(2)*SQRT(Z(3) + Z(4) + Z(5) + Z(6) + Z(19))
C
      VEC(1,1)=(Z(3) + Z(4) + Z(9) + Z(10) + Z(20))/2.0D0
      VEC(2,1)=Z(23)
      IF (ABS(AA(1)-AA(2)).LT.1.0D-3) THEN ! CIRCLE
         VEC(1,2)=-VEC(2,1)
         VEC(2,2)=VEC(1,1)
      ELSE
         VEC(1,2)=(Z(3) + Z(4) + Z(9) + Z(10) + Z(19))/2.0D0
         VEC(2,2)=Z(23)
      END IF
C
      DO II=1,2
         RR(II)=SQRT(VEC(1,II)*VEC(1,II)+VEC(2,II)*VEC(2,II))
      END DO
C
C     NORMALISE VECTOR-LENGTH
C
      IF (RR(1).LT.1.0D-10 .AND. RR(2).LT.1.0D-10) THEN
         WRITE(*,'(X,A,A,3(" [",F10.3,X,F10.3,"]"))') 
     &        MYNAME,'Singular axis.',
     &        XCM,YCM,XAM,YAM,XBM,YBM
         write(*,*) 
         IRC=634
         RETURN
      ELSE IF (RR(1).LT.1.0D-10) THEN
         VEC(1,2)=VEC(1,2)/RR(2)
         VEC(2,2)=VEC(2,2)/RR(2)
         VEC(1,1)=-VEC(2,2)
         VEC(2,1)=VEC(1,2)
      ELSE IF (RR(2).LT.1.0D-10) THEN
         VEC(1,1)=VEC(1,1)/RR(1)
         VEC(2,1)=VEC(2,1)/RR(1)
         VEC(1,2)=-VEC(2,1)
         VEC(2,2)=VEC(1,1)
      ELSE
         VEC(1,2)=VEC(1,2)/RR(2)
         VEC(2,2)=VEC(2,2)/RR(2)
         VEC(1,1)=VEC(1,1)/RR(1)
         VEC(2,1)=VEC(2,1)/RR(1)
      END IF
C
      RETURN
      END


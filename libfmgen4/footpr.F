      subroutine footpr(NX,NY,BPRO,
     &     XCM,YCM,XAM,YAM,XBM,YBM,
     &     MAXPOS,NRPOS,POS,
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
C     find satellite observation footprint and
C     returns list of grid points within elipse
C
C
      IMPLICIT NONE
      SAVE
C
      INTEGER  NX              ! dimension of grid in x position
      INTEGER  NY              ! dimension of grid in y position
      LOGICAL  BPRO            ! invalid input parameters?
      REAL   XCM,YCM            ! grid-position of center of elipse
      REAL   XAM,YAM            ! grid-position of semi-major axis extreme in elipse
      REAL   XBM,YBM            ! grid-position of semi-minor axis extreme in elipse
      INTEGER  MAXPOS          ! max number of allowed positions returned
      INTEGER  NRPOS           ! number of positions returned
      INTEGER  POS(2,MAXPOS)   ! list of grid-positions (x,y) within elipse
      INTEGER  IRC             ! error return code
C
      CHARACTER*8 MYNAME
      DATA MYNAME/'FOOTPR'/
      INTEGER  IXMIN,IXMAX,IYMIN,IYMAX,KODE,IY,IX
      REAL   Y,AA(2),VEC(2,2),YMIN,YMAX,XMIN,XMAX
C
      NRPOS=0
C
C     FIND YMAX AND YMIN
C
      CALL FOOTROT (XCM,YCM,XAM,YAM,XBM,YBM,AA,VEC,IRC)
      IF (IRC.NE.0) THEN
         WRITE(*,*) MYNAME,'Error return from FOOTROT.',IRC
         RETURN
      END IF
C
      KODE=-1
      CALL FOOT (KODE,
     &     AA,VEC,YCM,XCM,
     &     YMIN,YMAX,Y,XMIN,XMAX,
     &     IRC)
      IF (IRC.NE.0) THEN
         WRITE(*,*) MYNAME,'Error return from FOOT.',IRC
         RETURN
      END IF
C
C     LOOP OVER POINTS WITHIN ELIPSE
C
      KODE=0
      IYMAX=NINT(YMAX+0.5D0)
      IYMIN=NINT(YMIN-0.5D0)
C
      DO IY=IYMIN,IYMAX
         Y=REAL(IY)
         CALL FOOT (KODE,
     &        AA,VEC,YCM,XCM,
     &        YMIN,YMAX,Y,XMIN,XMAX,
     &        IRC)
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from FOOT.',IRC
            RETURN
         END IF
         IXMAX=NINT(XMAX+0.5D0)
         IXMIN=NINT(XMIN-0.5D0)
         DO IX=IXMIN,IXMAX
            IF (IX.LE.NX .AND. IX.GE.1.AND.
     &           IY.LE.NY.AND.IY.GE.1) THEN
               NRPOS=NRPOS+1
               IF (NRPOS.LE.MAXPOS) THEN
                  POS(1,NRPOS)=IX
                  POS(2,NRPOS)=IY
               END IF
            END IF
         END DO
      END DO
      IF (NRPOS.LE.0) THEN
         BPRO=.FALSE.
         RETURN
      END IF
      IF (NRPOS.GT.MAXPOS) THEN
         WRITE(*,*) MYNAME,'Too much grid in footprint.',
     &        NRPOS,'(',MAXPOS,')'
         IRC=934
         RETURN
      END IF
C
      KODE=1
      CALL FOOT (KODE,
     &     AA,VEC,YCM,XCM,
     &     YMIN,YMAX,Y,XMIN,XMAX,
     &     IRC)
      IF (IRC.NE.0) THEN
         WRITE(*,*) MYNAME,'Error return from FOOT.',IRC
         RETURN
      END IF
C
      IF (NRPOS.EQ.0) THEN
         WRITE(*,*) MYNAME,'Too small footprint.'
         IRC=384
         RETURN
      END IF
C
      RETURN
      END

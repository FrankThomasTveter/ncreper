      SUBROUTINE PROINK(KODE,TARGLEV, TARGPRES, TARGARR,
     &     CLEV, CPRES, CARR,
     &     NRLEV, PRES, ARR,
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
      INTEGER  KODE,TARGLEV
      REAL   TARGPRES(TARGLEV),TARGARR(TARGLEV)
      INTEGER  CLEV
      REAL   CPRES(CLEV), CARR(CLEV)
      INTEGER  NRLEV
      REAL   PRES(NRLEV), ARR(NRLEV)
      INTEGER  IRC
C
      CHARACTER*8 MYNAME
      DATA MYNAME /'PROINK'/
C
      INTEGER  MAXIND,MAXLEV
      PARAMETER (MAXIND = 10,MAXLEV = 100)
      INTEGER  NRIND(MAXLEV)
      LOGICAL  QCLIM(MAXLEV,MAXIND)
      REAL   ALPHA(MAXLEV,MAXIND)
      INTEGER  INDX(MAXLEV,MAXIND)
C
      INTEGER  II,JJ
C
      LOGICAL  FIRST
      DATA FIRST /.TRUE./
C
      IF (FIRST) THEN
C
C     DO SOME SIMPLE CHECKING OF RANGE...
C
         IF (TARGLEV.GT.MAXLEV) THEN
            WRITE(*,*) MYNAME,'Too large TARGLEV.',TARGLEV,MAXLEV
            IRC=737
            RETURN
         END IF
         IF (TARGPRES(1).GT.TARGPRES(TARGLEV)) THEN
            WRITE(*,*) MYNAME,'Inverted TARG pressure!'
            IRC=494
            RETURN
         END IF
         IF (CPRES(1).GT.CPRES(CLEV)) THEN
            WRITE(*,*) MYNAME,'Inverted CLIM pressure!'
            IRC=495
            RETURN
         END IF
         IF (PRES(1).GT.PRES(NRLEV)) THEN
            WRITE(*,*) MYNAME,'Inverted PROFILE pressure!'
            IRC=496
            RETURN
         END IF
         IF (TARGPRES(1).LT.CPRES(1) .OR.
     &        TARGPRES(TARGLEV).GT.CPRES(CLEV)) THEN
            WRITE(*,*) MYNAME,'TARG outside CLIM range!'
            IRC=497
            RETURN
         END IF
C
         FIRST=.FALSE.
      END IF
C
C     MAKE INDEX ARRAYS
C
      IF (KODE.EQ.-1) THEN
         CALL INDXARR (MAXIND,
     &        MAXLEV,TARGLEV, TARGPRES, ! target levels
     &        CLEV,CPRES,       ! climatology
     &        NRLEV,PRES,       ! profile
     &        NRIND,QCLIM,ALPHA,INDX,
     &        IRC)
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from INDXARR.',IRC
            RETURN
         END IF
      END IF
C
C     INTERPOLATE
C
      ARR(:)=0.0D0
      CARR(:)=0.0D0
C
      DO II=1,TARGLEV
         DO JJ=1,NRIND(II)
            IF (QCLIM(II,JJ)) THEN ! use climatology
               CARR(INDX(II,JJ))=CARR(INDX(II,JJ))+
     &              TARGARR(II)*ALPHA(II,JJ)
            ELSE                ! use profile
               ARR(INDX(II,JJ))=ARR(INDX(II,JJ))+
     &              TARGARR(II)*ALPHA(II,JJ)
            END IF
         END DO
      END DO
C     
      RETURN
      END
      

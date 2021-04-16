      SUBROUTINE GETARR(NRLEV,PRES,ARR,
     &        TARGPRES,TARGARR,IRC)
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
C
      SAVE
C
      INTEGER  NRLEV
      REAL   PRES(NRLEV), ARR(NRLEV), DP
      REAL      TARGPRES, TARGARR
      INTEGER  IRC
C
      INTEGER  LASTPOS, CPOS
      LOGICAL  REPERR
      DATA REPERR /.FALSE./
      DATA LASTPOS /1/
      CHARACTER*8 MYNAME
      DATA MYNAME /'GETARR'/
C
      LASTPOS=MIN(LASTPOS,NRLEV)
C
      IF (TARGPRES-PRES(NRLEV).GT.1.0D-6) THEN
C     extrapolate...
         TARGARR=ARR(NRLEV) + 
     &        (TARGPRES-PRES(NRLEV))*
     &        (ARR(NRLEV)-ARR(NRLEV-1))/
     &        (PRES(NRLEV)-PRES(NRLEV-1))
         IF (.NOT.REPERR) THEN
            WRITE(*,*) MYNAME,'PRESSURE OUT OF RANGE:',
     &           TARGPRES
            REPERR=.TRUE.
         END IF
         RETURN
      END IF
C
      IF (TARGPRES-PRES(1).LT.-1.0D-6) THEN
         TARGARR=ARR(1) + 
     &        (TARGPRES-PRES(1))*
     &        (ARR(1)-ARR(1+1))/
     &        (PRES(1)-PRES(1+1))
         IF (.NOT.REPERR) THEN
            WRITE(*,*) MYNAME,'PRESSURE OUT OF RANGE:',
     &           TARGPRES
            REPERR=.TRUE.
         END IF
         RETURN
      END IF
C
      CPOS=LASTPOS
C
      DO WHILE(TARGPRES.GT.PRES(CPOS))
         CPOS=CPOS+1
      END DO
C
      DO WHILE(TARGPRES.LT.PRES(CPOS))
         CPOS=CPOS-1
      END DO
C
      DP=LOG(PRES(CPOS+1)/PRES(CPOS))
      TARGARR=ARR(CPOS)+
     &     LOG(TARGPRES/PRES(CPOS))*
     &     (ARR(CPOS+1)-ARR(CPOS))/DP
C
      LASTPOS=CPOS
C
      RETURN
      END

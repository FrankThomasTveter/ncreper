      SUBROUTINE CLOSEST(NROBS,LONA,LATA,
     &     INDLON,INDLAT,IINDLON,IINDLAT,
     &     ID1A,C1,TPOS,EPOS,MINDIFF,IRC)
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
      INTEGER  NROBS
      REAL   LATA(NROBS), LONA(NROBS)
      INTEGER  INDLAT(NROBS), INDLON(NROBS), 
     &     IINDLAT(NROBS), IINDLON(NROBS)
      CHARACTER*1 ID1A(NROBS),C1
      INTEGER  TPOS, EPOS, IRC
C
      CHARACTER*10 MYNAME
      DATA MYNAME /'CLOSEST'/
      LOGICAL  BDONE
      INTEGER  CPOS
      REAL   TARGLON,TARGLAT,MINDIFF,DIFF
C
      TARGLON=LONA(TPOS)
      TARGLAT=LATA(TPOS)
C
      MINDIFF=180.0D0
C
C     GO FORWARD IN LATITUDE, CONTINUOUSLY RESETTING MINDIFF
C      
      EPOS=0
      CPOS=IINDLAT(TPOS)+1
      BDONE=(CPOS.GT.NROBS)
      DO WHILE (.NOT. BDONE)
         IF (C1.EQ.ID1A(INDLAT(CPOS))) THEN
            IF (ABS(TARGLAT-LATA(INDLAT(CPOS))).GT.MINDIFF) THEN
               BDONE=.TRUE.
            ELSE
               CALL DIST(TARGLON,TARGLAT,
     &              LONA(INDLAT(CPOS)),LATA(INDLAT(CPOS)),
     &              DIFF,IRC)
               IF (IRC.NE.0) THEN
                  WRITE(*,*) MYNAME,'Error return from DIST.',IRC
                  RETURN
               END IF
               IF (DIFF.LT.MINDIFF) THEN
                  EPOS=INDLAT(CPOS)
                  MINDIFF=DIFF
               END IF
               CPOS=CPOS+1
               BDONE=(CPOS.GT.NROBS)
            END IF
         ELSE
            CPOS=CPOS+1
            BDONE=(CPOS.GT.NROBS)
         END IF
      END DO
C
C     GO BACKWARDS IN LATITUDE, CONTINUOUSLY RESETTING MINDIFF
C
      CPOS=IINDLAT(TPOS)-1
      BDONE=(CPOS.LT.1)
      DO WHILE (.NOT. BDONE)
         IF (C1.EQ.ID1A(INDLAT(CPOS))) THEN
            IF (ABS(TARGLAT-LATA(INDLAT(CPOS))).GT.MINDIFF) THEN
               BDONE=.TRUE.
            ELSE
               CALL DIST(TARGLON,TARGLAT,
     &              LONA(INDLAT(CPOS)),LATA(INDLAT(CPOS)),
     &              DIFF,IRC)
               IF (IRC.NE.0) THEN
                  WRITE(*,*) MYNAME,'Error return from DIST.',IRC
                  RETURN
               END IF
               IF (DIFF.LT.MINDIFF) THEN
                  EPOS=INDLAT(CPOS)
                  MINDIFF=DIFF
               END IF
               CPOS=CPOS-1
               BDONE=(CPOS.LT.1)
            END IF
         ELSE
            CPOS=CPOS-1
            BDONE=(CPOS.LT.1)
         END IF
      END DO
C
C      WRITE(*,*) MYNAME,'Routine not implemented.'
C      IRC=999
C
      RETURN
      END


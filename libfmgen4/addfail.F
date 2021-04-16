      SUBROUTINE ADDFAIL(ERB250,ERR250,LENE,ECNT,NRERR,MAXERR,IRC)
C     ***********************************************************************
C     WRITTEN/MODIFIED BY:                                                  *
C     --------------------------------------------------------------------- *
C     |    NAME      |   DATE   |                 REASON                  | *
C     --------------------------------------------------------------------- *
C     | F. TVETER    | 18/07/95 | NEW                                     | *
C     |              |          |                                         | *
C     --------------------------------------------------------------------- *
C     ***********************************************************************
C     
      IMPLICIT NONE
C     
      INTEGER  MAXERR
      CHARACTER*250 ERB250, ERR250(MAXERR)
      INTEGER  LENE(MAXERR)
      INTEGER  ECNT(MAXERR)
      INTEGER  NRERR
      INTEGER  IRC
C     
      INTEGER  LENB, ECNTB, LENEB
      CHARACTER*250 ERR250B
      LOGICAL  FOUND,BDONE
      INTEGER  II, LENGTH
      EXTERNAL LENGTH
      CHARACTER*8 MYNAME
      DATA MYNAME /'ADDFAIL'/
C     
      CALL CHOP(ERB250,250)
      LENB=LENGTH(ERB250,250,40)
      FOUND = .FALSE.
      II=1
      BDONE=(FOUND.OR.II.GT.NRERR)
      DO WHILE (.NOT.BDONE)
         IF (ERR250(II)(1:LENE(II)).EQ.ERB250(1:LENB)) THEN
            FOUND = .TRUE.
         ELSE
            II=II+1
         END IF
         BDONE=(FOUND.OR.II.GT.NRERR)
      END DO
      IF (FOUND) THEN
         ECNT(II)=ECNT(II)+1
      ELSE
         IF (NRERR.GE.MAXERR) THEN
            WRITE(*,*) MYNAME,
     &           'TOO MANY BIT-FLAGS SET, IGNORING FLAG:'//
     &           ERB250(1:LENB)
         ELSE
            NRERR=MIN(NRERR+1,MAXERR)
            ERR250(NRERR)=ERB250
            LENE(NRERR)=LENB
            ECNT(NRERR)=1
C     SORT ARRAY
            II=NRERR
            BDONE=(II.LE.1)
            DO WHILE (.NOT.BDONE)
               IF (ERR250(II)(1:LENE(II)).LT.
     &              ERR250(II-1)(1:LENE(II-1))) THEN
C     FLIP
                  ERR250B=ERR250(II)
                  LENEB=LENE(II)
                  ECNTB=ECNT(II)
                  ERR250(II)=ERR250(II-1)
                  LENE(II)=LENE(II-1)
                  ECNT(II)=ECNT(II-1)
                  ERR250(II-1)=ERR250B
                  LENE(II-1)=LENEB
                  ECNT(II-1)=ECNTB
                  II=II-1
                  BDONE=(II.LE.1)
               ELSE
                  BDONE=.TRUE.
               END IF
            END DO
         END IF
      END IF
      RETURN
      END
      

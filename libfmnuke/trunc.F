      SUBROUTINE TRUNC(STRING,N)
C     ***********************************************************************
C     +                                                                     *
C     REMOVES TRAILING BLANKS FROM A STRING                                  *
C     +                                                                     *
CIO   STRING (C*N) = STRING TO BE EXAMINED                                  *
CI    N (I*4) = LENGTH OF STRING                                            *
C     +                                                                     *
C     VERSION                      : 18/07/95                               *
C     +                                                                     *
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
      SAVE
C
      INTEGER N
      CHARACTER*(*) STRING
      CHARACTER*1 C1
      INTEGER II
      LOGICAL BDONE
C
      LOGICAL FIRST
      DATA FIRST /.TRUE./
C
      IF (FIRST) THEN
         C1=CHAR(0)
         FIRST=.FALSE.
      ENDIF
C
      II=N
      BDONE=(II.LE.0)
      DO WHILE (.NOT. BDONE)
         IF (STRING(II:II).EQ.' ') THEN
            STRING(II:II)=C1
            II=II-1
            BDONE= (II.LE.0)
         ELSE IF (STRING(II:II).EQ.C1) THEN
            II=II-1
         ELSE
            BDONE=.TRUE.
         END IF
      END DO
C
      RETURN
      END

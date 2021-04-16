      INTEGER FUNCTION LENGTH(STRING,N,GUESS)
C     ***********************************************************************
C     +                                                                     *
C     GIVES THE LENGTH OF THE STRING                                        *
C     +                                                                     *
CI    STRING (C*N) = STRING TO BE EXAMINED                                  *
CI    N (I*4) = LENGTH OF STRING                                            *
CI    GUESS (I*4) = FIRST GUESS FOR LENGTH                                  *
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
      INTEGER N,LL,LENS,GUESS
      CHARACTER*(*) STRING
      CHARACTER*1 C1
      LOGICAL BL,FORW,DONE
C
      LOGICAL FIRST
      DATA FIRST /.TRUE./
C
      IF (FIRST) THEN
         C1=CHAR(0)
         FIRST=.FALSE.
      ENDIF
C
      DONE=.FALSE.
      LENS=0
      if (n.lt.1) return
C
      LL=GUESS
      IF (LL.LT.1) LL=1
      IF (LL.GT.N) LL=N
C
      BL=(STRING(LL:LL).EQ.C1)
      IF (.NOT.BL)LENS=LL
      FORW=(.NOT.BL)
C
      IF (FORW) THEN
         LL=LL+1
      ELSE
         LL=LL-1
      ENDIF
C
      DO WHILE (LL.GT.0.AND.LL.LE.N.AND..NOT.DONE)
         BL=(STRING(LL:LL).EQ.C1)
         IF (.NOT.BL)LENS=LL
         DONE=((BL.AND.FORW).OR.(.NOT.BL.AND..NOT.FORW))
         IF (FORW) THEN
            LL=LL+1
         ELSE
            LL=LL-1
         ENDIF
      ENDDO
C
      LENGTH=LENS
C
      RETURN
      END

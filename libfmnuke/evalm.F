      SUBROUTINE EVALM(STRING,N,CHANGED,IRC)
C     ***********************************************************************
C     +                                                                     *
C     EVALUATES AN EXPRESSION                                               *
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
      INTEGER N,IRC
      CHARACTER*(*) STRING
      LOGICAL CHANGED
C
      CHANGED=.FALSE.
C
C-----PERFORM ALL MATHEMATICAL OPERATIONS
C
      CALL PERFM(STRING,N,CHANGED,IRC)
C
C-----CLEAN UP
C
      CALL CHOP(STRING,N)
      CALL CLEANUP(STRING,N)
      CALL TRUNC0(STRING,N)
C
      RETURN
      END

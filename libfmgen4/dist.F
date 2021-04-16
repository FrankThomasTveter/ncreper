      SUBROUTINE DIST(LONA,LATA,
     &     LONB,LATB,DIFF,IRC)
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
      REAL   LONA,LATA,LONB,LATB,DIFF
      INTEGER  IRC
C
      REAL   CDIFF
C
      real sindeg, cosdeg,acosdeg
      external sindeg,cosdeg,acosdeg
C
      CHARACTER*8 MYNAME 
      DATA MYNAME /'DIST'/
C
      CDIFF=SINDEG(LATA)*SINDEG(LATB)+
     &     COSDEG(LATA)*COSDEG(LATB)*COSDEG(LONA-LONB)
C
      CDIFF=MAX(-1.0D0,MIN(1.0D0,CDIFF)) ! handle truncation errors
C
      DIFF=ABS(ACOSDEG(CDIFF))
C
      RETURN
      END

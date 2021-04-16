      SUBROUTINE CLW( PCAT, TCLW, SATZEN, BTEMPS, BTEMPSC, IRC) 
C
      IMPLICIT NONE
      SAVE
C
      INTEGER  PCAT
      REAL   TCLW ! OUTPUT: TOTAL CLOUD LIQUID WATER
      REAL   SATZEN ! SATELLITE ZENITH ANGLE
      REAL   BTEMPS(2) ! OBSERVED AMSUA-A CHANNEL ONE AND TWO
      REAL   BTEMPSC(2) ! CALCULATED AMSUA-A CHANNEL ONE AND TWO (NO CLUD)
      INTEGER  IRC
C
      real sindeg, cosdeg
      external sindeg,cosdeg
C
      REAL   COSQ, D0, D1, D2, TCLW_O, TCLW_C
      
      IF (PCAT.EQ.1) THEN       ! OVER SEA, USE NESDIS TCLW ALGORITHM
         COSQ=COSDEG(SATZEN)
         D0=8.24D0 - 2.622D0*COSQ + 1.846D0*COSQ*COSQ
         D1 = 0.754D0 
         D2 = - 2.265D0 
         TCLW_O=MAX(COSQ*(D0+
     &        D1*LOG(285.0D0-BTEMPS(1))+
     &        D2*LOG(285.0D0-BTEMPS(2)))-
     &        0.03,0.0D0)
         TCLW_C=MAX(COSQ*(D0+
     &        D1*LOG(285.0D0-BTEMPSC(1))+
     &        D2*LOG(285.0D0-BTEMPSC(2)))-
     &        0.03,0.0D0)
         TCLW=TCLW_O-TCLW_C
      ELSE                      ! GROSS ERROR CHECK
         TCLW=0.0D0
         IF (ABS(BTEMPS(1)-BTEMPSC(1)).GT.4.0D0) TCLW=10.0D0
      END IF
C     
      RETURN
      END

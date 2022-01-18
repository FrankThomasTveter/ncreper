*DECK BINOM
      FUNCTION BINOM (N, M)
C***BEGIN PROLOGUE  BINOM
C***PURPOSE  Compute the binomial coefficients.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C1
C***TYPE      SINGLE PRECISION (BINOM-S, DBINOM-D)
C***KEYWORDS  BINOMIAL COEFFICIENTS, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BINOM(N,M) calculates the binomial coefficient (N!)/((M!)*(N-M)!).
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNREL, R1MACH, R9LGMC, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BINOM
      LOGICAL FIRST
      SAVE SQ2PIL, BILNMX, FINTMX, FIRST
      DATA SQ2PIL / 0.9189385332 0467274E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BINOM
      IF (FIRST) THEN
         BILNMX = LOG (R1MACH(2))
         FINTMX = 0.9/R1MACH(3)
      ENDIF
      FIRST = .FALSE.
C
      IF (N .LT. 0 .OR. M .LT. 0) CALL XERMSG ('SLATEC', 'BINOM',
     +   'N OR M LT ZERO', 1, 2)
      IF (N .LT. M) CALL XERMSG ('SLATEC', 'BINOM', 'N LT M', 2, 2)
C
      K = MIN (M, N-M)
      IF (K.GT.20) GO TO 30
      IF (K*LOG(AMAX0(N,1)).GT.BILNMX) GO TO 30
C
      BINOM = 1.
      IF (K.EQ.0) RETURN
C
      DO 20 I=1,K
        BINOM = BINOM * REAL(N-I+1)/I
 20   CONTINUE
C
      IF (BINOM.LT.FINTMX) BINOM = AINT (BINOM+0.5)
      RETURN
C
C IF K.LT.9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM
 30   IF (K .LT. 9) CALL XERMSG ('SLATEC', 'BINOM',
     +   'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
C
      XN = N + 1
      XK = K + 1
      XNK = N - K + 1
C
      CORR = R9LGMC(XN) - R9LGMC(XK) - R9LGMC(XNK)
      BINOM = XK*LOG(XNK/XK) - XN*ALNREL(-(XK-1.)/XN)
     1  - 0.5*LOG(XN*XNK/XK) + 1.0 - SQ2PIL + CORR
C
      IF (BINOM .GT. BILNMX) CALL XERMSG ('SLATEC', 'BINOM',
     +   'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
C
      BINOM = EXP (BINOM)
      IF (BINOM.LT.FINTMX) BINOM = AINT (BINOM+0.5)
C
      RETURN
      END
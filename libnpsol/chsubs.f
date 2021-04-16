*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  chsubs.f
*
*     chcore   chfd     chfdlc   chfdls   chfgrd   chcjac   chfjac
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE CHCORE( DEBUG, DONE, FIRST, EPSA, EPSR, FX, X,
     $                   INFORM, ITER, ITMAX,
     $                   CDEST, FDEST, SDEST, ERRBND, F1,
     $                   F2, H, HOPT, HPHI )

      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      LOGICAL            DEBUG, DONE, FIRST

C***********************************************************************
C     CHCORE  implements algorithm  FD, the method described in
C     Gill, P.E., Murray, W., Saunders, M.A., and Wright, M. H.,
C     Computing Forward-Difference Intervals for Numerical Optimization,
C     Siam Journal on Scientific and Statistical Computing, vol. 4,
C     pp. 310-321, June 1983.
C     
C     The procedure is based on finding an interval (HPHI) that
C     produces an acceptable estimate of the second derivative, and
C     then using that estimate to compute an interval that should
C     produce a reasonable forward-difference approximation.
C     
C     One-sided difference estimates are used to ensure feasibility with
C     respect to an upper or lower bound on X.  If X is close to an 
C     upper bound, the trial intervals will be negative.  The final 
C     interval is always positive.
C     
C     CHCORE has been designed to use a reverse communication
C     control structure, i.e., all evaluations of the function occur
C     outside this routine.  The calling routine repeatedly calls  
C     CHCORE  after computing the indicated function values.
C     
C     CHCORE  is similar to subroutine FDCORE described in Report
C     SOL 83-6, Documentation of FDCORE and FDCALC, by P.E. Gill,
C     W. Murray,  M.A. Saunders, and M.H. Wright, Department of
C     Operations Research,  Stanford University, Stanford, California
C     94305, June 1983.
C     
C     Systems Optimization Laboratory, Stanford University.
C     Based on Fortran 66 Version 2.1 of  FDCORE  written June 1983.
C     Fortran 77 Version written 25-May-1985.
C     This version of  CHCORE  dated  11-February-1986.
C***********************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            CE1BIG, CE2BIG, TE2BIG, OVERFL
      SAVE               CDSAVE, FDSAVE, HSAVE, OLDH, RHO, SDSAVE
      SAVE               CE1BIG, CE2BIG, TE2BIG
      EXTERNAL           DDIV
      INTRINSIC          ABS   , MAX   , MIN  , SQRT

      PARAMETER         (BNDLO  =1.0D-3, BNDUP  =1.0D-1                )

      PARAMETER         (ZERO   =0.0D+0, SIXTH  =1.6D-1, FOURTH =2.5D-1)
      PARAMETER         (HALF   =5.0D-1, ONE    =1.0D+0, TWO    =2.0D+0)
      PARAMETER         (THREE  =3.0D+0, FOUR   =4.0D+0, TEN    =1.0D+1)

*     ------------------------------------------------------------------
*     Explanation of local variables...
*
*     BNDLO, BNDUP, and RHO control the logic of the routine.
*     BNDLO and BNDUP are the lower and upper bounds that define an
*     acceptable value of the bound on the relative condition error in
*     the second derivative estimate.
*
*     The scalar RHO is the factor by which the interval is multiplied
*     or divided, and also the multiple of the well-scaled interval
*     that is used as the initial trial interval.
*
*     All these values are discussed in the documentation.
*     ------------------------------------------------------------------

      ITER  = ITER + 1

*     Compute the forward-,  backward-,  central-  and second-order
*     difference estimates.

      FDEST  = DDIV  ( F1 - FX,     H, OVERFL )
      FDEST2 = DDIV  ( F2 - FX, TWO*H, OVERFL )

      OLDCD = CDEST
      CDEST = DDIV  ( FOUR*F1 - THREE*FX - F2, TWO*H, OVERFL )

      OLDSD = SDEST
      SDEST = DDIV  ( FX      - TWO*F1   + F2, H*H  , OVERFL )

*     Compute  FDCERR  and  SDCERR,  bounds on the relative condition
*     errors in the first and second derivative estimates.

      AFDMIN = MIN( ABS( FDEST ), ABS( FDEST2 ) )
      FDCERR = DDIV  ( EPSA, HALF*ABS( H )*AFDMIN, OVERFL )
      SDCERR = DDIV  ( EPSA, FOURTH*ABS( SDEST )*H*H, OVERFL )

      IF (DEBUG)
     $   WRITE (NOUT, 9000) ITER  , FX   , H,
     $                      F1    , FDEST,
     $                      F2    , FDEST2,
     $                      CDEST , SDEST,
     $                      FDCERR, SDCERR

*     ==================================================================
*     Select the correct case.
*     ==================================================================
      IF (FIRST) THEN
*        ---------------------------------------------------------------
*        First time through.
*        Check whether SDCERR lies in the acceptable range.
*        ------------------------------------------------------------
         FIRST  = .FALSE.
         DONE   = SDCERR .GE. BNDLO  .AND.  SDCERR .LE. BNDUP
         TE2BIG = SDCERR .LT. BNDLO
         CE2BIG = SDCERR .GT. BNDUP
         CE1BIG = FDCERR .GT. BNDUP

         IF (.NOT. CE1BIG) THEN
            HSAVE  = H
            FDSAVE = FDEST
            CDSAVE = CDEST
            SDSAVE = SDEST
         END IF

         RHO  = EPSR**(-SIXTH)/FOUR
         IF (TE2BIG) THEN

*           The truncation error may be too big  (same as saying
*           SDCERR is too small).  Decrease the trial interval.

            RHO    = TEN*RHO
            OLDH   = H
            H      = H / RHO
         ELSE IF (CE2BIG) THEN

*           SDCERR is too large.  Increase the trial interval.

            OLDH   = H
            H      = H*RHO
         END IF
      ELSE IF (CE2BIG) THEN
*        ---------------------------------------------------------------
*        During the last iteration,  the trial interval was
*        increased in order to decrease SDCERR.
*        ---------------------------------------------------------------
         IF (CE1BIG  .AND.  FDCERR .LE. BNDUP) THEN
            CE1BIG = .FALSE.
            HSAVE  = H
            FDSAVE = FDEST
            CDSAVE = CDEST
            SDSAVE = SDEST
         END IF

*        If SDCERR is small enough, accept H.  Otherwise,
*        increase H again.

         DONE   = SDCERR .LE. BNDUP
         IF (.NOT. DONE) THEN
            OLDH   = H
            H      = H*RHO
         END IF
      ELSE IF (TE2BIG) THEN
*        ---------------------------------------------------------------
*        During the last iteration,  the interval was decreased in order
*        to reduce the truncation error.
*        ---------------------------------------------------------------
         DONE   = SDCERR .GT. BNDUP
         IF (DONE) THEN

*           SDCERR has jumped from being too small to being too
*           large.  Accept the previous value of H.

            H     = OLDH
            SDEST = OLDSD
            CDEST = OLDCD
         ELSE

*           Test whether FDCERR is sufficiently small.

            IF (FDCERR .LE. BNDUP) THEN
               CE1BIG = .FALSE.
               HSAVE  = H
               FDSAVE = FDEST
               CDSAVE = CDEST
               SDSAVE = SDEST
            END IF

*           Check whether SDCERR is in range.

            DONE  = SDCERR .GE. BNDLO

            IF (.NOT. DONE) THEN

*              SDCERR is still too small, decrease H again.

               OLDH = H
               H    = H / RHO
            END IF
         END IF

      END IF

*     ==================================================================
*     We have either finished or have a new estimate of H.
*     ==================================================================
      IF (DONE) THEN

*        Sufficiently good second-derivative estimate found.
*        Compute the optimal interval.

         HPHI   = ABS( H )
         HOPT   = TWO * SQRT( EPSA ) / SQRT( ABS( SDEST ) )

*        ERR1 is the error bound on the forward-difference estimate
*        with the final value of H.  ERR2 is the difference of FDEST
*        and the central-difference estimate with HPHI.

         ERR1   = HOPT*ABS( SDEST )
         ERR2   = ABS( FDEST - CDEST )
         ERRBND = MAX( ERR1, ERR2 )

*        Set INFORM = 4  if the forward- and central-difference
*        estimates are not close.

         INFORM = 0
         IF (ERRBND .GT. HALF*ABS( FDEST )) INFORM = 4
      ELSE
*        ---------------------------------------------------------------
*        Check whether the maximum number of iterations has been
*        exceeded.  If not, exit.
*        ---------------------------------------------------------------
         DONE = ITER .GE. ITMAX
         IF (DONE) THEN
            IF (CE1BIG) THEN

*              FDCERR was never small.  Probably a constant function.

               INFORM = 1
               HPHI   = HOPT
               FDEST  = ZERO
               CDEST  = ZERO
               SDEST  = ZERO
               ERRBND = ZERO
            ELSE IF (CE2BIG) THEN

*              FDCERR was small,  but SDCERR was never small.
*              Probably a linear or odd function.

               INFORM = 2
               HPHI   = ABS( HSAVE )
               HOPT   = HPHI
               FDEST  = FDSAVE
               CDEST  = CDSAVE
               SDEST  = ZERO
               ERRBND = TWO*EPSA / HOPT
            ELSE

*              The only remaining case occurs when the second
*              derivative is changing too rapidly for an adequate
*              interval to be found (SDCERR remained small even
*              though H was decreased ITMAX times).

               INFORM = 3
               HPHI   = ABS( HSAVE )
               HOPT   = HPHI
               FDEST  = FDSAVE
               CDEST  = CDSAVE
               SDEST  = SDSAVE
               ERRBND = HOPT*ABS( SDEST )/TWO + TWO*EPSA/HOPT
            END IF
         END IF
      END IF

      IF (DEBUG) THEN
         WRITE (NOUT, 9001) CE1BIG, CE2BIG, TE2BIG
         IF (DONE)
     $      WRITE (NOUT, 9002) INFORM, HOPT, ERRBND
      END IF

      RETURN

 9000 FORMAT(/ ' //CHCORE//  ITN ', I3,
     $                             ' FX     H      ', 5X, 1P, 2D16.6
     $       / ' //CHCORE//  F1      FDEST         ', 5X, 1P, 2D16.6
     $       / ' //CHCORE//  F2      FDEST2        ', 5X, 1P, 2D16.6
     $       / ' //CHCORE//  CDEST   SDEST         ', 5X, 1P, 2D16.6
     $       / ' //CHCORE//  FDCERR  SDCERR        ', 5X, 1P, 2D16.6)
 9001 FORMAT(  ' //CHCORE//  CE1BIG  CE2BIG  TE2BIG', 5X, 3L2       )
 9002 FORMAT(  ' //CHCORE//  INFORM  HOPT    ERRBND', I5, 1P, 2D16.6)

*     End of  CHCORE.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE CHFD  ( INFORM, MSGLVL, LVLDER,
     $                   N, NCNLN, LDCJ, LDCJU,
     $                   BIGBND, EPSRF, FDNORM, OBJF,
     $                   OBJFUN, CONFUN, NEEDC,
     $                   BL, BU, C, C1, C2, CJAC, CJACU,
     $                   GRAD, GRADU, HFORWD, HCNTRL,
     $                   X, Y, W, LENW )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            NEEDC(*)
      DOUBLE PRECISION   BL(N), BU(N)
      DOUBLE PRECISION   C(*), C1(*), C2(*),
     $                   CJAC(LDCJ,*), CJACU(LDCJU,*)
      DOUBLE PRECISION   GRAD(N), GRADU(N)
      DOUBLE PRECISION   HFORWD(*), HCNTRL(*)
      DOUBLE PRECISION   X(N), Y(N), W(LENW)
      EXTERNAL           OBJFUN, CONFUN

************************************************************************
*  CHFD    computes difference intervals for the missing gradients of
*  F(x) and c(x).  Intervals are computed using a procedure that usually
*  requires about two function evaluations if the function is well
*  scaled.  Central-difference gradients are obtained as a by-product
*  of the algorithm.
*
*  On entry...
*     OBJF and C contain the problem functions at the point X.
*     An element of CJAC or GRAD not equal to RDUMMY signifies a known
*     gradient value.  Such values are not estimated by differencing.
*     CJACU and GRADU have dummy elements in the same positions as
*     CJAC and GRADU.
*
*  On exit...                          
*     CJAC and GRAD contain central-difference derivative estimates.
*     Elements of CJACU and GRADU are unaltered except for those
*     corresponding to constant derivatives, which are given the same
*     values as CJAC or GRAD.
*
*  The work array  W  is not used at present.
*
*  Systems Optimization Laboratory, Department of Operations Research,
*  Stanford University, Stanford, California 94305
*  Original version written 28-July-1985.
*  This version of CHFD   dated 14-July-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9

      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET

      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG

      LOGICAL            DEBUG , DONE  , FIRST , HEADNG, NEEDED
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      EXTERNAL           DNRM2
      PARAMETER         (RDUMMY =-11111.0              )
      PARAMETER         (FACTOR =0.97D+0               )
      PARAMETER         (ZERO   =0.0D+0, HALF   =0.5D+0, ONE   =1.0D+0)
      PARAMETER         (TWO    =2.0D+0, FOUR   =4.0D+0, TEN   =1.0D+1)

      INFORM = 0
      NEEDED = LVLDER .EQ. 0  .OR.  LVLDER .EQ. 2
     $                        .OR.  LVLDER .EQ. 1  .AND.  NCNLN .GT. 0
      IF (.NOT. NEEDED) RETURN

      DEBUG  = NPDBG  .AND.  INPDBG(5) .GT. 0
      IF (LFDSET .EQ. 0) THEN
         IF (MSGLVL .GT. 0) WRITE (NOUT, 1000)

         NSTATE = 0
         ITMAX  = 3
         MODE   = 0

         NCCNST = 0
         NFCNST = 0
         HEADNG = .TRUE.

         FDNORM = ZERO

*        ===============================================================
*        For each column of the Jacobian augmented by the transpose of
*        the objective gradient, rows IROW1 thru IROW2 are searched for
*        missing elements.
*        ===============================================================
         IROW1  = 1
         IROW2  = NCNLN + 1
         IF (LVLDER .EQ. 1) IROW2 = NCNLN
         IF (LVLDER .EQ. 2) IROW1 = NCNLN + 1

         BIGLOW = - BIGBND
         BIGUPP =   BIGBND

         IF (NCNLN  .GT. 0)
     $      CALL ILOAD ( NCNLN, (0), NEEDC, 1 )

         DO 600 J = 1, N
            XJ     = X(J)
            NFOUND = 0
            SUMSD  = ZERO
            SUMEPS = ZERO
            HFD    = ZERO
            HCD    = ZERO
            HMAX   = ZERO
            HMIN   = ONE / EPSPT3
            ERRMAX = ZERO
            ERRMIN = ZERO

            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J) .GT. BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J) .LT. BIGUPP) STEPBU = BU(J) - XJ

            SIGNH  = ONE
            IF (HALF*(STEPBL + STEPBU) .LT. ZERO) SIGNH =  - ONE

            DO 500 I = IROW1, IROW2

               IF (I .LE. NCNLN) THEN
                  TEST = CJACU(I,J)
               ELSE
                  TEST = GRADU(J)
               END IF

               IF (TEST .EQ. RDUMMY) THEN
*                 ======================================================
*                 Get the difference interval for this component.
*                 ======================================================
                  NFOUND = NFOUND + 1

                  IF (I .LE. NCNLN) THEN
                     NEEDC(I) = 1
                     FX       = C(I)
                     EPSA     = EPSRF*(ONE + ABS( C(I) ))
                  ELSE
                     FX       = OBJF
                     EPSA     = EPSRF*(ONE + ABS( FX ))
                  END IF

*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  ITER   = 0
                  HOPT   = TWO*(ONE + ABS( XJ ))*SQRT( EPSRF )
                  H      = SIGNH*TEN*HOPT
                  CDEST  = ZERO
                  SDEST  = ZERO
                  FIRST  = .TRUE.

*+                REPEAT
  400                X(J)  = XJ + H
                     IF (I .LE. NCNLN) THEN
                        CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                               NEEDC, X, C1, CJACU, NSTATE )
                        IF (MODE .LT. 0) GO TO 9999
                        F1 = C1(I)
                     ELSE
                        CALL OBJFUN( MODE, N, X, F1, GRADU, NSTATE )
                        IF (MODE .LT. 0) GO TO 9999
                     END IF

                     X(J)  = XJ + H + H
                    IF (I .LE. NCNLN) THEN
                       CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                              NEEDC, X, C1, CJACU, NSTATE )
                        IF (MODE .LT. 0) GO TO 9999
                        F2 = C1(I)
                     ELSE
                        CALL OBJFUN( MODE, N, X, F2, GRADU, NSTATE )
                        IF (MODE .LT. 0) GO TO 9999
                     END IF

                     CALL CHCORE( DEBUG, DONE, FIRST, EPSA, EPSRF,FX,XJ,
     $                            INFO, ITER, ITMAX,
     $                            CDEST, FDEST, SDEST, ERRBND, F1,
     $                            F2, H, HOPT, HPHI )
      
*+                UNTIL     DONE
                  IF (.NOT. DONE) GO TO 400

                  IF (I .LE. NCNLN) THEN
                     CJAC(I,J) = CDEST
                     IF (INFO .EQ. 1  .OR.  INFO .EQ. 2) THEN
                        NCCNST    =   NCCNST + 1
                        NCDIFF    =   NCDIFF - 1
                        CJACU(I,J) = - RDUMMY
                     END IF
                  ELSE
                     GRAD(J)   = CDEST
                     IF (INFO .EQ. 1  .OR.  INFO .EQ. 2) THEN
                        NFCNST    =   NFCNST + 1
                        NFDIFF    =   NFDIFF - 1
                        GRADU(J)  = - RDUMMY
                     END IF
                  END IF

                  SUMSD  = SUMSD  + ABS( SDEST )
                  SUMEPS = SUMEPS +      EPSA
                  IF (HOPT .GT. HMAX) THEN
                     HMAX   = HOPT
                     ERRMAX = ERRBND
                  END IF
                  IF (HOPT .LT. HMIN) THEN
                     HMIN   = HOPT
                     ERRMIN = ERRBND
                  END IF

                  IF (INFO .EQ. 0) HCD  = MAX ( HCD, HPHI )
               END IF
  500       CONTINUE

            IF (NFOUND .GT. 0) THEN
               IF (HMIN .GT. HMAX) THEN
                  HMIN   = HMAX
                  ERRMIN = ERRMAX
               END IF

               IF      (FOUR*SUMEPS .LT. HMIN*HMIN*SUMSD) THEN
                  HFD    = HMIN
                  ERRMAX = ERRMIN
               ELSE IF (FOUR*SUMEPS .GT. HMAX*HMAX*SUMSD) THEN
                  HFD    = HMAX
               ELSE
                  HFD    = TWO*SQRT( SUMEPS / SUMSD )
                  ERRMAX = TWO*SQRT( SUMEPS * SUMSD )
               END IF

               IF (HCD .EQ. ZERO) HCD = TEN*HFD

               IF (MSGLVL .GT. 0) THEN
                  IF (HEADNG) WRITE (NOUT, 1100)
                  WRITE (NOUT, 1200) J, XJ, HFD, HCD, ERRMAX
                  HEADNG = .FALSE.
               END IF
            END IF

            FDNORM    = MAX (FDNORM, HFD)
            HFORWD(J) = HFD / (ONE + ABS(XJ))
            HCNTRL(J) = HCD / (ONE + ABS(XJ))
            X(J)      = XJ
  600    CONTINUE

         IF (NCCNST + NFCNST .GT. 0) THEN

*           Check that the constants have been set properly by
*           evaluating the gradients at a strange (but feasible) point.

            D      =   ONE / N

            DO 710 J = 1, N
               XJ     =   X(J)
               STEPBL = - ONE
               STEPBU =   ONE
               IF (BL(J) .GT. BIGLOW)
     $            STEPBL = MAX( STEPBL, BL(J) - XJ )
               IF (BU(J) .LT. BIGUPP  .AND.  BU(J) .GT. BL(J))
     $            STEPBU = MIN( STEPBU, BU(J) - XJ )

               IF (HALF*(STEPBL + STEPBU) .LT. ZERO) THEN
                  Y(J) = XJ + D*STEPBL
               ELSE
                  Y(J) = XJ + D*STEPBU
               END IF

               D = FACTOR*D
  710       CONTINUE

            IF (NCNLN .GT. 0) THEN
               CALL ILOAD ( NCNLN, (1), NEEDC, 1 )
               CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                      NEEDC, Y, C2, CJACU, NSTATE )
               IF (MODE .LT. 0) GO TO 9999
            END IF

            CALL OBJFUN( MODE, N, Y, OBJF2, GRADU, NSTATE )
            IF (MODE .LT. 0) GO TO 9999

*           ------------------------------------------------------------
*           Loop over each of the components of  x.
*           ------------------------------------------------------------
            DO 800 J = 1, N
               YJ     = Y(J)
               DX     = HALF*(X(J) - YJ)
               Y(J)   = YJ + DX

               IF (NCNLN .GT. 0) THEN
                  NFOUND = 0
                  DO 720 I = 1, NCNLN
                     IF (CJACU(I,J) .EQ. - RDUMMY) THEN
                        NEEDC(I) = 1
                        NFOUND   = NFOUND + 1
                     ELSE
                        NEEDC(I) = 0
                     END IF
  720             CONTINUE

                  IF (NFOUND .GT. 0) THEN
                     CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                            NEEDC, Y, C1, CJACU, NSTATE )
                     IF (MODE .LT. 0) GO TO 9999

                     DO 730 I = 1, NCNLN
                        IF (NEEDC(I) .EQ. 1) THEN
                           CJDIFF = ( C1(I) -  C2(I) ) / DX
                           IF (CJDIFF .EQ. CJAC(I,J)) THEN
                              CJACU(I,J) = CJDIFF
                           ELSE
                              CJACU(I,J) = RDUMMY
                              NCCNST    = NCCNST - 1
                              NCDIFF    = NCDIFF + 1
                           END IF
                        END IF
  730                CONTINUE
                  END IF
               END IF

*              Now check the objective gradient component.

               IF (GRADU(J) .EQ. - RDUMMY) THEN

                  CALL OBJFUN( MODE, N, Y, F1, GRADU, NSTATE )
                  IF (MODE .LT. 0) GO TO 9999

                  GDIFF = (F1 - OBJF2)/DX
                  IF (GDIFF .EQ. GRAD(J)) THEN
                     GRADU(J) = GDIFF
                  ELSE
                     GRADU(J) = RDUMMY
                     NFDIFF   = NFDIFF + 1
                     NFCNST   = NFCNST - 1
                  END IF
               END IF

               Y(J)  = YJ
  800       CONTINUE

            IF (MSGLVL .GT. 0) THEN
               IF (LVLDER .LT. 2  .AND.  NCCNST .GT. 0)
     $            WRITE (NOUT, 1300) NCCNST
               IF (LVLDER .NE. 1  .AND.  NFCNST .GT. 0)
     $            WRITE (NOUT, 1400) NFCNST
            END IF

            IF (NCDIFF .EQ. 0  .AND.  LVLDER .LT. 2) THEN
               IF (LVLDER .EQ. 0) LVLDER = 2
               IF (LVLDER .EQ. 1) LVLDER = 3
               IF (MSGLVL .GT. 0) WRITE (NOUT, 1500) LVLDER
            END IF

            IF (NFDIFF .EQ. 0  .AND.  LVLDER .NE. 1) THEN
               IF (LVLDER .EQ. 0) LVLDER = 1
               IF (LVLDER .EQ. 2) LVLDER = 3
               IF (MSGLVL .GT. 0) WRITE (NOUT, 1600) LVLDER
            END IF
         END IF
      ELSE IF (LFDSET .EQ. 2) THEN

*        The user has supplied HFORWD and HCNTRL.
*        Check for wild values.

         DO 900 J = 1, N
            IF (HFORWD(J) .LE. ZERO) THEN
               if(msglvl.ne.0) WRITE (NOUT, 2000) J, HFORWD(J), EPSPT5
               HFORWD(J) = EPSPT5
            END IF
  900    CONTINUE
         DO 910 J = 1, N
            IF (HCNTRL(J) .LE. ZERO) THEN
               if(msglvl.ne.0) WRITE (NOUT, 2100) J, HCNTRL(J), EPSPT3
               HCNTRL(J) = EPSPT3
            END IF
  910    CONTINUE
      END IF

      RETURN

 9999 INFORM = MODE
      RETURN

 1000 FORMAT(//' Computation of the finite-difference intervals'
     $       / ' ----------------------------------------------' )
 1100 FORMAT(//'    J      X(J)   Forward DX(J)   Central DX(J) ',
     $         '     Error est.' /)
 1200 FORMAT(  I5, 1P, E10.2, E16.6, 2E16.6 )
 1300 FORMAT(/ I5,  '  constant constraint gradient elements assigned.')
 1400 FORMAT(/ I5,  '  constant  objective gradient elements assigned.')
 1500 FORMAT(//' All missing Jacobian elements are constants.  ',
     $         ' Derivative level increased to ', I4 )
 1600 FORMAT(//' All missing objective gradients are constants.  ',
     $         ' Derivative level increased to ', I4 )
 2000 FORMAT(' XXX  ', I4,'-th difference interval ',         1PE10.2,
     $       ' replaced by ', 1PE10.2 )
 2100 FORMAT(' XXX  ', I4,'-th central-difference interval ', 1PE10.2,
     $       ' replaced by ', 1PE10.2 )

*     End of  CHFD  .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chfdlc( inform, msglvl, lvlder, n, 
     $                   bigbnd, epsrf, fdnorm, objf,
     $                   objfun, bl, bu, 
     $                   grad, gradu, hforwd, hcntrl,
     $                   x, y, w, lenw )

      implicit           double precision(a-h,o-z)
      double precision   bl(n), bu(n)
      double precision   grad(n), gradu(n)
      double precision   hforwd(*), hcntrl(*)
      double precision   x(n), y(n), w(lenw)
      external           objfun

C***********************************************************************
C     CHFDLC  computes difference intervals for the missing gradients of
C     f(x).  Difference intervals are computed using a procedure that 
C     usually requires about two evaluations if  f(x)  is well scaled.
C     Central-difference gradients are obtained as a by-product of the
C     algorithm.
C
C     On entry...
C       OBJF  contains the objective function at the point X.
C       An element of  GRAD  not equal to  RDUMMY  signifies a user-
C       specified value.  Such values are not estimated by differencing.
C       GRAD  has dummy elements in the same positions as GRADU.
C
C     On exit...                          
C       GRAD  contain central-difference derivative estimates.
C       Elements of  GRADU  are unaltered except for those corresponding 
C       to constant derivatives, which are given values values found
C       by differencing.
C
C     The work array  W  is not used at present.
C
C     Original version written  5-July-1990.
C     This version of  CHFDLC  dated  5-Jul-1990.
C***********************************************************************
      common    /sol1cm/ nout
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            debug , done  , first , headng, needed
      intrinsic          abs   , max   , min   , sqrt
      external           dnrm2
      parameter         (rdummy =-11111.0              )
      parameter         (factor =0.97d+0               )
      parameter         (zero   =0.0d+0, half   =0.5d+0, one   =1.0d+0)
      parameter         (two    =2.0d+0, ten    =1.0d+1)

      inform = 0
      needed = lvlder .eq. 0  .or.  lvlder .eq. 2
      if (.not. needed) return

      debug  = npdbg  .and.  inpdbg(5) .gt. 0
      if (lfdset .eq. 0) then
         if (msglvl .gt. 0) write (nout, 1000)

         nstate = 0
         itmax  = 3
         mode   = 0

         nfcnst = 0
         headng = .true.

         fdnorm = zero

         biglow = - bigbnd
         bigupp =   bigbnd

         do 600, j = 1, n
            xj     = x(j)

            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

            signh  = one
            if (half*(stepbl + stepbu) .lt. zero) signh =  - one

            if (gradu(j) .eq. rdummy) then
*              =========================================================
*              This component needs a difference interval.
*              =========================================================
               epsa     = epsrf*(one + abs( objf ))
               iter   = 0
               hopt   = two*(one + abs( xj ))*sqrt( epsrf )
               h      = signh*ten*hopt
               cdest  = zero
               sdest  = zero
               first  = .true.

*+             repeat
  400             x(j)  = xj + h
                  call objfun( mode, n, x, f1, gradu, nstate )
                  if (mode .lt. 0) go to 9999

                  x(j)  = xj + h + h
                  call objfun( mode, n, x, f2, gradu, nstate )
                  if (mode .lt. 0) go to 9999

                  call chcore( debug, done, first, epsa, epsrf,
     $                         objf, xj, info, iter, itmax,
     $                         cdest, fdest, sdest, errbnd, f1,
     $                         f2, h, hopt, hphi )
*+             until     done
               if (.not. done) go to 400

               grad(j) = cdest
               if (info .eq. 1  .or.  info .eq. 2) then
                  nfcnst   =   nfcnst + 1
                  nfdiff   =   nfdiff - 1
                  gradu(j) = - rdummy
               end if

               if (info .eq. 0) then
                  hfd = hopt
                  hcd = hphi
               else 
                  hfd = two*(one + abs( xj ))*sqrt( epsrf )
                  hcd = ten*hfd
               end if

               if (msglvl .gt. 0) then
                  if (headng) write (nout, 1100)
                  write (nout, 1200) j, xj, hfd, hcd, errbnd
                  headng = .false.
               end if
            end if

            fdnorm    = max (fdnorm, hfd)
            hforwd(j) = hfd / (one + abs(xj))
            hcntrl(j) = hcd / (one + abs(xj))
            x(j)      = xj
  600    continue

         if (nfcnst .gt. 0) then

*           Check that the constant elements of  g  have been identified
*           correctly.  Re-evaluate the gradients at a strange feasible
*           point.

            d      =   one / n

            do 710, j = 1, n
               xj     =   x(j)
               stepbl = - one
               stepbu =   one
               if (bl(j) .gt. biglow)
     $            stepbl = max( stepbl, bl(j) - xj )
               if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $            stepbu = min( stepbu, bu(j) - xj )

               if (half*(stepbl + stepbu) .lt. zero) then
                  y(j) = xj + d*stepbl
               else
                  y(j) = xj + d*stepbu
               end if

               d = factor*d
  710       continue

            call objfun( mode, n, y, objf2, gradu, nstate )
            if (mode .lt. 0) go to 9999

*           ------------------------------------------------------------
*           Loop over each component of  x.
*           ------------------------------------------------------------
            do 800, j = 1, n
               yj     = y(j)
               dx     = half*(x(j) - yj)
               y(j)   = yj + dx

               if (gradu(j) .eq. - rdummy) then

                  call objfun( mode, n, y, f1, gradu, nstate )
                  if (mode .lt. 0) go to 9999

                  gdiff = (f1 - objf2)/dx
                  if (gdiff .eq. grad(j)) then
                     gradu(j) = gdiff
                  else
                     gradu(j) = rdummy
                     nfdiff   = nfdiff + 1
                     nfcnst   = nfcnst - 1
                  end if
               end if

               y(j)  = yj
  800       continue

            if (msglvl .gt. 0) then
               if (nfcnst .gt. 0)
     $            write (nout, 1400) nfcnst
            end if

            if (nfdiff .eq. 0) then
               if (lvlder .eq. 0) lvlder = 1
               if (lvlder .eq. 2) lvlder = 3
               if (msglvl .gt. 0) write (nout, 1600) lvlder
            end if
         end if
      else if (lfdset .eq. 2) then

*        The user has supplied  hforwd  and  hcntrl.
*        Check for wild values.

         do 900, j = 1, n
            if (hforwd(j) .le. zero) then
               if(msglvl.ne.0) write (nout, 2000) j, hforwd(j), epspt5
               hforwd(j) = epspt5
            end if
  900    continue
         do 910, j = 1, n
            if (hcntrl(j) .le. zero) then
               if(msglvl.ne.0) write (nout, 2100) j, hcntrl(j), epspt3
               hcntrl(j) = epspt3
            end if
  910    continue
      end if

      return

 9999 inform = mode

 1000 format(//' Computation of the finite-difference intervals'
     $       / ' ----------------------------------------------' )
 1100 format(//'    j      x(j)   Forward dx(j)   Central dx(j) ',
     $         '     Error est.' /)
 1200 format(  i5, 1p, e10.2, e16.6, 2e16.6 )
 1300 format(/ i5,  '  constant constraint gradient elements assigned.')
 1400 format(/ I5,  '  constant  objective gradient elements assigned.')
 1600 format(//' All missing objective gradients are constants.  ',
     $         ' Derivative level increased to ', i4 )
 2000 format(' XXX  ', i4,'-th difference interval ',         1pe10.2,
     $       ' replaced by ', 1pe10.2 )
 2100 format(' XXX  ', i4,'-th central-difference interval ', 1pe10.2,
     $       ' replaced by ', 1pe10.2 )

*     End of  CHFDLC.

      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE CHFDLS( INFORM, MSGLVL, LVLDER,
     $                   M, N, NCNLN, LDCJ, LDCJU, LDFJ, LDFJU,
     $                   BIGBND, EPSRF, FDNORM,
     $                   CONFUN, OBJFUN, NEEDC,
     $                   BL, BU, C, C1, C2, CJAC, CJACU,
     $                   F, F1, F2, FJAC, FJACU, HFORWD, HCNTRL,
     $                   X, Y, W, LENW )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            NEEDC(*)
      DOUBLE PRECISION   BL(N), BU(N)
      DOUBLE PRECISION   C(*), C1(*), C2(*),
     $                   CJAC(LDCJ,*), CJACU(LDCJU,*)
      DOUBLE PRECISION   F(*), F1(*), F2(*),
     $                   FJAC(LDFJ,*), FJACU(LDFJU,*)
      DOUBLE PRECISION   HFORWD(*), HCNTRL(*)
      DOUBLE PRECISION   X(N), Y(N), W(LENW)
      EXTERNAL           OBJFUN, CONFUN

C***********************************************************************
C     CHFDLS  computes difference intervals for missing Jacobian
C     elements of f(x) and c(x).  Intervals are computed using a
C     procedure that requires about two function evaluations if the
C     function is well scaled.  Central-difference gradients are
C     obtained as a by-product of the algorithm.
C
C     On entry...
C        F and C contain the problem functions at the point X.
C        An element of CJAC or FJAC not equal to RDUMMY signifies a known
C        gradient value.  Such values are not estimated by differencing.
C        CJACU and FJACU have dummy elements in the same positions as
C        CJAC and FJAC.
C
C     On exit...
C        CJAC and FJAC contain central-difference derivative estimates.
C        Elements of CJACU and FJACU are unaltered except for those
C        corresponding to constant derivatives, which are given the same
C        values as CJAC or FJAC.
C
C     The work array  W  is not used at present.
C
C     Systems Optimization Laboratory,
C     Department of Operations Research,
C     Stanford University, Stanford, California 94305
C     Original version written 10-May-1988.
C     This version of CHFDLS dated  8-Jun-1989.
C***********************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9

      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET

      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG

      LOGICAL            DEBUG , DONE  , FIRST , HEADNG, NEEDED
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      EXTERNAL           DNRM2
      PARAMETER         (RDUMMY =-11111.0              )
      PARAMETER         (FACTOR =0.97D+0               )
      PARAMETER         (ZERO   =0.0D+0, HALF   =0.5D+0, ONE   =1.0D+0)
      PARAMETER         (TWO    =2.0D+0, FOUR   =4.0D+0, TEN   =1.0D+1)

      INFORM = 0
      NEEDED = LVLDER .EQ. 0  .OR.  LVLDER .EQ. 2
     $                        .OR.  LVLDER .EQ. 1  .AND.  NCNLN .GT. 0
      IF (.NOT. NEEDED) RETURN

      DEBUG  = NPDBG  .AND.  INPDBG(5) .GT. 0
      IF (LFDSET .EQ. 0) THEN
         IF (MSGLVL .GT. 0) WRITE (NOUT, 1000)

         NSTATE = 0
         ITMAX  = 3
         MODE   = 0

         NCCNST = 0
         NFCNST = 0
         HEADNG = .TRUE.

         FDNORM = ZERO

*        ===============================================================
*        For each column of the matrix
*              ( CJAC  )
*              ( FJAC  ),
*        rows IROW1 thru IROW2 are searched for missing elements.
*        ===============================================================
         IROW1  = 1
         IROW2  = NCNLN + M
         IF (LVLDER .EQ. 1) IROW2 = NCNLN
         IF (LVLDER .EQ. 2) IROW1 = NCNLN + 1

         BIGLOW = - BIGBND
         BIGUPP =   BIGBND

         IF (NCNLN  .GT. 0)
     $      CALL ILOAD ( NCNLN, (0), NEEDC, 1 )

         DO 600, J = 1, N
            XJ     = X(J)
            NFOUND = 0
            SUMSD  = ZERO
            SUMEPS = ZERO
            HFD    = ZERO
            HCD    = ZERO
            HMAX   = ZERO
            HMIN   = ONE / EPSPT3
            ERRMAX = ZERO
            ERRMIN = ZERO

            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J) .GT. BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J) .LT. BIGUPP) STEPBU = BU(J) - XJ

            SIGNH  = ONE
            IF (HALF*(STEPBL + STEPBU) .LT. ZERO) SIGNH =  - ONE

            DO 500, IROW = IROW1, IROW2

               IF (IROW .LE. NCNLN) THEN
                  I    = IROW
                  TEST = CJACU(I,J)
               ELSE
                  I    = IROW - NCNLN
                  TEST = FJACU(I,J)
               END IF

               IF (TEST .EQ. RDUMMY) THEN
*                 ======================================================
*                 Get the difference interval for this component.
*                 ======================================================
                  NFOUND = NFOUND + 1

                  IF (IROW .LE. NCNLN) THEN
                     NEEDC(I) = 1
                     FX       = C(I)
                  ELSE
                     FX       = F(I)
                  END IF

                  EPSA     = EPSRF*(ONE + ABS( FX ))

*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  ITER   = 0
                  HOPT   = TWO*(ONE + ABS( XJ ))*SQRT( EPSRF )
                  H      = SIGNH*TEN*HOPT
                  CDEST  = ZERO
                  SDEST  = ZERO
                  FIRST  = .TRUE.

*+                REPEAT
  400                X(J)  = XJ + H
                     IF (IROW .LE. NCNLN) THEN
                        CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                               NEEDC, X, C1, CJACU, NSTATE )
                        IF (MODE .LT. 0) GO TO 9999
                        FFORW = C1(I)
                     ELSE
                        CALL OBJFUN( MODE, M    , N, LDFJU,
     $                                      X, F1, FJACU, NSTATE )
                        IF (MODE .LT. 0) GO TO 9999
                        FFORW = F1(I)
                     END IF

                     X(J)  = XJ + H + H
                    IF (IROW .LE. NCNLN) THEN
                       CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                              NEEDC, X, C1, CJACU, NSTATE )
                        IF (MODE .LT. 0) GO TO 9999
                        FBACK = C1(I)
                     ELSE
                        CALL OBJFUN( MODE, M    , N, LDFJU,
     $                                      X, F1, FJACU, NSTATE )
                        IF (MODE .LT. 0) GO TO 9999
                        FBACK = F1(I)
                     END IF

                     CALL CHCORE( DEBUG, DONE, FIRST, EPSA, EPSRF,FX,XJ,
     $                            INFO, ITER, ITMAX,
     $                            CDEST, FDEST, SDEST, ERRBND, FFORW,
     $                            FBACK, H, HOPT, HPHI )

*+                UNTIL     DONE
                  IF (.NOT. DONE) GO TO 400

                  IF (IROW .LE. NCNLN) THEN
                     CJAC(I,J) = CDEST
                     IF (INFO .EQ. 1  .OR.  INFO .EQ. 2) THEN
                        NCCNST     =   NCCNST + 1
                        NCDIFF     =   NCDIFF - 1
                        CJACU(I,J) = - RDUMMY
                     END IF
                  ELSE
                     FJAC(I,J) = CDEST
                     IF (INFO .EQ. 1  .OR.  INFO .EQ. 2) THEN
                        NFCNST     =   NFCNST + 1
                        NFDIFF     =   NFDIFF - 1
                        FJACU(I,J) = - RDUMMY
                     END IF
                  END IF

                  SUMSD  = SUMSD  + ABS( SDEST )
                  SUMEPS = SUMEPS +      EPSA
                  IF (HOPT .GT. HMAX) THEN
                     HMAX   = HOPT
                     ERRMAX = ERRBND
                  END IF
                  IF (HOPT .LT. HMIN) THEN
                     HMIN   = HOPT
                     ERRMIN = ERRBND
                  END IF

                  IF (INFO .EQ. 0) HCD  = MAX ( HCD, HPHI )
               END IF
  500       CONTINUE

            IF (NFOUND .GT. 0) THEN
               IF (HMIN .GT. HMAX) THEN
                  HMIN   = HMAX
                  ERRMIN = ERRMAX
               END IF

               IF      (FOUR*SUMEPS .LT. HMIN*HMIN*SUMSD) THEN
                  HFD    = HMIN
                  ERRMAX = ERRMIN
               ELSE IF (FOUR*SUMEPS .GT. HMAX*HMAX*SUMSD) THEN
                  HFD    = HMAX
               ELSE
                  HFD    = TWO*SQRT( SUMEPS / SUMSD )
                  ERRMAX = TWO*SQRT( SUMEPS * SUMSD )
               END IF

               IF (HCD .EQ. ZERO) HCD = TEN*HFD

               IF (MSGLVL .GT. 0) THEN
                  IF (HEADNG) WRITE (NOUT, 1100)
                  WRITE (NOUT, 1200) J, XJ, HFD, HCD, ERRMAX
                  HEADNG = .FALSE.
               END IF
            END IF

            FDNORM    = MAX (FDNORM, HFD)
            HFORWD(J) = HFD / (ONE + ABS(XJ))
            HCNTRL(J) = HCD / (ONE + ABS(XJ))
            X(J)      = XJ
  600    CONTINUE

         IF (NCCNST + NFCNST .GT. 0) THEN

*           Check that the constants have been set properly by
*           evaluating the gradients at a strange (but feasible) point.

            D      =   ONE / N

            DO 710, J = 1, N
               XJ     =   X(J)
               STEPBL = - ONE
               STEPBU =   ONE
               IF (BL(J) .GT. BIGLOW)
     $            STEPBL = MAX( STEPBL, BL(J) - XJ )
               IF (BU(J) .LT. BIGUPP  .AND.  BU(J) .GT. BL(J))
     $            STEPBU = MIN( STEPBU, BU(J) - XJ )

               IF (HALF*(STEPBL + STEPBU) .LT. ZERO) THEN
                  Y(J) = XJ + D*STEPBL
               ELSE
                  Y(J) = XJ + D*STEPBU
               END IF

               D = FACTOR*D
  710       CONTINUE

            IF (NCNLN .GT. 0) THEN
               CALL ILOAD ( NCNLN, (1), NEEDC, 1 )
               CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                      NEEDC, Y, C2, CJACU, NSTATE )
               IF (MODE .LT. 0) GO TO 9999
            END IF

            IF (M     .GT. 0) THEN
               CALL OBJFUN( MODE, M, N, LDFJU, Y, F2, FJACU, NSTATE )
               IF (MODE .LT. 0) GO TO 9999
            END IF

*           ------------------------------------------------------------
*           Loop over each of the components of  x.
*           ------------------------------------------------------------
            DO 800, J = 1, N
               YJ     = Y(J)
               DX     = HALF*(X(J) - YJ)
               Y(J)   = YJ + DX

               IF (NCNLN .GT. 0) THEN
                  NFOUND = 0
                  DO 720, I = 1, NCNLN
                     IF (CJACU(I,J) .EQ. - RDUMMY) THEN
                        NEEDC(I) = 1
                        NFOUND   = NFOUND + 1
                     ELSE
                        NEEDC(I) = 0
                     END IF
  720             CONTINUE

                  IF (NFOUND .GT. 0) THEN
                     CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                            NEEDC, Y, C1, CJACU, NSTATE )
                     IF (MODE .LT. 0) GO TO 9999

                     DO 730, I = 1, NCNLN
                        IF (NEEDC(I) .EQ. 1) THEN
                           CJDIFF = ( C1(I) -  C2(I) ) / DX
                           IF (CJDIFF .EQ. CJAC(I,J)) THEN
                              CJACU(I,J) = CJDIFF
                           ELSE
                              CJACU(I,J) = RDUMMY
                              NCCNST    = NCCNST - 1
                              NCDIFF    = NCDIFF + 1
                           END IF
                        END IF
  730                CONTINUE
                  END IF
               END IF

*              Now check the objective Jacobian.

               IF (M .GT. 0) THEN
                  NFOUND = 0
                  DO 740, I = 1, M
                     IF (FJACU(I,J) .EQ. - RDUMMY) THEN
                        NFOUND = NFOUND + 1
                     END IF
  740             CONTINUE

                  IF (NFOUND .GT. 0) THEN
                     CALL OBJFUN( MODE, M, N, LDFJU,
     $                            Y, F1, FJACU, NSTATE )
                     IF (MODE .LT. 0) GO TO 9999

                     DO 750, I = 1, M
                        IF (FJACU(I,J) .EQ. - RDUMMY) THEN
                           FJDIFF = ( F1(I) -  F2(I) ) / DX
                           IF (FJDIFF .EQ. FJAC(I,J)) THEN
                              FJACU(I,J) = FJDIFF
                           ELSE
                              FJACU(I,J) = RDUMMY
                              NFCNST     = NFCNST - 1
                              NFDIFF     = NFDIFF + 1
                           END IF
                        END IF
  750                CONTINUE
                  END IF
               END IF

               Y(J)  = YJ
  800       CONTINUE

            IF (MSGLVL .GT. 0) THEN
               IF (LVLDER .LT. 2  .AND.  NCCNST .GT. 0)
     $            WRITE (NOUT, 1300) NCCNST
               IF (LVLDER .NE. 1  .AND.  NFCNST .GT. 0)
     $            WRITE (NOUT, 1400) NFCNST
            END IF

            IF (NCDIFF .EQ. 0  .AND.  LVLDER .LT. 2) THEN
               IF (LVLDER .EQ. 0) LVLDER = 2
               IF (LVLDER .EQ. 1) LVLDER = 3
               IF (MSGLVL .GT. 0) WRITE (NOUT, 1500) LVLDER
            END IF

            IF (NFDIFF .EQ. 0  .AND.  LVLDER .NE. 1) THEN
               IF (LVLDER .EQ. 0) LVLDER = 1
               IF (LVLDER .EQ. 2) LVLDER = 3
               IF (MSGLVL .GT. 0) WRITE (NOUT, 1600) LVLDER
            END IF
         END IF
      ELSE IF (LFDSET .EQ. 2) THEN

*        The user has supplied HFORWD and HCNTRL.
*        Check for wild values.

         DO 900, J = 1, N
            IF (HFORWD(J) .LE. ZERO) THEN
               if(msglvl.ne.0) WRITE (NOUT, 2000) J, HFORWD(J), EPSPT5
               HFORWD(J) = EPSPT5
            END IF
  900    CONTINUE
         DO 910, J = 1, N
            IF (HCNTRL(J) .LE. ZERO) THEN
               if(msglvl.ne.0) WRITE (NOUT, 2100) J, HCNTRL(J), EPSPT3
               HCNTRL(J) = EPSPT3
            END IF
  910    CONTINUE
      END IF

      RETURN

 9999 INFORM = MODE
      RETURN

 1000 FORMAT(//' Computation of the finite-difference intervals'
     $       / ' ----------------------------------------------' )
 1100 FORMAT(//'    J      X(J)   Forward DX(J)   Central DX(J) ',
     $         '     Error est.' /)
 1200 FORMAT(  I5, 1P, E10.2, E16.6, 2E16.6 )
 1300 FORMAT(/ I5,  '  constant constraint gradient elements assigned.')
 1400 FORMAT(/ I5,  '  constant objective  gradient elements assigned.')
 1500 FORMAT(//' All missing constraint Jacobian elements',
     $         ' are constants.   Derivative level increased to ', I4 )
 1600 FORMAT(//' All missing objective  Jacobian elements',
     $         ' are constants.   Derivative level increased to ', I4 )
 2000 FORMAT(' XXX  ', I4,'-th difference interval ',         1PE10.2,
     $       ' replaced by ', 1PE10.2 )
 2100 FORMAT(' XXX  ', I4,'-th central-difference interval ', 1PE10.2,
     $       ' replaced by ', 1PE10.2 )

*     End of  CHFDLS.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        
      SUBROUTINE CHFGRD( INFORM, MSGLVL, N,
     $                   BIGBND, EPSRF, OKTOL, FDCHK, OBJF, XNORM,
     $                   OBJFUN,
     $                   BL, BU, GRAD, GRADU, DX, X, Y, W, LENW )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)

      DOUBLE PRECISION   BL(N), BU(N), GRAD(N), GRADU(N), DX(N)
      DOUBLE PRECISION   X(N), Y(N), W(LENW)
      EXTERNAL           OBJFUN

C***********************************************************************
C     CHFGRD  checks if the gradients of the objective function have
C     been coded correctly.
C
C     On input,  the value of the objective function at the point X is
C     stored in OBJF.  The corresponding gradient is stored in GRADU.
C     If any gradient component has not been specified,  it will have a
C     dummy value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods.  If this proves 
C     satisfactory and no further information is desired, CHFGRD is 
C     terminated.  Otherwise, the routine CHCORE is called to give 
C     optimal step-sizes and a forward-difference approximation to each 
C     component of the gradient for which a test is deemed necessary,
C     either by the program or the user.
C
C     Other inputs:
C
C           X         The n-dimensional point at which the
C                     gradient is to be verified.
C           EPSRF     The positive bound on the relative error
C                     associated with computing the function at
C                     the point x.
C           OKTOL     The desired relative accuracy which the
C                     components of the gradient should satisfy.
C
C     LVRFYC has the following meaning...
C
C       -1        do not perform any check.
C        0        do the cheap test only.
C        1 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  19-May-1985.
C     This version of CHFGRD dated  28-Jun-1989.  
C***********************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9

      COMMON    /SOL5NP/ LVRFYC, JVERFY(4)

      LOGICAL            NPDBG        
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG

      LOGICAL            CONST , DEBUG , DONE  , FIRST , HEADNG
      LOGICAL            NEEDED, OK
      CHARACTER*4        KEY   , LBAD  , LGOOD
      CHARACTER*18       RESULT(0:4)
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      EXTERNAL           DDOT
      PARAMETER         (RDUMMY =-11111.0              )
      PARAMETER         (ZERO   =0.0D+0, HALF  = 0.5D+0, POINT9 =0.9D+0)
      PARAMETER         (ONE    =1.0D+0, TWO   = 2.0D+0, TEN    =1.0D+1)
      PARAMETER         (LBAD   ='BAD?', LGOOD = '  OK')
      DATA               RESULT
     $                 / '                 ', 'Constant?      ',
     $                   'Linear or odd?   ', 'Too nonlinear?',
     $                   'Small derivative?'                   /

      INFORM = 0
      NEEDED = LVRFYC .EQ. 0  .OR.  LVRFYC .EQ. 1  .OR.  LVRFYC .EQ. 3
      IF (.NOT. NEEDED) RETURN
 
      if (msglvl.ne.0) WRITE (NOUT, 1000)
      DEBUG  = NPDBG  .AND.  INPDBG(5) .GT. 0
      NSTATE = 0

      BIGLOW = - BIGBND
      BIGUPP =   BIGBND

*     ==================================================================
*     Perform the cheap test.
*     ==================================================================
      H =     (ONE + XNORM)*FDCHK

      DXJ  = ONE / N
      DO 110 J = 1, N
         DX(J) =   DXJ
         DXJ   = - DXJ*POINT9
  110 CONTINUE

*     ------------------------------------------------------------------
*     Do not perturb X(J) if the  J-th  element is missing.
*     Compute the directional derivative.
*     ------------------------------------------------------------------
      NCHECK = 0
      DO 120 J = 1, N
         IF (GRAD(J) .EQ. RDUMMY) THEN
            DX(J) = ZERO
         ELSE
            NCHECK = NCHECK + 1

            XJ     =   X(J)                             
            STEPBL = - ONE       
            STEPBU =   ONE
            IF (BL(J) .GT. BIGLOW)
     $         STEPBL = MAX( STEPBL, BL(J) - XJ )
            IF (BU(J) .LT. BIGUPP  .AND.  BU(J) .GT. BL(J))
     $         STEPBU = MIN( STEPBU, BU(J) - XJ )

            IF (HALF*(STEPBL + STEPBU) .LT. ZERO) THEN
               DX(J) = DX(J)*STEPBL
            ELSE
               DX(J) = DX(J)*STEPBU
            END IF
         END IF
  120 CONTINUE

      IF (NCHECK .EQ. 0) THEN
         if (msglvl.ne.0) WRITE (NOUT, 3500)
         RETURN
      END IF
      GDX    = DDOT  ( N, GRADU, 1, DX, 1 )

*     ------------------------------------------------------------------
*     Make forward-difference approximation along  p.
*     ------------------------------------------------------------------
      CALL DCOPY ( N,     X, 1, Y, 1 )
      CALL DAXPY ( N, H, DX, 1, Y, 1 )
      
      MODE   = 0
      CALL OBJFUN( MODE, N, Y, OBJF1, GRADU, NSTATE )
      IF (MODE .LT. 0) GO TO 9999

      GDIFF =    (OBJF1 - OBJF) / H
      ERROR = ABS(GDIFF - GDX ) / (ABS(GDX) + ONE)

      OK    = ERROR .LE. OKTOL
      IF (OK) THEN
         IF (MSGLVL .GT. 0) WRITE (NOUT, 1100)
      ELSE
         if (msglvl.ne.0) WRITE (NOUT, 1200)
         IF (ERROR .GE. POINT9) INFORM = 1
      END IF

      IF (MSGLVL .GT. 0) WRITE (NOUT, 1300) GDX, GDIFF

*     ==================================================================
*     Component-wise check.
*     ==================================================================
      IF (LVRFYC .EQ. 1  .OR.  LVRFYC .EQ. 3) THEN
         HEADNG = .TRUE.
         ITMAX  = 3
         NWRONG = 0
         NGOOD  = 0
         JMAX   = 0
         EMAX   = ZERO
         NCHECK = 0
         J1     = JVERFY(1)
         J2     = JVERFY(2)

*        ---------------------------------------------------------------
*        Loop over each of the components of  x.
*        ---------------------------------------------------------------
         DO 500 J = J1, J2

            IF (GRAD(J) .NE. RDUMMY) THEN
*              ---------------------------------------------------------
*              Check this gradient component.
*              ---------------------------------------------------------
               NCHECK = NCHECK + 1
               GJ     = GRAD(J)
               GSIZE  = ABS( GJ )
               XJ     = X(J)
*              ---------------------------------------------------------
*              Find a finite-difference interval by iteration.
*              ---------------------------------------------------------
               ITER   = 0
               EPSA   = EPSRF*(ONE + ABS( OBJF ))
               CDEST  = ZERO
               SDEST  = ZERO
               FIRST  = .TRUE.

               STEPBL = BIGLOW
               STEPBU = BIGUPP
               IF (BL(J) .GT. BIGLOW) STEPBL = BL(J) - XJ
               IF (BU(J) .LT. BIGUPP) STEPBU = BU(J) - XJ

               HOPT   = TWO*(ONE + ABS( XJ ))*SQRT( EPSRF )
               H      = TEN*HOPT
               IF (HALF*(STEPBL + STEPBU) .LT. ZERO) H =  - H

*+             REPEAT
  400             X(J)  = XJ + H
                  CALL OBJFUN( MODE, N, X, F1, GRADU, NSTATE )
                  IF (MODE .LT. 0) GO TO 9999

                  X(J)  = XJ + H + H
                  CALL OBJFUN( MODE, N, X, F2, GRADU, NSTATE )
                  IF (MODE .LT. 0) GO TO 9999
                              
                  CALL CHCORE( DEBUG, DONE, FIRST, EPSA, EPSRF, OBJF,XJ,
     $                         INFO, ITER, ITMAX,
     $                         CDEST, FDEST, SDEST, ERRBND, F1,
     $                         F2, H, HOPT, HPHI )

*+             UNTIL     DONE
               IF (.NOT. DONE) GO TO 400

*              ---------------------------------------------------------
*              Exit for this variable.
*              ---------------------------------------------------------
               GDIFF = CDEST
               X(J)  = XJ

               ERROR = ABS(GDIFF - GJ) / (GSIZE + ONE)
               IF (ERROR .GE. EMAX) THEN
                  EMAX  = ERROR
                  JMAX  = J
               END IF

               OK =  ERROR .LE. OKTOL
               IF (OK) THEN
                  KEY    = LGOOD
                  NGOOD  = NGOOD  + 1
               ELSE
                  KEY    = LBAD
                  NWRONG = NWRONG + 1
               END IF

*              Zero components are not printed.

               CONST = OK .AND. INFO .EQ. 1 .AND. ABS(GJ) .LT. EPSPT8
               IF (.NOT. CONST) THEN
                  IF (HEADNG.and.msglvl.ne.0) WRITE (NOUT, 3000)
                  IF (OK) THEN
                     if (msglvl.ne.0) WRITE (NOUT, 3100) J, XJ, HOPT, 
     $                    GJ, GDIFF, KEY, ITER
                  ELSE
                     if (msglvl.ne.0) WRITE (NOUT, 3110) J, XJ, HOPT, 
     $                    GJ, GDIFF, KEY, ITER, RESULT(INFO)
                  END IF
                  HEADNG = .FALSE.
               END IF
            END IF
  500    CONTINUE

*        ===============================================================
*        Done.
*        ===============================================================
         INFORM = 0
         IF (NWRONG .EQ. 0) THEN
            if (msglvl.ne.0) WRITE (NOUT, 3200) NGOOD , NCHECK, J1, J2
         ELSE
            if (msglvl.ne.0) WRITE (NOUT, 3300) NWRONG, NCHECK, J1, J2
            IF (ERROR .GE. POINT9) INFORM = 1
         END IF
         if (msglvl.ne.0) WRITE (NOUT, 3400) EMAX, JMAX
      END IF

      CALL DCOPY ( N, GRAD, 1, GRADU, 1 )

      RETURN

 9999 INFORM = MODE
      RETURN

 1000 FORMAT(/// ' Verification of the objective gradients.'
     $       /   ' ----------------------------------------' )
 1100 FORMAT(/   ' The objective gradients seem to be ok.')
 1200 FORMAT(/   ' XXX  The objective gradients seem to be incorrect.')
 1300 FORMAT(/   ' Directional derivative of the objective', 1PE18.8/
     $           ' Difference approximation               ', 1PE18.8 )
 3000 FORMAT(// 4X, 'J', 4X, 'X(J)', 5X, 'DX(J)', 11X,
     $           'G(J)', 9X, '  Difference approxn  Itns' /)
 3100 FORMAT(  I5, 1P, 2E10.2,      2E18.8, 2X, A4, I6          )
 3110 FORMAT(  I5, 1P, 2E10.2,      2E18.8, 2X, A4, I6, 2X, A18 )
 3200 FORMAT(/ I7, '  Objective gradients out of the', I6,
     $             '  set in cols', I6, '  through', I6,
     $             '  seem to be ok.')
 3300 FORMAT(/   ' XXX  There seem to be', I6,
     $           '  incorrect objective gradients out of the', I6,
     $           '  set in cols', I6, '  through', I6 )
 3400 FORMAT(/   ' The largest relative error was', 1PE12.2,
     $           '   in element', I6 /)
 3500 FORMAT(/   ' No gradient elements assigned.' )

*     End of  CHFGRD.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE CHCJAC( INFORM, LVLDER, MSGLVL,
     $                   NCSET, N, NCNLN, LDCJ, LDCJU,
     $                   BIGBND, EPSRF, OKTOL, FDCHK, XNORM,
     $                   CONFUN, NEEDC,
     $                   BL, BU, C, C1, CJAC, CJACU, CJDX,
     $                   DX, ERR, X, Y, W, LENW )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            NEEDC(*)
      DOUBLE PRECISION   BL(N), BU(N), C(*), C1(*), CJDX(*),
     $                   CJAC(LDCJ,*), CJACU(LDCJU,*), ERR(*)
      DOUBLE PRECISION   DX(N), X(N), Y(N), W(LENW)
      EXTERNAL           CONFUN

C***********************************************************************
C     CHCJAC  checks if the gradients of the constraints have been coded
C     correctly.
C
C     On input,  the values of the constraints at the point X are stored
C     in C.  Their corresponding gradients are stored in CJACU.  If any
C     Jacobian component has not been specified,  it will have a dummy
C     value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods.  If this proves 
C     satisfactory and no further information is desired, CHCJAC is 
C     terminated.  Otherwise, CHCORE is called to give optimal stepsizes
C     and a central-difference approximation to each component of the 
C     Jacobian for which a test is deemed necessary, either by the 
C     program or the user.
C
C     LVRFYC has the following meaning...
C
C       -1        do not perform any check.
C        0        do the cheap test only.
C        2 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  19-May-1985.
C     This version of CHCJAC dated 28-Jun-1989.  
C***********************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9

      COMMON    /SOL5NP/ LVRFYC, JVERFY(4)

      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG

      LOGICAL            CONST , DEBUG , DONE  , FIRST , HEADNG
      LOGICAL            NEEDED, OK
      CHARACTER*4        KEY   , LBAD  , LGOOD
      CHARACTER*18       RESULT(0:4)
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      EXTERNAL           DDOT  , IDAMAX
      PARAMETER         (RDUMMY =-11111.0              )
      PARAMETER         (ZERO   =0.0D+0, HALF   =0.5D+0, POINT9 =0.9D+0)
      PARAMETER         (ONE    =1.0D+0, TWO    =2.0D+0, TEN    =1.0D+1)
      PARAMETER         (LBAD   ='BAD?', LGOOD  ='  OK')
      DATA               RESULT
     $                 / '                 ', 'Constant?      ',
     $                   'Linear or odd?   ', 'Too nonlinear?',
     $                   'Small derivative?'                   /

      INFORM = 0
      NEEDED = NCNLN  .GT. 0  .AND.
     $         LVRFYC .EQ. 0  .OR.   LVRFYC .EQ. 2  .OR.  LVRFYC .EQ. 3
      IF (.NOT. NEEDED) RETURN

      if (msglvl.ne.0) WRITE (NOUT, 1000)
      DEBUG  = NPDBG  .AND.  INPDBG(5) .GT. 0
      NSTATE = 0

      BIGLOW = - BIGBND
      BIGUPP =   BIGBND

*     ==================================================================
*     Perform the cheap test.
*     ==================================================================
      H = (ONE + XNORM)*FDCHK

      DXJ  = ONE / N
      DO 110 J = 1, N
         DX(J) =   DXJ
         DXJ   = - DXJ*POINT9
  110 CONTINUE

*     ------------------------------------------------------------------
*     Do not perturb  X(J)  if the  J-th  column contains any
*     unknown elements.  Compute the directional derivative for each
*     constraint gradient.
*     ------------------------------------------------------------------
      NCHECK = 0
      DO 140 J = 1, N
         DO 130 I = 1, NCNLN
            IF (CJAC(I,J) .EQ. RDUMMY) THEN
               DX(J) = ZERO      
               GO TO 140
            END IF           
  130    CONTINUE
         NCHECK = NCHECK + 1

         XJ     =   X(J)
         STEPBL = - ONE
         STEPBU =   ONE
         IF (BL(J) .GT. BIGLOW)
     $      STEPBL = MAX( STEPBL, BL(J) - XJ )
         IF (BU(J) .LT. BIGUPP  .AND.  BU(J) .GT. BL(J))
     $      STEPBU = MIN( STEPBU, BU(J) - XJ )

         IF (HALF*(STEPBL + STEPBU) .LT. ZERO) THEN
            DX(J) = DX(J)*STEPBL
         ELSE
            DX(J) = DX(J)*STEPBU
         END IF
  140 CONTINUE

      IF (NCHECK .EQ. 0) THEN
         if (msglvl.ne.0) WRITE (NOUT, 2300)
      ELSE

*        Compute  (Jacobian)*DX.

         CALL DGEMV ( 'Normal', NCNLN, N, ONE, CJACU, LDCJU,
     $                DX, 1, ZERO, CJDX, 1 )

*        ---------------------------------------------------------------
*        Make forward-difference approximation along DX.
*        ---------------------------------------------------------------
         CALL DCOPY ( N,     X, 1, Y, 1 )
         CALL DAXPY ( N, H, DX, 1, Y, 1 )

         CALL ILOAD ( NCNLN, (1), NEEDC, 1 )

         MODE   = 0
         CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                NEEDC, Y, C1, CJACU, NSTATE )
         IF (MODE .LT. 0) GO TO 9999

*        Set  ERR = (C1 - C)/H  - Jacobian*DX.  This should be small.

         DO 170 I = 1, NCNLN
            ERR(I) = (C1(I) - C(I)) / H  -  CJDX(I)
  170    CONTINUE            
         IMAX  = IDAMAX( NCNLN, ERR, 1 )
         EMAX  = ABS(ERR(IMAX)) / (ABS(CJDX(IMAX)) + ONE)

         IF (EMAX .LE. OKTOL) THEN
            IF (MSGLVL .GT. 0) WRITE (NOUT, 2000)
         ELSE
            !write(*,*) 'CHSUBS:',emax,oktol,point9,imax
            if (msglvl.ne.0) WRITE (NOUT, 2100)
            IF (EMAX .GE. POINT9) INFORM = 1
         END IF
         IF (MSGLVL .gt. 0) WRITE (NOUT, 2200) EMAX, IMAX
      END IF

*     ==================================================================
*     Component-wise check.
*     ==================================================================
      IF (LVRFYC .GE. 2) THEN
         IF (LVLDER .EQ. 3) THEN

*           Recompute the Jacobian to find the non-constant elements.

            CALL F06QHF( 'General', NCNLN, N, RDUMMY, RDUMMY, 
     $                   CJACU, LDCJU )

            CALL ILOAD ( NCNLN, (1), NEEDC, 1 )
            NSTATE = 0
            MODE   = 2

            CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                   NEEDC, X, C1, CJACU, NSTATE )
            IF (MODE .LT. 0) GO TO 9999

         END IF

         CALL ILOAD ( NCNLN, (0), NEEDC, 1 )

         ITMAX  =   3
         NCHECK =   0
         NWRONG =   0
         NGOOD  =   0
         COLMAX = - ONE
         JCOL   =   0
         IROW   =   0
         MODE   =   0
         J3     =   JVERFY(3)
         J4     =   JVERFY(4)

*        ---------------------------------------------------------------
*        Loop over each column.
*        ---------------------------------------------------------------
         DO 600 J = J3, J4

            CALL DLOAD ( NCNLN, ZERO, ERR, 1 )
            HEADNG = .TRUE.
            XJ     = X(J)

            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J) .GT. BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J) .LT. BIGUPP) STEPBU = BU(J) - XJ

            SIGNH  = ONE
            IF (HALF*(STEPBL + STEPBU) .LT. ZERO) SIGNH =  - ONE

            DO 500, I = 1, NCNLN
               EPSACI   = EPSRF*(ONE + ABS( C(I) ))

               IF (CJACU(I,J) .NE. RDUMMY) THEN
*                 ------------------------------------------------------
*                 Check this Jacobian element.
*                 ------------------------------------------------------
                  NCHECK   = NCHECK + 1
                  NEEDC(I) = 1

                  CIJ    = CJAC(I,J)
                  CJSIZE = ABS( CIJ )
*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  ITER   = 0
                  HOPT   = TWO*(ONE + ABS( XJ ))*SQRT( EPSRF )
                  H      = TEN*HOPT*SIGNH
                  CDEST  = ZERO
                  SDEST  = ZERO
                  FIRST  = .TRUE.

*+                REPEAT
  400                X(J)  = XJ + H
                     CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                            NEEDC, X, C1, CJACU, NSTATE )
                     IF (MODE .LT. 0) GO TO 9999
                     F1    = C1(I)

                     X(J)  = XJ + H + H
                     CALL CONFUN( MODE, NCNLN, N, LDCJU,
     $                            NEEDC, X, C1, CJACU, NSTATE )
                     IF (MODE .LT. 0) GO TO 9999
                     F2    = C1(I)
                                                     
                     CALL CHCORE( DEBUG,DONE,FIRST,EPSACI,EPSRF,C(I),XJ,
     $                            INFO, ITER, ITMAX,
     $                            CDEST, FDEST, SDEST, ERRBND, F1,
     $                            F2, H, HOPT, HPHI )

*+                UNTIL     DONE
                  IF (.NOT. DONE) GO TO 400

*                 ------------------------------------------------------
*                 Exit for this element.
*                 ------------------------------------------------------
                  CJDIFF   = CDEST
                  ERR(I)   = ABS(CJDIFF - CIJ) / (CJSIZE + ONE)

                  OK       = ERR(I) .LE. OKTOL
                  IF (OK) THEN
                     KEY    = LGOOD
                     NGOOD  = NGOOD  + 1
                  ELSE
                     KEY    = LBAD
                     NWRONG = NWRONG + 1
                  END IF

                  CONST = OK .AND. INFO       .EQ. 1
     $                       .AND. ABS( CIJ ) .LT. EPSPT8
                  IF (.NOT. CONST) THEN
                     IF (HEADNG) THEN
                        if (msglvl.ne.0) WRITE (NOUT, 4000)
                        IF (OK.and.msglvl.ne.0)
     $                     WRITE (NOUT, 4100)   J, XJ    , HOPT, I,
     $                                        CIJ, CJDIFF, KEY , ITER
                        IF (.NOT. OK.and.msglvl.ne.0)
     $                     WRITE (NOUT, 4110)   J, XJ    , HOPT, I,
     $                                        CIJ, CJDIFF, KEY , ITER,
     $                                        RESULT(INFO)
                        HEADNG = .FALSE.
                     ELSE
                        IF (OK.and.msglvl.ne.0)
     $                     WRITE (NOUT, 4200)              HOPT, I,
     $                                        CIJ, CJDIFF, KEY , ITER
                        IF (.NOT. OK.and.msglvl.ne.0)
     $                     WRITE (NOUT, 4210)              HOPT, I,
     $                                        CIJ, CJDIFF, KEY , ITER,
     $                                        RESULT(INFO)
                     END IF
                  END IF
                  NEEDC(I) = 0
               END IF
  500       CONTINUE

*           ------------------------------------------------------------
*           Finished with this column.
*           ------------------------------------------------------------
            IF (.NOT. HEADNG) THEN
               IMAX = IDAMAX( NCNLN, ERR, 1 )
               EMAX = ABS( ERR(IMAX) )

               IF (EMAX .GE. COLMAX) THEN
                  IROW   = IMAX
                  JCOL   = J
                  COLMAX = EMAX
               END IF
            END IF
            X(J) = XJ

  600    CONTINUE

         INFORM = 0
         IF (NCHECK .EQ. 0) THEN
            if (msglvl.ne.0) WRITE (NOUT, 4600) NCSET
         ELSE
            IF (NWRONG .EQ. 0) THEN
               if (msglvl.ne.0) WRITE (NOUT, 4300) NGOOD , NCHECK, J3, J4
            ELSE
               if (msglvl.ne.0) WRITE (NOUT, 4400) NWRONG, NCHECK, J3, J4
               IF (COLMAX .GE. POINT9) INFORM = 1
            END IF
            if (msglvl.ne.0) WRITE (NOUT, 4500) COLMAX, IROW, JCOL
         END IF

      END IF

*     Copy  ( constants + gradients + dummy values )  back into CJACU.
      
      CALL F06QFF( 'General', NCNLN, N, CJAC, LDCJ, CJACU, LDCJU )

      RETURN            

 9999 INFORM = MODE
      RETURN

 1000 FORMAT(/// ' Verification of the constraint gradients.'
     $       /   ' -----------------------------------------' )
 2000 FORMAT(/   ' The constraint Jacobian seems to be ok.')
 2100 FORMAT(/   ' XXX  The constraint Jacobian seems to be incorrect.')
 2200 FORMAT(/   ' The largest relative error was', 1PE12.2,
     $           '  in row', I5 /)
 2300 FORMAT(/   ' Every column contains a constant or',
     $           ' missing element.')
 4000 FORMAT(// ' Column    X(J)     DX(J)    Row   ',
     $          ' Jacobian Value      Difference Approxn  Itns'    )
 4100 FORMAT(/ I7,      1P, 2E10.2, I5, 2E18.8, 2X, A4, I6         )
 4110 FORMAT(/ I7,      1P, 2E10.2, I5, 2E18.8, 2X, A4, I6, 2X, A18)
 4200 FORMAT(  7X, 10X, 1P,  E10.2, I5, 2E18.8, 2X, A4, I6         )
 4210 FORMAT(  7X, 10X, 1P,  E10.2, I5, 2E18.8, 2X, A4, I6, 2X, A18)
 4300 FORMAT(/ I7, '  constraint Jacobian elements out of the', I6,
     $             '  set in cols', I6, '  through', I6,
     $             '  seem to be ok.')
 4400 FORMAT(/   ' XXX  There seem to be', I6,
     $           '  incorrect Jacobian elements out of the', I6,
     $           '  set in cols', I6, '  through', I6 )
 4500 FORMAT(/ ' The largest relative error was', 1PE12.2,
     $         '  in row', I5, ',  column', I5 /)
 4600 FORMAT(  ' All', I6, '   assigned Jacobian elements are',
     $         ' constant.' )

*     End of  CHCJAC.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE CHFJAC( INFORM, LVLDER, MSGLVL,
     $                   NFSET, M, N, LDFJ, LDFJU,
     $                   BIGBND, EPSRF, OKTOL, FDCHK, XNORM,
     $                   OBJFUN,
     $                   BL, BU, F, F1, FJAC, FJACU, FJDX,
     $                   DX, ERR, X, Y, W, LENW )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION   BL(N), BU(N), F(*), F1(*), FJDX(*),
     $                   FJAC(LDFJ,*), FJACU(LDFJU,*), ERR(*)
      DOUBLE PRECISION   DX(N), X(N), Y(N), W(LENW)
      EXTERNAL           OBJFUN

C***********************************************************************
C     CHFJAC  checks if the objective Jacobian matrix has been coded
C     correctly.
C
C     On input,  the values of the objective vector at the point X are
C     stored in F.  Their corresponding gradients are stored in FJACU.
C     If any Jacobian component has not been specified,  it will have a
C     dummy value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods.  If this proves
C     satisfactory and no further information is desired, CHFJAC is
C     terminated.  Otherwise, CHCORE is called to give optimal stepsizes
C     and a central-difference approximation to each component of the
C     Jacobian for which a test is deemed necessary, either by the
C     program or the user.
C
C     LVRFYC has the following meaning...
C
C       -1        do not perform any check.
C        0        do the cheap test only.
C        2 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  10-May-1988.
C     This version of CHFJAC dated  28-Jun-1989.
C***********************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9

      COMMON    /SOL5NP/ LVRFYC, JVERFY(4)

      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG

      LOGICAL            CONST , DEBUG , DONE  , FIRST , HEADNG
      LOGICAL            NEEDED, OK
      CHARACTER*4        KEY   , LBAD  , LGOOD
      CHARACTER*18       RESULT(0:4)
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      EXTERNAL           DDOT  , IDAMAX
      PARAMETER         (RDUMMY =-11111.0              )
      PARAMETER         (ZERO   =0.0D+0, HALF   =0.5D+0, POINT9 =0.9D+0)
      PARAMETER         (ONE    =1.0D+0, TWO    =2.0D+0, TEN    =1.0D+1)
      PARAMETER         (LBAD   ='BAD?', LGOOD  ='  OK')
      DATA               RESULT
     $                 / '                 ', 'Constant?      ',
     $                   'Linear or odd?   ', 'Too nonlinear?',
     $                   'Small derivative?'                   /

      INFORM = 0
      NEEDED = LVRFYC .EQ. 0  .OR.  LVRFYC .EQ. 1  .OR.  LVRFYC .EQ. 3
      IF (.NOT. NEEDED) RETURN

      if (msglvl.ne.0) WRITE (NOUT, 1000)
      DEBUG  = NPDBG  .AND.  INPDBG(5) .GT. 0
      NSTATE = 0

      BIGLOW = - BIGBND
      BIGUPP =   BIGBND

*     ==================================================================
*     Perform the cheap test.
*     ==================================================================
      H = (ONE + XNORM)*FDCHK

      DXJ  = ONE / N
      DO 110, J = 1, N
         DX(J) =   DXJ
         DXJ   = - DXJ*POINT9
  110 CONTINUE

*     ------------------------------------------------------------------
*     Do not perturb  X(J)  if the  J-th  column contains any
*     unknown elements.  Compute the directional derivative for each
*     objective gradient.
*     ------------------------------------------------------------------
      NCHECK = 0
      DO 140, J = 1, N
         DO 130, I = 1, M
            IF (FJAC(I,J) .EQ. RDUMMY) THEN
               DX(J) = ZERO
               GO TO 140
            END IF
  130    CONTINUE
         NCHECK = NCHECK + 1

         XJ     =   X(J)
         STEPBL = - ONE
         STEPBU =   ONE
         IF (BL(J) .GT. BIGLOW)
     $      STEPBL = MAX( STEPBL, BL(J) - XJ )
         IF (BU(J) .LT. BIGUPP  .AND.  BU(J) .GT. BL(J))
     $      STEPBU = MIN( STEPBU, BU(J) - XJ )

         IF (HALF*(STEPBL + STEPBU) .LT. ZERO) THEN
            DX(J) = DX(J)*STEPBL
         ELSE
            DX(J) = DX(J)*STEPBU
         END IF
  140 CONTINUE

      IF (NCHECK .EQ. 0) THEN
         if (msglvl.ne.0) WRITE (NOUT, 2300)
      ELSE

*        Compute  (Jacobian)*DX.

         CALL DGEMV ( 'Normal', M, N, ONE, FJACU, LDFJU,
     $                 DX, 1, ZERO, FJDX, 1 )

*        ---------------------------------------------------------------
*        Make forward-difference approximation along DX.
*        ---------------------------------------------------------------
         CALL DCOPY ( N,     X, 1, Y, 1 )
         CALL DAXPY ( N, H, DX, 1, Y, 1 )

         MODE   = 0
         CALL OBJFUN( MODE, M, N, LDFJU,
     $                Y, F1, FJACU, NSTATE )
         IF (MODE .LT. 0) GO TO 9999

*        Set  ERR = (F1 - F)/H  - Jacobian*DX.  This should be small.

         DO 170, I = 1, M
            ERR(I) = (F1(I) - F(I)) / H  -  FJDX(I)
  170    CONTINUE
         IMAX  = IDAMAX( M, ERR, 1 )
         EMAX  = ABS(ERR(IMAX)) / (ABS(FJDX(IMAX)) + ONE)

         IF (EMAX .LE. OKTOL) THEN
            IF (MSGLVL .GT. 0) WRITE (NOUT, 2000)
         ELSE
            if (msglvl.ne.0) WRITE (NOUT, 2100)
            IF (EMAX .GE. POINT9) INFORM = 1
         END IF
         IF (MSGLVL .gt. 0) WRITE (NOUT, 2200) EMAX, IMAX
      END IF

*     ==================================================================
*     Component-wise check.
*     ==================================================================
      IF (LVRFYC .GE. 2) THEN
         IF (LVLDER .EQ. 3) THEN

*           Recompute the Jacobian to find the non-constant elements.

            CALL F06QHF( 'General', M, N, RDUMMY, RDUMMY,
     $                   FJACU, LDFJU )

            NSTATE = 0
            MODE   = 2

            CALL OBJFUN( MODE, M, N, LDFJU,
     $                   X, F1, FJACU, NSTATE )
            IF (MODE .LT. 0) GO TO 9999

         END IF

         ITMAX  =   3
         NCHECK =   0
         NWRONG =   0
         NGOOD  =   0
         COLMAX = - ONE
         JCOL   =   0
         IROW   =   0
         MODE   =   0
         J3     =   JVERFY(3)
         J4     =   JVERFY(4)

*        ---------------------------------------------------------------
*        Loop over each column.
*        ---------------------------------------------------------------
         DO 600, J = J3, J4

            CALL DLOAD ( M, ZERO, ERR, 1 )
            HEADNG = .TRUE.
            XJ     = X(J)

            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J) .GT. BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J) .LT. BIGUPP) STEPBU = BU(J) - XJ   
                                                 
            SIGNH  = ONE
            IF (HALF*(STEPBL + STEPBU) .LT. ZERO) SIGNH =  - ONE

            DO 500, I = 1, M
               EPSAFI   = EPSRF*(ONE + ABS( F(I) ))

               IF (FJACU(I,J) .NE. RDUMMY) THEN
*                 ------------------------------------------------------
*                 Check this Jacobian element.
*                 ------------------------------------------------------
                  NCHECK   = NCHECK + 1

                  FIJ    = FJAC(I,J)
                  FJSIZE = ABS( FIJ )
*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  ITER   = 0
                  HOPT   = TWO*(ONE + ABS( XJ ))*SQRT( EPSRF )
                  H      = TEN*HOPT*SIGNH
                  CDEST  = ZERO
                  SDEST  = ZERO
                  FIRST  = .TRUE.

*+                REPEAT
  400                X(J)  = XJ + H
                     CALL OBJFUN( MODE, M, N, LDFJU,
     $                            X, F1, FJACU, NSTATE )
                     IF (MODE .LT. 0) GO TO 9999
                     FFORW = F1(I)

                     X(J)  = XJ + H + H
                     CALL OBJFUN( MODE, M, N, LDFJU,
     $                            X, F1, FJACU, NSTATE )
                     IF (MODE .LT. 0) GO TO 9999
                     FBACK = F1(I)

                     CALL CHCORE( DEBUG,DONE,FIRST,EPSAFI,EPSRF,F(I),XJ,
     $                            INFO, ITER, ITMAX,
     $                            CDEST, FDEST, SDEST, ERRBND, FFORW,
     $                            FBACK, H, HOPT, HPHI )

*+                UNTIL     DONE
                  IF (.NOT. DONE) GO TO 400

*                 ------------------------------------------------------
*                 Exit for this element.
*                 ------------------------------------------------------
                  FJDIFF   = CDEST
                  ERR(I)   = ABS(FJDIFF - FIJ) / (FJSIZE + ONE)

                  OK       = ERR(I) .LE. OKTOL
                  IF (OK) THEN
                     KEY    = LGOOD
                     NGOOD  = NGOOD  + 1
                  ELSE
                     KEY    = LBAD
                     NWRONG = NWRONG + 1
                  END IF

                  CONST = OK .AND. INFO       .EQ. 1
     $                       .AND. ABS( FIJ ) .LT. EPSPT8
                  IF (.NOT. CONST) THEN
                     IF (HEADNG) THEN
                        if (msglvl.ne.0) WRITE (NOUT, 4000)
                        IF (OK.and.msglvl.ne.0) WRITE (NOUT, 4100)   J, XJ    , HOPT, I,
     $                       FIJ, FJDIFF, KEY , ITER
                        IF (.NOT. OK.and.msglvl.ne.0) WRITE (NOUT, 4110)   J, XJ    , HOPT, I,
     $                                        FIJ, FJDIFF, KEY , ITER,
     $                                        RESULT(INFO)
                        HEADNG = .FALSE.
                     ELSE
                        IF (OK.and.msglvl.ne.0)
     $                     WRITE (NOUT, 4200)              HOPT, I,
     $                                        FIJ, FJDIFF, KEY , ITER
                        IF (.NOT. OK.and.msglvl.ne.0)
     $                     WRITE (NOUT, 4210)              HOPT, I,
     $                                        FIJ, FJDIFF, KEY , ITER,
     $                                        RESULT(INFO)
                     END IF
                  END IF
               END IF
  500       CONTINUE

*           ------------------------------------------------------------
*           Finished with this column.
*           ------------------------------------------------------------
            IF (.NOT. HEADNG) THEN
               IMAX = IDAMAX( M, ERR, 1 )
               EMAX = ABS( ERR(IMAX) )

               IF (EMAX .GE. COLMAX) THEN
                  IROW   = IMAX
                  JCOL   = J
                  COLMAX = EMAX
               END IF
            END IF
            X(J) = XJ

  600    CONTINUE

         INFORM = 0
         IF (NCHECK .EQ. 0) THEN
            if(msglvl.ne.0) WRITE (NOUT, 4600) NFSET
         ELSE
            IF (NWRONG .EQ. 0) THEN
               if(msglvl.ne.0) WRITE (NOUT, 4300) NGOOD , NCHECK, J3, J4
            ELSE
               if(msglvl.ne.0) WRITE (NOUT, 4400) NWRONG, NCHECK, J3, J4
               IF (COLMAX .GE. POINT9) INFORM = 1
            END IF
            if(msglvl.ne.0) WRITE (NOUT, 4500) COLMAX, IROW, JCOL
         END IF

      END IF

*     Copy  ( constants + gradients + dummy values )  back into FJACU.

      CALL F06QFF( 'General', M, N, FJAC, LDFJ, FJACU, LDFJU )

      RETURN

 9999 INFORM = MODE
      RETURN

 1000 FORMAT(/// ' Verification of the objective gradients.'
     $       /   ' ----------------------------------------' )
 2000 FORMAT(/   ' The objective Jacobian seems to be ok.')
 2100 FORMAT(/   ' XXX  The objective Jacobian seems to be incorrect.')
 2200 FORMAT(/   ' The largest relative error was', 1PE12.2,
     $           '  in row', I5 /)
 2300 FORMAT(/   ' Every column contains a constant or',
     $           ' missing element.')
 4000 FORMAT(// ' Column    X(J)     DX(J)    Row   ',
     $          ' Jacobian Value      Difference Approxn  Itns'    )
 4100 FORMAT(/ I7,      1P, 2E10.2, I5, 2E18.8, 2X, A4, I6         )
 4110 FORMAT(/ I7,      1P, 2E10.2, I5, 2E18.8, 2X, A4, I6, 2X, A18)
 4200 FORMAT(  7X, 10X, 1P,  E10.2, I5, 2E18.8, 2X, A4, I6         )
 4210 FORMAT(  7X, 10X, 1P,  E10.2, I5, 2E18.8, 2X, A4, I6, 2X, A18)
 4300 FORMAT(/ I7, '  Objective Jacobian elements out of the', I6,
     $             '  set in cols', I6, '  through', I6,
     $             '  seem to be ok.')
 4400 FORMAT(/   ' XXX  There seem to be', I6,
     $           '  incorrect objective Jacobian elements out of the',
     $             I6, '  set in cols', I6, '  through', I6 )
 4500 FORMAT(/ ' The largest relative error was', 1PE12.2,
     $         '  in row', I5, ',  column', I5 /)
 4600 FORMAT(  ' All', I6, '   assigned Jacobian elements are',
     $         ' constant.' )

*     End of  CHFJAC.

      END

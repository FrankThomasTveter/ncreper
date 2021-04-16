      SUBROUTINE GET_M (
     +     NRCHUSED,CHINDEX,CHINREL,CHINRELI,
     &     DBTODBTO,MM,ME,IRC)
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
#include "maxch.h"
#include "maxcat.h"
C     
      INTEGER NRCHUSED        ! NUMBER OF CHANNELS USED IN 1D-VAR
      INTEGER CHINDEX(MAXCH) ! INDEX OVER THE CHANNELS IN 2A-FILE
      INTEGER CHINREL(MAXCH) ! INDEX OVER THE CHANNELS USED IN 1D-VAR (CHINDEX SUBSET)
      INTEGER CHINRELI(MAXCH)
C
      REAL DBTODBTO(MAXCH,MAXCH)
      REAL MM(MAXCH,MAXCH)    ! eigenvectors
      REAL ME(MAXCH)          !eigenvalues
      INTEGER IRC
C
      CHARACTER*16 MYNAME
      DATA MYNAME /'GET_M'/
C
      INTEGER II,JJ
C
      INTEGER INFO,LWORK,N
      PARAMETER (LWORK=299)
      DOUBLE PRECISION EV(MAXCH),
     &     BUFF(NRCHUSED,NRCHUSED),WORK(LWORK)
C
      DO II=1,NRCHUSED
         if (CHINRELI(II).gt.maxch.or.CHINRELI(II).lt.1) then
            write(*,*) myname,'Invalid chinreli(',ii,')=',CHINRELI(II)
            irc=842
         end if
      end do
C
      if (irc.ne.0) return
C
      DO II=1,NRCHUSED
         DO JJ=1,NRCHUSED
            BUFF(II,JJ)=
     &           DBTODBTO(CHINRELI(II),CHINRELI(JJ))
         END DO
      END DO
C
C     GET EIGENVALUES AND VECTOR'S
C
      N=NRCHUSED
      CALL SSYEV('V','U',N,BUFF,N,
     &     EV,WORK,LWORK,INFO)
C
      IRC=INFO
      IF (IRC.NE.0) THEN
         WRITE(*,*) MYNAME,'Error return from DSYEV.',IRC
         RETURN
      END IF
C
      DO II=1,NRCHUSED
         ME(II)=EV(II)
         DO JJ=1,NRCHUSED
            MM(II,JJ)=BUFF(II,JJ)
         END DO
      END DO
C
      RETURN
      END




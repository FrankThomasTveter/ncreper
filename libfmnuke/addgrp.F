      SUBROUTINE ADDGRP(MXGRP,BHDRGRP,
     &     MAXGROUPS,NRGROUPS,DATAGRP,MXLVL,
     &     MXIND,NRIND,INDX,
     &     IRC)
C
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
      SAVE
C     HDRGRP(X,GRP)  CONTAINS GROUP NUMBER (1->) AND GROUP TYPE (&,%,$)
C     SUBGRP(GRP) CONTAINS FIRST SUB GROUP
C     INDX(1,SUBGRP(GRP)) CONTAINS INDEX TO NEXT SUB GROUP
C     INDX(2,SUBGRP(GRP)) CONTAINS SUB GROUP GRP
C     SUBGRP(INDX(2,SUBGRP(GRP))) CONTAINS GRP TO SUB-SUB GROUP...
C
      INTEGER MXGRP
      INTEGER NRGRP
      INTEGER PGRP
      INTEGER CGRP
      INTEGER BHDRGRP(2,0:MXGRP)
      INTEGER MAXGROUPS
      INTEGER NRGROUPS
      INTEGER DATAGRP(10,MAXGROUPS) ! hdr(2),x,find,lvl
      INTEGER MXLVL
      INTEGER LVLGRP(MAXGROUPS)
      INTEGER MXIND
      INTEGER NRIND
      INTEGER INDX(10,MXIND)
      INTEGER IRC
C
      CHARACTER*16 MYNAME
      DATA MYNAME /'ADDGRP'/
C
      logical first,bdone
      data first /.true./
      integer lvl, cind(0:maxgroups),curind,
     &     clev,curlev
      
C
      IF (NRGROUPS.EQ.0) THEN ! first is always super group
         PGRP=0
         LVL=-1
         MXLVL=-1
      ELSE                      ! find parent group position
         CALL FINDPOS(PGRP,LVL,
     &        MXGRP,BHDRGRP,
     &        MAXGROUPS,NRGROUPS,DATAGRP,
     &        MXIND,NRIND,INDX,
     &        IRC)
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from FINDPOS.',IRC
            RETURN
         END IF
      END IF
C
      LVL=LVL+1                 ! add the next levels
      BDONE=(LVL.GT.MXGRP)
      IF (.NOT.BDONE) BDONE=(BHDRGRP(2,LVL).EQ.0)
      DO WHILE (.NOT.BDONE)
         MXLVL=MAX(MXLVL,LVL)
         CALL ADDGROUP(PGRP,CGRP,LVL,
     &        MXGRP,BHDRGRP,
     &        MAXGROUPS,NRGROUPS,DATAGRP,
     &        MXIND,NRIND,INDX,
     &        IRC)
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from ADDGROUP.',IRC
            RETURN
         END IF
         PGRP=CGRP
         LVL=LVL+1              ! goto next level
         BDONE=(LVL.GT.MXGRP)
         IF (.NOT.BDONE) BDONE=(BHDRGRP(2,LVL).EQ.0)
      END DO
C     
      RETURN
C
      END

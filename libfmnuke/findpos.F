      subroutine findpos (PGRP,PLVL,
     &        MXGRP,BHDRGRP,
     &        MAXGROUPS,NRGROUPS,DATAGRP,
     &        MXIND,NRIND,INDX,
     &        IRC)
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
      implicit none
      save
C
      INTEGER PGRP            ! parent grp-position
      INTEGER PLVL             ! level of parent
      integer MXGRP           ! max levels of groups
      integer BHDRGRP(2,0:mxgrp) ! target grp identification
      INTEGER MAXGROUPS       ! MAX NUMBER GROUPS
      INTEGER NRGROUPS        ! NUMBER OF GROUPS
      INTEGER DATAGRP(10,MAXGROUPS) ! DATA FOR GROUPS
      INTEGER MXIND           ! MAX NUMBER OF INDEXES
      INTEGER NRIND           ! NUMBER OF INDEXES
      INTEGER INDX(10,MXIND)  ! INDEX ARRAY, 1=PARENT, 2=CHILD, 3=NEXT INDEX
      INTEGER IRC
C
      INTEGER CIND(0:MXGRP),POS,LVL
      LOGICAL BDONE
      CHARACTER*16 MYNAME
      DATA MYNAME /'FINDPOS'/
C     
      PGRP=1
      PLVL=0
      LVL=DATAGRP(5,1)          ! START AT TOP LEVEL
      CIND(LVL)=DATAGRP(4,1)    ! START AT TOP LEVEL
      BDONE=(LVL.GT.MXGRP)
      DO WHILE (.NOT. BDONE)
C     WRITE(*,*) MYNAME,'LVL, CIND:',LVL,CIND(LVL),
C     &        INDX(3,MAX(1,CIND(LVL)))
         IF (CIND(LVL).EQ.0) THEN
            IF (LVL.GT.0) THEN
               CIND(LVL)=0
               LVL=LVL-1
               CIND(LVL)=INDX(3,CIND(LVL))
            ELSE
               BDONE=.TRUE.
            END IF
         ELSE
            POS=INDX(2,CIND(LVL))
            IF (BHDRGRP(1,LVL+1).EQ. DATAGRP(1,POS).AND.
     &           BHDRGRP(2,LVL+1).EQ. DATAGRP(2,POS)) THEN ! SAME GRP
               PGRP=POS
               PLVL=LVL+1
               LVL=LVL+1
               CIND(LVL)=DATAGRP(4,POS)
            ELSE
               CIND(LVL)=INDX(3,CIND(LVL))
            END IF
         END IF
         IF (.NOT.BDONE) BDONE=(LVL.GT.MXGRP)
      END DO
C
      RETURN
      END

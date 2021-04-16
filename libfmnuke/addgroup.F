      subroutine addgroup (PGRP,CGRP,CLVL,
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
      INTEGER CGRP            ! child grp-position
      INTEGER CLVL            ! level of child
      integer MXGRP           ! max levels of groups
      integer BHDRGRP(2,0:mxgrp) ! target grp identification
      integer MAXGROUPS       ! max number groups
      integer NRGROUPS        ! number of groups
      integer DATAGRP(10,maxgroups) ! data for groups
      integer MXIND           ! max number of indexes
      integer NRIND           ! number of indexes
      integer INDX(10,mxind)  ! index array, 1=parent, 2=child, 3=next index
      integer IRC
C
      INTEGER CURIND,BIND
      CHARACTER*16 MYNAME
      DATA MYNAME /'ADDGROUP'/
C
      NRGROUPS=NRGROUPS+1
      IF (NRGROUPS.GT.MAXGROUPS) THEN
         WRITE(*,*) MYNAME,'Too many sub groups.',NRGROUPS
         IRC=945
         RETURN
      END IF
C     
      NRIND=NRIND+1
      IF (NRIND.GT.MXIND) THEN
         WRITE(*,*) MYNAME,'Too many group refs.'
         IRC=956
         RETURN
      END IF
C     
      DATAGRP(1,NRGROUPS)=BHDRGRP(1,CLVL)
      DATAGRP(2,NRGROUPS)=BHDRGRP(2,CLVL)
      DATAGRP(4,NRGROUPS)=0     ! no sub groups yet
      DATAGRP(6,NRGROUPS)=CLVL
      CGRP=NRGROUPS
C
      indx(1,nrind)=PGRP
      indx(2,nrind)=CGRP
      indx(3,nrind)=0
C     
      IF (PGRP.GT.0) THEN
         CURIND=DATAGRP(4,PGRP) ! first sub index
         IF (CURIND.EQ.0) THEN
            DATAGRP(4,PGRP)=NRIND
         ELSE
            BIND=CURIND
            DO WHILE (BIND.NE.0)
               CURIND=BIND
               BIND=INDX(3,BIND)
            END DO
            INDX(3,CURIND)=NRIND
         END IF
      END IF
C     
      RETURN
      END

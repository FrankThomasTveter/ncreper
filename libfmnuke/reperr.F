      SUBROUTINE REPERR (EGRP,CERR50,MXHDR,HDR250,
     &        MXGRP,NRGRP,HDRGRP,LFLGRP,
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
      IMPLICIT NONE
      SAVE
C
      INTEGER EGRP            ! group with error
      CHARACTER*50 CERR50       ! error message
      INTEGER MXHDR           ! max number of headers
      CHARACTER*250 HDR250(MXHDR) ! list of headers
      INTEGER MXGRP           ! MAX GROUPS
      INTEGER NRGRP           ! NUMBER OF GROUPS
      INTEGER HDRGRP(6,MXGRP) ! 1=parent, 2=next sibling, 3=first child,  4=type, 5=id/chdr, 6=lvl
      LOGICAL LFLGRP(MXGRP) ! is group set?
      INTEGER IRC
C
      CHARACTER*16 MYNAME
      DATA MYNAME /'REPERR'/
C
      INTEGER CIND(MXGRP),POS,LVL,TIND
      LOGICAL BDONE,BBDONE
      CHARACTER*250 PRE,PREF
      INTEGER LENH,LENP,LENC,LENGTH
      EXTERNAL LENGTH
C
      DATA PRE /'                                            '/
C     
      call chop(cerr50,50)
      LENC=LENGTH(CERR50,50,10)
      if (lenc.ne.0) WRITE(*,*) MYNAME(1:7),
     &     '======================================'
      LVL=1                     ! LEVEL (MUST HERE BE 1)
      CIND(LVL)=1               ! START WITH FIRST TOP LEVEL GROUP
      BDONE=.FALSE.
      DO WHILE (.NOT. BDONE)    ! LOOP OVER SUB GRP'S
         IF (CIND(LVL).EQ.0) THEN ! end of a chain reached
            IF (LVL.GT.1) THEN  ! move up
               CIND(LVL)=0
               LVL=LVL-1
               CIND(LVL)=HDRGRP(2,CIND(LVL)) ! move to next sibling
            ELSE
               BDONE=.TRUE.
            END IF
         ELSE
            POS=CIND(LVL)
            IF (POS.EQ.EGRP) THEN
               PREF=CERR50(1:LENC)//' ============'//
     &              '====================================='
               CALL CHOP(PREF,250)
               LENP=LENC+2*LVL
               PREF(LENP:LENP)='>'
            ELSE
               LENP=LENC+2*LVL
               PREF=PRE(1:LENP)
            END IF
            IF (LFLGRP(POS)) THEN
               PREF(LENP-1:LENP-1) = '*'
            END IF
            IF (HDRGRP(4,POS).EQ.1) THEN
               write (*,*) PREF(1:LENP)//
     &              '[all of the following:]'
            ELSE IF (HDRGRP(4,POS).EQ.2) THEN
               write (*,*) PREF(1:LENP)//
     &              '[up to one of the following:]'
            ELSE IF (HDRGRP(4,POS).EQ.3) THEN
               write (*,*) PREF(1:LENP)//'[possibly:]'
            ELSE IF (HDRGRP(4,POS).EQ.-1) THEN
               lenh=length(hdr250(-hdrgrp(5,POS)),250,10)
               write (*,*) PREF(1:LENP)//
     &              hdr250(-hdrgrp(5,POS))(1:LENH)
            END IF
            IF (HDRGRP(3,POS).NE.0) THEN
               LVL=LVL+1
               CIND(LVL)=HDRGRP(3,POS) ! first child
            ELSE
               CIND(LVL)=HDRGRP(2,CIND(LVL)) ! next sibling
            END IF
         END IF
      END DO
      if (lenc.ne.0) WRITE(*,*) MYNAME(1:7),
     &     '======================================'
C     
      return
      end

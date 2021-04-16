      SUBROUTINE STREP(STR250,NN,SRC100,REP100,LREP,IRC)
      IMPLICIT NONE
      SAVE
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
C     search for, and replace, key words in string with
C     given information...
C
      CHARACTER*250 STR250
      integer NN
      CHARACTER*100 SRC100(NN)
      CHARACTER*100 REP100(NN)
      logical lrep(nn)
      integer IRC
C
      logical BDONE
      CHARACTER*250 BUFF250
      CHARACTER*1000 BUFF1000
      integer II,JJ,LENGTH,LENS,LENR,LENI
C
      CHARACTER*16 MYNAME
      DATA MYNAME /'STREP'/
C
      DO II=1,NN
         CALL CHOP(SRC100(II),100)
         CALL CHOP(REP100(II),100)
         LENI=LENGTH(STR250,250,10)
         LENS=LENGTH(SRC100(II),100,2)
         LENR=LENGTH(REP100(II),100,2)
         if (lens.eq.0) then
            write(*,*) myname,'Invalid SRC at ',ii
            irc=922
            return
         end if
         BDONE=(LENI.EQ.0)
         DO WHILE (.NOT.BDONE)
            BUFF250=STR250
            BDONE=.TRUE.
            DO JJ=1,LENI-LENS+1
C     write(*,*) myname,'"',STR250(JJ:JJ+LENS-1),'"',
C     &              SRC100(II)(1:LENS),'"'
               IF (STR250(JJ:JJ+LENS-1).EQ.
     &              SRC100(II)(1:LENS)) THEN
                  BUFF1000=STR250(1:JJ-1)//
     &                 REP100(II)(1:LENR)//
     &                 STR250(JJ+LENS:LENI)
                  CALL CHOP(BUFF1000,251)
                  STR250=BUFF1000(1:250)
                  LENI=LENGTH(STR250,250,10)
C     write(*,*) myname,str250(1:leni)
                  lrep(ii)=.true.
                  BDONE=(BUFF250.eq.STR250)
               END IF
            END DO
         END DO
      END DO
C     
      RETURN
      END
C
C
C

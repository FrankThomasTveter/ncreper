      SUBROUTINE ARGUMENT(STRING,N,CPOS,BUFF,QQ,RR)
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
C     FINDS UP TO TWO ARGUMENTS OF A FUNCTION
C
      IMPLICIT NONE
      SAVE
      INTEGER N
      CHARACTER*(*) STRING
C
      CHARACTER*250 BUFF,BUFF1,BUFF2,BUFF3,NUKEHEAD
C
      INTEGER LENB1,LENB2,LENB3,LENGTH
      EXTERNAL LENGTH,NUKEHEAD
      LOGICAL QQ(2),QQ1,QQ2,ISINT,ISREAL
      EXTERNAL ISINT,ISREAL
      INTEGER CPOS,POS2(2),CODE,
     &     RI1,RI2
      REAL RR(2),RR1,RR2
C
      CHARACTER*8 MYNAME
      DATA MYNAME /'ARGUMENT'/
      INTEGER IRC
C
      LOGICAL ACTIVE,BDEB
      DATA ACTIVE /.FALSE./
C
      IF (.NOT.ACTIVE) CALL DEBUG(MYNAME,BDEB,ACTIVE)
C
      IF (BDEB) WRITE(*,*) MYNAME,'Debug: Routine starts.',IRC
C
C     FIND NEXT OBJECT
C
      CODE=+1
      CALL OBJECT(STRING,N,CPOS,POS2,CODE)
      CALL RESET(BUFF1,250)
      CALL CEQUAL(BUFF,250,STRING,POS2(1),POS2(2))
      CALL CHOP(BUFF,250)
      BUFF1=BUFF
      LENB1=LENGTH(BUFF1,250,LENB1)
C
C     INITIALISE
C
      QQ1=.FALSE.
      QQ2=.FALSE.
      RR1=0.0D0
      RR2=0.0D0
C
C     ANALYSE OBJECT
C
      IF (QQ(1)) THEN
C     FIND FIRST NUMERICAL ARGUMENT
         BUFF2=NUKEHEAD(BUFF1,250)
         BUFF3=BUFF1
         LENB2=LENGTH(BUFF2,250,LENB2)
         LENB3=LENGTH(BUFF3,250,LENB3)
         QQ1=ISREAL(BUFF2,RR1)
         IF (.NOT.QQ1) THEN
            QQ1=ISINT(BUFF2,RI1)
            IF (QQ1) RR1=DFLOAT(RI1)
         ENDIF
C
         IF (BDEB) WRITE(*,*) MYNAME,'Debug:',QQ1,BUFF2,RR1
C
         IF (QQ1.AND.QQ(2)) THEN
C     FIND SECOND ARGUMENT
            LENB3=LENGTH(BUFF3,250,LENB3)
            QQ2=ISREAL(BUFF3,RR2)
            IF (.NOT.QQ2) THEN
               QQ2=ISINT(BUFF3,RI2)
               IF (QQ2) RR2=DFLOAT(RI2)
            ENDIF
C
            IF (BDEB) WRITE(*,*) MYNAME,'Debug:',QQ2,BUFF3,RR2
C
         ELSEIF (QQ1) THEN
            QQ1=(LENB3.EQ.0)
         ENDIF
      ENDIF
C
C     ASSIGN NEW VALUES
C
      QQ(1)=QQ1
      QQ(2)=QQ2
      RR(1)=RR1
      RR(2)=RR2
C
      IF (BDEB) WRITE(*,*) MYNAME,'Debug: Routine ends.',IRC
C
      RETURN
      END

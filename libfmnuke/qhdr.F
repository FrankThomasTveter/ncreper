      SUBROUTINE QHDR(MAXUNIT,FILES,LUNIT,CUNIT,
     &     MAXHDR,NRHDR,HDR250,ACTHDR,NRLEN,LFLDAT,
     &     INTOUT,CNT,CODHDR,LINE,
     &     BUFF1000,LENB,
     &     PROCESS,BDONE,IRC)
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
      INTEGER MAXUNIT,CUNIT,LUNIT(MAXUNIT)
      CHARACTER*250 FILES(MAXUNIT)
      INTEGER MAXHDR,NRHDR
      CHARACTER*250 HDR250(MAXHDR)
      LOGICAL ACTHDR(MAXHDR)
      INTEGER NRLEN(MAXHDR)
      LOGICAL LFLDAT(MAXHDR)
      INTEGER INTOUT,CNT,CODHDR(MAXHDR),LINE
      CHARACTER*1000 BUFF1000
      INTEGER LENB
      LOGICAL PROCESS,BDONE
      INTEGER IRC
C
      INTEGER JJ,KK,LENF,LEND,LENGTH
      EXTERNAL LENGTH
      LOGICAL FOUND,ERROR
      DATA ERROR/.FALSE./
      CHARACTER*1000 DUMP
      CHARACTER*8 MYNAME
      DATA MYNAME /'QHDR'/
C
      JJ=1
      FOUND=.FALSE.
      DO WHILE(.NOT.FOUND .AND. JJ.LE.MAXHDR)
         FOUND=((NRLEN(JJ).LE.LENB).AND.(ACTHDR(JJ)))
         KK=1
         DO WHILE(KK.LE.NRLEN(JJ).AND.FOUND)
            IF(FOUND) 
     &           FOUND=(BUFF1000(KK:KK).EQ.HDR250(JJ)(KK:KK))
            KK=KK+1
         ENDDO
         IF (FOUND)THEN
            PROCESS=.FALSE.
C     
C     CHECK IF WE ARE INTERUPTING ANOTHER DATA BODY
C     
            IF (LINE.NE.-1) THEN
               IF (CODHDR(LINE).NE.-1 .AND. 
     &              CNT.NE.CODHDR(LINE)) THEN
                  ERROR=.TRUE.
                  BUFF1000='header in data body'
               END IF
            ENDIF
C     
C     PREPARE FOR DATABODY
C     
            LINE=JJ
            CNT=0
            INTOUT=LINE
            IF (CODHDR(LINE).EQ.0) THEN
               BDONE=.TRUE.
            ENDIF
         ELSE
            JJ=JJ+1
         ENDIF
      ENDDO
C
      IF (ERROR) THEN
         CALL CHOP(FILES(CUNIT),250)
         LENF=LENGTH(FILES(CUNIT),250,10)
         CALL CHOP(BUFF1000,250)
         LENB=LENGTH(BUFF1000,250,10)
         IF (LENF.GT.80)LENF=80
         WRITE(DUMP,*) 
     &        '"'//FILES(CUNIT)(1:LENF)//'", line ',
     &        LUNIT(CUNIT),': Error: '//BUFF1000(1:LENB)
         CALL CHOP(DUMP,200)
         LEND=LENGTH(DUMP,200,20)
         WRITE(*,*) DUMP(1:LEND)
         IRC=291
         RETURN
      ENDIF
C
      RETURN
      END


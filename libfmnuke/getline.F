      subroutine GETLINE(MAXUNIT,UNITA,FILES,LUNIT,CUNIT,
     &     CNT,MAXHDR,CODHDR,LINE,
     &     BUFF1000,LENB,ENDOFF,IRC)
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
      INTEGER MAXUNIT
      INTEGER UNITA(MAXUNIT),CUNIT,LUNIT(MAXUNIT)
      CHARACTER*250 FILES(MAXUNIT)
      INTEGER CNT,MAXHDR,CODHDR(MAXHDR),LINE
      CHARACTER*1000 BUFF1000
      integer lenb
      logical endoff
      integer irc
C
      INTEGER FTUNIT, LENGTH, LEND,LENF,IUNIT
      EXTERNAL FTUNIT, LENGTH
      LOGICAL BDONE
      CHARACTER*250 BUFF250
      CHARACTER*8 MYNAME
      DATA MYNAME /'GETLINE'/
C
      DATA IUNIT/-1/
C
      IF (IUNIT.NE.-1 .AND. IUNIT.LT.CUNIT) THEN
C     
C     OPEN THE NEW FILE
C     
         IF (CUNIT.GT.MAXUNIT) THEN
            WRITE(*,*) MYNAME,'Too deep file nesting.'
            IRC=275
            RETURN
         END IF
C     
         LUNIT(CUNIT)=0
         UNITA(CUNIT)=FTUNIT(IRC)
         IF (IRC.NE.0) RETURN
C     
         LENF=LENGTH(FILES(CUNIT),250,10)
         OPEN(UNIT=UNITA(CUNIT),FILE=FILES(CUNIT)(1:LENF),
     &        ACCESS='SEQUENTIAL',FORM='FORMATTED',
     &        STATUS='OLD',IOSTAT=IRC)
C     
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Unable to open file:'//
     &           FILES(CUNIT)(1:LENF)
            RETURN
         END IF
      END IF
C
      IUNIT=CUNIT
C
C     READ NEXT LINE
C
      BUFF1000=''
      LENB=0
      BDONE=.FALSE.
      DO WHILE (.NOT.BDONE)
         READ(UNITA(CUNIT),'(A250)',END=99,IOSTAT=IRC) BUFF250
         LUNIT(CUNIT)=LUNIT(CUNIT)+1
         IF (IRC.NE.0) THEN
            LENF=LENGTH(FILES(CUNIT),250,10)
            WRITE(*,*) 
     &           '"'//FILES(CUNIT)(1:LENF)//'", line ',
     &           LUNIT(CUNIT),': Error: Unable to read line.'
            RETURN
         END IF
C
         CALL CHOP (BUFF250,250)
         LEND=LENGTH(BUFF250,250,10)
C
C     IGNORE COMMENTS
C
         IF (LEND.NE.0 .AND. BUFF250(1:1).NE.'#') THEN
            BUFF1000=BUFF1000(1:LENB)//BUFF250(1:LEND)
            LENB=LENB+LEND
            IF (BUFF1000(LENB:LENB).EQ.'~') THEN
               BUFF1000(LENB:LENB)=' '
               LENB=LENB-1
               BDONE=.FALSE.
            ELSE
               BDONE=.TRUE.
            END IF
         END IF
C
         GOTO 250
 99      CONTINUE
         IRC=0                  ! NO ERROR, JUST EOF
         IF (CUNIT.GT.1) THEN
            CLOSE(UNITA(CUNIT),IOSTAT=IRC)
C
            CUNIT=CUNIT-1
            IUNIT=CUNIT
         ELSE
C
C     END OF FILE
C
            IF (LINE.NE.-1) THEN
               IF (CNT.LT.CODHDR(LINE)
     &              .AND. CODHDR(LINE).NE.-1) THEN
                  WRITE(*,*) MYNAME,'EOF interrupts data body.'
                  IRC=174
                  RETURN
               END IF
            END IF
            ENDOFF=.TRUE.
            BDONE=.TRUE.
         END IF
C
 250     CONTINUE
      END DO
      RETURN
C
      END

      

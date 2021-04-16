      SUBROUTINE NUKE(KODE,UNITI,HDR250,INTOUT,DAT250,NRDAT,IRC)
C     ***********************************************************************
C     +                                                                     *
C     INTERFACE SUBROUTINE TO FORTRAN-PROGRAM INPUT FILES                   *
C     +                                                                     *
CI    KODE (I*4) = INITIALISATION CODE, KODE=-1 REINITIALISES               *
CI    UNITI (I*4) = UNIT NUMBER FOR INPUT FILE (5=STANDARD INPUT)           *
CI    HDR250(250) (C*250) = HEADERS TO SEARCH FOR IN INPUT FILE             *
C     +   HEADER MUST BE ON FORM 'THIS IS THE HEADER [N]P', WHERE N         *
C     +   IS THE NUMBER OF LINES INCLUDED IN DATA LINE (CONCATENATED).      *
C     +   IF 'N'='*' THEN ONE DATA LINE IS RETURNED PER LINE IN THE FILE    *
C     +   (NOTE: THE HEADER IS IN THIS CASE ONLY NEEDED ONCE).              *
C     +   A 'P' AT THE END IMPLIES THAT THE DATA SHOULD BE PREPROCESSED.    *
CO    INTOUT(NRDAT) (I*4) = HEADER ID (LOC IN ARRAY) FOR DATA FOUND IN FILE *
CO    DAT250 (C*250) = DATA BODY                                     *
C     +      II=INTOUT(JJ) => HEADER HDR250(II) HAD DATABODY DAT250(JJ)     *
CI    NRDAT (I*4) = NUMBER OF DATA LINES (THE PROGRAM WILL PICK UP FROM     *
C     +    ITS LAST POSITION IN FILE IF IT IS CALLED SEVERAL TIMES).        *
C     +    IF THE ROUTINE REACHES END OF FILE (OR 'stop') BEFORE IT HAS     *
C     +    FILLED NRDAT LINES, THE LINE AFTER LAST LINE WILL HAVE INOUT=0   *
CIO   IRC (I*4) = ERROR RETURN CODE (IRC=0, ALL OK)                         *
C     +     A CALL WITH (IRC.NE.0) WILL INITIALISE THE ROUTINE              *
C     +                                                                     *
C     NOTE :                                                                *
C     +   o THE LINE '?' IN INPUT FILE WILL LIST ALL HEADERS AND COMMANDS   *
C     +   o IT IS POSSIBLE TO BYPASS THE MATHEMATICAL PREPROCESSOR          *
C     +     (IF IT FAILS TO WORK). IF THE MATHEMATICAL EXPRESSION IS        *
C     +     CONTAINED IN {}, IT WILL NOT BE EVALUATED ({} ARE DELETED), AN  *
C     +     ALTERNATIVE IS [] ([] ARE CONVERTED INTO ()).                   *
C     +                                                                     *
C     EXAMPLE OF CALL SEQUENCE :                                            *
C     +                                                                     *
C     INTEGER NRDAT,IRC,UNITI                                             *
C     PARAMETER (NRDAT=1)                                                   *
C     CHARACTER*250 HDR250(250),DAT250                               *
C     INTEGER INTOUT(NRDAT)                                               *
C     UNITI=5                                                               *
C     HDR250(1)='THIS IS A HEADER : [*]P'                                   *
C     CALL NUKE(KODE,UNITI,HDR250,INTOUT,DAT250,NRDAT,IRC)                       *
C     +                                                                     *
C     EXAMPLE OF INPUT FILE (THE 'C     ' AND '*' SHOULD BE IGNORED)        *
C     +                                                                     *
C     # first line (comment)                                                *
C     # test input file (yet another comment)                               *
C       THIS IS A HEADER :                                                  *
C        data-string data1                                                  *
C        1+2 +1.0D-4^-2 +23/1.01D2 +(1*(-1*(1*(-1*(1*(-2)))))) data2        *
C     # (...and another comment) The next line gives (some) help            *
C        ?                                                                  *
C     # last line (last comment)                                            *
C     +                                                                     *
C     THE DATA IN THE ABOVE FILE IS INTERPRETED AS :                        *
C     data-string data1                                                     *
C     3 +100000000. +0.22772277227723 +(-2) data2                           *
C     +                                                                     *
C     VERSION                      : 18/07/95                               *
C     +                                                                     *
C     WRITTEN/MODIFIED BY:                                                  *
C     --------------------------------------------------------------------- *
C     |    NAME      |   DATE   |                 REASON                  | *
C     --------------------------------------------------------------------- *
C     | F. TVETER    | 18/07/95 | NEW                                     | *
C     | S. EASTWOOD  | 29/02/08 | BUG IN LINE 1120 (in g77)               | *
C     |              |          |                                         | *
C     --------------------------------------------------------------------- *
C     ***********************************************************************
C
      IMPLICIT NONE
      SAVE
C
      INTEGER NRHDR
      PARAMETER (NRHDR=250)
C
C     INTERFACE VARIABLES
C
      INTEGER KODE,UNITI,NRDAT
      CHARACTER*250 HDR250(NRHDR),DAT250(NRDAT)
      INTEGER INTOUT(NRDAT),IRC
C
C     INTERNAL VARIABLES
C
      INTEGER NROPT
      PARAMETER (NROPT=5)
      LOGICAL PPDAT(NROPT,NRHDR)
      INTEGER NRLEN(NRHDR),CODHDR(NRHDR),MAXHDR,CURHDR
      LOGICAL ACTHDR(NRHDR)
C
      INTEGER MAXCOM,COMMAND,IRC1,
     &     LENGTH,LENS,LENB,LEND,LENF,LENH,LEN2,LEN3,MEAN
      PARAMETER (MAXCOM=12)
      CHARACTER*20 COM20(MAXCOM),HLP20(MAXCOM)
      CHARACTER*250 STRING,DUMP,JUNK,NUKEHEAD
      CHARACTER*1000 BUFF,BUFF2,BUFF3,PATHF
      CHARACTER*1 DEL(2),C1
      EXTERNAL FTUNIT,LENGTH
      INTEGER CMLEN(MAXCOM),system,putenv,FTUNIT,LEV
      EXTERNAL system,putenv
C
      INTEGER    II,JJ,KK,LINE,CNT,
     &     CPOS,POS(2),IFLVL
C
      CHARACTER*3 PRE
      CHARACTER*1 OPT,OPTL(NROPT)
      DATA OPTL /'V','F','M','L','R'/
C
      LOGICAL DONE,BDONE,SUCC,PROCESSED,ERROR,RDATA,
     &     SPLIT,FIRST,FOUND,ROCK,
     &     IFACT(10),IFDON(10),IFVAL,COMP,GO,DONCOM
      DATA FIRST /.TRUE./
C
C     AUXILIARY VARIABLES
C
      INTEGER MAXVAR
      PARAMETER (MAXVAR=50)
      INTEGER VARLEN(MAXVAR),NRVAR,CVAR
      CHARACTER*250 NAMVAR(MAXVAR),VALVAR(MAXVAR)
C
C     UNIT CONTROLLERS
C
      INTEGER MAXUNIT
      PARAMETER (MAXUNIT=10)
      INTEGER UNITA(MAXUNIT),CUNIT,LUNIT(MAXUNIT)
      CHARACTER*250 FILES(MAXUNIT)
C
C     EQUAL SIGN IDENTIFIERS
C
      INTEGER MAXEQN
      PARAMETER (MAXEQN=250)
      INTEGER POSEQN(MAXEQN),CEQN
C
      CHARACTER*8 MYNAME
      DATA MYNAME /'NUKE'/
      LOGICAL BDEB,ACTIVE
      DATA ACTIVE /.FALSE./
C
C     Debug System.
C
      IF (.NOT.ACTIVE) CALL DEBUG(MYNAME,BDEB,ACTIVE)
C
      IF (BDEB) WRITE(*,*) MYNAME,'Debug: Routine starts.',IRC
C
C     INITIALIZE (ALWAYS)
C
      ERROR=.FALSE.
      SUCC=.TRUE.
      CURHDR=0
      INTOUT(1)=0
C
C     CHECK IF WE SHOULD RE-INITIALISE
C
      IF (KODE.NE.0) THEN
         FIRST=.TRUE.
         KODE=0
      ENDIF
C
      IF (FIRST) THEN
C
         C1=CHAR(0)
C
C     DEFINE COMMANDS
C
         COM20(1)='if'
         HLP20(1)='(<expr>)'
         COM20(2)='elseif'
         HLP20(2)='(<expr>)'
         COM20(3)='else'
         CALL RESET(HLP20(3),20)
         COM20(4)='endif'
         CALL RESET(HLP20(4),20)
         COM20(5)='!'
         HLP20(5)='<UNIX command>'
         COM20(6)='include'
         HLP20(6)='(<file>)'
         COM20(7)='stop'
         CALL RESET(HLP20(7),20)
         COM20(8)='exit'
         HLP20(8)='<exit code>'
         COM20(9)='echo'
         HLP20(9)='<text>'
         COM20(10)='set'
         HLP20(10)='<var>=<expr[$(var)]>'
         COM20(11)='import'
         HLP20(11)='<var1> <var2> ...'
         COM20(12)='export'
         HLP20(12)='<var1> <var2> ...'
C
C     ANALYZE HEADERS
C
         LENH=1
         MEAN=0
         II=1
         DO WHILE(II.LE.NRHDR)
            CALL CHOP(HDR250(II),250)
            LENH=LENGTH(HDR250(II),250,LENH)
            IF (LENH.EQ.0) THEN
C     THE HEADER IS BLANK
               ACTHDR(II)=.FALSE.
            ELSE
               ACTHDR(II)=.TRUE.
               MAXHDR=II
C
C     FIND LENGTH OF HEADER, AND LOCATION OF []...
C
               DEL(1)='['
               DEL(2)=']'
               LEV=0
               POS(1)=1
               POS(2)=250
               IRC1=0
               CALL ITEM(HDR250(II),250,DEL,LEV,POS,IRC1)
               IF (IRC1.NE.0 .OR. LEV .EQ. -1) THEN
                  WRITE(*,*) MYNAME,
     &                 'Error: missing [] in header, ',II,'.'
                  IF(IRC.EQ.0)IRC=250
                  RETURN
               ENDIF
               NRLEN(II)=POS(1)-1
               MEAN=MEAN+NRLEN(II)
               IF (NRLEN(II).GT.0) THEN
                  IF (HDR250(II)(NRLEN(II):NRLEN(II)).EQ.' ')
     &                 NRLEN(II)=NRLEN(II)-1
               ELSE
                  WRITE(*,*) MYNAME,'Error: header too short, ',II,'.'
                  IF(IRC.EQ.0)IRC=101
               ENDIF
C
C     CHECK FOR PRE-PROCESSOR COMMANDS
C
               DO KK=1,NROPT
                  PPDAT(KK,II)=.FALSE.
               ENDDO
               DO JJ=POS(2)+1,LENH
                  OPT=HDR250(II)(JJ:JJ)
                  FOUND=.FALSE.
                  DO KK=1,NROPT
                     IF (OPT.EQ.OPTL(KK)) THEN
                        PPDAT(KK,II)=.TRUE.
                        FOUND=.TRUE.
                     ENDIF
                  ENDDO
                  IF (.NOT.FOUND) THEN
                     WRITE(*,*) MYNAME,
     &                    'Error: unknown header option, ',
     &                    OPT,'.'
                  ENDIF
               ENDDO
C
C     GET INFO CONCERNING THE NUMBER OF LINES IN THE DATA BODY
C
               READ(HDR250(II)(POS(1)+1:POS(2)-1),*,ERR=99) CODHDR(II)
 250           GOTO 101
 99            IF (HDR250(II)(POS(1)+1:POS(1)+1).EQ.'*'
     &              .OR.HDR250(II)(POS(1)+2:POS(1)+2).EQ.'*') THEN
                  CODHDR(II)=-1
               ELSE
                  CODHDR(II)=0
               ENDIF
 101           CONTINUE
C
            ENDIF
            II=II+1
         ENDDO
C
C     MEAN HEADER LENGTH
C
         IF (NRHDR.GT.0) MEAN=INT(MEAN/NRHDR)
C
C     ANALYZE COMMANDS
C
         DO II=1,MAXCOM
            CALL CHOP(COM20(II),20)
            CALL CHOP(HLP20(II),20)
C     FIND LENGTH OF COMMAND
            CMLEN(II)=LENGTH(COM20(II),20,3)
         ENDDO
C
C     INITIALIZE
C
         CNT=0
         SPLIT=.FALSE.
C
C
C     READ DATA FROM INPUT FILE..............................
C
C     INITIALIZE
C
         IFLVL=0
         CUNIT=1
         UNITA(CUNIT)=UNITI
         LUNIT(CUNIT)=0
         FILES(CUNIT)='input file'
         RDATA=.FALSE.
C
C     READ UNTIL WE FIND A '#! NUKE' STATEMENT IN INPUT FILE (=line 1)
C
         LUNIT(CUNIT)=1
         ROCK=.TRUE.
         DO WHILE(ROCK)
            READ(UNITA(CUNIT),'(A250)',ERR=50,END=50) STRING
            BUFF(1:250)=STRING
            CALL CHOP(BUFF,250)
            IF (BUFF(1:7).EQ.'#! NUKE' .OR.
     &           BUFF(1:6).EQ.'#!NUKE') THEN
               ROCK=.FALSE.
            ENDIF
            GOTO 55
 50         CUNIT=CUNIT-1
            WRITE(*,*) MYNAME,'Error: input file corrupted.'
            IRC=102
            ROCK=.FALSE.
 55         CONTINUE
         ENDDO
C
         FIRST=(IRC.NE.0)
      ENDIF
C
      PROCESSED=.FALSE.
      DONE=(CUNIT.LE.0.OR.IRC.NE.0)
C
      DO WHILE(.NOT.DONE)
C
C     CHECK IF INPUT IS SPLIT INTO SEVERAL LINES
C
         IF (.NOT. PROCESSED) THEN
C
            IF (SPLIT) THEN
C     LAST LINE WAS SPLIT
            ELSE
               CALL RESET(BUFF,1000)
               CPOS=1
            ENDIF
C
C     ADD TO BUFFER
C
            IF (CPOS.LT.900) THEN
               BUFF(CPOS:CPOS+250)=
     &              STRING(1:250)
               CALL CHOP(BUFF,CPOS+250)
               LENB=LENGTH(BUFF,CPOS+250,CPOS)
            ELSE
               PRE=': E'
               PROCESSED=.TRUE.
               ERROR=.TRUE.
               BUFF='buffer too small (recompile)'
            ENDIF
C
            SPLIT=.FALSE.
C
C     REDEFINE END OF LINE
C
            IF (BUFF(LENB:LENB).EQ.'&') THEN
               SPLIT=.TRUE.
               CPOS=LENB
               BUFF(LENB:LENB)=' '
               PRE=': &'
               PROCESSED=.TRUE.
            ENDIF
         ENDIF
C
C--------CHECK IF THE LINE IS A COMMENT
C
         IF (.NOT. PROCESSED) THEN
            IF (BUFF(1:1).EQ.'#' .OR.
     &           LENB.EQ.0) THEN
               PRE=': #'
               PROCESSED=.TRUE.
            ENDIF
         ENDIF
C
C--------CHECK IF THE LINE STARTS WITH A COMMAND
C
         IF (.NOT.PROCESSED) THEN
            COMMAND=0
            JJ=1
            FOUND=.FALSE.
            DO WHILE(.NOT.FOUND .AND. JJ.LE.MAXCOM)
               FOUND=(CMLEN(JJ).LE.LENB)
               KK=1
               DO WHILE(KK.LE.CMLEN(JJ).AND.FOUND)
                  IF(FOUND) FOUND=(BUFF(KK:KK).EQ.COM20(JJ)(KK:KK))
                  KK=KK+1
               ENDDO
               IF (FOUND)THEN
                  COMMAND=JJ
               ELSE
                  JJ=JJ+1
               ENDIF
            ENDDO
C
C--------PERFORM COMMAND (IF LINE STARTS WITH IT)
C
            IF (COMMAND.NE.0) THEN
C
               PRE=': C'
C
               DONCOM=.FALSE.
               IF (COMMAND.EQ.1) THEN
                  DONCOM=.TRUE.
                  DEL(1)='('
                  DEL(2)=')'
                  LEV=0
                  POS(1)=1
                  POS(2)=1000
                  CALL ITEM(BUFF,1000,DEL,LEV,POS,IRC)
                  CALL RESET(BUFF2,1000)
                  BUFF2=BUFF(1:POS(1)-1)
                  CALL CHOP(BUFF2,1000)
                  LEN2=LENGTH(BUFF2,1000,LEN2)
                  IF (IRC.NE.0) THEN
C     SYNTAX ERROR
                     PRE=': E'
                     PROCESSED=.TRUE.
                     ERROR=.TRUE.
                     BUFF='missing () in if statement.'
                  ELSEIF (COM20(1)(1:CMLEN(1)).EQ.BUFF2(1:LEN2).AND.
     &                    LENB.EQ.POS(2)) THEN
C     ANALYZE TEST VALUE
                     IF ((POS(2)-POS(1)).LE.1) THEN
                        IFVAL=.TRUE.
                     ELSE
                        CALL RESET(BUFF3,1000)
                        BUFF3=BUFF(POS(1)+1:POS(2)-1)
                        CALL EVAL(BUFF3,1000,NRVAR,NAMVAR,
     &                       VALVAR,VARLEN,PPDAT(1,1),IRC)
                        IF (IRC.NE.0) RETURN
                        CALL CHOP(BUFF3,1000)
                        LEN3=LENGTH(BUFF3,1000,1)
                        IFVAL=(BUFF3(1:LEN3).EQ.'T')
                     ENDIF
C
                     IFLVL=IFLVL+1
C
                     IF (IFVAL) THEN
                        IFACT(IFLVL)=.TRUE.
                        IFDON(IFLVL)=.TRUE.
                     ELSE
                        IFACT(IFLVL)=.FALSE.
                        IFDON(IFLVL)=.FALSE.
                     ENDIF
C
                  ELSE
C     SYNTAX ERROR
                     PRE=': E'
                     PROCESSED=.TRUE.
                     ERROR=.TRUE.
                     BUFF='strange if statement.'
                  ENDIF
               ELSEIF (COMMAND.EQ.2 .AND. IFLVL.GT.0) THEN
                  DONCOM=.TRUE.
                  IF (.NOT.IFDON(IFLVL)) THEN
                     DEL(1)='('
                     DEL(2)=')'
                     LEV=0
                     POS(1)=1
                     POS(2)=1000
                     CALL ITEM(BUFF,1000,DEL,LEV,POS,IRC)
                     CALL RESET(BUFF2,1000)
                     BUFF2=BUFF(1:POS(1)-1)
                     CALL CHOP(BUFF2,1000)
                     LEN2=LENGTH(BUFF2,1000,LEN2)
                     IF (IRC.NE.0) THEN
C     SYNTAX ERROR
                        PRE=': E'
                        PROCESSED=.TRUE.
                        ERROR=.TRUE.
                        BUFF='missing () in elseif statement.'
                     ELSEIF (COM20(1)(1:CMLEN(1))
     &                       .EQ.BUFF2(1:LEN2).AND.
     &                    LENB.EQ.POS(2)) THEN
C     ANALYZE TEST VALUE
                        IF ((POS(2)-POS(1)).LE.1) THEN
                           IFVAL=.TRUE.
                        ELSE
                           CALL RESET(BUFF3,1000)
                           BUFF3=BUFF(POS(1)+1:POS(2)-1)
                           CALL EVAL(BUFF3,1000,NRVAR,NAMVAR,
     &                          VALVAR,VARLEN,PPDAT(1,1),IRC)
                           IF (IRC.NE.0) RETURN
                           CALL CHOP(BUFF3,1000)
                           LEN3=LENGTH(BUFF3,1000,1)
                           IFVAL=(BUFF3(1:LEN3).EQ.'T')
                        ENDIF
C
                        IF (IFVAL) THEN
                           IFACT(IFLVL)=.TRUE.
                           IFDON(IFLVL)=.TRUE.
                        ELSE
                           IFACT(IFLVL)=.FALSE.
                           IFDON(IFLVL)=.FALSE.
                        ENDIF
C
                     ELSE
                        PRE=': E'
                        PROCESSED=.TRUE.
                        ERROR=.TRUE.
                        BUFF='strange elseif statement.'
                     ENDIF
                  ELSE
                     IFACT(IFLVL)=.FALSE.
                  ENDIF
               ELSEIF (COMMAND.EQ.3 .AND. IFLVL.GT.0) THEN
                  DONCOM=.TRUE.
                  IF (.NOT.IFDON(IFLVL)) THEN
                     IFACT(IFLVL)=.TRUE.
                     IFDON(IFLVL)=.TRUE.
                  ELSE
                     IFACT(IFLVL)=.FALSE.
                  ENDIF
               ELSEIF (COMMAND.EQ.4 .AND. IFLVL.GT.0) THEN
                  DONCOM=.TRUE.
                  IFACT(IFLVL)=.FALSE.
                  IFDON(IFLVL)=.FALSE.
                  IFLVL=IFLVL-1
               ENDIF
C
C     MAKE SURE WE REIGSTER NESTING ERRORS
C
               IF (.NOT.DONCOM .AND. COMMAND.LE.4) THEN
                  PRE=': E'
                  PROCESSED=.TRUE.
                  ERROR=.TRUE.
                  BUFF='elseif, else or endif out of place.'
               ENDIF
C
C     COMMANDS ONLY CHECKED IF WE ARE IN AN ACTIVE REGION
C
               GO = (IFLVL.EQ.0)
               IF (.NOT.GO) GO=(IFACT(IFLVL))
C
C     EVALUATE EXPRESSIONS (IN CASE OF ECHO OR INCLUDE)...
C
               IF (GO .AND. .NOT.DONCOM) THEN
                  IF (BDEB) WRITE(*,*) MYNAME,
     &                 'Debug: Calling EVAL.',IRC
                  CALL EVAL(BUFF,1000,NRVAR,NAMVAR,
     &                 VALVAR,VARLEN,PPDAT(1,1),IRC)
                  IF (IRC.NE.0) THEN
                     WRITE(*,*) MYNAME,'Error return from EVAL',IRC
                     RETURN
                  ENDIF
                  LENB=LENGTH(BUFF,1000,LENB)
               ENDIF
C
C               IF ((COMMAND.EQ.5.OR.COMMAND.EQ.6).AND. GO) THEN
C                  DONCOM=.TRUE.
C                  PRE=': E'
C                  PROCESSED=.TRUE.
C                  ERROR=.TRUE.
C                  LENB=LENGTH(BUFF,1000,LENB)
C                  BUFF='Permission denied: '//BUFF(1:LENB)
C               ELSEIF (COMMAND.EQ.5 .AND. GO) THEN
               IF (COMMAND.EQ.5 .AND. GO) THEN
C     EXECUTE UNIX COMMAND
                  DONCOM=.TRUE.
                  IRC1=system(BUFF(2:LENB))
                  IF (IRC1.NE.0) THEN
                     SUCC=.FALSE.
                     IRC1=0
                  ENDIF
               ELSEIF (COMMAND.EQ.6.AND.GO) THEN
C     INCLUDE FILE
                  DONCOM=.TRUE.
                  DEL(1)='('
                  DEL(2)=')'
                  LEV=0
                  POS(1)=1
                  POS(2)=LENB
                  CALL ITEM(BUFF,LENB,DEL,LEV,POS,IRC)
                  IF (IRC.NE.0) THEN
C     SYNTAX ERROR
                     PRE=': E'
                     PROCESSED=.TRUE.
                     ERROR=.TRUE.
                     BUFF='missing () in include statement.'
                  ELSE
                     PATHF=BUFF(POS(1)+1:POS(2)-1)
                     CALL CHOP(PATHF,LENB)
                     IF (CUNIT.LT.MAXUNIT) THEN
                        CUNIT=CUNIT+1
                        LUNIT(CUNIT)=0
                        UNITA(CUNIT)=FTUNIT(IRC)
                        IF (IRC.NE.0) THEN
                           WRITE(*,*) MYNAME,
     &                          'Error return from DEBUG.',IRC
                           RETURN
                        ENDIF
                        OPEN(UNIT=UNITA(CUNIT),FILE=PATHF,
     &                       ACCESS='SEQUENTIAL',FORM='FORMATTED',
     &                       STATUS='OLD',ERR=567)
                        FILES(CUNIT)=PATHF(1:250)
                        GOTO 568
 567                    CONTINUE
                        CUNIT=CUNIT-1
                        LENF=LENGTH(PATHF,250,10)
                        PRE=': E'
                        PROCESSED=.TRUE.
                        ERROR=.TRUE.
                        BUFF='corrupted file: '//PATHF(1:LENF)
 568                    CONTINUE
                     ELSE
                        LENF=LENGTH(PATHF,250,10)
                        PRE=': E'
                        PROCESSED=.TRUE.
                        ERROR=.TRUE.
                        BUFF='too deep file structure: '//
     &                       PATHF(1:LENF)
                     ENDIF
                  ENDIF
               ELSEIF (COMMAND.EQ.7.AND.GO) THEN
C     STOP COMMAND
                  GOTO 200
               ELSEIF (COMMAND.EQ.8.AND.GO) THEN
C     EXIT COMMAND
                  DONCOM=.TRUE.
                  JUNK=NUKEHEAD(BUFF,250)
                  READ(BUFF,*,ERR=10,END=10) IRC
 10               CALL exit(IRC)
               ELSEIF (COMMAND.EQ.9.AND.GO) THEN
C     WRITE COMMAND
                  DONCOM=.TRUE.
                  JUNK=NUKEHEAD(BUFF,250)
                  CALL CHOP(BUFF,LENB)
                  LENB=LENGTH(BUFF,LENB,LENB)
                  WRITE(*,*) BUFF(1:LENB)
               ELSEIF (COMMAND.EQ.10 .AND. GO) THEN
C     SET COMMAND
                  DONCOM=.TRUE.
C     IDENTIFY ALL EQUAL SIGNS
                  COMP=.FALSE.
                  CEQN=0
                  LENB=LENGTH(BUFF,1000,LENB)
                  DO II=1,LENB
                     IF (BUFF(II:II).EQ.'=') THEN
                        CEQN=CEQN+1
                        POSEQN(CEQN)=II
                     ENDIF
                  ENDDO
C
C     LOOP THROUGH ALL EQUATIONS
C
                  IF (CEQN.GT.0) THEN
C     FOUND AT LEAST ONE '='
                     DO KK=CEQN,1,-1
C
                        IF (KK.GT.1) THEN
                           POS(1)=POSEQN(KK-1)+1
                        ELSE
                           POS(1)=4
                        ENDIF
C
                        IF (KK.LT.CEQN) THEN
                           POS(2)=POSEQN(KK+1)-1
                        ELSE
                           POS(2)=1000
                        ENDIF
C
C     GET PART OF STRING BEFORE AND AFTER '='
C
                        CALL RESET(BUFF2,1000)
                        BUFF2=BUFF(POS(1):POSEQN(KK)-1)
                        IF (BDEB) WRITE(*,*) MYNAME,
     &                       'Debug: Calling EVAL.',IRC
                        CALL EVAL(BUFF2,1000,NRVAR,NAMVAR,
     &                       VALVAR,VARLEN,PPDAT(1,1),IRC)
                        IF (IRC.NE.0) RETURN
                        CALL CHOP(BUFF2,1000)
                        LEN2=LENGTH(BUFF2,250,LEN2)
                        CALL RESET(BUFF3,1000)
                        BUFF3=BUFF(POSEQN(KK)+1:POS(2))
                        IF (BDEB) WRITE(*,*) MYNAME,
     &                       'Debug: Calling EVAL.',IRC
                        CALL EVAL(BUFF3,1000,NRVAR,NAMVAR,
     &                       VALVAR,VARLEN,PPDAT(1,1),IRC)
                        IF (IRC.NE.0) RETURN
                        CALL CHOP(BUFF3,1000)
                        LEN3=LENGTH(BUFF3,250,LEN3)
C
C     FIND OUT IF VARIABLE EXISTS
C
                        II=0
                        FOUND=.FALSE.
                        DO WHILE (II.LT.NRVAR .AND. .NOT.FOUND)
                           II=II+1
                           JJ=0
                           FOUND=(VARLEN(II).EQ.LEN2.AND.LEN2.NE.0)
                           DO WHILE (JJ.LT.LEN2 .AND. FOUND)
                              JJ=JJ+1
                              IF (FOUND) FOUND=(BUFF2(JJ:JJ).EQ.
     &                             NAMVAR(II)(JJ:JJ))
                           ENDDO
                        ENDDO
                        CVAR=II
C
                        LEN3=MIN(LEN3,1)
                        IF (BUFF3(1:LEN3).EQ.'.'.OR.
     &                       BUFF3(1:LEN3).EQ.C1) THEN
                           IF (FOUND) THEN
C     REMOVE VARIABLE
                              NAMVAR(CVAR)=NAMVAR(NRVAR)
                              VALVAR(CVAR)=VALVAR(NRVAR)
                              VARLEN(CVAR)=VARLEN(NRVAR)
                              NRVAR=NRVAR-1
                           ELSE
                              PRE=': E'
                              PROCESSED=.TRUE.
                              ERROR=.TRUE.
                              BUFF='attempt to remove non-existent'//
     &                          ' variable.'
                           ENDIF
                        ELSE
C     ADD/UPDATE VARIABLE
                           IF (FOUND) THEN
                              VALVAR(CVAR)=BUFF3(1:250)
                           ELSE
                              IF (NRVAR.LT.MAXVAR) THEN
                                 NRVAR=NRVAR+1
                                 NAMVAR(NRVAR)=BUFF2(1:250)
                                 VALVAR(NRVAR)=BUFF3(1:250)
                                 VARLEN(NRVAR)=
     &                                LENGTH(NAMVAR(NRVAR),250,5)
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ELSEIF (COMMAND.EQ.11 .AND. GO) THEN
C     IMPORT
                  DONCOM=.TRUE.
C     FIRST REMOVE THE ACTUAL COMMAND...
                  BUFF2=NUKEHEAD(BUFF,250)
C     READ THE FIRST ARGUMENT...
                  BUFF2=NUKEHEAD(BUFF,250)
                  CALL CHOP(BUFF2,1000)
                  LEN2=LENGTH(BUFF2,250,LEN2)
                  BDONE = (LEN2.LE.0)
                  DO WHILE (.NOT. BDONE)
                     IF (BDEB) WRITE(*,*) MYNAME,
     &                    'Debug: Calling EVAL.',IRC
                     CALL EVAL(BUFF2,1000,NRVAR,NAMVAR,
     &                    VALVAR,VARLEN,PPDAT(1,1),IRC)
                     IF (IRC.NE.0) RETURN
                     CALL CHOP(BUFF2,1000)
                     LEN2=LENGTH(BUFF2,250,LEN2)
                     CALL RESET(BUFF3,1000)
C     IMPORT ENVIRONMENT VARIABLE FROM SHELL...
                     call getenv(buff2(1:len2),buff3)
                     IF (BDEB) WRITE(*,*) MYNAME,
     &                    'Debug: Calling EVAL.',IRC
C                     CALL EVAL(BUFF3,1000,NRVAR,NAMVAR,
C     &                    VALVAR,VARLEN,PPDAT(1,1),IRC)
                     IF (IRC.NE.0) RETURN
                     CALL CHOP(BUFF3,1000)
                     LEN3=LENGTH(BUFF3,250,LEN3)
C
C     FIND OUT IF VARIABLE EXISTS
C
                     II=0
                     FOUND=.FALSE.
                     DO WHILE (II.LT.NRVAR .AND. .NOT.FOUND)
                        II=II+1
                        JJ=0
                        FOUND=(VARLEN(II).EQ.LEN2.AND.LEN2.NE.0)
                        DO WHILE (JJ.LT.LEN2 .AND. FOUND)
                           JJ=JJ+1
                           IF (FOUND) FOUND=(BUFF2(JJ:JJ).EQ.
     &                          NAMVAR(II)(JJ:JJ))
                        ENDDO
                     ENDDO
                     CVAR=II
C
                     LEN3=MIN(LEN3,1)
                     IF (BUFF3(1:LEN3).EQ.'.'.OR.
     &                    BUFF3(1:LEN3).EQ.C1) THEN
                        IF (FOUND) THEN
C     REMOVE VARIABLE
                           NAMVAR(CVAR)=NAMVAR(NRVAR)
                           VALVAR(CVAR)=VALVAR(NRVAR)
                           VARLEN(CVAR)=VARLEN(NRVAR)
                           NRVAR=NRVAR-1
                        ELSE
                           PRE=': E'
                           PROCESSED=.TRUE.
                           ERROR=.TRUE.
                           BUFF='attempt to remove non-existent'//
     &                          ' variable.'
                        ENDIF
                     ELSE
C     ADD/UPDATE VARIABLE
                        IF (FOUND) THEN
                           VALVAR(CVAR)=BUFF3(1:250)
                        ELSE
                           IF (NRVAR.LT.MAXVAR) THEN
                              NRVAR=NRVAR+1
                              NAMVAR(NRVAR)=BUFF2(1:250)
                              VALVAR(NRVAR)=BUFF3(1:250)
                              VARLEN(NRVAR)=
     &                             LENGTH(NAMVAR(NRVAR),250,5)
                           ENDIF
                        ENDIF
                     ENDIF
C
C     READ THE NEXT ARGUMENT...
C
                     IF (ERROR) THEN
                        BDONE=.TRUE.
                     ELSE
                        BUFF2=NUKEHEAD(BUFF,250)
                        CALL CHOP(BUFF2,1000)
                        LEN2=LENGTH(BUFF2,250,LEN2)
                        BDONE = (LEN2.LE.0)
                     END IF
                  END DO
               ELSEIF (COMMAND.EQ.12 .AND. GO) THEN
C     EXPORT
                  DONCOM=.TRUE.
C     FIRST REMOVE THE ACTUAL COMMAND...
                  BUFF2=NUKEHEAD(BUFF,250)
C     READ THE FIRST ARGUMENT...
                  BUFF2=NUKEHEAD(BUFF,250)
                  CALL CHOP(BUFF2,1000)
                  LEN2=LENGTH(BUFF2,250,LEN2)
                  BDONE = (LEN2.LE.0)
                  DO WHILE (.NOT. BDONE)
                     IF (BDEB) WRITE(*,*) MYNAME,
     &                    'Debug: Calling EVAL.',IRC
                     CALL EVAL(BUFF2,1000,NRVAR,NAMVAR,
     &                    VALVAR,VARLEN,PPDAT(1,1),IRC)
                     IF (IRC.NE.0) RETURN
                     CALL CHOP(BUFF2,1000)
                     LEN2=LENGTH(BUFF2,250,LEN2)
C
C     FIND OUT IF VARIABLE EXISTS
C
                     II=0
                     FOUND=.FALSE.
                     DO WHILE (II.LT.NRVAR .AND. .NOT.FOUND)
                        II=II+1
                        JJ=0
                        FOUND=(VARLEN(II).EQ.LEN2.AND.LEN2.NE.0)
                        DO WHILE (JJ.LT.LEN2 .AND. FOUND)
                           JJ=JJ+1
                           IF (FOUND) FOUND=(BUFF2(JJ:JJ).EQ.
     &                          NAMVAR(II)(JJ:JJ))
                        ENDDO
                     ENDDO
                     CVAR=II
C
                     IF (FOUND) THEN
C     EXPORT VARIABLE TO SHELL...
                        BUFF3=VALVAR(CVAR)
                        IF (BDEB) WRITE(*,*) MYNAME,
     &                       'Debug: Calling EVAL.',IRC
                        CALL EVAL(BUFF3,1000,NRVAR,NAMVAR,
     &                       VALVAR,VARLEN,PPDAT(1,1),IRC)
                        IF (IRC.NE.0) RETURN
                        CALL CHOP(BUFF3,1000)
                        LEN3=LENGTH(BUFF3,250,LEN3)
                        BUFF3=BUFF2(1:LEN2)//'='//BUFF3(1:LEN3)
                        CALL CHOP(BUFF3,1000)
                        LEN3=LENGTH(BUFF3,250,LEN3)
                        IRC1=putenv(BUFF3(1:LEN3))
                        IF (IRC1.NE.0) THEN
                           SUCC=.FALSE.
                           IRC1=0
                        ENDIF
                     ELSE
                        PRE=': E'
                        PROCESSED=.TRUE.
                        ERROR=.TRUE.
                        BUFF='attempt to export non-existent'//
     &                       ' variable.'
                     ENDIF
C
C     READ THE NEXT ARGUMENT...
C
                     IF (ERROR) THEN
                        BDONE=.TRUE.
                     ELSE
                        BUFF2=NUKEHEAD(BUFF,250)
                        CALL CHOP(BUFF2,1000)
                        LEN2=LENGTH(BUFF2,250,LEN2)
                        BDONE = (LEN2.LE.0)
                     END IF
                  END DO
               ENDIF
               IF (DONCOM) PROCESSED=.TRUE.
            ENDIF
         ENDIF
C
C     MAKE SURE WE DISREGARD LINES IF WE ARE IN WRONG PART OF IF-TEST
C
         IF (IFLVL.GT.0.AND..NOT.PROCESSED.AND..NOT.GO) THEN
            PRE=': #'
            PROCESSED=.TRUE.
         ENDIF
C
C--------CHECK IF THE LINE IS A HELP STATEMENT
C
         IF (.NOT. PROCESSED) THEN
            IF (BUFF(1:1).EQ.'?') THEN
C     LIST ALL HEADERS
               WRITE(*,*) MYNAME,
     &              '--------LIST OF POSSIBLE HEADERS--------'
               DO II=1,NRHDR
                  LENS=LENGTH(HDR250(II),250,NRLEN(II)+3)
                  IF(ACTHDR(II))
     &                 WRITE(*,*) HDR250(II)(1:LENS)
               ENDDO
               WRITE(*,*) MYNAME,
     &              '--------LIST OF POSSIBLE COMMANDS-------'
               DO II=1,MAXCOM
                  LENS=LENGTH(COM20(II),20,3)
                  LENH=LENGTH(HLP20(II),20,3)
                  WRITE(*,*) COM20(II)(1:LENS)//' '//
     &                 HLP20(II)(1:LENH)
               ENDDO
               WRITE(*,*) MYNAME,
     &              '----------------------------------------'
               PRE=': ?'
               PROCESSED=.TRUE.
            ENDIF
         ENDIF
C
C--------CHECK IF THIS LINE IS A HEADER
C
         IF (.NOT.PROCESSED) THEN
            JJ=1
            FOUND=.FALSE.
            DO WHILE(.NOT.FOUND .AND. JJ.LE.MAXHDR)
               FOUND=((NRLEN(JJ).LE.LENB).AND.(ACTHDR(JJ)))
               KK=1
               DO WHILE(KK.LE.NRLEN(JJ).AND.FOUND)
                  IF(FOUND) FOUND=(BUFF(KK:KK).EQ.HDR250(JJ)(KK:KK))
                  KK=KK+1
               ENDDO
               IF (FOUND)THEN
C
C     CHECK IF WE ARE INTERUPTING ANOTHER DATA BODY
C
                  IF (RDATA) THEN
                     IF (CNT.NE.CODHDR(LINE)
     &                    .AND. CODHDR(LINE).NE.-1) THEN
                        PRE=': E'
                        PROCESSED=.TRUE.
                        ERROR=.TRUE.
                        BUFF='header in data body'
                     ENDIF
                  ENDIF
C
                  IF (.NOT.ERROR) THEN
C
C     PREPARE FOR DATABODY IF FIXED NUMBER OF DATA LINES
C
                     IF (CODHDR(JJ).NE.0) THEN
                        CNT=0
                        RDATA=.TRUE.
                     ELSE
                        RDATA=.FALSE.
                        CURHDR=CURHDR+1
                        INTOUT(CURHDR)=JJ
                        CALL RESET(DAT250(CURHDR),250)
                        IF (CURHDR.EQ.NRDAT) THEN
                           DONE=.TRUE.
                        ELSE
                           INTOUT(CURHDR+1)=0
                        ENDIF
                     ENDIF
                     LINE=JJ
                     PRE=': H'
                     PROCESSED=.TRUE.
                  ENDIF
               ELSE
                  JJ=JJ+1
               ENDIF
            ENDDO
         ENDIF
C
C--------CHECK IF THIS LINE SHOULD BE TREATED AS DATA BODY
C
         IF (.NOT.PROCESSED) THEN
            IF (RDATA) THEN
C
               PRE=': D'
               PROCESSED=.TRUE.
C
C     PRE-PROCESS LINE IF THIS IS REQUIRED
C
               IF (BDEB) WRITE(*,*) MYNAME,'Debug: Calling EVAL.',IRC
               CALL EVAL(BUFF,1000,NRVAR,NAMVAR,VALVAR,VARLEN,
     &              PPDAT(1,LINE),IRC)
               IF (IRC.NE.0) RETURN
               LENB=LENGTH(BUFF,1000,LENB)
C
               IF (LENB.EQ.0) THEN
C     DATA LINE IS EMPTY
                  WRITE(*,*) MYNAME,
     &                 'Warning: empty data line found under: '//
     &                 HDR250(LINE)(1:NRLEN(LINE))
               END IF
C
C     STORE DATA LINE
C
               IF (CODHDR(LINE).NE.0) THEN
C
                  CURHDR=CURHDR+1
                  INTOUT(CURHDR)=LINE
                  DAT250(CURHDR)=BUFF(1:250)
C
                  IF (CODHDR(LINE).GT.0) CNT=CNT+1
C
                  IF (CNT.EQ.CODHDR(LINE)) THEN
                     LINE=0
                     CNT=0
                     RDATA=.FALSE.
                  ENDIF
C
C     SET INDEX OF NEXT DATA BODY=0 (IN CASE THIS LINE IS THE LAST ONE)
C
                  IF (CURHDR.EQ.NRDAT) THEN
                     DONE=.TRUE.
                  ELSE
                     INTOUT(CURHDR+1)=0
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
C     DEBUG DUMP
C
         IF (BDEB) THEN
            CALL RESET(DUMP,250)
            IF (LENF.GT.80)LENF=80
            IF (LUNIT(CUNIT).EQ.0 .AND.CUNIT.GT.1) THEN
               CALL CHOP(FILES(CUNIT-1),250)
               LENF=LENGTH(FILES(CUNIT-1),250,10)
               WRITE(DUMP,*)
     &              '\"'//FILES(CUNIT-1)(1:LENF)//'\", line ',
     &              LUNIT(CUNIT-1),PRE,':'
            ELSE
               CALL CHOP(FILES(CUNIT),250)
               LENF=LENGTH(FILES(CUNIT),250,10)
               WRITE(DUMP,*) '\"'//FILES(CUNIT)(1:LENF)//'\", line ',
     &              LUNIT(CUNIT),PRE,':'
            ENDIF
            CALL CHOP(DUMP,250)
            LEND=LENGTH(DUMP,80,10)
            DUMP(LEND+1:LEND+1)=' '
            LENB=250-LEND-2
            WRITE(DUMP(LEND+2:250),*)STRING(1:LENB)
            CALL CHOP(DUMP,250)
            LEND=LENGTH(DUMP,250,LEND)
            WRITE(*,*) DUMP(1:LEND)
         ENDIF
C
C     IF WE'RE STILL NOT PROCESSED, SOMETHING IS WRONG
C
         IF (.NOT.PROCESSED .OR. ERROR) THEN
            CALL RESET(DUMP,250)
            CALL CHOP(FILES(CUNIT),250)
            LENF=LENGTH(FILES(CUNIT),250,10)
            IF (LENF.GT.80)LENF=80
            WRITE(DUMP,*) '\"'//FILES(CUNIT)(1:LENF)//'\", line ',
     &           LUNIT(CUNIT),': Error: '
            CALL CHOP(DUMP,250)
            LEND=LENGTH(DUMP,80,10)
            DUMP(LEND+1:LEND+1)=' '
            LENB=250-LEND-2
            WRITE(DUMP(LEND+2:250),*)BUFF(1:LENB)
            CALL CHOP(DUMP,250)
            LEND=LENGTH(DUMP,250,LEND)
            WRITE(*,*) DUMP(1:LEND)
            SUCC=.FALSE.
         ENDIF
C
         PROCESSED=.FALSE.
         ERROR=.FALSE.
C
 198     CONTINUE
         LUNIT(CUNIT)=LUNIT(CUNIT)+1
         READ(UNITA(CUNIT),'(A250)',ERR=200,END=200) STRING
         GOTO 205
 200     CONTINUE
         IF (UNITA(CUNIT).GT.30) CLOSE(UNITA(CUNIT))
C     CHECK IF END OF LINE INERRUPTED A DATA BODY
         IF (CUNIT.EQ.1.AND.CNT.NE.0) THEN
            BUFF='End of file interrupts data body.'
            CALL RESET(DUMP,250)
            CALL CHOP(FILES(CUNIT),250)
            LENF=LENGTH(FILES(CUNIT),250,10)
            IF (LENF.GT.80)LENF=80
            WRITE(DUMP,*) '\"'//FILES(CUNIT)(1:LENF)//'\", line ',
     &           LUNIT(CUNIT),': Error: '
            CALL CHOP(DUMP,250)
            LEND=LENGTH(DUMP,80,10)
            DUMP(LEND+1:LEND+1)=' '
            LENB=250-LEND-2
            WRITE(DUMP(LEND+2:250),*)BUFF(1:LENB)
            CALL CHOP(DUMP,250)
            LEND=LENGTH(DUMP,250,LEND)
            WRITE(*,*) DUMP(1:LEND)
            SUCC=.FALSE.
         ELSEIF (CUNIT.LE.0) THEN
            WRITE(*,*) MYNAME,'Internal error, file nesting.'
            IRC=103
            CALL exit(IRC)
         ENDIF
         CUNIT=CUNIT-1
C     READ NEXT LINE IN PREVIOUS FILE
         IF (CUNIT.GT.0) GOTO 198
C     END OF INPUT FILE HAS BEEN REACHED...
         IF (INTOUT(1).EQ.0 .OR. IRC.NE.0) THEN
            DONE=.TRUE.
            FIRST=.TRUE.
         END IF
 205     CONTINUE
C
      ENDDO
C
      IF (.NOT.SUCC .AND. IRC.EQ.0) IRC=104
C
CSE 998  FORMAT(X,A8,I5,A3,A40)
C
      IF (BDEB) WRITE(*,*) MYNAME,'Debug: Routine ends.',IRC
C
      RETURN
      END

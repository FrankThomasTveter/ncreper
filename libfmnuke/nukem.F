      SUBROUTINE NUKEM(KODE,UNITI,HDR250,INTOUT,DAT250,
     &     ENDOFF,IRC)
C     ***********************************************************************
C     +                                                                     *
C     INTERFACE SUBROUTINE TO FORTRAN-PROGRAM INPUT FILES                   *
C     +                                                                     *
CI    KODE (I*4) = INITIALISATION CODE,                                     *
C     +     KODE=-1 REINITIALISES                                           *
C     +     KODE=0 READ A DATA-BODY ITEM                                    *
C     +     KODE=-1 TERMINATE READING                                       *
CI    UNITI (I*4) = UNIT NUMBER FOR INPUT FILE (5=STANDARD INPUT)           *
CI    HDR250(250) (C*250) = HEADERS TO SEARCH FOR IN INPUT FILE             *
C     +   HEADER MUST BE ON FORM 'THIS IS THE HEADER [N]P', WHERE N         *
C     +   IS THE NUMBER OF LINES INCLUDED IN DATA LINE (CONCATENATED).      *
C     +   IF 'N'='*' THEN ONE DATA LINE IS RETURNED PER LINE IN THE FILE    *
C     +   (NOTE: THE HEADER IS IN THIS CASE ONLY NEEDED ONCE).              *
C     +   A 'P' AT THE END IMPLIES THAT THE DATA SHOULD BE PREPROCESSED.    *
CO    INTOUT (I*4) = HEADER ID (LOC IN ARRAY) FOR DATA FOUND IN FILE *
CO    DAT250 (C*250) = DATA BODY                                     *
C     +      II=INTOUT(JJ) => HEADER HDR250(II) HAD DATABODY DAT250(JJ)     *
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
C     INTEGER IRC,KODE,UNITI                                              *
C     CHARACTER*250 HDR250(250),DAT250                                      *
C     INTEGER INTOUT                                                      *
C     LOGICAL ENDOFF                                                      *
C     UNITI=5                                                               *
C     HDR250(1)='THIS IS A HEADER : [*]P'                                   *
C     KODE=-1                                                               *
C     CALL NUKEM(KODE,UNITI,HDR250,INTOUT,DAT250,ENDOFF,IRC)                *
C     KODE=0                                                                *
C     CALL NUKEM(KODE,UNITI,HDR250,INTOUT,DAT250,ENDOFF,IRC)                *
C     KODE=+1                                                               *
C     CALL NUKEM(KODE,UNITI,HDR250,INTOUT,DAT250,ENDOFF,IRC)                *
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
C     | S.EASTWOOD   | 10/12/07 | MINOR CHANGE IN DECLARATION ORDER.      | *
C     |              |          |                                         | *
C     --------------------------------------------------------------------- *
C     ***********************************************************************
C
      IMPLICIT NONE
      SAVE
C
      INTEGER MAXHDR
      PARAMETER (MAXHDR=100)
C
C     INTERFACE VARIABLES
C
      INTEGER KODE,UNITI
      CHARACTER*250 HDR250(MAXHDR),DAT250
      INTEGER INTOUT,IRC
      LOGICAL ENDOFF
C
C     INTERNAL VARIABLES
C
      LOGICAL LFLDAT(MAXHDR) ! 
      INTEGER NROPT, MAXGRP
      PARAMETER (NROPT=5,MAXGRP=250)
      LOGICAL PPDAT(NROPT,MAXHDR)
      LOGICAL ACTHDR(MAXHDR)
      INTEGER NRLEN(MAXHDR),CODHDR(MAXHDR),NRHDR,CURHDR
      INTEGER HDRGRP(4,MAXGRP),NRGRP
      LOGICAL LFLGRP(MAXGRP)   ! is group set?
      LOGICAL IENDOFF
C
      INTEGER MAXCOM,NRCOM,
     &     LENGTH,LENB,LENF,MEAN
      PARAMETER (MAXCOM=14)
      CHARACTER*20 COM20(MAXCOM),HLP20(MAXCOM)
      CHARACTER*250 NUKEHEAD
      CHARACTER*1000 BUFF1000,PATHF
      EXTERNAL FTUNIT,LENGTH
      INTEGER CMLEN(MAXCOM),system,putenv,FTUNIT
      EXTERNAL system,putenv
C
      INTEGER    LINE,CNT
C
      LOGICAL BDONE,PROCESS
C
C     AUXILIARY VARIABLES
C
      INTEGER LDUMP
      INTEGER MAXVAR,IFLVL
      PARAMETER (MAXVAR=50)
      INTEGER VARLEN(MAXVAR),NRVAR
      CHARACTER*250 NAMVAR(MAXVAR),VALVAR(MAXVAR)
C
C     UNIT CONTROLLERS
C
      INTEGER MAXUNIT
      PARAMETER (MAXUNIT=10)
      INTEGER UNITA(MAXUNIT),CUNIT,LUNIT(MAXUNIT),II
      CHARACTER*250 FILES(MAXUNIT)
C
C     EQUAL SIGN IDENTIFIERS
C
      CHARACTER*8 MYNAME
      DATA MYNAME /'NUKEM'/
      LOGICAL BDEB,ACTIVE,INITIALISED
      DATA ACTIVE /.FALSE./
      DATA INITIALISED /.FALSE./
C
C     Debug System.
C
      IF (.NOT.ACTIVE) CALL DEBUG(MYNAME,BDEB,ACTIVE)
C
      IF (BDEB) WRITE(*,*) MYNAME,'Debug: Routine starts.',IRC
C
      INTOUT=0
C
      IF (KODE.EQ.-1) THEN
         IF (INITIALISED) THEN
            WRITE(*,*) MYNAME,'ROUTINE IS ALLREADY INITIALISED'
            IRC=42
            RETURN
         END IF
C
C     INITIALISE
C
         LINE=-1
         KODE=0
         ENDOFF=.FALSE.
         NRGRP=0
         IFLVL=0
         DO II=1,MAXGRP
            LFLGRP(II)=.FALSE.
         END DO
         LDUMP=0
C
C     C1=CHAR(0)
C
C     DEFINE COMMANDS
C
         CALL GETCOM(MAXCOM,NRCOM,COM20,HLP20,CMLEN,IRC)
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from GETCOM.',IRC
            RETURN
         END IF
C     
C     ANALYZE HEADERS
C
         CALL ANAHDR(MAXHDR,NRHDR,HDR250,ACTHDR,NRLEN,
     &        CODHDR,PPDAT,MAXGRP,NRGRP,HDRGRP,
     &        IRC)
C
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from ANAHDR.',IRC
            RETURN
         END IF
C
         DO II=1,NRHDR
            LFLDAT(II)=.FALSE.
         END DO
C
C     GET FIRST LINE
C
         CNT=0
         CUNIT=1
         UNITA(CUNIT)=UNITI
         FILES(CUNIT)='input file'
         CALL CHOP(FILES(CUNIT),250)
         LENF=LENGTH(FILES(CUNIT),250,10)
C
         IENDOFF=.FALSE.
         CALL GETLINE(MAXUNIT,UNITA,FILES,LUNIT,CUNIT,
     &        CNT,MAXHDR,CODHDR,LINE,
     &        BUFF1000,LENB,IENDOFF,IRC)
C
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from GETLINE.',IRC
            RETURN
         END IF
C
         IF (IENDOFF) THEN
            IRC=266
            WRITE(*,*) MYNAME,
     &           'Use the "?" command in the input file '//
     &           'to list the possible headers.'
            WRITE(*,'(A)') '#'
            WRITE(*,'(A)') ' ?'
            WRITE(*,'(A)') '#'
            WRITE(*,*) MYNAME,'Input file is empty.',IRC
            RETURN
         END IF
C
         INITIALISED = .TRUE.
C
      ELSE IF (KODE.EQ.0) THEN  ! analyse next line
C
         IF (.NOT.INITIALISED) THEN
            WRITE(*,*) MYNAME,'ROUTINE NEEDS TO BE INITIALISED'
            IRC=44
            RETURN
         END IF
C
         INTOUT=0
         DAT250=''
         BDONE=IENDOFF
         DO WHILE (.NOT.BDONE)
C
            PROCESS=.TRUE.
C
C     CHECK FOR COMMAND
C
            CALL QCOM(MAXCOM,NRCOM,COM20,HLP20,CMLEN,PPDAT,
     &           MAXUNIT,UNITA,FILES,LUNIT,CUNIT,
     &           MAXHDR,NRHDR,HDR250,ACTHDR,NRLEN,
     &           MAXVAR,NRVAR,NAMVAR,VALVAR,VARLEN,IFLVL,
     &           MAXGRP,NRGRP,HDRGRP,LFLGRP,
     &           BUFF1000,LENB,PROCESS,BDONE,IENDOFF,LDUMP,IRC)
C
            IF (IRC.NE.0) THEN
               WRITE(*,*) MYNAME,'Error return from QCOM.',
     &              IRC
               RETURN
            END IF
C     
C     CHECK IF HAVE A HEADER
C
            IF (PROCESS) THEN
               CALL QHDR(MAXUNIT,FILES,LUNIT,CUNIT,
     &              MAXHDR,NRHDR,HDR250,ACTHDR,NRLEN,
     &              LFLDAT,INTOUT,CNT,CODHDR,LINE,
     &              BUFF1000,LENB,PROCESS,BDONE,IRC)
C
               IF (IRC.NE.0) THEN
                  WRITE(*,*) MYNAME,'Error return from QHDR.',
     &                 IRC
                  RETURN
               END IF
C
            END IF
C
C     CHECK IF WE SHOULD STORE AS A DATA-BODY
C
            IF (PROCESS) THEN
               CALL QBDY(MAXUNIT,FILES,LUNIT,CUNIT,
     &              MAXHDR,HDR250,NRLEN,PPDAT,INTOUT,
     &              CNT,CODHDR,LINE,DAT250,
     &              MAXVAR,NRVAR,NAMVAR,VALVAR,VARLEN,
     &              BUFF1000,LENB,PROCESS,BDONE,IRC)
C
               IF (IRC.NE.0) THEN
                  WRITE(*,*) MYNAME,'Error return from QBDY.',
     &                 IRC
                  RETURN
               END IF
C
            END IF
C
C     WRITE ERROR MESSAGE IF WE WERE NOT ABLE TO PROCESS LINE
C
            IF (PROCESS) THEN
               
               WRITE(*,*) MYNAME,
     &              '"'//FILES(CUNIT)(1:LENF)//'", line ',
     &              LUNIT(CUNIT),': Error: Unable to interpret:',
     &              BUFF1000(1:LENB)
               IRC=999
               RETURN
            END IF
C
C     GET NEXT LINE IF NOT EOF
C
            IF (.NOT.IENDOFF) THEN ! we may have encountered a "stop" statement
               CALL GETLINE(MAXUNIT,UNITA,FILES,LUNIT,CUNIT,
     &              CNT,MAXHDR,CODHDR,LINE,
     &              BUFF1000,LENB,IENDOFF,IRC)
C
               IF (IRC.NE.0) THEN
                  WRITE(*,*) MYNAME,'Error return from GETLINE.',
     &                 IRC
                  RETURN
               END IF
C
            END IF
C
            IF (IENDOFF) THEN   ! ok, wrap it up...
               BDONE=.TRUE.
               IF (LINE.NE.-1) THEN
                  IF (CNT.EQ.0 .AND. CODHDR(LINE).EQ.-1) THEN
                     LINE=-1
                     INTOUT=0
                  END IF
               ELSE
                  INTOUT=0
               END IF
            END IF
C     
         END DO
C
C     REPORT EOF IF WE ARE RETURNING DATA...
C
         ENDOFF=IENDOFF
C
         IF (INTOUT.NE.0) THEN
            LFLDAT(INTOUT)=.TRUE.
         END IF
C
         IF (INTOUT.NE.0 .AND. LDUMP.NE.0) THEN
            WRITE(LDUMP,*) HDR250(INTOUT)(1:NRLEN(INTOUT))
            IF (CODHDR(INTOUT).NE.0) THEN
               LENF=LENGTH(DAT250,250,10)
               WRITE(LDUMP,*) " "//DAT250(1:LENF)
            END IF
         END IF
C
      ELSE IF (KODE.EQ.1) THEN
         INITIALISED =.FALSE.
C     
C     CHECK THAT READING TERMINATED CORRECTLY
C     
C     
C     CHECK FOR EOF INSIDE FIXED FORMAT DATA-BODY...
C     
         IF (LINE.NE.-1) THEN
            IF (CODHDR(LINE).GT.0 .AND. CNT.NE.CODHDR(LINE)) THEN
               WRITE(*,*) MYNAME,
     &              'End of file interrupts data body.'
               IRC=345
               RETURN
            END IF
         END IF
C     
C     CHECK FOR NESTING ERROR IN IF-EXPRESSIONS...
C     
         IF (IFLVL.NE.0) THEN
            WRITE(*,*) MYNAME,
     &           'End of file interrupts if test.'
            IRC=345
            RETURN
         END IF
C     
C     CHECK FOR MISSING HEADERS...
C     
         CALL ANALFL(MAXHDR,NRHDR,HDR250,
     &        ACTHDR,LFLDAT,NRLEN,
     &        MAXGRP,NRGRP,HDRGRP,
     &        IRC)
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from ANALFL.',
     &           IRC
            RETURN
         END IF
C
         IF (LDUMP.NE.0) THEN
            CLOSE(UNIT=LDUMP,IOSTAT=IRC)
            IRC=0
         END IF
C     
      END IF
C
      IF (BDEB) WRITE(*,*) MYNAME,'Debug: Routine ends.',IRC
C
      RETURN
      END

      SUBROUTINE IOAAPPY_ASCII(FN250,KODE,
     +     TIME,LAT,LON,
     +     SHEIGHT,SATZEN,SATAZI,SATNAME,INSTRUMENT,
     +     NRCHAN,MAXCH,BTEMPS,OBSFAIL,
     +     BDONE,IRC)
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
C     I FN250 = NAME FOR DATA FILE
C     I  KODE = -1 OPEN DATA FILE TO READ, AND IDENTIFY SATNAME
C     +       =-10 OPEN DATA FILE TO WRITE
C     +       =  0 READ NEXT DATA ITEM
C     +       = 10 WRITE NEXT DATA ITEM
C     +       =  1 CLOSE DATA FILE (READ)
C     +       =  100 CLOSE DATA FILE (WRITE)
C     
C     IO TIME(6)  TIME OF OBS/PROFILE
C     +  TIME(1) = YEAR (YYYY)
C     +  TIME(2) = MONTH
C     +  TIME(3) = DAY
C     +  TIME(4) = HOUR
C     +  TIME(5) = MINUTE
C     +  TIME(6) = SEC * 10**4
C     O LAT = LATITUDE OF OBS
C     O LON = LONGTITUDE OF OBS
C     O SHEIGHT = SATELLITE HEIGHT
C     O SATZEN  = SATELLITE ZENITH ANGLE
C     O SATAZI  = SATELLITE AZIMUTH ANGLE
C     O SATNAM  = SATELLITE NAME (e.g. "noaa15")
C     O INSTRUMENT  = INSTRUMENT (e.g. "amsu-a"
C     O NRCHAN  = NUMBER OF SATELLITE CHANNELS
C     O MAXCH  = MAX VALUE FOR NRCHAN (PARAMETER)
C     O BTEMPS(MAXCH) BRIGHTNESS TEMPERATURE FOR EACH CHANNEL (K)
C     O OBSFAIL=FLAG SET WHEN OBSERVATION HAS ERRORS
C     O BDONE=FLAG SET WHEN LAST RECORD IS READ
C     O IRC=ERROR RETURN CODE, 0 = NO ERROR DETECTED
C     
      IMPLICIT NONE
      SAVE
C     
      CHARACTER*250 FN250
      INTEGER    KODE
      INTEGER  TIME(6)
      REAL   LAT
      REAL   LON
      REAL   SHEIGHT
      REAL   SATZEN
      REAL   SATAZI
      CHARACTER*8 SATNAME
      CHARACTER*12 INSTRUMENT
      INTEGER  NRCHAN,MAXCH
      REAL   BTEMPS(MAXCH)
      INTEGER    OBSFAIL,IRC
C     
      LOGICAL  OPENED,BDONE
      CHARACTER*100 FMT100,HDR100
      CHARACTER*14 MYNAME
      INTEGER  RECORD,IKODE,IIKODE,CREC,II,KK
C     
      INTEGER  LENA,LENS,LENGTH
      EXTERNAL LENGTH
C
      LOGICAL  OPENEDR, OPENEDW, FIRSTW
      EXTERNAL FTUNIT
      character*8 date8
      character*10 time10
      REAL   VERSION
      CHARACTER*24 IDB24,ID24
      INTEGER  UNITR, UNITW, FTUNIT,LENI,LENF,
     &     INRCHAN,INRLEV,INRPRED,IRC2,
     &     IVERSION
C     
      DATA MYNAME /'IOAAPPY_ASCII'/
      DATA OPENED /.FALSE./
C
      DATA OPENEDR, OPENEDW /.FALSE.,.FALSE./
C
!----------------------------------------
C     CHECK IF WE HAVE TO READ A NEW RECORD
!----------------------------------------
C
      IF (KODE.EQ.-1) THEN
C     ...check if a file is allready opened...
         IF (OPENEDR) THEN
            WRITE(*,*) MYNAME,'FILE IS ALLREADY OPENED'
            IRC=94
            RETURN
         END IF
C
         UNITR=FTUNIT(IRC)
         IF (IRC.NE.0) THEN
	    WRITE(*,*) MYNAME,'cannot open', irc
	    RETURN
	 ENDIF
C
         LENF=LENGTH(FN250,250,10)
         OPEN ( UNIT=UNITR, STATUS="OLD", FORM="FORMATTED",
     &        ACCESS="SEQUENTIAL",
     &        IOSTAT=IRC, FILE=FN250 )
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'UNABLE TO OPEN:'//FN250(1:LENF)
            RETURN
         END IF
C ---------------
C     READ HEADER
C ---------------
         WRITE(*,*) MYNAME,'Reading file:',FN250(1:LENF)
         READ(UNITR,'(A100)',IOSTAT=IRC) HDR100
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'UNABLE TO READ HEADER:'//FN250(1:LENF)
            RETURN
         END IF
         IF (HDR100(1:11).NE.'AAPP-ASCII ') THEN
            WRITE(*,*) MYNAME,' ATTEMPT TO USE FILE OTHER'//
     &           ' THAN AAPP-ASCII.'
            IRC=996
            RETURN
         END IF
         READ(HDR100(16:21),*,IOSTAT=IRC) VERSION
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'UNABLE TO READ NRCHAN:'//FN250(1:LENF)
            RETURN
         END IF
         IVERSION=INT(VERSION*100)
         IF (HDR100(12:15).NE.'Vrs:'.OR.IVERSION.NE.100) THEN
            WRITE(*,*) MYNAME,' UNKNOWN FILE VERSION:'//
     &           HDR100(12:21)
            IRC=997
            RETURN
         END IF
         ID24=HDR100(22:46)
         WRITE(*,*) MYNAME,'File Id:',ID24
         READ(HDR100(47:49),'(I2,X)',IOSTAT=IRC) 
     &        INRCHAN
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'UNABLE TO READ NRCHAN:'//FN250(1:LENF)
            RETURN
         END IF
         NRCHAN = INRCHAN
C
         WRITE(FMT100,'(A,I2.2,A)') 
     +        '(I4.4,X,4(I2.2,x),I7.7,X,2(F9.4,X),'//
     +        'F7.1,X,F7.2,X,F7.2,X,I5,X,',INRCHAN,'(F7.2))'
C     
         SATNAME=HDR100(50:57)
         CALL CHOP (SATNAME,8)
         LENS=LENGTH(SATNAME,8,4)
C
         INSTRUMENT=HDR100(59:70)
         CALL CHOP (INSTRUMENT,12)
         LENA=LENGTH(INSTRUMENT,12,4)
C
         WRITE(*,*) MYNAME,'Satellite:'//SATNAME(1:LENS)//
     &        ', Instrument:'//INSTRUMENT(1:LENA)
C
         IF (INRCHAN.GT.MAXCH) THEN
            WRITE(*,*) MYNAME,'Too many channels in file.',INRCHAN
            IRC=924
            RETURN
         END IF
C
         OPENEDR=.TRUE.
      ELSE IF (KODE.EQ.-10) THEN
C
!----------------------
!     OPEN FOR WRITING
!----------------------      
C     ...check if a file is allready opened...
         IF (OPENEDW) THEN
            WRITE(*,*) MYNAME,'FILE IS ALLREADY OPENED'
            IRC=95
            RETURN
         END IF
C
         FIRSTW=.TRUE.
         OPENEDW=.TRUE.
      ELSE IF (KODE.EQ.0) THEN
C
!-------------------------
C        READ FILE RECORDS
!-------------------------
C
         IF (.NOT.OPENEDR) THEN
            WRITE(*,*) MYNAME,'OPEN FILE BEFORE READING IT'
            IRC=96
            RETURN
         END IF
C
         READ(UNIT=UNITR,FMT=FMT100,IOSTAT=IRC2)
     +        (TIME(II),II=1,6),LAT,LON,
     +        SHEIGHT,SATZEN,SATAZI,OBSFAIL,
     +        (BTEMPS(II),II=1,INRCHAN)
         IF (IRC2.NE.0) THEN
            BDONE=.TRUE.
         END IF
C
      ELSE IF (KODE.EQ.10) THEN
C
!-------------------------
!        WRITE FILE RECORDS
!-------------------------
C
         IF (.NOT.OPENEDW) THEN
            WRITE(*,*) MYNAME,'OPEN FILE BEFORE WRITING TO IT'
            IRC=97
            RETURN
         END IF
C
C
         IF (FIRSTW) THEN
            UNITW=FTUNIT(IRC)
            IF (IRC.NE.0) THEN
	        WRITE(*,*) MYNAME,'cannot open to write',irc
	        RETURN
	    ENDIF
C     
C     OPEN ( UNIT=UNITW, STATUS="NEW", FORM="FORMATTED",
            LENF=LENGTH(FN250,250,10)
               OPEN ( UNIT=UNITW, STATUS="UNKNOWN", FORM="FORMATTED",
     &              ACCESS="SEQUENTIAL",
     &              IOSTAT=IRC, FILE=FN250 )
            IF (IRC.NE.0) THEN
               WRITE(*,*) MYNAME,'Unable to open:'//FN250(1:LENF)
               RETURN
            END IF
C     
C     WRITE HEADER...
C     
            WRITE(*,*) MYNAME,'Writing to file:',FN250(1:LENF)
            WRITE(HDR100(1:11),'(A11)')  'AAPP-ASCII '
            WRITE(HDR100(12:21),'(A10)') 'Vrs:  1.00'
            ID24=''             ! give new file a new id
            IDB24=ID24
            CALL CHOP(IDB24,24)
            LENI=LENGTH(IDB24,24,20)
            IF (LENI.NE.0) THEN
               HDR100(22:46)=ID24
            ELSE
               call date_and_time(DATE=date8,TIME=time10)
               WRITE(HDR100(22:46),'(A8,X,A10,5X)') date8,time10
            END IF
            INRCHAN = NRCHAN
C     
         WRITE(FMT100,'(A,I2.2,A)') 
     +        '(I4.4,X,4(I2.2,x),I7.7,X,2(F9.4,X),'//
     +           'F7.1,X,F7.2,X,F7.2,X,I5,X,',INRCHAN,'(F7.2))'
C     
            WRITE(HDR100(47:49),'(I2,X)') 
     &           INRCHAN
            IF (INRCHAN.GT.MAXCH) THEN
               WRITE(*,*) MYNAME,'Attempt to put too many channels'//
     &              ' in file.',INRCHAN
               IRC=924
               RETURN
            END IF
C     
            LENS=LENGTH(SATNAME,8,4)
            LENA=LENGTH(INSTRUMENT,12,4)
            WRITE(HDR100(50:57),'(A)') 
     &           SATNAME(1:LENS)
            WRITE(HDR100(58:58),'(X)')
            WRITE(HDR100(59:70),'(A)') 
     &           INSTRUMENT(1:LENA)
C     
            WRITE(UNITW,'(A)') HDR100(1:70)
            FIRSTW=.FALSE.
         END IF
C   
C    Write data to file
C     
C         WRITE(UNIT=UNITW,FMT=FMT100)
         WRITE(UNIT=UNITW,FMT=FMT100,IOSTAT=IRC)
     +        (TIME(II),II=1,6),LAT,LON,
     +        SHEIGHT,SATZEN,SATAZI,OBSFAIL,
     +        (BTEMPS(II),II=1,INRCHAN)
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'UNABLE TO WRITE TO:'//FN250(1:LENF)
            RETURN
         END IF
C
      ELSE IF (KODE.EQ.1) THEN
C      
         IF (OPENEDR) THEN
            CLOSE (UNIT=UNITR, IOSTAT=IRC)
         END IF
C
         OPENEDR=.FALSE.
C
      ELSE IF (KODE.EQ.100) THEN
C
         IF (OPENEDW) THEN
            CLOSE (UNIT=UNITW, IOSTAT=IRC)
         END IF
C
         OPENEDW=.FALSE.
C
      ELSE
         WRITE(*,*) MYNAME,'UNKNOWN KODE:',KODE
         IRC=999
         RETURN
      END IF
C
C     IDENTIFY THE NEXT OBS
C
C     RETURN OBS
C
      RETURN
      END

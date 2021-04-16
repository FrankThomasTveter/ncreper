      SUBROUTINE IO2AXBY5(FN250,KODE,
     +     SATNAME,INSTRUMENT,ID24,
     +     TIME,LAT,LON,SHEIGHT,SATZEN,SATAZI,
     +     MAXCH,NRCHAN,BTEMPS,
     +     LXB,BTEMPSC,BTEMPSB,NRPRED,PRED,
     +     GPRES,GWIND,GWINDIR,GTEMP,
     +     STYPE,FRACTL,PCAT,ICEC,
     +     MAXICE,NRICE,ICEP,MAXCL,NRCL,CLP,
     +     MAXLEV,NRLEV,PRES,TEMP,WATV,
     +     IRC)
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
C     I  FN250 = NAME FOR DATA FILE
C     I  KODE = -1 OPEN DATA FILE TO READ, AND IDENTIFY SATNAME
C     +       =-10 OPEN DATA FILE TO WRITE
C     +       =  0 READ NEXT DATA ITEM
C     +       = 10 WRITE NEXT DATA ITEM
C     +       =  1 CLOSE DATA FILE
C
C     OBS DATA
C
C     IO SATNAM  = SATELLITE ID
C     IO INSTRUMENT = INSTRUMENT CONFIGURATION
C     O ID24 = THE TIME STAMP OF THE FILE (RETURNED AT KODE=-1)
C     IO TIME(6)  TIME OF PROFILE
C     +  TIME(1) = YEAR (YYYY)
C     +  TIME(2) = MONTH
C     +  TIME(3) = DAY
C     +  TIME(4) = HOUR
C     +  TIME(5) = MINUTE
C     +  TIME(6) = SEC * 10^4
C     IO LAT = LATITUDE OF OBS/PROFILE
C     IO LON = LATITUDE OF OBS/PROFILE
C     IO SHEIGHT = SATELLITE HEIGHT
C     IO SATZEN  = SATELLITE ZENITH ANGLE
C     IO NRCHAN  = NUMBER OF BRIGHTNESS TEMPERATURE CHANNELS (MAX 100)
C     IO BTEMPS(NRCHAN) OBS BRIGHTNESS TEMPERATURE FOR EACH CHANNEL
C
C     CALCULATED BRIGHTNESS TEMPERATURE FROM BACKGROUND FIELD
C
C     IO LXB = BACKGROUND DATA HAS BEEN INCLUDED
C     IO BTEMPSC(NRCHAN) CALC BRIGHTNESS TEMPERATURE FOR EACH CHANNEL
C     IO BTEMPSB(NRCHAN) BIAS CORRECTED OBS BRIGHTNESS 
C     +  TEMPERATURE FOR EACH CHANNEL
C
C     PROFILE DATA
C
C     IO GPRES = GROUND PRESSURE
C     IO GWIND = WIND VELOCITY AT 10 METERS       (?)
C     IO GTEMP(2) = GROUND TEMPERATURE (0m, 2m)
C     IO STYPE = TYPE OF SURFACE (0=SEA, 1=LAND)
C     IO FRACTL = FRACTION OF LAND (0..1)
C     IO PCAT = SEA SURFACE TYPE CATEGORY (1=SEA ICE, 2=OPEN WATER)
C     IO ICEC  = ICECOVER (0..1)
C     IO NRICE  = NUMBER OF ICE PARAMETERS
C     IO ICEP(NRICE) = ICE PARAMETERS
C     IO NRCL = NUMBER OF CLOUD PARAMETERS
C     IO CLP(NRCL) = CLOUD PARAMETERS
C     IO NRLEV  = NUMBER OF LEVELS IN THE PROFILE (MAX 100)
C     IO PRES(NRLEV) = PRESSURE AT EACH LEVEL
C     IO TEMP(NRLEV) = AIR TEMPERATURE AT EACH LEVEL
C     IO WATV(NRLEV) = WATER VAPOUR AT EACH LEVEL
C
C     O IRC=ERROR RETURN CODE
C     +     = 0, NORMAL END
C
      IMPLICIT NONE
      SAVE
C
      INTEGER  MAXCH,MAXLEV,MAXICE,MAXCL
C     PARAMETER (MAXCH = 35, MAXLEV = 40,MAXICE=2,MAXCL=1)
C
      CHARACTER*250 FN250
      INTEGER    KODE
      CHARACTER*8 SATNAME
      CHARACTER*12 INSTRUMENT
      CHARACTER*24 ID24
      INTEGER  TIME(6)
      REAL   LAT
      REAL   LON
      REAL   SHEIGHT
      REAL   SATZEN,SATAZI
      INTEGER  NRCHAN
      REAL   BTEMPS(MAXCH)
      LOGICAL  LXB
      REAL   BTEMPSC(MAXCH)
      REAL   BTEMPSB(MAXCH)
      INTEGER  NRPRED
      REAL   PRED(NRPRED)
      REAL   GPRES
      REAL   GWIND,GWINDIR
      REAL   GTEMP(2)
      INTEGER  STYPE
      REAL   FRACTL
      INTEGER  PCAT
      REAL   ICEC
      INTEGER NRICE
      REAL   ICEP(MAXICE)
      INTEGER NRCL
      REAL   CLP(MAXCL)
      INTEGER  NRLEV
      REAL   PRES(MAXLEV)
      REAL   TEMP(MAXLEV)
      REAL   WATV(MAXLEV)
      INTEGER    IRC
C
      CHARACTER*24 IDB24
      INTEGER  UNITR, UNITW, FTUNIT,
     &     II,INRCHAN,INRLEV,INRPRED,
     &     IRC2,INRICE,INRCL,
     &     IVERSION
      CHARACTER*10 MYNAME
      LOGICAL  OPENEDR, OPENEDW, RLXB,WLXB
      EXTERNAL FTUNIT
      CHARACTER*100 HDR100
      character*8 date8
      character*10 time10
      REAL   VERSION
C
      logical firstw
C
      REAL   KAPPA,RPRES,FACT
      DATA KAPPA /0.286D0/,RPRES/1000.0D0/
C
      INTEGER  LENA,LENS,LENGTH,LENF,LENI,nrrd,nrwr
      EXTERNAL LENGTH
      logical bdeb
      data bdeb /.false./
C
      DATA MYNAME, OPENEDR, OPENEDW /'IO2AXBY5',.FALSE.,.FALSE./
C
C     CHECK IF WE HAVE TO READ A NEW RECORD
C
C      bdeb=.true.
      if (bdeb) write(*,*) myname,'Entered routine.',irc
      IF (KODE.EQ.-1) THEN
C     ...check if a file is allready opened...
         IF (OPENEDR) THEN
            WRITE(*,*) MYNAME,'FILE IS ALLREADY OPENED'
            IRC=94
            RETURN
         END IF
C
         UNITR=FTUNIT(IRC)
         IF (IRC.NE.0) RETURN
C
         LENF=LENGTH(FN250,250,10)
         OPEN ( UNIT=UNITR, STATUS="OLD", FORM="UNFORMATTED",
     &        ACCESS="SEQUENTIAL",
     &        IOSTAT=IRC, FILE=FN250 )
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'UNABLE TO OPEN:'//FN250(1:LENF)
            RETURN
         END IF
C
C     READ HEADER
C
         WRITE(*,*) MYNAME,'Reading file:',FN250(1:LENF)
         READ(UNITR,IOSTAT=IRC) HDR100
         IF (IRC.NE.0) RETURN
         IF (HDR100(1:10).NE.'Lvl:  2a  ') THEN
            WRITE(*,*) MYNAME,' ATTEMPT TO USE FILE OTHER'//
     &           ' THAN LEVEL 2A.'
            IRC=996
            RETURN
         END IF
         READ(HDR100(15:20),*,IOSTAT=IRC2) VERSION
         IVERSION=NINT(VERSION*100.0D0)
         IF (HDR100(11:14).NE.'Vrs:'
     &        .OR.(IVERSION.NE.501
     &        .and.IVERSION.NE.502)) THEN
            WRITE(*,*) MYNAME,' UNKNOWN FILE VERSION:'//
     &           HDR100(11:20)
            IRC=997
            RETURN
         END IF
         if (bdeb) write(*,*) myname,HDR100(11:20),version,iversion
         if (iversion.eq.501) then
            WRITE(*,*) MYNAME,'Setting T0m = T2m.'
         end if
         ID24=HDR100(21:45)
         WRITE(*,*) MYNAME,'File Id:',ID24
         READ(HDR100(46:65),'(5(I3,X))',IOSTAT=IRC) 
     &        INRCHAN,INRLEV,INRPRED,INRICE,INRCL
         NRCHAN = INRCHAN
         IF (IRC.NE.0) then
            write(*,*) myname,'Unable to read nrch/lev etc.'
            RETURN
         end if
         RLXB = (HDR100(92:92).EQ.'B')
         LXB=RLXB
         IF (RLXB) WRITE(*,*) MYNAME,'File includes '//
     &        ' background information.'
         IF (HDR100(91:91).EQ.'R') THEN
            WRITE(*,*) MYNAME,'File has been RE-RUN.'
         END IF
C
         READ(HDR100(70:77),'(A8)',IOSTAT=IRC) 
     &        SATNAME
         CALL CHOP (SATNAME,8)
         LENS=LENGTH(SATNAME,8,4)
C
         READ(HDR100(79:90),'(A12)',IOSTAT=IRC) 
     &        INSTRUMENT
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
         NRLEV = INRLEV
         IF (RLXB .AND. INRLEV.GT.MAXLEV) THEN
            WRITE(*,*) MYNAME,'Too many levels in file.',INRLEV
            IRC=925
            RETURN
         END IF
         IF (RLXB .AND. INRPRED.GT.NRPRED) THEN
            WRITE(*,*) MYNAME,'Wrong number of predictors in file.',
     &           INRPRED,NRPRED
            IRC=934
            RETURN
         END IF
         NRICE=INRICE
         IF (INRICE.GT.MAXICE) THEN
            WRITE(*,*) MYNAME,'Too much ice data in file.',inrice
            IRC=922
            RETURN
         END IF
         NRCL=INRCL
         IF (INRCL.GT.MAXCL) THEN
            WRITE(*,*) MYNAME,'Too much cloud data in file.',inrcl
            IRC=922
            RETURN
         END IF
C
         OPENEDR=.TRUE.
         nrrd=0
      ELSE IF (KODE.EQ.-10) THEN
C     ...check if a file is allready opened...
         IF (OPENEDW) THEN
            WRITE(*,*) MYNAME,'FILE IS ALLREADY OPENED'
            IRC=95
            RETURN
         END IF
C
         firstw=.true.
C
         OPENEDW=.TRUE.
         nrwr=0
      ELSE IF (KODE.EQ.0) THEN
         IF (.NOT.OPENEDR) THEN
            WRITE(*,*) MYNAME,'OPEN FILE BEFORE READING IT'
            IRC=96
            RETURN
         END IF
C
         READ(UNIT=UNITR,IOSTAT=IRC)
     +        (TIME(II),II=1,6),LAT,LON,
     +        SHEIGHT,SATZEN,SATAZI,
     +        (BTEMPS(II),II=1,INRCHAN)
         IF (IRC.NE.0) RETURN
         IF (RLXB) THEN
         READ(UNIT=UNITR,IOSTAT=IRC)
     +        (BTEMPSC(II),II=1,INRCHAN)
         IF (IRC.NE.0) RETURN
         READ(UNIT=UNITR,IOSTAT=IRC)
     +        (BTEMPSB(II),II=1,INRCHAN)
         IF (IRC.NE.0) RETURN
         READ(UNIT=UNITR,IOSTAT=IRC)
     +        (PRED(II),II=1,INRPRED)
         IF (IRC.NE.0) RETURN
         if (iversion.eq.501) then
            READ(UNIT=UNITR,IOSTAT=IRC)
     +           GPRES,GWIND,GWINDIR,GTEMP(2),
     +           STYPE,FRACTL,PCAT,ICEC
            gtemp(1)=gtemp(2) ! set t0m=t2m
         else
            READ(UNIT=UNITR,IOSTAT=IRC)
     +           GPRES,GWIND,GWINDIR,
     +           GTEMP(1),GTEMP(2),
     +           STYPE,FRACTL,PCAT,ICEC
         end if
         if (bdeb) write(*,*) 'R',gtemp(1),gtemp(2),pcat,icec
         IF (IRC.NE.0) RETURN
         READ(UNIT=UNITR,IOSTAT=IRC)
     +        (PRES(II),II=1,INRLEV),
     +        (TEMP(II),II=1,INRLEV),
     +        (WATV(II),II=1,INRLEV),
     +        (icep(II),II=1,inrice),
     +        (clp(II),II=1,inrcl)
         IF (IRC.NE.0) RETURN
         END IF
C     
C     version update
C     
         IF (INRPRED.EQ.6.AND.NRPRED.EQ.7) THEN
            PRED(7)=SATZEN
         END IF
         nrrd=nrrd+1
C
      ELSE IF (KODE.EQ.10) THEN
         IF (.NOT.OPENEDW) THEN
            WRITE(*,*) MYNAME,'OPEN FILE BEFORE WRITING TO IT'
            IRC=97
            RETURN
         END IF
C         
         if (firstw) then
            firstw=.false.
C
         UNITW=FTUNIT(IRC)
         IF (IRC.NE.0) RETURN
C
C     OPEN ( UNIT=UNITW, STATUS="NEW", FORM="UNFORMATTED",
         LENF=LENGTH(FN250,250,10)
         OPEN ( UNIT=UNITW, STATUS="UNKNOWN", FORM="UNFORMATTED",
     &        ACCESS="SEQUENTIAL",
     &        IOSTAT=IRC, FILE=FN250 )
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Unable to open:'//FN250(1:LENF)
            RETURN
         END IF
C     
C     WRITE HEADER...
C
         HDR100=''
         WLXB=LXB
         WRITE(*,*) MYNAME,'Writing to file:',FN250(1:LENF)
         WRITE(HDR100(1:10),'(A10)')  'Lvl:  2a  '
            WRITE(HDR100(11:20),'(A10)') 'Vrs:  5.02'
         IDB24=ID24
         CALL CHOP(IDB24,24)
         LENI=LENGTH(IDB24,24,20)
         IF (LENI.NE.0) THEN
            HDR100(21:45)=ID24
            HDR100(81:81)='R'
         ELSE
            call date_and_time(DATE=date8,TIME=time10)
            WRITE(HDR100(21:45),'(A,X,A,5X)') date8,time10
            HDR100(91:91)=' '
         END IF
         INRCHAN = NRCHAN
         INRLEV = NRLEV
         INRPRED = NRPRED
            INrICE = NRICE
            INRCL = NRCL
            WRITE(HDR100(46:65),'(5(I3,X))',iostat=irc) 
     &           INRCHAN, INRLEV, INRPRED,INRICE,INRCL
            IF (IRC.NE.0) then
               write(*,*) myname,'Unable to write nrch/lev etc.'
               RETURN
            end if
         IF (INRCHAN.GT.MAXCH) THEN
            WRITE(*,*) MYNAME,'Attempt to put too many channels'//
     &           ' in file.',INRCHAN
            IRC=924
            RETURN
         END IF
         IF (WLXB .AND. INRLEV.GT.MAXLEV) THEN
            WRITE(*,*) MYNAME,'Attempts to put too many levels'//
     &           ' in file.',INRLEV
            IRC=925
            RETURN
         END IF
         IF (WLXB .AND. INrICE.GT.MAXICE) THEN
            WRITE(*,*) MYNAME,'Attempts to put too much ice'//
     &           ' in file.',INrICE
            IRC=921
            RETURN
         END IF
         IF (WLXB .AND. INrCL.GT.MAXCL) THEN
            WRITE(*,*) MYNAME,'Attempts to put too much cloud'//
     &           ' in file.',INrCL
            IRC=921
            RETURN
         END IF
         IF (WLXB) THEN
            WRITE(HDR100(92:92),'(A1)') 'B'
         ELSE
            WRITE(HDR100(92:92),'(X)')
         END IF
         WRITE(HDR100(67:67),'(X)')
C     
         WRITE(HDR100(68:69),'(2X)')
         WRITE(HDR100(70:77),'(A8)') 
     &        SATNAME
         WRITE(HDR100(78:78),'(X)')
         WRITE(HDR100(79:90),'(A12)') 
     &        INSTRUMENT
C     
         WRITE(UNITW) HDR100
            
         END IF
C
         WRITE(UNIT=UNITW,IOSTAT=IRC)
     +        (TIME(II),II=1,6),LAT,LON,
     +        SHEIGHT,SATZEN,SATAZI,
     +        (BTEMPS(II),II=1,INRCHAN)
C
         IF (WLXB) THEN
            WRITE(UNIT=UNITW,IOSTAT=IRC)
     +           (BTEMPSC(II),II=1,INRCHAN)
            WRITE(UNIT=UNITW,IOSTAT=IRC)
     +           (BTEMPSB(II),II=1,INRCHAN)
            WRITE(UNIT=UNITW,IOSTAT=IRC)
     +           (PRED(II),II=1,INRPRED)
            WRITE(UNIT=UNITW,IOSTAT=IRC)
     +        GPRES,GWIND,GWINDIR,
     +        GTEMP(1),GTEMP(2),
     +        STYPE,FRACTL,PCAT,ICEC
            WRITE(UNIT=UNITW,IOSTAT=IRC)
     +           (PRES(II),II=1,INRLEV),
     +           (TEMP(II),II=1,INRLEV),
     +           (WATV(II),II=1,INRLEV),
     +        (icep(II),II=1,INRice),
     +        (clp(II),II=1,INrcl)
            if (bdeb) write(*,*) 'W',gtemp(1),gtemp(2),pcat,icec
         END IF
         nrwr=nrwr+1
C
      ELSE IF (KODE.EQ.1.OR.KODE.EQ.100) THEN
C
         IF (OPENEDR) THEN
            CLOSE (UNIT=UNITR, IOSTAT=IRC)
            write(*,*) myname,'Read obs:',nrrd
         END IF
C
         IF (OPENEDW) THEN
            CLOSE (UNIT=UNITW, IOSTAT=IRC)
            write(*,*) myname,'Wrote obs:',nrwr
         END IF
C
         OPENEDR=.FALSE.
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
      if (bdeb) write(*,*) myname,'Done.',irc
C
      RETURN
      END

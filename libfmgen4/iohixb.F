      SUBROUTINE IOHIXB(FN250,MAXF,KODE,
     +     TIME,LAT,LON,SATZEN,SATAZI,DELTA,REH,
     +     GPRES,GWIND,GWINDIR,GTEMP,GWATV,GOZONE,GWATL,
     +     SOILM,FRACTL,ALBEDO,ICEC,
     +     CLCOV,CLTOP,CLBOT,
     +     NRLEV,PRES,TEMP,WATV,OZONE,WATL,
     +     BPRO,IRC)
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
C     I  FN250(MAXF) = LIST OF BACKGROUND DATA FILES
C     I  MAXF = NUMBER OF DATA FILES
C     I  KODE = -1 OPEN DATA FILE AND IDENTIFY SATNAME
C     +       =  0 RETRIEVE NEXT SET OF BACKGROUND INFO
C     +       =  1 CLOSE DATA FILE
C     IO TIME(6)  TIME OF PROFILE
C     +  TIME(1) = YEAR (YYYY)
C     +  TIME(2) = MONTH
C     +  TIME(3) = DAY
C     +  TIME(4) = HOUR
C     +  TIME(5) = MINUTE
C     +  TIME(6) = SEC * 10**4
C     IO LAT = LATITUDE OF PROFILE
C     IO LON = LONGTITUDE OF PROFILE
C     IO SATZEN = SATELLITE ZENITH ANGLE (DEG)
C     IO SATAZI = SATELLITE AZIMUTH ANGLE (DEG)
C     IO DELTA  = SUB-SATELLITE INSTRUMENT RESOLUTION (DEG)
C     IO GPRES  = GROUND PRESSURE (MBAR)
C     IO GWIND  = WIND VELOCITY AT 10 METERS (M/S)
C     IO GTEMP  = GROUND TEMPERATURE (K)
C     IO GWATV  = GROUND WATER VAPOUR (KG/KG)
C     IO GOZONE = GROUND OZONE         (?,NOT USED?)
C     IO GWATL  = GROUND LIQUID WATER  (?,NOT USED?)
C     IO SOILM  = SOIL MOISTURE (KG/M**2, NOT USED?)
C     IO HEIGHT = LAND HEIGHT (M, NOT USED?)
C     IO FRACTL = FRACTION OF LAND (0..1)
C     IO ALBEDO = ALBEDO (0..1) 
C     IO ICEC  = ICECOVER (0..1)
C     IO CLCOV  = IR CLOUD COVER (0..1) (NOT USED)
C     IO CLTOP  = CLOUD TOP PRESSURE (MBAR, NOT USED)
C     IO CLBOT  = CLOUD BOTTOM PRESSURE (MBAR, NOT USED)
C     IO NRLEV  = NUMBER OF LEVELS IN THE PROFILE (<=MAXLEV)
C     IO PRES(NRLEV) = PRESSURE AT EACH LEVEL (MBAR)
C     IO TEMP(NRLEV) = AIR TEMPERATURE AT EACH LEVEL (K)
C     IO WATV(NRLEV) = WATER VAPOUR AT EACH LEVEL (KG/KG)
C     IO OZONE(NRLEV) = OZONE AT EACH LEVEL (?)
C     IO WATL(NRLEV) = LIQUID WATER AT EACH LEVEL (KG/M**2)
C     O BPRO=LOGICAL*4 FLAG THAT TELLS IF FIELDS ARE OK
C     O IRC=ERROR RETURN CODE, 0 = NO ERROR DETECTED
C
      IMPLICIT NONE
      SAVE
C
      INTEGER  MAXLEV
      PARAMETER (MAXLEV = 31)
C
C     INTERFACE PARAMETERS...
C
      INTEGER  MAXF
      CHARACTER*250 FN250(MAXF)
      INTEGER    KODE
      INTEGER  TIME(6)
      REAL   LAT
      REAL   LON
      REAL   SATZEN,SATAZI,DELTA,REH ! footprint geometry
      REAL   GPRES
      REAL   GWIND,GWINDIR
      REAL   GTEMP
      REAL   GWATV
      REAL   GOZONE
      REAL   GWATL,
     +     SOILM,HEIGHT,FRACTL,ALBEDO,ICEC
      REAL   CLCOV
      REAL   CLTOP
      REAL   CLBOT
      INTEGER  NRLEV
      REAL   PRES(MAXLEV)
      REAL   TEMP(MAXLEV)
      REAL   WATV(MAXLEV)
      REAL   OZONE(MAXLEV)
      REAL   WATL(MAXLEV)
      LOGICAL  BPRO
      INTEGER  IRC
C
C     OTHER DECLARATIONS...
C
      real cosdeg,sindeg,tandeg,atan2deg
      external cosdeg,sindeg,tandeg,atan2deg
      INTEGER  NX,NY,INX,INY,ILEV,NPAR,NTID,NPARS,
     &     LDATA,IPACK,IERROR,IKODE
      INTEGER  NT,LENF,NPAR2, LENGTH
      EXTERNAL LENGTH
C      PARAMETER ( INX=188, INY=152,ILEV=31 )
      PARAMETER ( INX=468, INY=378,ILEV=31 )
      PARAMETER ( NPAR = 4,NPARS = 6, NTID=2, NPAR2=4 )
      PARAMETER ( LDATA=20+INX*INY+50)
      INTEGER*2 IDATA(LDATA),IN(16)
      INTEGER  FTUNIT,ITIME(6),II,JJ,KK
      INTEGER    IUNIT,IPAR(NPAR),IPARS(NPARS),IPAR2(NPAR2)
      REAL*4 FELT(INX,INY,ILEV,NPAR,NTID), SURFACE(INX,INY,NPARS,NTID)
      REAL*4 FIXDAT(INX,INY,NPAR2)
      REAL   FELT_INT2(ILEV,NPAR,NTID),SURF_INT2(NPARS,NTID),
     &     FIX_INT(NPAR2)
      REAL   FELT_INT(ILEV,NPAR),SURF_INT(NPARS)
      REAL GRID(6)           ! GRID SPECIFICATIONS
C      REAL*4 GRID(6)           ! GRID SPECIFICATIONS
      REAL   KAPPA, A(ILEV), B(ILEV), CURPRES,CURTEMP,RPRES
      REAL   J_FELT, J_OBS     ! JULIAN DAYS
      REAL   SEC, DT, P, R
      INTEGER  IGTYPE         ! GRID TYPE
      CHARACTER*70 FN70
      INTEGER  MAXPOS,NRPOS
      PARAMETER (MAXPOS=10000)
      INTEGER  I,J, FTIME(NTID),POS(2,MAXPOS)
      REAL   X,Y,X1,Y1,X5,X6,X7,X8,
     &     XCM,YCM,XAM,YAM,XBM,YBM,XAMI,YAMI,XBMI,YBMI,
     &     CZ,CA,SA,CL,LATA,LATB,LONA,LONB
      REAL*4 XPOS,YPOS
      CHARACTER*8 MYNAME
      LOGICAL  OPENED, BDONE, FOUND,FIRST
      EXTERNAL FTUNIT
C
      LOGICAL  BDEB,ACTIVE
      DATA ACTIVE /.FALSE./
      DATA MYNAME,OPENED /'IOHIXB',.FALSE./
      DATA FTIME /6,12/
C
C
      INTEGER  CLEV,TYPE,NRN
      PARAMETER ( CLEV = 50 )
      REAL   CPRES(CLEV),COZONE(CLEV)
      REAL   TARGPRES(MAXLEV), TARGOZONE(MAXLEV),
     &     TPRES,TOZONE
C
C      DATA IOZONE /0.022068,0.015525,0.010460,0.0066280,
C     &     0.0078370,0.0086965,0.0099048,0.010763,0.011420,
C     &     0.011519,0.011184,0.0099956,0.0086073,0.0073708,
C     &     0.0062772,0.0033368,0.0024259,0.0017038,
C     &     0.0011434,0.00076103,0.00053896,0.00037990,
C     &     0.00029614,0.00015238,0.00010984,0.000078304,
C     &     0.000071910,0.000069046,0.000067485,
C     &     0.000065300,0.000063913,0.000059832,
C     &     0.000057019,0.000054246,0.000052504,
C     &     0.000047331,0.000042865,0.000038740,
C     &     0.000037067,0.000034964/
C
      IF (.NOT.ACTIVE) CALL DEBUG(MYNAME,BDEB,ACTIVE)
C
      IF (BDEB) WRITE(*,*) MYNAME,'Debug: Routine starts',KODE
C
C     CHECK IF WE HAVE TO READ A NEW RECORD
C
      IF (KODE.EQ.-1) THEN
C     ...CHECK IF A FILE IS ALLREADY OPENED...
         IF (OPENED) THEN
            WRITE(*,*) MYNAME,'FILE IS ALLREADY OPENED'
            IRC=48
            RETURN
         END IF
C
         NRLEV=ILEV
C
         IPARS(1)=8  ! PS
         IPARS(2)=30 ! TS
         IPARS(3)=31 ! T2M
         IPARS(4)=32 ! RH2M
         IPARS(5)=33 ! U10
         IPARS(6)=34 ! V10
C
         IPAR(1)=2  ! U
         IPAR(2)=3  ! V
         IPAR(3)=9  ! Q
         IPAR(4)=18 ! THETA
C
         IPAR2(1)=181  ! FRACTION OF LAND
         IPAR2(2)=101  ! TOPOGRAFI
C     IPAR2(3)=184  ! ALBEDO
         IPAR2(3)=191  ! FRACTION OF ICE
         IPAR2(4)=103  ! SST
C
         FIRST=.TRUE.
         OPENED=.TRUE.
      ELSE IF (KODE.EQ.0) THEN
C     ...CHECK IF FILE HAS BEEN OPENED...
         IF (.NOT.OPENED) THEN
            WRITE(*,*) MYNAME,'OPEN FILE BEFORE READING IT'
            IRC=55
            RETURN
         END IF
C
         IF (FIRST) THEN
C
            NRN=0
C
C     GET TIME OF OBS
C
            SEC = TIME(6)*1.0D-4
            CALL JD2000 (J_OBS, TIME(1),TIME(2),TIME(3),
     &           TIME(4),TIME(5),SEC) 
C     
C     OPEN FILE AND GET HEADER DATA
C     
            IPACK=1
C     
            IN( 1) = 88         ! PRODUSENT
C            IN( 2) = 2006       ! GRID-NUMBER
            IN( 2) = 2001       ! GRID-NUMBER
            IN( 3) = -32767
            IN( 4) = -32767
            IN( 5) = -32767
            IN( 9) = 2          ! TYPE,1=ANALYSIS,2=FORECAST
C     
            FOUND=.FALSE.
            II = 1
            BDONE = (II.GT.MAXF)
            DO WHILE (.NOT.BDONE)
C     
               IUNIT = FTUNIT(IRC)
               IF (IRC.NE.0) RETURN
C     
C  
               FN70 = FN250(II)
               CALL CHOP(FN70,70)
               LENF = LENGTH(FN70,70,10)
               CALL MRFELT(1,FN70,IUNIT,IN,IPACK,1,1.0,1.0,32,
     &              IDATA,IERROR)
               IF(IERROR.NE.0)THEN
                  WRITE(*,*) MYNAME, 'UNABLE TO OPEN FELT-FILE',
     &                 FN70(1:LENF)
               END IF
C     
               ITIME(1)=IDATA(5)
               ITIME(2)=IDATA(6)/100
               ITIME(3)=IDATA(6)-IDATA(6)/100*100
               ITIME(4)=IDATA(7)/100
               ITIME(5)=0
               ITIME(6)=0
               SEC = 0
               CALL JD2000 (J_FELT, ITIME(1),ITIME(2),ITIME(3),
     &              ITIME(4),ITIME(5),SEC)
C     
               IF (J_OBS-J_FELT .LT. 0.5 .AND. J_OBS-J_FELT 
     &              .GE. 0.25 ) THEN
                  FOUND=.TRUE.
               ELSE
                  II = II+1
C     CLOSE FILE
                  CALL MRFELT(3,FN70,IUNIT,IN,IPACK,
     &                 1,1.0,1.0,1,1,IERROR)
                  IERROR=0
               ENDIF
               
               BDONE=(II.GT.MAXF.OR.FOUND)
            END DO
C     
            IF (.NOT. FOUND) THEN
               WRITE(*,*) MYNAME, 'UNABLE TO FIND CORRECT FELT-FILE'
               IRC = 49
               RETURN
            ENDIF
C     
            WRITE(*,*) MYNAME,'Using felt file:'//FN70(1:LENF)
C     
C     PARAMETER LOOP,LEVEL DATA
            DO NT = 1,NTID      ! +6, +12
               IN(10) = FTIME(NT) ! FORECAST LENGHT
               IN(11) = 10      ! VERTICAL COORDINATE,1:PRESSURE,10:ETA
               DO JJ = 1,ILEV
C     IN(13) = 100
                  IN(13) = JJ
                  IN(14) = -32767 ! LEVEL 2
                  DO II=1,NPAR
C     IN(12)=308
                     IN(12)=IPAR(II) ! PARAMETER: U,V,Q,TETA, ! CL
                     CALL MRFELT(2,FN70,IUNIT,IN,IPACK,INX*INY,
     +                    FELT(1,1,JJ,II,NT),1.0,LDATA,IDATA,IERROR)
C     
                     IF(IERROR.NE.0)THEN
                        WRITE(*,*) MYNAME,'Error return from MRFELT.',
     &                       IERROR
                        IRC=IERROR
                        RETURN
                     END IF
C     
                  END DO
C     
C     ASSUME THAT PARAMETER WAS = 18 FOR CURRENT IDATA...
C     ...READ A AND B PARAMETERS FOR ETA <-> PRESSURE CONVERTION
C     
                  IF (JJ.NE.IDATA(7)) THEN
                     WRITE(*,*) MYNAME,'SYSTEM ERROR',JJ,IDATA(7)
                     IRC=51
                     RETURN
                  END IF
C     
                  A(JJ) = IDATA(8)*0.1D0
                  B(JJ) = IDATA(19)*0.0001D0
C     
               END DO
C     
               IN(11) = 2       ! VERTICAL COORDINATE,1:PRESSURE,10:ETA
C     
               IN(13) = 1000    ! SURFACE LEVEL
               IN(14) = -32767  ! LEVEL 2
C     
C     STOP 'BDONE...'
C     
C     PARAMETER LOOP,SURFACE DATA
               DO II=1,NPARS
                  IN(12)=IPARS(II) ! PARAMETER: U,V,Q,TETA, ! CL
                  CALL MRFELT(2,FN70,IUNIT,IN,IPACK,INX*INY,
     +                 SURFACE(1,1,II,NT),1.0,LDATA,IDATA,IERROR)
                  IF(IERROR.NE.0)THEN
                     IRC=52
                     RETURN
                  END IF
C     
               END DO
               
C     
            END DO
C     SURFACE PARAMETERS
            IN( 9) = 4          ! TYPE: PARAMETER FELT
            IN(10) = 0          ! FORECAST LENGHT
            IN(11) = 2          ! VERTICAL COORDINATE,1:PRESSURE,10:ETA
            IN(13) = 1000       ! SURFACE LEVEL
            IN(14) = -32767     ! LEVEL 2
C     
C     STOP 'BDONE...'
C     
C     PARAMETER LOOP,SURFACE DATA
            DO II=1,NPAR2
               IN(12)=IPAR2(II) ! PARAMETER: LAND FRAC.,TOPO, ALBEDO, SNOW, ICE
               CALL MRFELT(2,FN70,IUNIT,IN,IPACK,INX*INY,
     +              FIXDAT(1,1,II),1.0,LDATA,IDATA,IERROR)
               IF(IERROR.NE.0)THEN
                  IRC=723
                  RETURN
               END IF
C     
            ENDDO
C
C     CLOSE FILE
            CALL MRFELT(3,FN70,IUNIT,IN,IPACK,1,1.0,1.0,1,1,IERROR)
C     
C     FIND GRID PARAMETERS
C     
            CALL GRIDPAR( +1, LDATA,IDATA, IGTYPE,NX,NY,GRID, IERROR)
            IF(IERROR.NE.0)THEN
               IRC=53
               RETURN
            END IF
C     
            IF(NX.NE.INX .OR. NY.NE.INY) THEN
               WRITE(6,FMT=999) NX,NY,INX,INY
 999           FORMAT(X,'Bad input:',I4,X,I4,' (',I4,X,I4,')')
               IRC=54
               RETURN
            END IF
            FIRST=.FALSE.
         END IF
C
         CA=COSDEG(-SATAZI)
         SA=SINDEG(-SATAZI)
         CZ=SINDEG(SATZEN)/(REH*SIN(TANDEG(SATZEN)/REH))
C     CZ=COSDEG(SATZEN)
         CL=COSDEG(LAT)
         IF (CZ.LT.1.0D-10) THEN
            WRITE(*,*) MYNAME,'Invalid satellite azimuth angle.',
     &           SATAZI
            IRC=374
            RETURN
         END IF
         LATA=LAT+DELTA*CA/(CZ*CZ)
         LONA=LON+DELTA*SA/(CZ*CZ*CL)
         LATB=LAT+DELTA*SA/CZ
         LONB=LON-DELTA*CA/(CZ*CL)
C
C     GET OBS-POSITION IN GRID
C
C     SEMI-MAJOR AXIS (LAND MASK)
         XPOS = LATA
         YPOS = LONA
         CALL GC2XY(1,1,XPOS,YPOS,IGTYPE,GRID,IERROR)
         IF(IERROR.NE.0)THEN
            IRC=56
            RETURN
         END IF
         XAM=XPOS
         YAM=YPOS
C
C     SEMI-MINOR AXIS (LAND MASK)
         XPOS = LATB
         YPOS = LONB
         CALL GC2XY(1,1,XPOS,YPOS,IGTYPE,GRID,IERROR)
         IF(IERROR.NE.0)THEN
            IRC=56
            RETURN
         END IF
         XBM=XPOS
         YBM=YPOS
C
C     SEMI-MAJOR AXIS (ICE MASK)
         XPOS = 2.0D0*LATA-LAT
         YPOS = 2.0D0*LONA-LON
         CALL GC2XY(1,1,XPOS,YPOS,IGTYPE,GRID,IERROR)
         IF(IERROR.NE.0)THEN
            IRC=56
            RETURN
         END IF
         XAMI=XPOS
         YAMI=YPOS
C
C     SEMI-MINOR AXIS (ICE MASK)
         XPOS = 2.0D0*LATB-LAT
         YPOS = 2.0D0*LONB-LON
         CALL GC2XY(1,1,XPOS,YPOS,IGTYPE,GRID,IERROR)
         IF(IERROR.NE.0)THEN
            IRC=56
            RETURN
         END IF
         XBMI=XPOS
         YBMI=YPOS
C
C     CENTER POSITION
         XPOS = LAT
         YPOS = LON
         CALL GC2XY(1,1,XPOS,YPOS,IGTYPE,GRID,IERROR)
         IF(IERROR.NE.0)THEN
            IRC=56
            RETURN
         END IF
         XCM=XPOS
         YCM=YPOS
C
C     CHECK IF POSITION IS OUTSIDE GRID
C
         IF (XPOS.LT.0 .OR. XPOS.GT.NX .OR.
     &        YPOS.LT.0 .OR. YPOS.GT.NY) THEN
            BPRO=.FALSE.
            RETURN
         ELSE
            BPRO=.TRUE.
         END IF
C
C     GET FOOTPRINT
C
         CALL FOOTPR (NX,NY,BPRO,
     &     XCM,YCM,XAM,YAM,XBM,YBM,
     &     MAXPOS,NRPOS,POS,
     &     IRC)
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from FOOTPR.',IRC
            RETURN
         END IF
C
         IF (.NOT.BPRO) RETURN
         II=1                   ! LAND
         FIX_INT(II)=0.0D0
         DO KK=1,NRPOS
            FIX_INT(II)=FIX_INT(II)+
     &           FIXDAT(POS(1,KK),POS(2,KK),II)
         END DO
         FIX_INT(II)=FIX_INT(II)/DFLOAT(NRPOS)
C
C     GET FOOTPRINT
C
         CALL FOOTPR (NX,NY,BPRO,
     &     XCM,YCM,XAMI,YAMI,XBMI,YBMI,
     &     MAXPOS,NRPOS,POS,
     &     IRC)
         IF (IRC.NE.0) THEN
            WRITE(*,*) MYNAME,'Error return from FOOTPR.',IRC
            RETURN
         END IF
C
         IF (.NOT.BPRO) RETURN
         II=3                   ! ICE
         FIX_INT(II)=0.0D0
         DO KK=1,NRPOS
            FIX_INT(II)=FIX_INT(II)+
     &           FIXDAT(POS(1,KK),POS(2,KK),II)
         END DO
         FIX_INT(II)=FIX_INT(II)/DFLOAT(NRPOS)
C     
C     BILINEAER INETRPOLASJON
C
         X=XPOS
         Y=YPOS
         J=Y
         Y1=Y-J
         I=X
         X1=X-I
         X5=1-X1-Y1+X1*Y1
         X6=X1-X1*Y1
         X7=Y1-X1*Y1
         X8=X1*Y1
C
         IF (I+3.GT.NX .OR. I-1.LT.1 .OR.
     &        J+3.GT.NY .OR. J-1.LT.1) THEN
            BPRO=.FALSE.
            RETURN
         END IF
C
C     LOOPS ON PARAMETER AND LEVELS
C
         DO NT = 1,NTID
            DO II = 1,NPAR
               DO JJ = 1,ILEV
                  FELT_INT2(JJ,II,NT) =
     &                 X5*FELT(I,J,JJ,II,NT)+X6*FELT(I+1,J,JJ,II,NT)
     &                 +X7*FELT(I,J+1,JJ,II,NT)+
     &                 X8*FELT(I+1,J+1,JJ,II,NT)
C     IF (FELT_INT2(JJ,II,NT).LT.0.0D0) THEN
C     WRITE(*,*) MYNAME,FELT_INT2(JJ,II,NT)
C     END IF
               ENDDO
            ENDDO
            DO II = 1,NPARS
               SURF_INT2(II,NT) =
     &              X5*SURFACE(I,J,II,NT)+X6*SURFACE(I+1,J,II,NT)
     &              +X7*SURFACE(I,J+1,II,NT)+X8*SURFACE(I+1,J+1,II,NT)
            ENDDO
         ENDDO
C
C     NEED A WIDE LAND MASK...
C
         DO II=1,NPAR2
            IF (II.EQ.1.OR.II.EQ.3) THEN ! wide land/ice mask
C     THIS IS DONE ELSEWHERE
            ELSE
               FIX_INT(II)=
     &              X5*FIXDAT(I,J,II)+X6*FIXDAT(I+1,J,II)
     &              +X7*FIXDAT(I,J+1,II)+X8*FIXDAT(I+1,J+1,II)
            END IF
         END DO
C     
C     FIND TIME INTERPOLATION FACTOR
C         
         DT=(J_OBS - J_FELT - FTIME(1)/24.0D0)/
     &        (FTIME(2) - FTIME(1)/24.0)
         IF (DT.LT.0 .OR. DT.GT.1.0) THEN
            WRITE(*,*) MYNAME,'SYSTEM ERROR, INVALID TIME INT. FACTOR',
     &           DT
            IRC=724
            RETURN
         END IF
C
C     LINEAR INTERPOLATION
C
         DO II=1,NPAR
            DO JJ=1,ILEV
               FELT_INT(JJ,II)=DT*FELT_INT2(JJ,II,2)+
     &              (1.0-DT)*FELT_INT2(JJ,II,1)
C     IF (FELT_INT(JJ,II).LT.0.0D0) THEN
C     WRITE(*,*) MYNAME,FELT_INT(JJ,II)
C     END IF
            END DO
         END DO
         DO II=1,NPARS
            SURF_INT(II)=DT*SURF_INT2(II,2)+(1.0-DT)*SURF_INT2(II,1)
         END DO
C
C     ASSIGN VALUE TO OUTPUT VARIABLES
C
         IF (SURF_INT(1).GT.500) THEN
            GPRES = SURF_INT(1)
            GWIND = SQRT(SURF_INT(5)*SURF_INT(5)+
     &           SURF_INT(6)*SURF_INT(6))
            GWINDIR=ATAN2DEG(SURF_INT(6),SURF_INT(5))
C     
            SOILM  = 0.0D0
            FRACTL = MIN(MAX(0.0D0,FIX_INT(1)*0.01D0),1.0D0)
            HEIGHT = FIX_INT(2)
            ALBEDO = 0.0D0
C            ALBEDO = FIX_INT(3)
            GTEMP  = FRACTL*SURF_INT(2) + (1.0D0-FRACTL)*FIX_INT(4)
            ICEC   = MIN(MAX(0.0D0,FIX_INT(3)*0.01D0),1.0D0)
C     
            CLCOV = 0.0
            CLTOP = 200
            CLBOT = 500
C
            R = SURF_INT(4)/100.0D0 ! relative humidity (0..1)
            P = GPRES*100.0 ! pressure in pascal
            CALL RQ(P,GWATV,R,GTEMP,1)
            GOZONE = 0.0D0 ! should eventually be close to 0.00000003496
            GWATL = 0.0D0
C     
            KAPPA = 0.286D0
            RPRES = 1000.0D0
C     
            DO II=1,ILEV
               CURPRES = A(II)+B(II)*GPRES
               CURTEMP = FELT_INT(II,4)*(CURPRES/RPRES)**(KAPPA)
               PRES(II)  = CURPRES
               TARGPRES(II)  = CURPRES
               TEMP(II)  = CURTEMP
               IF ( FELT_INT(II,3).LT.0.0D0) THEN
C     WRITE(*,*) MYNAME,'Negative Watv:',FELT_INT(II,3)
                  nrn=nrn+1
               END IF
               WATV(II)  = MAX( FELT_INT(II,3),0.0D0)
               WATL(II)  = 0.0
            END DO
C
C     CHECK FOR SUPER-SATURATION
C
            CALL SATURATE_FT(KODE,1,GPRES,GTEMP,GWATV,IRC)
            IF (IRC.NE.0) THEN
               WRITE(*,*) MYNAME,'Error return from SATURATE_FT.',IRC
               RETURN
            END IF
C
            CALL SATURATE_FT(KODE,NRLEV,PRES,TEMP,WATV,IRC)
            IF (IRC.NE.0) THEN
               WRITE(*,*) MYNAME,'Error return from SATURATE_FT.',IRC
               RETURN
            END IF
C
C     GET OZONE CLIMATOLOGY...
C
            TYPE=0              ! pressure levels (mBar)
            CALL GETCLIM(TYPE,TIME,LAT,LON,CPRES,IRC)
C     
            TYPE=7              ! ozone (Kg/Kg)
            CALL GETCLIM(TYPE,TIME,LAT,LON,COZONE,IRC)
C     
            TPRES=GPRES
            CALL GETARR(CLEV,CPRES,COZONE,
     &           TPRES,TOZONE,IRC)
            IF (IRC.NE.0) THEN
               WRITE(*,*) MYNAME,'Error return from GETARR.',IRC
               RETURN
            END IF
            GOZONE=MAX(TOZONE,0.0D0)
C     
C     INTERPOLATE OZONE TO THE HIRLAM PRESSURE LEVELS
C     
            IKODE=-1
            CALL PROINT (IKODE,ILEV, TARGPRES, TARGOZONE,
     &           CLEV, CPRES, COZONE,
     &           CLEV, CPRES, COZONE,
     &           IRC)
            IF (IRC.NE.0) THEN
               WRITE(*,*) MYNAME,'Error return from PROINT',IRC
               RETURN
            END IF
C     
            DO II=1,ILEV
               OZONE(II) = TARGOZONE(II)
            END DO
C
         ELSE
            IRC=997
            RETURN
         END IF
C
      ELSE IF (KODE.EQ.1) THEN
C
         OPENED=.FALSE.
         write(*,*) myname,'Number of negative watv:', nrn
C
      ELSE
         WRITE(*,*) MYNAME,'UNKNOWN KODE:',KODE
         IRC=999
         RETURN
      END IF
C
      IF (BDEB) WRITE(*,*) MYNAME,'Debug: Routine ends',IRC
C
      RETURN
      END

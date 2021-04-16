      SUBROUTINE MAPGRID (NRLON, NRLAT, OGTYPE,OG,GRID, 
     &     NROBS, LONA, LATA, NRIND, IND,
     &     SATZENA,SATAZIA,DELTA,REH,UNDEF,
     &     MAXLONB,MINLONB,MAXLATB,MINLATB,
     &     IRC)
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
C     MAKES A GRID INDEX MAP (TO THE OBSERVATIONS)...
C
      IMPLICIT NONE
      SAVE
C
      INTEGER  NRLON,NRLAT,GRID(NRLON,NRLAT),NROBS,
     &     NRIND, IND(2,NRIND)
      REAL   LONA(NROBS), LATA(NROBS), 
     &     LONC,LATC,SATZENA(NROBS),SATAZIA(NROBS),DELTA,REH
      INTEGER  UNDEF
      REAL   MAXLONB,MINLONB,MAXLATB,MINLATB
      INTEGER  IRC
      real cosdeg,sindeg,tandeg,atan2deg
      external cosdeg,sindeg,tandeg,atan2deg
C
      INTEGER  NROBSB,ONE,igtype,ogtype,ierror
      DATA ONE /1/
      REAL LONB(NROBS),LATB(NROBS),IG(6),OG(6),XPOS,YPOS
      REAL   MAXLONA,MINLONA,MAXLATA,MINLATA
C
      INTEGER  NDELTA,POSA(2,NROBS),
     &     II,JJ,KK,LL,MM,LS,LE,JS,JE,
     &     MAXPOS,NRPOS
      INTEGER, allocatable::  POS(:,:)
      REAL   DIFF,DDLON,DDLAT,DLON,DLAT,DLONA,DLATA,
     &     XCM,YCM,XAM,YAM,XBM,YBM,CA,SA,CZ,CL,DD
      LOGICAL  BPRO
      CHARACTER*8 MYNAME
      DATA MYNAME /'MAPGRID'/
C
C      write(*,*)myname,'Entering.'
C     
C     ASSIGN INTERNAL TABLES
C
      NROBSB=NROBS
      IF (NROBSB.LE.0) THEN
         WRITE(*,*) MYNAME,'Too few observations!',NROBSB
         IRC=293
         RETURN
      END IF
C
C     allocate work-array
C
      maxpos=2*max(1000,(nrlon*nrlat)/nrobs)
      allocate(pos(2,maxpos),stat=irc)
      if (irc.ne.0) then
         write(*,*) myname,'Unable to allocate:',maxpos
         return
      end if
C
      DO KK=1,NROBSB
         LONB(KK)=LONA(KK)
         LATB(KK)=LATA(KK)
         IF (KK.EQ.1) THEN
            MAXLONA=LONA(KK)
            MAXLATA=LATA(KK)
            MINLONA=LONA(KK)
            MINLATA=LATA(KK)
         ELSE
            MAXLONA=MAX(MAXLONA,LONA(KK))
            MAXLATA=MAX(MAXLATA,LATA(KK))
            MINLONA=MIN(MINLONA,LONA(KK))
            MINLATA=MIN(MINLATA,LATA(KK))
         END IF
      END DO
C
      DLONA=MAXLONA-MINLONA
      DLATA=MAXLATA-MINLATA
C
      IF (DLATA.LT.1.0D-10.OR.DLONA.LT.1.0D-10) THEN
         WRITE(*,*) MYNAME,'Obs. area too small.'
         IRC=932
         RETURN
      END IF
C
C     ROTATE COORDINATES
C
      IG(1)=1.0
      IG(2)=1.0
      IG(3)=1.0
      IG(4)=1.0
      IG(5)=0.0
      IG(6)=0.0
      IGTYPE=2                  ! CONVERT FROM GEOGRAPHIC LAT,LON
C
      CALL XYCONVERT(NROBSB,LONB(1),LATB(1),IGTYPE,IG,
     +     OGTYPE,OG,IERROR)
      IF (IERROR.NE.0) THEN
         WRITE(*,*) MYNAME,'Error return from XYCONVERT.',IERROR
         IRC=IERROR
         RETURN
      END IF
C      write(*,*)myname,'POS:',lona(1),lata(1),lonb(1),latb(1)
C
C     DETERMINE MIN/MAX IN TRANSFORMED COORDINATE SYSTEM
C
      DO II=1,NROBS
         IF (II.EQ.1) THEN
            MAXLONB=LONB(II)
            MAXLATB=LATB(II)
            MINLONB=LONB(II)
            MINLATB=LATB(II)
         ELSE
            MAXLONB=MAX(MAXLONB,LONB(II))
            MAXLATB=MAX(MAXLATB,LATB(II))
            MINLONB=MIN(MINLONB,LONB(II))
            MINLATB=MIN(MINLATB,LATB(II))
         END IF
      END DO
C
      DLON=MAXLONB-MINLONB
      DLAT=MAXLATB-MINLATB
C
      MAXLONB=MAXLONB+DELTA*(DLAT/DLATA)
      MAXLATB=MAXLATB+DELTA*(DLAT/DLATA)
      MINLONB=MINLONB-DELTA*(DLAT/DLATA)
      MINLATB=MINLATB-DELTA*(DLAT/DLATA)
C
      IF (OGTYPE.EQ.9) THEN ! don't readjust coordinates
         MINLATB=MINLATB-DLAT/1000.0
         MAXLATB=MAXLATB+DLAT/1000.0
         MINLONB=MINLONB+DLON/1000.0
         MAXLONB=MAXLONB-DLAT/1000.0
      ELSEIF (OGTYPE.NE.1.AND.OGTYPE.NE.4) THEN ! don't readjust polar stereogr.
         DIFF=DLON-DLAT
         IF (DIFF.GT.0.0) THEN
            MINLATB=MAX(-90.0,MINLATB-DIFF/2.0)
            MAXLATB=MIN(90.0,MAXLATB+DIFF/2.0)
         ELSE
            MINLONB=MAX(MINLONB+DIFF/2.0,-180.0)
            MAXLONB=MIN(MAXLONB-DIFF/2.0,180.0)
         END IF
      ELSE
         DIFF=DLON-DLAT
         IF (DIFF.GT.0) THEN
            MINLATB=MINLATB-DIFF/2.0
            MAXLATB=MAXLATB+DIFF/2.0
         ELSE
            MINLONB=MINLONB+DIFF/2.0
            MAXLONB=MAXLONB-DIFF/2.0
         END IF
      END IF
      NDELTA=NRLAT*(DELTA/DLATA)
C
      DLON=MAXLONB-MINLONB
      DLAT=MAXLATB-MINLATB
C     
      DDLON=DLON/NRLON
      DDLAT=DLAT/NRLAT
C
C     SET GRID TO UNDEFINED
C
      GRID(:,:)=UNDEF
C
C     LOOP THROUGH POSITIONS, FILL ADJACENT PIXELS
C
      DO MM=1,NRIND
C
C     MARK CENTER-POINTS
C
         DO II=IND(1,MM),IND(2,MM)
            if (ii.lt.1.or.ii.gt.nrobs) then
               write(*,*)myname,'Error:',ii,nrobs
            end if
            POSA(1,II)=1+(NRLON-1)*(LONB(II)-MINLONB)/DLON
            POSA(2,II)=1+(NRLAT-1)*(LATB(II)-MINLATB)/DLAT
            IF (POSA(1,II).LT.1 .OR. POSA(1,II).GT.NRLON) THEN
               WRITE(*,*) MYNAME,'System error.',
     &              POSA(1,II),MINLONB,LONB(II),DLON,NRLON
               posa(1,ii)=max(min(POSA(1,II),nrlon),1)
            END IF
            IF (POSA(2,II).LT.1 .OR. POSA(2,II).GT.NRLAT) THEN
               WRITE(*,*) MYNAME,'Another system error.',
     &              POSA(2,II),MINLATB,LATB(II),DLAT,NRLAT
               posa(2,ii)=max(min(POSA(2,II),nrlon),1)
            END IF
            IF (GRID(POSA(1,II),POSA(2,II)).EQ.0) 
     &           GRID(POSA(1,II),POSA(2,II))=II
         END DO
C
C     MAKE FOOTPRINT
C
         DO II=IND(1,MM),IND(2,MM)
            IF (GRID(POSA(1,II),POSA(2,II)).EQ.II) THEN
C
C     GET OBS-POSITION IN GRID
C
C
               CA=COSDEG(-SATAZIA(II))
               SA=SINDEG(-SATAZIA(II))
               CZ=SINDEG(SATZENA(II))/(REH*SIN(TANDEG(SATZENA(II))/REH))
C     CZ=COSD(SATZENA(II))
               CL=COSDEG(LATA(II))
               IF (CZ.LT.1.0D-10) THEN
                  WRITE(*,*) MYNAME,
     &                 'Invalid satellite azimuth angle.',
     &                 SATAZIA(II)
                  IRC=374
                  RETURN
               END IF
C     SEMI-MAJOR AXIS
               XPOS = LONA(II)+DELTA*SA/(CZ*CZ*CL)
               YPOS = LATA(II)+DELTA*CA/(CZ*CZ)
               CALL XYCONVERT(ONE,XPOS,YPOS,IGTYPE,IG,
     +              OGTYPE,OG,IERROR)
               IF (IERROR.NE.0) THEN
                  WRITE(*,*) MYNAME,'Error return from XYCONVERT.',
     &                       IERROR
                  IRC=IERROR
                  RETURN
               END IF
               XAM=1+(NRLON-1)*(XPOS-MINLONB)/DLON
               YAM=1+(NRLAT-1)*(YPOS-MINLATB)/DLAT
C     
C     SEMI-MINOR AXIS
               XPOS = LONA(II)-DELTA*CA/(CZ*CL)
               YPOS = LATA(II)+DELTA*SA/CZ
               CALL XYCONVERT(ONE,XPOS,YPOS,IGTYPE,IG,
     +              OGTYPE,OG,IERROR)
               IF (IERROR.NE.0) THEN
                  WRITE(*,*) MYNAME,'Error return from XYCONVERT.',
     &                       IERROR
                  IRC=IERROR
                  RETURN
               END IF
               XBM=1+(NRLON-1)*(XPOS-MINLONB)/DLON
               YBM=1+(NRLAT-1)*(YPOS-MINLATB)/DLAT
C
               XCM=1+(NRLON-1)*(LONB(II)-MINLONB)/DLON
               YCM=1+(NRLAT-1)*(LATB(II)-MINLATB)/DLAT
C
C     GET FOOTPRINT
C
               BPRO=.TRUE.
               CALL FOOTPR (NRLON,NRLAT,BPRO,
     &              XCM,YCM,XAM,YAM,XBM,YBM,
     &              MAXPOS,NRPOS,POS,
     &              IRC)
               IF (IRC.NE.0) THEN
                  WRITE(*,*) MYNAME,'Error return from FOOTPR.',IRC
                  RETURN
               END IF
               IF (BPRO) THEN
                  DO LL=1,NRPOS
                     IF (GRID(POS(1,LL),POS(2,LL)).EQ.0) THEN
                        GRID(POS(1,LL),POS(2,LL))=II
                     ELSE IF (GRID(POS(1,LL),POS(2,LL))
     &                       .GE.IND(1,MM) .AND.
     &                       GRID(POS(1,LL),POS(2,LL))
     &                       .LE.IND(2,MM)) THEN
C     USE OBS CLOSEST TO PIXEL
                        DD=(POS(1,LL)-POSA(1,II))*
     &                       (POS(1,LL)-POSA(1,II))*DDLON*DDLON+
     &                       (POS(2,LL)-POSA(2,II))*
     &                       (POS(2,LL)-POSA(2,II))*DDLAT*DDLAT -
     &                       (POS(1,LL)-POSA(1,GRID(POS(1,LL),
     &                       POS(2,LL))))*
     &                       (POS(1,LL)-POSA(1,GRID(POS(1,LL),
     &                       POS(2,LL))))*DDLON*DDLON -
     &                       (POS(2,LL)-POSA(2,GRID(POS(1,LL),
     &                       POS(2,LL))))*
     &                       (POS(2,LL)-POSA(2,GRID(POS(1,LL),
     &                       POS(2,LL))))*DDLAT*DDLAT
                        IF (DD.LT.0) THEN
                           GRID(POS(1,LL),POS(2,LL))=II
                        END IF
                     END IF
                  END DO
               ELSE
C     OUT OF RANGE
               END IF
            END IF
         END DO
      END DO
C
      if (allocated(pos)) deallocate(pos)
      RETURN
      END


      SUBROUTINE INDXARR (MAXIND,
     &     MAXLEV,TARGLEV, TARGPRES,   ! target levels
     &     CLEV,CPRES,          ! climatology
     &     NRLEV,PRES,          ! profile
     &     NRIND,QCLIM,ALPHA,INDX,
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
      IMPLICIT NONE
      SAVE
C
      INTEGER  MAXLEV,TARGLEV
      REAL   TARGPRES(TARGLEV)
      INTEGER  CLEV
      REAL   CPRES(CLEV)
      INTEGER  NRLEV
      REAL   PRES(NRLEV)
      INTEGER  IRC
C
      CHARACTER*8 MYNAME
      DATA MYNAME /'INDXARR'/
C
      INTEGER  MAXIND
      INTEGER  NRIND(MAXLEV)
      LOGICAL  QCLIM(MAXLEV,MAXIND)
      REAL   ALPHA(MAXLEV,MAXIND)
      INTEGER  INDX(MAXLEV,MAXIND)
C
      INTEGER  CPOS,II
      REAL   DPT,DP,FACT,FX
      DATA FX/10.0D0/
      LOGICAL  BDONE, FERR
      DATA FERR/.TRUE./
C
      CPOS=1
C
      DO II=1,TARGLEV
C     
         IF (TARGPRES(II).GT.PRES(NRLEV).OR.
     &        TARGPRES(II).LT.PRES(1)) THEN ! USE CLIMATOLOGY
C     
C     FIND POS IN CLIMATOLOGY ARRAY
C         
            IF (TARGPRES(II).LT.CPRES(CLEV).AND.
     &        TARGPRES(II).GT.CPRES(1)) THEN ! USE CLIMATOLOGY
               
               CPOS=MAX(MIN(CPOS,CLEV),1)
               BDONE=(TARGPRES(II).LT.CPRES(CPOS))
               DO WHILE(.NOT.BDONE)
                  CPOS=CPOS+1
                  IF (CPOS.GT.CLEV) THEN
                     WRITE(*,*) MYNAME,'System error.'
                     IRC=989
                     RETURN
                  ELSE
                     BDONE = (TARGPRES(II).LT.CPRES(CPOS))
                  END IF
               END DO
C     
               CPOS=MAX(MIN(CPOS,CLEV),1)
               BDONE=(TARGPRES(II).GT.CPRES(CPOS))
               DO WHILE(.NOT.BDONE)
                  CPOS=CPOS-1
                  IF (CPOS.LT.1) THEN
                     WRITE(*,*) MYNAME,'System error.'
                     IRC=990
                     RETURN
                  ELSE
                     BDONE = (TARGPRES(II).GT.CPRES(CPOS))
                  END IF
               END DO
C     
               DPT=CPRES(CPOS+1)-CPRES(CPOS)
               DP=TARGPRES(II)-CPRES(CPOS)
C
               IF (TARGPRES(II).GT.PRES(NRLEV)) THEN
                  FACT=EXP(-FX*
     &                 (TARGPRES(II)-PRES(NRLEV))/
     &                 (CPRES(CPOS+1)))
               ELSE
                  FACT=EXP(FX*
     &                 (TARGPRES(II)-PRES(1))/
     &                 (CPRES(CPOS)))
               END IF
C     
               NRIND(II)=3
               IF(NRIND(II).GT.MAXIND) THEN
                  WRITE(*,*) MYNAME,'Maxind error.',MAXIND
                  IRC=993
                  RETURN
               END IF
C     
               QCLIM(II,1)=.TRUE.
               ALPHA(II,1)=(1.0D0-FACT)*DP/DPT ! weight
               INDX(II,1)=CPOS+1 ! array index 
C     
               QCLIM(II,2)=.TRUE.
               ALPHA(II,2)=(1.0D0-FACT)*(1.0D0-(DP/DPT))
               INDX(II,2)=CPOS
               
                  
               IF (TARGPRES(II).GT.PRES(NRLEV)) THEN
                  QCLIM(II,3)=.FALSE.
                  ALPHA(II,3)=FACT
                  INDX(II,3)=NRLEV
               ELSE
                  QCLIM(II,3)=.FALSE.
                  ALPHA(II,3)=FACT
                  INDX(II,3)=1
               END IF
C
            ELSE
C
C     OUTSIDE CLIMATOLOGY TOO!!
C
               WRITE(*,*) MYNAME,'Pressure out of bounds.',
     &              TARGPRES(II)
               IRC=234
               RETURN
C
            END IF
C
         ELSE
C
C     HANDLE THE BORDERLINE
C
            IF (ABS(TARGPRES(II)-PRES(NRLEV)).LT.1.0D-5) THEN ! max pressure
               CPOS=NRLEV
            ELSE IF (ABS(TARGPRES(II)-PRES(1)).LT.1.0D-5) THEN ! min pressure
               CPOS=2
            ELSE
C     
C     FIND POS IN PRES ARRAY
C     
               CPOS=MAX(MIN(CPOS,NRLEV),1)
               BDONE=(TARGPRES(II).LE.PRES(CPOS))
               DO WHILE(.NOT.BDONE)
                  CPOS=CPOS+1
                  IF (CPOS.GT.NRLEV) THEN
                     WRITE(*,*) MYNAME,'System error.'
                     IRC=991
                     RETURN
                  ELSE
                     BDONE = (TARGPRES(II).LE.PRES(CPOS))
                  END IF
               END DO
C     
               CPOS=MAX(MIN(CPOS,NRLEV),1)
               BDONE=(TARGPRES(II).GT.PRES(CPOS)+1.0D-5)
               DO WHILE(.NOT.BDONE)
                  CPOS=CPOS-1
                  IF (CPOS.LT.1) THEN
                     WRITE(*,*) MYNAME,'System error.'
                     IRC=992
                     RETURN
                  ELSE
                     BDONE = (TARGPRES(II).GT.PRES(CPOS))
                  END IF
               END DO
            END IF
C     
            DPT=PRES(CPOS+1)-PRES(CPOS)
            DP=TARGPRES(II)-PRES(CPOS)
C     
            NRIND(II)=2
            IF(NRIND(II).GT.MAXIND) THEN
               WRITE(*,*) MYNAME,'Maxind error.',MAXIND
               IRC=993
               RETURN
            END IF
C     
            QCLIM(II,1)=.FALSE.
            ALPHA(II,1)=DP/DPT  ! weight
            INDX(II,1)=CPOS+1   ! array index 
C     
            QCLIM(II,2)=.FALSE.
            ALPHA(II,2)=1.0D0-DP/DPT
            INDX(II,2)=CPOS
C     
         END IF
      END DO
C     
C     WRITE(*,*) MYNAME,'Routine not implemented.'
C     IRC=999
C
      RETURN
      END

      SUBROUTINE SATURATE_FT (KODE,NRLEV, PRES,TEMP,WATV,IRC)
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
      INTEGER  KODE,NRLEV
      REAL   PRES(NRLEV), TEMP(NRLEV), WATV(NRLEV)
      INTEGER  IRC
C
      REAL   SVP,AA1, AA2,SATLIM,MEANCORR,STDCORR
      EXTERNAL SVP
      INTEGER  II,NRCORR
      DATA SATLIM /1.01D0/
      LOGICAL  INITIALISED
C
      CHARACTER*12 MYNAME
      DATA MYNAME/'SATURATE_FT'/
      DATA INITIALISED /.FALSE./
C
      IF (KODE.EQ.-1) THEN
         IF (INITIALISED) THEN
            WRITE(*,*) MYNAME,'ROUTINE IS ALLREADY INITIALISED'
            IRC=36
            RETURN
         END IF
C
         NRCORR=0
         MEANCORR=0.0D0
         STDCORR=0.0D0
         INITIALISED=.TRUE.
C
      ELSE IF (KODE.EQ.0) THEN
         IF (.NOT.INITIALISED) THEN
            WRITE(*,*) MYNAME,'INITIALISE BEFORE USING ROUTINE'
            IRC=41
            RETURN
         END IF
C         
         DO II=1,NRLEV
            AA1=SVP(TEMP(II))*SATLIM ! CURRENT TEMP
            IF (AA1.LT.PRES(II)) THEN
               AA2=0.622*AA1/(PRES(II)-AA1) ! TARPRES(II) = CURRENT PRESSURE
               IF (WATV(II).GT.AA2) THEN
                  NRCORR=NRCORR+1
                  MEANCORR=MEANCORR+(WATV(II)-AA2)
                  STDCORR=STDCORR+(WATV(II)-AA2)**2
                  WATV(II)=AA2
               ENDIF
            ENDIF
         END DO
C
      ELSE IF (KODE.EQ.1) THEN
C
         IF (INITIALISED) THEN
            IF (NRCORR.NE.0) THEN
               MEANCORR=MEANCORR/DFLOAT(NRCORR)
               STDCORR=DSQRT(MAX(0.0D0,
     &              STDCORR/DFLOAT(NRCORR)-
     &              MEANCORR*MEANCORR))
               WRITE(*,999) MYNAME,NRCORR,MEANCORR,STDCORR
 999           format(X,A,'Number of saturated Watv:',I8,
     &              ' (Mn:',E10.2,' Std:',E10.2,')')
            ELSE
               WRITE(*,*) MYNAME,'No Watv were saturated.'
            END IF
            INITIALISED=.FALSE.
         END IF
      END IF
C
      RETURN
      END

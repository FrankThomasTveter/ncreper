      SUBROUTINE GET_IWATV ( IWATV, NRLEV, PRES, WATV, GPRES,
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
C     INTEGRATE TOTAL WATER VAPOUR CONTENT (KG/M**2)
C
      IMPLICIT NONE
      SAVE
C
      REAL   IWATV
      INTEGER  NRLEV
      REAL   PRES(NRLEV), WATV(NRLEV), GPRES
      INTEGER  IRC
C
      REAL   PTOP, PBOT,G
      PARAMETER (G=9.81D0)
      INTEGER  II
C      
      IWATV=0.0D0
      PTOP=0.0D0
      II=1
      DO WHILE (II.LT.NRLEV .AND. PRES(II).LT.GPRES)
         PBOT=0.5D0*(PRES(II)+PRES(II+1))
         IWATV=IWATV+g*WATV(II)*(PBOT-PTOP)
         PTOP=PBOT
         II=II+1
      END DO
      PBOT=GPRES
      IWATV=IWATV+g*WATV(II)*(PBOT-PTOP)
C
      RETURN
      END


      SUBROUTINE GINTUV(IUNDEF,UNDEF,NX,NY,U,V,GRIDIN,GRIDOT)
C
C        Rotation of wind components after interpolation.
C        Undefined values allowed.
C        Use GINTER for each component and then GINTUV.
C
C  Input:
C        IUNDEF: 0 = no undefined values in the fields
C                1 = undefined values in the fields
C        UNDEF:      undefined value (input and output)
C        GRIDIN(4):  XP,YP,AN,FI for input  grid, i.e. input to GINTER
C        GRIDOT(4):  XP,YP,AN,FI for output grid
C        NX,NY:      dimensions of fields
C        U(NX,NY):   unrotated x component field
C        V(NX,NY):   unrotated y component field
C
C  Output:
C        U(NX,NY):   x component field
C        V(NX,NY):   y component field
C
C-------------------------------------------------------------------
C  DNMI/FoU   18.11.1992   Anstein Foss
C-------------------------------------------------------------------
C
      INTEGER IUNDEF,NX,NY
      REAL    UNDEF
      REAL    U(NX,NY),V(NX,NY),GRIDIN(4),GRIDOT(4)
C
      FI =GRIDOT(4)
      FIB=GRIDIN(4)
C
C..check if rotation is required
      IF(ABS(FI-FIB).LT.0.001) RETURN
C
      ALFA=(FI-FIB)*3.1415927/180.
C
      DUU=+COS(ALFA)
      DVU=+SIN(ALFA)
      DUV=-SIN(ALFA)
      DVV=+COS(ALFA)
C
      IF(IUNDEF.EQ.0) THEN
        DO J=1,NY
          DO I=1,NX
            UI=U(I,J)
            VI=V(I,J)
            U(I,J)=DUU*UI+DVU*VI
            V(I,J)=DUV*UI+DVV*VI
          END DO
        END DO
      ELSE
        UD =UNDEF
        UDL=UNDEF*0.9
        DO J=1,NY
          DO I=1,NX
            IF(U(I,J).LT.UDL .AND. V(I,J).LT.UDL) THEN
              UI=U(I,J)
              VI=V(I,J)
              U(I,J)=DUU*UI+DVU*VI
              V(I,J)=DUV*UI+DVV*VI
            ELSE
              U(I,J)=UD
              V(I,J)=UD
            END IF
          END DO
        END DO
      END IF
C
      RETURN
      END

      SUBROUTINE XYTOXY(GRID1,GRID2,CXY)
C
C        Calculate constants for polarstereographic grid transformation.
C
C        For calculations from 'GRID1' to 'GRID2'
C
C        (X1,Y1) - position in 'GRID1'
C        (X2,Y2) - position in 'GRID2'
C
C        X2 = CXY(1) + CXY(2)*X1 + CXY(3)*Y1
C        Y2 = CXY(4) + CXY(5)*X1 + CXY(6)*Y1
C
C        INPUT:  GRID1(1:4) - XP,YP,AN,FI
C                GRID2(1:4) - XP,YP,AN,FI
C        OUTPUT: CXY(1:6)
C
C-------------------------------------------------------------------
C  DNMI/FoU   xx.xx.1990   Anstein Foss
C-------------------------------------------------------------------
C
      REAL      GRID1(4),GRID2(4),CXY(6)
C
      REAL      ALFA,CALFA,SALFA,S
C
      ALFA=(GRID1(4)-GRID2(4))*3.1415927/180.
      CALFA=COS(ALFA)
      SALFA=SIN(ALFA)
      S=GRID2(3)/GRID1(3)
      CXY(2)=+S*CALFA
      CXY(3)=-S*SALFA
      CXY(5)=+S*SALFA
      CXY(6)=+S*CALFA
      CXY(1)=GRID2(1)-CXY(2)*GRID1(1)-CXY(3)*GRID1(2)
      CXY(4)=GRID2(2)-CXY(5)*GRID1(1)-CXY(6)*GRID1(2)
C
      RETURN
      END

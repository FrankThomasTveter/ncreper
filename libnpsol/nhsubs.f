*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  nhsubs.f:  Auxilliary routines developed at Norsk Hydro
*
*     npnout
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE NPNOUT( NOUT )

      INTEGER NOUT

C***********************************************************************
C     NPNOUT   sets the logical file unit for NPSOL output
C***********************************************************************

      DOUBLEPRECISION WMACH
      COMMON /SOLMCH/ WMACH(15)
      SAVE   /SOLMCH/
C
      WMACH(11) = NOUT
C
      RETURN
      END

      subroutine rlunit(lrunit)
c
c  PURPOSE: Return the unit length in bytes used for RECL= when
c           opening a direct access file.
c           According to the ANSI FORTRAN 77 standard this unit
c           is machine dependant.
c
c-----------------------------------------------------------------------
c  DNMI/FoU  09.10.1992  Anstein Foss
c-----------------------------------------------------------------------
c
c..output:
      integer lrunit
c
c.SGI and DEC: record length in unit 32 bit words (= 4 bytes)
      lrunit=4
c
c.SUN: record length in unit bytes
cc    lrunit=1
c
      return
      end

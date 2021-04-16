      subroutine rlunit(lrunit)
c
c  PURPOSE: Return the unit length in bytes used for RECL= when
c           opening a direct access file.
c           According to the ANSI FORTRAN 77 standard this unit
c           is machine dependant.
c
c-----------------------------------------------------------------------
c  DNMI/FoU  09.10.1992  Anstein Foss
c  DNMI/FoU  29.11.1995  Anstein Foss ... IBM/RS 6000 xl fortran
c-----------------------------------------------------------------------
c
c..output:
      integer lrunit
c
c.SGI and DEC: record length in unit 32 bit words (= 4 bytes)
cc    lrunit=4
c
c.SUN and IBM RS/6000: record length in unit bytes
      lrunit=1
c
      return
      end

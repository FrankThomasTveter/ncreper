      subroutine termchar(chr)
c
c  PURPOSE: Return a character used to terminate text strings
c           for internal free format read.
c
c-----------------------------------------------------------------------
c  DNMI/FoU  05.12.1995  Anstein Foss
c-----------------------------------------------------------------------
c
c..output:
      character*1 chr
c
c..SGI, DEC, SUN:
cc    chr=char(0)
c
c..IBM RS/6000 xl fortran:
      chr=' '
c
      return
      end

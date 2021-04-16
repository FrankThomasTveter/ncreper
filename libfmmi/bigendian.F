      logical function bigendian()
c
c  Returns true if machine is big endian, false if little endian.
c  (SGI:true  Linux/Intel:false)
c
c-----------------------------------------------------------------------
c  met.no/FoU  22.03.2003  Anstein Foss ... first version
c-----------------------------------------------------------------------
c
      implicit none
c
      integer ia,ib
      integer*2   i2
      character*2 c2
      equivalence (i2,c2)
c
      ia=ichar('a')
      ib=ichar('b')
      i2=ia*256+ib
c
      bigendian= (c2.eq.'ab')
c
      return
      end

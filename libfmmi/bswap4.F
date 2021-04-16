      subroutine bswap4(ldata,idata)
c
c       swap bytes in integer*4 words
c
c-----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.199x  Anstein Foss
c-----------------------------------------------------------------------
c
      implicit none
c
      integer   ldata
      integer*4 idata(ldata)
c
      integer i,iswap,iswapo
c
      do i=1,ldata
        iswap=idata(i)
        iswapo=ior(iand(ishft(iswap,-24),255),
     +             iand(ishft(iswap,-8),65280))
        iswapo=ior(iswapo,ishft(iand(iswap,65280),+8))
        iswapo=ior(iswapo,ishft(iand(iswap,255),+24))
        idata(i)=iswapo
      end do
c
      return
      end

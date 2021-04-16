      subroutine bswap2(ldata,idata)
c
c       swap bytes in integer*2 words
c
c-----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.198x  Anstein Foss
c-----------------------------------------------------------------------
c
      implicit none
c
      integer   ldata
      integer*2 idata(ldata)
c
      integer i,iswap
c
      do i=1,ldata
        iswap=idata(i)
        iswap=ior(iand(ishft(iswap,-8),255),
     +            ishft(iand(iswap,255),+8))
        idata(i)=iswap
      end do
c
      return
      end

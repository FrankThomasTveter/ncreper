      subroutine ebcdic(ntext,text,nchar)
c
c  NAME:
c     ebcdic
c
c  PURPOSE:
c     Conversion of characters in a text string from standard ascii
c     to ebcdic for ibm mainframe.
c
c  SYNOPSIS:
c     subroutine ebcdic(ntext,text,nchar)
c     integer       ntext,nchar
c     character*(*) text(ntext)
c
c  INPUT:
c     ntext  - no. of text strings
c     nchar  - no. of characters in each text string
c	       (-1 = all characters in each string)
c
c  INPUT/OUTPUT:
c     text   - the text strings
c
c-----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.198x  Anstein Foss
c  DNMI/FoU  06.04.1995  Anstein Foss ... new library version
c-----------------------------------------------------------------------
c
      implicit none
c
      integer       ntext,nchar
      character*(*) text(ntext)
c
      integer nch,n,k,i
      integer a2e(0:255)
c
      data (a2e(i),i=000,009)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=010,019)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=020,029)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=030,039)/064,064,064,079,127,123,091,108,080,125/
c                                          |   "   #   $   %   &   '
      data (a2e(i),i=040,049)/077,093,092,078,107,096,075,097,240,241/
c                              (   )   *   +   ,   -   .   /   0   1
      data (a2e(i),i=050,059)/242,243,244,245,246,247,248,249,122,094/
c                              2   3   4   5   6   7   8   9   :   ;
      data (a2e(i),i=060,069)/076,126,110,111,124,193,194,195,196,197/
c                              <   =   >   ?   @   a   b   c   d   e
      data (a2e(i),i=070,079)/198,199,200,201,209,210,211,212,213,214/
c                              f   g   h   i   j   k   l   m   n   o
      data (a2e(i),i=080,089)/215,216,217,226,227,228,229,230,231,232/
c                              p   q   r   s   t   u   v   w   x   y
      data (a2e(i),i=090,099)/233,123,124,091,095,109,121,129,130,131/
c                              z   #   @   $   ^   _   `   a   b   c
      data (a2e(i),i=100,109)/132,133,134,135,136,137,145,146,147,148/
c                              d   e   f   g   h   i   j   k   l   m
      data (a2e(i),i=110,119)/149,150,151,152,153,162,163,164,165,166/
c                              n   o   p   q   r   s   t   u   v   w
      data (a2e(i),i=120,129)/167,168,169,192,106,208,161,064,064,064/
c                              x   y   z   {   º   }   ~
      data (a2e(i),i=130,139)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=140,149)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=150,159)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=160,169)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=170,179)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=180,189)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=190,199)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=200,209)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=210,219)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=220,229)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=230,239)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=240,249)/064,064,064,064,064,064,064,064,064,064/
c
      data (a2e(i),i=250,255)/064,064,064,064,064,064/
c
      nch=nchar
      if(nch.lt.0) nch=len(text(1))
c
      do n=1,ntext
	do k=1,nch
	  text(n)(k:k)=char(a2e(ichar(text(n)(k:k))))
	end do
      end do
c
      return
      end

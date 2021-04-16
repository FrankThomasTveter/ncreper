      subroutine utmgeo ( rnorth, reast, izone, glat, glon )
c
c  NAME:
c     utmgeo
c
c  SYNOPSIS:
c     subroutine utmgeo ( rnorth, reast, izone, glat, glon )
c     integer          izone
c     double precision rnorth, reast, glat, glon
c
c-----------------------------------------------------------------------
c     transformasjon fra plane utm til geografisk - koord
c-----------------------------------------------------------------------
c     algoritmen baserer seg p]    1) holsen : geodesi 2 ,kap. 4
c                                  2) ngo    : publikasjon nr 16
c     det er tatt med  b-koeff. til 8. ledd i rekken
c     programmerer :   aage hovden
c     spr]k        :   vax-11 fortran
c     sentralmeridian, glz=9 gr. for sone 32
c
c     kl|dd p} av arne edmund ofstad, lance, 3. oktober, 1985
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c  DNMI/FoU  18.05.1995  Anstein Foss ... library version
c  DNMI/FoU  29.11.1995  Anstein Foss ... IBM/RS 6000 xl fortran
c-----------------------------------------------------------------------
c

      implicit none
      double precision  rnorth, reast, glat, glon,
     <  pi, a, b, e2, em2, f,
     <  zs, glz, ba1, ba2, ba3, ba4, ba5, ba6, e4, e6, e8, e10, fak, g,
     <  fi1, fi2, fi4, fi6, fi8, fi10, diff, dfi, dy, w, rn, cfi, sfi,
     <  t, t2, t4, t6, t8, t10, eta2, eta4, eta6, eta8, eta10, eta12,
     <  d1, d2, d3, d4, d5, b1, b2, b3, b4, b5, b6, b7, b8, dq, dbr,
     <  dlra,rho

      integer izone

      parameter (pi  = 3.1415926535898d0)
      parameter (a   = 6378137.d0)
      parameter (b   = 6356752.3142d0)
      parameter (e2  = 1.0d0 - (b/a)*(b/a))
      parameter (em2 = (a/b)*(a/b) - 1.0d0)
      parameter (zs  = 0.9996d0)
      parameter (rho = 57.295779513082d0)

!      parameter  ( pi  = 3.14159265358979323846264338327950d0 )
!      parameter  ( a   = 6378388.d0         )
!      parameter  ( f   = 1.d0 / 297.d0      )
!      parameter  ( b   = ( 1.d0 - f ) * a   )
!      parameter  ( e2  = 1.d0 - (b/a)*(b/a) )
!      parameter  ( em2 = (a/b)*(a/b) - 1.d0 )
!      parameter  ( zs  = 0.9996d0           )
!      parameter  ( rho = 57.2957795130823208767981548141052d0)

c-------beregning av meridianbuelengde-koeffisienter-----------------

      glz = (6.d0*dble(izone)-183.d0)/rho

      ba1 = (1.d0+e2*(0.75d0+e2*(45.d0/64.d0+e2*(175.d0/256.d0+e2
     1      *(11025.d0/16384.d0+e2*43659.d0/65536.d0)))))

      ba2 = 0.5d0*e2*(.75d0+e2*(15.d0/16.d0+e2*(525.d0/512.d0+e2*
     1      (2205.d0/2048.d0+e2*72765.d0/65536.d0))))

      e4  = e2 * e2
      ba3 = 0.25d0*e4*(15.d0/64.d0+e2*(105.d0/256.d0+e2*(2205.d0/
     1      4096.d0+e2*10395.d0/16384.d0)))

      e6  = e4 * e2
      ba4 = 1.d0/6.d0*e6*(35.d0/512.d0+e2*(315.d0/2048.d0
     1      +e2*31185.d0/131072.d0))

      e8  = e2 * e6
      ba5 = 0.125d0*e8*(315.d0/16384.d0+e2*3465.d0/65536.d0)

      e10 = e8 * e2
      ba6 = 0.1d0*e10*693.d0/131072.d0

ccc----beregning av fi1-----------------------------------------

      fak = 1.d0 / ( a*zs*(1.d0-e2)*ba1)
      fi1 = rnorth*fak

c-------iterasjonsloop start-------------------------------------

      diff = 1.d0
      dfi  = 0.d0

      do while ( dabs(diff) .ge. 0.001d0 )

        fi1=fi1-dfi

        fi2 =fi1+fi1
        fi4 =fi2+fi2
        fi6 =fi4+fi2
        fi8 =fi6+fi2
        fi10=fi8+fi2

        g = a*zs*(1.d0-e2)*(ba1*fi1-ba2*dsin(fi2)+ba3*dsin(fi4)
     1         -ba4*dsin(fi6)+ba5*dsin(fi8)-ba6*dsin(fi10))

        diff=g-rnorth
        dfi=diff*fak

      end do

ccc----iterasjonsloop slutt-------------------------------------

      dy=reast-500000.d0

c-------bredde-avhengige st\rrelser------------------------------

      w   =dsqrt( 1.d0 - e2*dsin(fi1)*dsin(fi1) )
      rn  =zs*a/w
      eta2=em2*dcos(fi1)*dcos(fi1)
      t   =dtan(fi1)
      cfi =dcos(fi1)
      sfi =dsin(fi1)

c------------------------------------------------------------------------

      t2 =t*t
      t4 =t2*t2
      t6 =t4*t2
      t8 =t6*t2
      t10=t8*t2

      eta4 =eta2*eta2
      eta6 =eta2*eta4
      eta8 =eta2*eta6
      eta10=eta4*eta6
      eta12=eta6*eta6

c------- d-koeffisienter skissert i ngo publikasjon nr 16 -----------------

      d1= (1.d0+eta2)*cfi
      d2= - 0.5d0*cfi**2*t*(1.d0+4.d0*eta2+3.d0*eta4)
      d3= - 1.d0/6.d0*cfi**3*(1.d0-t2+5.d0*eta2-1.3d1*eta2*t2+
     1     7.d0*eta4-2.7d1*eta4*t2)
      d4= 1.d0/24.d0*cfi**4*t*(5.d0-t2+5.6d1*eta2-4.d1*eta2*t2)
      d5= 1.d0/120.d0*cfi**5*(5.d0-1.8d1*t2-t4)

c------- b-koeffisienter-----------------------------------------------

      b1=1.d0 / (rn*cfi)

      b2=   t / (2.d0*rn**2*cfi)

      b3= (1.d0+2.d0*t2+eta2) / (cfi*6.d0)

      b4= 1.d0/24.d0*t*(5.d0+6.d0*t2+eta2-4.d0*eta4) / (cfi)

      b5= 1.d0/120.d0*(5.d0+2.8d1*t2+2.4d1*t4+eta2*(6.d0+8.d0*t2)+
     1    eta4*(1-3.d0+4.d0*t2)+eta6*(-4.d0+24.d0*t2)) / cfi

      b6=1.d0/720.d0*t*(61.d0+180.d0*t2+120.d0*t4+46.d0*eta2+
     1   48.d0*t2*eta2+eta4*(-3.d0-36.d0*t2)+eta6*(100.d0-96.d0*t2)
     2   +eta8*(88.d0-192.d0*t2))/ cfi

      b7=1.d0/5040.d0*(61.d0+662.d0*t2+1320.d0*t4+720.d0*t6+eta2
     1   *(107.d0+440.d0*t2+336.d0*t4)+eta4*(43.d0-234.d0*t2
     2   -192.d0*t4)+eta6*(97.d0-772.d0*t2+408.d0*t4)+
     3   eta8*(188.d0-2392.d0*t2+1536.d0*t4)+(88.d0-1632.d0
     4   *t2+1920.d0*t4))/cfi

      b8=1.d0/40320.d0*t*(1385.d0+7266.d0*t2+10920.d0*t4+5040.d0*t6+
     1   eta2*(1731.d0+4416.d0*t2+2688.d0*t4)+eta4*(-573.d0-1830.d0*t2
     2   -1536.d0*t4)+eta6*(-2927.d0+5052.d0*t2+744.d0*t4)+eta8*(
     3   -8808.d0+27456.d0*t2-7872.d0*t4)+eta10*(-11472.d0+53952.d0
     4   *t2-24960.d0*t4)+eta12*(-4672.d0+30528.d0*t2-23040.d0
     5   *t4)) / cfi

c-------delta isometrisk bredde---------------------------------------

      dq= -b2*dy**2 +b4*(dy/rn)**4 -b6*(dy/rn)**6 +b8*(dy/rn)**8

c-------delta geografisk bredde---------------------------------------

      dbr= d1*dq +d2*dq**2 +d3*dq**3 +d4*dq**4 +d5*dq**5

c-------delta geografisk lengde---------------------------------------

      dlra= b1*dy -b3*(dy/rn)**3 +b5*(dy/rn)**5 -b7*(dy/rn)**7

c-------geografisk bredde og lengde----------------------------------

      glat = fi1 + dbr
      glon = glz + dlra

      return
      end

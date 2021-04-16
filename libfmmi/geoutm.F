      subroutine geoutm ( glat, glon, izone, rnorth, reast )
c
c  NAME:
c     geoutm
c
c  SYNOPSIS:
c     subroutine geoutm ( glat, glon, izone, rnorth, reast )
c     integer          izone
c     double precision glat, glon, rnorth, reast
c
c-----------------------------------------------------------------------
c     transformasjon fra geografisk til plane utm - koord
c-----------------------------------------------------------------------
c     algoritmen baserer seg p]    1) holsen : geodesi 2 ,kap. 4
c                                  2) ngo    : publikasjon nr 16
c     det er tatt med a- og b-koeff. til 8. ledd i rekken
c     programmerer :   aage hovden
c     spr]k        :   vax-11 fortran
c     sentralmeridian utledes av innlest sone
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c  DNMI/FoU  18.05.1995  Anstein Foss ... library version
c-----------------------------------------------------------------------
c

      implicit double precision (a-h,o-z)
      implicit integer          (i-n)

      parameter (a   = 6378388.d0)
      parameter (zs  = 0.9996d0)
      parameter (b   = 6356911.946127946d0)
      parameter (pi  = 3.1415926535898d0)
      parameter (rho = 57.295779513082d0)
      parameter (e2  = 0.006722670062316885d0)
      parameter (em2 = 0.0067681702377508769d0)

c-------sonen bestemmer  sentralmeridian (glz)  ------------------

      glz= dfloat (izone*6-183) / rho

c-------bredde-avhengige st\rrelser------------------------------

      w=dsqrt( 1.d0 - e2*dsin(glat)*dsin(glat) )
      rc=a*zs/w
      rm=a*(1.d0-e2)/(w**3)
      eta2=em2*dcos(glat)*dcos(glat)
      t=dtan(glat)

c-------beregning av meridianbuelengde-----------------------------

      ba1=(1.d0+e2*(0.75d0+e2*(45.d0/64.d0+e2*(175.d0/256.d0+e2
     1    *(11025.d0/16384.d0+e2*43659.d0/65536.d0)))))

      ba2=0.5d0*e2*(.75d0+e2*(15.d0/16.d0+e2*(525.d0/512.d0+e2*
     1    (2205.d0/2048.d0+e2*72765.d0/65536.d0))))

      e4 =e2*e2
      ba3=0.25d0*e4*(15.d0/64.d0+e2*(105.d0/256.d0+e2*(2205.d0/
     1    4096.d0+e2*10395.d0/16384.d0)))

      e6 =e4*e2
      ba4=1.d0/6.d0*e6*(35.d0/512.d0+e2*(315.d0/2048.d0
     1    +e2*31185.d0/131072.d0))

      e8 =e2*e6
      ba5=.125d0*e8*(315.d0/16384.d0+e2*3465.d0/65536.d0)

      e10=e8*e2
      ba6=0.1d0*e10*693.d0/131072.d0

c-----------------------------------------------------------------------

      fi2 =glat+glat
      fi4 =fi2+fi2
      fi6 =fi4+fi2
      fi8 =fi6+fi2
      fi10=fi8+fi2

      g=a*zs*(1.d0-e2)*(ba1*glat-ba2*dsin(fi2)+ba3*dsin(fi4)
     1       -ba4*dsin(fi6)+ba5*dsin(fi8)-ba6*dsin(fi10))

c-------dl=avstand i rad. mellom punkt og akse--------------------------

      dl= glon - glz

c-------eksponenter av t og eta som brukes i rekke-utviklingen----------

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

c------- a-koeffisienter skissert i ngo publikasjon nr 16 --------------

      a01=   rc*dcos(glat)

      a02= - rc*dcos(glat)**2*t/2.d0

      a03= - (rc*dcos(glat)**3)/(6.d0) * (1.d0-t2+eta2)

      a04=    dcos(glat)**4*rc*t/24.d0 * (5.d0-t2+9.d0*eta2+
     1                                 4.d0*eta4)

      a05=  1.d0/120.d0 * rc * dcos(glat)**5 *
     1     ( 5.d0-18.d0*t2+t4+14.d0*eta2-58.d0*t2*eta2
     2     + eta4*(13.d0-64.d0*t2) + eta6*(4.d0-24.d0*t2) )

      a06=  1.d0/720.d0 * rc * dcos(glat)**6 * t *
     1     (-61.d0+58.d0*t2-t4-270.d0*eta2+330.d0*eta2*t2 +
     2     eta4*(-445.d0+680.d0*t2) + eta6*(-324.d0+600.d0*t2) +
     3     eta8*(-88.d0+192.d0*t2) )

      a07=  1.d0/5040.d0 * rc * dcos(glat)**7 *
     1     (-61.d0 +    479.d0 *t2-179.d0*t4  + t6 +
     2     eta2 * (-331.d0+3298.d0*t2-1771.d0*t4) +
     3     eta4 * (-715.d0+8655.d0*t2-6080.d0*t4) +
     4     eta6 * (-769.d0+10964.d0*t2+9480.d0*t4) +
     5     eta8 * (-412.d0+6760.d0*t2-6912.d0*t4) +
     6     eta10* (-88.d0+1632.d0*t2-1920.d0*t4) )

      a08=  1.d0/40320.d0 * rc * dcos(glat)**8 * t *
     1     (1385.d0-3111.d0*t2+543.d0*t4-t6 +
     2     eta2 * ( 10899.d0-  32802.d0*t2 +9219.d0*t4 ) +
     3     eta4 * ( 34419.d0-129087.d0 *t2+ 49644.d0* t4 ) +
     4     eta6 * ( 56385.d0-252084.d0*t2+121800.d0 * t4 ) +
     5     eta8 * ( 50856.d0-263088.d0*t2+151872.d0* t4 ) +
     6     eta10* ( 24048.d0-140928.d0*t2+94080.d0*t4 ) +
     7     eta12* ( 4672.d0- 30528.d0* t2 + 23040.d0* t4 ) )

c-------delta bredde og lengde i meter---------------------------------

      dbm= (-a02 * dl**2 + a04 * dl**4 - a06*dl**6 + a08*dl**8 )

      dlm= ( a01* dl - a03 * dl**3 + a05 * dl**5 - a07 * dl**7 )

c-------utm koordinater------------------------------------------------

      rnorth = g         + dbm
      reast  = 500000.d0 + dlm

      return
      end

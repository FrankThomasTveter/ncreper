      subroutine SetEarthRadius(ir)
      real ir
      REAL :: ecc, Prd          ! Ellipsoid eccentricity, 
!     Earth equatorial radius (m)
      common /lambert/ecc,prd

      prd=ir
      return
      END SUBROUTINE SetEarthRadius

      SUBROUTINE mvGeod2LambConf(latin1, latin2, 
     &     lon0, lat0, x_lat,x_lon,x_ee, x_nn,
     &     irc)

      IMPLICIT none

!
!     PURPOSE
!     Given Geodetic Coordinates,
!     projection parameters and eccentricity,
!     Computes Lambert Conformal Coordinates
!     
!     INPUT
!     latin1, latin2    The two standard parallels
!     lon0, lat0        Long and lat of the reference point (deg)
!     x_lon, x_lat      Geographical coordinates (deg)
!     
!     OUTPUT
!     x_ee, x_mm        Easting and Northing in projected plane
!     irc               Error return code (0=ok)
!
!     ... from Metview Projection.c to Fortran90 ...
!     
      
      REAL :: latin1, latin2, lon0, lat0
      REAL :: x_ee, x_nn 

      REAL :: x_lon, x_lat
      integer :: irc
      logical, parameter :: bdeb=.false.
      
!     Other vars for computation
      
      REAL :: ecc, Prd          ! Ellipsoid eccentricity, 
!     Earth equatorial radius (m)
      common /lambert/ecc,prd
      REAL :: ff, nn, rho, rho_0, theta
      REAL :: latin1r, latin2r, lon0r, lat0r
      REAL :: x_lonr, x_latr
      REAL, PARAMETER :: pi = 3.14159265358979323846
      REAL :: m1, m2, aux0, aux1, aux2, t1, t2, t0
      REAL :: equad
      REAL :: aux, t, flt
      real ir
      
!     Degrees to radians
      
      latin1r = latin1*pi /180.0D0
      latin2r = latin2*pi /180.0D0
      lon0r = lon0*pi /180.0D0
      lat0r = lat0*pi /180.0D0
      x_lonr = x_lon*pi /180.0D0
      x_latr = x_lat*pi /180.0D0
      
      if (bdeb) then
         write (*,*) '================================='//
     &        '======================'
         write (*,*) '1st parallel (deg)' , latin1 , ' (rad) ' , latin1r
         write (*,*) '2nd parallel (deg)' , latin2 , ' (rad) ' , latin2r
         write (*,*) 'Lon ref (deg)' , lon0 , ' (rad) ' , lon0r
         write (*,*) 'Lat ref (deg)' , lat0 , ' (rad) ' , lat0r 
      
      write(*,*) '====================================='//
     &     '=================='
      write(*,*) ' '
      
      end if
!
      if (prd.le.0.0D0) Prd=6367.47D3
      ecc=0.0D0

C      flt=1.0D0/297.0D0
C      Prd=6378.16D3
C      ecc=sqrt(2.0D0*flt - flt*flt) ! ecc**2 = 2 flt - flt**2
      
!     Computations
      
      equad = ecc**2
      m1 = cos(latin1r)/sqrt(1.0D0-equad *(sin(latin1r)**2)) 
      m2 = cos(latin2r)/sqrt(1.0D0-equad *(sin(latin2r)**2)) 
      aux1 = sqrt((1.0D0-ecc*sin(latin1r))/(1.0D0+ecc*sin(latin1r)))
      aux2 = sqrt((1.0D0-ecc*sin(latin2r))/(1.0D0+ecc*sin(latin2r)))
      aux0 = sqrt((1.0D0-ecc*sin(lat0r))/(1.0D0+ecc*sin(lat0r)))
      t1 =((1.0D0-tan(latin1r/2.0D0))/(1.0D0+tan(latin1r/2.0D0)))/
     &     (aux1**ecc)
      t2 =((1.0D0-tan(latin2r/2.0D0))/(1.0D0+tan(latin2r/2.0D0)))/
     &     (aux2**ecc)
      t0 =((1.0D0-tan(lat0r/2.0D0))/(1.0D0+tan(lat0r/2.0D0)))/
     &     (aux0**ecc)
      

      IF ( latin1r .eq. latin2r ) THEN
         nn = sin (latin1r)
      ELSE
         nn = (log (m1) - log (m2)) / (log (t1) - log (t2))
      ENDIF

C      nn =(log(m1)-log(m2))/(log(t1)-log(t2))

      ff = m1/(nn *(t1**nn))
      rho_0 = Prd*ff *(t0**nn)
      
      aux = sqrt((1.0D0-ecc*sin(x_latr))/(1.0D0+ecc*sin(x_latr)))
      t =((1.0D0-tan(x_latr/2.0D0))/(1.0D0+tan(x_latr/2.0D0)))/
     &     (aux**ecc)
      rho = Prd*ff *(t**nn)
      theta = nn *(x_lonr-lon0r)
      
      x_ee = rho*sin(theta)
      x_nn = rho_0-rho*cos(theta)
      return
      END SUBROUTINE mvGeod2LambConf
      
      SUBROUTINE mvLambConf2Geod(latin1, latin2, 
     &     lon0, lat0, x_ee, x_nn,x_lon, x_lat)
      
      
      IMPLICIT none
      
!     
!     PURPOSE
!     Given Lambert Conformal Coordinates,
!     projection parameters and ellipsoid eccentricity,
!     Computes Coordinates
!     
!     INPUT
!     latin1, latin2    The two standard parallels
!     lon0, lat0        Long and lat of the reference point(deg)
!     x_ee, x_nn        Easting and Northing in projected plane
!     ecc, Prd	     Ellipsoid eccentricity, Earth equatorial radius(m)
!     
!     OUTPUT
!     x_lon, x_lat        Geographical coordinates(deg)
!     
!     ... from Metview Projection.c to Fortran90 ...
!     
      
      REAL :: latin1, latin2, lon0, lat0
      REAL :: x_ee, x_nn 
      REAL :: ecc, Prd
      common /lambert/ecc,prd
      REAL :: x_lon, x_lat
      
!     Other vars for computation
      
      REAL :: ff, nn, rho, rho_0, theta
      REAL :: latin1r, latin2r, lon0r, lat0r
      REAL :: x_lonr, x_latr
      REAL, PARAMETER :: pi   = 3.14159265358979323846
      REAL, PARAMETER :: pi_2 = 1.57079632679489661923
      REAL :: m1, m2, aux0, aux1, aux2, t1, t2, t0
      REAL :: equad
      REAL :: aux, t
      REAL :: aux3, aux4, aux5, xx
      
!     Degrees to radians
!
      if (prd.le.0.0D0) Prd=6367.47D3
      ecc=0.0D0
      
      latin1r = latin1*pi/180.0D0
      latin2r = latin2*pi/180.0D0
      lon0r = lon0*pi/180.0D0
      lat0r = lat0*pi/180.0D0
      
!     Computations
      
      equad = ecc**2
      m1 = cos(latin1r)/sqrt(1.0D0-equad *(sin(latin1r)**2))
      m2 = cos(latin2r)/sqrt(1.0D0-equad *(sin(latin2r)**2))
      aux1 = sqrt((1.0D0-ecc*sin(latin1r))/(1.0D0+ecc*sin(latin1r)))
      aux2 = sqrt((1.0D0-ecc*sin(latin2r))/(1.0D0+ecc*sin(latin2r)))
      aux0 = sqrt((1.0D0-ecc*sin(lat0r))/(1.0D0+ecc*sin(lat0r)))
      t1 =((1.0D0-tan(latin1r/2.0D0))/(1.0D0+tan(latin1r/2.0D0)))/
     &     (aux1**ecc)
      t2 =((1.0D0-tan(latin2r/2.0D0))/(1.0D0+tan(latin2r/2.0D0)))/
     &     (aux2**ecc)
      t0 =((1.0D0-tan(lat0r/2.0D0))/(1.0D0+tan(lat0r/2.0D0)))/
     &     (aux0**ecc)
      
C     write(*,*) 'lamb:',m1,m2,aux1,aux2,aux0,t1,t2,t0


      IF ( latin1r .eq. latin2r ) THEN
         nn = sin (latin1r)
      ELSE
         nn = (log (m1) - log (m2)) / (log (t1) - log (t2))
      ENDIF

C      nn =(log(m1)-log(m2))/(log(t1)-log(t2))

      ff = m1/(nn *(t1**nn))
      rho_0 = Prd*ff *(t0**nn)
      
C      write(*,'(X,5(F10.4,X))') m1,m2,nn,t1,t2


      rho = sign(sqrt(x_ee*x_ee +(rho_0-x_nn) *(rho_0-x_nn)), nn)
      theta = atan(x_ee/(rho_0-x_nn))
      t =(rho/(Prd*ff)) **(1.0D0/nn)
      xx = pi_2-2.0D0*atan(t)
      aux3 = equad/2.0D0+5.0D0*equad*equad/24.0D0+
     &     equad*equad*equad/12.0D0
      aux4 = 7.0D0*equad*equad/48.0D0+29.0D0 *equad*equad*equad/240.0D0
      aux5 =( 7.0D0*equad*equad* equad/120.0D0)*sin(12.0D0*atan(t))
      
      x_lonr = theta/nn+lon0r
      x_latr = xx+aux3*sin(4.0D0*atan(t))-aux4*sin(8.0D0*atan(t))+aux5 
      
!     Radians to Degrees
      
      x_lon = x_lonr*180.0D0/pi
      x_lat = x_latr*180.0D0/pi
      
      
      END SUBROUTINE mvLambConf2Geod
      

      real*8 function degtor(x)
      implicit none
      real*8 x
      real*8 pi
      parameter (pi=3.14159265359)
      degtor=x*pi/180.
      end

      real*8 function rtodeg(x)
      implicit none
      real*8 x
      real*8 pi
      parameter (pi=3.14159265359)
      rtodeg=x*180./pi
      end
      
      real*8 function sindeg(x)
      implicit none
      real*8 x,degtor
      sindeg=sin(degtor(x))
      end

      real*8 function cosdeg(x)
      implicit none
      real*8 x,degtor
      cosdeg=cos(degtor(x))
      end

      real*8 function tandeg(x)
      implicit none
      real*8 x,degtor
      tandeg=tan(degtor(x))
      end

      real*8 function asindeg(x)
      implicit none
      real*8 x,rtodeg
      asindeg=rtodeg(asin(x))
      end

      real*8 function acosdeg(x)
      implicit none
      real*8 x,rtodeg
      acosdeg=rtodeg(acos(x))
      end

      real*8 function atandeg(x)
      implicit none
      real*8 x,rtodeg
      atandeg=rtodeg(atan(x))
      end

      real*8 function atan2deg(y,x)
      implicit none
      real*8 y,x,rtodeg
      atan2deg=rtodeg(atan2(y,x))
      end


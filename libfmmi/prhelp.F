      subroutine prhelp(iunit,chelp)
c
c  PURPOSE: Print help lines in 'program.input' files.
c
c  METHOD:  Help lines start with the string 'chelp'.
c           The file is first rewinded.
c           The program finds the first help line and displays this
c           and the following lines until a not help line is found.
c           (Max. length of a help line is 128 characters)
c
c  INPUT:   iunit  - integer       - file unit (to read from)
c           chelp  - character*(*) - comment line identifier
c
c-----------------------------------------------------------------------
c  DNMI/FoU  12.11.1992  Anstein Foss
c  DNMI/FoU  29.11.1995  Anstein Foss ... IBM/RS 6000 xl fortran
c-----------------------------------------------------------------------
c
c..input:
      integer       iunit
      character*(*) chelp
c
c..local:
      character*128  cline
c
      lcline=128
c
      lchelp=index(chelp,' ')
      if(lchelp.eq.0) then
        lchelp=len(chelp)
      elseif(lchelp.gt.1) then
        lchelp=lchelp-1
      end if
c
      if(lchelp.gt.lcline) lchelp=lcline
      cline=' '
c
      rewind(iunit)
c
      do while (cline(1:lchelp).ne.chelp(1:lchelp))
        read(iunit,fmt='(a)',iostat=ios,err=80,end=90) cline
      end do
c
      write(6,*)
      do while (cline(1:lchelp).eq.chelp(1:lchelp))
	nchr=lenstr(cline,0)
        if(nchr.gt.lchelp) then
c..don't print the blanks at the end of the line
          write(6,*) cline(lchelp+1:nchr)
        else
          write(6,*)
        end if
        read(iunit,fmt='(a)',iostat=ios,err=80,end=80) cline
      end do
c
   80 write(6,*)
      return
c
   90 write(6,*)
      write(6,*) '  Sorry. No help found.'
      write(6,*)
      return
      end

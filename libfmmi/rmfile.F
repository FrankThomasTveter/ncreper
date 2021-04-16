      subroutine rmfile(filnam,iprint,ierror)
c
c  PURPOSE: Remove file.
c           Useful before 'creating' a direct access file in cases
c           where it already exists (and the existing file is larger
c           than the new file will be).
c           The file will not be removed if the input file name
c           is a link.
c
c  INPUT:   filnam    - character*(*)  - the file path and name
c           iprint    - integer        - 0 = do not print message
c                                        1 = print message
c
c  OUTPUT:  ierror    - integer        - 0 = no error (always)
c
c  METHOD:
c    Using fortran library functions (SGI): stat,lstat,system.
c    ( call system('rm -f filnam') ).
c
c  BUGS:
c    Max lenght of file path and name is 256 characters.
c    Always returning 'no error'.
c
c  NOTE:
c    If problems when porting this routine to another machine:
c    Just make it a dummy routine, and possibly remove (output) files
c    before running the programs.
c
c-----------------------------------------------------------------------
c  DNMI/FoU  04.01.1993  Anstein Foss
c  DNMI/FoU  03.05.2000  Anstein Foss ... very simple linux version
c-----------------------------------------------------------------------
c
c..input/output:
      integer       iprint,ierror
      character*(*) filnam
c
c..local:
      integer       i,l
      character*262 cmd
c
      l=index(filnam,' ')-1
      if(l.lt.1) l=len(filnam)
c
         if(iprint.eq.1) then
            write(6,*) 'Remove file:'
            write(6,*)  filnam(1:l)
          end if
          cmd = 'rm -f '//filnam
          call system(cmd)
c
      ierror=0
c
      return
      end

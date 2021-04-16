      subroutine rcomnt(iunit,cexit,comnt,nlines,ierror)
c
c  PURPOSE: Read comment lines in 'program.input' files.
c
c  METHOD:  Comment lines are read until a line starting with
c           the string 'cexit' is found. If 'comnt' is not space,
c           comment lines before the last must be identified with
c           the string 'comnt'.
c
c  INPUT:   iunit  - integer       - file unit (to read from)
c           cexit  - character*(*) - last comment line identifier,
c                                    no spaces allowed,
c                                    max 16 characters
c           comnt  - character*(*) - comment line identifier,
c                                    no check on comment lines if
c                                    first character is space,
c                                    max 16 characters
c           nlines - integer       - no. of lines already read
c
c  OUTPUT:  nlines - integer       - no. of lines read (updated)
c           ierror - integer       - 0 = no error
c                                    1 = error (syntax or read)
c
c-----------------------------------------------------------------------
c  DNMI/FoU  14.10.1992  Anstein Foss
c-----------------------------------------------------------------------
c
c..input:
      integer       iunit
      character*(*) cexit,comnt
c
c..output:
      integer       nlines
      integer       ierror
c
c..local:
      character*16  cline
c
      lcline=16
c
      lcexit=index(cexit,' ')
      if(lcexit.eq.0) then
        lcexit=len(cexit)
      elseif(lcexit.gt.1) then
        lcexit=lcexit-1
      endif
c
      lcomnt=index(comnt,' ')
      if(lcomnt.eq.0) then
        lcomnt=len(comnt)
      else
        lcomnt=lcomnt-1
      endif
c
      if(lcexit.gt.lcline) lcexit=lcline
      if(lcomnt.gt.lcline) lcomnt=lcline
c
      if(lcomnt.eq.0) then
c
        do while (.true.)
          nlines=nlines+1
          read(iunit,fmt='(a)',iostat=ios,err=90,end=90) cline
          if(cline(1:lcexit).eq.cexit(1:lcexit)) goto 80
        end do
c
      else
c
        do while (.true.)
          nlines=nlines+1
          read(iunit,fmt='(a)',iostat=ios,err=90,end=90) cline
          if(cline(1:lcexit).eq.cexit(1:lcexit)) goto 80
          if(cline(1:lcomnt).ne.comnt(1:lcomnt)) goto 90
        end do
c
      end if
c
   80 ierror=0
      return
c
   90 ierror=1
      return
      end

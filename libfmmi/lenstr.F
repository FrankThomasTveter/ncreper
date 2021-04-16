      integer function lenstr(text,minlen)
c
c  NAME:
c     lenstr
c
c  PURPOSE:
c     Return actual length of text string, excl. trailing blanks.
c
c  SYNOPSIS:
c     integer function lenstr(text,minlen)
c     character*(*) text
c     integer       minlen
c
c  INPUT:
c     text   - text string
c     minlen - the minimum length returned (see NOTES below)
c
c  OUTPUT:
c     lenstr - the string length
c
c  NOTES:
c     Use minlen=1 if you use lenstr to print a text string without
c     trailing blanks. E.g. filenames declared with a large maximum
c     length which hardly ever is longer than one line, but always
c     is printed with trailing blanks on two ore more lines.
c     Example:  write(6,*) filename(1:lenstr(filename,1))
c     (filename(1:0) may abort the program if filename=' ')
c
c-----------------------------------------------------------------------
c  DNMI/FoU  15.03.1995  Anstein Foss
c-----------------------------------------------------------------------
c
      character*(*) text
      integer       minlen
c
      lt=len(text)
      l=0
      do k=1,lt
        if(text(k:k).ne.' ') l=k
      end do
c
      if(l.lt.minlen) l=minlen
c
      lenstr=l
c
      return
      end

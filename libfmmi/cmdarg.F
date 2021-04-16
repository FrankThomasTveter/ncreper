      subroutine cmdarg(nopt,copt,iopt,iopts,margs,nargs,cargs,
     +                  mispec,ispec,mrspec,rspec,mcspec,cspec,
     +                                           ierror,nerror)
c
************* MAN_BEGIN ************************************************
*
*  PURPOSE:
*     Search for specifications in command line arguments.
*     All options must start with a - character.
*     Each option consists of a one character identifier and
*     may be followed by further specifications as
*     integer(s), real(s), character or character string.
*     Command line arguments not preceeded by an option identifier
*     are also recorded.
*     No error messages are printed.
*     The intention is to make fortran programs work the UNIX way.
*     At least programs requiring a modest amount of user input.
*
*  SYNOPSIS:
*     integer       nopt
*     character*1   copt(nopt)
*     integer       iopt(nopt)
*     integer       iopts(2,nopt)
*     integer       margs
*     integer       nargs
*     character*(*) cargs(margs)
*     integer       mispec
*     integer       ispec(mispec)
*     integer       mrspec
*     real          rspec(mrspec)
*     integer       mcspec
*     character*(*) cspec(mcspec)
*     integer       ierror
*     integer       nerror
*
*     call cmdarg(nopt,copt,iopt,iopts,margs,nargs,cargs,
*    +            mispec,ispec,mrspec,rspec,mcspec,cspec,
*    +                                     ierror,nerror)
*
*  INPUT DESCRIPTION:
*     nopt        - no. of options to look for
*     copt        - option identifier (one character)
*     iopt(nopt)  - option type, 'arguments' following the identifier:
*                                   0 = none
*                                   1 = one or more integers
*                                   2 = one or more reals
*                                   3 = one character
*                                   4 = one character string
*                                  (other values will be treated as 0)
*     margs       - max. no. of not option type command line arguments
*                   (e.g. file names not preceeded by option identifier)
*     mispec      - max. no. of integers to be returned
*     mrspec      - max. no. of reals to be returned
*     mcspec      - max. no. of character strings to be returned
*                   (note that length of each string also is important)
*
*  OUTPUT DESCRIPTION:
*     iopts(2,nopt) - iopts(1,n) = 0 : option not found
*                     iopts(1,n) > 0 : option found,
*                                      for option type 1, 2 and 4 this
*                                      is the number of integers, reals
*                                      or characters specified.
*                     iopts(2,n)     : pointer to the ispec, rspec or
*                                      cspec arrays, if more than one
*                                      integer or real those will be
*                                      found in the following array
*                                      elements (while all characters
*                                      will be stored as one string)
*     nargs         - no. of not option type command line arguments
*                     found
*     cargs( )      - the not option type command line arguments
*     ispec( )      - the integers found
*     rspec( )      - the reals
*     cspec( )      - the characters and character strings found
*     ierror        - 0 = no error
*                     1 = no command line arguments
*                     2 = error detected in a command line argument
*                     3 = option found more than once.
*                     4 = error due to too small maximum sizes
*                         (array size or character string length)
*                         or possibly too much input on command line.
*                     5 = error due to too long command line argument
*     nerror        - command line argument no. with (first) error
*                     (possibly caused by an error in previous argument)
*
*  METHOD:
*     If a command line argument ends with an option identifier
*     requiring an 'argument', the next command line argument is used.
*     Integer and real numbers:
*       One text string input, i.e. no spaces.
*       Separator between each number is the , character.
*       String is terminated by a not acceptable number character or
*       end of command line argument.
*       Acceptable characters are -+0123456789 and . for reals.
*       A real number doesn't need the . character.
*       Integers can maximum have 9 digits.
*     Character strings:
*       Terminated by end of command line argument.
*
*  EXAMPLES:
*     Accepted input (simplest understandable way first):
*     Option type = 0 : -d
*     Option type = 1 : -d +25     -d +25    -abcd +25    -abcd+25x
*     Option type = 1 : -d 1,2,3   -d1,2,3   -abcd1,2,3x
*     Option type = 2 : -d -1.5    -d-1.5    -abcd-1.5x
*     Option type = 2 : -d -1.5,2.5          -d-1.5,2.5
*     Option type = 3 : -dx        -d x      -dxabc       -abcd x
*     Option type = 4 : -dtext     -d text   -abcdtext    -abcd text
*
*  LIMITATIONS:
*     Each option identifier cannot be found more than once.
*     Maximum length of one command line argument is 127 characters.
*     The routine terminates search when an error is detected.
*
*  CALLS:
*     Fortran System library: iargc
*                             getarg
*
*  CODED AND UPDATED:
*                     DNMI/FoU  09.10.1992  Anstein Foss
*                     DNMI/FoU  16.11.1994  Anstein Foss
*                     DNMI/FoU  30.11.1994  Anstein Foss
*
************* MAN_END **************************************************
c
c..input/output:
      integer       nopt
      character*1   copt(nopt)
      integer       iopt(nopt)
      integer       iopts(2,nopt)
      integer       margs
      integer       nargs
      character*(*) cargs(margs)
      integer       mispec
      integer       ispec(mispec)
      integer       mrspec
      real          rspec(mrspec)
      integer       mcspec
      character*(*) cspec(mcspec)
      integer       ierror
      integer       nerror
c
c..local:
      integer       narg,iarg,nispec,nrspec,ncspec
      integer       lcargs,lcspec,lcarg
      integer       n,no,k,k1,k2
      character*128 carg,cread
      logical       search
c
c..functions:
      integer       iargc
c
      do n=1,nopt
        iopts(1,n)=0
        iopts(2,n)=0
      end do
      nargs=0
c
      narg=iargc()
      if(narg.lt.1) goto 910
c
      lcargs=len(cargs(1))
      lcspec=len(cspec(1))
c
      nispec=0
      nrspec=0
      ncspec=0
      iarg=0
c
      do while (iarg.lt.narg)
c
        iarg=iarg+1
        call getarg(iarg,carg)
        lcarg=index(carg,' ')-1
c..the space is not empty if no space is found,
c..but search for numbers terminates by checking the next character
c..(wich must exist)
        if(lcarg.lt.0) goto 950
c
        if(carg(1:1).ne.'-') then
c..not option type argument
          nargs=nargs+1
          if(lcarg.gt.lcargs .or. nargs.gt.margs) goto 940
          cargs(nargs)=carg
        else
c..start of option
          k=2
c
          do while (k.le.lcarg)
c
            no=0
c..find option identifier
            do n=1,nopt
              if(carg(k:k).eq.copt(n)) no=n
            end do
            if(no.eq.0) goto 910
c..check that it's not set before
            if(iopts(1,no).ne.0) goto 930
c
            if(iopt(no).lt.1 .or. iopt(no).gt.4) then
c..option without further specifications
              iopts(1,no)=1
              k=k+1
            else
c..option with further specifications, possibly get next argument
              k=k+1
              if(k.gt.lcarg) then
                if(iarg.eq.narg) goto 920
                iarg=iarg+1
                call getarg(iarg,carg)
                lcarg=index(carg,' ')-1
                if(lcarg.lt.0) goto 950
                k=1
              end if
              if(iopt(no).eq.1) then
c..integer(s)
                iopts(2,no)=nispec+1
                search=.true.
                do while (search)
                  k1=k
                  if(carg(k:k).eq.'-' .or. carg(k:k).eq.'+') k=k+1
                  k2=k
                  do while (k.le.lcarg .and. carg(k:k).ge.'0'
     +                                 .and. carg(k:k).le.'9')
                    k=k+1
                  end do
                  if(k.eq.k2 .or. k-k2.gt.9) goto 920
                  nispec=nispec+1
                  if(nispec.gt.mispec) goto 940
                  cread=carg(k1:k-1)
                  read(cread,*,iostat=ios,err=920) ispec(nispec)
                  iopts(1,no)=iopts(1,no)+1
                  if(k.le.lcarg .and. carg(k:k).eq.',') then
                    k=k+1
                  else
                    search=.false.
                  end if
                end do
              elseif(iopt(no).eq.2) then
c..real(s)
                iopts(2,no)=nrspec+1
                search=.true.
                do while (search)
                  k1=k
                  if(carg(k:k).eq.'-' .or. carg(k:k).eq.'+') k=k+1
                  k2=k
                  do while (k.le.lcarg .and. carg(k:k).ge.'0'
     +                                 .and. carg(k:k).le.'9')
                    k=k+1
                  end do
                  if(k.le.lcarg .and. carg(k:k).eq.'.') then
                    k=k+1
                    do while (k.le.lcarg .and. carg(k:k).ge.'0'
     +                                   .and. carg(k:k).le.'9')
                      k=k+1
                    end do
                  end if
                  if(k.eq.k2) goto 920
                  nrspec=nrspec+1
                  if(nrspec.gt.mrspec) goto 940
                  cread=carg(k1:k-1)
                  read(cread,*,iostat=ios,err=920) rspec(nrspec)
                  iopts(1,no)=iopts(1,no)+1
                  if(k.le.lcarg .and. carg(k:k).eq.',') then
                    k=k+1
                  else
                    search=.false.
                  end if
                end do
              elseif(iopt(no).eq.3) then
c..one character
                ncspec=ncspec+1
                if(ncspec.gt.mcspec) goto 940
                cspec(ncspec)=carg(k:k)
                iopts(1,no)=1
                iopts(2,no)=ncspec
                k=k+1
              elseif(iopt(no).eq.4) then
c..character string
                ncspec=ncspec+1
                if(ncspec.gt.mcspec .or. lcarg-k+1.gt.lcspec) goto 940
                cspec(ncspec)=carg(k:lcarg)
                iopts(1,no)=lcarg-k+1
                iopts(2,no)=ncspec
                k=lcarg+1
              end if
            end if
c
          end do
c
        end if
c
      end do
c
      ierror=0
      nerror=0
      return
c
  910 ierror=1
      nerror=0
      return
c
  920 ierror=2
      nerror=iarg
      return
c
  930 ierror=3
      nerror=iarg
      return
c
  940 ierror=4
      nerror=iarg
      return
c
  950 ierror=5
      nerror=iarg
      return
c
      end

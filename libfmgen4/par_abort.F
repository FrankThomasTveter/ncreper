      subroutine par_abort(id,len)
c
c     KSM: Abort to take care of problems with calling 
c     abort on MPI systems
c     
      implicit none
      character(*) id
      integer len
#ifdef MPILIB
#include "mpif.h"
      integer errormpi
#endif

      write (6,*)'ABORTING ',id(1:len)
#ifdef NEC
      stop
#else
#ifdef MPILIB 
      call mpi_abort(MPI_COMM_WORLD,errormpi)
#else
      call abort
#endif
#endif
      end

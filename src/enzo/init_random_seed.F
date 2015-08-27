#include "fortran.def"

c=======================================================================
c/////////////////////  SUBROUTINE F_ERROR  \\\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine init_random_seed(in_seed)
c
c     Initialize FORTRAN RNG stream
c
c
c     INPUTS:
c
c     in_seed    - Integer seed for the RNG.
c=======================================================================

      implicit none
#include "fortran_types.def"

      INTG_PREC in_seed

      integer*4 i, n
      integer*4, allocatable :: seed(:)

      call random_seed(size=n)
      allocate(seed(n))
      do i = 1, n
c        Offset by one to avoid the poor RNG stream generated
c        when the seed is set to 0.
         seed(i) = in_seed+1
      end do
      call random_seed(put=seed)
      return
      end

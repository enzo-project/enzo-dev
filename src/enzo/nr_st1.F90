      subroutine nr_st1(x, n1, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: n1, idir
      CMPLX_PREC :: x(n1)

      INTG_PREC :: n(3)
      R_PREC :: factor

      factor = 1.0_RKIND/REAL(n1,RKIND)

      n(1) = n1
      n(2) = 1
      n(3) = 1

      call fourn(x, n, 1_IKIND, idir)

      if( idir == 1 ) then
        x = x * factor
      end if

      return
      end

#include "error.def"
#include "fortran.def"

      subroutine wrapper1d(x, rank, n1, n2, n3, dir, method)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, dir
      CMPLX_PREC :: x(n1)
      external :: method

      INTG_PREC :: n(3)

      if( rank /= 1 ) then
        write(0,*) '1D wrapper rank != 1'
        ERROR_MESSAGE
      end if

      if( n2 /= 1 ) then
        write(0,*) '1D wrapper dim2 > 1'
        ERROR_MESSAGE
      end if

      if( n3 /= 1 ) then
        write(0,*) '1D wrapper dim3 > 1'
        ERROR_MESSAGE
      end if

      n(1) = n1
      n(2) = 1
      n(3) = 1

      call fftwrap1d( x(1), n, dir, method )

      return

      end


      subroutine fftwrap1d( a, n, dir, method )

      implicit none
#include "fortran_types.def"

      CMPLX_PREC :: a(*)
      INTG_PREC :: n(3)
      INTG_PREC :: dir
      external :: method

      call method(a, n(1), dir)

      return

      end

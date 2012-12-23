      subroutine nr_3d(x, rank, n1, n2, n3, dir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, dir
      CMPLX_PREC :: x(n1,n2,n3)

      INTG_PREC :: n(3)
      R_PREC :: factor

!     write(*,*) 'NR_3D ',rank,n1,n2,n3,dir

      if( rank /= 3 ) then
        write(0,*) 'NR_3D rank != 3'
        call stop_all_cpus
      end if

      n(1) = n1
      n(2) = n2
      n(3) = n3

      factor = 1._RKIND/REAL(n1*n2*n3,RKIND)

      if( dir == -1 ) then
        call fourn(x, n, rank, dir)
      else
        call fourn(x, n, rank, dir)
        x = x * factor
      end if

      return
      end

#include "error.def"

      subroutine s66_1d(x, rank, n1, n2, n3, dir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, dir
      CMPLX_PREC :: x(n1)

      INTG_PREC :: n(3)

!     write(*,*) 'S66_1D ',rank,n1,n2,n3,dir

      if( rank /= 1 ) then
        write(0,*) 'S66_1D rank != 1'
        ERROR_MESSAGE
      end if

      if( n2 /= 1 ) then
        write(0,*) 'S66_1D dim2 > 1'
        ERROR_MESSAGE
      end if

      if( n3 /= 1 ) then
         write(0,*) 'S66_1D dim3 > 1'
         ERROR_MESSAGE
      end if

      n(1) = n1
      n(2) = 1
      n(3) = 1

      call sf66( x(1), n1, dir )

      return

      end


      subroutine sf66( a, n, dir )

      implicit none
#include "fortran_types.def"

      INTG_PREC :: n, dir
      R_PREC :: a(2*n)
      INTG_PREC :: i

      if( dir == -1 ) then
        call fft66(a(1),a(2),n,n,n,-2)
!       do i = 1,2*n
!       a(i) = a(i) * sqrt(REAL(n,RKIND))
!       end do
      else
        call fft66(a(1),a(2),n,n,n,2)
        do i = 1,2*n
        a(i) = a(i) / REAL(n,RKIND)
        end do
      end if

      return
      end

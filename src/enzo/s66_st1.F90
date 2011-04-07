      subroutine s66_st1( a, n, dir )

      implicit none
#include "fortran_types.def"

      integer :: n, dir
      R_PREC :: a(2*n)
      integer :: i

      if( dir == -1 ) then
        call fft66(a(1),a(2),n,n,n,-2)
!       do i = 1,2*n
!       a(i) = a(i) * sqrt(REAL(n))
!       end do
      else
        call fft66(a(1),a(2),n,n,n,2)
        do i = 1,2*n
        a(i) = a(i) / REAL(n)
        end do
      end if

      return
      end

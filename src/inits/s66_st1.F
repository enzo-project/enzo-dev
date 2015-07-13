      subroutine s66_st1( a, n, dir )

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: n, dir
      R_PREC :: a(2*n)
      INTG_PREC :: i

      if( dir == -1 ) then
        call fft66(a(1),a(2),n,n,n,-2_IKIND)
!       do i = 1,2*n
!       a(i) = a(i) * sqrt(REAL(n,RKIND))
!       end do
      else
        call fft66(a(1),a(2),n,n,n,2_IKIND)
        do i = 1,2*n
        a(i) = a(i) / REAL(n,RKIND)
        end do
      end if

      return
      end

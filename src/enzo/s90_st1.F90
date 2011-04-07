#include "fortran.def"
      subroutine s90_st1( a, n, dir )

      use singleton

      implicit none
#include "fortran_types.def"

      INTG_PREC :: n, dir
      CMPLX_PREC :: a(n)

      INTG_PREC :: i
      INTG_PREC :: sn(1)

      sn(1) = n

      if( dir == -1 ) then
        call fftn(a, sn, inv=.false. )
        do i = 1,n
        a(i) = a(i) * sqrt(REAL(n,RKIND))
        end do
      else
        call fftn(a, sn, inv=.true. )
        do i = 1,n
        a(i) = a(i) / sqrt(REAL(n,RKIND))
        end do
      end if

      return
      end

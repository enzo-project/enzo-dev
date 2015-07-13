#include "../enzo/fortran.def"
      subroutine s90_st1( a, n, dir )

      use singleton

      implicit none
#include "../enzo/fortran_types.def"

      INTEGER(fftintk) :: n, dir
      COMPLEX(fftkind) :: a(n)

      INTEGER(fftintk) :: i
      INTEGER(fftintk) :: sn(1)

      sn(1) = n

      if( dir == -1 ) then
        call fftn(a, sn, inv=.false._fftlogk )
        do i = 1,n
        a(i) = a(i) * sqrt(REAL(n,fftkind))
        end do
      else
        call fftn(a, sn, inv=.true._fftlogk )
        do i = 1,n
        a(i) = a(i) / sqrt(REAL(n,fftkind))
        end do
      end if

      return
      end

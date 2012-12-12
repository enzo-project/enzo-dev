      subroutine ffte_st1(x, n1, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: n1, idir
      CMPLX_PREC :: x(n1)

      R_PREC :: factor
      R_PREC :: scale
      CMPLX_PREC, allocatable :: work(:)

      INTG_PREC :: nwork, jdir
      INTG_PREC :: m1

      m1 = n1
      nwork = n1*2
      jdir = idir

      allocate( work(nwork) )

      call ffte_zfft1d(x, m1, 0_IKIND, work)
      call ffte_zfft1d(x, m1, jdir,    work)

      deallocate( work )

!     factor = 1.0/REAL(n1,RKIND)

!     if( jdir == 1 ) then
!       scale = 1.0_RKIND
!     else
!       scale = factor
!     end if

      return
      end

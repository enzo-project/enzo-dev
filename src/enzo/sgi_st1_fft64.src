#include "error.def"

#ifdef ALTIX

#ifdef CONFIG_BFLOAT_4


      subroutine sgi_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      real*4 :: factor
      real*4 :: scale

      integer*4 :: jdir
      integer*4 :: m1, i0, i1

      real*4 :: work(2*n1)
      real*4 :: table(2*n1+256)
      integer*4 :: isys(0:1)

      isys(0) = 1
      m1 = n1
!     SGI is backwards: inverse = 1
      jdir = -idir
      i0 = 0
      i1 = idir

      factor = 1.0e-00/real(n1)

      if( idir == -1 ) then
        scale = 1.0e-00
      else
        scale = factor
      end if

      call ccfft(i0, m1, scale, x, x, table, work, isys)
      call ccfft(i1, m1, scale, x, x, table, work, isys)

      return
      end

#endif



#ifdef CONFIG_BFLOAT_8

      subroutine sgi_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      real*8 :: factor
      real*8 :: scale

      integer*4 :: jdir
      integer*4 :: m1, i0, i1

      real*8 :: work(2*n1)
      real*8 :: table(2*n1+256)
      integer*4 :: isys(0:1)

      isys(0) = 1
      m1 = n1
!     SGI is backwards: inverse = 1
      jdir = -idir
      i0 = 0
      i1 = idir

      factor = 1.0d-00/real(n1)

      if( idir == -1 ) then
        scale = 1.0d-00
      else
        scale = factor
      end if

!     write(0,'("ZZFFT idir, jdir, scale, factor", 2i4, 2(1pe12.5))')
!    &  idir, jdir, scale, factor

      call zzfft(i0, m1, scale, x, x, table, work, isys)
      call zzfft(i1, m1, scale, x, x, table, work, isys)

      return
      end

#endif

#else

      subroutine sgi_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      write(0,'("SGI stride 1 FFT error")')
      ERROR_MESSAGE
      
      return
      end
      
#endif

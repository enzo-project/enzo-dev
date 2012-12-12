#include "../enzo/fortran.def"
#ifdef CRAYX1

#ifdef CONFIG_BFLOAT_4

      subroutine cray_st1(x, n1, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: n1, idir
      complex :: x(n1)

      real*4 :: factor
      real*4 :: scale

      real*8 :: table(100+8*n1)
      real*4 :: work(8*n1)

      integer*4 :: jdir
      integer*4 :: m1, i0

      m1 = n1
      i0 = 0
      jdir = idir

      factor = 1.0/real(n1)

      if( jdir == -1 ) then
        scale = 1.0
      else
        scale = factor
      end if

      call ccfft(  i0, m1, scale, x, x, table, work, i0)
      call ccfft(jdir, m1, scale, x, x, table, work, i0)

      return
      end

#endif

#ifdef CONFIG_BFLOAT_8

      subroutine cray_st1(x, n1, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: n1, idir
      complex :: x(n1)

      real*8 :: factor
      real*8 :: scale

      real*8 :: table(100+8*n1)
      real*8 :: work(8*n1)

      integer*4 :: jdir
      integer*4 :: m1, i0

      m1 = n1
      jdir = idir
      i0 = 0

      factor = 1.d0/REAL(n1,8)

      if( jdir == -1 ) then
        scale = 1.0d-00
      else
        scale = factor
      end if

      call zzfft(  i0, m1, scale, x, x, table, work, i0)
      call zzfft(jdir, m1, scale, x, x, table, work, i0)

      return
      end

#endif

#else

      subroutine cray_st1(x, n1, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: n1, idir
      complex :: x(n1)

      write(0,'("Dummy Cray X1 1D FFT - error")')
      call stop_all_cpus

      return
      end

#endif

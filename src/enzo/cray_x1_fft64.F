#include "error.def"
#include "fortran.def"
#ifdef CRAYX1

#ifdef CONFIG_BFLOAT_4

      subroutine cray_1d(x, rank, n1, n2, n3, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1)

      REAL*4 :: factor
      REAL*4 :: scale

      REAL*8 :: table(100+8*n1)
      REAL*4 :: work(8*n1)

      integer*4 :: jdir
      integer*4 :: m1, m2, m3, i0

      m1 = n1
      i0 = 0
      jdir = idir

      factor = 1.0/REAL(n1,RKIND)

      if( jdir == -1 ) then
        scale = 1.0_RKIND
      else
        scale = factor
      end if

      call ccfft(  i0, m1, scale, x, x, table, work, i0)
      call ccfft(jdir, m1, scale, x, x, table, work, i0)

      return
      end




      subroutine cray_2d(x, rank, n1, n2, n3, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1,n2)

      REAL*4 :: factor
      REAL*4 :: scale

      REAL*8 :: table(100+n1+n2)
      REAL*4 :: work(2048*max(n1,n2))

      integer*4 :: jdir
      integer*4 :: m1, m2, m3, i0

      m1 = n1
      m2 = n2
      jdir = idir
      i0 = 0

      factor = 1.0/(REAL(n1,RKIND)*REAL(n2,RKIND))

      if( jdir == -1 ) then
        scale = 1.0_RKIND
      else
        scale = factor
      end if

      call ccfft2d(  i0, m1, m2, scale, x, m1, x, m1, table, work, i0)
      call ccfft2d(jdir, m1, m2, scale, x, m1, x, m1, table, work, i0)

      return
      end




      subroutine cray_3d(x, rank, n1, n2, n3, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1,n2,n3)

      REAL*4 :: factor
      REAL*4 :: scale

      REAL*8 :: table(100+n1+n2+n3)
      REAL*4 :: work(2048*max(n1,n2,n3))

      integer*4 :: jdir
      integer*4 :: m1, m2, m3, i0

      m1 = n1
      m2 = n2
      m3 = n3
      jdir = idir
      i0 = 0

      factor = 1.0/(REAL(n1,RKIND)*REAL(n2,RKIND)*REAL(n3,RKIND))

      if( jdir == -1 ) then
        scale = 1.0_RKIND
      else
        scale = factor
      end if

      call ccfft3d(  i0, m1, m2, m3, scale, x, m1, m2, x, m1, m2,
     &            table, work, i0)
      call ccfft3d(jdir, m1, m2, m3, scale, x, m1, m2, x, m1, m2,
     &            table, work, i0)

      return
      end

#endif

#ifdef CONFIG_BFLOAT_8

      subroutine cray_1d(x, rank, n1, n2, n3, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1)

      REAL*8 :: factor
      REAL*8 :: scale

      REAL*8 :: table(100+8*n1)
      REAL*8 :: work(8*n1)

      integer*4 :: jdir
      integer*4 :: m1, m2, m3, i0

      m1 = n1
      jdir = idir
      i0 = 0

      factor = 1.0d-00/REAL(n1,RKIND)

      if( jdir == -1 ) then
        scale = 1.0_RKIND
      else
        scale = factor
      end if

      call zzfft(  i0, m1, scale, x, x, table, work, i0)
      call zzfft(jdir, m1, scale, x, x, table, work, i0)

      return
      end




      subroutine cray_2d(x, rank, n1, n2, n3, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1,n2)

      REAL*8 :: factor
      REAL*8 :: scale

      REAL*8 :: table(100+2*(n1+n2))
      REAL*8 :: work(2048*max(n1,n2))

      integer*4 :: jdir
      integer*4 :: m1, m2, m3, i0

      m1 = n1
      m2 = n2
      jdir = idir
      i0 = 0

      factor = 1.0d-00/(REAL(n1,RKIND)*REAL(n2,RKIND))

      if( jdir == -1 ) then
        scale = 1.0_RKIND
      else
        scale = factor
      end if

      call zzfft2d(  i0, m1, m2, scale, x, m1, x, m1, table, work, i0)
      call zzfft2d(jdir, m1, m2, scale, x, m1, x, m1, table, work, i0)

      return
      end




      subroutine cray_3d(x, rank, n1, n2, n3, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1,n2,n3)

      REAL*8 :: factor
      REAL*8 :: scale

      REAL*8 :: table(100+2*(n1+n2+n3))
      REAL*8 :: work(2048*max(n1,n2,n3))

      integer*4 :: jdir
      integer*4 :: m1, m2, m3, i0

!     for isys=1
!     R_PREC :: work(4*ncpus*max(n1*n2, n2*n3, n3*n1))

      m1 = n1
      m2 = n2
      m3 = n3
      jdir = idir
      i0 = 0

      factor = 1.0d-00/(REAL(n1,RKIND)*REAL(n2,RKIND)*REAL(n3,RKIND))

      if( jdir == -1 ) then
        scale = 1.0_RKIND
      else
        scale = factor
      end if

      call zzfft3d(  i0, m1, m2, m3, scale, x, m1, m2, x, m1, m2,
     &            table, work, i0)
      call zzfft3d(jdir, m1, m2, m3, scale, x, m1, m2, x, m1, m2,
     &            table, work, i0)

      return
      end

#endif

#else

      subroutine cray_1d(x, rank, n1, n2, n3, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC rank, n1, n2, n3, idir
      CMPLX_PREC x(n1)

      write(0,'("Dummy Cray X1 1D FFT - error")')
      ERROR_MESSAGE

      return
      end

      subroutine cray_2d(x, rank, n1, n2, n3, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC rank, n1, n2, n3, idir
      CMPLX_PREC x(n1,n2)

      write(0,'("Dummy Cray X1 2D FFT - error")')
      ERROR_MESSAGE

      return
      end

      subroutine cray_3d( x, rank, n1, n2, n3, idir )

      implicit none
#include "fortran_types.def"

      INTG_PREC rank, n1, n2, n3, idir
      CMPLX_PREC x(n1,n2,n3)

      write(0,'("Dummy Cray X1 3D FFT - error")')
      ERROR_MESSAGE

      return
      end

#endif

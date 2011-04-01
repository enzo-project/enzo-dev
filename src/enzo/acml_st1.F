#include "error.def"

#ifdef XT3

#ifdef CONFIG_BFLOAT_4

      subroutine acml_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      real*4 :: factor
      real*4 :: scale
      complex*8, allocatable :: work(:)

      integer*4 :: nwork, jdir
      integer*4 :: m1, info, i

      m1 = n1
      nwork = 5*n1+100
      jdir = idir

      allocate( work(nwork) )

      call cfft1d(   0, m1, x, work, info)
      call cfft1d(jdir, m1, x, work, info)

      deallocate( work )

      if( idir == -1 ) then
        do i = 1, n1
          x(i) = x(i) * sqrt(real(n1))
        end do
      else
        do i = 1, n1
          x(i) = x(i) / sqrt(real(n1))
        end do
      end if

      return
      end

#endif

#ifdef CONFIG_BFLOAT_8

      subroutine acml_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      real*8 :: factor
      real*8 :: scale
      complex*16, allocatable :: work(:)

      integer*4 :: nwork, jdir
      integer*4 :: m1, info, i

      m1 = n1
      nwork = 5*n1+100
      jdir = idir

      allocate( work(nwork) )

      call zfft1d(   0, m1, x, work, info)
      if( info .ne. 0 ) write(0,'("Info1 = ",i4)') info
      call zfft1d(jdir, m1, x, work, info)
      if( info .ne. 0 ) write(0,'("Info2 = ",i4)') info


      deallocate( work )

      if( idir == -1 ) then
        do i = 1, n1
          x(i) = x(i) * sqrt(real(n1))
        end do
      else
        do i = 1, n1
          x(i) = x(i) / sqrt(real(n1))
        end do
      end if

      return
      end

#endif

#else

      subroutine acml_st1(x, n1, idir)

      implicit none

      integer n1, idir
      complex x(n1)

      write(0,'("ACML stride 1 FFT error")')
      ERROR_MESSAGE

      return
      end

#endif

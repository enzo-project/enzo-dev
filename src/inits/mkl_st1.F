#include "../enzo/fortran.def"


! Complex to CMPLX_PREC

! isign = 0   initialize coeffts
! isign = -1  forward normal
! isign = +1  inverse normal
! isign = -2  forward normal in, bit-rev out
! isign = +2  inverse input bit-rev, out normal

! cfft1d( r, n, isign, wsave )
! zfft1d( r, n, isign, wsave )

! r(n)  CMPLX_PREC / double CMPLX_PREC
! n     INTG_PREC must be power of 2
! wsave CMPLX_PREC / double CMPLX_PREC  array((3*n)/2)


#ifdef MKL

#ifdef CONFIG_BFLOAT_4

      subroutine mkl_st1(x, n1, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: n1, idir
      CMPLX_PREC :: x(n1)

      integer*8 :: power_of_2

      REAL*4 :: factor
      REAL*4 :: scale
      complex*8, allocatable :: wsave(:)

      integer*4 :: isign, n

      if( power_of_2(n1) .ne. 0 ) then
        write(0,'("Non-power-of-2 in mkl_st1 call")')
        stop 'BAD_FFT'
      end if

      n = n1

      allocate( wsave(3*n/2) ) 
      
      isign = 0
      call cfft1d( x, n, isign, wsave )
      isign = idir
      call cfft1d( x, n, isign, wsave )

      deallocate( wsave )

      return
      end

#endif

#ifdef CONFIG_BFLOAT_8

      subroutine mkl_st1(x, n1, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: n1, idir
      CMPLX_PREC :: x(n1)

      integer*8 :: power_of_2

      REAL*8 :: factor
      REAL*8 :: scale
      complex*16, allocatable :: wsave(:)

      integer*4 :: isign, n

      if( power_of_2(n1) .ne. 0 ) then
        write(0,'("Non-power-of-2 in mkl_st1 call")')
        stop 'BAD_FFT'
      end if

      n = n1

      allocate( wsave(3*n/2) ) 
      
      isign = 0
      call zfft1d( x, n, isign, wsave )
      isign = idir
      call zfft1d( x, n, isign, wsave )

      deallocate( wsave )

      return
      end

#endif

#else

      subroutine mkl_st1(x, n1, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: n1, idir
      CMPLX_PREC :: x(n1)

      write(0,'("MKL stride 1 FFT error")')
      call stop_all_cpus

      return
      end

#endif

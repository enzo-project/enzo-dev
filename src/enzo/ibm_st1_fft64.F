#include "error.def"
#include "fortran.def"

#ifdef SP2

#ifdef CONFIG_BFLOAT_4

      subroutine ibm_st1(x, n1, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: n1, idir
      CMPLX_PREC :: x(n1)

      REAL*4 :: factor
      REAL*4 :: scale
      REAL*8, allocatable :: aux1(:), aux2(:)

      integer*4 :: inc1x, inc2x, inc1y, inc2y
      integer*4 :: naux1, naux2, jdir, mult
      integer*4 :: m1, i0, i1

      m1 = n1
      jdir = -idir
      mult = 1
      inc1x = 1
      inc2x = 0
      inc1y = 1
      inc2y = 0
      i0 = 0
      i1 = 1

! SCFT

      if( n1 <= 8192 ) then
        naux1 = 20000
      else
        naux1 = 20000+(1.14)*n1
      end if

      if( n1 <= 8192 ) then
        naux2 = 20000
      else
        naux2 = 20000+(1.14)*n1
      end if 

      allocate ( aux1(naux1) )
      allocate ( aux2(naux2) )

      factor = 1._RKIND/REAL(n1,RKIND)

      if( jdir == 1 ) then
        scale = 1._RKIND
      else
        scale = factor
      end if

      call scft(i1,x,i1,i0,x,i1,i0,m1,i1,jdir,scale,
     &          aux1,naux1,aux2,naux2)
      call scft(i0,x,i1,i0,x,i1,i0,m1,i1,jdir,scale,
     &          aux1,naux1,aux2,naux2)

      deallocate ( aux1 )
      deallocate ( aux2 )

      return
      end

#endif




#ifdef CONFIG_BFLOAT_8

      subroutine ibm_st1(x, n1, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: n1, idir
      CMPLX_PREC :: x(n1)

      REAL*8 :: factor
      REAL*8 :: scale
      REAL*8, allocatable :: aux1(:), aux2(:)

      integer*4 :: inc1x, inc2x, inc1y, inc2y
      integer*4 :: naux1, naux2, jdir, mult
      integer*4 :: m1, i0, i1

      m1 = n1
      jdir = -idir
      mult = 1
      inc1x = 1
      inc2x = 0
      inc1y = 1
      inc2y = 0
      i0 = 0
      i1 = 1

! DCFT

      if( n1 <= 2048 ) then
        naux1 = 20000
      else
        naux1 = 20000+(2.28)*n1
      end if

      if( n1 <= 2048 ) then
        naux2 = 20000
      else
        naux2 = 20000+(2.28)*n1
      end if

      allocate ( aux1(naux1) )
      allocate ( aux2(naux2) )

      factor = 1._RKIND/REAL(n1,RKIND)

      if( jdir == 1 ) then
        scale = 1._RKIND
      else
        scale = factor
      end if

      call dcft(i1,x,i1,i0,x,i1,i0,m1,i1,jdir,scale,
     &          aux1,naux1,aux2,naux2)
      call dcft(i0,x,i1,i0,x,i1,i0,m1,i1,jdir,scale,
     &          aux1,naux1,aux2,naux2)

      deallocate ( aux1 )
      deallocate ( aux2 )

      return
      end

#endif

#else

      subroutine ibm_st1(x, n1, idir)

      implicit none
#include "fortran_types.def"

      INTG_PREC n1, idir
      CMPLX_PREC x(n1)

      write(0,'("IBM stride 1 FFT error")')
      ERROR_MESSAGE

      return
      end

#endif

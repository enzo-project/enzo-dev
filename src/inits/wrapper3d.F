#include "../enzo/fortran.def"

      subroutine wrapper3d(x, rank, n1, n2, n3, dir, method)

!     All stride-1, cache-friendly permuted wrapper
!     Written by: Robert Harkness
!     Date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

!     Arguments

      INTG_PREC :: rank, n1, n2, n3, dir
      CMPLX_PREC :: x(n1,n2,n3)
      external :: method

!     Locals

      CMPLX_PREC, allocatable :: y(:,:,:), z(:,:,:)
      INTG_PREC :: n(3)
      INTG_PREC :: i,j,k

      if( rank /= 3 ) then
        write(0,*) '3D wrapper rank != 3'
        call stop_all_cpus
      end if

      n(1) = n1
      n(2) = 1
      n(3) = 1

      do k=1,n3
      do j=1,n2
      call fftwrap3d( x(1,j,k), n, dir, method )
      end do
      end do

      allocate( y(n2,n3,n1) )

      call rotate3d(x,n1,n2,n3,y)

      n(1) = n2
      do i=1,n1
      do k=1,n3
      call fftwrap3d( y(1,k,i), n, dir, method )
      end do
      end do

      allocate( z(n3,n1,n2) )

      call rotate3d(y,n2,n3,n1,z)

      deallocate( y)

      n(1) = n3
      do j=1,n2
      do i=1,n1
      call fftwrap3d( z(1,i,j), n, dir, method )
      end do
      end do

      call rotate3d(z,n3,n1,n2,x)

      deallocate( z )

      return

      end


      subroutine fftwrap3d( a, n, dir, method )

      implicit none
#include "../enzo/fortran_types.def"

      CMPLX_PREC :: a(*)
      INTG_PREC :: n(3)
      INTG_PREC :: dir
      external :: method

      call method(a, n(1), dir)

      return

      end

#include "error.def"

      subroutine wrapper3d(x, rank, n1, n2, n3, dir, method)

!     All stride-1, cache-friendly permuted wrapper
!     Written by: Robert Harkness
!     Date:       November, 2003

      implicit none

!     Arguments

      integer :: rank, n1, n2, n3, dir
      complex :: x(n1,n2,n3)
      external :: method

!     Locals

      complex, allocatable :: y(:,:,:), z(:,:,:)
      integer :: n(3)
      integer :: i,j,k

      if( rank /= 3 ) then
        write(0,*) '3D wrapper rank != 3'
        ERROR_MESSAGE
      end if

      n(1) = n1
      n(2) = 1
      n(3) = 1

!$omp parallel do private(j) schedule(static)
      do k=1,n3
      do j=1,n2
      call fftwrap3d( x(1,j,k), n, dir, method )
      end do
      end do

      allocate( y(n2,n3,n1) )

      call rotate3d(x,n1,n2,n3,y)

      n(1) = n2
!$omp parallel do private(k) schedule(static)
      do i=1,n1
      do k=1,n3
      call fftwrap3d( y(1,k,i), n, dir, method )
      end do
      end do

      allocate( z(n3,n1,n2) )

      call rotate3d(y,n2,n3,n1,z)

      deallocate( y)

      n(1) = n3
!$omp parallel do private(i) schedule(static)
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

      complex :: a(*)
      integer :: n(3)
      integer :: dir
      external :: method

      call method(a, n(1), dir)

      return

      end

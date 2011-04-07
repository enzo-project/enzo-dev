#include "error.def"

      subroutine s66_3d(x, rank, n1, n2, n3, dir)

      implicit none

      INTG_PREC :: rank, n1, n2, n3, dir
      CMPLX_PREC :: x(n1,n2,n3)

      CMPLX_PREC, allocatable :: y(:,:,:), z(:,:,:)
      INTG_PREC :: n(3)
      INTG_PREC :: i,j,k

!     write(*,*) 'S66_3D ',rank,n1,n2,n3,dir

      if( rank /= 3 ) then
        write(0,*) 'S66_3D rank != 3'
        ERROR_MESSAGE
      end if

      n(1) = n1
      n(2) = 1
      n(3) = 1

      do k=1,n3
      do j=1,n2
      call sf66( x(1,j,k), n1, dir )
      end do
      end do

      allocate( y(n2,n3,n1) )

!     do i=1,n1
!     do k=1,n3
!     do j=1,n2
!     y(j,k,i) = x(i,j,k)
!     end do
!     end do
!     end do

      call rotate3d(x,n1,n2,n3,y)

!     deallocate( x )

      n(1) = n2
      do i=1,n1
      do k=1,n3
      call sf66( y(1,k,i), n2, dir )
      end do
      end do

      allocate( z(n3,n1,n2) )

!     do j=1,n2
!     do i=1,n1
!     do k=1,n3
!     z(k,i,j) = y(j,k,i)
!     end do
!     end do
!     end do

      call rotate3d(y,n2,n3,n1,z)

      deallocate( y)

      n(1) = n3
      do j=1,n2
      do i=1,n1
      call sf66( z(1,i,j), n3, dir )
      end do
      end do

!     allocate( x(n1,n2,n3) )

!     do k=1,n3
!     do j=1,n2
!     do i=1,n1
!     x(i,j,k) = z(k,i,j)
!     end do
!     end do
!     end do

      call rotate3d(z,n3,n1,n2,x)

      deallocate( z )

      return

      end


      subroutine sf66( a, n, dir )

      implicit none

      INTG_PREC :: n, dir
      R_PREC :: a(2*n)
      INTG_PREC :: i

      if( dir == -1 ) then
        call fft66(a(1),a(2),n,n,n,-2)
!       do i = 1,2*n
!       a(i) = a(i) * sqrt(REAL(n,RKIND))
!       end do
      else
        call fft66(a(1),a(2),n,n,n,2)
        do i = 1,2*n
        a(i) = a(i) / REAL(n,RKIND)
        end do
      end if

      return
      end

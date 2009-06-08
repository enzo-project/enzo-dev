#include "error.def"

      subroutine s90_2d(x, rank, n1, n2, n3, dir)

      implicit none

      integer :: rank, n1, n2, n3, dir
      complex :: x(n1,n2)

      complex, allocatable :: y(:,:)
      integer :: n(3)
      integer :: i,j

!     write(*,*) 'S90_2D ',rank,n1,n2,n3,dir

      if( rank /= 2 ) then
        write(0,*) 'S90_2D rank != 2'
        ERROR_MESSAGE
      end if

      if( n3 /= 1 ) then
        write(0,*) 'S90_2D dim3 > 1'
        ERROR_MESSAGE
      end if

      n(1) = n1
      n(2) = 1
      n(3) = 1

      do j=1,n2
      call sf90( x(1,j), n1, dir )
      end do

      allocate( y(n2,n1) )
      call rotate2d(x,n1,n2,y)

      n(1) = n2
      do i=1,n1
      call sf90( y(1,i), n2, dir )
      end do

      call rotate2d(y,n2,n1,x)
      deallocate( y )

      return

      end


      subroutine sf90( a, n, dir )

      use singleton

      implicit none

      integer :: n, dir
      complex :: a(n)

      integer :: i
      integer :: sn(1)

      sn(1) = n

      if( dir == -1 ) then
        call fftn(a, sn, inv=.false. )
        do i = 1,n
        a(i) = a(i) * sqrt(real(n))
        end do
      else
        call fftn(a, sn, inv=.true. )
        do i = 1,n
        a(i) = a(i) / sqrt(real(n))
        end do
      end if

      return
      end

#include "error.def"

      subroutine s90_1d(x, rank, n1, n2, n3, dir)

      implicit none

      integer :: rank, n1, n2, n3, dir
      complex :: x(n1)

      integer :: n(3)

!     write(*,*) 'S90_1D ',rank,n1,n2,n3,dir

      if( rank /= 1 ) then
        write(0,*) 'S90_1D rank != 1'
        ERROR_MESSAGE
      end if

      if( n2 /= 1 ) then
        write(0,*) 'S90_1D dim2 > 1'
        ERROR_MESSAGE
      end if

      if( n3 /= 1 ) then
        write(0,*) 'S90_1D dim3 > 1'
        ERROR_MESSAGE
      end if

      n(1) = n1
      n(2) = 1
      n(3) = 1

      call sf90( x(1), n1, dir )

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

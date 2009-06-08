#include "error.def"

      subroutine nr_3d(x, rank, n1, n2, n3, dir)

      implicit none

      integer :: rank, n1, n2, n3, dir
      complex :: x(n1,n2,n3)

      integer :: n(3)
      real :: factor

!     write(*,*) 'NR_3D ',rank,n1,n2,n3,dir

      if( rank /= 3 ) then
        write(0,*) 'NR_3D rank != 3'
        ERROR_MESSAGE
      end if

      n(1) = n1
      n(2) = n2
      n(3) = n3

      factor = 1.0/real(n1*n2*n3)

      if( dir == -1 ) then
        call fourn(x, n, rank, dir)
      else
        call fourn(x, n, rank, dir)
        x = x * factor
      end if

      return
      end

      subroutine nr_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      integer :: n(3)
      real :: factor

      factor = 1.0/real(n1)

      n(1) = n1
      n(2) = 1
      n(3) = 1

      call fourn(x, n, 1, idir)

      if( idir == 1 ) then
        x = x * factor
      end if

      return
      end

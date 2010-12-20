      subroutine s66_st1( a, n, dir )

      implicit none

      integer :: n, dir
      real :: a(2*n)
      integer :: i

      if( dir == -1 ) then
        call fft66(a(1),a(2),n,n,n,-2)
!       do i = 1,2*n
!       a(i) = a(i) * sqrt(real(n))
!       end do
      else
        call fft66(a(1),a(2),n,n,n,2)
        do i = 1,2*n
        a(i) = a(i) / real(n)
        end do
      end if

      return
      end

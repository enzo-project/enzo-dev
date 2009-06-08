      subroutine s90_st1( a, n, dir )

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

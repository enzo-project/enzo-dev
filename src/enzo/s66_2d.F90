#include "error.def"

      subroutine s66_2d(x, rank, n1, n2, n3, dir)

      implicit none

      INTG_PREC :: rank, n1, n2, n3, dir
      CMPLX_PREC :: x(n1,n2)

      CMPLX_PREC, allocatable :: y(:,:)
      INTG_PREC :: n(3)
      INTG_PREC :: i,j

!     write(*,*) 'S66_2D ',rank,n1,n2,n3,dir

      if( rank /= 2 ) then
        write(0,*) 'S66_2D rank != 2'
        ERROR_MESSAGE
      end if

      if( n3 /= 1 ) then
        write(0,*) 'S66_2D dim3 > 1'
        ERROR_MESSAGE
      end if

      n(1) = n1
      n(2) = 1
      n(3) = 1

      do j=1,n2
      call sf66( x(1,j), n1, dir )
      end do

      allocate( y(n2,n1) )
      call rotate2d(x,n1,n2,y)

      n(1) = n2
      do i=1,n1
      call sf66( y(1,i), n2, dir )
      end do

      call rotate2d(y,n2,n1,x)
      deallocate(y)

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

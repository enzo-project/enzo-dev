#include "error.def"

      subroutine s90_2d(x, rank, n1, n2, n3, dir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, dir
      CMPLX_PREC :: x(n1,n2)

      CMPLX_PREC, allocatable :: y(:,:)
      INTG_PREC :: n(3)
      INTG_PREC :: i,j

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

      INTG_PREC :: n, dir
      CMPLX_PREC :: a(n)

      INTG_PREC :: i
      INTG_PREC :: sn(1)

      sn(1) = n

      if( dir == -1 ) then
        call fftn(a, sn, inv=.false._fftlogk )
        do i = 1,n
        a(i) = a(i) * sqrt(REAL(n,fftkind))
        end do
      else
        call fftn(a, sn, inv=.true._fftlogk )
        do i = 1,n
        a(i) = a(i) / sqrt(REAL(n,fftkind))
        end do
      end if

      return
      end

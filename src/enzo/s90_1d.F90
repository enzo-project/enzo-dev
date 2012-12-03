#include "error.def"

      subroutine s90_1d(x, rank, n1, n2, n3, dir)

      implicit none
#include "fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, dir
      CMPLX_PREC :: x(n1)

      INTG_PREC :: n(3)

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

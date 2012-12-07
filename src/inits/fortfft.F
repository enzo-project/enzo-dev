#include "../enzo/fortran.def"
      subroutine fortfft(x, rank, dim1, dim2, dim3, dir)

!  written by: Robert Harkness
!  date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

!     Arguments

      INTG_PREC :: rank, dim1, dim2, dim3, dir
      CMPLX_PREC :: x(dim1,dim2,dim3)

!     Locals

      INTG_PREC :: method
      INTG_PREC :: powers(3)

!     External function

      INTG_PREC :: power_of_2


      method = 4

      powers(1) = power_of_2(dim1)
      powers(2) = power_of_2(dim2)
      powers(3) = power_of_2(dim3)

      write(*,'("Fortran FFT method ",i1,":",4i8)') 
     &   method, dim1, dim2, dim3, dir

      if( method == 1 )
     &   call nr_3d(x, rank, dim1, dim2, dim3, dir)
      if( method == 2 )
     &   call nrx_3d(x, rank, dim1, dim2, dim3, dir)
      if( method == 3 )
     &   call s90_3d(x, rank, dim1, dim2, dim3, dir)
      if( method == 4 )
     &   call s66_3d(x, rank, dim1, dim2, dim3, dir)

      write(*,'("Fortran FFT complete")')

      return
      end

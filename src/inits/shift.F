#include "../enzo/fortran.def"
!=======================================================================
!////////////////////////  SUBROUTINE SHIFT  \\\\\\\\\\\\\\\\\\\\\\\\\\\

      subroutine shift(field, nx, ny, nz, in, jn, kn, temp)

!  SHIFTS A 3D FIELD BY 1/2 CELL USING CIC_INTERPOLATE
!
!     written by: Greg Bryan
!     date:       June, 1997
!     modified:   Robert Harkness
!     date:       November, 2003
!
!  PURPOSE: 
!
!  INPUTS:
!        i,j,kn      = real dimensions of green
!        nx,ny,nz    = active dimensions of green
!        temp        = temp field
!
!  Outputs:
!        field       = input & output field

      implicit none
#include "../enzo/fortran_types.def"

!  Arguments

      INTG_PREC :: in, jn, kn, nx, ny, nz
      R_PREC    :: field(in, jn, kn), temp(in,jn,kn)

!  Locals

      INTG_PREC :: i, j, k, j1, k1
      R_PREC    :: coef


!     Shift

      coef = 1._RKIND/8._RKIND

      do k=1, nz
         k1 = k+1
         if (k .eq. nz) k1 = 1
         do j=1, ny
            j1 = j+1
            if (j .eq. ny) j1 = 1
            do i=1, nx
               temp(i,j,k) = coef*(
     &             field(i,j ,k) + field(i,j ,k1) +
     &             field(i,j1,k) + field(i,j1,k1) )
            enddo
            do i=1, nx-1
               temp(i,j,k) = temp(i,j,k) + coef*(
     &             field(i+1,j ,k) + field(i+1,j ,k1) +
     &             field(i+1,j1,k) + field(i+1,j1,k1) )
            enddo
            temp(nx,j,k) = temp(nx,j,k) + coef*(
     &             field(1,j ,k) + field(1,j ,k1) +
     &             field(1,j1,k) + field(1,j1,k1) )
         enddo
      enddo

!     Copy back

      do k=1, nz
         do j=1, ny
            do i=1, nx
               field(i,j,k) = temp(i,j,k)
            enddo
         enddo
      enddo

      return
      end

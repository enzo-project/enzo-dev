#include "../enzo/fortran.def"
!=======================================================================
!////////////////////////  SUBROUTINE ADJFFT  \\\\\\\\\\\\\\\\\\\\\\\\\\

      subroutine adjfft(array, nxz, nyz, nzz, in, jn)

!  ADJUSTS A COMPLEX ARRAY TO SATISFY THE CONJUGATE RELATIONS FOR REAL ARRAYS
!
!  written by: Greg Bryan
!  date:       January, 1993
!  modified:   Robert Harkness
!  date:       November, 2003
!
!  INPUTS: array(in,jn,*) - 3d array to be adjusted
!          in,jn        - first and second declared dimension of array
!          nxz,nyz,nzz  - size (within array) of appropriate 3d array
!
!  OUTPUTS: array(in,jn,*)

      implicit NONE
#include "../enzo/fortran_types.def"

!     Arguments

      INTG_PREC :: nxz, nyz, nzz, in, jn
      R_PREC :: array(in, jn, *)

!     Locals

      INTG_PREC :: j, k

!-----------------------------------------------------------------------

!  Adjust the array to satisfy the complex conjugation relations for real
!  arrays.

!  1) adjust corners

      array(    2,       1,       1) = 0.0_RKIND
      array(nxz+2,       1,       1) = 0.0_RKIND
      array(    2, nyz/2+1,       1) = 0.0_RKIND
      array(    2,       1, nzz/2+1) = 0.0_RKIND
      array(nxz+2, nyz/2+1,       1) = 0.0_RKIND
      array(nxz+2,       1, nzz/2+1) = 0.0_RKIND
      array(    2, nyz/2+1, nzz/2+1) = 0.0_RKIND
      array(nxz+2, nyz/2+1, nzz/2+1) = 0.0_RKIND

!  2) adjust faces

      do j = 2, nyz/2
         array(    1, j,       1) =  array(    1, nyz+2-j,       1)
         array(    2, j,       1) = -array(    2, nyz+2-j,       1)
         array(nxz+1, j,       1) =  array(nxz+1, nyz+2-j,       1)
         array(nxz+2, j,       1) = -array(nxz+2, nyz+2-j,       1)
         array(    1, j, nzz/2+1) =  array(    1, nyz+2-j, nzz/2+1)
         array(    2, j, nzz/2+1) = -array(    2, nyz+2-j, nzz/2+1)
         array(nxz+1, j, nzz/2+1) =  array(nxz+1, nyz+2-j, nzz/2+1)
         array(nxz+2, j, nzz/2+1) = -array(nxz+2, nyz+2-j, nzz/2+1)
         do k = 2, nzz/2
            array(    1, j,       k) =  array(    1, nyz+2-j, nzz+2-k)
            array(    2, j,       k) = -array(    2, nyz+2-j, nzz+2-k)
            array(nxz+1, j,       k) =  array(nxz+1, nyz+2-j, nzz+2-k)
            array(nxz+2, j,       k) = -array(nxz+2, nyz+2-j, nzz+2-k)
            array(    1, j, nzz+2-k) =  array(    1, nyz+2-j,       k)
            array(    2, j, nzz+2-k) = -array(    2, nyz+2-j,       k)
            array(nxz+1, j, nzz+2-k) =  array(nxz+1, nyz+2-j,       k)
            array(nxz+2, j, nzz+2-k) = -array(nxz+2, nyz+2-j,       k)
         enddo
      enddo

!  3) adjust sides

      do k = 2, nzz/2
         array(    1,       1, k) =  array(    1,       1, nzz+2-k)
         array(    2,       1, k) = -array(    2,       1, nzz+2-k)
         array(nxz+1,       1, k) =  array(nxz+1,       1, nzz+2-k)
         array(nxz+2,       1, k) = -array(nxz+2,       1, nzz+2-k)
         array(    1, nyz/2+1, k) =  array(    1, nyz/2+1, nzz+2-k)
         array(    2, nyz/2+1, k) = -array(    2, nyz/2+1, nzz+2-k)
         array(nxz+1, nyz/2+1, k) =  array(nxz+1, nyz/2+1, nzz+2-k)
         array(nxz+2, nyz/2+1, k) = -array(nxz+2, nyz/2+1, nzz+2-k)
      enddo

      return
      end

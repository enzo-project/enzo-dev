#include "fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE CALC_EIGEN  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine calc_eigen(dslice, cs, is, ie, idim, lem, rem)
c
c  COMPUTES THE LEFT AND RIGHT EIGENMATRICES FOR ADIABATIC HYDRODYNAMICS
c
c  written by: John Wise
c  date:       June, 2010
c  modified1:  
c
c  PURPOSE:  To project the primitive variables into the characteristic
c    variables and vice-versa, these matrices are needed.  We use these
c    variables to calculate the differences during interpolation.  This
c    makes the interpolation total variation diminishing (LeVeque 2002).
c
c  INPUT:
c    dslice - extracted 2d slice of the density, d
c    cs     - sound speed
c    
c  OUTPUT:
c    lem,rem  - left and right eigenmatrices
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer idim, is, ie
      real dslice(idim), cs(idim)
      real lem(idim,5,5), rem(idim,5,5)
c
c     Locals
c
      integer i, n, m
      real csi, csqi, csq
c
      do m = 1,5
         do n = 1,5
            do i = 1, idim
               lem(i,n,m) = 0.d0
               rem(i,n,m) = 0.d0
            enddo
         enddo
      enddo

      do i = 1, idim
         csi = 1.0/cs(i)
         csq = cs(i)**2
         csqi = 1.d0/csq
c
c     Left eigenvectors
c         
         lem(i,2,1) = -0.5d0*dslice(i)*csi
         lem(i,5,1) = 0.5d0*csqi
         lem(i,1,2) = 1.d0
         lem(i,5,2) = -csqi
         lem(i,3,3) = 1.d0
         lem(i,4,4) = 1.0
         lem(i,2,5) = -lem(i,2,1)
         lem(i,5,5) = lem(i,5,1)
c
c     Right eigenvectors
c
         rem(i,1,1) = 1.d0
         rem(i,1,2) = -cs(i)/dslice(i)
         rem(i,1,5) = csq
         rem(i,2,1) = 1.d0
         rem(i,3,3) = 1.d0
         rem(i,4,4) = 1.d0
         rem(i,5,1) = 1.d0
         rem(i,5,2) = -rem(i,1,2)
         rem(i,5,5) = csq

      enddo

      end

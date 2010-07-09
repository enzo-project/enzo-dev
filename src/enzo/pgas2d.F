c=======================================================================
c////////////////////////  SUBROUTINE PGAS2D  \\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine pgas2d(dslice, eslice, pslice, uslice, vslice,
     &                  wslice, idim, jdim, i1, i2, j1, j2,
     &                  gamma, pmin)
c
c  COMPUTES GAS PRESSURE ON A SLICE
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  July, 1994; by Greg Bryan; changed to a slice
c
c  PURPOSE: Computes gas pressure needed for interpolations for PPM
c           solver.  EOS used is p=(gamma-1)*d*[E-(u**2+v**2+w**2)/2],
c          (eqn 2.1)
c          Only works for gamma law gas with gamma.gt.0
c          All quantities are zone centered (defined over nxz,nyz,nzz).
c
c  INPUTS:
c     dslice  - density slice
c     eslice  - total specific energy slice
c     gamma   - gamma law parameter
c     i1,i2   - starting and ending addresses for dimension 1
c     idim    - declared leading dimension of slices
c     j1,j2   - starting and ending addresses for dimension 2
c     jdim    - declared second dimension of slices
c     pmin    - minimum allowed pressure
c     uslice  - velocity-1 slice
c     vslice  - velocity-2 slice
c     wslice  - velocity-3 slice
c
c  OUTPUT ARGUMENTS: 
c     pslice  - pressure slice
c
c  EXTERNALS: 
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer i1, i2, idim, j1, j2, jdim
      real    gamma, pmin
      real dslice(idim,jdim), eslice(idim,jdim), pslice(idim,jdim),
     &     uslice(idim,jdim), vslice(idim,jdim), wslice(idim,jdim)
c
c  locals
c
      integer i, j
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c  compute the pressure
c
      do j=j1,j2
         do i=i1,i2
            pslice(i,j) =
     &           (gamma-1.0)*dslice(i,j)*(
     &              eslice(i,j) - 0.5*(uslice(i,j)**2 +
     &                                 vslice(i,j)**2 +
     &                                 wslice(i,j)**2  )
     &                                   )
            if (pslice(i,j) .lt. pmin) pslice(i,j) = pmin
         enddo
      enddo
c
      return
      end

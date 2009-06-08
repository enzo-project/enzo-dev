c=======================================================================
c//////////////////////  SUBROUTINE PGAS2D_DUAL  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine pgas2d_dual(dslice, eslice, geslice, pslice, 
     &                       uslice, vslice, wslice, eta1, eta2,
     &                       idim, jdim, i1, i2, j1, j2, gamma, pmin)
c
c  COMPUTES GAS PRESSURE ON A SLICE (DUAL ENERGY VERSION)
c
c  written by: Greg Bryan
c  date:       March, 1994
c  modified1:
c
c  PURPOSE: Computes gas pressure needed for interpolations for PPM
c           solver.  EOS used is p=(gamma-1)*d*[E-(u**2+v**2+w**2)/2],
c           (eqn 2.1) Only works for gamma law gas with gamma > 0.
c           This particular version uses the dual energy formalism
c           (see PPM gravity method paper) to determine the pressure
c           from a combination of the total energy and the gas energy.
c
c  INPUTS:
c     dslice  - density slice
c     eslice  - total specific energy slice
c     eta1    - selection parameter for gas energy (typically ~0.001)
c     eta2    - selection parameter for total energy (typically ~0.1)
c     gamma   - gamma law parameter
c     geslice - specific gas energy slice
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
      real    eta1, eta2, gamma, pmin
      real dslice(idim,jdim), eslice(idim,jdim), pslice(idim,jdim),
     &     uslice(idim,jdim), vslice(idim,jdim), wslice(idim,jdim),
     &    geslice(idim,jdim)
c
c  locals
c
      integer i, j, im1, ip1
      real    demax, ge1, ge2, ke
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c  compute the pressure
c
      do j=j1,j2
         do i=i1,i2
c
c           Compute the specific energy from the total energy
c
            ke = 0.5*(uslice(i,j)**2 + vslice(i,j)**2 + wslice(i,j)**2)
            ge1 = eslice(i,j) - ke
c
c           Find the maximum nearby total energy (not specific)
c
            im1 = max(i-1,i1)
            ip1 = min(i+1,i2)
            demax = max(dslice(i  ,j)*eslice(i  ,j),
     &                  dslice(im1,j)*eslice(im1,j),
     &                  dslice(ip1,j)*eslice(ip1,j))
c
c           If the ratio of the gas energy to the max nearby total energy is
c             is > eta2 then use the gas energy computed from the total energy
c             to update the specific energy (total energy is better).
c
            if (geslice(i,j) .le. 0.0) write(6,*) 
     &            'pga:',i,j,geslice(i,j),demax
            if (ge1*dslice(i,j)/demax .gt. eta2) geslice(i,j) = ge1
            if (geslice(i,j) .le. 0.0) write(6,*) 
     &            'pgb:',i,j,geslice(i,j),dslice(i,j),demax,eslice(i,j)
c
c           If the ratio of the specific gas energy to specific total energy 
c             (in this cell) is < eta1 then use the gas energy to update
c             the specific energy (gas energy is better).
c
            if (ge1/eslice(i,j) .gt. eta1) then
               ge2 = ge1
            else
               ge2 = geslice(i,j)
            endif
c
c           If pressure is below the minimum, set it to the minimum
c
            ge2 = max(ge2, pmin/((gamma - 1.0)*dslice(i,j)))
c
c           Update the total energy
c
            eslice(i,j) = eslice(i,j) - ge1 + ge2
c
c           Compute the pressure with the total energy derived gas energy
c
            pslice(i,j) = (gamma - 1.0)*dslice(i,j)*ge2
c
         enddo
      enddo
c
      return
      end

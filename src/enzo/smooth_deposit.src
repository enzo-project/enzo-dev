#include "error.def"

c=======================================================================
c//////////////////////  SUBROUTINE SMOOTH_DEPOSIT  \\\\\\\\\\\\\\\\\\\\
c
      subroutine smooth_deposit(posx, posy, posz, ndim, npositions, 
     &                      mass, field, leftedge, 
     &                      dim1, dim2, dim3, cellsize, rsmooth)
#ifndef CONFIG_PFLOAT_16
#include "fortran.def"
c
c  PERFORMS 3D SMOOTHED INTERPOLATION FROM FIELD TO SUMFIELD
c
c  written by: Greg Bryan
c  date:       January, 2000
c  modified1:
c
c  PURPOSE: This routine performs a three-dimension, second-order
c           interpolation from field to sumfield (without clearing sumfield
c           first) at the positions specified.
c
c  INPUTS:
c     ndim       - dimensionality
c     cellsize   - the cell size of field
c     dim1,2,3   - real dimensions of field
c     leftedge   - the left edge(s) of field
c     npositions - number of particles
c     posx,y,z   - particle positions
c     sumfield   - 1D field (length npositions) of masses
c
c  OUTPUT ARGUMENTS: 
c     field      - field to be deposited to
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
      integer dim1, dim2, dim3, npositions, ndim
      REALSUB posx(npositions), posy(npositions), posz(npositions),
     &        leftedge(3)
      real    mass(npositions), field(dim1, dim2, dim3)
      real    cellsize, rsmooth
c
c  locals
c     
      integer i, j, k, n
      REALSUB xpos, ypos, zpos, rsmsqr, rsqr
      real    coef, pi
      parameter (pi = 3.14159)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      rsmsqr = rsmooth**2
      coef   = 3.0/pi*(cellsize/rsmooth)**3
c
      if (ndim .ne. 3) then
         write(6,*) 'SMOOTH_DEPOSIT: only ndim=3 supported.'
         ERROR_MESSAGE
      endif
c
c     3D

      if (ndim .eq. 3) then
c
c        loop over grid
c
         do k=1, dim3
            zpos = leftedge(3) + (float(k) - 0.5)*cellsize
            do j=1, dim2
               ypos = leftedge(2) + (float(j) - 0.5)*cellsize
               do i=1, dim1
                  xpos = leftedge(1) + (float(i) - 0.5)*cellsize
c
c                 Loop over particles
c
                  do n=1, npositions
c
c                     Compute distance from particle to cell center
c
                     rsqr = (posx(n) - xpos)**2 +
     &                      (posy(n) - ypos)**2 +
     &                      (posz(n) - zpos)**2
c
                     if (rsqr .lt. rsmsqr)
     &                    field(i,j,k) = field(i,j,k) +
     &                         mass(n)*coef*(1.0 - sqrt(rsqr)/rsmooth)
c
                  enddo
c
               enddo
            enddo
         enddo
c
      endif
c
      return
#endif
      end

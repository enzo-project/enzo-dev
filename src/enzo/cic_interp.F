c=======================================================================
c//////////////////////  SUBROUTINE CIC_INTERP  \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cic_interp(posx, posy, posz, ndim, npositions, 
     &                      sumfield, field, leftedge, 
     &                      dim1, dim2, dim3, cellsize)
#ifndef CONFIG_PFLOAT_16
#include "fortran.def"
c
c  PERFORMS 1/2/3D CLOUD-IN-CELL INTERPOLATION FROM FIELD TO SUMFIELD
c
c  written by: Greg Bryan
c  date:       January, 1998
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
c     field      - field to interpolate from
c     leftedge   - the left edge(s) of field
c     npositions - number of particles
c     posx,y,z   - particle positions
c     sumfield   - 1D field (length npositions) to interpolate into
c
c  OUTPUT ARGUMENTS: 
c     sumfield   - 1D field (length npositions) to interpolate into
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
     &        leftedge(3), cellsize
      real    sumfield(npositions), field(dim1, dim2, dim3)
c
c  locals
c
      integer i1, j1, k1, n
      real    xpos, ypos, zpos, dx, dy, dz, fact
      REALSUB half, edge1, edge2, edge3
c
c  constants (half is slightly more than 1/2 to make sure rounding is
c             performed in the correct direction)
c
      parameter (half = 0.5001)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      fact = 1.0/cellsize
      edge1 = real(dim1) - half
      edge2 = real(dim2) - half
      edge3 = real(dim3) - half
c
c     1D
c
      if (ndim .eq. 1) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
c
c           Convert this into an integer index
c
            i1  = int(xpos + 0.5)
c
c           Compute the weights
c
            dx = real(i1) + 0.5 - xpos
c
c           Interpolate from field into sumfield
c
            sumfield(n) = sumfield(n) +
     &              field(i1  ,1,1)*     dx  +
     &              field(i1+1,1,1)*(1.0-dx)
c
         enddo
c
      endif
c
c     2D
c
      if (ndim .eq. 2) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
            ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
c
c           Convert this into an integer index
c
            i1  = int(xpos + 0.5)
            j1  = int(ypos + 0.5)
c
c           Compute the weights
c
            dx = real(i1) + 0.5 - xpos
            dy = real(j1) + 0.5 - ypos
c
c           Interpolate from field into sumfield
c
            sumfield(n) = sumfield(n) +
     &              field(i1  ,j1  ,1)*     dx *     dy  +
     &              field(i1+1,j1  ,1)*(1.0-dx)*     dy  +
     &              field(i1  ,j1+1,1)*     dx *(1.0-dy) +
     &              field(i1+1,j1+1,1)*(1.0-dx)*(1.0-dy)
c
         enddo
c
      endif
c
c     3D
c
      if (ndim .eq. 3) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
            ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
            zpos = min(max((posz(n) - leftedge(3))*fact, half), edge3)
c
c           Convert this into an integer index
c
            i1  = int(xpos + 0.5)
            j1  = int(ypos + 0.5)
            k1  = int(zpos + 0.5)
c
c           Compute the weights
c
            dx = real(i1) + 0.5 - xpos
            dy = real(j1) + 0.5 - ypos
            dz = real(k1) + 0.5 - zpos
c
c           Interpolate from field into sumfield
c 
    
            sumfield(n) = sumfield(n) +
     &              field(i1  ,j1  ,k1  )*     dx *     dy *     dz  +
     &              field(i1+1,j1  ,k1  )*(1.0-dx)*     dy *     dz  +
     &              field(i1  ,j1+1,k1  )*     dx *(1.0-dy)*     dz  +
     &              field(i1+1,j1+1,k1  )*(1.0-dx)*(1.0-dy)*     dz  +
     &              field(i1  ,j1  ,k1+1)*     dx *     dy *(1.0-dz) +
     &              field(i1+1,j1  ,k1+1)*(1.0-dx)*     dy *(1.0-dz) +
     &              field(i1  ,j1+1,k1+1)*     dx *(1.0-dy)*(1.0-dz) +
     &              field(i1+1,j1+1,k1+1)*(1.0-dx)*(1.0-dy)*(1.0-dz)
c

         enddo
c
      endif
c

      return
#endif
      end

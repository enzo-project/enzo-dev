#include "fortran.def"
#include "error.def"
c=======================================================================
c//////////////////////  SUBROUTINE NGP_DEPOSIT  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine ngp_deposit(posx, posy, posz, ndim, npositions, 
     &                      mass, field, leftedge, 
     &                      dim1, dim2, dim3, cellsize)
#ifndef CONFIG_PFLOAT_16
c
c  PERFORMS 1/2/3D NEAREST-GRID-POINT INTERPOLATION FROM FIELD TO SUMFIELD
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:
c
c  PURPOSE: This routine performs a three-dimension, first-order
c           deposition to field.
c
c  INPUTS:
c     ndim       - dimensionality
c     cellsize   - the cell size of field
c     dim1,2,3   - real dimensions of field
c     leftedge   - the left edge(s) of field
c     npositions - number of particles
c     posx,y,z   - particle positions
c     mass       - 1D field (length npositions) of masses
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
#include "fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC dim1, dim2, dim3, npositions, ndim
      P_PREC posx(npositions), posy(npositions), posz(npositions),
     &        leftedge(3)
      R_PREC    mass(npositions), field(dim1, dim2, dim3), cellsize
c
c  locals
c
      INTG_PREC i1, j1, k1, n
      R_PREC    xpos, ypos, zpos
      P_PREC    edge1, edge2, edge3, fact, half
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c

!     write(0,*) npositions, leftedge, dim1, dim2, dim3, cellsize

      fact = 1._PKIND/cellsize
      half = 0.5001_PKIND
      edge1 = REAL(dim1,PKIND) - half
      edge2 = REAL(dim2,PKIND) - half
      edge3 = REAL(dim3,PKIND) - half
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
c           Convert this into an INTG_PREC index
c
            i1  = int(xpos,IKIND) + 1
c     
c           Add mass to field
c     
            field(i1, 1, 1) = field(i1,1,1) + mass(n)
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
c           Convert this into an INTG_PREC index
c
            i1  = int(xpos,IKIND) + 1
            j1  = int(ypos,IKIND) + 1
c     
c           Add mass to field
c     
            field(i1, j1, 1) = field(i1,j1,1) + mass(n)
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
c           Convert this into an INTG_PREC index
c
            i1  = int(xpos,IKIND) + 1
            j1  = int(ypos,IKIND) + 1
            k1  = int(zpos,IKIND) + 1
c     
c           Add mass to field
c     
            field(i1, j1, k1) = field(i1,j1,k1) + mass(n)
c
         enddo
      endif
c
      return
#endif
      end

#include "fortran.def"
#include "error.def"
c=======================================================================
c//////////////////////  SUBROUTINE CIC_DEPOSIT  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cic_deposit(posx, posy, posz, ndim, npositions, 
     &                      mass, field, leftedge, 
     &                      dim1, dim2, dim3, cellsize, cloudsize)
#ifndef CONFIG_PFLOAT_16
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
c           The CIC cloud size is set by cloudsize, which should be
c           equal to or less than cellsize.  If cloudsize is equal to
c           cellsize, then this is the typical CIC deposit, but if it
c           is less than cellsize, then the effective cloud size is
c           smaller and the overlap calculations are adjusted.
c
c  INPUTS:
c     ndim       - dimensionality
c     cellsize   - the cell size of field
c     cloudsize  - size of the particle "cloud" (must be <= cellsize)
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
#include "fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC dim1, dim2, dim3, npositions, ndim
      P_PREC posx(npositions), posy(npositions), posz(npositions),
     &        leftedge(3)
      R_PREC    mass(npositions), field(dim1, dim2, dim3), cellsize,
     &          cloudsize
c
c  locals
c
      INTG_PREC iii, jjj, kkk
      INTG_PREC i1, j1, k1, n
      R_PREC    xpos, ypos, zpos, dx, dy, dz
      P_PREC edge1, edge2, edge3, fact, shift, half, refine
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c

!     write(0,*) npositions, leftedge, dim1, dim2, dim3, cellsize

      fact = 1._PKIND/cellsize
      refine = cellsize/cloudsize
      half = 0.5001_PKIND/refine
      shift = 0.5_PKIND/refine
      edge1 = REAL(dim1,PKIND) - half
      edge2 = REAL(dim2,PKIND) - half
      edge3 = REAL(dim3,PKIND) - half
c
      if (cloudsize .gt. cellsize) then
         write(6,*) "cloudsize > cellsize in cic_deposit!"
         write(6,*) cloudsize, cellsize
         ERROR_MESSAGE
      endif
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
            i1  = int(xpos - shift,IKIND) + 1
c
c           Compute the weights
c
            dx = min((REAL(i1,RKIND)+shift-xpos)*refine, 1.0_RKIND)
c
c           Interpolate from field into sumfield
c
            field(i1  ,1,1) = field(i1  ,1,1) + mass(n)*dx
            field(i1+1,1,1) = field(i1+1,1,1) + mass(n)*(1._RKIND-dx)
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
            i1  = int(xpos - shift,IKIND) + 1
            j1  = int(ypos - shift,IKIND) + 1
c
c           Compute the weights
c
            dx = min((REAL(i1,RKIND)+shift-xpos)*refine, 1.0_RKIND)
            dy = min((REAL(j1,RKIND)+shift-ypos)*refine, 1.0_RKIND)
c
c           Interpolate from field into sumfield
c
            field(i1  ,j1  ,1) = field(i1  ,j1  ,1) +
     &                           mass(n)*      dx *      dy
            field(i1+1,j1  ,1) = field(i1+1,j1  ,1) +
     &                           mass(n)*(1._RKIND-dx)*      dy
            field(i1  ,j1+1,1) = field(i1  ,j1+1,1) + 
     &                           mass(n)*      dx *(1._RKIND-dy)
            field(i1+1,j1+1,1) = field(i1+1,j1+1,1) +
     &                           mass(n)*(1._RKIND-dx)*(1._RKIND-dy)
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
c     Compute the position of the central cell
c     
            xpos = min(max((posx(n) - leftedge(1))*fact, half),edge1)
            ypos = min(max((posy(n) - leftedge(2))*fact, half),edge2)
            zpos = min(max((posz(n) - leftedge(3))*fact, half),edge3)
c     
c     Convert this into an INTG_PREC index
c     
            i1  = int(xpos - shift,IKIND) + 1
            j1  = int(ypos - shift,IKIND) + 1
            k1  = int(zpos - shift,IKIND) + 1
c     
c     Compute the weights
c     
            dx = min((REAL(i1,RKIND)+shift-xpos)*refine, 1.0_RKIND)
            dy = min((REAL(j1,RKIND)+shift-ypos)*refine, 1.0_RKIND)
            dz = min((REAL(k1,RKIND)+shift-zpos)*refine, 1.0_RKIND)
c     
c     Interpolate from field into sumfield
c     
            field(i1  ,j1  ,k1  ) = field(i1  ,j1  ,k1  ) +
     &              mass(n)*     dx *     dy *    dz
            field(i1+1,j1  ,k1  ) = field(i1+1,j1  ,k1  ) +
     &              mass(n)*(1._RKIND-dx)*     dy *    dz
            field(i1  ,j1+1,k1  ) = field(i1  ,j1+1,k1  ) + 
     &              mass(n)*     dx *(1._RKIND-dy)*    dz
            field(i1+1,j1+1,k1  ) = field(i1+1,j1+1,k1  ) +
     &              mass(n)*(1._RKIND-dx)*(1._RKIND-dy)*    dz
            field(i1  ,j1  ,k1+1) = field(i1  ,j1  ,k1+1) +
     &              mass(n)*     dx *     dy *(1._RKIND-dz)
            field(i1+1,j1  ,k1+1) = field(i1+1,j1  ,k1+1) +
     &              mass(n)*(1._RKIND-dx)*     dy *(1._RKIND-dz)
            field(i1  ,j1+1,k1+1) = field(i1  ,j1+1,k1+1) + 
     &              mass(n)*     dx *(1._RKIND-dy)*(1._RKIND-dz)
            field(i1+1,j1+1,k1+1) = field(i1+1,j1+1,k1+1) +
     &              mass(n)*(1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)
c     
         enddo
      
      endif
c
      return
#endif
      end

#include "fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE NGPINTERP  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine ngpinterp(parent, dim1, dim2, dim3, ndim,
     &                    start1, start2, start3,
     &                    end1, end2, end3, 
     &                    refine1, refine2, refine3, grid,
     &                    gdim1, gdim2, gdim3, 
     &                    gstart1, gstart2, gstart3)
c
c  PERFORMS A 3D NGP-LIKE INTERPOLATION FROM THE FIELD PARENT TO GRID
c
c     written by: Greg Bryan
c     date:       May, 1995
c     modified1:  
c
c  PURPOSE:  This routine takes the field parent and interpolates it using
c     a second order accurate scheme.
c     NOTE: There is a restriction.  The interpolation must be done in
c        blocks of the parent grid.
c
c  INPUTS:
c     ndim        - rank of fields
c     parent      - parent field
c     dim1,2,3    - declared dimension of parent
c     start1,2,3  - starting index in parent in units of grid (one based)
c     end1,2,3    - ending index in parent in units of grid (one based)
c     refine1,2,3 - INTG_PREC refinement factors
c     gdim1,2,3   - declared dimension of grid
c     gstart1,2,3 - starting offset of refined region in grid (one based)
c
c  OUTPUTS:
c     grid        - grid with refined region
c
c  LOCALS:
c
c  EXTERNALS:
c
c  LOCALS:
c-----------------------------------------------------------------------
      implicit NONE
#include "fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  arguments
c
      INTG_PREC dim1, dim2, dim3, start1, start2, start3, ndim,
     &        end1, end2, end3, refine1, refine2, refine3, 
     &        gdim1, gdim2, gdim3, gstart1, gstart2, gstart3
      R_PREC    parent(dim1, dim2, dim3), grid(gdim1, gdim2, gdim3)
c
c  locals
c
      INTG_PREC i, j, k, iparent, jparent, kparent
      R_PREC    dx, dy, dz, xpos, ypos, zpos
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c  Loop over area to be refined and interpolate
c
c     1d)
c
      if (ndim .eq. 1) then
c
         do i = start1, end1
            iparent = (i-1)/refine1 + 1
            grid(i-start1+gstart1, 1, 1) = parent(iparent  , 1, 1)
         enddo
c
      endif
c
c     2d)
c
      if (ndim .eq. 2) then
c
        do j = start2, end2
          jparent = (j-1)/refine2 + 1
c
          do i = start1, end1
            iparent = (i-1)/refine1 + 1
c
            grid(i-start1+gstart1, j-start2+gstart2, 1) =
     &           parent(iparent  , jparent  , 1)
c
          enddo
        enddo
c
      endif
c
c     3d)
c
      if (ndim .eq. 3) then
      do k = start3, end3
        kparent = (k-1)/refine3 + 1
c
        do j = start2, end2
          jparent = (j-1)/refine2 + 1
c
          do i = start1, end1
            iparent = (i-1)/refine1 + 1
c
            grid(i-start1+gstart1, j-start2+gstart2, k-start3+gstart3) =
     &                         parent(iparent  , jparent  , kparent  )
c
          enddo
        enddo
      enddo
      endif
c
      return
      end

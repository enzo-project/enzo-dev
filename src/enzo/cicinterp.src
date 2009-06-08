#include "fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE CICINTERP  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cicinterp(parent, dim1, dim2, dim3, ndim,
     &                    start1, start2, start3,
     &                    end1, end2, end3, 
     &                    refine1, refine2, refine3, grid,
     &                    gdim1, gdim2, gdim3, 
     &                    gstart1, gstart2, gstart3, ishift)
c
c  PERFORMS A 3D CIC-LIKE INTERPOLATION FROM THE FIELD PARENT TO GRID
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
c     refine1,2,3 - integer refinement factors
c     gdim1,2,3   - declared dimension of grid
c     gstart1,2,3 - starting offset of refined region in grid (one based)
c     ishift      - specified face-centered (0 - none, 1/2/3 x,y,z)
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
c
c-----------------------------------------------------------------------
c
c  arguments
c
      integer dim1, dim2, dim3, start1, start2, start3, ndim,
     &        end1, end2, end3, refine1, refine2, refine3, 
     &        gdim1, gdim2, gdim3, gstart1, gstart2, gstart3, ishift
      real    parent(dim1, dim2, dim3), grid(gdim1, gdim2, gdim3)
c
c  locals
c
      integer i, j, k, iparent, jparent, kparent
      real    dx, dy, dz, xpos, ypos, zpos, xshift, yshift, zshift
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
      xshift = 0.5
      yshift = 0.5
      zshift = 0.5
      if (ishift .eq. 1) xshift = 1.0
      if (ishift .eq. 2) yshift = 1.0
      if (ishift .eq. 3) zshift = 1.0
c
c  Loop over area to be refined and interpolate
c
c     1d)
c
      if (ndim .eq. 1) then
c
         do i = start1, end1
            xpos    = (real(i)-xshift)/real(refine1) + 1.0
            iparent = int(xpos-(1.0-xshift))
            dx      = real(iparent) - xpos + 2.0 - xshift
c
            grid(i-start1+gstart1, 1, 1) =
     &      dx  * parent(iparent  , 1, 1) +
     &   (1-dx) * parent(iparent+1, 1, 1)
         enddo
c
      endif
c
c     2d)
c
      if (ndim .eq. 2) then
c
        do j = start2, end2
          ypos    = (real(j)-yshift)/real(refine2) + 1.0
          jparent = int(ypos-(1.0-yshift))
          dy      = real(jparent) - ypos + 2.0 - yshift
c
          do i = start1, end1
            xpos    = (real(i)-xshift)/real(refine1) + 1.0
            iparent = int(xpos-(1.0-xshift))
            dx      = real(iparent) - xpos + 2.0 - xshift
c
            grid(i-start1+gstart1, j-start2+gstart2, 1) =
     &      dx *   dy * parent(iparent  , jparent  , 1) +
     &   (1-dx)*   dy * parent(iparent+1, jparent  , 1) +
     &      dx *(1-dy)* parent(iparent  , jparent+1, 1) +
     &   (1-dx)*(1-dy)* parent(iparent+1, jparent+1, 1)
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
        zpos    = (real(k)-zshift)/real(refine3) + 1.0
        kparent = int(zpos-(1.0-zshift))
        dz      = real(kparent) - zpos + 2.0 - zshift
c
        do j = start2, end2
          ypos    = (real(j)-yshift)/real(refine2) + 1.0
          jparent = int(ypos-(1.0-yshift))
          dy      = real(jparent) - ypos + 2.0 - yshift
c
          do i = start1, end1
            xpos    = (real(i)-xshift)/real(refine1) + 1.0
            iparent = int(xpos-(1.0-xshift))
            dx      = real(iparent) - xpos + 2.0 - xshift
c
            grid(i-start1+gstart1, j-start2+gstart2, k-start3+gstart3) =
     &      dx *   dy *   dz * parent(iparent  , jparent  , kparent  ) +
     &   (1-dx)*   dy *   dz * parent(iparent+1, jparent  , kparent  ) +
     &      dx *(1-dy)*   dz * parent(iparent  , jparent+1, kparent  ) +
     &   (1-dx)*(1-dy)*   dz * parent(iparent+1, jparent+1, kparent  ) +
     &      dx *   dy *(1-dz)* parent(iparent  , jparent  , kparent+1) +
     &   (1-dx)*   dy *(1-dz)* parent(iparent+1, jparent  , kparent+1) +
     &      dx *(1-dy)*(1-dz)* parent(iparent  , jparent+1, kparent+1) +
     &   (1-dx)*(1-dy)*(1-dz)* parent(iparent+1, jparent+1, kparent+1)
c
          enddo
        enddo
      enddo
      endif
c
      return
      end

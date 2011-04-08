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
c     refine1,2,3 - INTG_PREC refinement factors
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
#include "fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  arguments
c
      INTG_PREC dim1, dim2, dim3, start1, start2, start3, ndim,
     &        end1, end2, end3, refine1, refine2, refine3, 
     &        gdim1, gdim2, gdim3, gstart1, gstart2, gstart3, ishift
      R_PREC    parent(dim1, dim2, dim3), grid(gdim1, gdim2, gdim3)
c
c  locals
c
      INTG_PREC i, j, k, iparent, jparent, kparent
      R_PREC    dx, dy, dz, xpos, ypos, zpos, xshift, yshift, zshift
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
      xshift = 0.5_RKIND
      yshift = 0.5_RKIND
      zshift = 0.5_RKIND
      if (ishift .eq. 1) xshift = 1._RKIND
      if (ishift .eq. 2) yshift = 1._RKIND
      if (ishift .eq. 3) zshift = 1._RKIND
c
c  Loop over area to be refined and interpolate
c
c     1d)
c
      if (ndim .eq. 1) then
c
         do i = start1, end1
            xpos    = (REAL(i,RKIND)-xshift)/REAL(refine1,RKIND) 
     &           + 1._RKIND
            iparent = int(xpos-(1._RKIND-xshift),IKIND)
            dx      = REAL(iparent,RKIND) - xpos + 2._RKIND - xshift
c
            grid(i-start1+gstart1, 1, 1) =
     &                 dx  * parent(iparent  , 1, 1) +
     &           (1._RKIND-dx) * parent(iparent+1, 1, 1)
         enddo
c
      endif
c
c     2d)
c
      if (ndim .eq. 2) then
c
        do j = start2, end2
          ypos    = (REAL(j,RKIND)-yshift)/REAL(refine2,RKIND) 
     &          + 1._RKIND
          jparent = int(ypos-(1._RKIND-yshift),IKIND)
          dy      = REAL(jparent,RKIND) - ypos + 2._RKIND - yshift
c
          do i = start1, end1
            xpos = (REAL(i,RKIND)-xshift)/REAL(refine1,RKIND) 
     &            + 1._RKIND
            iparent = int(xpos-(1._RKIND-xshift),IKIND)
            dx      = REAL(iparent,RKIND) - xpos + 2._RKIND - xshift
c
            grid(i-start1+gstart1, j-start2+gstart2, 1) =
     &         dx *      dy * parent(iparent  , jparent  , 1) +
     &   (1._RKIND-dx)*      dy * parent(iparent+1, jparent  , 1) +
     &         dx *(1._RKIND-dy)* parent(iparent  , jparent+1, 1) +
     &   (1._RKIND-dx)*(1._RKIND-dy)* parent(iparent+1, jparent+1, 1)
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
        zpos    = (REAL(k,RKIND)-zshift)/REAL(refine3,RKIND) + 1._RKIND
        kparent = int(zpos-(1._RKIND-zshift),IKIND)
        dz      = REAL(kparent,RKIND) - zpos + 2._RKIND - zshift
c
        do j = start2, end2
          ypos    = (REAL(j,RKIND)-yshift)/REAL(refine2,RKIND) 
     &          + 1._RKIND
          jparent = int(ypos-(1._RKIND-yshift),IKIND)
          dy      = REAL(jparent,RKIND) - ypos + 2._RKIND - yshift
c
          do i = start1, end1
            xpos    = (REAL(i,RKIND)-xshift)/REAL(refine1,RKIND) 
     &            + 1._RKIND
            iparent = int(xpos-(1._RKIND-xshift),IKIND)
            dx      = REAL(iparent,RKIND) - xpos + 2._RKIND - xshift
c
            grid(i-start1+gstart1, j-start2+gstart2, k-start3+gstart3) =
     &           dx *   dy *   dz * 
     &           parent(iparent  , jparent  , kparent  ) +
     &           (1._RKIND-dx)*   dy *   dz * 
     &           parent(iparent+1, jparent  , kparent  ) +
     &           dx *(1._RKIND-dy)*   dz * 
     &           parent(iparent  , jparent+1, kparent  ) +
     &           (1._RKIND-dx)*(1._RKIND-dy)*   dz * 
     &           parent(iparent+1, jparent+1, kparent  ) +
     &           dx *   dy *(1._RKIND-dz)* 
     &           parent(iparent  , jparent  , kparent+1) +
     &           (1._RKIND-dx)*   dy *(1._RKIND-dz)* 
     &           parent(iparent+1, jparent  , kparent+1) +
     &           dx *(1._RKIND-dy)*(1._RKIND-dz)* 
     &           parent(iparent  , jparent+1, kparent+1) +
     &           (1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)* 
     &           parent(iparent+1, jparent+1, kparent+1)
c
          enddo
        enddo
      enddo
      endif
c
      return
      end

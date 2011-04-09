#include "fortran.def"
#include "error.def"

c=======================================================================
c///////////////////////  SUBROUTINE TSCINT1D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine tscint1d(parent, dim1, start1, end1, refine1,
     &                    grid, gdim1, gstart1)
c
c  PERFORMS A 1D TSC-LIKE INTERPOLATION FROM THE FIELD PARENT TO GRID
c
c     written by: Greg Bryan
c     date:       November, 1994
c     modified1:  January, 1995 to make the scheme monotonic (GB)
c
c  PURPOSE:  This routine takes the field parent and interpolates it using
c     a varient of the third-order triangular-shaped-cloud interpolation
c     scheme.
c     NOTE: There is a restriction.  The interpolation must be done in by
c        blocks of the parent grid.
c
c  INPUTS:
c     parent    - parent field
c     dim1      - declared dimension of parent
c     start1    - starting index in parent in units of grid (one based)
c     end1      - ending index in parent in units of grid (one based)
c     refine1   - INTG_PREC refinement factors
c     gdim1     - declared dimension of grid
c     gstart1   - starting offset of refined region in grid (one based)
c
c  OUTPUTS:
c     grid      - grid with refined region
c
c  EXTERNALS:
c
c  LOCALS:
c     fraction  - a work array of size gdim1/refine1
c-----------------------------------------------------------------------
      implicit NONE
#include "fortran_types.def"
c
      INTG_PREC ijkn, MAX_REFINE
      parameter (ijkn = MAX_ANY_SINGLE_DIRECTION)
      parameter (MAX_REFINE = 32)
c
c-----------------------------------------------------------------------
c
c  arguments
c
      INTG_PREC dim1, start1, end1, refine1, gdim1, gstart1
      R_PREC    parent(dim1), grid(gdim1)
c
c  locals
c
      INTG_PREC i, j, k
      R_PREC  coefm(MAX_REFINE), coefp(MAX_REFINE), coef0(MAX_REFINE),
     &        del1left, del2left, del1right, del2right, factor,
     &        fraction(ijkn), sum
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c  error check
c
      if (refine1 .gt. MAX_REFINE) then
         write (6,*) "Error: refine1 > MAX_REFINE"
         ERROR_MESSAGE
      endif
c
c  Compute constants
c
      do i = 0, refine1-1
         coefm(i+1) = REAL((refine1-i)**3 - (refine1-i-1)**3,RKIND)/
     &                REAL(6*refine1**3                     ,RKIND)
         coefp(i+1) = REAL(3*i**2 + 3*i + 1,RKIND)                 /
     &                REAL(6*refine1**3                     ,RKIND)
         coef0(i+1) = 1._RKIND/refine1 - coefm(i+1) - coefp(i+1)
      enddo
c
c  Loop over area to be refined and interpolate
c
      do i = start1, end1
         j = (i-1)/refine1 + 1
         k = i - (j-1)*refine1
         grid(i-start1+gstart1) = coefm(k)*parent(j-1) + 
     &                            coefp(k)*parent(j+1) + 
     &                            coef0(k)*parent(j  )
      enddo
c
c  Correct to make sure the interpolation is conservatiove
c
      do i = start1, end1, refine1
         sum = 0._RKIND
         do j = 0, refine1-1
            sum = sum + grid(i-start1+gstart1+j)
         enddo
         factor = (parent((i-1)/refine1+1) - sum)/REAL(refine1,RKIND)
         do j = 0, refine1-1
            grid(i-start1+gstart1+j) = grid(i-start1+gstart1+j) +
     &               factor
c     &               parent((i-1)/refine1+1)/(sum + tiny)

         enddo
      enddo
c
c  Correct to ensure monotonicity.  This is done by comparing the ends
c    to their neighbour and if they are local minima or maxima (relative
c    to the neighbour and the average value), then decrease the entire
c    profile as needed, but do it in such a way as to keep the average value
c    the same.
c
      do i = start1, end1, refine1
         if (i .eq. start1) then
            del1left = parent((i-1)/refine1)/REAL(refine1,RKIND) -
     &                 grid(i-start1+gstart1)
         else
            del1left = grid(i-start1+gstart1-1) - grid(i-start1+gstart1)
         endif
         del2left = parent((i-1)/refine1+1)/REAL(refine1,RKIND) - 
     &              grid(i-start1+gstart1)
         if (abs(del2left) .lt. tiny) 
     &        del2left = sign(REAL(tiny,RKIND), del2left)
         fraction((i-start1)/refine1 + 1) = 
     &             min(max(del1left/del2left, 0._RKIND), 1._RKIND)
      enddo
c
c       repeat for right side of each distribution
c
      do i = start1, end1, refine1
         if (i .eq. end1-refine1+1) then
            del1right = parent((i-1)/refine1+2)/REAL(refine1,RKIND) -
     &                 grid(i-start1+gstart1+refine1-1)
         else
            del1right = grid(i-start1+gstart1+refine1) - 
     &                 grid(i-start1+gstart1+refine1-1)
         endif
         del2right = parent((i-1)/refine1+1)/REAL(refine1,RKIND) - 
     &              grid(i-start1+gstart1+refine1-1)
         if (abs(del2right) .lt. tiny) del2right = 
     &        sign(REAL(tiny,RKIND), del2right)
         fraction((i-start1)/refine1 + 1) = max(
     &                    min(max(del1right/del2right, 0._RKIND), 
     &                    1._RKIND), fraction((i-start1)/refine1 + 1) )
      enddo
c
c   Modify distributions based on the computed fraction decrease
c
      do i = start1, end1
            grid(i-start1+gstart1) = 
     &           parent((i-1)/refine1+1)/REAL(refine1,RKIND)*
     &                        fraction((i-start1)/refine1 + 1) +
     &           grid(i-start1+gstart1)             *
     &                 (1._RKIND - fraction((i-start1)/refine1 + 1) )
      enddo
c
c     Remove the irritating division by the refinement volume
c
      do i = gstart1, gstart1 + end1-start1
         grid(i) = grid(i) * refine1
      enddo
c
      return
      end

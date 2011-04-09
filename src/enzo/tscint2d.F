#include "fortran.def"
#include "error.def"

c=======================================================================
c///////////////////////  SUBROUTINE TSCINT2D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine tscint2d(parent, dim1, dim2, start1, start2, 
     &                    end1, end2, refine1, refine2, grid,
     &                    gdim1, gdim2, gstart1, gstart2)
c
c  PERFORMS A 2D TSC-LIKE INTERPOLATION FROM THE FIELD PARENT TO GRID
c
c     written by: Greg Bryan
c     date:       February, 1995
c     modified1:  
c
c  PURPOSE:  This routine takes the field parent and interpolates it using
c     a varient of the third-order triangular-shaped-cloud interpolation
c     scheme.
c     NOTE: There is a restriction.  The interpolation must be done in
c        blocks of the parent grid.
c     We also assume refine1 = refine2.  This is relatively easy to left,
c     but so far has not been required.
c
c  INPUTS:
c     parent    - parent field
c     dim1,2    - declared dimension of parent
c     start1,2  - starting index in parent in units of grid (one based)
c     end1,2    - ending index in parent in units of grid (one based)
c     refine1,2 - INTG_PREC refinement factors
c     gdim1,2   - declared dimension of grid
c     gstart1,2 - starting offset of refined region in grid (one based)
c
c  OUTPUTS:
c     grid      - grid with refined region
c
c  LOCALS:
c     fraction  - a work array of size gdim1/refine1
c
c  EXTERNALS:
c
c  LOCALS:
c-----------------------------------------------------------------------
      implicit NONE
#include "fortran_types.def"
c
      INTG_PREC ijkn, MAX_REFINE, num_iter
      parameter (ijkn = MAX_ANY_SINGLE_DIRECTION, num_iter = 2)
      parameter (MAX_REFINE = 32)
c
c-----------------------------------------------------------------------
c
c  arguments
c
      INTG_PREC dim1, dim2, start1, start2, end1, end2, 
     &        refine1, refine2, gdim1, gdim2, gstart1, gstart2
      R_PREC    parent(dim1, dim2), grid(gdim1, gdim2)
c
c  locals
c
      INTG_PREC i, icoef, igrid, iparent, isub, iter,
     &        j, jcoef, jgrid, jparent, jsub, k
      R_PREC coefm1(MAX_REFINE), coefp1(MAX_REFINE), coef01(MAX_REFINE),
     &       coefm2(MAX_REFINE), coefp2(MAX_REFINE), coef02(MAX_REFINE),
     &       del1left, del2left, del1right, del2right, dV, factor, sum
      R_PREC   temp(ijkn), fraction(ijkn)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c  error check
c
      if (refine1 .gt. MAX_REFINE .or. refine2 .gt. MAX_REFINE) then
         write (6,*) "tscint2d: Error: refine1 or 2 > MAX_REFINE"
         ERROR_MESSAGE
      endif
c
c  Compute constants
c
      do i = 0, refine1-1
         coefm1(i+1) = REAL((refine1-i)**3 - (refine1-i-1)**3,RKIND)/
     &                 REAL(6*refine1**3                     ,RKIND)
         coefp1(i+1) = REAL(3*i**2 + 3*i + 1,RKIND)                 /
     &                 REAL(6*refine1**3                     ,RKIND)
         coef01(i+1) = 1._RKIND/refine1 - coefm1(i+1) - coefp1(i+1)
      enddo
c
      do i = 0, refine2-1
         coefm2(i+1) = REAL((refine2-i)**3 - (refine2-i-1)**3,RKIND)/
     &                 REAL(6*refine2**3                     ,RKIND)
         coefp2(i+1) = REAL(3*i**2 + 3*i + 1,RKIND)                 /
     &                 REAL(6*refine2**3                     ,RKIND)
         coef02(i+1) = 1._RKIND/refine2 - coefm2(i+1) - coefp2(i+1)
      enddo
c
c  Loop over area to be refined and interpolate
c
      do j = start2, end2
         jparent = (j-1)/refine2 + 1
         jcoef   = j - (jparent-1)*refine2
         do i = start1, end1
            iparent = (i-1)/refine1 + 1
            icoef   = i - (iparent-1)*refine1
            grid(i-start1+gstart1, j-start2+gstart2) =
     &       coefm1(icoef)*coefm2(jcoef)*parent(iparent-1, jparent-1) +        
     &       coefm1(icoef)*coef02(jcoef)*parent(iparent-1, jparent  ) +        
     &       coefm1(icoef)*coefp2(jcoef)*parent(iparent-1, jparent+1) +        
     &       coef01(icoef)*coefm2(jcoef)*parent(iparent  , jparent-1) +        
     &       coef01(icoef)*coef02(jcoef)*parent(iparent  , jparent  ) +        
     &       coef01(icoef)*coefp2(jcoef)*parent(iparent  , jparent+1) +        
     &       coefp1(icoef)*coefm2(jcoef)*parent(iparent+1, jparent-1) +
     &       coefp1(icoef)*coef02(jcoef)*parent(iparent+1, jparent  ) +
     &       coefp1(icoef)*coefp2(jcoef)*parent(iparent+1, jparent+1)
         enddo
      enddo
c
c  Correct to make sure the interpolation is conservative
c
      do j = start2, end2, refine2
         do i = start1, end1, refine1
            sum = 0._RKIND
            do jsub = 0, refine2-1
               do isub = 0, refine1-1
                  sum = sum + grid(i-start1+gstart1+isub, 
     &                             j-start2+gstart2+jsub )
               enddo
            enddo
            factor = (parent((i-1)/refine1+1, (j-1)/refine2+1) - sum)/
     &                REAL(refine1*refine2,RKIND)
            do jsub = 0, refine2-1
               do isub = 0, refine1-1
c                  grid(i-start1+gstart1+isub, j-start2+gstart2+jsub) = 
c     &             grid(i-start1+gstart1+isub, j-start2+gstart2+jsub)*
c     &             parent((i-1)/refine1+1, (j-1)/refine2+1)/
c     &             (sum + tiny)
                  grid(i-start1+gstart1+isub, j-start2+gstart2+jsub) = 
     &             grid(i-start1+gstart1+isub, j-start2+gstart2+jsub) +
     &             factor
               enddo
            enddo
         enddo
      enddo
c
c  Correct to ensure monotonicity.  This is done by comparing the ends
c    to their neighbour and if they are local minima or maxima (relative
c    to the neighbour and the average value), then decrease the entire
c    profile as needed, but do it in such a way as to keep the average value
c    the same.
c
      dV = 1.0/REAL(refine1*refine2,RKIND)
      do iter = 1, num_iter
c
c  1) Do it first in the i direction
c
      do j = start2, end2
         jparent = (j-1)/refine2+1
         jgrid   = j-start2+gstart2
c
c        Compute the average along a single line of the interpolated values
c           a) zero temp space
c           b) average over one line of the interpolated grid in i direction
c
         do i = 1, (end1-start1+1)/refine1
            temp(i) = 0._RKIND
         enddo
         do i = start1, end1
            temp((i-start1)/refine1 + 1) = temp((i-start1)/refine1 + 1)
     &         + grid(i-start1+gstart1, jgrid)/REAL(refine1,RKIND)
         enddo
c
c        Check to see if the end of the distribution is a local maxima
c         (compared to the nearest interpolated grid point to the left
c          and the mean value along this interpolated line -- computed above).
c
         do i = start1, end1, refine1
            if (i .eq. start1) then
               del1left = parent((i-1)/refine1, jparent)*dV -
     &                    grid(i-start1+gstart1  , jgrid)
            else
               del1left = grid(i-start1+gstart1-1, jgrid) - 
     &                    grid(i-start1+gstart1  , jgrid)
            endif
            del2left = temp((i-start1)/refine1+1) -
     &                 grid(i-start1+gstart1, jgrid)
            if (abs(del2left) .lt. tiny) 
     &           del2left = sign(REAL(tiny,RKIND), del2left)
            fraction((i-start1)/refine1 + 1) = 
     &           min(max(del1left/del2left, 0._RKIND), 1._RKIND)
         enddo
c
c        repeat for right side of each distribution
c
         do i = start1, end1, refine1
            if (i .eq. end1-refine1+1) then
               del1right = parent((i-1)/refine1+2, jparent)*dV -
     &                     grid(i-start1+gstart1+refine1-1, jgrid)
            else
               del1right = grid(i-start1+gstart1+refine1  , jgrid) - 
     &                     grid(i-start1+gstart1+refine1-1, jgrid)
            endif
            del2right = temp((i-start1)/refine1+1) -
     &           grid(i-start1+gstart1+refine1-1, jgrid)
            if (abs(del2right) .lt. tiny) del2right = 
     &           sign(REAL(tiny,RKIND), del2right)
            fraction((i-start1)/refine1 + 1) = max(
     &                    min(max(del1right/del2right, 0._RKIND), 
     &                    1._RKIND), fraction((i-start1)/refine1 + 1) )
         enddo
c
c        Modify distributions based on the computed fraction decrease
c
         do i = start1, end1
            grid(i-start1+gstart1, jgrid) = 
     &         temp((i-start1)/refine1+1) 
     &           *        fraction((i-start1)/refine1+1) +
     &         grid(i-start1+gstart1, jgrid)
     &           * (1._RKIND - fraction((i-start1)/refine1+1) )
         enddo
      enddo
c
c  2) Repeat for the j-direction
c
      do i = start1, end1
         iparent = (i-1)/refine1+1
         igrid   = i-start1+gstart1
c
c        Compute the average along a single line of the interpolated values
c
         do j = 1, (end2-start2+1)/refine2
            temp(j) = 0._RKIND
         enddo
         do j = start2, end2
            temp((j-start2)/refine2 + 1) = temp((j-start2)/refine2 + 1)
     &         + grid(igrid, j-start2+gstart2)/REAL(refine2,RKIND)
         enddo
c
c        Check to see if the end of the distribution is a local maxima
c         (compared to the nearest interpolated grid point to the left
c          and the mean value along this interpolated line -- computed above).
c
         do j = start2, end2, refine2
            if (j .eq. start2) then
               del1left = parent(iparent, (j-1)/refine2)*dV -
     &                    grid(igrid, j-start2+gstart2  )
            else
               del1left = grid(igrid, j-start2+gstart2-1) - 
     &                    grid(igrid, j-start2+gstart2  )
            endif
            del2left = temp((j-start2)/refine2+1) -
     &                 grid(igrid, j-start2+gstart2)
            if (abs(del2left) .lt. tiny) 
     &           del2left = sign(REAL(tiny,RKIND), del2left)
            fraction((j-start2)/refine2 + 1) = 
     &           min(max(del1left/del2left, 0._RKIND), 1._RKIND)
         enddo
c
c        repeat for right side of each distribution
c
         do j = start2, end2, refine2
            if (j .eq. end2-refine2+1) then
               del1right = parent(iparent, (j-1)/refine2+2)*dV -
     &                     grid(igrid, j-start2+gstart2+refine2-1)
            else
               del1right = grid(igrid, j-start2+gstart2+refine2  ) - 
     &                     grid(igrid, j-start2+gstart2+refine2-1)
            endif
            del2right = temp((j-start2)/refine2+1) -
     &           grid(igrid, j-start2+gstart2+refine2-1)
            if (abs(del2right) .lt. tiny) del2right = 
     &           sign(REAL(tiny,RKIND), del2right)
            fraction((j-start2)/refine2 + 1) = max(
     &                    min(max(del1right/del2right, 0._RKIND), 
     &                    1._RKIND), fraction((j-start2)/refine2 + 1) )
         enddo
c
c        Modify distributions based on the computed fraction decrease
c
         do j = start2, end2
            grid(igrid, j-start2+gstart2) = 
     &         temp((j-start2)/refine2+1) 
     &           *        fraction((j-start2)/refine2+1) +
     &         grid(igrid, j-start2+gstart2)
     &           * (1._RKIND - fraction((j-start2)/refine2+1))
         enddo
      enddo
c
c     Next iteration
c
      enddo
c
c     Remove the irritating division by the refinement volume
c
      do j = gstart2, gstart2 + end2-start2
         do i = gstart1, gstart1 + end1-start1
            grid(i,j) = grid(i,j) * refine1*refine2
         enddo
      enddo
c
      return
      end

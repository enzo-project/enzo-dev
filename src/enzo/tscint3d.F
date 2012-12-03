#include "fortran.def"
#include "error.def"

c=======================================================================
c///////////////////////  SUBROUTINE TSCINT2D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine tscint3d(parent, dim1, dim2, dim3, 
     &                    start1, start2, start3,
     &                    end1, end2, end3, 
     &                    refine1, refine2, refine3, grid,
     &                    gdim1, gdim2, gdim3, 
     &                    gstart1, gstart2, gstart3)
c
c  PERFORMS A 3D TSC-LIKE INTERPOLATION FROM THE FIELD PARENT TO GRID
c
c     written by: Greg Bryan
c     date:       May, 1995
c     modified1:  
c
c  PURPOSE:  This routine takes the field parent and interpolates it using
c     a varient of the third-order triangular-shaped-cloud interpolation
c     scheme.
c     NOTE: There is a restriction.  The interpolation must be done in
c        blocks of the parent grid.
c
c  INPUTS:
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
c     temp        - temporary work field of size max(gdim#/refine#)
c     fraction    - a work array of size max(gdim#/refine#)
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
      INTG_PREC dim1, dim2, dim3, start1, start2, start3, 
     &        end1, end2, end3, refine1, refine2, refine3, 
     &        gdim1, gdim2, gdim3, gstart1, gstart2, gstart3
      R_PREC    parent(dim1, dim2, dim3), grid(gdim1, gdim2, gdim3)
c
c  locals
c
      INTG_PREC i, i1, i2, icoef, igrid, iparent, isub, iter,
     &        j, j1, j2, jcoef, jgrid, jparent, jsub, 
     &        k, k1, k2, kcoef, kgrid, kparent, ksub, 
     &        wdim1, wdim2, wdim3
      R_PREC coefm1(MAX_REFINE), coefp1(MAX_REFINE), coef01(MAX_REFINE),
     &       coefm2(MAX_REFINE), coefp2(MAX_REFINE), coef02(MAX_REFINE),
     &       coefm3(MAX_REFINE), coefp3(MAX_REFINE), coef03(MAX_REFINE),
     &       del1left, del2left, del1next, del2next, 
     &       del1right, del2right, dV, factor, sum,
     &       w01, w02, w03, wm1, wm2, wm3, wp1, wp2, wp3
      R_PREC   temp(ijkn), fraction(ijkn)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c  error check
c
      if (refine1 .gt. MAX_REFINE .or. refine2 .gt. MAX_REFINE .or.
     &    refine2 .gt. MAX_REFINE) then
         write (6,*) "tscint3d: Error: refine1,2 or 3 > MAX_REFINE"
         ERROR_MESSAGE
      endif
c
c  Set work dimensions
c
      wdim1 = (end1-start1+1)/refine1
      wdim2 = (end2-start2+1)/refine2
      wdim3 = (end3-start3+1)/refine3
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
      do i = 0, refine3-1
         coefm3(i+1) = REAL((refine3-i)**3 - (refine3-i-1)**3,RKIND)/
     &                 REAL(6*refine3**3                     ,RKIND)
         coefp3(i+1) = REAL(3*i**2 + 3*i + 1,RKIND)                 /
     &                 REAL(6*refine3**3                     ,RKIND)
         coef03(i+1) = 1._RKIND/refine3 - coefm3(i+1) - coefp3(i+1)
      enddo
c
c  Loop over area to be refined and interpolate
c
      do k = start3, end3
        kparent = (k-1)/refine3 + 1
        kcoef   = k - (kparent-1)*refine3
        wm3     = coefm3(kcoef)
        w03     = coef03(kcoef)
        wp3     = coefp3(kcoef)
c
        do j = start2, end2
          jparent = (j-1)/refine2 + 1
          jcoef   = j - (jparent-1)*refine2
          wm2     = coefm2(jcoef)
          w02     = coef02(jcoef)
          wp2     = coefp2(jcoef)
c
          do i = start1, end1
            iparent = (i-1)/refine1 + 1
            icoef   = i - (iparent-1)*refine1
            wm1     = coefm1(icoef)
            w01     = coef01(icoef)
            wp1     = coefp1(icoef)
c
            grid(i-start1+gstart1, j-start2+gstart2, k-start3+gstart3) =
     &       wm1*wm2*wm3* parent(iparent-1, jparent-1, kparent-1) +        
     &       w01*wm2*wm3* parent(iparent  , jparent-1, kparent-1) +        
     &       wp1*wm2*wm3* parent(iparent+1, jparent-1, kparent-1) +
     &       wm1*w02*wm3* parent(iparent-1, jparent  , kparent-1) +        
     &       w01*w02*wm3* parent(iparent  , jparent  , kparent-1) +        
     &       wp1*w02*wm3* parent(iparent+1, jparent  , kparent-1) +
     &       wm1*wp2*wm3* parent(iparent-1, jparent+1, kparent-1) +        
     &       w01*wp2*wm3* parent(iparent  , jparent+1, kparent-1) +        
     &       wp1*wp2*wm3* parent(iparent+1, jparent+1, kparent-1) +

     &       wm1*wm2*w03* parent(iparent-1, jparent-1, kparent  ) +        
     &       w01*wm2*w03* parent(iparent  , jparent-1, kparent  ) +        
     &       wp1*wm2*w03* parent(iparent+1, jparent-1, kparent  ) +
     &       wm1*w02*w03* parent(iparent-1, jparent  , kparent  ) +        
     &       w01*w02*w03* parent(iparent  , jparent  , kparent  ) +        
     &       wp1*w02*w03* parent(iparent+1, jparent  , kparent  ) +
     &       wm1*wp2*w03* parent(iparent-1, jparent+1, kparent  ) +        
     &       w01*wp2*w03* parent(iparent  , jparent+1, kparent  ) +        
     &       wp1*wp2*w03* parent(iparent+1, jparent+1, kparent  ) +

     &       wm1*wm2*wp3* parent(iparent-1, jparent-1, kparent+1) +        
     &       w01*wm2*wp3* parent(iparent  , jparent-1, kparent+1) +        
     &       wp1*wm2*wp3* parent(iparent+1, jparent-1, kparent+1) +
     &       wm1*w02*wp3* parent(iparent-1, jparent  , kparent+1) +        
     &       w01*w02*wp3* parent(iparent  , jparent  , kparent+1) +        
     &       wp1*w02*wp3* parent(iparent+1, jparent  , kparent+1) +
     &       wm1*wp2*wp3* parent(iparent-1, jparent+1, kparent+1) +        
     &       w01*wp2*wp3* parent(iparent  , jparent+1, kparent+1) +        
     &       wp1*wp2*wp3* parent(iparent+1, jparent+1, kparent+1)
c
          enddo
        enddo
      enddo
c
c  Correct to make sure the interpolation is conservative.
c    We do this by compute the interpolated sum and then adding the difference
c    between that value and the top grid value.
c
      do k = 1, wdim3
         do j = 1, wdim2
c     
c           Clear the work field
c
            do i = 1, wdim1
               temp(i) = 0._RKIND
            enddo
c
c           Compute the sum over each parent zone
c
            do ksub = gstart3, gstart3+refine3-1
               do jsub = gstart2, gstart2+refine2-1
                  do isub = gstart1, gstart1+refine1-1
                     do i = 1, wdim1
                        temp(i) = temp(i) + grid((i-1)*refine1 + isub,
     &                                           (j-1)*refine2 + jsub,
     &                                           (k-1)*refine3 + ksub)
                     enddo
                  enddo
               enddo
            enddo
c
c           Compute difference between subgrid sum and parent zone value
c
            do i = 1, wdim1
               temp(i) = (parent((start1-1)/refine1+i, 
     &                           (start2-1)/refine2+j, 
     &                           (start3-1)/refine3+k) - temp(i))/
     &                   REAL(refine1*refine2*refine3)
            enddo
c
c           Add difference to each subgrid value
c
            do ksub = gstart3, gstart3+refine3-1
               do jsub = gstart2, gstart2+refine2-1
                  do isub = gstart1, gstart1+refine1-1
                     do i = 1, wdim1
                        grid((i-1)*refine1 + isub, 
     &                       (j-1)*refine2 + jsub,
     &                       (k-1)*refine3 + ksub) = 
     &                  grid((i-1)*refine1 + isub, 
     &                       (j-1)*refine2 + jsub,
     &                       (k-1)*refine3 + ksub) + temp(i)
                     enddo
                  enddo
               enddo
            enddo
c
c        Next j,k
c
         enddo
      enddo
c
c  Correct to ensure monotonicity.  This is done by comparing the ends
c    to their neighbour and if they are local minima or maxima (relative
c    to the neighbour and the average value), then decrease the entire
c    profile as needed, but do it in such a way as to keep the average value
c    the same.
c
      dV = 1._RKIND/REAL(refine1*refine2*refine3)
      do iter = 1, num_iter
c
c  1) Do it first in the i direction
c
      do k = start3, end3
       kparent = (k-1)/refine3+1
       kgrid   = k-start3+gstart3
       do j = start2, end2
         jparent = (j-1)/refine2+1
         jgrid   = j-start2+gstart2
c
c        Compute the average along a single line of the interpolated values
c           a) zero temp space
c           b) average over one line of the interpolated grid in i direction
c
         do i = 1, wdim1
            temp(i) = 0._RKIND
         enddo
c
         do i1 = gstart1, gstart1+refine1-1
            do i2 = 1, wdim1
               temp(i2) = temp(i2) + 
     &              grid((i2-1)*refine1+i1, jgrid, kgrid)/REAL(refine1)
            enddo
         enddo
c
c        Check to see if the end of the distribution is a local maxima
c         (compared to the nearest interpolated grid point to the left
c          and the mean value along this interpolated line -- computed above).
c
         del1next = parent((start1-1)/refine1, jparent, kparent)*dV
         do i = 1, wdim1
            del1left = del1next - 
     &                 grid((i-1)*refine1+gstart1  , jgrid, kgrid)
            del1next = grid( i   *refine1+gstart1-1, jgrid, kgrid)
            del2left = temp(i) - 
     &                 grid((i-1)*refine1+gstart1  , jgrid, kgrid)
c
            if (abs(del2left) .lt. tiny) 
     &           del2left = sign(REAL(tiny,RKIND), del2left)
            fraction(i) = min(max(del1left/del2left, 0._RKIND),1._RKIND)
         enddo
c
c        repeat for right side of each distribution
c
         del1next = parent(end1/refine1+1, jparent, kparent)*dV
         do i = wdim1, 1, -1
            del1right = del1next -
     &                  grid( i   *refine1+gstart1-1, jgrid, kgrid)
            del1next  = grid((i-1)*refine1+gstart1  , jgrid, kgrid)
            del2right = temp(i) -
     &                  grid( i   *refine1+gstart1-1, jgrid, kgrid)
c
            if (abs(del2right) .lt. tiny) 
     &           del2right = sign(REAL(tiny,RKIND), del2right)
            fraction(i) =  max(max(min(del1right/del2right, 0._RKIND),
     &           1._RKIND), fraction(i)                             )
         enddo
c
c        Modify distributions based on the computed fraction decrease
c
c         do i = start1, end1
c            grid(i-start1+gstart1, jgrid, kgrid) = 
c     &         temp((i-start1)/refine1+1) 
c     &                        *        fraction((i-start1)/refine1+1) +
c     &         grid(i-start1+gstart1, jgrid, kgrid)
c     &                        * (1._RKIND - fraction((i-start1)/refine1+1))
c         enddo
c
         do i1 = gstart1, gstart1+refine1-1
            do i2 = 1, wdim1
               grid((i2-1)*refine1+i1, jgrid, kgrid) = 
     &                     fraction(i2)  * temp(i2) +
     &              (1._RKIND - fraction(i2)) * grid((i2-1)*refine1+i1, 
     &                                          jgrid, kgrid)
            enddo
         enddo
c
       enddo
      enddo
c
c  2) Repeat for the j-direction
c
      do k = start3, end3
       kparent = (k-1)/refine3+1
       kgrid   = k-start3+gstart3
       do i = start1, end1
         iparent = (i-1)/refine1+1
         igrid   = i-start1+gstart1
c
c        Compute the average along a single line of the interpolated values
c
         do j = 1, wdim2
            temp(j) = 0._RKIND
         enddo
c
         do j1 = gstart2, gstart2+refine2-1
            do j2 = 1, wdim2
               temp(j2) = temp(j2) + 
     &              grid(igrid, (j2-1)*refine2+j1, kgrid)/REAL(refine2)
            enddo
         enddo
c
c        Check to see if the end of the distribution is a local maxima
c         (compared to the nearest interpolated grid point to the left
c          and the mean value along this interpolated line -- computed above).
c
         del1next = parent(iparent, (start2-1)/refine2, kparent)*dV
         do j = 1, wdim2
            del1left = del1next - 
     &                 grid(igrid, (j-1)*refine2+gstart2  , kgrid)
            del1next = grid(igrid,  j   *refine2+gstart2-1, kgrid)
            del2left = temp(j) - 
     &                 grid(igrid, (j-1)*refine2+gstart2  , kgrid)
c
            if (abs(del2left) .lt. tiny) 
     &           del2left = sign(REAL(tiny,RKIND), del2left)
            fraction(j) = min(max(del1left/del2left, 0._RKIND),1._RKIND)
         enddo
c
c        repeat for right side of each distribution
c
         del1next = parent(iparent, end2/refine2+1, kparent)*dV
         do j = wdim2, 1, -1
            del1right = del1next -
     &                  grid(igrid,  j   *refine2+gstart2-1, kgrid)
            del1next  = grid(igrid, (j-1)*refine2+gstart2  , kgrid)
            del2right = temp(j) -
     &                  grid(igrid,  j   *refine2+gstart2-1, kgrid)
c
            if (abs(del2right) .lt. tiny) 
     &           del2right = sign(REAL(tiny,RKIND), del2right)
            fraction(j) =  max(max(min(del1right/del2right, 0._RKIND),
     &                         1._RKIND), fraction(j)               )
         enddo
c
c        Modify distributions based on the computed fraction decrease
c
         do j1 = gstart2, gstart2+refine2-1
            do j2 = 1, wdim2
               grid(igrid, (j2-1)*refine2+j1, kgrid) = 
     &                     fraction(j2)  * temp(j2) +
     &              (1._RKIND - fraction(j2)) * grid(igrid, 
     &                                         (j2-1)*refine2+j1, kgrid)
            enddo
         enddo
c
       enddo
      enddo
c
c  3) Repeat for the k-direction
c
      do i = start1, end1
        iparent = (i-1)/refine1+1
        igrid   = i-start1+gstart1
        do j = start2, end2
          jparent = (j-1)/refine2+1
          jgrid   = j-start2+gstart2
c
c        Compute the average along a single line of the interpolated values
c
         do k = 1, wdim3
            temp(k) = 0._RKIND
         enddo
c
         do k1 = gstart3, gstart3+refine3-1
            do k2 = 1, wdim3
               temp(k2) = temp(k2) + 
     &              grid(igrid, jgrid, (k2-1)*refine3+k1)/REAL(refine3)
            enddo
         enddo
c
c        Check to see if the end of the distribution is a local maxima
c         (compared to the nearest interpolated grid point to the left
c          and the mean value along this interpolated line -- computed above).
c
         del1next = parent(iparent, jparent, (start3-1)/refine3)*dV
         do k = 1, wdim3
            del1left = del1next - 
     &                 grid(igrid, jgrid, (k-1)*refine3+gstart3  )
            del1next = grid(igrid, jgrid,  k   *refine3+gstart3-1)
            del2left = temp(k) - 
     &                 grid(igrid, jgrid, (k-1)*refine3+gstart3  )
c
            if (abs(del2left) .lt. tiny) 
     &           del2left = sign(REAL(tiny,RKIND), del2left)
            fraction(k) = min(max(del1left/del2left, 0._RKIND),1._RKIND)
         enddo
c
c        repeat for right side of each distribution
c
         del1next = parent(iparent, jparent, end3/refine3+1)*dV
         do k = wdim3, 1, -1
            del1right = del1next -
     &                  grid(igrid, jgrid,  k   *refine3+gstart3-1)
            del1next  = grid(igrid, jgrid, (k-1)*refine3+gstart3  )
            del2right = temp(k) -
     &                  grid(igrid, jgrid,  k   *refine3+gstart3-1)
c
            if (abs(del2right) .lt. tiny) 
     &           del2right = sign(REAL(tiny,RKIND), del2right)
            fraction(k) =  max(max(min(del1right/del2right, 0._RKIND),
     &                         1._RKIND), fraction(k)                )
         enddo
c
c        Modify distributions based on the computed fraction decrease
c
         do k1 = gstart3, gstart3+refine3-1
            do k2 = 1, wdim3
               grid(igrid, jgrid, (k2-1)*refine3+k1) = 
     &                     fraction(k2)  * temp(k2) +
     &              (1._RKIND - fraction(k2)) * grid(igrid, jgrid,
     &                                         (k2-1)*refine3+k1)
            enddo
         enddo
c
       enddo
      enddo
c
c     Next iteration
c
      enddo
c
c     Remove the irritating division by the refinement volume
c
      do k = gstart3, gstart3 + end3-start3
         do j = gstart2, gstart2 + end2-start2
            do i = gstart1, gstart1 + end1-start1
               grid(i,j,k) = grid(i,j,k) * refine1*refine2*refine3
            enddo
         enddo
      enddo
c
      return
      end

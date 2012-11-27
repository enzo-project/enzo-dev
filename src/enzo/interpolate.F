#include "error.def"
c=======================================================================
c//////////////////////  SUBROUTINE INTERPOLATE  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine interpolate(ndim, parent, pdims, pstart, pend, refine,
     &                       grid, gdims, gstart, work, imethod,
     &                       iposflag, ierror)
c
c  PERFORMS A TSC-LIKE INTERPOLATION FROM THE FIELD PARENT TO GRID
c
c     written by: Greg Bryan
c     date:       November, 1994
c     modified1:
c
c  PURPOSE:  This routine takes the field parent and interpolates it using
c     a varient of the third-order triangular-shaped-cloud interpolation
c     scheme.
c     NOTE: There is a restriction.  The interpolation must be done in by
c        blocks of the parent grid.
c
c  INPUTS:
c     gdims(3)  - declared dimension of grid
c     gstart(3) - starting offset of refined region in grid
c     imethod   - interpolation method (0 = ThirdOrderA, 1 = SecondOrderA,
c                                       2 = SecondOrderB)
c     iposflag  - flag indicating what kind of positivity/monotonicty is
c                 required (only used if imethod = 2) 
c                 (0 - none, 1 - positivity, 2 - monotonicity)
c     ndim      - rank of fields
c     parent    - parent field
c     pdims(3)  - declared dimension of parent
c     pend(3)   - ending dimension in parent in index units of grid
c     pstart(3) - starting dimension in parent in index units of grid
c     refine(3) - INTG_PREC refinement factors
c     work      - array of size Multiply[i=1,ndim][(gdims(i)/refine(i) + 1)] 
c
c  OUTPUTS:
c     grid      - grid with refined region
c
c  EXTERNALS:
c
c  LOCALS:
c-----------------------------------------------------------------------
      implicit NONE
#include "fortran_types.def"
c
c  arguements
c
      INTG_PREC imethod, iposflag, ndim, ierror
      INTG_PREC pdims(ndim), pstart(ndim), pend(ndim), refine(ndim),
     &        gdims(ndim), gstart(ndim)
      R_PREC    grid(1), parent(1), work(1)
c
c  locals
c
      INTG_PREC i
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c  Error check for number of dimensions
c
      if (ndim .lt. 1 .or. ndim .gt. 3) then
         write (6,*) "INTERPOLATE: ERROR ndim =", ndim
         ERROR_MESSAGE
      endif
c
c  Error check for interpolation method
c
      if (imethod .lt. 0 .or. imethod .gt. 4) then
         write (6,*) "INTERPOLATE: ERROR imethod =", imethod
         ERROR_MESSAGE
      endif
c
c  Error check for bounds (must have one parent cell on each side)
c
      do i = 1, ndim
         if (pstart(i) .lt. refine(i) .or. 
     &       pend(i) .ge. (pdims(i)-1)*refine(i)) then
            write(6,1030) pstart(i), pend(i), pdims(i), i
 1030       format("INTERPOLATE: ERROR pstart =", i5, "  pend =", i5,
     &             "  pdims =", i5, "  i =", i2)
            ERROR_MESSAGE
         endif
      enddo
c
c  Error check for self-consistent interpolation indexes
c
      do i = 1, ndim
         if (int(pstart(i)/refine(i),IKIND)*refine(i) /= pstart(i)) then
            write(6,*) "INTERPOLATE: ERROR pstart(i) =", pstart(i),
     &                 " refine(i) =", refine(i), " i =", i
            ERROR_MESSAGE
         endif
      enddo
c
c  A) ThirdOrderA
c    Do the appropriate interpolation (add one to convert to fortran 1 based)
c
      if (imethod .eq. 0) then
c
       if (ndim .eq. 1) call tscint1d(parent, pdims(1), 
     &        pstart(1)+1_IKIND, pend(1)+1_IKIND, refine(1),
     &        grid, gdims(1), gstart(1)+1_IKIND)
c
       if (ndim .eq. 2) call tscint2d(parent, pdims(1), pdims(2),
     &      pstart(1)+1_IKIND, pstart(2)+1_IKIND,
     &      pend(1)+1_IKIND, pend(2)+1_IKIND,
     &      refine(1), refine(2),
     &      grid, gdims(1), gdims(2),
     &      gstart(1)+1_IKIND, gstart(2)+1_IKIND)
c
       if (ndim .eq. 3) call tscint3d(parent, 
     &      pdims(1), pdims(2), pdims(3),
     &      pstart(1)+1_IKIND, pstart(2)+1_IKIND, pstart(3)+1_IKIND,
     &      pend(1)+1_IKIND, pend(2)+1_IKIND, pend(3)+1_IKIND,
     &      refine(1), refine(2), refine(3),
     &      grid, gdims(1), gdims(2), gdims(3),
     &      gstart(1)+1_IKIND, gstart(2)+1_IKIND, gstart(3)+1_IKIND)
c
      endif
c
c  B) SecondOrderA
c
      if (imethod .eq. 1 .or. 
     &    (imethod .eq. 2 .and. iposflag .eq. 2)) then
c
         if (ndim .eq. 1) call interp1d(parent, work, 
     &                            pdims(1),
     &                            pstart(1)+1_IKIND,
     &                            pend(1)+1_IKIND,
     &                            refine(1),
     &                            grid, gdims(1),
     &                            gstart(1)+1_IKIND,
     &                            gdims(1)/refine(1)+1_IKIND)
c
         if (ndim .eq. 2) call interp2d(parent, work, 
     &                            pdims(1), pdims(2),
     &                            pstart(1)+1_IKIND, pstart(2)+1_IKIND,
     &                            pend(1)+1_IKIND, pend(2)+1_IKIND,
     &                            refine(1), refine(2),
     &                            grid, gdims(1), gdims(2),
     &                            gstart(1)+1_IKIND, gstart(2)+1_IKIND,
     &                            gdims(1)/refine(1)+1_IKIND, 
     &                            gdims(2)/refine(2)+1_IKIND)
c
         if (ndim .eq. 3) call interp3d(parent, work, 
     &        pdims(1), pdims(2), pdims(3),
     &        pstart(1)+1_IKIND, pstart(2)+1_IKIND, pstart(3)+1_IKIND,
     &        pend(1)+1_IKIND, pend(2)+1_IKIND, pend(3)+1_IKIND,
     &        refine(1), refine(2), refine(3),
     &        grid, gdims(1), gdims(2), gdims(3),
     &        gstart(1)+1_IKIND, gstart(2)+1_IKIND, gstart(3)+1_IKIND,
     &        gdims(1)/refine(1)+1_IKIND, 
     &        gdims(2)/refine(2)+1_IKIND,
     &        gdims(3)/refine(3)+1_IKIND, ierror)
c
      endif
c
c  C) SecondOrderB
c
      if (imethod .eq. 2 .and. iposflag .ne. 2) then
c
#if 0
         if (ndim .eq. 1) call interp1d(parent, work, 
     &        pdims(1),
     &        pstart(1)+1_IKIND,
     &        pend(1)+1_IKIND,
     &        refine(1),
     &        grid, gdims(1),
     &        gstart(1)+1_IKIND,
     &        gdims(1)/refine(1)+1_IKIND)
c
         if (ndim .eq. 2) call interp2d(parent, work, 
     &        pdims(1), pdims(2),
     &        pstart(1)+1_IKIND, pstart(2)+1_IKIND,
     &        pend(1)+1_IKIND, pend(2)+1_IKIND,
     &        refine(1), refine(2),
     &        grid, gdims(1), gdims(2),
     &        gstart(1)+1_IKIND, gstart(2)+1_IKIND,
     &        gdims(1)/refine(1)+1_IKIND, 
     &        gdims(2)/refine(2)+1_IKIND)
#endif
         if (ndim .ne. 3) then
            write(6,*) 'INTERPOLATE: imethod=2 only for 3d'
            ERROR_MESSAGE
         endif
c
         if (ndim .eq. 3) call int_lin3d(parent, work, 
     &        pdims(1), pdims(2), pdims(3),
     &        pstart(1)+1_IKIND, pstart(2)+1_IKIND, pstart(3)+1_IKIND,
     &        pend(1)+1_IKIND, pend(2)+1_IKIND, pend(3)+1_IKIND,
     &        refine(1), refine(2), refine(3),
     &        grid, gdims(1), gdims(2), gdims(3),
     &        gstart(1)+1_IKIND, gstart(2)+1_IKIND, gstart(3)+1_IKIND,
     &        gdims(1)/refine(1)+1_IKIND, 
     &        gdims(2)/refine(2)+1_IKIND,
     &        gdims(3)/refine(3)+1_IKIND, iposflag)
c
      endif
c
c  D) SecondOrderC
c
      if (imethod .eq. 3) then
c
         call cicinterp(parent, pdims(1), pdims(2), pdims(3), ndim,
     &        pstart(1)+1_IKIND, pstart(2)+1_IKIND, pstart(3)+1_IKIND,
     &        pend(1)+1_IKIND, pend(2)+1_IKIND, pend(3)+1_IKIND,
     &        refine(1), refine(2), refine(3),
     &        grid, gdims(1), gdims(2), gdims(3),
     &        gstart(1)+1_IKIND, gstart(2)+1_IKIND, gstart(3)+1_IKIND, 
     &        iposflag)
c
      endif
c
c  E) FirstOrderA
c
      if (imethod .eq. 4) then
c
         call ngpinterp(parent, pdims(1), pdims(2), pdims(3), ndim,
     &        pstart(1)+1_IKIND, pstart(2)+1_IKIND, pstart(3)+1_IKIND,
     &        pend(1)+1_IKIND, pend(2)+1_IKIND, pend(3)+1_IKIND,
     &        refine(1), refine(2), refine(3),
     &        grid, gdims(1), gdims(2), gdims(3),
     &        gstart(1)+1_IKIND, gstart(2)+1_IKIND, gstart(3)+1_IKIND)
c
      endif
c
      return
      end

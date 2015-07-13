#include "error.def"
c=======================================================================
c///////////////////////  SUBROUTINE INTERP1D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine interp1d(parent, work, dim1, 
     &                    start1, end1, refine1, grid,
     &                    gdim1, gstart1, wdim1)
c
c  PERFORMS A 1D LINEAR, MONOTONIC, CONSERVATIVE INTERPOLATION
c
c     written by: Greg Bryan
c     date:       August, 1995
c     modified1:  
c
c  PURPOSE:  This routine takes the field parent and interpolates it using
c     a linear, monotonic, conservative interpolation technique.
c
c     NOTE: There is a restriction.  The interpolation must be done in
c        blocks of the parent grid.
c
c  INPUTS:
c     parent      - parent field
c     dim1        - declared dimension of parent
c     start1      - starting index in parent in units of grid (one based)
c     end1        - ending index in parent in units of grid (one based)
c     refine1     - integer refinement factors
c     gdim1       - declared dimension of grid
c     gstart1     - starting offset of refined region in grid (one based)
c     wdim1       - work dimensions
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
      INTG_PREC MAX_REFINE
      parameter (MAX_REFINE = 32)
c
c-----------------------------------------------------------------------
c
c  arguments
c
      INTG_PREC dim1, start1, end1, refine1, gdim1, gstart1, wdim1
      R_PREC    parent(dim1), grid(gdim1), work(wdim1)
c
c  locals
c
      INTG_PREC i, i1, pstart1, pdim1
      R_PREC    delf0, fbar, fx, one
      parameter (one=1._RKIND)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c     Error check
c
      if (refine1 .gt. MAX_REFINE) then
         write (6,*) "interp1d: Error: refine1 > MAX_REFINE"
         ERROR_MESSAGE
      endif
c
c     Set work dimensions (number of coarse cell in refinedment in each dim)
c
      pdim1 = (end1-start1+1)/refine1
c
c     Set parent start indexes (note: zero-based)
c
      pstart1 = (start1-1)/refine1
c
c     Interpolate values from cell centers to cell corners
c
      do i = 1, pdim1+1
         work(i) = 0.5_RKIND*(
     &           parent(pstart1+i-1) +
     &           parent(pstart1+i  )
     &                 )
      enddo
c
c     Loop over coarse cells and compute the interpolation coefficients.
c
      do i = 1, pdim1
c
c        Store the average value of this cell in fbar for easy reference
c
         fbar  = parent(pstart1+i)
c
c        Compute the minimum delta f across the cell
c
         delf0 =     min(abs(fbar - work(i  )),
     &                   abs(fbar - work(i+1)) )*
     &            sign(one, (fbar - work(i  )))
         if ((work(i+1) - fbar) * (fbar - work(i)) .lt. 0.0) 
     &        delf0=0.0
c
c        Now, find the the interpolation coefficients.
c
         fx = 2._RKIND*delf0
c
c        And interpolate to the subgrid.
c
         do i1 = 0, refine1-1
            grid(gstart1+(i-1)*refine1+i1) = (fbar +
     &            (REAL(i1)+0.5_RKIND - 
     &             0.5_RKIND*REAL(refine1)) / 
     &             REAL(refine1)*fx)
         enddo
c
      enddo
c
      return
      end

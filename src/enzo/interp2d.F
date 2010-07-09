#include "error.def"
c=======================================================================
c///////////////////////  SUBROUTINE INTERP2D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine interp2d(parent, work, dim1, dim2,
     &                    start1, start2,
     &                    end1, end2, 
     &                    refine1, refine2, grid,
     &                    gdim1, gdim2, 
     &                    gstart1, gstart2, 
     &                    wdim1, wdim2)
c
c  PERFORMS A 2D LINEAR, MONOTONIC, CONSERVATIVE INTERPOLATION
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
c     dim1,2      - declared dimension of parent
c     start1,2    - starting index in parent in units of grid (one based)
c     end1,2      - ending index in parent in units of grid (one based)
c     refine1,2   - integer refinement factors
c     gdim1,2     - declared dimension of grid
c     gstart1,2   - starting offset of refined region in grid (one based)
c     wdim1,2     - work dimensions
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
c
      integer MAX_REFINE
      parameter (MAX_REFINE = 32)
c
c-----------------------------------------------------------------------
c
c  arguments
c
      integer dim1, dim2, start1, start2,
     &        end1, end2, refine1, refine2,
     &        gdim1, gdim2, gstart1, gstart2, 
     &        wdim1, wdim2
      real    parent(dim1, dim2), grid(gdim1, gdim2),
     &        work(wdim1, wdim2)
c
c  locals
c
      integer i, i1, j, j1, pstart1, pstart2, pdim1, pdim2
      real    delf0, delf1, fbar, fx, fy
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c     Error check
c
      if (refine1 .gt. MAX_REFINE .or. refine2 .gt. MAX_REFINE) then
         write (6,*) "interp2d: Error: refine1 or 2 > MAX_REFINE"
         ERROR_MESSAGE
      endif
c
c     Set work dimensions (number of coarse cell in refinedment in each dim)
c
      pdim1 = (end1-start1+1)/refine1
      pdim2 = (end2-start2+1)/refine2
c
c     Set parent start indexes (note: zero-based)
c
      pstart1 = (start1-1)/refine1
      pstart2 = (start2-1)/refine2
c
c     Interpolate values from cell centers to cell corners
c
      do j = 1, pdim2+1
         do i = 1, pdim1+1
            work(i, j) = 0.25*(
     &           parent(pstart1+i-1, pstart2+j-1) +
     &           parent(pstart1+i-1, pstart2+j  ) +
     &           parent(pstart1+i  , pstart2+j-1) +
     &           parent(pstart1+i  , pstart2+j  )
     &                        )
         enddo
      enddo
c
c     Loop over coarse cells and compute the interpolation coefficients.
c
      do j = 1, pdim2
         do i = 1, pdim1
c
c           Store the average value of this cell in fbar for easy reference
c
            fbar  = parent(pstart1+i, pstart2+j)
c
c            if (fbar .lt. 0.0) write(6,*) 'a',fbar, pstart1+i, pstart2+j
c
c           Compute the minimum delta f across both diagonals
c
            delf0 =     min(abs(fbar - work(i  , j  )),
     &                      abs(fbar - work(i+1, j+1)) )*
     &               sign(1.0, (fbar - work(i  , j  )))
            if ((work(i+1, j+1) - fbar) * (fbar - work(i, j)) .lt. 0.0)
     &         delf0 = 0.0
c
            delf1 =     min(abs(fbar - work(i+1, j  )),
     &                      abs(fbar - work(i  , j+1)) )*
     &               sign(1.0, (fbar - work(i+1, j  )))
            if ((work(i, j+1) - fbar) * (fbar - work(i+1, j)) .lt. 0.0)
     &         delf1 = 0.0
c
c           Now, find the the interpolation coefficients.
c
            fx = delf0 - delf1
            fy = delf0 + delf1
c
c              And interpolate to the subgrid.
c
            do j1 = 0, refine2-1
               do i1 = 0, refine1-1
                  grid(gstart1+(i-1)*refine1+i1,
     &                 gstart2+(j-1)*refine2+j1) = (fbar +
     &               (real(i1)+0.5-0.5*real(refine1))/real(refine1)*fx +
     &               (real(j1)+0.5-0.5*real(refine2))/real(refine2)*fy )
c                  if (grid(gstart1+(i-1)*refine1+i1,
c     &                     gstart2+(j-1)*refine2+j1) .lt. 0.0)
c     &                write(6,*) 'b',i,j,i1,j1
               enddo
            enddo
c
         enddo
      enddo
c
      return
      end

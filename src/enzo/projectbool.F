#include "fortran.def"
c=======================================================================
c////////////////////////  SUBROUTINE PROJECT  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine projectbool(rank, i1, i2, i3, iline, pdim, field, line)
c
c  PROJECTS A 2 OR 3D BOOLEAN FIELD TO A LINE
c
c  written by: Greg Bryan
c  date:       October, 1995
c  modified1:
c
c  PURPOSE:
c      A good fortran compiler will pull the loops (below) apart and
c        produce fast (albeit long) code.
c
c  INPUTS:
c     field    - 2 or 3D field
c     i1,i2,i3 - dimensions of field
c     iline    - dimension of line
c     pdim     - direction of projection (0,1 or 2)
c     rank     - rank of field
c
c  OUTPUT ARGUMENTS: 
c     line     - projected line
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
      INTG_PREC i1, i2, i3, iline, pdim, rank
      INTG_PREC line(iline)
      logical*1 field(i1, i2, i3)
c
c  locals
c
      INTG_PREC i, j, k
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     Clear line
c
      do i = 1, iline
         line(i) = 0
      enddo
c
c     Do projection
c
      do k = 1, i3
         do j = 1, i2
            do i = 1, i1
               if (field(i,j,k)) then
                  if (pdim .eq. 0) line(i) = line(i) + 1
                  if (pdim .eq. 1) line(j) = line(j) + 1
                  if (pdim .eq. 2) line(k) = line(k) + 1
               endif
            enddo
         enddo
      enddo
c
      return
      end

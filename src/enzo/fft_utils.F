c=======================================================================
c/////////////////////////  SUBROUTINE COPY3DFT  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine copy3dft(source, dest, sdim1, sdim2, sdim3, 
     &                  ddim1, ddim2, ddim3,
     &                  sstart1, sstart2, sstart3, 
     &                  dstart1, dstart2, dstart3)
c
c  COPIES PART OF FIELD FROM TO FIELD TO (FORWARD TRANSPOSE)
c
c  written by: Greg Bryan
c  date:       October, 1995
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source, dest - source and destination fields
c     sdim1-3      - source dimension
c     ddim1-3      - destination dimension
c     sstart1-3    - start of source field in global coordinates
c     dstart1-3    - start of destination field in global coordinates
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
      INTG_PREC ddim1, ddim2, ddim3, dstart1, dstart2, dstart3, 
     &        sdim1, sdim2, sdim3, sstart1, sstart2, sstart3
      R_PREC  source(sdim1, sdim2, sdim3), dest(ddim3*2, ddim2, ddim1/2)
c
c  locals
c
      INTG_PREC i, j, k, end1, end2, end3, start1, start2, start3
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     determine the overlap area in global coordinates
c
      start1 = max(sstart1, dstart1)
      start2 = max(sstart2, dstart2)
      start3 = max(sstart3, dstart3)
c
      end1 = min(sstart1+sdim1, dstart1+ddim1) - 1
      end2 = min(sstart2+sdim2, dstart2+ddim2) - 1
      end3 = min(sstart3+sdim3, dstart3+ddim3) - 1
c
c      write(6,*) 'start',sstart1, sstart2, sstart3
c      write(6,*) dstart1, dstart2, dstart3
c      write(6,*) 'dim',sdim1, sdim2, sdim3
c      write(6,*) ddim1, ddim2, ddim3
c      write(6,*) 'se',start1, start2, start3
c      write(6,*) end1, end2, end3
c
      do k = start3, end3
         do j = start2, end2
            do i = start1, end1, 2
               dest  ((k-dstart3)*2+1, j-dstart2+1, (i-dstart1)/2+1) =
     &         source( i-sstart1   +1, j-sstart2+1,  k-sstart3   +1)
               dest  ((k-dstart3)*2+2, j-dstart2+1, (i-dstart1)/2+1) =
     &         source( i-sstart1   +2, j-sstart2+1,  k-sstart3   +1)
            enddo
         enddo
      enddo
c
      return
      end


c=======================================================================
c/////////////////////////  SUBROUTINE COPY3DRT  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine copy3drt(source, dest, sdim1, sdim2, sdim3, 
     &                  ddim1, ddim2, ddim3,
     &                  sstart1, sstart2, sstart3, 
     &                  dstart1, dstart2, dstart3)
c
c  COPIES PART OF FIELD FROM TO FIELD TO (FORWARD TRANSPOSE)
c
c  written by: Greg Bryan
c  date:       October, 1995
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source, dest - source and destination fields
c     sdim1-3      - source dimension
c     ddim1-3      - destination dimension
c     sstart1-3    - start of source field in global coordinates
c     dstart1-3    - start of destination field in global coordinates
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
      INTG_PREC ddim1, ddim2, ddim3, dstart1, dstart2, dstart3, 
     &        sdim1, sdim2, sdim3, sstart1, sstart2, sstart3
      R_PREC  source(sdim3*2, sdim2, sdim1/2), dest(ddim1, ddim2, ddim3)
c
c  locals
c
      INTG_PREC i, j, k, end1, end2, end3, start1, start2, start3
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     determine the overlap area in global coordinates
c
      start1 = max(sstart1, dstart1)
      start2 = max(sstart2, dstart2)
      start3 = max(sstart3, dstart3)
c
      end1 = min(sstart1+sdim1, dstart1+ddim1) - 1
      end2 = min(sstart2+sdim2, dstart2+ddim2) - 1
      end3 = min(sstart3+sdim3, dstart3+ddim3) - 1
c
      do k = start3, end3
         do j = start2, end2
            do i = start1, end1, 2
               dest  ( i-dstart1   +1, j-dstart2+1,  k-dstart3   +1) =
     &         source((k-sstart3)*2+1, j-sstart2+1, (i-sstart1)/2+1)
               dest  ( i-dstart1   +2, j-dstart2+1,  k-dstart3   +1) =
     &         source((k-sstart3)*2+2, j-sstart2+1, (i-sstart1)/2+1)
            enddo
         enddo
      enddo
c
      return
      end

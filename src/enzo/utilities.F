#include "fortran.def"
c=======================================================================
c/////////////////////////  SUBROUTINE COPY3D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine copy3d(source, dest, sdim1, sdim2, sdim3, 
     &                  ddim1, ddim2, ddim3,
     &                  sstart1, sstart2, sstart3, 
     &                  dstart1, dstart2, dstart3)
c
c  COPIES PART OF FIELD FROM TO FIELD TO
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
     &          sdim1, sdim2, sdim3, sstart1, sstart2, sstart3
      R_PREC    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
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
      end1 = min(sstart1+sdim1, dstart1+ddim1)  - 1
      end2 = min(sstart2+sdim2, dstart2+ddim2)  - 1
      end3 = min(sstart3+sdim3, dstart3+ddim3)  - 1
c
      do k = start3, end3
         do j = start2, end2
            do i = start1, end1
               dest  (i-dstart1+1, j-dstart2+1, k-dstart3+1) =
     &         source(i-sstart1+1, j-sstart2+1, k-sstart3+1)
            enddo
         enddo
      enddo
c
      return
      end

c=======================================================================
c///////////////////////  SUBROUTINE COPY3DREL  \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine copy3drel(source, dest,
     &                     dim1, dim2, dim3,
     &                     sdim1, sdim2, sdim3,
     &                     ddim1, ddim2, ddim3,
     &                     sstart1, sstart2, sstart3, 
     &                     dstart1, dstart2, dstart3)
c
c  COPIES PART OF FIELD FROM TO FIELD TO
c
c  written by: Greg Bryan
c  date:       October, 1995
c  modified1:  John Wise (May, 2010) for positions relative to 
c              their grids
c
c  PURPOSE:
c
c  INPUTS:
c     source, dest - source and destination fields
c     dim1-3       - dimension of copied region
c     sdim1-3      - source dimension
c     ddim1-3      - destinaation dimension
c     sstart1-3    - start of source field in local grid coordinates
c     dstart1-3    - start of destination field in local grid coordinates
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
     &          sdim1, sdim2, sdim3, sstart1, sstart2, sstart3,
     &          dim1, dim2, dim3
      R_PREC    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
c
c  locals
c
      INTG_PREC i, j, k
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     determine the overlap area in global coordinates
c
      do k = 1, dim3
         do j = 1, dim2
            do i = 1, dim1
               dest  (i+dstart1, j+dstart2, k+dstart3) =
     &         source(i+sstart1, j+sstart2, k+sstart3)
            enddo
         enddo
      enddo
c
      return
      end


c=======================================================================
c////////////////////////  SUBROUTINE COPY3DINT  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine copy3dint(source, dest, sdim1, sdim2, sdim3, 
     &                     ddim1, ddim2, ddim3,
     &                     sstart1, sstart2, sstart3, 
     &                     dstart1, dstart2, dstart3)
c
c  COPIES PART OF FIELD (INTG_PREC VERSION)
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
     &          sdim1, sdim2, sdim3, sstart1, sstart2, sstart3
      INTG_PREC source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
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
            do i = start1, end1
               dest  (i-dstart1+1, j-dstart2+1, k-dstart3+1) =
     &         source(i-sstart1+1, j-sstart2+1, k-sstart3+1)
            enddo
         enddo
      enddo
c
      return
      end

c=======================================================================
c////////////////////////  SUBROUTINE COPY3DBOOL  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine copy3dbool(source, dest, sdim1, sdim2, sdim3, 
     &                      ddim1, ddim2, ddim3,
     &                      sstart1, sstart2, sstart3, 
     &                      dstart1, dstart2, dstart3)
c
c  COPIES PART OF FIELD (BOOLEAN VERSION)
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
     &          sdim1, sdim2, sdim3, sstart1, sstart2, sstart3
      logical*1 source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
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
            do i = start1, end1
               dest  (i-dstart1+1, j-dstart2+1, k-dstart3+1) =
     &         source(i-sstart1+1, j-sstart2+1, k-sstart3+1)
            enddo
         enddo
      enddo
c
      return
      end



c=======================================================================
c/////////////////////////  SUBROUTINE MULT3D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine mult3d(source, dest, sdim1, sdim2, sdim3, 
     &                  ddim1, ddim2, ddim3,
     &                  sstart1, sstart2, sstart3, 
     &                  dstart1, dstart2, dstart3)
c
c  PERFORMS THE FIELDS OPERATION: DEST = DEST * SOURCE
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
     &          sdim1, sdim2, sdim3, sstart1, sstart2, sstart3
      R_PREC    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
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
            do i = start1, end1
               dest  (i-dstart1+1, j-dstart2+1, k-dstart3+1) =
     &         dest  (i-dstart1+1, j-dstart2+1, k-dstart3+1) *
     &         source(i-sstart1+1, j-sstart2+1, k-sstart3+1)
            enddo
         enddo
      enddo
c
      return
      end


c=======================================================================
c/////////////////////////  SUBROUTINE DIV3D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine div3d(source, dest, sdim1, sdim2, sdim3, 
     &                  ddim1, ddim2, ddim3,
     &                  sstart1, sstart2, sstart3, 
     &                  dstart1, dstart2, dstart3,
     &                  rstart1, rstart2, rstart3,
     &                  rend1, rend2, rend3)
c
c  PERFORMS THE FIELDS OPERATION: DEST = DEST / SOURCE
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
c     rstart1-3    - start of region to operate on
c     rend1-3      - end (inclusive) of region to operate on
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
     &          sdim1, sdim2, sdim3, sstart1, sstart2, sstart3,
     &          rstart1, rstart2, rstart3, rend1, rend2, rend3
      R_PREC    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
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
      start1 = max(sstart1, dstart1, rstart1)
      start2 = max(sstart2, dstart2, rstart2)
      start3 = max(sstart3, dstart3, rstart3)
c
      end1 = min(sstart1+sdim1-1, dstart1+ddim1-1, rend1)
      end2 = min(sstart2+sdim2-1, dstart2+ddim2-1, rend2)
      end3 = min(sstart3+sdim3-1, dstart3+ddim3-1, rend3)
c
      do k = start3, end3
         do j = start2, end2
            do i = start1, end1
               dest  (i-dstart1+1, j-dstart2+1, k-dstart3+1) =
     &         dest  (i-dstart1+1, j-dstart2+1, k-dstart3+1) /
     &         source(i-sstart1+1, j-sstart2+1, k-sstart3+1)
            enddo
         enddo
      enddo
c
      return
      end



c=======================================================================
c///////////////////////  SUBROUTINE COMBINDE3D  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine combine3d(source1, weight1, source2, weight2,
     &                     dest, sdim1, sdim2, sdim3, 
     &                     ddim1, ddim2, ddim3,
     &                     sstart1, sstart2, sstart3, 
     &                     dstart1, dstart2, dstart3,
     &                     ivel_shift, irefine)
c
c  COPIES PART OF FIELD FROM TO FIELD TO
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
      INTG_PREC ddim1, ddim2, ddim3, dstart1, dstart2, dstart3, irefine,
     &        sdim1, sdim2, sdim3, sstart1, sstart2, sstart3, ivel_shift
      R_PREC  source1(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3),
     &        source2(sdim1, sdim2, sdim3), weight1, weight2
c
c  locals
c
      INTG_PREC i, j, k, end1, end2, end3, start1, start2, start3
      R_PREC    fact1, fact2
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
      if (ivel_shift .eq. 0) then
       do k = start3, end3
          do j = start2, end2
             do i = start1, end1
                dest   (i-dstart1+1, j-dstart2+1, k-dstart3+1) =
     &          source1(i-sstart1+1, j-sstart2+1, k-sstart3+1)*weight1 +
     &          source2(i-sstart1+1, j-sstart2+1, k-sstart3+1)*weight2
             enddo
          enddo
       enddo
      endif
c
      if (ivel_shift .eq. 1) then
       do k = start3, end3
          do j = start2, end2
             fact1 = (REAL(irefine,RKIND)+1.0)/(2.0*REAL(irefine,RKIND))
             fact2 = (REAL(irefine,RKIND)-1.0)/(2.0*REAL(irefine,RKIND))
             do i = start1, end1
                dest   (i-dstart1+1, j-dstart2+1, k-dstart3+1) = fact1*(
     &          source1(i-sstart1+1, j-sstart2+1, k-sstart3+1)*weight1 +
     &          source2(i-sstart1+1, j-sstart2+1, k-sstart3+1)*weight2 )
     &                                                         + fact2*(
     &          source1(i-sstart1+2, j-sstart2+1, k-sstart3+1)*weight1 +
     &          source2(i-sstart1+2, j-sstart2+1, k-sstart3+1)*weight2 )
             enddo
          enddo
       enddo
      endif
c
      if (ivel_shift .eq. 2) then
       do k = start3, end3
          do j = start2, end2
             fact1 = (REAL(irefine,RKIND)+1.0)/(2.0*REAL(irefine,RKIND))
             fact2 = (REAL(irefine,RKIND)-1.0)/(2.0*REAL(irefine,RKIND))
             do i = start1, end1
                dest   (i-dstart1+1, j-dstart2+1, k-dstart3+1) = fact1*(
     &          source1(i-sstart1+1, j-sstart2+1, k-sstart3+1)*weight1 +
     &          source2(i-sstart1+1, j-sstart2+1, k-sstart3+1)*weight2 )
     &                                                         + fact2*(
     &          source1(i-sstart1+1, j-sstart2+2, k-sstart3+1)*weight1 +
     &          source2(i-sstart1+1, j-sstart2+2, k-sstart3+1)*weight2 )
             enddo
          enddo
       enddo
      endif
c
      if (ivel_shift .eq. 3) then
       do k = start3, end3
          do j = start2, end2
             fact1 = (REAL(irefine,RKIND)+1.0)/(2.0*REAL(irefine,RKIND))
             fact2 = (REAL(irefine,RKIND)-1.0)/(2.0*REAL(irefine,RKIND))
             do i = start1, end1
                dest   (i-dstart1+1, j-dstart2+1, k-dstart3+1) = fact1*(
     &          source1(i-sstart1+1, j-sstart2+1, k-sstart3+1)*weight1 +
     &          source2(i-sstart1+1, j-sstart2+1, k-sstart3+1)*weight2 )
     &                                                         + fact2*(
     &          source1(i-sstart1+2, j-sstart2+1, k-sstart3+2)*weight1 +
     &          source2(i-sstart1+2, j-sstart2+1, k-sstart3+2)*weight2 )
             enddo
          enddo
       enddo
      endif
c
      return
      end


c=======================================================================
c/////////////////////////  SUBROUTINE COPY1D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine copy1d(source, dest, sdim1, ddim1, sstart1, dstart1)
c
c  COPIES PART OF FIELD FROM TO FIELD TO
c
c  written by: Greg Bryan
c  date:       October, 1995
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source, dest - source and destination fields
c     sdim1      - source dimension
c     ddim1      - destination dimension
c     sstart1    - start of source field in global coordinates
c     dstart1    - start of destination field in global coordinates
c
c  OUTPUT ARGUMENTS: 
c     dest       - destination field
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
      INTG_PREC ddim1, dstart1, sdim1, sstart1
      R_PREC    source(sdim1), dest(ddim1)
c
c  locals
c
      INTG_PREC i, end1, start1
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     determine the overlap area in global coordinates
c
      start1 = max(sstart1, dstart1)
c
      end1 = min(sstart1+sdim1, dstart1+ddim1) - 1
c
      do i = start1, end1
         dest  (i-dstart1+1) = source(i-sstart1+1)
      enddo
c
      return
      end


c=======================================================================
c/////////////////////////  SUBROUTINE SET1D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine set1d(dest, ddim1, value)
c
c  SETS A FIELD TO A CONSTANT VALUE
c
c  written by: Greg Bryan
c  date:       October, 1995
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     ddim1      - destination dimension
c     value      - constant
c
c  OUTPUT ARGUMENTS: 
c     dest       - destination field
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
      INTG_PREC ddim1, value
      R_PREC    dest(ddim1)
c
c  locals
c
      INTG_PREC i
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      do i = 1, ddim1
         dest(i) = value
      enddo
c
      return
      end


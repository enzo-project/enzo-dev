#include "fortran.def"
c=======================================================================
c/////////////////////////  SUBROUTINE SMOOTH  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine smooth(source1, source2, source3, ndim, 
     &                  sdim1, sdim2, sdim3, nsmooth)
c
c  SMOOTH FIELD
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source       - source field
c     sdim1-3      - source dimension
c     ndim         - rank of fields
c
c  OUTPUT ARGUMENTS: 
c     source       - 
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
      INTG_PREC sdim1, sdim2, sdim3, ndim, nsmooth
      R_PREC    source1(sdim1, sdim2, sdim3), 
     &        source2(sdim1, sdim2, sdim3), 
     &        source3(sdim1, sdim2, sdim3)
c
c  locals
c
      INTG_PREC i, j, k, n, i1, j1, k1
      R_PREC    temp(MAX_ANY_SINGLE_DIRECTION), fact
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      fact = 1._RKIND/(2**nsmooth + 1._RKIND)
c
c     1) x-direction
c
      if (ndim .ge. 1) then
         do k=1, sdim3
            do j=1, sdim2
               do i=1, sdim1
                  temp(i) = 0._RKIND
               enddo
               do n=-nsmooth,+nsmooth
                  do i=1, sdim1
                     i1 = min(max(i+n, 2), sdim1)
                     temp(i) = temp(i) + source1(i1,j,k)
                  enddo
               enddo
               do i=2, sdim1
                  source1(i,j,k) = temp(i)*fact
               enddo
            enddo
         enddo
      endif
c
c     2) y-direction
c
      if (ndim .ge. 2) then
         do k=1, sdim3
            do i=1, sdim1
               do j=1, sdim2
                  temp(j) = 0._RKIND
               enddo
               do n=-nsmooth,+nsmooth
                  do j=1, sdim2
                     j1 = min(max(j+n, 2), sdim2)
                     temp(j) = temp(j) + source2(i,j1,k)
                  enddo
               enddo
               do j=2, sdim2
                  source2(i,j,k) = temp(j)*fact
               enddo
            enddo
         enddo
      endif
c
c     3) z-direction
c
      if (ndim .ge. 3) then
         do j=1, sdim2
            do i=1, sdim1
               do k=1, sdim3
                  temp(k) = 0._RKIND
               enddo
               do n=-nsmooth,+nsmooth
                  do k=1, sdim3
                     k1 = min(max(k+n, 2), sdim3)
                     temp(k) = temp(k) + source3(i,j,k1)
                  enddo
               enddo
               do k=2, sdim3
                  source3(i,j,k) = temp(k)*fact
               enddo
            enddo
         enddo
      endif
c
      return
      end


c=======================================================================
c/////////////////////////  SUBROUTINE SMOOTH2  \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine smooth2(source, dest, ndim, sdim1, sdim2, sdim3)
c
c  SMOOTH FIELD
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source       - source field
c     dest         - destination field
c     sdim1-3      - source dimension
c     ndim         - rank of fields
c
c  OUTPUT ARGUMENTS: 
c     source       - 
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
      INTG_PREC sdim1, sdim2, sdim3, ndim
      R_PREC    source(sdim1, sdim2, sdim3), dest(sdim1, sdim2, sdim3)
c
c  locals
c
      INTG_PREC i, j, k
      R_PREC    fact
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     1) 1d
c
      if (ndim .eq. 1) then
         fact = 1._RKIND/3._RKIND
         do i=2, sdim1-1
            dest(i,1,1) = fact*(source(i-1,1,1) + source(i,1,1) +
     &                          source(i+1,1,1))
         enddo
         dest(    1,1,1) = source(    1,1,1)
         dest(sdim1,1,1) = source(sdim1,1,1)
      endif
c
c     2) 2d
c
      if (ndim .eq. 2) then
         fact = 1._RKIND/5._RKIND
         do j=2, sdim2-1
            do i=2, sdim1-1
               dest(i,j,1) = fact*(source(i-1,j,1) + source(i,j,1) +
     &                             source(i+1,j,1) + source(i,j-1,1) +
     &                                               source(i,j+1,1))
            enddo
            dest(    1,j,1) = source(    1,j,1)
            dest(sdim1,j,1) = source(sdim1,j,1)
         enddo
         do i=1, sdim1
            dest(i,    1,1) = source(i,    1,1)
            dest(i,sdim2,1) = source(i,sdim2,1)
         enddo
      endif
c
c     3) 3d
c
      if (ndim .eq. 3) then
         fact = 1._RKIND/7._RKIND
         do k=2, sdim3-1
            do j=2, sdim2-1
               do i=2, sdim1-1
                  dest(i,j,k) = fact*(source(i,j,k) +
     &                       source(i-1,j,k) + source(i+1,j,k) +
     &                       source(i,j-1,k) + source(i,j+1,k) +
     &                       source(i,j,k-1) + source(i,j,k+1) )
               enddo
               dest(    1,j,k) = source(    1,j,k)
               dest(sdim1,j,k) = source(sdim1,j,k)
            enddo
            do i=1, sdim1
               dest(i,    1,k) = source(i,    1,k)
               dest(i,sdim2,k) = source(i,sdim2,k)
            enddo
         enddo
         do j=1, sdim2
            do i=1, sdim1
               dest(i,j,    1) = source(i,j,    1)
               dest(i,j,sdim3) = source(i,j,sdim3)
            enddo
         enddo
      endif
c
      return
      end

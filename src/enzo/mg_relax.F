c=======================================================================
c//////////////////////////  SUBROUTINE MG_RELAX  \\\\\\\\\\\\\\\\\\\\\\
c
      subroutine mg_relax(solution, rhs, ndim, dim1, dim2, dim3)
c
c  MULTIGRID: RELAX SOLUTION WITH DIFFERENCED POISSON OPERATOR
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:  Oliver Hahn, 
c  date:       February, 2010
c
c  PURPOSE:
c
c  INPUTS:
c     solution     - solution field
c     rhs          - right hand side
c     dim1-3       - dimensions
c     ndim         - rank of fields
c
c  OUTPUT ARGUMENTS: 
c     solution     - solution field
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
      INTG_PREC ndim, dim1, dim2, dim3
      R_PREC    solution(dim1, dim2, dim3), rhs(dim1, dim2, dim3)
c
c  locals
c
      INTG_PREC i, j, k, ipass, istart, jstart, kstart
      R_PREC    h1, h2, h3, coef1, coef2, coef3
      
#if defined(GRAVITY_4S) && !defined(GRAVITY_6S)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 4th order 13-point Laplacian --
c

c
c     Precompute some things
c
      h1 = 1._RKIND/REAL(dim1-1,RKIND)
      if (ndim .ge. 2) h2 = h1/REAL(dim2-1,RKIND)
      if (ndim .ge. 3) h3 = h2/REAL(dim3-1,RKIND)
      coef1 = 12._RKIND/30._RKIND
      coef2 = 12._RKIND/60._RKIND
      coef3 = 12._RKIND/90._RKIND
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do ipass=1, 2
            do i=ipass+2, dim1-2, 2
               solution(i,1,1) = coef1*((
     &                  - 1._RKIND * (solution(i+2,1,1)+
     &                           solution(i-2,1,1)) +
     &                  +16._RKIND * (solution(i+1,1,1)+
     &                           solution(i-1,1,1)) )
     &                  /12._RKIND-h1*rhs(i,1,1))
            enddo
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         jstart = 1
         do ipass=1, 2
            istart = jstart
            do j=3, dim2-2
               do i=istart+2, dim1-2, 2
                 solution(i,j,1) = coef2*((
     &                     - 1._RKIND * (solution(i+2,j,1)+
     &                              solution(i-2,j,1)+
     &                              solution(i,j+2,1)+
     &                              solution(i,j-2,1)) +
     &                     +16._RKIND * (solution(i+1,j,1)+
     &                              solution(i-1,j,1)+
     &				    solution(i,j+1,1)+
     &                              solution(i,j-1,1)) )
     &                     /12._RKIND-h2*rhs(i,j,1))
               enddo
               istart = 3-istart
            enddo
            jstart = 3-jstart
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         kstart = 1
         do ipass=1, 2
            jstart = kstart
            do k=3, dim3-2
               istart = jstart
               do j=3, dim2-2
                  do i=istart+2, dim1-2, 2
                    solution(i,j,k) = coef3*((
     &                   - 1._RKIND * (solution(i+2,j,k)+
     &                            solution(i-2,j,k)+
     &                            solution(i,j+2,k)+
     &                            solution(i,j-2,k)+
     &                            solution(i,j,k+2)+
     &                            solution(i,j,k-2))
     &                   +16._RKIND * (solution(i+1,j,k)+
     &                            solution(i-1,j,k)+
     &				  solution(i,j+1,k)+
     &                            solution(i,j-1,k)+
     &                            solution(i,j,k+1)+
     &                            solution(i,j,k-1)) )
     &                   /12._RKIND-h3*rhs(i,j,k))
                  enddo
                  istart = 3-istart
               enddo
               jstart = 3-jstart
            enddo
            kstart = 3-kstart
         enddo
      endif


#elif defined(GRAVITY_6S)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 6th order 19-point Laplacian --
c

c
c     Precompute some things
c
      h1 = 1._RKIND/REAL(dim1-1,RKIND)
      if (ndim .ge. 2) h2 = h1/REAL(dim2-1,RKIND)
      if (ndim .ge. 3) h3 = h2/REAL(dim3-1,RKIND)
      coef1 = 1080._RKIND/2720._RKIND
      coef2 = 1080._RKIND/5440._RKIND
      coef3 = 1080._RKIND/8160._RKIND
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do ipass=1, 2
            do i=ipass+3, dim1-3, 2
               solution(i,1,1) = coef1*((
     &                    1._RKIND * (solution(i+3,1,1)+
     &                           solution(i-3,1,1)) +
     &                  -96._RKIND * (solution(i+2,1,1)+
     &                           solution(i-2,1,1)) +
     &                +1455._RKIND * (solution(i+1,1,1)+
     &                           solution(i-1,1,1)) )
     &                /1080._RKIND-h1*rhs(i,1,1))
            enddo
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         jstart = 1
         do ipass=1, 2
            istart = jstart
            do j=4, dim2-3
               do i=istart+3, dim1-3, 2
                 solution(i,j,1) = coef2*((
     &                       1._RKIND * (solution(i+3,j,1)+
     &                              solution(i-3,j,1)+
     &                              solution(i,j+3,1)+
     &                              solution(i,j-3,1))
     &                     -96._RKIND * (solution(i+2,j,1)+
     &                              solution(i-2,j,1)+
     &                              solution(i,j+2,1)+
     &                              solution(i,j-2,1))
     &                   +1455._RKIND * (solution(i+1,j,1)+
     &                              solution(i-1,j,1)+
     &		                    solution(i,j+1,1)+
     &                              solution(i,j-1,1)))
     &                   /1080._RKIND-h2*rhs(i,j,1))
               enddo
               istart = 3-istart
            enddo
            jstart = 3-jstart
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         kstart = 1
         do ipass=1, 2
            jstart = kstart
            do k=4, dim3-3
               istart = jstart
               do j=4, dim2-3
                  do i=istart+3, dim1-3, 2
                    solution(i,j,k) = coef3*((
     &                     1._RKIND * (solution(i+3,j,k)+
     &                            solution(i-3,j,k)+
     &                            solution(i,j+3,k)+
     &                            solution(i,j-3,k)+
     &                            solution(i,j,k+3)+
     &                            solution(i,j,k-3))
     &                   -96._RKIND * (solution(i+2,j,k)+
     &                            solution(i-2,j,k)+
     &                            solution(i,j+2,k)+
     &                            solution(i,j-2,k)+
     &                            solution(i,j,k+2)+
     &                            solution(i,j,k-2))
     &                 +1455._RKIND * (solution(i+1,j,k)+
     &                            solution(i-1,j,k)+
     &				  solution(i,j+1,k)+
     &                            solution(i,j-1,k)+
     &                            solution(i,j,k+1)+
     &                            solution(i,j,k-1)) )
     &                 /1080._RKIND-h3*rhs(i,j,k))
                  enddo
                  istart = 3-istart
               enddo
               jstart = 3-jstart
            enddo
            kstart = 3-kstart
         enddo
      endif


#else
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 2nd order 7-point Laplacian --
c

c
c     Precompute some things
c
      h1 = 1._RKIND/REAL(dim1-1,RKIND)
      if (ndim .ge. 2) h2 = h1/REAL(dim2-1,RKIND)
      if (ndim .ge. 3) h3 = h2/REAL(dim3-1,RKIND)
      coef1 = 1._RKIND/2._RKIND
      coef2 = 1._RKIND/4._RKIND
      coef3 = 1._RKIND/6._RKIND
      
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do ipass=1, 2
            do i=ipass+1, dim1-1, 2
               solution(i,1,1) = coef1*(
     &                           solution(i-1,1,1)+solution(i+1,1,1) -
     &                           h1*rhs(i,1,1))
            enddo
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         jstart = 1
         do ipass=1, 2
            istart = jstart
            do j=2, dim2-1
               do i=istart+1, dim1-1, 2
                  solution(i,j,1) = coef2*(
     &                    solution(i-1,j  ,1) + solution(i+1,j  ,1) +
     &                    solution(i  ,j-1,1) + solution(i  ,j+1,1) -
     &                    h2*rhs(i,j,1))
               enddo
               istart = 3-istart
            enddo
            jstart = 3-jstart
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         kstart = 1
#define CACHE_AWARE
#ifdef CACHE_AWARE
! Cache aware Gauss-Seidel method: Data is passed through cache only once.         
         ! Red cells (even-even, odd-odd) on the first slab
         jstart = kstart
         istart = jstart
         do j = 2, dim2-1
            do i = istart+1, dim1-1, 2
               solution(i,j,2) = coef3*(
     &              solution(i-1,j  ,2) + solution(i+1,j  ,2) +
     &              solution(i  ,j-1,2) + solution(i  ,j+1,2) +
     &              solution(i  ,j  ,1) + solution(i  ,j  ,3) -
     &              h3*rhs(i,j,2))
            enddo
            istart = 3-istart
         enddo
         jstart = 3-jstart

         ! Fused loop with red/black cells
         do k = 3, dim3-1
            istart = jstart
            do j = 2, dim2-1
               do i = istart+1, dim1-1, 2
                  ! Red on this slab
                  solution(i,j,k) = coef3*(
     &                 solution(i-1,j  ,k  ) + solution(i+1,j  ,k  ) +
     &                 solution(i  ,j-1,k  ) + solution(i  ,j+1,k  ) +
     &                 solution(i  ,j  ,k-1) + solution(i  ,j  ,k+1) -
     &                 h3*rhs(i,j,k))
                  ! Black on previous slab
                  solution(i,j,k-1) = coef3*(
     &                 solution(i-1,j,k-1) + solution(i+1,j,k-1) +
     &                 solution(i,j-1,k-1) + solution(i,j+1,k-1) +
     &                 solution(i,j  ,k-2) + solution(i,j  ,k) -
     &                 h3*rhs(i,j,k-1))
               enddo
               istart = 3-istart
            enddo
            jstart = 3-jstart
         enddo

         ! Black cells (even-odd, odd-even) on the last slab
         istart = jstart
         do j = 2, dim2-1
            do i = istart+1, dim1-1, 2
               solution(i,j,dim3-1) = coef3*(
     &            solution(i-1,j,dim3-1) + solution(i+1,j,dim3-1)+
     &            solution(i,j-1,dim3-1) + solution(i,j+1,dim3-1)+
     &            solution(i,j  ,dim3-2) + solution(i,j  ,dim3) -
     &            h3*rhs(i,j,dim3-1))
            enddo
            istart = 3-istart
         enddo
#else
! Usual method: data is passed through cache twice
         do ipass=1, 2
            jstart = kstart
            do k=2, dim3-1
               istart = jstart
               do j=2, dim2-1
                  do i=istart+1, dim1-1, 2
                     solution(i,j,k) = coef3*(
     &                  solution(i-1,j  ,k  ) + solution(i+1,j  ,k  ) +
     &                  solution(i  ,j-1,k  ) + solution(i  ,j+1,k  ) +
     &                  solution(i  ,j  ,k-1) + solution(i  ,j  ,k+1) -
     &                  h3*rhs(i,j,k))
                  enddo
                  istart = 3-istart
               enddo
               jstart = 3-jstart
            enddo
            kstart = 3-kstart
         enddo
#endif /* CACHE_AWARE */
      endif
#endif
c
      return
      end

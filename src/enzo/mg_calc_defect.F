c=======================================================================
c////////////////////////  SUBROUTINE MG_CALC_DEFECT  \\\\\\\\\\\\\\\\\\
c
      subroutine mg_calc_defect(solution, rhs, defect, ndim, 
     &                          dim1, dim2, dim3, norm)
c
c  MULTIGRID: CALCULATE (NEGATIVE) DEFECT: -(Lu - f)
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
c     defect       - negative defect = -(Lu - f)
c     norm         - L2 norm of defect
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
      R_PREC    solution(dim1, dim2, dim3), rhs(dim1, dim2, dim3),
     &        defect(dim1, dim2, dim3), norm
      REAL*8  sum
c
c  locals
c
      INTG_PREC i, j, k
      R_PREC    h1, h2, h3
      
#if defined(GRAVITY_4S) && !defined(GRAVITY_6S)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 4th order 13-point Laplacian --
c

c
c     Precompute some things
c
      h1 = -REAL(dim1-1,RKIND)/12._RKIND
      h2 = h1*REAL(dim2-1,RKIND)
      h3 = h2*REAL(dim3-1,RKIND)
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=3, dim1-2
            defect(i,1,1) = h1*(-1._RKIND*(solution(i+2,1,1)+
     &                                solution(i-2,1,1))+
     &                          16._RKIND*(solution(i+1,1,1)+
     &                                solution(i-1,1,1))-
     &                          30._RKIND* solution(i,1,1)) + rhs(i,1,1)
         enddo
         defect(1,1,1) = 0._RKIND
	 	 defect(2,1,1) = 0._RKIND
         defect(dim1,1,1) = 0._RKIND
	 	 defect(dim1-1,1,1)=0._RKIND
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=3, dim2-2
            do i=3, dim1-2
          	  defect(i,j,1) = h2*
     &		     	   (-1._RKIND*(solution(i+2,j,1)+
     &                            solution(i-2,j,1)+
     &                            solution(i,j+2,1)+
     &                            solution(i,j-2,1))+
     &                      16._RKIND*(solution(i+1,j,1)+
     &                            solution(i-1,j,1)+
     &                            solution(i,j+1,1)+
     &                            solution(i,j-1,1))-
     &                      60._RKIND* solution(i,j,1)) + rhs(i,j,1)
      
            enddo
            defect(1,j,1) = 0._RKIND
	    	defect(2,j,1) = 0._RKIND
            defect(dim1,j,1) = 0._RKIND
	    	defect(dim1-1,j,1) = 0._RKIND
         enddo
         do i=1,dim1
            defect(i,1,1) = 0._RKIND
	    	defect(i,2,1) = 0._RKIND
            defect(i,dim2,1) = 0._RKIND
            defect(i,dim2-1,1) = 0._RKIND
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         do k=3, dim3-2
            do j=3, dim2-2
               do i=3, dim1-2
            	 defect(i,j,k) = h3*
     &			   (-1._RKIND*(solution(i+2,j,k)+
     &                            solution(i-2,j,k)+
     &                            solution(i,j+2,k)+
     &                            solution(i,j-2,k)+
     &			          solution(i,j,k+2)+
     &                            solution(i,j,k-2))+
     &                      16._RKIND*(solution(i+1,j,k)+ 
     &                            solution(i-1,j,k)+
     &                            solution(i,j+1,k)+
     &                            solution(i,j-1,k)+
     &		                  solution(i,j,k+1)+
     &                            solution(i,j,k-1))-
     &                      90._RKIND* solution(i,j,k)) + rhs(i,j,k)
        
               enddo
               defect(1,j,k) = 0._RKIND
	       defect(2,j,k) = 0._RKIND
               defect(dim1,j,k) = 0._RKIND
               defect(dim1-1,j,k) = 0._RKIND
            enddo
            do i=1, dim1
               defect(i,1,k) = 0._RKIND
               defect(i,2,k) = 0._RKIND
               defect(i,dim2,k) = 0._RKIND
               defect(i,dim2-1,k) = 0._RKIND
            enddo
         enddo
         do j=1, dim2
            do i=1, dim1
               defect(i,j,1) = 0._RKIND
               defect(i,j,2) = 0._RKIND
               defect(i,j,dim3) = 0._RKIND
               defect(i,j,dim3-1) = 0._RKIND
            enddo
         enddo
      endif
c
c     Compute norm
c
      sum = 0._RKIND
      do k=1, dim3
         do j=1, dim2
            do i=1, dim1
               sum = sum + defect(i,j,k)**2
            enddo
         enddo
      enddo
      norm = sqrt(sum)/(dim1*dim2*dim3)
      
#elif defined(GRAVITY_6S)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 6th order 19-point Laplacian --
c

c
c     Precompute some things
c
      h1 = -REAL(dim1-1,RKIND)/1080._RKIND
      h2 = h1*REAL(dim2-1,RKIND)
      h3 = h2*REAL(dim3-1,RKIND)
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=4, dim1-3
            defect(i,1,1) = h1*( 1._RKIND*(solution(i+3,1,1)+
     &                                solution(i-3,1,1))
     &                         -96._RKIND*(solution(i+2,1,1)+
     &                                solution(i-2,1,1))
     &                       +1455._RKIND*(solution(i+1,1,1)+
     &                                solution(i-1,1,1))
     &                       -2720._RKIND* solution(i,1,1)) + rhs(i,1,1)
         enddo
         defect(1,1,1) = 0._RKIND
         defect(2,1,1) = 0._RKIND
         defect(3,1,1) = 0._RKIND
         defect(dim1,1,1) = 0._RKIND
         defect(dim1-1,1,1) = 0._RKIND
         defect(dim1-2,1,1) = 0._RKIND
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=4, dim2-3
            do i=4, dim1-3
          	  defect(i,j,1) = h2*(
     &		    	     1._RKIND*(solution(i+3,j,1)+
     &                            solution(i-3,j,1)+
     &                            solution(i,j+3,1)+
     &                            solution(i,j-3,1))
     &                     -96._RKIND*(solution(i+2,j,1)+
     &                            solution(i-2,j,1)+
     &                            solution(i,j+2,1)+
     &                            solution(i,j-2,1))
     &                   +1455._RKIND*(solution(i+1,j,1)+
     &                            solution(i-1,j,1)+
     &                            solution(i,j+1,1)+
     &                            solution(i,j-1,1))
     &                   -5440._RKIND* solution(i,j,1)) + rhs(i,j,1)
      
            enddo
            defect(1,j,1) = 0._RKIND
            defect(2,j,1) = 0._RKIND
            defect(3,j,1) = 0._RKIND
            defect(dim1,j,1) = 0._RKIND
            defect(dim1-1,j,1) = 0._RKIND
            defect(dim1-2,j,1) = 0._RKIND
         enddo
         do i=1,dim1
            defect(i,1,1) = 0._RKIND
            defect(i,2,1) = 0._RKIND
            defect(i,3,1) = 0._RKIND
            defect(i,dim2,1) = 0._RKIND
            defect(i,dim2-1,1) = 0._RKIND
            defect(i,dim2-2,1) = 0._RKIND
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         do k=4, dim3-3
            do j=4, dim2-3
               do i=4, dim1-3
                  defect(i,j,k) = h3*(
     &		             1._RKIND*(solution(i+3,j,k)+
     &                            solution(i-3,j,k)+
     &                            solution(i,j+3,k)+
     &                            solution(i,j-3,k)+
     &		                  solution(i,j,k+3)+
     &                            solution(i,j,k-3))
     &                     -96._RKIND*(solution(i+2,j,k)+ 
     &                            solution(i-2,j,k)+
     &                            solution(i,j+2,k)+
     &                            solution(i,j-2,k)+
     &		                  solution(i,j,k+2)+
     &                            solution(i,j,k-2))
     &                   +1455._RKIND*(solution(i+1,j,k)+ 
     &                            solution(i-1,j,k)+
     &                            solution(i,j+1,k)+
     &                            solution(i,j-1,k)+
     &		                  solution(i,j,k+1)+
     &                            solution(i,j,k-1))
     &                   -8160._RKIND* solution(i,j,k)) + rhs(i,j,k)
  
               enddo
               defect(1,j,k) = 0._RKIND
               defect(2,j,k) = 0._RKIND
               defect(3,j,k) = 0._RKIND
               defect(dim1,j,k) = 0._RKIND
               defect(dim1-1,j,k) = 0._RKIND
               defect(dim1-2,j,k) = 0._RKIND
            enddo
            do i=1, dim1
               defect(i,1,k) = 0._RKIND
               defect(i,2,k) = 0._RKIND
               defect(i,3,k) = 0._RKIND
               defect(i,dim2,k) = 0._RKIND
               defect(i,dim2-1,k) = 0._RKIND
               defect(i,dim2-2,k) = 0._RKIND
            enddo
         enddo
         do j=1, dim2
            do i=1, dim1
               defect(i,j,1) = 0._RKIND
               defect(i,j,2) = 0._RKIND
               defect(i,j,3) = 0._RKIND
               defect(i,j,dim3) = 0._RKIND
               defect(i,j,dim3-1) = 0._RKIND
               defect(i,j,dim3-2) = 0._RKIND
            enddo
         enddo
      endif
c
c     Compute norm
c
      sum = 0._RKIND
      do k=1, dim3
         do j=1, dim2
            do i=1, dim1
               sum = sum + defect(i,j,k)**2
            enddo
         enddo
      enddo
      norm = sqrt(sum)/(dim1*dim2*dim3)
c

#else
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 2nd order 7-point Laplacian --
c


c
c     Precompute some things
c
      h1 = -REAL(dim1-1,RKIND)
      h2 = h1*REAL(dim2-1,RKIND)
      h3 = h2*REAL(dim3-1,RKIND)
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=2, dim1-1
            defect(i,1,1) = h1*(solution(i-1,1,1) + solution(i+1,1,1) -
     &                           2._RKIND*solution(i,1,1)) + rhs(i,1,1)
         enddo
         defect(1,1,1) = 0._RKIND
         defect(dim1,1,1) = 0._RKIND
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=2, dim2-1
            do i=2, dim1-1
               defect(i,j,1) = h2*(
     &                    solution(i-1,j  ,1) + solution(i+1,j  ,1) +
     &                    solution(i  ,j-1,1) + solution(i  ,j+1,1) -
     &                    4._RKIND*solution(i,j,1)) + rhs(i,j,1)
            enddo
            defect(1,j,1) = 0._RKIND
            defect(dim1,j,1) = 0._RKIND
         enddo
         do i=1,dim1
            defect(i,1,1) = 0._RKIND
            defect(i,dim2,1) = 0._RKIND
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         do k=2, dim3-1
            do j=2, dim2-1
               do i=2, dim1-1
                  defect(i,j,k) = h3*(
     &                 solution(i-1,j  ,k  ) + solution(i+1,j  ,k  ) +
     &                 solution(i  ,j-1,k  ) + solution(i  ,j+1,k  ) +
     &                 solution(i  ,j  ,k-1) + solution(i  ,j  ,k+1) -
     &                 6._RKIND*solution(i,j,k)) + rhs(i,j,k)
               enddo
               defect(1,j,k) = 0._RKIND
               defect(dim1,j,k) = 0._RKIND
            enddo
            do i=1, dim1
               defect(i,1,k) = 0._RKIND
               defect(i,dim2,k) = 0._RKIND
            enddo
         enddo
         do j=1, dim2
            do i=1, dim1
               defect(i,j,1) = 0._RKIND
               defect(i,j,dim3) = 0._RKIND
            enddo
         enddo
      endif
c
c     Compute norm
c
      sum = 0._RKIND
      do k=1, dim3
         do j=1, dim2
            do i=1, dim1
               sum = sum + defect(i,j,k)**2
            enddo
         enddo
      enddo
      norm = sqrt(sum)/(dim1*dim2*dim3)
c

#endif

      return
      end

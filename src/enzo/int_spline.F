#include "fortran.def"
#include "error.def"

c=======================================================================
c/////////////////////////  SUBROUTINE INT_SPLINE \\\\\\\\\\\\\\\\\\\\\\
c
      subroutine int_spline(source, dest, ndim, sdim1, sdim2, sdim3,
     &                   ddim1, ddim2, ddim3, start1, start2, start3,
     &                   refine1, refine2, refine3, temp1, temp2)
c
c  SPLINE INTERPOLATE FROM SOURCE TO DEST
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
c     ddim1-3      - destination dimension
c     ndim         - rank of fields
c     start1-3     - dest start index in destination cells
c     refine1-3    - refinement factors
c     temp1,2      - temporary space
c
c  OUTPUT ARGUMENTS: 
c     dest         - prolonged field
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
      INTG_PREC ijkn
      parameter (ijkn = MAX_ANY_SINGLE_DIRECTION)
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim,
     &        start1, start2, start3, refine1, refine2, refine3
      R_PREC    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3),
     &        temp1(ddim1, ddim2/refine2+2), 
     &        temp2(ddim1, ddim2, ddim3/refine3+2)
c
c  locals
c
      INTG_PREC i, j, k, i1, j1, k1, i2, j2, k2
      R_PREC work1(ijkn), work2(ijkn), work3(ijkn), work4(ijkn)
c     
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     Precompute some things
c
      i1 = int(start1/refine1,IKIND)
      i2 = int((start1+ddim1)/refine1,IKIND)+1
      j1 = int(start2/refine2,IKIND)
      j2 = int((start2+ddim2)/refine2,IKIND)+1
      k1 = int(start3/refine3,IKIND)
      k2 = int((start3+ddim3)/refine3,IKIND)+1
c
      if ((i2-i1-1)*refine1 .ne. ddim1) then
         write(6,*) i2-i1-1, refine1, ddim1
         ERROR_MESSAGE
      endif
c
c     a) 1D
c
      if (ndim .eq. 1) then
         call spline_make(source(i1,1,1), i2-i1+1_IKIND, work1, work4)
         call spline_eval(source(i1,1,1), i2-i1+1_IKIND, work1,
     &                    refine1, dest(1,1,1))
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=j1, j2
            call spline_make(source(i1,j,1), i2-i1+1_IKIND, work1,work4)
            call spline_eval(source(i1,j,1), i2-i1+1_IKIND, work1,
     &                       refine1, temp1(1,j-j1+1))
         enddo
c
         do i=1, ddim1
            do j=1, j2-j1+1
               work1(j) = temp1(i,j)
            enddo
            call spline_make(work1, j2-j1+1_IKIND, work2, work4)
            call spline_eval(work1, j2-j1+1_IKIND, work2, refine2,work3)
     &                          
            do j=1, ddim2
               dest(i,j,1) = work3(j)
            enddo
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
c
c        Loop over k-slices and generate splined slices
c
         do k=k1, k2
            do j=j1, j2
               call spline_make(source(i1,j,k), i2-i1+1_IKIND, work1, 
     &              work4)
               call spline_eval(source(i1,j,k), i2-i1+1_IKIND, work1,
     &                           refine1, temp1(1,j-j1+1))
            enddo
c
            do i=1, ddim1
               do j=1, j2-j1+1
                  work1(j) = temp1(i,j)
               enddo
               call spline_make(work1, j2-j1+1_IKIND, work2, work4)
               call spline_eval(work1, j2-j1+1_IKIND, work2, refine2, 
     &              work3)
     &                          
               do j=1, ddim2
                  temp2(i,j,k-k1+1) = work3(j)
               enddo
            enddo
         enddo
c
         do j=1, ddim2
            do i=1, ddim1
               do k=1, k2-k1+1
                  work1(k) = temp2(i,j,k)
               enddo
               call spline_make(work1, k2-k1+1_IKIND, work2, work4)
               call spline_eval(work1, k2-k1+1_IKIND, work2, refine3, 
     &              work3)
               do k=1, ddim3
                  dest(i,j,k) = work3(k)
               enddo
            enddo
         enddo
c
      endif
c
      return
      end



      subroutine spline_eval(y, n, y2, refine, result)
c
      implicit none
#include "fortran_types.def"
c
      INTG_PREC i, j, n, refine
      R_PREC y(n), y2(n), result((n-2)*refine), a, b, step, x
c
      step = 1._RKIND/refine
      x = 1._RKIND - 0.5_RKIND*step
      do i=1, (n-2)*refine
         j = int(x,IKIND) + 1
         a = REAL(j,RKIND) - x
         b = 1._RKIND - a
         result(i) = a*y(j) + b*y(j+1) + ((a**3 - a)*y2(j) +
     &               (b**3-b)*y2(j+1))/6._RKIND
         x = x + step
c         write(7,1000) i,a*y(j) + b*y(j+1),
c     &        result(i),y(j),y(j+1),y2(j),y2(j+1)
      enddo
c 1000 format(i5,1p,10(1x,e12.3))
c
      return
      end



      subroutine spline_make(y, n, y2, work)
c
c     This routine assumes that all y points are deltax=1 apart
c
      implicit none
#include "fortran_types.def"
c
      INTG_PREC i, n
      R_PREC y(n), y2(n), work(n), p
c
c     hardwire for natural derivatives (zero 2nd deriv at ends)
c
      y2(1) = 0._RKIND
      work(1) = 0._RKIND
c
      do i=2, n-1
         p = 1._RKIND/(0.5_RKIND*y2(i-1) + 2._RKIND)
         y2(i) = -0.5_RKIND*p
         work(i) = (3._RKIND*((y(i+1)-y(i)) - (y(i)-y(i-1))) - 
     &              0.5_RKIND*work(i-1))*p
      enddo
c
      y2(n) = 0._RKIND
c
c     back substitution
c
      do i=n-1, 1, -1
         y2(i) = y2(i) * y2(i+1) + work(i)
      enddo
c
      return
      end

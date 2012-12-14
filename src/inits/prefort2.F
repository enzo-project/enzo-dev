c=======================================================================
c///////////////////////  SUBROUTINE PREFORT  \\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine prefort(buffer, tempbuffer, dim1, dim2, dim3,
     &                   rdim1, rdim2, rdim3, rank, isign)
c
c  PREPARES THE REAL/CMPLX_PREC BUFFER FOR A CALL TO FORTRAN CMPLX_PREC FFT
c
c     written by: Robert Harkness
c     date:       May, 2003
c     modified1:  
c
c  PURPOSE:
c
c  INPUTS:
c     buffer      - real or CMPLX_PREC array to be transformed
c     tempbuffer  - temporary buffer of size 2*dim1*dim2*dim3
c     dim1,2,3    - dimension of transform
c     rdim1,2,3   - declared dimension of buffer
c     rank        - number of dimensions to be used in transform
c     isign       - direction of transform
c
c  OUTPUTS:
c     buffer      - transformed array
c
c  LOCALS:
c
c  EXTERNALS:
c
c  LOCALS:
c-----------------------------------------------------------------------
      implicit NONE
#include "../enzo/fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  arguments
c
      INTG_PREC dim1, dim2, dim3, isign, rank, rdim1, rdim2, rdim3
      R_PREC    buffer(rdim1, rdim2, rdim3), 
     &        tempbuffer(2*dim1, dim2, dim3)
c
c  locals
c
      INTG_PREC i, i1, j, j1, k, k1, dims(3)
      R_PREC    factor
c
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c     Error check
c

      write(0,'("Pre-Fortran FFT ",8i8)')
     & rdim1, rdim2, rdim3, dim1, dim2, dim3, rank, isign

      if (rdim1 .lt. dim1+2) then
         write (0,*) 'PREFORT: rdim1 < dim1+2:', rdim1, dim1+2
         stop
      endif
c
c     Copy dims into INTG_PREC array for fourn
c
      dims(1) = dim1
      dims(2) = dim2
      dims(3) = dim3
c
c     1) Forward direction
c
      if (isign .eq. -1) then
c
c        Copy real data into CMPLX_PREC temporary
c
         do k = 1, dim3
            do j = 1, dim2
               do i = 1, dim1
                  tempbuffer(i*2-1, j, k) = buffer(i, j, k)
                  tempbuffer(i*2  , j, k) = 0._RKIND
               enddo
            enddo
         enddo
c
c        Carry out the transform
c

!        call nr_fftx(tempbuffer, rank, dim1, dim2, dim3, isign)
!        call nrot_fftx(tempbuffer, rank, dim1, dim2, dim3, isign)
!        call sing90_fftx(tempbuffer, rank, dim1, dim2, dim3, isign)
!        call sing66_fftx(tempbuffer, rank, dim1, dim2, dim3, isign)

         call fortfft(tempbuffer, rank, dim1, dim2, dim3, isign)

c
c        Copy unique part of CMPLX_PREC temporary back to buffer.
c
         do k = 1, dim3
            do j = 1, dim2
               do i = 1, dim1+2
                  buffer(i, j, k) = tempbuffer(i, j, k)
               enddo
            enddo
         enddo
c
c     2) Inverse transform
c
      else
c
c        Recreate the other half of the CMPLX_PREC data.
c
         do k = 1, dim3
            do j = 1, dim2
c
c              First copy the 1/2 that is already there.
c
               do i = 1, dim1+2
                  tempbuffer(i, j, k) = buffer(i, j, k)
               enddo
c
c              Use conjugate relations to make up other 1/2.
c
               k1 = dim3 + 2 - k
               if (k .eq. 1) k1 = 1
               j1 = dim2 + 2 - j
               if (j .eq. 1) j1 = 1
c
               do i = dim1/2+2, dim1
                  i1 = dim1 + 2 - i
                  tempbuffer(i*2-1, j, k) =   buffer(i1*2-1, j1, k1)
                  tempbuffer(i*2  , j, k) = - buffer(i1*2  , j1, k1)
               enddo
c
            enddo
         enddo
c
c        Perform transform
c

!        call nr_fftx(tempbuffer, rank, dim1, dim2, dim3, isign)
!        call nrot_fftx(tempbuffer, rank, dim1, dim2, dim3, isign)
!        call sing66_fftx(tempbuffer, rank, dim1, dim2, dim3, isign)
!        call sing90_fftx(tempbuffer, rank, dim1, dim2, dim3, isign)

         call fortfft(tempbuffer, rank, dim1, dim2, dim3, isign)

c
c        Copy real part of tempbuffer to buffer & scale.
c

!        factor = 1.0/REAL(dim1*dim2*dim3,RKIND)

         do k = 1, dim3
            do j = 1, dim2
               do i = 1, dim1
                  buffer(i, j, k) = tempbuffer(i*2-1, j, k) ! *factor
               enddo
            enddo
         enddo
c
      endif
c
c     Done
c
      write(0,'("Exit Pre-Fortran")')

      return
      end

#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE FOURN  \\\\\\\\\\\\\\\\\\\\\\\\\\\
c
        subroutine fourn(data, nn, ndim, isign)
c
c  PERFORM A CMPLX_PREC-CMPLX_PREC N-DIMENSIONAL FFT
c
c   written by: Numerical Recipes
c   data:
c   modified1:
c
c  PURPOSE:  Replaces DATA by its NDIM-dimensional discrete Fourier transform
c      if ISIGN is 1.  NN is an INTG_PREC array of length NDIM, containing
c      the lengths of each dimension (number of CMPLX_PREC values), which
c      MUST all be powers of 2.  DATA is a real array of length twice the
c      product of these lengths, in which the data are stored as in a
c      multidimensional CMPLX_PREC FORTRAN array.  If ISIGN in input as -1,
c      DATA is replaced by its inverse transform times the product of the
c      lengths of all dimensions.
c
c  INPUTS:
c    nn    - an INTG_PREC array of length NDIM, containing the lengths of 
c          each dimension
c    data  - data to be transformed
c    ndim  - number of dimensions
c    isign - +1 is forward, -1 is inverse
c
c  OUTPUTS:
c    data  - transformed array
c
c-----------------------------------------------------------------------
      implicit NONE
#include "fortran_types.def"
c
c  Arguments
c
      R_PREC    data(*)
      INTG_PREC isign, ndim, nn(ndim)
c
c  External functions
c
      integer*8 power_of_2
c
c  Locals
c
      INTG_PREC i1, i2, i3, i2rev, i3rev, ibit, idim,
     &        ifp1, ifp2, ip1, ip2, ip3, k2, n, nprev, nrem, ntot
      INTG_PREC ibad
      R_PREC    tempi, tempr
c     real*8 twopi, theta, wr, wi, wpr, wpi, wtemp
      R_PREC    twopi, theta, wr, wi, wpr, wpi, wtemp
      parameter (twopi=6.2831853071796d0)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c     Compute total number of CMPLX_PREC values
c

!     write(*,*) 'FOURN:', nn, ndim, isign

      ibad = 0

      do idim = 1, ndim
      if( power_of_2(nn(idim)) .ne. 0 ) then
        write(*,'("Index ",i1," is not a power of 2")') idim
        ibad = 1
      end if
      end do

      if( ibad .ne. 0 ) then
         ERROR_MESSAGE
      end if

      ntot = 1
      do idim=1, ndim
         ntot = ntot*nn(idim)
      enddo
      nprev=1
c
c     Main loop over the dimensions.
c
      do idim=1, ndim
         n     = nn(idim)
         nrem  = ntot/(n*nprev)
         ip1   = 2*nprev
         ip2   = ip1*n
         ip3   = ip2*nrem
         i2rev = 1
c
c        This is the bit reversal section of the routine.
c
         do i2=1,ip2,ip1
            if (i2.lt.i2rev) then
               do i1=i2, i2+ip1-2, 2
                  do i3=i1, ip3, ip2
                     i3rev = i2rev + i3 - i2
                     tempr = data(i3  )
                     tempi = data(i3+1)
                     data(i3     ) = data(i3rev  )
                     data(i3+1   ) = data(i3rev+1)
                     data(i3rev  ) = tempr
                     data(i3rev+1) = tempi
                  enddo
               enddo
            endif
c
            ibit = ip2/2
 1          if ((ibit.ge.ip1) .and. (i2rev.gt.ibit)) then
               i2rev = i2rev-ibit
               ibit  = ibit/2
               go to 1
            end if
            i2rev = i2rev + ibit
         enddo
c
c        Here begins the Danielson-Lanczos section of the routine.
c
          ifp1 = ip1
 2        if (ifp1 .lt. ip2) then
             ifp2 = 2*ifp1
c
c         Initialize for the trigonometric recurrence.
c
             theta = isign*twopi/(ifp2/ip1)
             wpr   = -2._RKIND*sin(0.5_RKIND*theta)**2
             wpi   = sin(theta)
             wr    = 1._RKIND
             wi    = 0._RKIND
             do i3 = 1, ifp1, ip1
                do i1 = i3, i3+ip1-2, 2
                   do i2 = i1, ip3, ifp2
c
c                     Danielson-Lanczos formula.
c
                      k2 = i2 + ifp1
                      tempr = wr * data(k2  ) - wi*data(k2+1)
                      tempi = wr * data(k2+1) + wi*data(k2  )
                      data(k2  ) = data(i2  ) - tempr
                      data(k2+1) = data(i2+1) - tempi
                      data(i2  ) = data(i2  ) + tempr
                      data(i2+1) = data(i2+1) + tempi
                   enddo
                enddo
c
c              Trigonometric recurrence.
c
                wtemp = wr
                wr    = wr*wpr - wi*wpi    + wr
                wi    = wi*wpr + wtemp*wpi + wi
             enddo
c
             ifp1 = ifp2
             go to 2
          end if
          nprev = n*nprev
c
       enddo
c
       return
       end

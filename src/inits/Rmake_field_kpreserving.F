#include "../enzo/fortran.def"
!=======================================================================
!//////////////////////  SUBROUTINE MAKE_FIELD  \\\\\\\\\\\\\\\\\\\\\\\\

      subroutine make_field_kpreserving(field, nx, ny, nz, 
     &                      in, jn, kn, itype, iseed, box,
     &                      PSTable, PSMin, PSStep, kfcutoff)

!  COMPUTES RANDOM GAUSSIAN FIELD FROM SPECIFIED POWER SPECTRUM
!
!  written by: Greg Bryan
!  date:       June, 1997
!  modified:   Robert Harkness
!  date:       November, 2003
!
!  PURPOSE: 
!
!  INPUTS:
!        i,j,kn      = real dimensions of green
!        nx,ny,nz    = active dimensions of green
!        itype       = field type (0 - density, 1/2/3 - x/y/z displacement)
!        iseed       = random number seed (negative)
!        box         = size
!        PSTable     = Table of precomputed PS values
!        PSMin       = minimum x value in PSTable
!        PSStep      = x step in PSTable
!        kfcutoff    = high k filter (sharp) in units of the fundamental
!
!  Outputs:
!        field       = gaussian random field
!
!  LOCALS:
!        num_dim     = number of dimensions to be used for force law
!        nx,y,zmid   = midpoint (+1) of the grid in each axis
!        nx,y,zd2    = number of grid points divided by 2 for each axis

      implicit NONE
#include "../enzo/fortran_types.def"

!     Arguments

      INTG_PREC :: in, jn, kn, nx, ny, nz, nxmax, nymax, nzmax, 
     &           itype, iseed, kfcutoff
      R_PREC ::    field(in, jn, kn), box, 
     &           PSMin, PSPart, PSStep, PSTable(1)

!     Locals

      INTG_PREC :: i, i1, j, j1, n, n1
      R_PREC ::    dummy, twopi, kcutoffsq, dk
      CMPLX_PREC :: z

      INTG_PREC :: long_seed

!     External functions

      R_PREC ::    ran1

!     Set constants

      twopi  = 8.0_RKIND*atan(1.0_RKIND)
      dk     = twopi/box
      kcutoffsq = 1.0e30_RKIND

      if (kfcutoff .gt. 0) kcutoffsq = (kfcutoff*dk)**2

!     Initialize random # generator with random seed

      long_seed = iseed
      n = iseed
!     dummy = ran1(n)
      call enzo_seed(long_seed)

!     Loop over k-box sizes, so that we fill k-space from low-k outwards

      do n=1,nx/2

         do i=-n+1, n
            do j=-n+1, n

               i1 = mod(i+nx,nx)+1
               j1 = mod(j+nx,nx)+1
               n1 = mod(1_IKIND-n+nx,nx)+1

!              1) +i plane

               call processk(n,i,j, dk, PSMin, PSStep, PSTable, 
     &                       itype, z, kcutoffsq, box)

               field((n+1)*2-1,i1,j1) = REAL(z,RKIND)
               field((n+1)*2  ,i1,j1) = imag(z)

!              2) +j and -j plane
!                 (the i .ne. n is to avoid overlapping with (1))

               if (i .ge. 0 .and. i .ne. n) then

                  call processk(i,n,j, dk, PSMin, PSStep, PSTable, 
     &                          itype, z, kcutoffsq, box)

                  field(i1*2-1,n+1,j1) = REAL(z,RKIND)
                  field(i1*2  ,n+1,j1) = imag(z)

                  call processk(i,1_IKIND-n,j, dk, PSMin, PSStep, 
     &                          PSTable, itype, z, kcutoffsq, box)

                  field(i1*2-1,n1,j1) = REAL(z,RKIND)
                  field(i1*2  ,n1,j1) = imag(z)

               endif

!              3) +k and -k plane
!                 (the logic involving j is to avoid overlapping with (2))

               if (i .ge. 0 .and. i .ne. n .and. 
     &             j .ne. -n+1 .and. j .ne. n) then

                  call processk(i,j,n, dk, PSMin, PSStep, PSTable,
     &                          itype, z, kcutoffsq, box)

                  field(i1*2-1,j1,n+1) = REAL(z,RKIND)
                  field(i1*2  ,j1,n+1) = imag(z)

                  call processk(i,j,1_IKIND-n, dk, PSMin, PSStep, 
     &                          PSTable, itype, z, kcutoffsq, box)

                  field(i1*2-1,j1,n1) = REAL(z,RKIND)
                  field(i1*2  ,j1,n1) = imag(z)

               endif

            enddo
         enddo

      enddo

      do i=1, in
         do j=1, jn
            do n=1, kn
               field(i,j,n) = field(i,j,n) * REAL(nx*ny*nz,RKIND)
            enddo
         enddo
      enddo

!     Clear the zero wavenumber position

      field(1,1,1) = 0.0_RKIND
      field(2,1,1) = 0.0_RKIND

!     Adjust the field to satisfy the conjugate relations that
!     are implied by a zero imaginary part.

      call adjfft(field, nx, ny, nz, in, jn)

      return
      end


c===================================================================

      subroutine processk(i, j, k, dk, PSMin, PSStep, PSTable, 
     &                    itype, z, kcutoffsq, box)

      implicit none
#include "../enzo/fortran_types.def"

!     Parameter

      R_PREC, parameter :: twopi = 2._RKIND*3.14159265358979324_RKIND

!     Arguments

      INTG_PREC :: i, j, k, itype
      R_PREC ::    dk, PSMin, PSStep, PSTable(*)
      R_PREC ::    kcutoffsq, box
      CMPLX_PREC :: z

!     Locals

      INTG_PREC :: index
      R_PREC :: psval, kdir, klog, ang, amp, kmodsq
      R_PREC :: ranf_min

!     External function

      R_PREC :: ran1
      R_PREC :: enzo_ranf



!     Define table lookup function
 
      R_PREC ::    Table1, Table2, Step, Min, Tablex, TableLookUp
      INTG_PREC :: Tablei

      TableLookUp(Table1, Table2, Step, Min, Tablei, Tablex) = 
     &            Table1 + (Tablex - REAL(Tablei-1,RKIND)*Step - Min) 
     &            / Step * (Table2 - Table1)


      kmodsq = max(i**2 + j**2 + k**2, 1)*dk**2
      klog   = 0.5_RKIND*log(kmodsq)
      index = int((klog - PSMin)/PSStep)
      psval = TableLookUp(PSTable(index), PSTable(index+1),
     &                    PSStep, PSMin, index, klog)
      psval = psval * dk**3

      if (kmodsq .gt. kcutoffsq) psval = 0.0_RKIND

!     Generate a complex number with random phase and amplitude
!     Gaussian distributed with a mean of sqrt(psval) with the
!     Box-Muller method.  Note we have supressed a factor of
!     sqrt(2) since we must also divide by this factor to account
!     for the dreary fact that we are really generating two random
!     fields (if we were doing a complex-to-complex transform
!     this would show up when we discarded the perfectly
!     good imaginary component of the transformed field).  whew.

      ranf_min = 1.e-37_RKIND

      ang = twopi*enzo_ranf()
      amp = sqrt(-log(max(enzo_ranf(),ranf_min)) * psval)
      z   = cmplx(cos(ang), sin(ang)) * amp

!     Process this on the basis of itype:
!      0)   density field - just leave it be.
!      1-3) displacement field - multiply by vec(k)/k^2
!           (and then convert from Mpc to fraction of box).

      if (itype .ne. 0) then
         if (itype .eq. 1) kdir = REAL(i,RKIND)*dk
         if (itype .eq. 2) kdir = REAL(j,RKIND)*dk
         if (itype .eq. 3) kdir = REAL(k,RKIND)*dk
         z = z * cmplx(0._RKIND,1._RKIND) * kdir / (kmodsq * box)
      endif

      return
      end

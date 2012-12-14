#include "../enzo/fortran.def"
!=======================================================================
!//////////////////////  SUBROUTINE MAKE_FIELD  \\\\\\\\\\\\\\\\\\\\\\\\

      subroutine make_field(field, nx, ny, nz, nxmax, nymax, nzmax,
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
!        nx,y,zmax   = dimensions of k field (for random number countinf)
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
     &           PSMin, PSPart, PSStep, PSTable(*)

!     Locals

      INTG_PREC :: i, ii, j, jj, k, kk, index, nxmid, nymid, nzmid
      R_PREC ::    ang, amp, d3k, dummy, kmodsq, kx, ky, kz, kdir,
     &           klog, psval, twopi, kcutoffsq, dkx, dky, dkz
      CMPLX_PREC :: z

      R_PREC ::    ranf_min
      INTG_PREC :: long_seed

!     External function

      R_PREC ::    ran1, enzo_ranf

!  Define table lookup function

      R_PREC ::    Table1, Table2, Step, Min, Tablex, TableLookUp
      INTG_PREC :: Tablei

      TableLookUp(Table1, Table2, Step, Min, Tablei, Tablex) = 
     &            Table1 + (Tablex - REAL(Tablei-1,RKIND)*Step - Min) 
     &            / Step * (Table2 - Table1)


!     Set constants

      twopi  = 8.0_RKIND*atan(1.0_RKIND)
      dkx    = twopi/box
      dky    = twopi/(ny*box/nx)
      dkz    = twopi/(nz*box/nx)
      d3k    = (twopi/box)**3
      kcutoffsq = 1.0e30_RKIND

      if (kfcutoff .gt. 0) kcutoffsq = (kfcutoff*dkx)**2

!     Initialize random # generator with random seed

      long_seed = iseed

!     dummy = ran1(iseed)

      call enzo_seed(long_seed)

!     Set maximum wavenumbers in the positive direction

      nxmid  = max(nx/2_IKIND + 1_IKIND, 1_IKIND)
      nymid  = max(ny/2_IKIND + 1_IKIND, 1_IKIND)
      nzmid  = max(nz/2_IKIND + 1_IKIND, 1_IKIND)

!     Loop over k

      do k = 1, nzmax
         kk = k-1
         if (k .gt. nzmax/2+1) kk = kk - nzmax
         kz    = REAL(kk,RKIND)*dkz

!        Loop over j

         do j = 1, nymax
            jj = j-1
            if (j .gt. nymax/2+1) jj = jj - nymax
            ky    = REAL(jj,RKIND)*dky
 
!           If this j corresponds to a wavenumber vector that is not
!           inside field then just skip over it (& eat the random nums).

            if (jj .ge. nymid .or. jj .lt. -ny/2+1 .or.
     &          kk .ge. nzmid .or. kk .lt. -nz/2+1     ) then

               do i=1,nxmax/2+1
                  dummy = enzo_ranf()
                  dummy = enzo_ranf()
               enddo

               goto 100

            endif

!           Loop over i

            do i = 1, nxmid
               ii    = i-1
               kx    = REAL(ii,RKIND)*dkx

!              Compute kmod and lookup value in table
!              (and convert from power per mode).

               kmodsq  = (kx**2 + ky**2 + kz**2)
               if (kmodsq .eq. 0) kmodsq = 1.0_RKIND
               klog   = 0.5_RKIND*log(kmodsq)
               index = int((klog - PSMin)/PSStep,IKIND)
               psval = TableLookUp(PSTable(index), PSTable(index+1),
     &                             PSStep, PSMin, index, klog)
               psval = psval * d3k
               if (kmodsq .gt. kcutoffsq) psval = 0.0_RKIND

!              Generate a complex number with random phase and amplitude
!              Gaussian distributed with a mean of sqrt(psval) with the
!              Box-Muller method.  Note we have supressed a factor of
!              sqrt(2) since we must also divide by this factor to account
!              for the dreary fact that we are really generating two random
!              fields (if we were doing a complex-to-complex transform
!              this would show up when we discarded the perfectly
!              good imaginary component of the transformed field).  whew.

               ranf_min = 1.e-37_RKIND

               ang = twopi*enzo_ranf()
               amp = sqrt(-log(max(enzo_ranf(),ranf_min)) * psval)
               z   = cmplx(cos(ang), sin(ang)) * amp

!              Process this on the basis of itype:
!                0)   density field - just leave it be.
!                1-3) displacement field - multiply by vec(k)/k^2
!                     (and then convert from Mpc to fraction of box).

               if (itype .ne. 0) then
                  if (itype .eq. 1) kdir = kx
                  if (itype .eq. 2) kdir = ky
                  if (itype .eq. 3) kdir = kz
                  z = z * cmplx(0._RKIND,1._RKIND) * kdir
     &                 / (kmodsq**2 * box)
               endif

!              Set the complex field

               field(i*2-1,j,k) = REAL(z,RKIND)
               field(i*2  ,j,k) = IMAG(z)

            enddo

!           Now loop over the rest of the kx space to use up the
!           allotted number of random numbers

            do i=nxmid+1, nxmax/2+1
               dummy = enzo_ranf()
               dummy = enzo_ranf()
            enddo

 100        continue

         enddo
      enddo

!     Clear the zero wavenumber position

      field(1,1,1) = 0.0_RKIND
      field(2,1,1) = 0.0_RKIND

!     Adjust the field to satisfy that conjugate relations that
!     are implied by a zero imaginary part.

      call adjfft(field, nx, ny, nz, in, jn)

      return
      end

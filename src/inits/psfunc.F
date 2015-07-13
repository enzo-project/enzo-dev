#include "../enzo/fortran.def"
!=======================================================================
! Set R_PREC to match what's in fortran_types.def (since it cannot be 
! included before the first function declaration)
#ifdef CONFIG_BFLOAT_4
#define R_PREC real*4
#endif

#ifdef CONFIG_BFLOAT_8
#define R_PREC real*8
#endif
!=======================================================================


!=======================================================================
!//////////////////////////  FUNCTION PS_FUNC  \\\\\\\\\\\\\\\\\\\\\\\\\

      R_PREC function psfunc(k, ispectrum, omega0, hub, omega_nu,
     &                  omega_lam, psindex, omegab0, z, gamma,
     &                  psnorm, growth_factor, kcutoff)

!  written by: Greg Bryan
!  date:       May, 1996
!  modified:   Robert Harkness
!  date:       November, 2003
!
!  PURPOSE: Generalized power spectrum function
!
!  INPUT: 
!    k        - wavenumber (in Mpc^-1)
!    (also from cosmo_parm, including:
!    ipower_spectrum - indicates which function to use)
!
!  OUTPUT:
!    ps_func  - spectral density (Mpc^3)

      implicit none
#include "../enzo/fortran_types.def"

!     Arguments

      R_PREC :: k
      INTG_PREC :: ispectrum
      R_PREC :: omega0, hub, omega_nu, omega_lam, psindex,
     &        omegab0, z, gamma
      R_PREC :: psnorm, growth_factor, kcutoff

!     1) CDM declarations

      R_PREC, parameter :: T0 = 2.726_RKIND
      R_PREC :: q, result

!     2) CHDM declarations

      R_PREC, parameter :: d1 = 0.004321_RKIND,
     &                   d2 = 2.217e-6_RKIND,
     &                   d3 = 11.63_RKIND,
     &                   d4 = 3.317_RKIND
      R_PREC :: pscold, ptotal, gamma_nu, x, x0

!     4) Read a spectrum declarations

      INTG_PREC, parameter :: max_points=400
      INTG_PREC :: i, npoints
      LOGIC_PREC :: found, first_call
      R_PREC :: kvalue(max_points), power(max_points), dummy
      character (len = 80) :: filename, line

      data first_call /.true._LKIND/
      save npoints, kvalue, power, first_call


!     Below the spectral cutoff (if present), set power to zero

      if (k .lt. kcutoff) then
         psfunc = 0.0_RKIND
         return
      endif

      if (ispectrum .eq. 1) then

!  1) CDM Power-Spectrum Bardeen et al 1986 augmented by:
!       Peacock and Dodds 1994 (MNRAS 267, 1020), 
!       Sugiyama 1995 (ApJS 100, 281)

!     Compute shape parameter

         gamma = omega0*hub*(T0/2.7_RKIND)**(-2)*
     &     exp(-omegab0*(1._RKIND + sqrt(2._RKIND*hub)/omega0))

!    &     exp(-2._RKIND*omegab0)  ! old

         q = k/(gamma*hub)

!     Compute z=0 spectrum

         result = psnorm * k**psindex * 
     &         (log(1._RKIND+2.34_RKIND*q)/(2.34_RKIND*q))**2 *
     &         (1._RKIND + 3.89_RKIND*q + (16.1_RKIND*q)**2 + 
     &         (5.46_RKIND*q)**3 + (6.71_RKIND*q)**4)**(-0.5_RKIND)

!     Use growth factor to linearly scale to requested z

         psfunc = result * growth_factor**2
         return

      elseif (ispectrum .eq. 2) then

!  2) CHDM Power-Spectrum From Ma (1996) -- 1 massive neutrino

!     First, compute cold spectrum as above
     
         gamma = omega0*hub*(T0/2.7_RKIND)**(-2)*
     &     exp(-omegab0*(1._RKIND + sqrt(2._RKIND*hub)/omega0))

         q = k/(gamma*hub)
         result = psnorm * k**psindex * 
     &         (log(1._RKIND+2.34_RKIND*q)/(2.34_RKIND*q))**2 *
     &         (1._RKIND + 3.89_RKIND*q + (16.1_RKIND*q)**2 +
     &         (5.46_RKIND*q)**3 + (6.71_RKIND*q)**4)**(-0.5_RKIND)
         pscold = result * growth_factor**2

!     Now modify with Ma eq. 11 (does this assume EdS?)

         gamma_nu = sqrt(1._RKIND/(1._RKIND+z)) * omega_nu * hub**2
         x  = k/gamma_nu
         x0 = k/ (omega_nu * hub**2)
         ptotal = pscold*((1._RKIND + d1*x**(d4/2._RKIND) + d2*x**d4) /
     &        (1._RKIND + d3*x0**d4)       )**(omega_nu**1.05_RKIND)
         psfunc = ptotal
         return

      elseif (ispectrum .eq. 3) then

!  3) Power-law spectrum (scale-free)

         result = psnorm * k**psindex
         psfunc = result * growth_factor**2
         return

      elseif (ispectrum .eq. 4) then

!  4) Read a power spectrum from a file
!     Spectrum must power ordered by increasing k values
!     (format: dummy, k, P(k) on each line, comment symbol is #)

!     This version ignores the normalization and z, etc.

!        If this is the first call call, prompt for name and read spectrum

         if (first_call) then

            write(*,*) 'Enter ps name:'
            read(*,'(a80)') filename

            open(1, file=filename, status='old')

!           Loop until done (ignore lines starting with #)

            npoints = 0

            do while (.true.)
               read(1, '(a80)', end=100) line
               if (line(1:1) .ne. '#') then
                  npoints = npoints + 1
                  read(line, *) dummy, kvalue(npoints), power(npoints)
!                 write(0,*) npoints, kvalue(npoints), power(npoints)
               endif
            enddo

 100        continue

!           Close file and set first_call to false

            close(1)

            first_call = .false._LKIND

         endif

!        Search for given k value in spectrum

         found = .false._LKIND
         i = 0

         do while (.not. found .and. i .lt. npoints)
            i = i + 1
            if (k .lt. kvalue(i)) found = .true._LKIND
         enddo

!        Linearly interpolate or set power to zero if outside spectrum

         if (i .eq. 1 .or. .not. found) then
            result = 0._RKIND
         else
            result = power(i-1) + (power(i)    - power(i-1)) *
     &                           (k           - kvalue(i-1))/
     &                           (kvalue(i)   - kvalue(i-1))
         endif

         psfunc = result
         return

      elseif (ispectrum .eq. 5) then

!  5) CHDM Power-Spectrum From Ma (1996) -- 2 neutrinos, equal masses
!     (I just doubled the free-streaming length)

!     First, compute cold spectrum as above
     
         gamma = omega0*hub*(T0/2.7_RKIND)**(-2)*
     &     exp(-omegab0*(1._RKIND + sqrt(2._RKIND*hub)/omega0))

         q = k/(gamma*hub)
         result = psnorm * k**psindex * 
     &         (log(1._RKIND+2.34_RKIND*q)/(2.34_RKIND*q))**2 *
     &         (1._RKIND + 3.89_RKIND*q + (16.1_RKIND*q)**2 + 
     &         (5.46_RKIND*q)**3 + (6.71_RKIND*q)**4)**(-0.5_RKIND)
         pscold = result * growth_factor**2

!     Now modify with Ma eq. 11 (does this assume EdS?)
!       (just half gamma_nu to reflect the doubled velocity)

         gamma_nu = sqrt(1._RKIND/(1._RKIND+z)) * omega_nu 
     &        * hub**2 / 2._RKIND
         x  = k/gamma_nu
         x0 = k/(omega_nu * hub**2 / 2._RKIND)
         ptotal = pscold*((1._RKIND + d1*x**(d4/2._RKIND) + d2*x**d4) /
     &        (1._RKIND + d3*x0**d4)       )**(omega_nu**1.05_RKIND)

         psfunc = ptotal
         return

      elseif (ispectrum .eq. 6) then

!  6) CDM-like power spectrum with a fixed shape parameter.
!
!  BWO: This is included for historical interest.  The spectrum is
!  the same as the bardeen power spectrum (type #1) but with a free
!  parameter, gamma, that allows the spectral shape to be modified 
!  independently of the cosmological parameters.  This is now
!  HIGHLY DEPRECATED (i.e., wrong)

!     Compute shape parameter

!        gamma = 0.25   (now set in common block)

         q = k/(gamma*hub)

!     Compute z=0 spectrum

         result = psnorm * k**psindex * 
     &         (log(1._RKIND+2.34_RKIND*q)/(2.34_RKIND*q))**2 *
     &         (1._RKIND + 3.89_RKIND*q + (16.1_RKIND*q)**2 + 
     &         (5.46_RKIND*q)**3 +  (6.71_RKIND*q)**4)**(-0.5_RKIND)

!     use growth factor to linearly scale to requested z

         psfunc = result * growth_factor**2
         return

      else

         write(0,*) 'PSFUNC: unrecognized ispectrum:', ispectrum
         stop

      endif

      end

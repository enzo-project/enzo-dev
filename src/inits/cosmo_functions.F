#include "../enzo/fortran.def"
!=======================================================================
! Set the current radiation density (incl. massless neutrinos) for h=1
!define OMEGA0_RAD 4.22e-5
#define OMEGA0_RAD 0.0
!=======================================================================

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
!//////////////////////  FUNCTION COMPUTE_TIME  \\\\\\\\\\\\\\\\\\\\\\\\

      R_PREC function calc_time(aye_temp)

!  COMPUTES THE TIME GIVEN THE EXPANSION FACTOR
!
!  written by: Greg Bryan
!  date:       February, 1992
!  modified:   Robert Harkness
!  date:       November, 2003
!
!  PURPOSE:  
!
!  INPUTS:
!
!  OUTPUTS:

      implicit none
#include "../enzo/fortran_types.def"

#include "cosmo.h"

!     Argument

      R_PREC :: aye_temp

!     Locals

      R_PREC :: time_temp

!     Externals

      R_PREC :: dtda, midpnt
      external :: dtda, midpnt


      call qromo(dtda, 0.0_RKIND, aye_temp, time_temp, midpnt)

      calc_time = time_temp

      return
      end


!=======================================================================
!/////////////////////////  FUNCTION DADT  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      R_PREC function calc_ayed(aye_temp)

!  modified:   Robert Harkness
!  date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

!     Argument

      R_PREC :: aye_temp

!     External function

      R_PREC :: dtda


      calc_ayed = 1.0_RKIND/dtda(aye_temp)

      return
      end


!=======================================================================
!////////////////////////  FUNCTION D2A/DT2  \\\\\\\\\\\\\\\\\\\\\\\\\\\

      R_PREC function calc_ayedd(aye_temp)

!  modified:   Robert Harkness
!  date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

#include "cosmo.h"

!     Argument

      R_PREC :: aye_temp

!     Locals

      R_PREC :: omega0_rad, omega0_mrad

      omega0_rad = OMEGA0_RAD*hub**(-2)
      omega0_mrad = omega0 - omega0_rad

      calc_ayedd = aye_temp*(lam0 
     &     - 0.5_RKIND*omega0_mrad/(aye_temp*uaye)**3
     &     - omega0_rad /(aye_temp*uaye)**4)

!     Convert to code units (the factors of H_0 have already been cancelled)

     &     /(sqrt(1.5_RKIND*omega0)*(1.0_RKIND+zri)**1.5_RKIND)**2

      return
      end


!=======================================================================
!/////////////////////////  FUNCTION DTDA  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      R_PREC function dtda(aye_temp)

!  COMPUTES THE VALUE OF DT/DA, GIVEN A

!  written by: Greg Bryan
!  date:       February, 1992
!  modified:   Robert Harkness
!  date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

#include "cosmo.h"

!     Argument

      R_PREC :: aye_temp

!     Locals

      R_PREC :: at2, omega0_rad, omega0_mrad


!     We include the (small) effect of radiation (Peebles 6.80)

      omega0_rad = OMEGA0_RAD*hub**(-2)
      omega0_mrad = omega0 - omega0_rad

!     Convert aye from code units

      at2 = aye_temp*uaye

!     Compute dt/da (Peebles1993, p. 312)

      dtda = 1._RKIND/sqrt(omega0_mrad/at2 + omega0_rad/at2**2 +
     &                lam0*at2**2 + 1._RKIND - lam0 - omega0)

!     Convert to code units (the factors of H_0 have already been cancelled)

     &       *sqrt(1.5_RKIND*omega0)*(1._RKIND+zri)**1.5_RKIND * uaye

      return
      end


!=======================================================================
!////////////////////////  FUNCTION CALC_AYE  \\\\\\\\\\\\\\\\\\\\\\\\\\

      R_PREC function calc_aye(time_temp)

!  COMPUTES THE EXPANSION FACTOR GIVEN THE TIME
!
!  written by: Greg Bryan
!  date:       February, 1992
!  modified:   Robert Harkness
!  date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

#include "cosmo.h"

!     Argument

      R_PREC :: time_temp

!     Locals

      R_PREC :: aye_temp, aye_old, calc_time, dtda, tfinal, tfromfinal
      INTG_PREC :: i

!     Parameters

      INTG_PREC, parameter :: niter = 10
      R_PREC, parameter :: tolerance = 1.e-5_RKIND


!     Make an initial guess based on Taylor expansion (i.e. use q0)

      tfinal = calc_time(1.0_RKIND+zri)
      tfromfinal = sqrt(2._RKIND/3._RKIND/omega0) * 
     &     (1._RKIND+zri)**(-1.5_RKIND) * (tfinal - time_temp)
      aye_temp = (1._RKIND+zri)*(1._RKIND - tfromfinal - 
     &     0.5_RKIND*(0.5_RKIND*omega0 - lam0)*tfromfinal**2)

!     Do a Newton-Raphson iteration

      do i = 1, niter
         aye_old = aye_temp
         aye_temp = aye_old + 1._RKIND/dtda(aye_old) *
     &                        (time_temp - calc_time(aye_old))
         if (abs(aye_old-aye_temp)/aye_temp .lt. tolerance) goto 100
      enddo

      write(0,*) 'NR in calc_aye failed.'
      stop

!     Done

 100  continue

      calc_aye = aye_temp

      return
      end


!=======================================================================
!/////////////////////////  FUNCTION CALC_F  \\\\\\\\\\\\\\\\\\\\\\\\\\\

      R_PREC function calc_f(aye_temp)

!  COMPUTES THE FUNCTION D LOG (DELTA_PLUS) / D LOG (AYE)
!
!  written by: Greg Bryan
!  date:       February, 1995
!  modified:   Robert Harkness
!  date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

#include "cosmo.h"

!     Argument

      R_PREC :: aye_temp

!     Locals

      R_PREC :: ayed, ayedd, sum, at2

      R_PREC :: fhelper, calc_ayed, calc_ayedd, midpnt

      external :: fhelper, midpnt


!     We calculate f(z) through PEEBLES93, eq 13.81:
!
!     f = ayedd*aye/ayed^2 - 1 + 1 / (a^2 E^3 G)
!
!        where G = int(z to infifinty) (1+z)/E(z)^3 dz
!                = int(0 to a) da/(E(a) * a)^3
!              E = da/dt / H_0 

!     Compute G (using the usual a convention)

      at2 = aye_temp*uaye

      call qromo(fhelper, 0.0_RKIND, at2, sum, midpnt)

      ayed = calc_ayed(aye_temp)
      ayedd = calc_ayedd(aye_temp)

!     Note that f is dimensionless, specifically all the unusual aye
!     convention (a=1 at z=zri) cancels for the first term.  The third
!     term is computed entirely in the usual convention.

      calc_f = ayedd*aye_temp/ayed**2 - 1.0_RKIND + at2*fhelper(at2)/sum

      return
      end


!=======================================================================
!/////////////////////////  FUNCTION FFUNC  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      R_PREC function fhelper(at2)

!  modified:   Robert Harkness
!  date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

#include "cosmo.h"

!     returns the function 1/(a * E)^3  
!       (where, here only, a is not aye, i.e. it obeys the usual convention
!        of a = 1 at z = 0)

!     Argument

      R_PREC :: at2

!     Locals

      R_PREC :: E, omega0_rad, omega0_mrad


      omega0_rad = OMEGA0_RAD*hub**(-2)
      omega0_mrad = omega0 - omega0_rad

      E = sqrt(omega0_mrad/at2**3 + (1.0_RKIND - omega0 - lam0)/at2**2 + 
     &         lam0 + omega0_rad/at2**4)

      fhelper = 1.0_RKIND/(at2*E)**3

      return
      end


!=======================================================================
!///////////////////////  FUNCTION CALC_GROWTH  \\\\\\\\\\\\\\\\\\\\\\\\

      R_PREC function calc_growth(z1)

!  COMPUTES THE GROWTH FUNCTION D(z), NORMALIZED AT HIGH REDSHIFT
!
!  written by: Greg Bryan
!  date:       February, 1995
!  modified:   Robert Harkness
!  date:       November, 2003

!  PURPOSE:
!    Note: this function is _not_ normalized to D(z=0) = 1
!
!  INPUTS:
!    z      - redshift function is to be evaluated at
!    omega0 - matter density ratio at z=0
!    omega_lam - \Lambda/(3*H_0^2) at z=0
!
!  OUTPUTS:
!    calc_growth - linear growth, normalized to 1/(1+z) at z=infinity


      implicit none
#include "../enzo/fortran_types.def"

#include "cosmo.h"

!     Argument

      R_PREC :: z1

!     Locals

      R_PREC :: a, sum, E

!     Externals

      R_PREC :: fhelper, midpnt
      external :: fhelper, midpnt


!     We calculate D(z) through PEEBLES93, eq 13.78:
!
!     D(z) = E(z) * G(z)
!
!        where G = 5*omega0/2 * int(z to infinity) (1+z)/E(z)^3 dz
!                = 5*omega0/2 * int(0 to a) da/(E(a) * a)^3
!              E = da/dt / H_0

!     Compute G, missing the prefactor

      a = 1.0_RKIND/(1.0_RKIND+z1)

      call qromo(fhelper, 1.e-10_RKIND, a, sum, midpnt)

      E = sqrt(omega0/a**3 + (1._RKIND - omega0 - lam0)/a**2 + lam0)

      calc_growth = 5.0_RKIND*omega0/2.0_RKIND * sum * E

      return
      end


!=======================================================================
!/////////////////////  SUBROUTINE SET_COMMON  \\\\\\\\\\\\\\\\\\\\\\\\\

      subroutine set_common(lam0_in, omega0_in, zri_in, hub_in)

!  modified:   Robert Harkness
!  date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

#include "cosmo.h"

!     Arguments

      R_PREC :: lam0_in, omega0_in, zri_in, hub_in


      lam0   = lam0_in
      omega0 = omega0_in
      zri    = zri_in
      hub    = hub_in
      uaye   = 1.0_RKIND/(1.0_RKIND + zri)

      return
      end

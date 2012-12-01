#include "fortran.def"
!=======================================================================
!
! Copyright 2006 Daniel R. Reynolds
! Copyright 2006 Laboratory for Computational Astrophysics
! Copyright 2006 Regents of the University of California
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine gFLDProblem_GasEnergySource(ecsrc, time, Era, eca, nHIa,     &
     nHeIa, nHeIIa, Temp, rhoa, eha, vxa, vya, vza, a, Model, ProbType, &
     Nchem, HFrac, aUnits, DenUnits, VelUnits, LenUnits, TimeUnits,     &
     ErUnits, ecUnits, NiUnits, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,     &
     NGzl, NGzr, x0L, x0R, x1L, x1R, x2L, x2R, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2007
!
!  PURPOSE: Computes any external sources to the gas energy correction 
!           equation.
!
!  INPUTS:
!     time       - simulation time for evaluation
!     Era        - radiation energy density array
!     eca        - specific gas energy correction array
!     nHIa       - Hydrogen I density array
!     nHeIa      - Helium I density array
!     nHeIIa     - Helium II density array
!     Temp       - fluid temperature array
!     rhoa       - fluid density array
!     eha        - specific gas energy array
!     vxa        - x-directional gas velocity array
!     vya        - y-directional gas velocity array
!     vza        - z-directional gas velocity array
!     a          - cosmological expansion parameter
!     Model      - flag denoting physical model to use
!     ProbType   - flag denoting problem to run
!     Nchem      - number of chemical species
!     HFrac      - percentage of mass composed of Hydrogen
!     *Units     - variable scaling constants
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     x*L/x*R    - left/right subdomain boundaries (comoving, no ghosts)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     ecsrc      - array of gas energy sources
!     ier        - success/failure flag (1->success, 0->failure)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
  implicit none
#include "fortran_types.def"
!#ifdef USE_MPI
!#include "mpif.h"
!#endif

!--------------
! argument declarations
  INTG_PREC, intent(in) :: Nchem, Model, ProbType
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC,  intent(in) :: a
  R_PREC, intent(in) :: time, HFrac
  R_PREC, intent(in) :: aUnits, DenUnits, VelUnits, LenUnits, TimeUnits, &
       ErUnits, ecUnits, NiUnits
  R_PREC, intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  R_PREC, intent(in),                                        &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) &
       :: Era, eca, nHIa, nHeIa, nHeIIa, Temp, rhoa, eha, vxa, vya, vza
  R_PREC, intent(out) ::                                     &
       ecsrc(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  
!--------------
! locals
  INTG_PREC :: i, j, k

!=======================================================================

!!$  write(*,*) 'Entering gFLDProblem::GasEnergySource routine'

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  ecsrc = 0._RKIND


  if (ProbType == 416) then
!!$     if (time/TimeUnits > 0.005) then
        ecsrc = 1.e4_RKIND
!!$     end if
  end if


  ! compute gas energy sources based on physical model
  !   [no models use this feature at the current time]

  return
end subroutine gFLDProblem_GasEnergySource
!=======================================================================

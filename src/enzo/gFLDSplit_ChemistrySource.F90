#include "fortran.def"
!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine gFLDSplit_ChemistrySource(HIsrc, HeIsrc, HeIIsrc, time, a,   &
     ProbType, Nchem, HFrac, aUnits, LenUnits, TimeUnits, NiUnits, Nx,  &
     Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, x0L, x0R, x1L, x1R,    &
     x2L, x2R, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!
!  PURPOSE: Computes any external sources to the chemical ionization 
!           equations
!
!  INPUTS:
!     time       - simulation time for evaluation
!     a          - cosmological expansion parameter
!     ProbType   - flag denoting problem to use
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
!     HIsrc      - array of HI sources
!     HeIsrc     - array of HeI sources
!     HeIIsrc    - array of HeII sources
!     ier        - success/failure flag (1->success, 0->failure)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
  implicit none
#include "fortran_types.def"

!--------------
! argument declarations
  INTG_PREC, intent(in) :: Nchem, ProbType
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in) :: a
  R_PREC,    intent(in) :: time, HFrac
  R_PREC,    intent(in) :: aUnits, LenUnits, TimeUnits, NiUnits
  R_PREC,    intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(out) &
       :: HIsrc, HeIsrc, HeIIsrc
  
!--------------
! locals
  INTG_PREC :: i, j, k

!=======================================================================

!!$  write(*,*) 'Entering gFLDSplit::ChemistrySource routine'

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  if (Nchem == 1) then
     HIsrc = 0._RKIND
  endif
  if (Nchem == 3) then
     HIsrc   = 0._RKIND
     HeIsrc  = 0._RKIND
     HeIIsrc = 0._RKIND
  endif


  ! compute chemistry sources based on ProbType
  !   [no problems use this feature at the current time]

  return
end subroutine gFLDSplit_ChemistrySource
!=======================================================================

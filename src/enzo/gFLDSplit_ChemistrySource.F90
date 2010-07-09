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
#include "fortran.def"

!--------------
! argument declarations
  integer, intent(in) :: Nchem, ProbType
  integer, intent(in) :: Nx, NGxl, NGxr
  integer, intent(in) :: Ny, NGyl, NGyr
  integer, intent(in) :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in) :: a
  real,    intent(in) :: time, HFrac
  real,    intent(in) :: aUnits, LenUnits, TimeUnits, NiUnits
  real,    intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(out) &
       :: HIsrc, HeIsrc, HeIIsrc
  
!--------------
! locals
  integer :: i, j, k

!=======================================================================

!!$  write(*,*) 'Entering gFLDSplit::ChemistrySource routine'

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  if (Nchem == 1) then
     HIsrc = 0.d0
  endif
  if (Nchem == 3) then
     HIsrc   = 0.d0
     HeIsrc  = 0.d0
     HeIIsrc = 0.d0
  endif


  ! compute chemistry sources based on ProbType
  !   [no problems use this feature at the current time]

  return
end subroutine gFLDSplit_ChemistrySource
!=======================================================================

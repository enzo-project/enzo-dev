!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine gFLDSplit_GasEnergySource(ecsrc, time, a, ProbType, aUnits, &
     VelUnits, LenUnits, TimeUnits, ecUnits, Nx, Ny, Nz, NGxl, NGxr,   &
     NGyl, NGyr, NGzl, NGzr, x0L, x0R, x1L, x1R, x2L, x2R, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!
!  PURPOSE: Computes any external sources to the gas energy correction 
!           equation.
!
!  INPUTS:
!     time       - simulation time for evaluation
!     a          - cosmological expansion parameter
!     ProbType   - flag denoting problem to run
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
#include "fortran.def"

!--------------
! argument declarations
  integer, intent(in) :: ProbType
  integer, intent(in) :: Nx, NGxl, NGxr
  integer, intent(in) :: Ny, NGyl, NGyr
  integer, intent(in) :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in) :: a
  real,    intent(in) :: time
  real,    intent(in) :: aUnits, VelUnits, LenUnits, TimeUnits, ecUnits
  real,    intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  real,    intent(out) :: ecsrc(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  
!--------------
! locals
  integer :: i, j, k

!=======================================================================

!!$  write(*,*) 'Entering gFLDSplit::GasEnergySource routine'

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  ecsrc = 0.d0

  return
end subroutine gFLDSplit_GasEnergySource
!=======================================================================

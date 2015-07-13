#include "fortran.def"
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
#include "fortran_types.def"

!--------------
! argument declarations
  INTG_PREC, intent(in) :: ProbType
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in) :: a
  R_PREC,    intent(in) :: time
  R_PREC,    intent(in) :: aUnits, VelUnits, LenUnits, TimeUnits, ecUnits
  R_PREC,    intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  R_PREC,    intent(out) :: ecsrc(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  
!--------------
! locals
  INTG_PREC :: i, j, k

!=======================================================================

!!$  write(*,*) 'Entering gFLDSplit::GasEnergySource routine'

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  ecsrc = 0._RKIND

  return
end subroutine gFLDSplit_GasEnergySource
!=======================================================================

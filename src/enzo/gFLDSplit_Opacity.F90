#include "fortran.def"
#include "phys_const.def"
!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine gFLDSplit_Opacity(kappaE, time, rho, n_HI, n_HeI, n_HeII, &
     a, Model, IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII,  &
     IsEsHeIInu, EmC0, EmC1, EmC2, aUnits, DenUnits, LenUnits,       &
     TimeUnits, NiUnits, Nchem, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,  &
     NGzl, NGzr, x0L, x0R, x1L, x1R, x2L, x2R, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!
!  PURPOSE: Computes the local opacities over the domain 
!           (NOTE: these are normalized by NiUnits, and must be 
!            converted to the appropriate redshift-dependent values 
!            in physics routines via kappaE_true = kappaE*NiUnits).
!
!  INPUTS:
!     time       - simulation time for evaluation
!     rho        - gas density
!     n_HI       - proportional density of HI species
!     n_HeI      - proportional density of HeI species
!     n_HeII     - proportional density of HeII species
!     a          - cosmological expansion parameter
!     Model      - flag denoting physical model to use
!     IsE        - int_{nu0_HI}^{inf} sigE dnu
!     IsEsHI     - int_{nu0_HI}^{inf} sigE*sigHI dnu
!     IsEsHInu   - int_{nu0_HI}^{inf} sigE*sigHI/nu dnu
!     IsEsHeI    - int_{nu0_HeI}^{inf} sigE*sigHeI dnu
!     IsEsHeInu  - int_{nu0_HeI}^{inf} sigE*sigHeI/nu dnu
!     IsEsHeII   - int_{nu0_HeII}^{inf} sigE*sigHeII dnu
!     IsEsHeIInu - int_{nu0_HeII}^{inf} sigE*sigHeII/nu dnu
!     EmC0-EmC2  - input parameters for Energy mean opacity
!     *Units     - variable scaling constants
!     Nchem      - number of chemistry species in problem {0, 1, 3}
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     x*L/x*R    - left/right subdomain boundaries (comoving, no ghosts)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     kappaE     - Opacity coefficient (Energy mean)
!     ier        - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in) :: Model, Nchem
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in) :: a
  R_PREC,    intent(in) :: time, IsE, IsEsHI, IsEsHInu, IsEsHeI
  R_PREC,    intent(in) :: IsEsHeInu, IsEsHeII, IsEsHeIInu
  R_PREC,    intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  R_PREC,    intent(in) :: EmC0, EmC1, EmC2
  R_PREC,    intent(in) :: aUnits, DenUnits, LenUnits, TimeUnits, NiUnits
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) &
       :: rho, n_HI, n_HeI, n_HeII
  R_PREC, intent(out) :: kappaE(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)

  !--------------
  ! locals
  INTG_PREC :: i, j, k
  R_PREC    :: mp, HIconst, HeIconst, HeIIconst
  R_PREC    :: rhoval, Tval

  !=======================================================================

!!$  write(*,*) 'Entering gFLDSplit::Opacity routine'

  ! initialize outputs to zero, flag to success
  kappaE = 0._RKIND
  ier = 1

  ! set shortcut values
  mp = mass_h                    ! mass of a proton [g]

  ! compute opacity shortcuts
  HIconst   = IsEsHI/IsE
  HeIconst  = IsEsHeI/IsE/4._RKIND
  HeIIconst = IsEsHeII/IsE/4._RKIND

  ! compute opacity over domain depending on number of chemical species 

  !   power-law test case (ZEUS-like opacities)
  if (Model == 10) then
     do k=1-NGzl,Nz+NGzr,1
        do j=1-NGyl,Ny+NGyr,1
           do i=1-NGxl,Nx+NGxr,1
              rhoval = rho(i,j,k)*DenUnits
              kappaE(i,j,k) = EmC0 * (rhoval/EmC1)**EmC2 / NiUnits
           enddo
        enddo
     enddo

  !   chemistry-based opacity (default)
  else

     ! Hydrogen only
     if (Nchem == 1) then
        do k=1-NGzl,Nz+NGzr,1
           do j=1-NGyl,Ny+NGyr,1
              do i=1-NGxl,Nx+NGxr,1
                 kappaE(i,j,k) = n_HI(i,j,k)*HIconst
              enddo
           enddo
        enddo

     ! Hydrogen plus Helium
     else if (Nchem == 3) then
        do k=1-NGzl,Nz+NGzr,1
           do j=1-NGyl,Ny+NGyr,1
              do i=1-NGxl,Nx+NGxr,1
                 kappaE(i,j,k) = n_HI(i,j,k)  *HIconst   &
                               + n_HeI(i,j,k) *HeIconst  &
                               + n_HeII(i,j,k)*HeIIconst
              enddo
           enddo
        enddo

     else
        write(0,*) 'gFLDSplit_Opacity ERROR: illegal Nchem =',Nchem, &
             ', Model =',Model,' requires Nchem = {1, 3}'
        ier = 0
     endif

  end if  ! Model

  return

end subroutine gFLDSplit_Opacity
!=======================================================================

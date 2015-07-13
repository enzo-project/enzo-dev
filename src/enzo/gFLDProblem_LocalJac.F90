#include "fortran.def"
#include "phys_const.def"
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
subroutine gFLDProblem_LocalJac(Erjac_Er, Erjac_ec, Erjac_HI, Erjac_HeI,  &
     Erjac_HeII, ecjac_Er, ecjac_ec, ecjac_HI, ecjac_HeI, ecjac_HeII,     &
     HIjac_Er, HIjac_ec, HIjac_HI, HIjac_HeI, HIjac_HeII, HeIjac_Er,      &
     HeIjac_ec, HeIjac_HI, HeIjac_HeI, HeIjac_HeII, HeIIjac_Er,           &
     HeIIjac_ec, HeIIjac_HI, HeIIjac_HeI, HeIIjac_HeII, time, Era, eca,   &
     n_HIa, n_HeIa, n_HeIIa, eha, rhoa, vx, vy, vz, Nchem, HFrac,         &
     Model, ESpectrum, ProbType, DualEnergy, a, adot, CompA, Comp_xray,   &
     Comp_temp, IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII,      &
     IsEsHeIInu, PmC0, PmC1, PmC2, PmC3, PmC4, EmC0, EmC1, EmC2, EmC3,    &
     EmC4, gamma, NTempBins, TempStart, TempEnd, k1Tb, k2Tb, k3Tb, k4Tb,  &
     k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb,    &
     ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, piHI,    &
     piHeI, piHeII, aUnits, DenUnits, VelUnits, LenUnits, ErUnits,        &
     ecUnits, NiUnits, ecScale, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, &
     NGyr, NGzl, NGzr, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August, 2006
  !  modified1:  September 19, 2007, by John Hayes; minor mods on behalf of
  !              ProbTypes 204 & 205.
  !
  !  PURPOSE: Computes the spatially-local components of the chemical 
  !           Jacobian for the Gray FLD problem.
  !
  !  INPUTS:
  !     time       - simulation time for evaluation
  !     Era        - Gray radiation energy density
  !     eca        - fluid energy correction array
  !     n_HIa      - number density of Hydrogen I species
  !     n_HeIa     - number density of Helium I species
  !     n_HeIIa    - number density of Helium II species
  !     eha        - total fluid energy array
  !     rhoa       - fluid density, assumed partitioned into either 
  !     vx         - fluid velocity (x-direction)
  !     vy         - fluid velocity (y-direction)
  !     vz         - fluid velocity (z-direction)
  !     Nchem      - number of chemical species (allowed: 0, 1, 3)
  !     HFrac      - percentage of mass composed of Hydrogen
  !     Model      - flag denoting physical model to use
  !     ESpectrum  - radiation spectrum choice
  !                       1 -> 1e5 black body spectrum
  !                       0 -> power law spectrum
  !                      -1 -> monochromatic 
  !     ProbType   - flag denoting problem type (kluge)
  !     DualEnergy - flag denoting dual energy formalism
  !     a          - cosmological expansion parameter
  !     adot       - da/dt
  !     CompA      - Compton cooling coefficient 1 (multiplier)
  !     Comp_xray  - X-ray Compton heating coefficient
  !     Comp_temp  - X-ray Compton heating temperature 
  !     IsE        - int_{nu0_HI}^{inf} sigE dnu
  !     IsEsHI     - int_{nu0_HI}^{inf} sigE*sigHI dnu
  !     IsEsHInu   - int_{nu0_HI}^{inf} sigE*sigHI/nu dnu
  !     IsEsHeI    - int_{nu0_HeI}^{inf} sigE*sigHeI dnu
  !     IsEsHeInu  - int_{nu0_HeI}^{inf} sigE*sigHeI/nu dnu
  !     IsEsHeII   - int_{nu0_HeII}^{inf} sigE*sigHeII dnu
  !     IsEsHeIInu - int_{nu0_HeII}^{inf} sigE*sigHeII/nu dnu
  !     PmC0-PmC4  - input parameters for Planck mean opacity
  !     EmC0-EmC4  - input parameters for Energy mean opacity
  !     gamma      - ideal gas parameter (ratio of specific heats)
  !     *Units     - variable scaling constants
  !     dx,dy,dz   - mesh spacings in each direction
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: all jacobian components are scaled to correspond 
  !                    with the internally-stored Enzo units for the 
  !                    relevant time-dependent equation
  !     Erjac_Er     - local radiation energy Jacobian w.r.t. Er
  !     Erjac_ec     - local radiation energy Jacobian w.r.t. ec
  !     Erjac_HI     - local radiation energy Jacobian w.r.t. nHI
  !     Erjac_HeI    - local radiation energy Jacobian w.r.t. nHeI
  !     Erjac_HeII   - local radiation energy Jacobian w.r.t. nHeII
  !     ecjac_Er     - local gas energy Jacobian w.r.t. Er
  !     ecjac_ec     - local gas energy Jacobian w.r.t. ec
  !     ecjac_HI     - local gas energy Jacobian w.r.t. nHI
  !     ecjac_HeI    - local gas energy Jacobian w.r.t. nHeI
  !     ecjac_HeII   - local gas energy Jacobian w.r.t. nHeII
  !     HIjac_Er     - local HI chemistry Jacobian w.r.t. Er
  !     HIjac_ec     - local HI chemistry Jacobian w.r.t. ec
  !     HIjac_HI     - local HI chemistry Jacobian w.r.t. nHI
  !     HIjac_HeI    - local HI chemistry Jacobian w.r.t. nHeI
  !     HIjac_HeII   - local HI chemistry Jacobian w.r.t. nHeII
  !     HeIjac_Er    - local HeI chemistry Jacobian w.r.t. Er
  !     HeIjac_ec    - local HeI chemistry Jacobian w.r.t. ec
  !     HeIjac_HI    - local HeI chemistry Jacobian w.r.t. nHI
  !     HeIjac_HeI   - local HeI chemistry Jacobian w.r.t. nHeI
  !     HeIjac_HeII  - local HeI chemistry Jacobian w.r.t. nHeII
  !     HeIIjac_Er   - local HeII chemistry Jacobian w.r.t. Er
  !     HeIIjac_ec   - local HeII chemistry Jacobian w.r.t. ec
  !     HeIIjac_HI   - local HeII chemistry Jacobian w.r.t. nHI
  !     HeIIjac_HeI  - local HeII chemistry Jacobian w.r.t. nHeI
  !     HeIIjac_HeII - local HeII chemistry Jacobian w.r.t. nHeII
  !     ier          - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in) :: Nchem, Model, ESpectrum, ProbType, DualEnergy
  INTG_PREC, intent(in) :: NTempBins
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in) :: a, adot
  R_PREC,    intent(in) :: dx, dy, dz
  R_PREC,    intent(in) :: time, HFrac, gamma, IsE, IsEsHI, IsEsHInu
  R_PREC,    intent(in) :: IsEsHeI, IsEsHeInu, IsEsHeII, IsEsHeIInu
  R_PREC,    intent(in) :: CompA, Comp_xray, Comp_temp
  R_PREC,    intent(in) :: TempStart, TempEnd, piHI, piHeI, piHeII
  R_PREC,    intent(in) :: PmC0, PmC1, PmC2, PmC3, PmC4
  R_PREC,    intent(in) :: EmC0, EmC1, EmC2, EmC3, EmC4
  R_PREC,    intent(in) :: aUnits, DenUnits, VelUnits, LenUnits,  &
       ErUnits, ecUnits, NiUnits, ecScale
  R_PREC, intent(in),                                             &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)  &
       :: Era, eca, n_HIa, n_HeIa, n_HeIIa, eha, rhoa, vx, vy, vz
  R_PREC, intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb, k4Tb,      &
       k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, &
       ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb
  R_PREC, intent(out),                                           &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       Erjac_Er, Erjac_ec, Erjac_HI, Erjac_HeI, Erjac_HeII,       &
       ecjac_Er, ecjac_ec, ecjac_HI, ecjac_HeI, ecjac_HeII,       &
       HIjac_Er, HIjac_ec, HIjac_HI, HIjac_HeI, HIjac_HeII,       &
       HeIjac_Er, HeIjac_ec, HeIjac_HI, HeIjac_HeI, HeIjac_HeII,  &
       HeIIjac_Er, HeIIjac_ec, HeIIjac_HI, HeIIjac_HeI, HeIIjac_HeII

  !--------------
  ! locals
  INTG_PREC :: i, j, k
  R_PREC :: dxi2, dyi2, dzi2, DivV, GradRhoDotV, aval

  !=======================================================================

  ! set flag to success
  ier = 1

  ! check that chemistry constants make sense
  if ((Nchem /= 0) .and. (Nchem /= 1) .and. (Nchem /= 3)) then
     write(*,*) 'Chemistry ERROR: illegal value, Nchem = ',Nchem, &
          ',  Nchem must be one of {0, 1, 3}.  Returning!'
     ier = 0
     return
  endif
  if ((HFrac < 0._RKIND) .or. (HFrac > 1._RKIND)) then
     write(*,*) 'Chemistry ERROR: illegal value, HFrac = ',HFrac, &
          ',  value must be in the inteval [0,1].  Returning!'
     ier = 0
     return     
  endif

  ! set shortcut values
  aval = a*aunits
  dxi2 = 0.5_RKIND*a/(dx*LenUnits)   ! convert to proper units during transition
  dyi2 = 0.5_RKIND*a/(dy*LenUnits)
  dzi2 = 0.5_RKIND*a/(dz*LenUnits)

  !  First compute correction-based adjustments alone
  !  (as these involve derivatives, use loop for appropriate dimension)
  ecjac_ec = 0._RKIND
!!$  if (Ny == 1) then
!!$     do k=1,Nz,1
!!$        do j=1,Ny,1
!!$           do i=1,Nx,1
!!$
!!$              ! velocity divergence
!!$              DivV = dxi2*(vx(i+1,j,k)-vx(i-1,j,k)) * VelUnits
!!$
!!$              ! (grad density).dot.(velocity)/density
!!$              GradRhoDotV = dxi2*(rhoa(i+1,j,k)-rhoa(i-1,j,k)) &
!!$                   * vx(i,j,k) / rhoa(i,j,k) * VelUnits
!!$
!!$              ! put it together
!!$              ecjac_ec(i,j,k) = ecUnits/aval*(DivV - (gamma-1._RKIND)*GradRhoDotV)
!!$
!!$           enddo
!!$        enddo
!!$     enddo
!!$     ! 2D model
!!$  else if (Nz == 1) then
!!$     do k=1,Nz,1
!!$        do j=1,Ny,1
!!$           do i=1,Nx,1
!!$
!!$              ! velocity divergence
!!$              DivV = (dxi2*(vx(i+1,j,k)-vx(i-1,j,k))  &
!!$                    + dyi2*(vy(i,j+1,k)-vy(i,j-1,k))) * VelUnits
!!$
!!$              ! (grad density).dot.(velocity)/density
!!$              GradRhoDotV = (dxi2*(rhoa(i+1,j,k)-rhoa(i-1,j,k))*vx(i,j,k)  &
!!$                   + dyi2*(rhoa(i,j+1,k)-rhoa(i,j-1,k))*vy(i,j,k)) &
!!$                   / rhoa(i,j,k)*VelUnits
!!$
!!$              ! put it together
!!$              ecjac_ec(i,j,k) = ecUnits/aval*(DivV - (gamma-1._RKIND)*GradRhoDotV)
!!$
!!$           enddo
!!$        enddo
!!$     enddo
!!$     ! 3D model
!!$  else
!!$     do k=1,Nz,1
!!$        do j=1,Ny,1
!!$           do i=1,Nx,1
!!$
!!$              ! velocity divergence
!!$              DivV = (dxi2*(vx(i+1,j,k)-vx(i-1,j,k))  &
!!$                    + dyi2*(vy(i,j+1,k)-vy(i,j-1,k))  &
!!$                    + dzi2*(vz(i,j,k+1)-vz(i,j,k-1))) * VelUnits
!!$
!!$              ! (grad density).dot.(velocity)/density
!!$              GradRhoDotV = (dxi2*(rhoa(i+1,j,k)-rhoa(i-1,j,k))*vx(i,j,k)  &
!!$                           + dyi2*(rhoa(i,j+1,k)-rhoa(i,j-1,k))*vy(i,j,k)  &
!!$                           + dzi2*(rhoa(i,j,k+1)-rhoa(i,j,k-1))*vz(i,j,k)) &
!!$                          / rhoa(i,j,k) * VelUnits
!!$
!!$              ! put it together
!!$              ecjac_ec(i,j,k) = ecUnits/aval*(DivV - (gamma-1._RKIND)*GradRhoDotV)
!!$
!!$           enddo
!!$        enddo
!!$     enddo
!!$  endif


  ! call Model-specific local Jacobian routine
  if (Model == 10) then    ! standard Zeus-like model

     call gFLDProblem_LocalJac10(Erjac_Er, Erjac_ec, ecjac_Er, ecjac_ec, &
          time, Era, eca, eha, rhoa, vx, vy, vz, Nchem, Model, ESpectrum, &
          ProbType, DualEnergy, a, adot, PmC0, PmC1, PmC2, PmC3, PmC4,    &
          EmC0, EmC1, EmC2, EmC3, EmC4, gamma, DenUnits, VelUnits,        &
          ErUnits, ecUnits, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)

  else if (Model == 4) then   ! isothermal case-B Hydrogen recombination
     
     call gFLDProblem_LocalJac4(Erjac_Er, Erjac_ec, Erjac_HI, ecjac_Er, &
          ecjac_ec, ecjac_HI, HIjac_Er, HIjac_ec, HIjac_HI, time, Era,   &
          eca, n_HIa, eha, rhoa, vx, vy, vz, Nchem, HFrac, Model,        &
          ESpectrum, DualEnergy, a, adot, IsE, IsEsHI, IsEsHInu, gamma,  &
          NTempBins, TempStart, TempEnd, k1Tb, k2Tb, DenUnits, VelUnits, &
          ErUnits, ecUnits, NiUnits, ecScale, Nx, Ny, Nz, NGxl, NGxr,    &
          NGyl, NGyr, NGzl, NGzr, ier)

  else if (Model == 5) then   ! isothermal point-source emissivity

     call gFLDProblem_LocalJac5(Erjac_Er, Erjac_ec, Erjac_HI, ecjac_Er,   &
          ecjac_ec, ecjac_HI, HIjac_Er, HIjac_ec, HIjac_HI, Era, n_HIa,    &
          Nchem, Model, ESpectrum, IsE, IsEsHI, a, adot, ErUnits, NiUnits, &
          Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
     
  else if (Model == 2) then   ! case-A cosmological problems, with emissivity

     call gFLDProblem_LocalJac2(Erjac_Er, Erjac_ec, Erjac_HI, Erjac_HeI,     &
          Erjac_HeII, ecjac_Er, ecjac_ec, ecjac_HI, ecjac_HeI, ecjac_HeII,    &
          HIjac_Er, HIjac_ec, HIjac_HI, HIjac_HeI, HIjac_HeII, HeIjac_Er,     &
          HeIjac_ec, HeIjac_HI, HeIjac_HeI, HeIjac_HeII, HeIIjac_Er,          &
          HeIIjac_ec, HeIIjac_HI, HeIIjac_HeI, HeIIjac_HeII, time, Era, eca,  &
          n_HIa, n_HeIa, n_HeIIa, eha, rhoa, vx, vy, vz, Nchem, HFrac, Model, &
          ESpectrum, DualEnergy, a, adot, CompA, Comp_xray, Comp_temp, IsE,   &
          IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII, IsEsHeIInu, gamma,  &
          NTempBins, TempStart, TempEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb,  &
          ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb,     &
          reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, piHI, piHeI,      &
          piHeII, aUnits, DenUnits, VelUnits, ErUnits, ecUnits, NiUnits, Nx,  &
          Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
     
  else if (Model == 1) then   ! case-B cosmological problems

     call gFLDProblem_LocalJac1(Erjac_Er, Erjac_ec, Erjac_HI, Erjac_HeI,     &
          Erjac_HeII, ecjac_Er, ecjac_ec, ecjac_HI, ecjac_HeI, ecjac_HeII,    &
          HIjac_Er, HIjac_ec, HIjac_HI, HIjac_HeI, HIjac_HeII, HeIjac_Er,     &
          HeIjac_ec, HeIjac_HI, HeIjac_HeI, HeIjac_HeII, HeIIjac_Er,          &
          HeIIjac_ec, HeIIjac_HI, HeIIjac_HeI, HeIIjac_HeII, time, Era, eca,  &
          n_HIa, n_HeIa, n_HeIIa, eha, rhoa, vx, vy, vz, Nchem, HFrac, Model, &
          ESpectrum, DualEnergy, a, adot, CompA, Comp_xray, Comp_temp, IsE,   &
          IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII, IsEsHeIInu, gamma,  &
          NTempBins, TempStart, TempEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb,  &
          ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb,     &
          reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, piHI, piHeI,      &
          piHeII, aUnits, DenUnits, VelUnits, ErUnits, ecUnits, NiUnits, Nx,  &
          Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)

  !  we have not yet implemented analytical Jacobians for other models
  else

     write(0,*) 'gFLDProblem_LocalJac: Model =',Model,' not yet implemented!'

  endif ! Model

  return
end subroutine gFLDProblem_LocalJac
!=======================================================================





subroutine gFLDProblem_LocalJac1(Erjac_Er, Erjac_ec, Erjac_HI, Erjac_HeI, &
     Erjac_HeII, ecjac_Er, ecjac_ec, ecjac_HI, ecjac_HeI, ecjac_HeII,      &
     HIjac_Er, HIjac_ec, HIjac_HI, HIjac_HeI, HIjac_HeII, HeIjac_Er,       &
     HeIjac_ec, HeIjac_HI, HeIjac_HeI, HeIjac_HeII, HeIIjac_Er,            &
     HeIIjac_ec, HeIIjac_HI, HeIIjac_HeI, HeIIjac_HeII, time, Era, eca,    &
     n_HIa, n_HeIa, n_HeIIa, eha, rhoa, vx, vy, vz, Nchem, HFrac, Model,   &
     ESpectrum, DualEnergy, a, adot, CompA, Comp_xray, Comp_temp, IsE,     &
     IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII, IsEsHeIInu, gamma,    &
     NTempBins, TempStart, TempEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb,    &
     ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb,       &
     reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, piHI, piHeI,        &
     piHeII, aUnits, DenUnits, VelUnits, ErUnits, ecUnits, NiUnits, Nx,    &
     Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August, 2006
  !
  !  PURPOSE: Computes the spatially-local components of the chemical 
  !           Jacobian for the Gray FLD problem (Model 1).
  !
  !  INPUTS:
  !     time       - simulation time for evaluation
  !     Era        - Gray radiation energy density
  !     eca        - fluid energy correction array
  !     n_HIa      - number density of Hydrogen I species
  !     n_HeIa     - number density of Helium I species
  !     n_HeIIa    - number density of Helium II species
  !     eha        - total fluid energy array
  !     rhoa       - fluid density, assumed partitioned into either 
  !     vx         - fluid velocity (x-direction)
  !     vy         - fluid velocity (y-direction)
  !     vz         - fluid velocity (z-direction)
  !     Nchem      - number of chemical species (allowed: 0, 1, 3)
  !     HFrac      - percentage of mass composed of Hydrogen
  !     Model      - flag denoting physical model to use
  !     ESpectrum  - radiation spectrum choice
  !                       1 -> 1e5 black body spectrum
  !                       0 -> power law spectrum
  !                      -1 -> monochromatic 
  !     DualEnergy - flag denoting dual energy formalism
  !     a          - cosmological expansion parameter
  !     adot       - da/dt
  !     CompA      - Compton cooling coefficient 1 (multiplier)
  !     Comp_xray  - X-ray Compton heating coefficient
  !     Comp_temp  - X-ray Compton heating temperature 
  !     IsE        - int_{nu0_HI}^{inf} sigE dnu
  !     IsEsHI     - int_{nu0_HI}^{inf} sigE*sigHI dnu
  !     IsEsHInu   - int_{nu0_HI}^{inf} sigE*sigHI/nu dnu
  !     IsEsHeI    - int_{nu0_HeI}^{inf} sigE*sigHeI dnu
  !     IsEsHeInu  - int_{nu0_HeI}^{inf} sigE*sigHeI/nu dnu
  !     IsEsHeII   - int_{nu0_HeII}^{inf} sigE*sigHeII dnu
  !     IsEsHeIInu - int_{nu0_HeII}^{inf} sigE*sigHeII/nu dnu
  !     gamma      - ideal gas parameter (ratio of specific heats)
  !     *Units     - variable scaling constants
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: 
  !     Erjac_Er     - local radiation energy Jacobian w.r.t. Er
  !     Erjac_ec     - local radiation energy Jacobian w.r.t. ec
  !     Erjac_HI     - local radiation energy Jacobian w.r.t. nHI
  !     Erjac_HeI    - local radiation energy Jacobian w.r.t. nHeI
  !     Erjac_HeII   - local radiation energy Jacobian w.r.t. nHeII
  !     ecjac_Er     - local gas energy Jacobian w.r.t. Er
  !     ecjac_ec     - local gas energy Jacobian w.r.t. ec
  !     ecjac_HI     - local gas energy Jacobian w.r.t. nHI
  !     ecjac_HeI    - local gas energy Jacobian w.r.t. nHeI
  !     ecjac_HeII   - local gas energy Jacobian w.r.t. nHeII
  !     HIjac_Er     - local HI chemistry Jacobian w.r.t. Er
  !     HIjac_ec     - local HI chemistry Jacobian w.r.t. ec
  !     HIjac_HI     - local HI chemistry Jacobian w.r.t. nHI
  !     HIjac_HeI    - local HI chemistry Jacobian w.r.t. nHeI
  !     HIjac_HeII   - local HI chemistry Jacobian w.r.t. nHeII
  !     HeIjac_Er    - local HeI chemistry Jacobian w.r.t. Er
  !     HeIjac_ec    - local HeI chemistry Jacobian w.r.t. ec
  !     HeIjac_HI    - local HeI chemistry Jacobian w.r.t. nHI
  !     HeIjac_HeI   - local HeI chemistry Jacobian w.r.t. nHeI
  !     HeIjac_HeII  - local HeI chemistry Jacobian w.r.t. nHeII
  !     HeIIjac_Er   - local HeII chemistry Jacobian w.r.t. Er
  !     HeIIjac_ec   - local HeII chemistry Jacobian w.r.t. ec
  !     HeIIjac_HI   - local HeII chemistry Jacobian w.r.t. nHI
  !     HeIIjac_HeI  - local HeII chemistry Jacobian w.r.t. nHeI
  !     HeIIjac_HeII - local HeII chemistry Jacobian w.r.t. nHeII
  !     ier          - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in) :: Nchem, Model, ESpectrum, DualEnergy
  INTG_PREC, intent(in) :: NTempBins
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in) :: a, adot
  R_PREC,    intent(in) :: time, HFrac, gamma, IsE, IsEsHI, IsEsHInu
  R_PREC,    intent(in) :: IsEsHeI, IsEsHeInu, IsEsHeII, IsEsHeIInu
  R_PREC,    intent(in) :: TempStart, TempEnd, piHI, piHeI, piHeII
  R_PREC,    intent(in) :: CompA, Comp_xray, Comp_temp
  R_PREC,    intent(in) :: aUnits, DenUnits, VelUnits, ErUnits, ecUnits, NiUnits
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) :: &
       Era, eca, n_HIa, n_HeIa, n_HeIIa, eha, rhoa, vx, vy, vz
  R_PREC, dimension(NTempBins), intent(in) :: k1Tb, k2Tb, k3Tb, k4Tb,      &
       k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, &
       ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(out) :: &
       Erjac_Er, Erjac_ec, Erjac_HI, Erjac_HeI, Erjac_HeII,    &
       ecjac_Er, ecjac_ec, ecjac_HI, ecjac_HeI, ecjac_HeII,       &
       HIjac_Er, HIjac_ec, HIjac_HI, HIjac_HeI, HIjac_HeII,       &
       HeIjac_Er, HeIjac_ec, HeIjac_HI, HeIjac_HeI, HeIjac_HeII,  &
       HeIIjac_Er, HeIIjac_ec, HeIIjac_HI, HeIIjac_HeI, HeIIjac_HeII

  !--------------
  ! locals
  INTG_PREC :: i, j, k, Tidx, Tidxp
  R_PREC :: lTempS, lTempE, dlTemp, lTemp, Tl, Tr, Tfac, Tfac2
  R_PREC :: k1, k2, k3, k4, k5, k6
  R_PREC :: dk1_dT, dk2_dT, dk3_dT, dk4_dT, dk5_dT, dk6_dT
  R_PREC :: aval, grey
  R_PREC :: T, dT_dec, dT_dHI, dT_dHeI, dT_dHeII, min_temp
  R_PREC :: lamT, dlamT_dT
  R_PREC :: pi, c, hp, kb, StBz, eint, Er, rho, zr, afac, mp, KEconst
  R_PREC :: mol_weight, dmu_dHI, dmu_dHeI, dmu_dHeII
  R_PREC :: comp1, comp2
  R_PREC :: HIconst, HeIconst, HeIIconst, nu0_HI, nu0_HeI, nu0_HeII
  R_PREC :: nH, nHI, dnHI_dHI, nHII, dnHII_dHI
  R_PREC :: ne, dne_dHI, dne_dHeI, dne_dHeII
  R_PREC :: nHe, nHeI, dnHeI_dHeI, nHeII, dnHeII_dHeII
  R_PREC :: nHeIII, dnHeIII_dHeI, dnHeIII_dHeII 
  R_PREC :: kappaE, dkappaE_dHI, dkappaE_dHeI, dkappaE_dHeII
  R_PREC :: eta, deta_dec, deta_dHI, deta_dHeI, deta_dHeII
  R_PREC :: brem, dbrem_dT, ceHI, dceHI_dT, ceHeI, dceHeI_dT
  R_PREC :: ceHeII, dceHeII_dT, ciHI, dciHI_dT, ciHeI, dciHeI_dT
  R_PREC :: ciHeIS, dciHeIS_dT, ciHeII, dciHeII_dT, reHII, dreHII_dT
  R_PREC :: reHeII1, dreHeII1_dT, reHeII2, dreHeII2_dT, reHeIII, dreHeIII_dT
  R_PREC :: dG_dEr, dG_dHI, dG_dHeI, dG_dHeII
  R_PREC :: dLambda_dec, dLambda_dHI, dLambda_dHeI, dLambda_dHeII
  R_PREC :: G_HI, dGHI_dEr, G_HeI, dGHeI_dEr, G_HeII, dGHeII_dEr

  !=======================================================================

  ! set flag to success
  ier = 1

  ! initialize outputs to have all zero values
  Erjac_Er = 0._RKIND
  Erjac_ec = 0._RKIND
  ecjac_Er = 0._RKIND
  if (Nchem == 1) then
     Erjac_HI = 0._RKIND

     ecjac_HI = 0._RKIND

     HIjac_Er = 0._RKIND
     HIjac_ec = 0._RKIND
     HIjac_HI = 0._RKIND
  endif
  if (Nchem == 3) then
     Erjac_HI   = 0._RKIND
     Erjac_HeI  = 0._RKIND
     Erjac_HeII = 0._RKIND

     ecjac_HI   = 0._RKIND
     ecjac_HeI  = 0._RKIND
     ecjac_HeII = 0._RKIND

     HIjac_Er   = 0._RKIND
     HIjac_ec   = 0._RKIND
     HIjac_HI   = 0._RKIND
     HIjac_HeI  = 0._RKIND
     HIjac_HeII = 0._RKIND

     HeIjac_Er   = 0._RKIND
     HeIjac_ec   = 0._RKIND
     HeIjac_HI   = 0._RKIND
     HeIjac_HeI  = 0._RKIND
     HeIjac_HeII = 0._RKIND

     HeIIjac_Er   = 0._RKIND
     HeIIjac_ec   = 0._RKIND
     HeIIjac_HI   = 0._RKIND
     HeIIjac_HeI  = 0._RKIND
     HeIIjac_HeII = 0._RKIND
  endif


  ! initialize constants
  aval = a*aunits
  pi = pi_val
  afac = adot/a                ! adot/a
  mp = mass_h                  ! mass of a proton [g]
  c  = c_light                 ! speed of light    [cm/s]
  hp = hplanck                 ! Planck's constant [ergs*s]
  kb = kboltz                  ! boltzmann constant [erg/K]
  StBz  = 5.6704e-5_RKIND      ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  min_temp = 1.d-1             ! minimum temperature [K]
  nu0_HI = 13.6*ev2erg/hp      ! ionization frequency of HI
  nu0_HeI = 24.6*ev2erg/hp     ! ionization frequency of HeI
  nu0_HeII = 54.4*ev2erg/hp    ! ionization frequency of HeII
  zr = 1._RKIND/aval - 1._RKIND        ! cosmological redshift
  if (aval .ne. 1._RKIND) then        ! Compton cooling coefficients
     Comp1 = CompA*(1._RKIND + zr)**4
     Comp2 = 2.73_RKIND*(1._RKIND + zr)
  else
     Comp1 = 0._RKIND
     Comp2 = 0._RKIND
  endif
  if (DualEnergy == 1) then
     KEconst = 0._RKIND
  else
     KEconst = 0.5_RKIND
  endif
  grey = 1._RKIND
  if (ESpectrum == -1)  grey = 0._RKIND

  !   lookup table constants
  lTempS = log(TempStart)
  lTempE = log(TempEnd)
  dlTemp = (lTempE - lTempS)/(1._RKIND*NTempBins - 1._RKIND)

  ! set shortcut values
  HIconst   = c*(IsEsHI - nu0_HI*IsEsHInu)/IsE
  HeIconst  = c*(IsEsHeI - nu0_HeI*IsEsHeInu)/IsE
  HeIIconst = c*(IsEsHeII - nu0_HeII*IsEsHeIInu)/IsE

  if (Model /= 1) then
     write(0,*) 'gFLDProblem_LocalJac1: incorrect Model =',Model
     return
  endif

  ! Hydrogen only case
  if (Nchem == 1) then
     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1

              ! set shortcut values for this spatial location
              ! converting densities from comoving to proper
              Er  = Era(i,j,k)*ErUnits
              rho = rhoa(i,j,k)*DenUnits
              nHI = n_HIa(i,j,k)*NiUnits
              nH = Hfrac*rho/mp
              nHII = max(nH - nHI, 0._RKIND)
              ne  = nHII
              dnHI_dHI = NiUnits
              dnHII_dHI = -dnHI_dHI
              dne_dHI = dnHII_dHI

              ! compute temperature and shortcuts
              !!!! This assumes ALL density is H, otherwise mol_weight it wrong !!!!
              eint = eca(i,j,k)*ecUnits + VelUnits*VelUnits*(eha(i,j,k)  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              mol_weight = rho/mp/(nHI + nHII + ne)
              dmu_dHI = -rho/mp*(dnHI_dHI + dnHII_dHI + dne_dHI) &
                   / (nHI + nHII + ne)**2
              T = (gamma-1._RKIND)*mol_weight*mp*eint/kb
              T = max(T,min_temp)
              dT_dec = ecUnits*(gamma-1._RKIND)*mol_weight*mp/kb
              dT_dHI = (gamma-1._RKIND)*mp*eint/kb*dmu_dHI
              lamT = 3.15614e5_RKIND/T
              dlamT_dT = -3.15614e5_RKIND/T/T

              ! look up rates, derivatives
              lTemp = min(max(log(T), lTempS), lTempE)
              Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
              Tidxp = Tidx+1
              Tl = lTempS + (Tidx-1)*dlTemp
              Tr = lTempS +  Tidx*dlTemp
              Tfac = (lTemp - Tl)/(Tr - Tl)
              Tfac2 = 1._RKIND/(Tr - Tl)/T
              k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
              dk1_dT = (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac2
              ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
              dceHI_dT = (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac2
              ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
              dciHI_dT = (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac2
              reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
              dreHII_dT = (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac2
              brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac
              dbrem_dT = (bremTb(Tidxp) - bremTb(Tidx))*Tfac2

              ! compute opacities and derivatives
              kappaE = nHI*IsEsHI/IsE
              dkappaE_dHI = IsEsHI/IsE*dnHI_dHI

              ! compute case B Hydrogen recombination coefficient 
              ! [Hui & Gnedin, 1997: RI^B_{HII}]
              !   (still need this because table holds case A coefficient!)
              k2 = 2.753e-14_RKIND*lamT**(1.5_RKIND) *                         &
                   (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-2.242_RKIND)
              dk2_dT = 2.753e-14_RKIND*dlamT_dT*sqrt(lamT)*(               &
                   1.5_RKIND*(1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-2.242_RKIND) &
                   - lamT*(2.242_RKIND*0.407_RKIND/2.74_RKIND)*                  &
                   (lamT/2.74_RKIND)**(-0.593_RKIND)*                        &
                   (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-3.242_RKIND))

              ! compute Hydrogen photoionization rate & derivs
              G_HI = c*Er/hp*IsEsHInu/IsE
              dGHI_dEr = c/hp*IsEsHInu/IsE*ErUnits

              ! compute fluid cooling rate.  Terms (in order):
              !    Collisional Excitations
              !    Collisional Ionizations
              !    Recombinations
              !    Compton cooling or heating 
              !    X-ray Compton heating
              !    Bremsstrahlung
              !  Lambda = ne/rho*(ceHI*nHI + ciHI*nHI
              !     + reHII*nHII + Comp1*(T-Comp2)
              !     + Comp_xray*(T-Comp_temp) + brem*nHII)
              dLambda_dec = ne/rho*dT_dec*( nHI*(dceHI_dT + dciHI_dT) &
                   + nHII*(dreHII_dT + dbrem_dT) + Comp1 + Comp_xray )
              dLambda_dHI = dne_dHI/rho*((ceHI+ciHI)*nHI   &
                     + (reHII+brem)*nHII                   &
                     + Comp1*(T-Comp2)                     &
                     + Comp_xray*(T-Comp_temp))            &
                   + ne/rho*( (ceHI + ciHI)*dnHI_dHI       &
                     + (dceHI_dT + dciHI_dT)*dT_dHI*nHI    &
                     + (reHII + brem)*dnHII_dHI            &
                     + (dreHII_dT + dbrem_dT)*dT_dHI*nHII  &
                     + (Comp1 + Comp_xray)*dT_dHI)

              ! compute fluid heating rate
              !   G = nHI/rho*Er*HIconst
              dG_dEr = ErUnits*nHI/rho*HIconst
              dG_dHI = dnHI_dHI/rho*Er*HIconst

              ! compute emissivity 
              ! [WE HAVE NO FORMULA FOR ETA IN CASE-B CASE!!]
              eta      = 0._RKIND
              deta_dec = 0._RKIND
              deta_dHI = 0._RKIND

              ! put it all together
              !   rhs_Er = (src_Er + 4._RKIND*pi*eta 
              !            - (afac*grey + c*kappaE)*Er)/ErUnits
              Erjac_Er(i,j,k) = -afac*grey - c*kappaE
              Erjac_ec(i,j,k) =  4._RKIND*pi*deta_dec/ErUnits
              Erjac_HI(i,j,k) = (4._RKIND*pi*deta_dHI - c*Er*dkappaE_dHI)/ErUnits

              !   rhs_HI = (src_HI + k2*ne*nHII - nHI*(k1*ne + G_HI))/NiUnits
              HIjac_Er(i,j,k) = -nHI*dGHI_dEr/NiUnits
              HIjac_ec(i,j,k) = ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dec/NiUnits
              HIjac_HI(i,j,k) = (dk2_dT*dT_dHI*ne*nHII + k2*dne_dHI*nHII  &
                   + k2*ne*dnHII_dHI - dnHI_dHI*(k1*ne + G_HI)            &
                   - nHI*(dk1_dT*dT_dHI*ne + k1*dne_dHI))/NiUnits

              !   rhs_ec = (src_ec + G - Lambda - 2*afac*ec)/ecUnits
              ecjac_ec(i,j,k) = (ecjac_ec(i,j,k) - 2._RKIND*afac*ecUnits &
                   - dLambda_dec)/ecUnits
              ecjac_Er(i,j,k) = dG_dEr/ecUnits
              ecjac_HI(i,j,k) = (dG_dHI - dLambda_dHI)/ecUnits

           enddo
        enddo
     enddo

     ! Hydrogen + Helium case
  else if (Nchem == 3) then

     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1

              ! set shortcut values for this spatial location
              ! converting densities from comoving to proper
              Er  = Era(i,j,k)*ErUnits
              rho = rhoa(i,j,k)*DenUnits
              nHI = n_HIa(i,j,k)*NiUnits
              nH = Hfrac*rho/mp
              nHII = max(nH - nHI, 0._RKIND)
              nHe    = (1._RKIND-HFrac)*rho/4._RKIND/mp
              nHeI   = n_HeIa(i,j,k)*NiUnits/4._RKIND
              nHeII  = n_HeIIa(i,j,k)*NiUnits/4._RKIND
              nHeIII = max(nHe - nHeI - nHeII, 0._RKIND)
              ne     = nHII + 0.25_RKIND*nHeII + 0.5_RKIND*nHeIII
              dnHI_dHI = NiUnits
              dnHII_dHI = -dnHI_dHI
              dnHeI_dHeI = NiUnits/4._RKIND
              dnHeII_dHeII = NiUnits/4._RKIND
              dnHeIII_dHeI = -dnHeI_dHeI
              dnHeIII_dHeII = -dnHeII_dHeII
              dne_dHI = dnHII_dHI
              dne_dHeI = 0.5_RKIND*dnHeIII_dHeI
              dne_dHeII = 0.25_RKIND*dnHeII_dHeII + 0.5_RKIND*dnHeIII_dHeII

              ! compute temperature and shortcuts
              eint = eca(i,j,k)*ecUnits + VelUnits*VelUnits*(eha(i,j,k)  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              mol_weight = rho/mp/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)
              dmu_dHI = -rho/mp/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)**2 &
                   * (dnHI_dHI+dnHII_dHI+dne_dHI)
              dmu_dHeI = -rho/mp/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)**2 &
                   * (0.25_RKIND*(dnHeI_dHeI+dnHeIII_dHeI) + dne_dHeI)
              dmu_dHeII = -rho/mp/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)**2 &
                   * (0.25_RKIND*(dnHeII_dHeII+dnHeIII_dHeII) + dne_dHeII)
              T = (gamma-1._RKIND)*mol_weight*mp*eint/kb
              T = max(T,min_temp)
              dT_dec = ecUnits*(gamma-1._RKIND)*mol_weight*mp/kb
              dT_dHI = (gamma-1._RKIND)*mp*eint/kb*dmu_dHI
              dT_dHeI = (gamma-1._RKIND)*mp*eint/kb*dmu_dHeI
              dT_dHeII = (gamma-1._RKIND)*mp*eint/kb*dmu_dHeII
              lamT = 3.15614e5_RKIND/T
              dlamT_dT = -3.15614e5_RKIND/T/T

              ! look up rates
              lTemp = min(max(log(T), lTempS), lTempE)
              Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
              Tidxp = Tidx+1
              Tl = lTempS + (Tidx-1)*dlTemp
              Tr = lTempS +  Tidx*dlTemp
              Tfac = (lTemp - Tl)/(Tr - Tl)
              Tfac2 = 1._RKIND/(Tr - Tl)/T
              k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
              dk1_dT = (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac2
              k2 = k2Tb(Tidx) + (k2Tb(Tidxp) - k2Tb(Tidx))*Tfac
              dk2_dT = (k2Tb(Tidxp) - k2Tb(Tidx))*Tfac2
              k3 = k3Tb(Tidx) + (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac
              dk3_dT = (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac2
              k4 = k4Tb(Tidx) + (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac
              dk4_dT = (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac2
              k5 = k5Tb(Tidx) + (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac
              dk5_dT = (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac2
              k6 = k6Tb(Tidx) + (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac
              dk6_dT = (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac2
              ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
              dceHI_dT = (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac2
              ceHeI = ceHeITb(Tidx) + (ceHeITb(Tidxp) - ceHeITb(Tidx))*Tfac
              dceHeI_dT = (ceHeITb(Tidxp) - ceHeITb(Tidx))*Tfac2
              ceHeII = ceHeIITb(Tidx) + (ceHeIITb(Tidxp) - ceHeIITb(Tidx))*Tfac
              dceHeII_dT = (ceHeIITb(Tidxp) - ceHeIITb(Tidx))*Tfac2
              ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
              dciHI_dT = (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac2
              ciHeI = ciHeITb(Tidx) + (ciHeITb(Tidxp) - ciHeITb(Tidx))*Tfac
              dciHeI_dT = (ciHeITb(Tidxp) - ciHeITb(Tidx))*Tfac2
              ciHeII = ciHeIITb(Tidx) + (ciHeIITb(Tidxp) - ciHeIITb(Tidx))*Tfac
              dciHeII_dT = (ciHeIITb(Tidxp) - ciHeIITb(Tidx))*Tfac2
              ciHeIS = ciHeISTb(Tidx) + (ciHeISTb(Tidxp) - ciHeISTb(Tidx))*Tfac
              dciHeIS_dT = (ciHeISTb(Tidxp) - ciHeISTb(Tidx))*Tfac2
              reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
              dreHII_dT = (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac2
              reHeII1 = reHeII1Tb(Tidx) + (reHeII1Tb(Tidxp) - reHeII1Tb(Tidx))*Tfac
              dreHeII1_dT = (reHeII1Tb(Tidxp) - reHeII1Tb(Tidx))*Tfac2
              reHeII2 = reHeII2Tb(Tidx) + (reHeII2Tb(Tidxp) - reHeII2Tb(Tidx))*Tfac
              dreHeII2_dT = (reHeII2Tb(Tidxp) - reHeII2Tb(Tidx))*Tfac2
              reHeIII = reHeIIITb(Tidx) + (reHeIIITb(Tidxp) - reHeIIITb(Tidx))*Tfac
              dreHeIII_dT = (reHeIIITb(Tidxp) - reHeIIITb(Tidx))*Tfac2
              brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac
              dbrem_dT = (bremTb(Tidxp) - bremTb(Tidx))*Tfac2

              ! compute opacities and derivatives
              kappaE = nHI*IsEsHI/IsE + nHeI*IsEsHeI/IsE + nHeII*IsEsHeII/IsE
              dkappaE_dHI = IsEsHI/IsE*dnHI_dHI
              dkappaE_dHeI = IsEsHeI/IsE*dnHeI_dHeI
              dkappaE_dHeII = IsEsHeII/IsE*dnHeII_dHeII

              ! compute case B Hydrogen recombination coefficient 
              ! [Hui & Gnedin, 1997: RI^B_{HII}]
              !   (still need this because table holds case A coefficient!)
              k2 = 2.753e-14_RKIND*lamT**(1.5_RKIND) *                         &
                   (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-2.242_RKIND)
              dk2_dT = 2.753e-14_RKIND*dlamT_dT*sqrt(lamT)*(               &
                   1.5_RKIND*(1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-2.242_RKIND) &
                   - lamT*(2.242_RKIND*0.407_RKIND/2.74_RKIND)*                  &
                   (lamT/2.74_RKIND)**(-0.593_RKIND)*                        &
                   (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-3.242_RKIND))

              ! compute Hydrogen photoionization rate & derivs
              G_HI = c*Er/hp*IsEsHInu/IsE
              dGHI_dEr = c/hp*IsEsHInu/IsE*ErUnits

              ! compute Helium photoionization rates & derivs
              G_HeI  = c*Er/hp*IsEsHeInu/IsE
              G_HeII = c*Er/hp*IsEsHeIInu/IsE
              dGHeI_dEr = c/hp*IsEsHeInu/IsE*ErUnits
              dGHeII_dEr = c/hp*IsEsHeIInu/IsE*ErUnits

              ! compute fluid cooling rate.  Terms (in order):
              !    Collisional Excitations
              !    Collisional Ionizations
              !    Recombinations
              !    Compton cooling or heating 
              !    X-ray Compton heating
              !    Bremsstrahlung
              !  Lambda = ne/rho*(ceHI*nHI + ciHI*nHI + reHII*nHII 
              !     + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII
              !     + reHeII1 + reHeII2 + brem)*nHeII/4._RKIND
              !     + (ciHeI*nHeI + reHeIII*nHeIII)/4._RKIND
              !     + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp)
              !     + brem*(nHII+nHeIII))
              dLambda_dec = ne/rho*dT_dec*(dceHI_dT*nHI + dciHI_dT*nHI &
                   + (dceHeI_dT*ne + dciHeIS_dT*ne + dceHeII_dT        &
                   + dciHeII_dT + dreHeII1_dT + dreHeII2_dT            &
                   + dbrem_dT)*nHeII/4._RKIND + dciHeI_dT*nHeI/4._RKIND        &
                   + dreHeIII_dT*nHeIII/4._RKIND + dreHII_dT*nHII          &
                   + Comp1 + Comp_xray + dbrem_dT*(nHII+nHeIII))
              dLambda_dHI = ne/rho*(                                 &
                   (ceHI + ciHI)*dnHI_dHI + (reHII+brem)*dnHII_dHI + &
                   dT_dHI*(dceHI_dT*nHI + dciHI_dT*nHI               &
                   + (dceHeI_dT*ne + dciHeIS_dT*ne + dceHeII_dT      &
                   + dciHeII_dT + dreHeII1_dT + dreHeII2_dT          &
                   + dbrem_dT)*nHeII/4._RKIND + dciHeI_dT*nHeI/4._RKIND      &
                   + dreHeIII_dT*nHeIII/4._RKIND + dreHII_dT*nHII        &
                   + Comp1 + Comp_xray + dbrem_dT*(nHII+nHeIII))) +  &
                   dne_dHI/rho*(ceHI*nHI + reHII*nHII + ciHI*nHI     &
                   + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII         &
                   + reHeII1 + reHeII2 + brem)*nHeII/4._RKIND            &
                   + (reHeIII*nHeIII + ciHeI*nHeI)/4._RKIND              &
                   + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp)       &
                   + brem*(nHII+nHeIII))
              dLambda_dHeI = ne/rho*(                                &
                   + (ciHeI*dnHeI_dHeI + reHeIII*dnHeIII_dHeI)/4._RKIND  &
                   + brem*dnHeIII_dHeI +                             &
                   dT_dHeI*(dceHI_dT*nHI + dciHI_dT*nHI              &
                   + (dceHeI_dT*ne + dciHeIS_dT*ne + dceHeII_dT      &
                   + dciHeII_dT + dreHeII1_dT + dreHeII2_dT          &
                   + dbrem_dT)*nHeII/4._RKIND + dciHeI_dT*nHeI/4._RKIND      &
                   + dreHeIII_dT*nHeIII/4._RKIND + dreHII_dT*nHII        &
                   + Comp1 + Comp_xray + dbrem_dT*(nHII+nHeIII))) +  &
                   dne_dHeI/rho*(ceHI*nHI + reHII*nHII + ciHI*nHI    &
                   + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII         &
                   + reHeII1 + reHeII2 + brem)*nHeII/4._RKIND            &
                   + (reHeIII*nHeIII + ciHeI*nHeI)/4._RKIND              &
                   + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp)       &
                   + brem*(nHII+nHeIII))
              dLambda_dHeII = ne/rho*(                               &
                   + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII         &
                   + reHeII1 + reHeII2 + brem)*dnHeII_dHeII/4._RKIND     &
                   + (reHeIII/4._RKIND + brem)*dnHeIII_dHeII +           &
                   dT_dHeII*(dceHI_dT*nHI + dciHI_dT*nHI             &
                   + (dceHeI_dT*ne + dciHeIS_dT*ne + dceHeII_dT      &
                   + dciHeII_dT + dreHeII1_dT + dreHeII2_dT          &
                   + dbrem_dT)*nHeII/4._RKIND + dciHeI_dT*nHeI/4._RKIND      &
                   + dreHeIII_dT*nHeIII/4._RKIND + dreHII_dT*nHII        &
                   + Comp1 + Comp_xray + dbrem_dT*(nHII+nHeIII))) +  &
                   dne_dHeII/rho*(ceHI*nHI + reHII*nHII + ciHI*nHI   &
                   + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII         &
                   + reHeII1 + reHeII2 + brem)*nHeII/4._RKIND            &
                   + (reHeIII*nHeIII + ciHeI*nHeI)/4._RKIND              &
                   + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp)       &
                   + brem*(nHII+nHeIII))

              ! compute fluid heating rate
              !   G = Er/rho*(nHI*HIconst + nHeI*HeIconst + nHeII*HeIIconst)
              dG_dEr = ErUnits/rho*(nHI*HIconst + nHeI*HeIconst + nHeII*HeIIconst)
              dG_dHI = Er/rho*HIconst*dnHI_dHI
              dG_dHeI = Er/rho*HeIconst*dnHeI_dHeI
              dG_dHeII = Er/rho*HeIIconst*dnHeII_dHeII

              ! compute emissivity 
              ! [WE HAVE NO FORMULA FOR ETA IN CASE-B CASE!!]
              eta        = 0._RKIND
              deta_dec   = 0._RKIND
              deta_dHI   = 0._RKIND
              deta_dHeI  = 0._RKIND
              deta_dHeII = 0._RKIND

              ! put it all together
              !   rhs_Er = (src_Er + 4._RKIND*pi*eta 
              !             - (afac*grey + c*kappaE)*Er)/ErUnits
              Erjac_Er(i,j,k)   = -(afac*grey + c*kappaE)
              Erjac_ec(i,j,k)   =  (4._RKIND*pi*deta_dec)/ErUnits
              Erjac_HI(i,j,k)   =  (4._RKIND*pi*deta_dHI-c*Er*dkappaE_dHI)/ErUnits
              Erjac_HeI(i,j,k)  =  (4._RKIND*pi*deta_dHeI-c*Er*dkappaE_dHeI)/ErUnits
              Erjac_HeII(i,j,k) =  (4._RKIND*pi*deta_dHeII-c*Er*dkappaE_dHeII)/ErUnits

              !   rhs_ec = (src_ec + G - Lambda - 2._RKIND*afac*ec)/ecUnits
              ecjac_ec(i,j,k) = (ecJac_ec(i,j,k) - 2._RKIND*afac*ecUnits &
                   - dLambda_dec)/ecUnits
              ecjac_Er(i,j,k) = dG_dEr/ecUnits
              ecjac_HI(i,j,k) = (dG_dHI - dLambda_dHI)/ecUnits
              ecjac_HeI(i,j,k) = (dG_dHeI - dLambda_dHeI)/ecUnits
              ecjac_HeII(i,j,k) = (dG_dHeII - dLambda_dHeII)/ecUnits

              !   rhs_HI = (src_HI + ne*(k2*nHII - k1*nHI) - nHI*G_HI)/NiUnits
              HIjac_Er(i,j,k) = -nHI*dGHI_dEr/NiUnits
              HIjac_ec(i,j,k) = ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dec/NiUnits
              HIjac_HI(i,j,k) = (ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dHI  &
                   + (k2*nHII - k1*nHI)*dne_dHI + k2*ne*dnHII_dHI      &
                   - dnHI_dHI*k1*ne - dnHI_dHI*G_HI)/NiUnits
              HIjac_HeI = (ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dHeI       &
                   + (k2*nHII - k1*nHI)*dne_dHeI)/NiUnits
              HIjac_HeII = (ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dHeII     &
                   + (k2*nHII - k1*nHI)*dne_dHeII)/NiUnits

              !  rhs_HeI = (src_HeI - nHeI*G_HeI
              !            + ne*(k4*nHeII - k3*nHeI))/NiUnits
              HeIjac_Er(i,j,k) = -nHeI*dGHeI_dEr/NiUnits
              HeIjac_ec(i,j,k) = ne*(nHeII*dk4_dT - nHeI*dk3_dT)           &
                   *dT_dec/NiUnits
              HeIjac_HI(i,j,k) = (ne*(nHeII*dk4_dT - nHeI*dk3_dT)*dT_dHI   &
                   + (k4*nHeII - k3*nHeI)*dne_dHI)/NiUnits
              HeIjac_HeI(i,j,k) = (ne*(nHeII*dk4_dT - nHeI*dk3_dT)*dT_dHeI &
                   + (k4*nHeII - k3*nHeI)*dne_dHeI - ne*k3*dnHeI_dHeI      &
                   - G_HeI*dnHeI_dHeI)/NiUnits
              HeIjac_HeII(i,j,k) = (ne*(nHeII*dk4_dT-nHeI*dk3_dT)*dT_dHeII &
                   + (k4*nHeII - k3*nHeI)*dne_dHeII                        &
                   + ne*k4*dnHeII_dHeII)/NiUnits

              ! rhs_HeII = (src_HeII + nHeI*G_HeI - nHeII*G_HeII
              !            + ne*(k3*nHeI + k6*nHeIII 
              !                  - k4*nHeII + k5*nHeII))/NiUnits
              HeIIjac_Er(i,j,k) = -nHeII*dGHeII_dEr*4._RKIND/NiUnits
              HeIIjac_ec(i,j,k) = ne*(nHeI*dk3_dT + nHeIII*dk6_dT    &
                   - nHeII*dk4_dT + nHeII*dk5_dT)*dT_dec/NiUnits
              HeIIjac_HI(i,j,k) = (ne*(nHeI*dk3_dT + nHeIII*dk6_dT   &
                   - nHeII*dk4_dT + nHeII*dk5_dT)*dT_dHI             &
                   + (k3*nHeI + k6*nHeIII - k4*nHeII                 &
                   + k5*nHeII)*dne_dHI)/NiUnits
              HeIIjac_HeI(i,j,k) = (ne*(nHeI*dk3_dT + nHeIII*dk6_dT  &
                   - nHeII*dk4_dT + nHeII*dk5_dT)*dT_dHeI            &
                   + (k3*nHeI + k6*nHeIII - k4*nHeII                 &
                   + ne*(k3*dnHeI_dHeI + k6*dnHeIII_dHeI)            &
                   + k5*nHeII)*dne_dHeI)/NiUnits
              HeIIjac_HeII(i,j,k) = (ne*(nHeI*dk3_dT + nHeIII*dk6_dT &
                   - nHeII*dk4_dT + nHeII*dk5_dT)*dT_dHeII           &
                   + (k3*nHeI + k6*nHeIII - k4*nHeII                 &
                   + k5*nHeII)*dne_dHeII                             &
                   + ne*(k6*dnHeIII_dHeII                            &
                   - (k4-k5)*dnHeII_dHeII)                           &
                   - G_HeII*dnHeII_dHeII)/NiUnits

           enddo
        enddo
     enddo

  else
     write(0,*) 'gFLDProblem_LocalJac1: illegal Nchem =',Nchem
  endif  ! Nchem

  return
end subroutine gFLDProblem_LocalJac1
!=======================================================================





subroutine gFLDProblem_LocalJac2(Erjac_Er, Erjac_ec, Erjac_HI, Erjac_HeI,   &
     Erjac_HeII, ecjac_Er, ecjac_ec, ecjac_HI, ecjac_HeI, ecjac_HeII,       &
     HIjac_Er, HIjac_ec, HIjac_HI, HIjac_HeI, HIjac_HeII, HeIjac_Er,        &
     HeIjac_ec, HeIjac_HI, HeIjac_HeI, HeIjac_HeII, HeIIjac_Er, HeIIjac_ec, &
     HeIIjac_HI, HeIIjac_HeI, HeIIjac_HeII, time, Era, eca, n_HIa, n_HeIa,  &
     n_HeIIa, eha, rhoa, vx, vy, vz, Nchem, HFrac, Model, ESpectrum,        &
     DualEnergy, a, adot, CompA, Comp_xray, Comp_temp, IsE, IsEsHI,         &
     IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII, IsEsHeIInu, gamma, NTempBins,  &
     TempStart, TempEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb,        &
     ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb, reHIITb,       &
     reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, piHI, piHeI, piHeII, aUnits,  &
     DenUnits, VelUnits, ErUnits, ecUnits, NiUnits, Nx, Ny, Nz, NGxl, NGxr, &
     NGyl, NGyr, NGzl, NGzr, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August, 2006
  !
  !  PURPOSE: Computes the spatially-local components of the chemical 
  !           Jacobian for the Gray FLD problem (Model 2).
  !
  !  INPUTS:
  !     time       - simulation time for evaluation
  !     Era        - Gray radiation energy density
  !     eca        - fluid energy correction array
  !     n_HIa      - number density of Hydrogen I species
  !     n_HeIa     - number density of Helium I species
  !     n_HeIIa    - number density of Helium II species
  !     eha        - total fluid energy array
  !     rhoa       - fluid density, assumed partitioned into either 
  !     vx         - fluid velocity (x-direction)
  !     vy         - fluid velocity (y-direction)
  !     vz         - fluid velocity (z-direction)
  !     Nchem      - number of chemical species (allowed: 0, 1, 3)
  !     HFrac      - percentage of mass composed of Hydrogen
  !     Model      - flag denoting physical model to use
  !     ESpectrum  - radiation spectrum choice
  !                       1 -> 1e5 black body spectrum
  !                       0 -> power law spectrum
  !                      -1 -> monochromatic 
  !     DualEnergy - flag denoting dual energy formalism
  !     a          - cosmological expansion parameter
  !     adot       - da/dt
  !     CompA      - Compton cooling coefficient 1 (multiplier)
  !     Comp_xray  - X-ray Compton heating coefficient
  !     Comp_temp  - X-ray Compton heating temperature 
  !     IsE        - int_{nu0_HI}^{inf} sigE dnu
  !     IsEsHI     - int_{nu0_HI}^{inf} sigE*sigHI dnu
  !     IsEsHInu   - int_{nu0_HI}^{inf} sigE*sigHI/nu dnu
  !     IsEsHeI    - int_{nu0_HeI}^{inf} sigE*sigHeI dnu
  !     IsEsHeInu  - int_{nu0_HeI}^{inf} sigE*sigHeI/nu dnu
  !     IsEsHeII   - int_{nu0_HeII}^{inf} sigE*sigHeII dnu
  !     IsEsHeIInu - int_{nu0_HeII}^{inf} sigE*sigHeII/nu dnu
  !     gamma      - ideal gas parameter (ratio of specific heats)
  !     *Units     - variable scaling constants
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: 
  !     Erjac_Er     - local radiation energy Jacobian w.r.t. Er
  !     Erjac_ec     - local radiation energy Jacobian w.r.t. ec
  !     Erjac_HI     - local radiation energy Jacobian w.r.t. nHI
  !     Erjac_HeI    - local radiation energy Jacobian w.r.t. nHeI
  !     Erjac_HeII   - local radiation energy Jacobian w.r.t. nHeII
  !     ecjac_Er     - local gas energy Jacobian w.r.t. Er
  !     ecjac_ec     - local gas energy Jacobian w.r.t. ec
  !     ecjac_HI     - local gas energy Jacobian w.r.t. nHI
  !     ecjac_HeI    - local gas energy Jacobian w.r.t. nHeI
  !     ecjac_HeII   - local gas energy Jacobian w.r.t. nHeII
  !     HIjac_Er     - local HI chemistry Jacobian w.r.t. Er
  !     HIjac_ec     - local HI chemistry Jacobian w.r.t. ec
  !     HIjac_HI     - local HI chemistry Jacobian w.r.t. nHI
  !     HIjac_HeI    - local HI chemistry Jacobian w.r.t. nHeI
  !     HIjac_HeII   - local HI chemistry Jacobian w.r.t. nHeII
  !     HeIjac_Er    - local HeI chemistry Jacobian w.r.t. Er
  !     HeIjac_ec    - local HeI chemistry Jacobian w.r.t. ec
  !     HeIjac_HI    - local HeI chemistry Jacobian w.r.t. nHI
  !     HeIjac_HeI   - local HeI chemistry Jacobian w.r.t. nHeI
  !     HeIjac_HeII  - local HeI chemistry Jacobian w.r.t. nHeII
  !     HeIIjac_Er   - local HeII chemistry Jacobian w.r.t. Er
  !     HeIIjac_ec   - local HeII chemistry Jacobian w.r.t. ec
  !     HeIIjac_HI   - local HeII chemistry Jacobian w.r.t. nHI
  !     HeIIjac_HeI  - local HeII chemistry Jacobian w.r.t. nHeI
  !     HeIIjac_HeII - local HeII chemistry Jacobian w.r.t. nHeII
  !     ier          - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in) :: Nchem, Model, ESpectrum, DualEnergy
  INTG_PREC, intent(in) :: NTempBins
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in) :: a, adot
  R_PREC, intent(in) :: time, HFrac, gamma, IsE, IsEsHI, IsEsHInu
  R_PREC, intent(in) :: IsEsHeI, IsEsHeInu, IsEsHeII, IsEsHeIInu
  R_PREC, intent(in) :: TempStart, TempEnd, piHI, piHeI, piHeII
  R_PREC, intent(in) :: CompA, Comp_xray, Comp_temp
  R_PREC, intent(in) :: aUnits, DenUnits, VelUnits, ErUnits, ecUnits, NiUnits
  R_PREC, intent(in), dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       Era, eca, n_HIa, n_HeIa, n_HeIIa, eha, rhoa, vx, vy, vz
  R_PREC, intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb, k4Tb, k5Tb,   &
       k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb,          &
       ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(out) :: &
       Erjac_Er, Erjac_ec, Erjac_HI, Erjac_HeI, Erjac_HeII,    &
       ecjac_Er, ecjac_ec, ecjac_HI, ecjac_HeI, ecjac_HeII,       &
       HIjac_Er, HIjac_ec, HIjac_HI, HIjac_HeI, HIjac_HeII,       &
       HeIjac_Er, HeIjac_ec, HeIjac_HI, HeIjac_HeI, HeIjac_HeII,  &
       HeIIjac_Er, HeIIjac_ec, HeIIjac_HI, HeIIjac_HeI, HeIIjac_HeII

  !--------------
  ! locals
  INTG_PREC :: i, j, k, Tidx, Tidxp
  R_PREC :: lTempS, lTempE, dlTemp, lTemp, Tl, Tr, Tfac, Tfac2
  R_PREC :: k1, k2, k3, k4, k5, k6
  R_PREC :: dk1_dT, dk2_dT, dk3_dT, dk4_dT, dk5_dT, dk6_dT
  R_PREC :: aval, grey
  R_PREC :: T, dT_dec, dT_dHI, dT_dHeI, dT_dHeII, min_temp
  R_PREC :: pi, c, hp, kb, StBz, eint, Er, rho, zr, afac, mp, KEconst
  R_PREC :: mol_weight, dmu_dHI, dmu_dHeI, dmu_dHeII
  R_PREC :: alpha, beta, comp1, comp2
  R_PREC :: HIconst, HeIconst, HeIIconst, nu0_HI, nu0_HeI, nu0_HeII
  R_PREC :: nH, nHI, dnHI_dHI, nHII, dnHII_dHI
  R_PREC :: ne, dne_dHI, dne_dHeI, dne_dHeII
  R_PREC :: nHe, nHeI, dnHeI_dHeI, nHeII, dnHeII_dHeII
  R_PREC :: nHeIII, dnHeIII_dHeI, dnHeIII_dHeII 
  R_PREC :: kappaE, dkappaE_dHI, dkappaE_dHeI, dkappaE_dHeII
  R_PREC :: eta, deta_dec, deta_dHI, deta_dHeI, deta_dHeII
  R_PREC :: brem, dbrem_dT, ceHI, dceHI_dT, ceHeI, dceHeI_dT
  R_PREC :: ceHeII, dceHeII_dT, ciHI, dciHI_dT, ciHeI, dciHeI_dT
  R_PREC :: ciHeIS, dciHeIS_dT, ciHeII, dciHeII_dT, reHII, dreHII_dT
  R_PREC :: reHeII1, dreHeII1_dT, reHeII2, dreHeII2_dT, reHeIII, dreHeIII_dT
  R_PREC :: dG_dEr, dG_dHI, dG_dHeI, dG_dHeII
  R_PREC :: dLambda_dec, dLambda_dHI, dLambda_dHeI, dLambda_dHeII
  R_PREC :: G_HI, dGHI_dEr, G_HeI, dGHeI_dEr, G_HeII, dGHeII_dEr

  !=======================================================================

  ! set flag to success
  ier = 1

  ! initialize outputs to have all zero values
  Erjac_Er = 0._RKIND
  Erjac_ec = 0._RKIND
  ecjac_Er = 0._RKIND
  if (Nchem == 1) then
     Erjac_HI = 0._RKIND

     ecjac_HI = 0._RKIND

     HIjac_Er = 0._RKIND
     HIjac_ec = 0._RKIND
     HIjac_HI = 0._RKIND
  endif
  if (Nchem == 3) then
     Erjac_HI   = 0._RKIND
     Erjac_HeI  = 0._RKIND
     Erjac_HeII = 0._RKIND

     ecjac_HI   = 0._RKIND
     ecjac_HeI  = 0._RKIND
     ecjac_HeII = 0._RKIND

     HIjac_Er   = 0._RKIND
     HIjac_ec   = 0._RKIND
     HIjac_HI   = 0._RKIND
     HIjac_HeI  = 0._RKIND
     HIjac_HeII = 0._RKIND

     HeIjac_Er   = 0._RKIND
     HeIjac_ec   = 0._RKIND
     HeIjac_HI   = 0._RKIND
     HeIjac_HeI  = 0._RKIND
     HeIjac_HeII = 0._RKIND

     HeIIjac_Er   = 0._RKIND
     HeIIjac_ec   = 0._RKIND
     HeIIjac_HI   = 0._RKIND
     HeIIjac_HeI  = 0._RKIND
     HeIIjac_HeII = 0._RKIND
  endif


  ! initialize constants
  aval = a*aunits
  pi = pi_val
  afac = adot/a                ! adot/a
  mp = mass_h                  ! mass of a proton [g]
  c  = c_light                 ! speed of light    [cm/s]
  hp = hplanck                 ! Planck's constant [ergs*s]
  kb = kboltz                  ! boltzmann constant [erg/K]
  StBz  = 5.6704e-5_RKIND      ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  min_temp = 0.1_RKIND         ! minimum temperature [K]
  alpha = -0.4910454_RKIND     ! exponent in emissivity fitting
  beta  = 2.17748887e-24_RKIND ! scaling in emissivity fitting
  nu0_HI = 13.6*ev2erg/hp      ! ionization frequency of HI
  nu0_HeI = 24.6*ev2erg/hp     ! ionization frequency of HeI
  nu0_HeII = 54.4*ev2erg/hp    ! ionization frequency of HeII
  zr = 1._RKIND/aval - 1._RKIND        ! cosmological redshift
  if (aval .ne. 1._RKIND) then        ! Compton cooling coefficients
     Comp1 = CompA*(1._RKIND + zr)**4
     Comp2 = 2.73_RKIND*(1._RKIND + zr)
  else
     Comp1 = 0._RKIND
     Comp2 = 0._RKIND
  endif
  if (DualEnergy == 1) then
     KEconst = 0._RKIND
  else
     KEconst = 0.5_RKIND
  endif
  grey = 1._RKIND
  if (ESpectrum == -1)  grey = 0._RKIND

  !   lookup table constants
  lTempS = log(TempStart)
  lTempE = log(TempEnd)
  dlTemp = (lTempE - lTempS)/(1._RKIND*NTempBins - 1._RKIND)

  ! set shortcut values
  HIconst   = c*(IsEsHI - nu0_HI*IsEsHInu)/IsE     ! matches RadHydro Document
  HeIconst  = c*(IsEsHeI - nu0_HeI*IsEsHeInu)/IsE
  HeIIconst = c*(IsEsHeII - nu0_HeII*IsEsHeIInu)/IsE

  if (Model /= 2) then
     write(0,*) 'gFLDProblem_LocalJac2: incorrect Model =',Model
  endif
  
  
  ! Hydrogen only case
  if (Nchem == 1) then
     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1

              ! set shortcut values for this spatial location
              ! converting densities from comoving to proper
              Er  = Era(i,j,k)*ErUnits
              rho = rhoa(i,j,k)*DenUnits
              nHI = n_HIa(i,j,k)*NiUnits
              nH = Hfrac*rho/mp
              nHII = max(nH - nHI, 0._RKIND)
              ne  = nHII
              dnHI_dHI = NiUnits
              dnHII_dHI = -dnHI_dHI
              dne_dHI = dnHII_dHI

              ! compute temperature and shortcuts
              !!!! This assumes ALL density is H, otherwise mol_weight it wrong !!!!
              eint = eca(i,j,k)*ecUnits + VelUnits*VelUnits*(eha(i,j,k)  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              mol_weight = rho/mp/(nHI + nHII + ne)
              dmu_dHI = -rho/mp*(dnHI_dHI + dnHII_dHI + dne_dHI) &
                   / (nHI + nHII + ne)**2
              T = (gamma-1._RKIND)*mol_weight*mp*eint/kb
              T = max(T,min_temp)
              dT_dec = ecUnits*(gamma-1._RKIND)*mol_weight*mp/kb
              dT_dHI = (gamma-1._RKIND)*mp*eint/kb*dmu_dHI

              ! look up rates, derivatives
              lTemp = min(max(log(T), lTempS), lTempE)
              Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
              Tidxp = Tidx+1
              Tl = lTempS + (Tidx-1)*dlTemp
              Tr = lTempS +  Tidx*dlTemp
              Tfac = (lTemp - Tl)/(Tr - Tl)
              Tfac2 = 1._RKIND/(Tr - Tl)/T
              k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
              dk1_dT = (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac2
              ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
              dceHI_dT = (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac2
              ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
              dciHI_dT = (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac2
              reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
              dreHII_dT = (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac2
              brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac
              dbrem_dT = (bremTb(Tidxp) - bremTb(Tidx))*Tfac2

              ! compute opacities and derivatives
              kappaE = nHI*IsEsHI/IsE
              dkappaE_dHI = IsEsHI/IsE*dnHI_dHI

              ! compute Hydrogen photoionization rate & derivs
              G_HI = c*Er/hp*IsEsHInu/IsE
              dGHI_dEr = c/hp*IsEsHInu/IsE*ErUnits

              ! compute fluid cooling rate.  Terms (in order):
              !    Collisional Excitations
              !    Collisional Ionizations
              !    Recombinations
              !    Compton cooling or heating 
              !    X-ray Compton heating
              !    Bremsstrahlung
              !  Lambda = ne/rho*(ceHI*nHI + ciHI*nHI
              !      + reHII*nHII + Comp1*(T-Comp2)
              !      + Comp_xray*(T-Comp_temp) + brem*nHII)
              dLambda_dec = ne/rho*dT_dec*(dceHI_dT*nHI       &
                   + dciHI_dT*nHI + dreHII_dT*nHII            &
                   + Comp1 + Comp_xray + dbrem_dT*nHII)
              dLambda_dHI = dne_dHI/rho*(ceHI*nHI + ciHI*nHI      &
                   + reHII*nHII + Comp1*(T-Comp2)                 &
                   + Comp_xray*(T-Comp_temp) + brem*nHII)         &
                   + ne/rho*(dceHI_dT*dT_dHI*nHI + ceHI*dnHI_dHI  &
                   + dciHI_dT*dT_dHI*nHI + ciHI*dnHI_dHI          &
                   + dreHII_dT*dT_dHI*nHII + reHII*dnHII_dHI      &
                   + Comp1*dT_dHI + Comp_xray*dT_dHI              &
                   + dbrem_dT*dT_dHI*nHII + brem*dnHII_dHI)

              ! compute fluid heating rate
              !   G = nHI/rho*Er*HIconst
              dG_dEr = nHI/rho*ErUnits*HIconst
              dG_dHI = dnHI_dHI/rho*Er*HIconst

              ! compute emissivity 
              eta = nHII*ne*beta*T**alpha
              deta_dec = dT_dec*alpha*nHII*ne*beta*T**(alpha-1._RKIND)
              deta_dHI = dT_dHI*alpha*nHII*ne*beta*T**(alpha-1._RKIND) &
                   + (dnHII_dHI*ne + nHII*dne_dHI)*beta*T**alpha

              ! put it all together
              !   rhs_Er = (src_Er + 4._RKIND*pi*eta 
              !            - (afac*grey + c*kappaE)*Er)/ErUnits
              Erjac_Er(i,j,k) = -(afac*grey + c*kappaE)
              Erjac_ec(i,j,k) =  (4._RKIND*pi*deta_dec)/ErUnits
              Erjac_HI(i,j,k) =  (4._RKIND*pi*deta_dHI - c*Er*dkappaE_dHI)/ErUnits

              !   rhs_HI = (src_HI + k2*ne*nHII - nHI*(k1*ne + G_HI))/NiUnits
              HIjac_Er(i,j,k) = -nHI*dGHI_dEr/NiUnits
              HIjac_ec(i,j,k) = ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dec/NiUnits
              HIjac_HI(i,j,k) = (dk2_dT*dT_dHI*ne*nHII + k2*dne_dHI*nHII &
                   + k2*ne*dnHII_dHI - dnHI_dHI*k1*ne - nHI*k1*dne_dHI   &
                   - nHI*dk1_dT*dT_dHI*ne - dnHI_dHI*G_HI)/NiUnits

              !   rhs_ec = (src_ec + G - Lambda - 2._RKIND*afac*ec)/ecUnits
              ecjac_ec(i,j,k) = (ecjac_ec(i,j,k) - 2._RKIND*afac*ecUnits &
                   - dLambda_dec)/ecUnits
              ecjac_Er(i,j,k) = dG_dEr/ecUnits
              ecjac_HI(i,j,k) = (dG_dHI - dLambda_dHI)/ecUnits

           enddo
        enddo
     enddo

     ! Hydrogen + Helium case
  else if (Nchem == 3) then

     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1

              ! set shortcut values for this spatial location
              ! converting densities from comoving to proper
              Er  = Era(i,j,k)*ErUnits
              rho = rhoa(i,j,k)*DenUnits
              nHI = n_HIa(i,j,k)*NiUnits
              nH = Hfrac*rho/mp
              nHII = max(nH - nHI, 0._RKIND)
              nHe    = (1._RKIND-HFrac)*rho/4._RKIND/mp
              nHeI   = n_HeIa(i,j,k)*NiUnits/4._RKIND
              nHeII  = n_HeIIa(i,j,k)*NiUnits/4._RKIND
              nHeIII = max(nHe - nHeI - nHeII, 0._RKIND)
              ne     = nHII + 0.25_RKIND*nHeII + 0.5_RKIND*nHeIII
              dnHI_dHI = NiUnits
              dnHII_dHI = -dnHI_dHI
              dnHeI_dHeI = NiUnits/4._RKIND
              dnHeII_dHeII = NiUnits/4._RKIND
              dnHeIII_dHeI = -dnHeI_dHeI
              dnHeIII_dHeII = -dnHeII_dHeII
              dne_dHI = dnHII_dHI
              dne_dHeI = 0.5_RKIND*dnHeIII_dHeI
              dne_dHeII = 0.25_RKIND*dnHeII_dHeII + 0.5_RKIND*dnHeIII_dHeII

              ! compute temperature and shortcuts
              eint = eca(i,j,k)*ecUnits + VelUnits*VelUnits*(eha(i,j,k)  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              mol_weight = rho/mp/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)
              dmu_dHI = -rho/mp/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)**2 &
                   * (dnHI_dHI+dnHII_dHI+dne_dHI)
              dmu_dHeI = -rho/mp/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)**2 &
                   * (0.25_RKIND*(dnHeI_dHeI+dnHeIII_dHeI) + dne_dHeI)
              dmu_dHeII = -rho/mp/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)**2 &
                   * (0.25_RKIND*(dnHeII_dHeII+dnHeIII_dHeII) + dne_dHeII)
              T = (gamma-1._RKIND)*mol_weight*mp*eint/kb
              T = max(T,min_temp)
              dT_dec = ecUnits*(gamma-1._RKIND)*mol_weight*mp/kb
              dT_dHI = (gamma-1._RKIND)*mp*eint/kb*dmu_dHI
              dT_dHeI = (gamma-1._RKIND)*mp*eint/kb*dmu_dHeI
              dT_dHeII = (gamma-1._RKIND)*mp*eint/kb*dmu_dHeII

              ! look up rates
              lTemp = min(max(log(T), lTempS), lTempE)
              Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
              Tidxp = Tidx+1
              Tl = lTempS + (Tidx-1)*dlTemp
              Tr = lTempS +  Tidx*dlTemp
              Tfac = (lTemp - Tl)/(Tr - Tl)
              Tfac2 = 1._RKIND/(Tr - Tl)/T
              k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
              dk1_dT = (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac2
              k2 = k2Tb(Tidx) + (k2Tb(Tidxp) - k2Tb(Tidx))*Tfac
              dk2_dT = (k2Tb(Tidxp) - k2Tb(Tidx))*Tfac2
              k3 = k3Tb(Tidx) + (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac
              dk3_dT = (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac2
              k4 = k4Tb(Tidx) + (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac
              dk4_dT = (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac2
              k5 = k5Tb(Tidx) + (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac
              dk5_dT = (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac2
              k6 = k6Tb(Tidx) + (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac
              dk6_dT = (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac2
              ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
              dceHI_dT = (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac2
              ceHeI = ceHeITb(Tidx) + (ceHeITb(Tidxp) - ceHeITb(Tidx))*Tfac
              dceHeI_dT = (ceHeITb(Tidxp) - ceHeITb(Tidx))*Tfac2
              ceHeII = ceHeIITb(Tidx) + (ceHeIITb(Tidxp) - ceHeIITb(Tidx))*Tfac
              dceHeII_dT = (ceHeIITb(Tidxp) - ceHeIITb(Tidx))*Tfac2
              ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
              dciHI_dT = (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac2
              ciHeI = ciHeITb(Tidx) + (ciHeITb(Tidxp) - ciHeITb(Tidx))*Tfac
              dciHeI_dT = (ciHeITb(Tidxp) - ciHeITb(Tidx))*Tfac2
              ciHeII = ciHeIITb(Tidx) + (ciHeIITb(Tidxp) - ciHeIITb(Tidx))*Tfac
              dciHeII_dT = (ciHeIITb(Tidxp) - ciHeIITb(Tidx))*Tfac2
              ciHeIS = ciHeISTb(Tidx) + (ciHeISTb(Tidxp) - ciHeISTb(Tidx))*Tfac
              dciHeIS_dT = (ciHeISTb(Tidxp) - ciHeISTb(Tidx))*Tfac2
              reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
              dreHII_dT = (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac2
              reHeII1 = reHeII1Tb(Tidx) + (reHeII1Tb(Tidxp) - reHeII1Tb(Tidx))*Tfac
              dreHeII1_dT = (reHeII1Tb(Tidxp) - reHeII1Tb(Tidx))*Tfac2
              reHeII2 = reHeII2Tb(Tidx) + (reHeII2Tb(Tidxp) - reHeII2Tb(Tidx))*Tfac
              dreHeII2_dT = (reHeII2Tb(Tidxp) - reHeII2Tb(Tidx))*Tfac2
              reHeIII = reHeIIITb(Tidx) + (reHeIIITb(Tidxp) - reHeIIITb(Tidx))*Tfac
              dreHeIII_dT = (reHeIIITb(Tidxp) - reHeIIITb(Tidx))*Tfac2
              brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac
              dbrem_dT = (bremTb(Tidxp) - bremTb(Tidx))*Tfac2

              ! compute opacities and derivatives
              kappaE = nHI*IsEsHI/IsE+nHeI*IsEsHeI/IsE+nHeII*IsEsHeII/IsE
              dkappaE_dHI = IsEsHI/IsE*dnHI_dHI
              dkappaE_dHeI = IsEsHeI/IsE*dnHeI_dHeI
              dkappaE_dHeII = IsEsHeII/IsE*dnHeII_dHeII

              ! compute Hydrogen photoionization rate & derivs
              G_HI = c*Er/hp*IsEsHInu/IsE
              dGHI_dEr = c/hp*IsEsHInu/IsE*ErUnits

              ! compute Helium photoionization rates & derivs
              G_HeI  = c*Er/hp*IsEsHeInu/IsE
              G_HeII = c*Er/hp*IsEsHeIInu/IsE
              dGHeI_dEr = c/hp*IsEsHeInu/IsE*ErUnits
              dGHeII_dEr = c/hp*IsEsHeIInu/IsE*ErUnits

              ! compute fluid cooling rate.  Terms (in order):
              !    Collisional Excitations
              !    Collisional Ionizations
              !    Recombinations
              !    Compton cooling or heating 
              !    X-ray Compton heating
              !    Bremsstrahlung
              !  Lambda = ne/rho*(ceHI*nHI + ciHI*nHI + reHII*nHII
              !     + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII
              !     + reHeII1 + reHeII2 + brem)*nHeII/4._RKIND
              !     + (ciHeI*nHeI + reHeIII*nHeIII)/4._RKIND
              !     + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp)
              !     + brem*(nHII+nHeIII))
              dLambda_dec = ne/rho*dT_dec*(dceHI_dT*nHI + dciHI_dT*nHI &
                   + (dceHeI_dT*ne + dciHeIS_dT*ne + dceHeII_dT        &
                   + dciHeII_dT + dreHeII1_dT + dreHeII2_dT            &
                   + dbrem_dT)*nHeII/4._RKIND + dciHeI_dT*nHeI/4._RKIND        &
                   + dreHeIII_dT*nHeIII/4._RKIND + dreHII_dT*nHII          &
                   + Comp1 + Comp_xray + dbrem_dT*(nHII+nHeIII))
              dLambda_dHI = ne/rho*(                                 &
                   (ceHI + ciHI)*dnHI_dHI + (reHII+brem)*dnHII_dHI + &
                   dT_dHI*(dceHI_dT*nHI + dciHI_dT*nHI               &
                   + (dceHeI_dT*ne + dciHeIS_dT*ne + dceHeII_dT      &
                   + dciHeII_dT + dreHeII1_dT + dreHeII2_dT          &
                   + dbrem_dT)*nHeII/4._RKIND + dciHeI_dT*nHeI/4._RKIND      &
                   + dreHeIII_dT*nHeIII/4._RKIND + dreHII_dT*nHII        &
                   + Comp1 + Comp_xray + dbrem_dT*(nHII+nHeIII))) +  &
                   dne_dHI/rho*(ceHI*nHI + reHII*nHII + ciHI*nHI     &
                   + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII         &
                   + reHeII1 + reHeII2 + brem)*nHeII/4._RKIND            &
                   + (reHeIII*nHeIII + ciHeI*nHeI)/4._RKIND              &
                   + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp)       &
                   + brem*(nHII+nHeIII))
              dLambda_dHeI = ne/rho*(                                &
                   + (ciHeI*dnHeI_dHeI + reHeIII*dnHeIII_dHeI)/4._RKIND  &
                   + brem*dnHeIII_dHeI +                             &
                   dT_dHeI*(dceHI_dT*nHI + dciHI_dT*nHI              &
                   + (dceHeI_dT*ne + dciHeIS_dT*ne + dceHeII_dT      &
                   + dciHeII_dT + dreHeII1_dT + dreHeII2_dT          &
                   + dbrem_dT)*nHeII/4._RKIND + dciHeI_dT*nHeI/4._RKIND      &
                   + dreHeIII_dT*nHeIII/4._RKIND + dreHII_dT*nHII        &
                   + Comp1 + Comp_xray + dbrem_dT*(nHII+nHeIII))) +  &
                   dne_dHeI/rho*(ceHI*nHI + reHII*nHII + ciHI*nHI    &
                   + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII         &
                   + reHeII1 + reHeII2 + brem)*nHeII/4._RKIND            &
                   + (reHeIII*nHeIII + ciHeI*nHeI)/4._RKIND              &
                   + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp)       &
                   + brem*(nHII+nHeIII))
              dLambda_dHeII = ne/rho*(                               &
                   + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII         &
                   + reHeII1 + reHeII2 + brem)*dnHeII_dHeII/4._RKIND     &
                   + (reHeIII/4._RKIND + brem)*dnHeIII_dHeII +           &
                   dT_dHeII*(dceHI_dT*nHI + dciHI_dT*nHI             &
                   + (dceHeI_dT*ne + dciHeIS_dT*ne + dceHeII_dT      &
                   + dciHeII_dT + dreHeII1_dT + dreHeII2_dT          &
                   + dbrem_dT)*nHeII/4._RKIND + dciHeI_dT*nHeI/4._RKIND      &
                   + dreHeIII_dT*nHeIII/4._RKIND + dreHII_dT*nHII        &
                   + Comp1 + Comp_xray + dbrem_dT*(nHII+nHeIII))) +  &
                   dne_dHeII/rho*(ceHI*nHI + reHII*nHII + ciHI*nHI   &
                   + (ceHeI*ne + ciHeIS*ne + ceHeII + ciHeII         &
                   + reHeII1 + reHeII2 + brem)*nHeII/4._RKIND            &
                   + (reHeIII*nHeIII + ciHeI*nHeI)/4._RKIND              &
                   + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp)       &
                   + brem*(nHII+nHeIII))

              ! compute fluid heating rate
              !   G = Er/rho*(nHI*HIconst + nHeI*HeIconst + nHeII*HeIIconst)
              dG_dEr = ErUnits/rho*(nHI*HIconst + nHeI*HeIconst + nHeII*HeIIconst)
              dG_dHI = Er/rho*HIconst*dnHI_dHI
              dG_dHeI = Er/rho*HeIconst*dnHeI_dHeI
              dG_dHeII = Er/rho*HeIIconst*dnHeII_dHeII

              ! compute emissivity 
              ! [WE HAVE NO FORMULA FOR ETA IN MULTISPECIES CASE!!]
              eta        = 0._RKIND
              deta_dec   = 0._RKIND
              deta_dHI   = 0._RKIND
              deta_dHeI  = 0._RKIND
              deta_dHeII = 0._RKIND

              ! put it all together
              !   rhs_Er = (src_Er + 4._RKIND*pi*eta - (afac*grey + c*kappaE)*Er)/ErUnits
              Erjac_Er(i,j,k)   = -(afac*grey + c*kappaE)
              Erjac_ec(i,j,k)   =  (4._RKIND*pi*deta_dec)/ErUnits
              Erjac_HI(i,j,k)   =  (4._RKIND*pi*deta_dHI-c*Er*dkappaE_dHI)/ErUnits
              Erjac_HeI(i,j,k)  =  (4._RKIND*pi*deta_dHeI-c*Er*dkappaE_dHeI)/ErUnits
              Erjac_HeII(i,j,k) =  (4._RKIND*pi*deta_dHeII-c*Er*dkappaE_dHeII)/ErUnits

              !   rhs_ec = (src_ec + G - Lambda - 2._RKIND*afac*ec)/ecUnits
              ecjac_ec(i,j,k) = (ecjac_ec(i,j,k) - 2._RKIND*afac*ecUnits &
                   - dLambda_dec)/ecUnits
              ecjac_Er(i,j,k) = dG_dEr/ecUnits
              ecjac_HI(i,j,k) = (dG_dHI - dLambda_dHI)/ecUnits
              ecjac_HeI(i,j,k) = (dG_dHeI - dLambda_dHeI)/ecUnits
              ecjac_HeII(i,j,k) = (dG_dHeII - dLambda_dHeII)/ecUnits

              !   rhs_HI = (src_HI + ne*(k2*nHII - k1*nHI) - nHI*G_HI)/NiUnits
              HIjac_Er(i,j,k) = -nHI*dGHI_dEr/NiUnits
              HIjac_ec(i,j,k) = ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dec/NiUnits
              HIjac_HI(i,j,k) = (ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dHI &
                   + (k2*nHII - k1*nHI)*dne_dHI + k2*ne*dnHII_dHI     &
                   - dnHI_dHI*k1*ne - dnHI_dHI*G_HI)/NiUnits
              HIjac_HeI = (ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dHeI      &
                   + (k2*nHII - k1*nHI)*dne_dHeI)/NiUnits
              HIjac_HeII = (ne*(nHII*dk2_dT - nHI*dk1_dT)*dT_dHeII    &
                   + (k2*nHII - k1*nHI)*dne_dHeII)/NiUnits

              !  rhs_HeI = (src_HeI - nHeI*G_HeI 
              !            + ne*(k4*nHeII - k3*nHeI))/NiUnits
              HeIjac_Er(i,j,k) = -nHeI*dGHeI_dEr/NiUnits
              HeIjac_ec(i,j,k) = ne*(nHeII*dk4_dT - nHeI*dk3_dT)*dT_dec/NiUnits
              HeIjac_HI(i,j,k) = (ne*(nHeII*dk4_dT - nHeI*dk3_dT)*dT_dHI   &
                   + (k4*nHeII - k3*nHeI)*dne_dHI)/NiUnits
              HeIjac_HeI(i,j,k) = (ne*(nHeII*dk4_dT - nHeI*dk3_dT)*dT_dHeI &
                   + (k4*nHeII - k3*nHeI)*dne_dHeI - ne*k3*dnHeI_dHeI      &
                   - G_HeI*dnHeI_dHeI)/NiUnits
              HeIjac_HeII(i,j,k) = (ne*(nHeII*dk4_dT-nHeI*dk3_dT)*dT_dHeII &
                   + (k4*nHeII - k3*nHeI)*dne_dHeII                        &
                   + ne*k4*dnHeII_dHeII)/NiUnits

              ! rhs_HeII = (src_HeII + nHeI*G_HeI - nHeII*G_HeII
              !            + ne*(k3*nHeI + k6*nHeIII 
              !            - k4*nHeII + k5*nHeII))/NiUnits
              HeIIjac_Er(i,j,k) = -nHeII*dGHeII_dEr/NiUnits
              HeIIjac_ec(i,j,k) = ne*(nHeI*dk3_dT + nHeIII*dk6_dT    &
                   - nHeII*dk4_dT + nHeII*dk5_dT)*dT_dec/NiUnits
              HeIIjac_HI(i,j,k) = (ne*(nHeI*dk3_dT + nHeIII*dk6_dT   &
                   - nHeII*dk4_dT + nHeII*dk5_dT)*dT_dHI             &
                   + (k3*nHeI + k6*nHeIII - k4*nHeII                 &
                   + k5*nHeII)*dne_dHI)/NiUnits
              HeIIjac_HeI(i,j,k) = (ne*(nHeI*dk3_dT + nHeIII*dk6_dT  &
                   - nHeII*dk4_dT + nHeII*dk5_dT)*dT_dHeI            &
                   + (k3*nHeI + k6*nHeIII - k4*nHeII                 &
                   + ne*(k3*dnHeI_dHeI + k6*dnHeIII_dHeI)            &
                   + k5*nHeII)*dne_dHeI)/NiUnits
              HeIIjac_HeII(i,j,k) = (ne*(nHeI*dk3_dT + nHeIII*dk6_dT &
                   - nHeII*dk4_dT + nHeII*dk5_dT)*dT_dHeII           &
                   + (k3*nHeI + k6*nHeIII - k4*nHeII                 &
                   + k5*nHeII)*dne_dHeII + ne*(k6*dnHeIII_dHeII      &
                   - (k4-k5)*dnHeII_dHeII) - G_HeII*dnHeII_dHeII)/NiUnits

           enddo
        enddo
     enddo

  else
     write(0,*) 'gFLDProblem_LocalJac2: illegal Nchem =',Nchem
  endif  ! Nchem

  return
end subroutine gFLDProblem_LocalJac2
!=======================================================================





subroutine gFLDProblem_LocalJac4(Erjac_Er, Erjac_ec, Erjac_HI, ecjac_Er,  &
     ecjac_ec, ecjac_HI, HIjac_Er, HIjac_ec, HIjac_HI, time, Era, eca,    &
     n_HIa, eha, rhoa, vx, vy, vz, Nchem, HFrac, Model, ESpectrum,        &
     DualEnergy, a, adot, IsE, IsEsHI, IsEsHInu, gamma, NTempBins,        &
     TempStart, TempEnd, k1Tb, k2Tb, DenUnits, VelUnits, ErUnits,         &
     ecUnits, NiUnits, ecScale, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, &
     NGzr, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August, 2006
  !
  !  PURPOSE: Computes the spatially-local components of the chemical 
  !           Jacobian for the Gray FLD problem (Model 4).
  !
  !  INPUTS:
  !     time       - simulation time for evaluation
  !     Era        - Gray radiation energy density
  !     eca        - fluid energy correction array
  !     n_HIa      - number density of Hydrogen I species
  !     n_HeIa     - number density of Helium I species
  !     n_HeIIa    - number density of Helium II species
  !     eha        - total fluid energy array
  !     rhoa       - fluid density, assumed partitioned into either 
  !     vx         - fluid velocity (x-direction)
  !     vy         - fluid velocity (y-direction)
  !     vz         - fluid velocity (z-direction)
  !     Nchem      - number of chemical species (allowed: 0, 1, 3)
  !     HFrac      - percentage of mass composed of Hydrogen
  !     Model      - flag denoting physical model to use
  !     ESpectrum  - radiation spectrum choice
  !                       1 -> 1e5 black body spectrum
  !                       0 -> power law spectrum
  !                      -1 -> monochromatic 
  !     DualEnergy - flag denoting dual energy formalism
  !     a          - cosmological expansion parameter
  !     adot       - da/dt
  !     IsE        - int_{nu0_HI}^{inf} sigE dnu
  !     IsEsHI     - int_{nu0_HI}^{inf} sigE*sigHI dnu
  !     IsEsHInu   - int_{nu0_HI}^{inf} sigE*sigHI/nu dnu
  !     gamma      - ideal gas parameter (ratio of specific heats)
  !     *Units     - variable scaling constants
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: 
  !     Erjac_Er     - local radiation energy Jacobian w.r.t. Er
  !     Erjac_ec     - local radiation energy Jacobian w.r.t. ec
  !     Erjac_HI     - local radiation energy Jacobian w.r.t. nHI
  !     ecjac_Er     - local gas energy Jacobian w.r.t. Er
  !     ecjac_ec     - local gas energy Jacobian w.r.t. ec
  !     ecjac_HI     - local gas energy Jacobian w.r.t. nHI
  !     HIjac_Er     - local HI chemistry Jacobian w.r.t. Er
  !     HIjac_ec     - local HI chemistry Jacobian w.r.t. ec
  !     HIjac_HI     - local HI chemistry Jacobian w.r.t. nHI
  !     ier          - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: Nchem, Model, ESpectrum, DualEnergy
  INTG_PREC, intent(in)  :: NTempBins
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a, adot
  R_PREC,    intent(in) :: time, HFrac, gamma, IsE, IsEsHI, IsEsHInu
  R_PREC,    intent(in) :: TempStart, TempEnd
  R_PREC,    intent(in) :: DenUnits, VelUnits, ErUnits, ecUnits, NiUnits, ecScale
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) :: &
       Era, eca, n_HIa, eha, rhoa, vx, vy, vz
  R_PREC, intent(in), dimension(NTempBins) :: k1Tb, k2Tb
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(out) :: &
       Erjac_Er, Erjac_ec, Erjac_HI, ecjac_Er, ecjac_ec, ecjac_HI, &
       HIjac_Er, HIjac_ec, HIjac_HI

  !--------------
  ! locals
  INTG_PREC :: i, j, k, Tidx, Tidxp
  R_PREC :: lTempS, lTempE, dlTemp, lTemp, Tl, Tr, Tfac, k1, k2
  R_PREC :: afac, grey, T, min_temp, lamT
  R_PREC :: c, hp, kb, StBz, eint, Er, rho, mp, KEconst
  R_PREC :: mol_weight, nH, nHI, nHII, ne
  R_PREC :: kappaE, dkappaE_dHI, G_HI, dGHI_dEr
  R_PREC :: dnHI_dHI, dnHII_dHI, dne_dHI

  !=======================================================================

  ! set flag to success
  ier = 1

  if (Model /= 4) then
     write(0,*) 'gFLDProblem_LocalJac4: illegal Model =',Model
  endif

  ! ensure that we're running with proper chemistry
  if (Nchem /= 1) then
     write(0,*) 'gFLDProblem_LocalJac4: Nchem = ',Nchem,' and Model =', &
          Model,' do not match!'
     return
  endif

  ! initialize outputs to have all zero values
  Erjac_Er = 0._RKIND
  Erjac_ec = 0._RKIND
  Erjac_HI = 0._RKIND

  ecjac_Er = 0._RKIND
  ecjac_HI = 0._RKIND

  HIjac_Er = 0._RKIND
  HIjac_ec = 0._RKIND
  HIjac_HI = 0._RKIND


  ! initialize constants
  afac = adot/a        ! adot/a
  mp = mass_h          ! mass of a proton [g]
  c  = c_light         ! speed of light [cm/s]
  hp = hplanck         ! Planck's constant [ergs*s]
  kb = kboltz          ! boltzmann constant [erg/K]
  StBz  = 5.6704e-5_RKIND  ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  min_temp = 0.1_RKIND     ! minimum temperature [K]
  mol_weight = 0.6_RKIND   ! mean molecular mass
  if (DualEnergy == 1) then
     KEconst = 0._RKIND
  else
     KEconst = 0.5_RKIND
  endif
  grey = 1._RKIND
  if (ESpectrum == -1)  grey = 0._RKIND

  !   lookup table constants
  lTempS = log(TempStart)
  lTempE = log(TempEnd)
  dlTemp = (lTempE - lTempS)/(1._RKIND*NTempBins - 1._RKIND)

  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           ! set shortcut values for this spatial location
           ! converting densities from comoving to proper
           Er  = Era(i,j,k)*ErUnits
           rho = rhoa(i,j,k)*DenUnits
           nHI = n_HIa(i,j,k)*NiUnits
           nH = Hfrac*rho/mp
           nHII = max(nH - nHI, 0._RKIND)
           ne = nHII
           dnHI_dHI = NiUnits
           dnHII_dHI = -dnHI_dHI
           dne_dHI = dnHII_dHI

           ! compute temperature and shortcuts
           !*** For Model 4 with cosmology, the temperature is held in ecScale ***!
           if (adot == 0._RKIND) then
              eint = eca(i,j,k)*ecUnits + VelUnits*VelUnits*(eha(i,j,k)  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              T = (gamma-1._RKIND)*mol_weight*mp*eint/kb
              T = max(T,min_temp)
           else
              T = ecScale
           endif
           lamT = 3.15614e5_RKIND/T

           ! look up rates
           lTemp = min(max(log(T), lTempS), lTempE)
           Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
           Tidxp = Tidx+1
           Tl = lTempS + (Tidx-1)*dlTemp
           Tr = lTempS +  Tidx*dlTemp
           Tfac = (lTemp - Tl)/(Tr - Tl)
           k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac

           ! compute opacities and derivatives
           kappaE = nHI*IsEsHI/IsE
           dkappaE_dHI = IsEsHI/IsE*NiUnits

           ! compute case B Hydrogen recombination coefficient 
           ! [Hui & Gnedin, 1997: RI^B_{HII}]
           !   (still need this because table holds case A coefficient!)
           k2 = 2.753e-14_RKIND*lamT**(1.5_RKIND) *                 &
                (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-2.242_RKIND)

           ! compute Hydrogen photoionization rate & derivs
           G_HI = c*Er/hp*IsEsHInu/IsE
           dGHI_dEr = c/hp*IsEsHInu/IsE*ErUnits

           ! put it all together
           !   rhs_Er = -(afac*grey + c*kappaE)*Er/ErUnits
           Erjac_Er(i,j,k) = -afac*grey - c*kappaE 
           Erjac_ec(i,j,k) = 0._RKIND
           Erjac_HI(i,j,k) = -c*Er*dkappaE_dHI/ErUnits

           !   rhs_HI = (k2*ne*nHII - nHI*(k1*ne + G_HI))/NiUnits
           HIjac_Er(i,j,k) = -nHI*dGHI_dEr/NiUnits
           HIjac_ec(i,j,k) = 0._RKIND
           HIjac_HI(i,j,k) = (k2*dne_dHI*nHII + k2*ne*dnHII_dHI &
                - dnHI_dHI*(k1*ne + G_HI) - nHI*k1*dne_dHI)/NiUnits

           ecjac_ec(i,j,k) = 0._RKIND
           ecjac_Er(i,j,k) = 0._RKIND
           ecjac_HI(i,j,k) = 0._RKIND

        enddo
     enddo
  enddo

  return
end subroutine gFLDProblem_LocalJac4
!=======================================================================





subroutine gFLDProblem_LocalJac5(Erjac_Er, Erjac_ec, Erjac_HI, ecjac_Er,  &
     ecjac_ec, ecjac_HI, HIjac_Er, HIjac_ec, HIjac_HI, Era, n_HIa, Nchem, &
     Model, ESpectrum, IsE, IsEsHI, a, adot, ErUnits, NiUnits, Nx, Ny,    &
     Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August, 2006
  !
  !  PURPOSE: Computes the spatially-local components of the chemical 
  !           Jacobian for the Gray FLD problem (Model 5).
  !
  !  INPUTS:
  !     Era        - Gray radiation energy density
  !     n_HIa      - number density of Hydrogen I species
  !     Nchem      - number of chemical species (allowed: 0, 1, 3)
  !     Model      - flag denoting physical model to use
  !     ESpectrum  - radiation spectrum choice
  !                       1 -> 1e5 black body spectrum
  !                       0 -> power law spectrum
  !                      -1 -> monochromatic 
  !     IsE        - int_{nu0_HI}^{inf} sigE dnu
  !     IsEsHI     - int_{nu0_HI}^{inf} sigE*sigHI dnu
  !     *Units     - variable scaling constants
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: 
  !     Erjac_Er     - local radiation energy Jacobian w.r.t. Er
  !     Erjac_ec     - local radiation energy Jacobian w.r.t. ec
  !     Erjac_HI     - local radiation energy Jacobian w.r.t. nHI
  !     ecjac_Er     - local gas energy Jacobian w.r.t. Er
  !     ecjac_ec     - local gas energy Jacobian w.r.t. ec
  !     ecjac_HI     - local gas energy Jacobian w.r.t. nHI
  !     HIjac_Er     - local HI chemistry Jacobian w.r.t. Er
  !     HIjac_ec     - local HI chemistry Jacobian w.r.t. ec
  !     HIjac_HI     - local HI chemistry Jacobian w.r.t. nHI
  !     ier          - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in) :: Nchem, Model, ESpectrum
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC,  intent(in) :: a, adot
  R_PREC, intent(in) :: ErUnits, NiUnits
  R_PREC, intent(in) :: IsE, IsEsHI
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) :: &
       Era, n_HIa
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(out) :: &
       Erjac_Er, Erjac_ec, Erjac_HI, ecjac_Er, ecjac_ec, ecjac_HI, &
       HIjac_Er, HIjac_ec, HIjac_HI

  !--------------
  ! locals
  INTG_PREC :: i, j, k
  R_PREC :: c, nHI, Er, kappaE, afac, grey

  !=======================================================================

  ! set flag to success
  ier = 1

  ! ensure that we're running with proper chemistry
  if (Nchem /= 1) then
     write(0,*) 'gFLDProblem_LocalJac5: Nchem = ',Nchem,' and Model =', &
          Model,' do not match!'
     return
  endif

  ! initialize outputs to have all zero values
  Erjac_Er = 0._RKIND
  Erjac_ec = 0._RKIND
  Erjac_HI = 0._RKIND

  ecjac_Er = 0._RKIND
  ecjac_HI = 0._RKIND

  HIjac_Er = 0._RKIND
  HIjac_ec = 0._RKIND
  HIjac_HI = 0._RKIND


  ! initialize constants
  c  = c_light         ! speed of light [cm/s]
  afac = adot/a        ! adot/a
  grey = 1._RKIND
  if (ESpectrum == -1)  grey = 0._RKIND

  if (Model /= 5) then
     write(0,*) 'gFLDProblem_LocalJac5: illegal Model =',Model
  endif
  
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           ! set shortcut values for this spatial location
           ! converting densities from comoving to proper
           Er  = Era(i,j,k)*ErUnits
           nHI = n_HIa(i,j,k)*NiUnits

           ! compute opacities and derivatives
           kappaE = nHI*IsEsHI/IsE

           ! put it all together
           ! rhs_Er = -(afac*grey + c*kappaE)*Er/ErUnits
           Erjac_Er(i,j,k) = -afac*grey - c*kappaE
           Erjac_ec(i,j,k) = 0._RKIND
           Erjac_HI(i,j,k) = 0._RKIND

           HIjac_Er(i,j,k) = 0._RKIND
           HIjac_ec(i,j,k) = 0._RKIND
           HIjac_HI(i,j,k) = 0._RKIND

           ecjac_ec(i,j,k) = 0._RKIND
           ecjac_Er(i,j,k) = 0._RKIND
           ecjac_HI(i,j,k) = 0._RKIND

        enddo
     enddo
  enddo

  return
end subroutine gFLDProblem_LocalJac5
!=======================================================================





subroutine gFLDProblem_LocalJac10(Erjac_Er, Erjac_ec, ecjac_Er, ecjac_ec, &
     time, Era, eca, eha, rhoa, vx, vy, vz, Nchem, Model, ESpectrum,      &
     ProbType, DualEnergy, a, adot, PmC0, PmC1, PmC2, PmC3, PmC4, EmC0,   &
     EmC1, EmC2, EmC3, EmC4, gamma, DenUnits, VelUnits, ErUnits,          &
     ecUnits, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August, 2006
  !  modified1:  September 19, 2007, by John Hayes; minor mods on behalf of
  !              ProbTypes 204 & 205.
  !
  !  PURPOSE: Computes the spatially-local components of the chemical 
  !           Jacobian for the Gray FLD problem (Model 10).
  !
  !  INPUTS:
  !     time       - simulation time for evaluation
  !     Era        - Gray radiation energy density
  !     eca        - fluid energy correction array
  !     eha        - total fluid energy array
  !     rhoa       - fluid density, assumed partitioned into either 
  !     vx         - fluid velocity (x-direction)
  !     vy         - fluid velocity (y-direction)
  !     vz         - fluid velocity (z-direction)
  !     Nchem      - number of chemical species (allowed: 0, 1, 3)
  !     Model      - flag denoting physical model to use
  !     ESpectrum  - radiation spectrum choice
  !                       1 -> 1e5 black body spectrum
  !                       0 -> power law spectrum
  !                      -1 -> monochromatic 
  !     ProbType   - flag denoting problem type (kluge)
  !     DualEnergy - flag denoting dual energy formalism
  !     a          - cosmological expansion parameter
  !     adot       - da/dt
  !     PmC0-PmC4  - input parameters for Planck mean opacity
  !     EmC0-EmC4  - input parameters for Energy mean opacity
  !     gamma      - ideal gas parameter (ratio of specific heats)
  !     *Units     - variable scaling constants
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: 
  !     Erjac_Er     - local radiation energy Jacobian w.r.t. Er
  !     Erjac_ec     - local radiation energy Jacobian w.r.t. ec
  !     ecjac_Er     - local gas energy Jacobian w.r.t. Er
  !     ecjac_ec     - local gas energy Jacobian w.r.t. ec
  !     ier          - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: Nchem, Model, ESpectrum, ProbType, DualEnergy
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a, adot
  R_PREC,    intent(in) :: time, gamma
  R_PREC,    intent(in) :: PmC0, PmC1, PmC2, PmC3, PmC4
  R_PREC,    intent(in) :: EmC0, EmC1, EmC2, EmC3, EmC4
  R_PREC,    intent(in) :: DenUnits, VelUnits, ErUnits, ecUnits
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) :: &
       Era, eca, eha, rhoa, vx, vy, vz
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(out) :: &
       Erjac_Er, Erjac_ec, ecjac_Er, ecjac_ec

  !--------------
  ! locals
  INTG_PREC :: i, j, k
  R_PREC :: T, dT_dec, min_temp, mol_weight
  R_PREC :: pi, c, hp, kb, rc, StBz, eint, Er, rho, afac, grey, mp, KEconst
  R_PREC :: kappaP, dkappaP_dec, kappaE, dkappaE_dec
  R_PREC :: eta, deta_dec, Cv

  !=======================================================================

  ! set flag to success
  ier = 1

  ! initialize outputs to have all zero values
  Erjac_Er = 0._RKIND
  Erjac_ec = 0._RKIND
  ecjac_Er = 0._RKIND


  ! initialize constants
  pi = pi_val
  afac = adot/a                      ! adot/a
  mp = mass_h                        ! mass of a proton [g]
  c  = c_light                       ! speed of light    [cm/s]
  hp = hplanck                       ! Planck's constant [ergs*s]
  kb = kboltz                        ! boltzmann constant [erg/K]
  rc = 7.56e-15_RKIND                ! radiation constant [erg/cm^3/K^4]
  StBz  = 5.6704e-5_RKIND            ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  min_temp = 0.1_RKIND               ! minimum temperature [K]
  Cv    = 2.218056e12_RKIND * kb / ev2erg  ! used for problem 205
  if (ProbType == 405) then
     mol_weight = ev2erg/(gamma-1._RKIND)/2.218056e12_RKIND/mp  ! mean molecular mass
  else if (ProbType == 404) then
     mol_weight = 0.5_RKIND
  else
     mol_weight = 0.6_RKIND   ! mean molecular mass
  endif
  if (DualEnergy == 1) then
     KEconst = 0._RKIND
  else
     KEconst = 0.5_RKIND
  endif
  grey = 1._RKIND
  if (ESpectrum == -1)  grey = 0._RKIND

  if (Model /= 10) then
     write(0,*) 'gFLDProblem_LocalJac10: illegal Model =',Model
  endif

  ! ensure that we're running without chemistry
  if (Nchem /= 0) then
     write(0,*) 'gFLDProblem_LocalJac10: Nchem = ',Nchem,' and Model =', &
          Model,' do not match!'
     return
  endif
  
  ! iterate over the domain
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1
           
           ! set shortcut values for this spatial location
           ! converting densities from comoving to proper
           rho = rhoa(i,j,k)*DenUnits
           Er  = Era(i,j,k)*ErUnits
           
           ! compute temperature and derivatives
           eint = eca(i,j,k)*ecUnits + VelUnits*VelUnits*(eha(i,j,k)  &
                - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
           if ( ProbType /= 405 ) then
              T = (gamma-1._RKIND)*mol_weight*mp*eint/kb
              dT_dec = ecUnits*(gamma-1._RKIND)*mol_weight*mp/kb
           else
              T = eint/Cv
              dT_dec = ecUnits/Cv
           endif
           
           T = max(T,min_temp)
           
           ! compute opacities and derivatives
           kappaP = PmC0 * (rho/PmC1)**PmC2 * (T/PmC3)**PmC4
           kappaE = EmC0 * (rho/EmC1)**EmC2 * (T/EmC3)**EmC4
           dkappaP_dec = PmC0*PmC4 * (rho/PmC1)**PmC2 &
                * (T/PmC3)**(PmC4-1._RKIND)*dT_dec/PmC3
           dkappaE_dec = EmC0*EmC4 * (rho/EmC1)**EmC2 &
                * (T/EmC3)**(EmC4-1._RKIND)*dT_dec/EmC3
           
           ! compute emissivity and derivatives
           eta = kappaP*StBz/pi*T**4
           deta_dec = StBz/pi*T**4*dkappaP_dec &
                + 4._RKIND*kappaP*StBz/pi*T**3*dT_dec
           
           ! put it all together
           !   rhs_ec = (-2._RKIND*afac*ec + (c*kappaE*Er - 4._RKIND*pi*eta)/rho)/ecUnits
           ecjac_ec(i,j,k) = (ecjac_ec(i,j,k) - 2._RKIND*afac*ecUnits  &
                + (c*Er*dkappaE_dec - 4._RKIND*pi*deta_dec)/rho)/ecUnits
           ecjac_Er(i,j,k) = c*kappaE/rho*ErUnits/ecUnits
           
           !   rhs_Er = (4._RKIND*pi*eta - (afac*grey + c*kappaE)*Er)/ErUnits
           Erjac_Er(i,j,k) = -(afac*grey + c*kappaE)
           Erjac_ec(i,j,k) =  (4._RKIND*pi*deta_dec - c*Er*dkappaE_dec)/ErUnits
           
        enddo
     enddo
  enddo
  
  
  return
end subroutine gFLDProblem_LocalJac10
!=======================================================================

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
subroutine gFLDSplit_AnalyticChemistry(Er, ec, HI, HeI, HeII, Er0, ec0,  &
     HI0, HeI0, HeII0, dt, vx, vy, vz, rho, eh, src_ec, src_HI, src_HeI, &
     src_HeII, kappa, gamma, HFrac, Model, PType, DualEnergy, a, adot,   &
     CompA, CompXray, CompTemp, IsE, IsEsHI, IsEsHInu, IsEsHeI,          &
     IsEsHeInu, IsEsHeII, IsEsHeIInu, NTempBins, TempStart, TempEnd,     &
     k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb,      &
     ciHITb, ciHeITb, ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, &
     reHeIIITb, bremTb, aUnits, DenUnits, VelUnits, LenUnits, TimeUnits, &
     ErUnits, ecUnits, NiUnits, ecScale, Nchem, Nx, Ny, Nz, NGxl, NGxr,  &
     NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!
!  PURPOSE: Driver to call the correct routine depending on 'model'
!
!=======================================================================
  implicit none
#include "fortran_types.def"
    
!--------------
! argument declarations
  INTG_PREC, intent(in)  :: Model, PType, Nchem, NTempBins, DualEnergy
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a, adot
  R_PREC,    intent(in) :: dt, gamma, HFrac, TempStart, TempEnd
  R_PREC,    intent(in) :: CompA, CompXray, CompTemp
  R_PREC,    intent(in) :: DenUnits, VelUnits, TimeUnits, LenUnits,   &
       ErUnits, ecUnits, NiUnits, aUnits, ecScale
  R_PREC,    intent(in) :: IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, &
       IsEsHeII, IsEsHeIInu
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr),    &
       intent(in) :: Er, ec0, Er0, HI0, HeI0, HeII0, vx, vy, vz,    &
       rho, eh, src_ec, src_HI, src_HeI, src_HeII, kappa
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr),    &
       intent(out) :: ec, HI, HeI, HeII
  R_PREC, dimension(NTempBins), intent(in) :: k1Tb, k2Tb, k3Tb, k4Tb, &
       k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,      &
       ciHeIITb, ciHeISTb, reHIITb, reHeII1Tb, reHeII2Tb,           &
       reHeIIITb, bremTb

!=======================================================================

  ! initialize success/fail flag to fail
  ier = 0

  ! call the appropriate 'model'-specific routine
  if (Model == 1) then

     call gFLDSplit_AnalyticChemistry1(Er, ec, HI, HeI, HeII, Er0, ec0,  &
          HI0, HeI0, HeII0, dt, vx, vy, vz, rho, eh, src_ec, src_HI,     &
          src_HeI, src_HeII, gamma, HFrac, Model, DualEnergy, a, adot,   &
          CompA, CompXray, CompTemp, IsE, IsEsHI, IsEsHInu, IsEsHeI,     &
          IsEsHeInu, IsEsHeII, IsEsHeIInu, NTempBins, TempStart,         &
          TempEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb,  &
          ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb, reHIITb,        &
          reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUnits, DenUnits,     &
          VelUnits, LenUnits, TimeUnits, ErUnits, ecUnits, NiUnits,      &
          Nchem, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)

  else if (Model == 4) then
     
     call gFLDSplit_AnalyticChemistry4(Er, HI, Er0, HI0, dt, vx, vy, vz, &
          rho, eh, src_HI, gamma, HFrac, Model, DualEnergy, a, adot,     &
          IsE, IsEsHI, IsEsHInu, NTempBins, TempStart, TempEnd, k1Tb,    &
          DenUnits, VelUnits, ErUnits, NiUnits, ecScale, Nchem, Nx, Ny,  &
          Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)     

  else if (Model == 10) then
     
     call gFLDSplit_AnalyticChemistry10(Er, ec, Er0, ec0, dt, vx, vy,    &
          vz, rho, eh, src_ec, kappa, gamma, Model, PType, DualEnergy,   &
          a, adot, aUnits, DenUnits, VelUnits, LenUnits, TimeUnits,      &
          ErUnits, ecUnits, NiUnits, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, &
          NGzl, NGzr, ier)     

  else
     
     write(*,*) 'gFLDSplit_AnalyticChemistry: cannot process Model =',Model
     write(*,*) '  returning with an error'

  endif

  ! exit subroutine
  return

end subroutine gFLDSplit_AnalyticChemistry
!=======================================================================






subroutine gFLDSplit_AnalyticChemistry1(Er, ec, HI, HeI, HeII, Er0, ec0, &
     HI0, HeI0, HeII0, dt, vx, vy, vz, rho, eh, src_ec, src_HI, src_HeI, &
     src_HeII, gamma, HFrac, Model, DualEnergy, a, adot, CompA,          &
     CompXray, CompTemp, IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu,      &
     IsEsHeII, IsEsHeIInu, NTempBins, TempStart, TempEnd, k1Tb, k2Tb,    &
     k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, &
     ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb,       &
     bremTb, aUnits, DenUnits, VelUnits, LenUnits, TimeUnits, ErUnits,   &
     ecUnits, NiUnits, Nchem, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl,  &
     NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!
!  PURPOSE: Computes the solutions to the gas energy correction and 
!           chemistry equations, using a quasi-steady-state approximation 
!           and the resulting analytical solution of the relevant ODEs.
!
!  INPUTS:
!     Er         - radiation energy array (new)
!     Er0        - radiation energy array (old)
!     ec0        - fluid energy correction array (old)
!     HI0        - Hydrogen I number density array (old)
!     HeI0       - Helium I number density array (old)
!     HeII0      - Helium II number density array (old)
!     dt         - time step size
!     vx,vy,vz   - velocity arrays in each direction
!     rho        - density array
!     eh         - total fluid energy array
!     src_ec     - source function values for gas energy correction eq.
!     src_HI     - source function values for HI eq.
!     src_HeI    - source function values for HeI eq.
!     src_HeII   - source function values for HeII eq.
!     gamma      - constant in ideal gas law
!     HFrac      - percentage of mass composed of Hydrogen
!     Model      - flag denoting physical model to use
!     DualEnergy - flag denoting dual energy formalism
!     a          - cosmological expansion parameter
!     adot       - da/dt
!     CompA      - Compton cooling coefficient 1 (multiplier)
!     CompXray   - X-ray Compton heating coefficient
!     CompTemp   - X-ray Compton heating temperature 
!     IsE        - int_{nu0}^{inf} sigE dnu
!     IsEsHI     - int_{nu0}^{inf} sigE*sigHI dnu
!     IsEsHInu   - int_{nu0}^{inf} sigE*sigHI/nu dnu
!     IsEsHeI    - int_{nu0_HeI}^{inf} sigE*sigHeI dnu
!     IsEsHeInu  - int_{nu0_HeI}^{inf} sigE*sigHeI/nu dnu
!     IsEsHeII   - int_{nu0_HeII}^{inf} sigE*sigHeII dnu
!     IsEsHeIInu - int_{nu0_HeII}^{inf} sigE*sigHeII/nu dnu
!     NTempBins  - number of temperature bins for rate tables
!     TempStart  - starting temperature for rate table
!     TempEnd    - ending temperature for rate table
!     k*Tb       - chemistry rate tables
!     c*Tb, r*Tb - heating/cooling rate tables
!     *Units     - variable scaling constants
!     Nchem      - number of chemical species in simulation
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     ec      - fluid energy correction array (new)
!     HI      - Hydrogen I number density array (new)
!     HeI     - Helium I number density array (new)
!     HeII    - Helium II number density array (new)
!     ier     - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: Model, Nchem, NTempBins, DualEnergy
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a, adot
  R_PREC,    intent(in) :: dt, gamma, HFrac, TempStart, TempEnd
  R_PREC,    intent(in) :: CompA, CompXray, CompTemp
  R_PREC,    intent(in) :: DenUnits, VelUnits, TimeUnits, LenUnits,   &
       ErUnits, ecUnits, NiUnits, aUnits
  R_PREC,    intent(in) :: IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, &
       IsEsHeII, IsEsHeIInu
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr),    &
       intent(in) :: Er, ec0, Er0, HI0, HeI0, HeII0, vx, vy, vz,    &
       rho, eh, src_ec, src_HI, src_HeI, src_HeII
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr),    &
       intent(out) :: ec, HI, HeI, HeII
  R_PREC, dimension(NTempBins), intent(in) :: k1Tb, k2Tb, k3Tb, k4Tb, &
       k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,      &
       ciHeIITb, ciHeISTb, reHIITb, reHeII1Tb, reHeII2Tb,           &
       reHeIIITb, bremTb

!--------------
! locals
  INTG_PREC :: i, j, k, l, sweeps, sweeps2
  REAL*8  :: zr, Comp1, Comp2, CompX, CompT, KEconst, eint, gam
  REAL*8  :: HydFrac, deltat, deltax, int1, int2, int3, int4, int5, int6, int7
  REAL*8  :: TStart, TEnd, aUn, dUn, vUn, lUn, tUn, rUn, eUn, nUn
  REAL*8  :: res_Er, res_ec, res_HI, res_HeI, res_HeII, aval, dadt
  REAL*8  :: Ernew, ecnew, HInew, HeInew, HeIInew, change, lam, lam2, FPtol
  REAL*8  :: Erold, ecold, HIold, HeIold, HeIIold
  REAL*8  :: rhoval, ecsrc, HIsrc, HeIsrc, HeIIsrc
  
!=======================================================================

  ! initialize success/fail flag to success
  ier = 1

  ! we only have this enabled for Model 1 (case B HII recomb rate), 
  ! with chemistry
  if ((Nchem == 0) .or. (Model /= 1)) then
     write(*,*) 'AnalyticChemistry1 ERROR: only implemented for Model 1 w/ chem', &
          'called with Model =',Model,' and Nchem = ',Nchem
     ier = 0
     return
  endif

  ! compute shortcuts (need to compute these in double precision)
  gam = gamma
  HydFrac = HFrac
  deltat = dt
  deltax = 1.d0   ! send in a dummy value since we aren't using radiation
  int1 = IsE
  int2 = IsEsHI
  int3 = IsEsHInu
  int4 = IsEsHeI
  int5 = IsEsHeInu
  int6 = IsEsHeII
  int7 = IsEsHeIInu
  TStart = TempStart
  TEnd = TempEnd
  aUn = aUnits
  dUn = DenUnits
  vUn = VelUnits
  lUn = LenUnits
  tUn = TimeUnits
  rUn = ErUnits
  eUn = ecUnits
  nUn = NiUnits
  aval = a
  dadt = adot
  Comp1 = 0.d0
  Comp2 = 0.d0
  CompX = CompXray
  CompT = CompTemp
  zr = 1.d0/(a*aUnits) - 1.d0
  if (a*aUnits /= 1.d0) then
     Comp1 = CompA*(1.d0 + zr)**4
     Comp2 = 2.73d0*(1.d0 + zr)
  endif
  KEconst = 0.5d0
  if (DualEnergy == 1)  KEconst = 0.d0

  ! set fixed-point iteration parameters
  FPtol = 1.0d-8
  sweeps = 50
  sweeps2 = 10000
  lam = 1.d0
  lam2 = 0.1d0

  ! perform iteration based on chemistry

  !    no chemistry
  if (Nchem == 0) then

     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              Ernew = Er(i,j,k)
              ecnew = ec(i,j,k)
              Erold = Er0(i,j,k)
              ecold = ec0(i,j,k)
              rhoval = rho(i,j,k)
              ecsrc = src_ec(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! perform coarse fixed-point iteration to find analytical solution
              do l=0,sweeps
                 
                 ! call the local residual routine
                 call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,     &
                      res_HeI, res_HeII, Ernew, ecnew, 0.d0, 0.d0, 0.d0,        &
                      Erold, ecold, 0.d0, 0.d0, 0.d0, deltat, deltax, rhoval,   &
                      eint, 0.d0, ecsrc, 0.d0, 0.d0, 0.d0, gam, HydFrac,        &
                      Model, 0_IKIND, aval, dadt, Comp1, Comp2, CompX, CompT,   &
                      int1, int2, int3, int4, int5, int6, int7, NTempBins,      &
                      TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, &
                      ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb,   &
                      reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUn,    &
                      dUn, vUn, lUn, tUn, rUn, eUn, nUn, Nchem, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! update the current guesses
                 ecnew = ecnew - lam*res_ec
                 
                 ! compute the FP change in iteration, check for convergence
                 change = abs(res_ec)
                 if (change < FPtol)  exit
                 
              end do
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l=0,sweeps2
                    
                    ! call the local residual routine
                    call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,     &
                         res_HeI, res_HeII, Ernew, ecnew, 0.d0, 0.d0, 0.d0,        &
                         Erold, ecold, 0.d0, 0.d0, 0.d0, deltat, deltax, rhoval,   &
                         eint, 0.d0, ecsrc, 0.d0, 0.d0, 0.d0, gam, HydFrac,        &
                         Model, 0_IKIND, aval, dadt, Comp1, Comp2, CompX, CompT,   &
                         int1, int2, int3, int4, int5, int6, int7, NTempBins,      &
                         TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, &
                         ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb,   &
                         reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUn,    &
                         dUn, vUn, lUn, tUn, rUn, eUn, nUn, Nchem, ier)
                    
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! update the current guesses
                    ecnew = ecnew - lam2*res_ec
                    
                    ! compute the FP change in iteration, check for convergence
                    change = abs(res_ec)
                    if (change < FPtol)  exit
                    
                 end do
                 
              end if
              
              ! fill the relevant results
              ec(i,j,k) = ecnew
              
           end do
        end do
     end do

  !    Hydrogen chemistry
  elseif (Nchem == 1) then

     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              Ernew = Er(i,j,k)
              ecnew = ec(i,j,k)
              HInew = HI(i,j,k)
              Erold = Er0(i,j,k)
              ecold = ec0(i,j,k)
              HIold = HI0(i,j,k)
              rhoval = rho(i,j,k)
              ecsrc = src_ec(i,j,k)
              HIsrc = src_HI(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! perform coarse fixed-point iteration to find analytical solution
              do l=0,sweeps
                 
                 ! call the local residual routine
                 call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,     &
                      res_HeI, res_HeII, Ernew, ecnew, HInew, 0.d0, 0.d0,       &
                      Erold, ecold, HIold, 0.d0, 0.d0, deltat, deltax, rhoval,  &
                      eint, 0.d0, ecsrc, HIsrc, 0.d0, 0.d0, gam, HydFrac,       &
                      Model, 0_IKIND, aval, dadt, Comp1, Comp2, CompX, CompT,   &
                      int1, int2, int3, int4, int5, int6, int7, NTempBins,      &
                      TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, &
                      ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb,   &
                      reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUn,    &
                      dUn, vUn, lUn, tUn, rUn, eUn, nUn, Nchem, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! update the current guesses
                 ecnew = ecnew - lam*res_ec
                 HInew = HInew - lam*res_HI
                 
                 ! compute the FP change in iteration, check for convergence
                 change = max(abs(res_ec),abs(res_HI))
                 if (change < FPtol)  exit
                 
              end do
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l=0,sweeps2
                    
                    ! call the local residual routine
                    call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,     &
                         res_HeI, res_HeII, Ernew, ecnew, HInew, 0.d0, 0.d0,       &
                         Erold, ecold, HIold, 0.d0, 0.d0, deltat, deltax, rhoval,  &
                         eint, 0.d0, ecsrc, HIsrc, 0.d0, 0.d0, gam, HydFrac,       &
                         Model, 0_IKIND, aval, dadt, Comp1, Comp2, CompX, CompT,   &
                         int1, int2, int3, int4, int5, int6, int7, NTempBins,      &
                         TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, &
                         ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb,   &
                         reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUn,    &
                         dUn, vUn, lUn, tUn, rUn, eUn, nUn, Nchem, ier)
                    
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! update the current guesses
                    ecnew = ecnew - lam2*res_ec
                    HInew = HInew - lam2*res_HI
                    
                    ! compute the FP change in iteration, check for convergence
                    change = max(abs(res_ec),abs(res_HI))
                    if (change < FPtol)  exit
                    
                 end do
                 
              end if
              
              ! fill the relevant results
              ec(i,j,k) = ecnew
              HI(i,j,k) = HInew
              
           end do
        end do
     end do


  !    Hydrogen + Helium chemistry
  else

     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              Ernew   = Er(i,j,k)
              ecnew   = ec(i,j,k)
              HInew   = HI(i,j,k)
              HeInew  = HeI(i,j,k)
              HeIInew = HeII(i,j,k)
              Erold   = Er0(i,j,k)
              ecold   = ec0(i,j,k)
              HIold   = HI0(i,j,k)
              HeIold  = HeI0(i,j,k)
              HeIIold = HeII0(i,j,k)
              rhoval  = rho(i,j,k)
              ecsrc   = src_ec(i,j,k)
              HIsrc   = src_HI(i,j,k)
              HeIsrc  = src_HeI(i,j,k)
              HeIIsrc = src_HeII(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! perform coarse fixed-point iteration to find analytical solution
              do l=0,sweeps
                 
                 ! call the local residual routine
                 call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,     &
                      res_HeI, res_HeII, Ernew, ecnew, HInew, HeInew, HeIInew,  &
                      Erold, ecold, HIold, HeIold, HeIIold, deltat, deltax,     &
                      rhoval, eint, 0.d0, ecsrc, HIsrc, HeIsrc, HeIIsrc, gam,   &
                      HydFrac, Model, 0_IKIND, aval, dadt, Comp1, Comp2, CompX, &
                      CompT, int1, int2, int3, int4, int5, int6, int7,          &
                      NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb,    &
                      k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,         &
                      ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb,        &
                      reHeIIITb, bremTb, aUn, dUn, vUn, lUn, tUn, rUn, eUn,     &
                      nUn, Nchem, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! update the current guesses
                 ecnew   = ecnew   - lam*res_ec
                 HInew   = HInew   - lam*res_HI
                 HeInew  = HeInew  - lam*res_HeI
                 HeIInew = HeIInew - lam*res_HeII

                 ! compute the FP change in iteration, check for convergence
                 change = max(abs(res_ec),abs(res_HI))
                 change = max(change,abs(res_HeI))
                 change = max(change,abs(res_HeII))
                 if (change < FPtol)  exit
                 
              end do
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l=0,sweeps2
                    
                    ! call the local residual routine
                    call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,    &
                         res_HeI, res_HeII, Ernew, ecnew, HInew, HeInew, HeIInew, &
                         Erold, ecold, HIold, HeIold, HeIIold, deltat, deltax,    &
                         rhoval, eint, 0.d0, ecsrc, HIsrc, HeIsrc, HeIIsrc, gam,  &
                         HydFrac, Model, 0_IKIND, aval, dadt, Comp1, Comp2, CompX,&
                         CompT, int1, int2, int3, int4, int5, int6, int7,         &
                         NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb,   &
                         k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,        &
                         ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb,       &
                         reHeIIITb, bremTb, aUn, dUn, vUn, lUn, tUn, rUn, eUn,    &
                         nUn, Nchem, ier)
                    
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! update the current guesses
                    ecnew   = ecnew   - lam2*res_ec
                    HInew   = HInew   - lam2*res_HI
                    HeInew  = HeInew  - lam2*res_HeI
                    HeIInew = HeIInew - lam2*res_HeII
                    
                    ! compute the FP change in iteration, check for convergence
                    change = max(abs(res_ec),abs(res_HI))
                    change = max(change,abs(res_HeI))
                    change = max(change,abs(res_HeII))
                    if (change < FPtol)  exit
                    
                 end do
                 
              end if
              
              ! fill the relevant results
              ec(i,j,k)   = ecnew
              HI(i,j,k)   = HInew
              HeI(i,j,k)  = HeInew
              HeII(i,j,k) = HeIInew
              
           end do
        end do
     end do

  end if

  ! exit subroutine
  return

end subroutine gFLDSplit_AnalyticChemistry1
!=======================================================================






subroutine gFLDSplit_AnalyticChemistry4(Er, HI, Er0, HI0, dt, vx, vy, vz, &
     rho, eh, src_HI, gamma, HFrac, Model, DualEnergy, a, adot, IsE,      &
     IsEsHI, IsEsHInu, NTempBins, TempStart, TempEnd, k1Tb, DenUnits,     &
     VelUnits, ErUnits, NiUnits, ecScale, Nchem, Nx, Ny, Nz, NGxl, NGxr,  &
     NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       October 2009
!
!  PURPOSE: Computes the solution to the chemistry equations, using a 
!           quasi-steady-state approximation and the resulting 
!           analytical solution of the relevant ODEs.
!
!  INPUTS:
!     Er         - radiation energy array (new)
!     Er0        - radiation energy array (old)
!     ec0        - fluid energy correction array (old)
!     HI0        - Hydrogen I number density array (old)
!     HeI0       - Helium I number density array (old)
!     HeII0      - Helium II number density array (old)
!     dt         - time step size
!     vx,vy,vz   - velocity arrays in each direction
!     rho        - density array
!     eh         - total fluid energy array
!     src_HI     - source function values for HI eq.
!     src_HeI    - source function values for HeI eq.
!     src_HeII   - source function values for HeII eq.
!     gamma      - constant in ideal gas law
!     HFrac      - percentage of mass composed of Hydrogen
!     Model      - flag denoting physical model to use
!     DualEnergy - flag denoting dual energy formalism
!     a          - cosmological expansion parameter
!     adot       - da/dt
!     CompA      - Compton cooling coefficient 1 (multiplier)
!     CompXray   - X-ray Compton heating coefficient
!     CompTemp   - X-ray Compton heating temperature 
!     IsE        - int_{nu0}^{inf} sigE dnu
!     IsEsHI     - int_{nu0}^{inf} sigE*sigHI dnu
!     IsEsHInu   - int_{nu0}^{inf} sigE*sigHI/nu dnu
!     IsEsHeI    - int_{nu0_HeI}^{inf} sigE*sigHeI dnu
!     IsEsHeInu  - int_{nu0_HeI}^{inf} sigE*sigHeI/nu dnu
!     IsEsHeII   - int_{nu0_HeII}^{inf} sigE*sigHeII dnu
!     IsEsHeIInu - int_{nu0_HeII}^{inf} sigE*sigHeII/nu dnu
!     NTempBins  - number of temperature bins for rate tables
!     TempStart  - starting temperature for rate table
!     TempEnd    - ending temperature for rate table
!     k*Tb       - chemistry rate tables
!     c*Tb, r*Tb - heating/cooling rate tables
!     *Units     - variable scaling constants
!     ecScale    - holds the temperature for Model 4 (if cosmology enabled)
!     Nchem      - number of chemical species in simulation
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     HI      - Hydrogen I number density array (new)
!     HeI     - Helium I number density array (new)
!     HeII    - Helium II number density array (new)
!     ier     - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: Model, Nchem, NTempBins, DualEnergy
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a, adot
  R_PREC,    intent(in) :: dt, gamma, HFrac, TempStart, TempEnd
  R_PREC,    intent(in) :: DenUnits, VelUnits, ErUnits, NiUnits, ecScale
  R_PREC,    intent(in) :: IsE, IsEsHI, IsEsHInu
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in) :: Er, Er0, HI0, vx, vy, vz, rho, eh, src_HI
  R_PREC,    intent(out) :: HI(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  R_PREC,    intent(in) :: k1Tb(NTempBins)
       

!--------------
! locals
  INTG_PREC :: i, j, k, Tidx, Tidxp
  REAL*8  :: KEconst, eint, HydFrac, deltat, int1, int2, int3
  REAL*8  :: dUn, rUn, nUn, res_Er, res_HI, aval, dadt, Ernew, HInew
  REAL*8  :: Erold, HIold, rhoval, HIsrc, T, lamT, lTemp
  REAL*8  :: lTempS, lTempE, dlTemp, Tl, Tr, Tfac, k1, k2
  
!=======================================================================

  ! initialize success/fail flag to success
  ier = 1

  ! we only have this enabled for Model 1 (case B HII recomb rate), 
  ! with chemistry
  if ((Nchem /= 1) .or. (Model /= 4)) then
     write(*,*) 'AnalyticChemistry4 ERROR: only implemented for Nchem = 1', &
          'called with Model =',Model,' and Nchem = ',Nchem
     ier = 0
     return
  endif

  ! compute shortcuts (need to compute these in double precision)
  HydFrac = HFrac
  deltat = dt
  int1 = IsE
  int2 = IsEsHI
  int3 = IsEsHInu
  dUn = DenUnits
  rUn = ErUnits
  nUn = NiUnits
  aval = a
  dadt = adot
  KEconst = 0.5d0
  if (DualEnergy == 1)  KEconst = 0.d0

  ! first compute the fixed temperature (isothermal model)
  if (adot == 0.d0) then
     eint = VelUnits*VelUnits*(eh(1,1,1)                        &
          - KEconst*(vx(1,1,1)**2 + vy(1,1,1)**2 + vz(1,1,1)**2))
     T = eint*(gamma-1.d0)*0.6d0*mass_h/kboltz
     T = max(T,1.d0)
  else
     T = ecScale
  endif
  
  ! now compute the Hydrogen reaction rates
  lamT = 3.15614d5/T
  lTempS = log(TempStart)
  lTempE = log(TempEnd)
  dlTemp = (lTempE - lTempS)/(1.d0*NTempBins - 1.d0)
  lTemp = min(max(log(T), lTempS), lTempE)
  Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
  Tidxp = Tidx+1
  Tl = lTempS + (Tidx-1)*dlTemp
  Tr = lTempS +  Tidx*dlTemp
  Tfac = (lTemp - Tl)/(Tr - Tl)
  k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
  
  ! compute case B Hydrogen recombination coefficient 
  ! [Hui & Gnedin, 1997: RI^B_{HII}]
  !   (still need this because table holds case A coefficient!)
  k2 = 2.753d-14*lamT**(1.5d0) *                 &
       (1.d0+(lamT/2.74d0)**(0.407d0))**(-2.242d0)
  
  ! iterate over the domain
  do k=1,Nz
     do j=1,Ny
        do i=1,Nx
           
           ! get shortcut values at this location
           Ernew = Er(i,j,k)
           HInew = HI(i,j,k)
           Erold = Er0(i,j,k)
           HIold = HI0(i,j,k)
           rhoval = rho(i,j,k)
           HIsrc = src_HI(i,j,k)
           
           ! call the local residual routine
           call gFLDProblem_AnalyticLocResid4(res_Er, res_HI, Ernew, HInew,  &
                Erold, HIold, deltat, rhoval, 0.d0, HIsrc, HydFrac, 0_IKIND, &
                aval, dadt, int1, int2, int3, k1, k2, dUn, rUn, nUn, ier)
           
           ! check the error flag
           if (ier /= 1)  return
           
           ! fill the relevant results
           HI(i,j,k) = HInew - res_HI
           
        end do
     end do
  end do
  
  ! exit subroutine
  return

end subroutine gFLDSplit_AnalyticChemistry4
!=======================================================================






subroutine gFLDSplit_AnalyticChemistry10(Er, ec, Er0, ec0, dt, vx, vy,  &
     vz, rho, eh, src_ec, kappa, gamma, Model, PType, DualEnergy, a,    &
     adot, aUnits, DenUnits, VelUnits, LenUnits, TimeUnits, ErUnits,    &
     ecUnits, NiUnits, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!
!  PURPOSE: Computes the solutions to the gas energy correction and 
!           chemistry equations, using a quasi-steady-state approximation 
!           and the resulting analytical solution of the relevant ODEs.
!
!  INPUTS:
!     Er         - radiation energy array (new)
!     Er0        - radiation energy array (old)
!     ec0        - fluid energy correction array (old)
!     dt         - time step size
!     vx,vy,vz   - velocity arrays in each direction
!     rho        - density array
!     eh         - total fluid energy array
!     src_ec     - source function values for gas energy correction eq.
!     kappa      - opacity
!     gamma      - constant in ideal gas law
!     Model      - flag denoting physical model to use
!     PType      - problem type
!     DualEnergy - flag denoting dual energy formalism
!     a          - cosmological expansion parameter
!     adot       - da/dt
!     *Units     - variable scaling constants
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     ec      - fluid energy correction array (new)
!     ier     - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: Model, PType, DualEnergy
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a, adot
  R_PREC,    intent(in) :: dt, gamma
  R_PREC,    intent(in) :: DenUnits, VelUnits, TimeUnits, LenUnits, &
       ErUnits, ecUnits, aUnits, NiUnits
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr),  &
       intent(in) :: Er, ec0, Er0, vx, vy, vz, rho, eh, src_ec, kappa
  R_PREC, intent(out) :: ec(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)

!--------------
! locals
  INTG_PREC :: i, j, k, l
  REAL*8  :: zr, KEconst, eint, gam, deltat, aUn, dUn, vUn, lUn, tUn, rUn
  REAL*8  :: eUn, nUn, kap, res_Er, res_ec, aval, dadt, Ernew, ecnew
  REAL*8  :: Erold, ecold, rhoval, ecsrc
  
!=======================================================================


!!$  print *, 'entering AnalyticChemistry10'

  ! initialize success/fail flag to success
  ier = 1

  ! compute shortcuts (need to compute these in double precision)
  gam = gamma
  deltat = dt
  aUn = aUnits
  dUn = DenUnits
  vUn = VelUnits
  lUn = LenUnits
  tUn = TimeUnits
  rUn = ErUnits
  eUn = ecUnits
  nUn = NiUnits
  aval = a
  dadt = adot
  KEconst = 0.5d0
  if (DualEnergy == 1)  KEconst = 0.d0

  ! iterate over the domain
  do k=1,Nz
     do j=1,Ny
        do i=1,Nx
           
           ! get shortcut values at this location
           Ernew = Er(i,j,k)
           ecnew = ec(i,j,k)
           Erold = Er0(i,j,k)
           ecold = ec0(i,j,k)
           rhoval = rho(i,j,k)
           ecsrc = src_ec(i,j,k)
           kap = kappa(i,j,k)
           eint = vUn*vUn*(eh(i,j,k)                                  &
                - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
           
           ! call the local residual routine
           call gFLDProblem_AnalyticLTEResid(res_Er, res_ec, Ernew,      &
                ecnew, Erold, ecold, deltat, rhoval, eint, 0.d0,         &
                ecsrc, kap, gam, Model, PType, 0_IKIND, aval, dadt, aUn, &
                dUn, lUn, tUn, rUn, eUn, nUn, ier)
           
           ! check the error flag
           if (ier /= 1)  return
           
           ! update the current guesses
           ec(i,j,k) = ecnew - res_ec
             
        end do
     end do
  end do

  ! exit subroutine
  return

end subroutine gFLDSplit_AnalyticChemistry10
!=======================================================================






subroutine gFLDSplit_AnalyticInitGuess(Er, ec, HI, HeI, HeII, dt, vx,    &
     vy, vz, rho, eh, src_Er, src_ec, src_HI, src_HeI, src_HeII, gamma,  &
     HFrac, Model, ESpectrum, DualEnergy, a, adot, CompA, CompXray,      &
     CompTemp, IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII,      &
     IsEsHeIInu, NTempBins, TempStart, TempEnd, k1Tb, k2Tb, k3Tb, k4Tb,  &
     k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb,   &
     ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUnits, &
     DenUnits, VelUnits, LenUnits, TimeUnits, ErUnits, ecUnits, NiUnits, &
     ecScale, Nchem, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,     &
     NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!
!  PURPOSE: Computes an initial guess to the time-evolved PDE system 
!     using a Gauss-Seidel iteration on the local ODE reaction system.
!
!  INPUTS:
!     Er         - radiation energy (old time on input, new on output)
!     ec         - fluid energy correction array
!     HI         - Hydrogen I number density array
!     HeI        - Helium I number density array
!     HeII       - Helium II number density array
!     dt         - time step size
!     vx,vy,vz   - velocity arrays in each direction
!     rho        - density array
!     eh         - total fluid energy array
!     src_Er     - source function values for radiation equation
!     src_ec     - source function values for gas energy correction eq.
!     src_HI     - source function values for HI eq.
!     src_HeI    - source function values for HeI eq.
!     src_HeII   - source function values for HeII eq.
!     gamma      - constant in ideal gas law
!     HFrac      - percentage of mass composed of Hydrogen
!     Model      - flag denoting physical model to use
!     ESpectrum  - flag denoting choice of radiation energy spectrum
!     DualEnergy - flag denoting dual energy formalism
!     a          - cosmological expansion parameter
!     adot       - da/dt
!     CompA      - Compton cooling coefficient 1 (multiplier)
!     CompXray   - X-ray Compton heating coefficient
!     CompTemp   - X-ray Compton heating temperature 
!     IsE        - int_{nu0}^{inf} sigE dnu
!     IsEsHI     - int_{nu0}^{inf} sigE*sigHI dnu
!     IsEsHInu   - int_{nu0}^{inf} sigE*sigHI/nu dnu
!     IsEsHeI    - int_{nu0_HeI}^{inf} sigE*sigHeI dnu
!     IsEsHeInu  - int_{nu0_HeI}^{inf} sigE*sigHeI/nu dnu
!     IsEsHeII   - int_{nu0_HeII}^{inf} sigE*sigHeII dnu
!     IsEsHeIInu - int_{nu0_HeII}^{inf} sigE*sigHeII/nu dnu
!     NTempBins  - number of temperature bins for rate tables
!     TempStart  - starting temperature for rate table
!     TempEnd    - ending temperature for rate table
!     k*Tb       - chemistry rate tables
!     c*Tb, r*Tb - heating/cooling rate tables
!     *Units     - variable scaling constants
!     Nchem      - number of chemical species in simulation
!     Nx,Ny,Nz   - active mesh size in each direction
!     dx,dy,dz   - spatial mesh cell size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     Er   - guess at radiation energy density
!     ec   - guess at gas energy correction 
!     HI   - guess at HI density (unfilled if Nchem < 1)
!     HeI  - guess at HeI density (unfilled if Nchem < 3)
!     HeII - guess at HeII density (unfilled if Nchem < 3)
!     ier  - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: Model, Nchem, NTempBins, DualEnergy, ESpectrum
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a, adot
  R_PREC,    intent(in) :: dt, dx, dy, dz, gamma, HFrac, TempStart, TempEnd
  R_PREC,    intent(in) :: CompA, CompXray, CompTemp
  R_PREC,    intent(in) :: DenUnits, VelUnits, TimeUnits, LenUnits,   &
       ErUnits, ecUnits, NiUnits, aUnits, ecScale
  R_PREC,    intent(in) :: IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, &
       IsEsHeII, IsEsHeIInu
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) &
       :: vx, vy, vz, rho, eh, src_Er, src_ec, src_HI, src_HeI, src_HeII
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) ::  &
       ec, Er, HI, HeI, HeII
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) ::  &
       ecres, HIres, HeIres, HeIIres
  R_PREC, dimension(NTempBins), intent(in) :: k1Tb, k2Tb, k3Tb, k4Tb, &
       k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,      &
       ciHeIITb, ciHeISTb, reHIITb, reHeII1Tb, reHeII2Tb,           &
       reHeIIITb, bremTb

!--------------
! locals
  INTG_PREC :: ijk, i, j, k, l, l2, sweeps, sweeps2, lmax, lmin, lsum, Tidx, Tidxp
  REAL*8  :: chmin, chmax, chsum, eint, gam
  REAL*8  :: FPtol, lam, lam2, change, zr, Comp1, Comp2, CompX, CompT, KEconst
  REAL*8  :: HydFrac, deltat, deltax, int1, int2, int3, int4, int5, int6, int7
  REAL*8  :: TStart, TEnd, aUn, dUn, vUn, lUn, tUn, rUn, eUn, nUn, rhoval
  REAL*8  :: res_Er, res_ec, res_HI, res_HeI, res_HeII, aval, dadt
  REAL*8  :: Ernew, ecnew, HInew, HeInew, HeIInew
  REAL*8  :: Erold, ecold, HIold, HeIold, HeIIold
  REAL*8  :: Ersrc, ecsrc, HIsrc, HeIsrc, HeIIsrc
  REAL*8  :: T, lamT, lTemp, lTempS, lTempE, dlTemp, Tl, Tr, Tfac, k1, k2
  
!=======================================================================

  ! initialize success/fail flag to success
  ier = 1

  ! we only have this enabled for Model 1 (case B HII recomb rate), 
  ! with chemistry 
  if ((Nchem == 0) .or. ((Model /= 1) .and. (Model /= 4))) then
     write(*,*) 'AnalyticInitGuess ERROR: only implemented for Models 1 & 4', &
          ' w/ chem, called with Model =',Model,' and Nchem = ',Nchem
     ier = 0
     return
  endif

  ! compute shortcuts (need to compute these in double precision)
  gam = gamma
  HydFrac = HFrac
  deltat = dt
  deltax = (dx*dy*dz)**(1.d0/3.d0)
  int1 = IsE
  int2 = IsEsHI
  int3 = IsEsHInu
  int4 = IsEsHeI
  int5 = IsEsHeInu
  int6 = IsEsHeII
  int7 = IsEsHeIInu
  TStart = TempStart
  TEnd = TempEnd
  aUn = aUnits
  dUn = DenUnits
  vUn = VelUnits
  lUn = LenUnits
  tUn = TimeUnits
  rUn = ErUnits
  eUn = ecUnits
  nUn = NiUnits
  aval = a
  dadt = adot
  Comp1 = 0.d0
  Comp2 = 0.d0
  CompX = CompXray
  CompT = CompTemp
  zr = 1.d0/(a*aUnits) - 1.d0
  if (a*aUnits /= 1.d0) then
     Comp1 = CompA*(1.d0 + zr)**4
     Comp2 = 2.73d0*(1.d0 + zr)
  endif
  KEconst = 0.5d0
  if (DualEnergy == 1)  KEconst = 0.d0

  ! set some parameters for the Gauss-Seidel iteration
  FPtol = 1.0d-8
  sweeps = 50
  sweeps2 = 10000
  lam = 1.d0
  lam2 = 0.1d0


  ! perform method based on chemistry/model

  if (Model == 4) then

     ! first compute the fixed temperature (isothermal model)
     if (adot == 0.d0) then
        eint = vUn*vUn*(eh(1,1,1)                                  &
             - KEconst*(vx(1,1,1)**2 + vy(1,1,1)**2 + vz(1,1,1)**2))
        T = eint*(gamma-1.d0)*0.6d0*mass_h/kboltz
        T = max(T,1.d0)
     else
        T = ecScale
     endif

     ! now compute the Hydrogen reaction rates
     lamT = 3.15614d5/T
     lTempS = log(TempStart)
     lTempE = log(TempEnd)
     dlTemp = (lTempE - lTempS)/(1.d0*NTempBins - 1.d0)
     lTemp = min(max(log(T), lTempS), lTempE)
     Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
     Tidxp = Tidx+1
     Tl = lTempS + (Tidx-1)*dlTemp
     Tr = lTempS +  Tidx*dlTemp
     Tfac = (lTemp - Tl)/(Tr - Tl)
     k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
     
     ! compute case B Hydrogen recombination coefficient 
     ! [Hui & Gnedin, 1997: RI^B_{HII}]
     !   (still need this because table holds case A coefficient!)
     k2 = 2.753d-14*lamT**(1.5d0) *                 &
          (1.d0+(lamT/2.74d0)**(0.407d0))**(-2.242d0)

     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              Erold = Er(i,j,k)
              HIold = HI(i,j,k)
              Ernew = Erold
              HInew = HIold
              rhoval = rho(i,j,k)
              Ersrc = src_Er(i,j,k)
              HIsrc = src_HI(i,j,k)
              
              ! perform coarse fixed-point iteration
              do l=0,sweeps
                 
                 ! call the local residual routine
                 call gFLDProblem_AnalyticLocResid4(res_Er, res_HI, Ernew,   &
                      HInew, Erold, HIold, deltat, rhoval, Ersrc, HIsrc,     &
                      HydFrac, ESpectrum, aval, dadt, int1, int2, int3, k1,  &
                      k2, dUn, rUn, nUn, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! compute FP change in this iteration and check for convergence
                 change = max(abs(res_Er),abs(res_HI))
                 if (change < FPtol)  exit
                 
                 ! update the current guesses
                 Ernew = Ernew - lam*res_Er
                 HInew = HInew - lam*res_HI
                 
              end do  ! coarse fixed-point iteration
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l2=0,sweeps2
                    
                    ! call the local residual routine
                    call gFLDProblem_AnalyticLocResid4(res_Er, res_HI, Ernew,   &
                         HInew, Erold, HIold, deltat, rhoval, Ersrc, HIsrc,     &
                         HydFrac, ESpectrum, aval, dadt, int1, int2, int3, k1,  &
                         k2, dUn, rUn, nUn, ier)
                 
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! compute FP change in iteration and check for convergence
                    change = max(abs(res_Er),abs(res_HI))
                    if (change < FPtol)  exit
                 
                    ! update the current guesses
                    Ernew = Ernew - lam2*res_Er
                    HInew = HInew - lam2*res_HI
                 
                 end do  ! fine fixed-point iteration
                 
              end if

              ! update the outputs with the current values
              Er(i,j,k) = Ernew
              HI(i,j,k) = HInew

           end do
        end do
     end do

  else     ! model /= 4

  !    no chemistry
  if (Nchem == 0) then
  
     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              Erold = Er(i,j,k)
              ecold = ec(i,j,k)
              Ernew = Erold
              ecnew = ecold
              rhoval = rho(i,j,k)
              Ersrc = src_Er(i,j,k)
              ecsrc = src_ec(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                 &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! perform coarse fixed-point iteration
              do l=0,sweeps
                 
                 ! call the local residual routine
                 call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,   &
                      res_HeI, res_HeII, Ernew, ecnew, 0.d0, 0.d0, 0.d0,      &
                      Erold, ecold, 0.d0, 0.d0, 0.d0, deltat, deltax,         &
                      rhoval, eint, Ersrc, ecsrc, 0.d0, 0.d0, 0.d0, gam,      &
                      HydFrac, Model, ESpectrum, aval, dadt, Comp1, Comp2,    &
                      CompX, CompT, int1, int2, int3, int4, int5, int6, int7, &
                      NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb,  &
                      k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,       &
                      ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb,      &
                      reHeIIITb, bremTb, aUn, dUn, vUn, lUn, tUn, rUn, eUn,   &
                      nUn, Nchem, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! compute FP change in this iteration and check for convergence
                 change = max(abs(res_Er),abs(res_ec))
                 if (change < FPtol)  exit
                 
                 ! update the current guesses
                 Ernew = Ernew - lam*res_Er
                 ecnew = ecnew - lam*res_ec
                 
              end do  ! coarse fixed-point iteration
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l2=0,sweeps2
                    
                    ! call the local residual routine
                    call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,   &
                         res_HeI, res_HeII, Ernew, ecnew, 0.d0, 0.d0, 0.d0,      &
                         Erold, ecold, 0.d0, 0.d0, 0.d0, deltat, deltax,         &
                         rhoval, eint, Ersrc, ecsrc, 0.d0, 0.d0, 0.d0, gam,      &
                         HydFrac, Model, ESpectrum, aval, dadt, Comp1, Comp2,    &
                         CompX, CompT, int1, int2, int3, int4, int5, int6, int7, &
                         NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb,  &
                         k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,       &
                         ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb,      &
                         reHeIIITb, bremTb, aUn, dUn, vUn, lUn, tUn, rUn, eUn,   &
                         nUn, Nchem, ier)
                 
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! compute FP change in iteration and check for convergence
                    change = max(abs(res_Er),abs(res_ec))
                    if (change < FPtol)  exit
                 
                    ! update the current guesses
                    Ernew = Ernew - lam2*res_Er
                    ecnew = ecnew - lam2*res_ec
                 
                 end do  ! fine fixed-point iteration
                 
              end if

              ! update the outputs with the current values
              Er(i,j,k) = Ernew
              ec(i,j,k) = ecnew

           end do
        end do
     end do
     
  !    Hydrogen chemistry
  else if (Nchem == 1) then
  
     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              Erold = Er(i,j,k)
              ecold = ec(i,j,k)
              HIold = HI(i,j,k)
              Ernew = Erold
              ecnew = ecold
              HInew = HIold
              rhoval = rho(i,j,k)
              Ersrc = src_Er(i,j,k)
              ecsrc = src_ec(i,j,k)
              HIsrc = src_HI(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                 &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! initialize iteration counters
              l = 0
              l2 = 0

              ! perform coarse fixed-point iteration
              do l=0,sweeps
                 
                 ! call the local residual routine
                 call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,   &
                      res_HeI, res_HeII, Ernew, ecnew, HInew, 0.d0, 0.d0,     &
                      Erold, ecold, HIold, 0.d0, 0.d0, deltat, deltax,        &
                      rhoval, eint, Ersrc, ecsrc, HIsrc, 0.d0, 0.d0, gam,     &
                      HydFrac, Model, ESpectrum, aval, dadt, Comp1, Comp2,    &
                      CompX, CompT, int1, int2, int3, int4, int5, int6, int7, &
                      NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb,  &
                      k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,       &
                      ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb,      &
                      reHeIIITb, bremTb, aUn, dUn, vUn, lUn, tUn, rUn, eUn,   &
                      nUn, Nchem, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! compute FP change in this iteration and check for convergence
                 change = max(abs(res_Er),abs(res_ec))
                 change = max(change,abs(res_HI))
                 if (change < FPtol)  exit
                 
                 ! update the current guesses
                 Ernew = Ernew - lam*res_Er
                 ecnew = ecnew - lam*res_ec
                 HInew = HInew - lam*res_HI
                 
              end do  ! coarse fixed-point iteration
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l2=0,sweeps2
                    
                    ! call the local residual routine
                    call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,   &
                         res_HeI, res_HeII, Ernew, ecnew, HInew, 0.d0, 0.d0,     &
                         Erold, ecold, HIold, 0.d0, 0.d0, deltat, deltax,        &
                         rhoval, eint, Ersrc, ecsrc, HIsrc, 0.d0, 0.d0, gam,     &
                         HydFrac, Model, ESpectrum, aval, dadt, Comp1, Comp2,    &
                         CompX, CompT, int1, int2, int3, int4, int5, int6, int7, &
                         NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb,  &
                         k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,       &
                         ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb,      &
                         reHeIIITb, bremTb, aUn, dUn, vUn, lUn, tUn, rUn, eUn,   &
                         nUn, Nchem, ier)
                 
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! compute FP change in iteration and check for convergence
                    change = max(abs(res_Er),abs(res_ec))
                    change = max(change,abs(res_HI))
                    if (change < FPtol)  exit
                 
                    ! update the current guesses
                    Ernew = Ernew - lam2*res_Er
                    ecnew = ecnew - lam2*res_ec
                    HInew = HInew - lam2*res_HI
                 
                 end do  ! fine fixed-point iteration
                 
              end if

              ! update the outputs with the current values
              Er(i,j,k) = Ernew
              ec(i,j,k) = ecnew
              HI(i,j,k) = HInew

           end do
        end do
     end do
     
  !    Hydrogen + Helium chemistry
  else
  
     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              Erold   = Er(i,j,k)
              ecold   = ec(i,j,k)
              HIold   = HI(i,j,k)
              HeIold  = HeI(i,j,k)
              HeIIold = HeII(i,j,k)
              Ernew   = Erold
              ecnew   = ecold
              HInew   = HIold
              HeInew  = HeIold
              HeIInew = HeIIold
              rhoval  = rho(i,j,k)
              Ersrc   = src_Er(i,j,k)
              ecsrc   = src_ec(i,j,k)
              HIsrc   = src_HI(i,j,k)
              HeIsrc  = src_HeI(i,j,k)
              HeIIsrc = src_HeII(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                 &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! initialize iteration counters
              l = 0
              l2 = 0
              
              ! perform coarse fixed-point iteration
              do l=0,sweeps

                 ! call the local residual routine
                 call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,   &
                      res_HeI, res_HeII, Ernew, ecnew, HInew, HeInew,         &
                      HeIInew, Erold, ecold, HIold, HeIold, HeIIold, deltat,  &
                      deltax, rhoval, eint, Ersrc, ecsrc, HIsrc, HeIsrc,      &
                      HeIIsrc, gam, HydFrac, Model, ESpectrum, aval, dadt,    &
                      Comp1, Comp2, CompX, CompT, int1, int2, int3, int4,     &
                      int5, int6, int7, NTempBins, TStart, TEnd, k1Tb, k2Tb,  &
                      k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb,      &
                      ciHITb, ciHeITb, ciHeISTb, ciHeIITb, reHIITb,           &
                      reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUn, dUn, vUn, &
                      lUn, tUn, rUn, eUn, nUn, Nchem, ier)
              
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! compute FP change in iteration and check for convergence
                 change = max(abs(res_Er),abs(res_ec))
                 change = max(change,abs(res_HI))
                 change = max(change,abs(res_HeI))
                 change = max(change,abs(res_HeII))
                 if (change < FPtol)  exit
              
                 ! update the current guesses
                 Ernew   = Ernew   - lam*res_Er
                 ecnew   = ecnew   - lam*res_ec
                 HInew   = HInew   - lam*res_HI
                 HeInew  = HeInew  - lam*res_HeI
                 HeIInew = HeIInew - lam*res_HeII
                 
              end do  ! coarse fixed-point iteration
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l2=0,sweeps2
                    
                    ! call the local residual routine
                    call gFLDProblem_AnalyticLocResid1(res_Er, res_ec, res_HI,   &
                         res_HeI, res_HeII, Ernew, ecnew, HInew, HeInew,         &
                         HeIInew, Erold, ecold, HIold, HeIold, HeIIold, deltat,  &
                         deltax, rhoval, eint, Ersrc, ecsrc, HIsrc, HeIsrc,      &
                         HeIIsrc, gam, HydFrac, Model, ESpectrum, aval, dadt,    &
                         Comp1, Comp2, CompX, CompT, int1, int2, int3, int4,     &
                         int5, int6, int7, NTempBins, TStart, TEnd, k1Tb, k2Tb,  &
                         k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb,      &
                         ciHITb, ciHeITb, ciHeISTb, ciHeIITb, reHIITb,           &
                         reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUn, dUn, vUn, &
                         lUn, tUn, rUn, eUn, nUn, Nchem, ier)
                    
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! compute FP change in iteration and check for convergence
                    change = max(abs(res_Er),abs(res_ec))
                    change = max(change,abs(res_HI))
                    change = max(change,abs(res_HeI))
                    change = max(change,abs(res_HeII))
                    if (change < FPtol)  exit
                    
                    ! update the current guesses
                    Ernew   = Ernew   - lam2*res_Er
                    ecnew   = ecnew   - lam2*res_ec
                    HInew   = HInew   - lam2*res_HI
                    HeInew  = HeInew  - lam2*res_HeI
                    HeIInew = HeIInew - lam2*res_HeII
                    
                 end do  ! fine fixed-point iteration
                 
              end if

              ! update the outputs with the current values
              Er(i,j,k)   = Ernew
              ec(i,j,k)   = ecnew
              HI(i,j,k)   = HInew
              HeI(i,j,k)  = HeInew
              HeII(i,j,k) = HeIInew
              
           end do
        end do
     end do
     
  end if  ! nchem
  end if  ! model

  ! exit subroutine
  return

end subroutine gFLDSplit_AnalyticInitGuess
!=======================================================================






subroutine gFLDProblem_AnalyticLTEResid(Erres, ecres, Er, ec, Er0, ec0,    &
     dt, rho, eint, src_Er, src_ec, kappa, gamma, Model, PType, ESpectrum, &
     a, adot, aUnits, DenUnits, LenUnits, TimeUnits, ErUnits, ecUnits,     &
     NiUnits, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       May 2009
  !
  !  PURPOSE: Computes the deviation from the analytical solutions to a 
  !           quasi-steady-state approximation of the (Er,ec) system in 
  !           local thermodynamic equilibrium (LTE).
  !
  !  INPUTS:
  !     Er         - radiation energy 
  !     ec         - fluid energy correction
  !     Er0        - old radiation energy 
  !     ec0        - old fluid energy correction
  !     dt         - time step size
  !     rho        - density 
  !     eint       - fluid internal energy
  !     src_Er     - source function values for radiation equation
  !     src_ec     - source function values for gas energy correction eq.
  !     kappa      - opacity 
  !     gamma      - constant in ideal gas law
  !     Model      - flag denoting physical model to use
  !     PType      - problem type
  !     ESpectrum  - flag denoting choice of radiation energy spectrum
  !     a          - cosmological expansion parameter
  !     adot       - da/dt
  !     *Units     - variable scaling constants
  !
  !  OUTPUT ARGUMENTS: 
  !     Erres   - radiation energy density residual
  !     ecres   - gas energy correction residual
  !     HIres   - Hydrogen I number density residual
  !     HeIres  - Helium I number density residual
  !     HeIIres - Helium II number density residual
  !     ier     - success/failure flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: Model, PType, ESpectrum
  INTG_PREC, intent(out) :: ier
  REAL*8,  intent(in) :: a, adot
  REAL*8,  intent(in) :: dt, gamma, kappa
  REAL*8,  intent(in) :: DenUnits, TimeUnits, LenUnits, &
       ErUnits, ecUnits, aUnits, NiUnits
  REAL*8,  intent(in) :: Er, ec, Er0, ec0, rho, eint, src_Er, src_ec
  REAL*8,  intent(out) :: Erres, ecres

  !--------------
  ! locals
  REAL*8 :: afac, StBz, c, kb, mp, min_temp, min_rad, grey, eta
  REAL*8 :: ecval, Erval, rhoval, Eranal, ecanal
  REAL*8 :: T, P1, Q1, P2, Q2

  !=======================================================================

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  ecres = 0.d0
  Erres = 0.d0

  ! initialize constants
  afac = adot/a
  StBz  = 5.6704d-5            ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  c  = c_light                 ! speed of light [cm/s]
  kb = kboltz                  ! boltzmann constant [erg/K]
  mp = mass_h                  ! Mass of a proton [g]
  min_temp = 1.d0              ! minimum temperature [K]
  min_rad  = 0.d0              ! minimum radiation density [ergs/cm^3]
  grey = 1.d0                  ! grey vs monochromatic coeff for eqns
  if (ESpectrum == -1)  grey = 0.d0


  ! get shortcut values time-centered variables
  Erval = (Er0 + Er)*0.5d0*ErUnits
  ecval = (ec0 + ec)*0.5d0*ecUnits
  rhoval = rho*DenUnits

  ! compute temperature and ODE terms
  if (PType == 405) then   ! Lowrie & Edwards radiating shock
     T = (eint+ecval)/2.218056d12/kb*ev2erg
  else
     if (PType == 404) then
        T = (gamma-1.0)*0.5d0*mp*(eint+ecval)/kb
     else
        T = (gamma-1.0)*0.6d0*mp*(eint+ecval)/kb
     endif
  endif
  T = max(T, min_temp)
  eta = 4.d0*kappa*NiUnits*StBz*T**4

  ! compute radiation ODE rates
  P1 = grey*afac + c*kappa*NiUnits
  Q1 = (src_Er + eta)/ErUnits

  ! compute gas ODE rates
  P2 = 2.d0*afac
  Q2 = (src_ec + (c*kappa*NiUnits*Erval - eta)/rhoval)/ecUnits
  
  ! compute quasi-steady-state solution for Er, place in Eranal
  if (abs(P1) < 1.0d-14) then
     Eranal = Er0 + Q1*dt
  else
     if (P1*dt > 7.0d2) then
        Eranal = Q1/P1
     else
        Eranal = (Er0 - Q1/P1)*exp(-P1*dt) + Q1/P1
     end if
  end if
  

  ! compute quasi-steady-state solution for ec, place in ecanal
  if (abs(P2) < 1.0d-14) then
     ecanal = ec0 + Q2*dt
  else
     if (P2*dt > 7.0d2) then
        ecanal = Q2/P2
     else
        ecanal = (ec0 - Q2/P2)*exp(-P2*dt) + Q2/P2
     end if
  end if

  ! enforce bounds
  Eranal = max(Eranal, min_rad/ErUnits)
  ecanal = max(ecanal, -eint/ecUnits)

  ! compute residuals
  Erres = Er - Eranal
  ecres = ec - ecanal

  ! check some things
  !    this statement checks if Eranal = NaN
  if (Eranal /= Eranal) then
     print *,'NaN encountered in AnalyticLTEResid (Er)!!'
     print '(2(A,es12.5))', '   P =',P1,', Q =',Q1
     ier = 0
     return
  end if
  
  !    this statement checks if ecanal = NaN
  if (ecanal /= ecanal) then
     print *,'NaN encountered in AnalyticLTEResid (ec)!!'
     print '(2(A,es12.5))', '   P =',P2,', Q =',Q2
     ier = 0
     return
  end if
  
  ! exit subroutine
  return

end subroutine gFLDProblem_AnalyticLTEResid
!=======================================================================

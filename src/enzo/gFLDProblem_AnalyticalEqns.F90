#include "fortran.def"
#include "phys_const.def"
!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
! Copyright 2009 Laboratory for Computational Astrophysics
! Copyright 2009 Regents of the University of California
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine gFLDProblem_AnalyticResid(ecres, HIres, HeIres, HeIIres, Er,  &
     ec, HI, HeI, HeII, Er0, ec0, HI0, HeI0, HeII0, dt, vx, vy, vz, rho, &
     eh, src_Er, src_ec, src_HI, src_HeI, src_HeII, gamma, HFrac, Model, &
     DualEnergy, a, adot, CompA, Comp_xray, Comp_temp, IsE, IsEsHI,      &
     IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII, IsEsHeIInu, NTempBins,      &
     TempStart, TempEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb,     &
     ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb, reHIITb,    &
     reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUnits, DenUnits,          &
     VelUnits, LenUnits, TimeUnits, ErUnits, ecUnits, NiUnits, ecScale,  &
     Nchem, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       May 2009
!
!  PURPOSE: Computes the residuals
!                ecres = ec - ec_analytical
!                HIres = HI - HI_analytical
!                HeIres = HeI - HeI_analytical
!                HeIIres = HeII - HeII_analytical
!     where the *_analytical come from using a quasi-steady-state
!     approximation and the resulting analytical solution of the 
!     relevant ODEs.
!
!  INPUTS:
!     Er         - radiation energy array
!     ec         - fluid energy correction array
!     HI         - Hydrogen I number density array
!     HeI        - Helium I number density array
!     HeII       - Helium II number density array
!     Er0        - radiation energy array (old)
!     ec0        - fluid energy correction array (old)
!     HI0        - Hydrogen I number density array (old)
!     HeI0       - Helium I number density array (old)
!     HeII0      - Helium II number density array (old)
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
!     DualEnergy - flag denoting dual energy formalism
!     a          - cosmological expansion parameter
!     adot       - da/dt
!     CompA      - Compton cooling coefficient 1 (multiplier)
!     Comp_xray  - X-ray Compton heating coefficient
!     Comp_temp  - X-ray Compton heating temperature 
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
!     ecres   - gas energy correction residual 
!     HIres   - HI density residual (unfilled if Nchem < 1)
!     HeIres  - HeI density residual (unfilled if Nchem < 3)
!     HeIIres - HeII density residual (unfilled if Nchem < 3)
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
  R_PREC, intent(in) :: dt, gamma, HFrac, TempStart, TempEnd
  R_PREC, intent(in) :: CompA, Comp_xray, Comp_temp
  R_PREC, intent(in) :: DenUnits, VelUnits, TimeUnits, LenUnits,    &
       ErUnits, ecUnits, NiUnits, aUnits, ecScale
  R_PREC, intent(in) :: IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu,  &
       IsEsHeII, IsEsHeIInu
  R_PREC, intent(in), dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       ec, Er, HI, HeI, HeII, ec0, Er0, HI0, HeI0, HeII0, vx, vy, &
       vz, rho, eh, src_Er, src_ec, src_HI, src_HeI, src_HeII
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) ::  &
       ecres, HIres, HeIres, HeIIres
  R_PREC, intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, &
       k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeIITb,        &
       ciHeISTb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb

!--------------
! locals
  INTG_PREC :: i, j, k, l, sweeps, sweeps2, Tidx, Tidxp
  REAL*8 :: zr, Comp1, Comp2, CompX, CompT, KEconst, eint, gam
  REAL*8 :: HydFrac, deltat, deltax, int1, int2, int3, int4, int5, int6, int7
  REAL*8 :: TStart, TEnd, aUn, dUn, vUn, lUn, tUn, rUn, eUn, nUn
  REAL*8 :: res_Er, res_ec, res_HI, res_HeI, res_HeII, aval, dadt
  REAL*8 :: Ernew, ecnew, HInew, HeInew, HeIInew, change, lam, lam2, FPtol
  REAL*8 :: Erold, ecold, HIold, HeIold, HeIIold
  REAL*8 :: rhoval, ecsrc, HIsrc, HeIsrc, HeIIsrc
  REAL*8 :: T, lamT, lTemp, lTempS, lTempE, dlTemp, Tl, Tr, Tfac, k1, k2
  
!=======================================================================

  ! initialize success/fail flag to success
  ier = 1

  ! initialize outputs to have all zero values
  ecres = 0._RKIND
  if (Nchem > 0)  HIres = 0._RKIND
  if (Nchem > 1) then
     HeIres = 0._RKIND
     HeIIres = 0._RKIND
  end if

  ! we only have this enabled for Models 1 & 4 (case B HII recomb rate), 
  ! with chemistry
  if ((Nchem == 0) .or. ((Model /= 1) .and. (Model /= 4))) then
     write(*,*) 'AnalyticResid ERROR: only implemented for Models 1 & 4', &
          ' w/ chem, called with Model =',Model,' and Nchem = ',Nchem
     ier = 0
     return
  endif

  ! compute shortcuts
  gam = gamma
  HydFrac = HFrac
  deltat = dt
  deltax = 1.d0   ! send in a dummy value, since residual does not use Er
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
  CompX = Comp_xray
  CompT = Comp_temp
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


  ! perform iteration based on chemistry/model

  if (Model == 4) then

     ! first compute the fixed temperature (isothermal model)
     if (adot == 0.d0) then
        eint = vUn*vUn*(eh(1,1,1)                                  &
             - KEconst*(vx(1,1,1)**2 + vy(1,1,1)**2 + vz(1,1,1)**2))
        T = eint*(gamma-1.d0)*0.6d0*mass_h/kboltz
        T = max(1.d0*T,1.d0)
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
                 
              ! fill the relevant residuals
              HIres(i,j,k) = res_HI
              
           end do
        end do
     end do

  else     ! model /= 4

  !    Hydrogen chemistry
  if (Nchem == 1) then

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
              
              ! fill the relevant residuals
              ecres(i,j,k) = ec(i,j,k) - ecnew
              HIres(i,j,k) = HI(i,j,k) - HInew
              
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
              
              ! fill the relevant residuals
              ecres(i,j,k)   = ec(i,j,k)   - ecnew
              HIres(i,j,k)   = HI(i,j,k)   - HInew
              HeIres(i,j,k)  = HeI(i,j,k)  - HeInew
              HeIIres(i,j,k) = HeII(i,j,k) - HeIInew
              
           end do
        end do
     end do

  end if  ! nchem
  end if  ! model

  ! exit subroutine
  return

end subroutine gFLDProblem_AnalyticResid
!=======================================================================






subroutine gFLDProblem_AnalyticInitGuess(Er, ec, HI, HeI, HeII, dt, vx,  &
     vy, vz, rho, eh, src_Er, src_ec, src_HI, src_HeI, src_HeII, gamma,  &
     HFrac, Model, ESpectrum, DualEnergy, a, adot, CompA, Comp_xray,     &
     Comp_temp, IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII,     &
     IsEsHeIInu, NTempBins, TempStart, TempEnd, k1Tb, k2Tb, k3Tb, k4Tb,  &
     k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb,   &
     ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, aUnits, &
     DenUnits, VelUnits, LenUnits, TimeUnits, ErUnits, ecUnits, NiUnits, &
     ecScale, Nchem, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,     &
     NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       May 2009
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
!     Comp_xray  - X-ray Compton heating coefficient
!     Comp_temp  - X-ray Compton heating temperature 
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
  R_PREC, intent(in) :: dt, dx, dy, dz, gamma, HFrac, TempStart, TempEnd
  R_PREC, intent(in) :: CompA, Comp_xray, Comp_temp
  R_PREC, intent(in) :: DenUnits, VelUnits, TimeUnits, LenUnits,     &
       ErUnits, ecUnits, NiUnits, aUnits, ecScale
  R_PREC, intent(in) :: IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu,   &
       IsEsHeII, IsEsHeIInu
  R_PREC, intent(in),                                                &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) ::  &
       vx, vy, vz, rho, eh, src_Er, src_ec, src_HI, src_HeI, src_HeII
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       ec, Er, HI, HeI, HeII
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       ecres, HIres, HeIres, HeIIres
  R_PREC, intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb,      &
       k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb,        &
       ciHeITb, ciHeIITb, ciHeISTb, reHIITb, reHeII1Tb, reHeII2Tb, &
       reHeIIITb, bremTb

!--------------
! locals
  INTG_PREC :: i, j, k, l, l2, sweeps, sweeps2, lmax, lmin, lsum, Tidx, Tidxp
  REAL*8 :: chmin, chmax, chsum, eint, gam
  REAL*8 :: FPtol, lam, lam2, change, zr, Comp1, Comp2, CompX, CompT, KEconst
  REAL*8 :: HydFrac, deltat, deltax, int1, int2, int3, int4, int5, int6, int7
  REAL*8 :: TStart, TEnd, aUn, dUn, vUn, lUn, tUn, rUn, eUn, nUn, rhoval
  REAL*8 :: res_Er, res_ec, res_HI, res_HeI, res_HeII, aval, dadt
  REAL*8 :: Ernew, ecnew, HInew, HeInew, HeIInew
  REAL*8 :: Erold, ecold, HIold, HeIold, HeIIold
  REAL*8 :: Ersrc, ecsrc, HIsrc, HeIsrc, HeIIsrc
  REAL*8 :: T, lamT, lTemp, lTempS, lTempE, dlTemp, Tl, Tr, Tfac, k1, k2
  
!=======================================================================

  ! initialize success/fail flag to success
  ier = 1

  ! we only have this enabled for Models 1 & 4 (case B HII recomb rate), 
  ! with chemistry 
  if ((Nchem == 0) .or. ((Model /= 1) .and. (Model /= 4))) then
     write(*,*) 'AnalyticInitGuess ERROR: only implemented for Models 1 & 4', &
          ' w/ chem, called with Model =',Model,' and Nchem = ',Nchem
     ier = 0
     return
  endif

  ! compute shortcuts
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
  CompX = Comp_xray
  CompT = Comp_temp
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

  !    Hydrogen chemistry
  if (Nchem == 1) then
  
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

end subroutine gFLDProblem_AnalyticInitGuess
!=======================================================================






subroutine gFLDProblem_AnalyticLocResid1(Erres, ecres, HIres, HeIres,  &
     HeIIres, Er, ec, HI, HeI, HeII, Er0, ec0, HI0, HeI0, HeII0, dt,   &
     dx, rho, eint, src_Er, src_ec, src_HI, src_HeI, src_HeII, gamma,  &
     HFrac, Model, ESpectrum, a, adot, Comp1, Comp2, Comp_xray,        &
     Comp_temp, IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII,   &
     IsEsHeIInu, NTempBins, TempStart, TempEnd, k1Tb, k2Tb, k3Tb,      &
     k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,     &
     ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb,     &
     bremTb, aUnits, DenUnits, VelUnits, LenUnits, TimeUnits, ErUnits, &
     ecUnits, NiUnits, Nchem, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       May 2009
  !
  !  PURPOSE: 
  !
  !  INPUTS:
  !     Er         - radiation energy array
  !     ec         - fluid energy correction array
  !     HI         - Hydrogen I number density array
  !     HeI        - Helium I number density array
  !     HeII       - Helium II number density array
  !     Er0        - old radiation energy array
  !     ec0        - old fluid energy correction array
  !     HI0        - old Hydrogen I number density array
  !     HeI0       - old Helium I number density array
  !     HeII0      - old Helium II number density array
  !     dt         - time step size
  !     dx         - spatial mesh size
  !     vx,vy,vz   - velocity arrays in each direction
  !     rho        - density array
  !     eing       - fluid internal energy
  !     src_Er     - source function values for radiation equation
  !     src_ec     - source function values for gas energy correction eq.
  !     src_HI     - source function values for HI eq.
  !     src_HeI    - source function values for HeI eq.
  !     src_HeII   - source function values for HeII eq.
  !     gamma      - constant in ideal gas law
  !     HFrac      - percentage of mass composed of Hydrogen
  !     Model      - flag denoting physical model to use
  !     ESpectrum  - flag denoting choice of radiation energy spectrum
  !     a          - cosmological expansion parameter
  !     adot       - da/dt
  !     Comp1      - Compton cooling coefficient 1
  !     Comp2      - Compton cooling coefficient 2
  !     Comp_xray  - X-ray Compton heating coefficient
  !     Comp_temp  - X-ray Compton heating temperature 
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
  INTG_PREC, intent(in)  :: Model, Nchem, NTempBins, ESpectrum
  INTG_PREC, intent(out) :: ier
  REAL*8, intent(in) :: a, adot
  REAL*8, intent(in) :: dt, dx, gamma, HFrac, TempStart, TempEnd
  REAL*8, intent(in) :: Comp1, Comp2, Comp_xray, Comp_temp
  REAL*8, intent(in) :: DenUnits, VelUnits, TimeUnits, LenUnits,    &
       ErUnits, ecUnits, NiUnits, aUnits
  REAL*8, intent(in) :: IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu,  &
       IsEsHeII, IsEsHeIInu
  REAL*8, intent(in) :: Er, ec, HI, HeI, HeII, Er0, ec0, HI0, HeI0, &
       HeII0, rho, eint, src_Er, src_ec, src_HI, src_HeI, src_HeII
  REAL*8, intent(out) :: Erres, ecres, HIres, HeIres, HeIIres
  R_PREC, intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb, k4Tb, &
       k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,      &
       ciHeIITb, ciHeISTb, reHIITb, reHeII1Tb, reHeII2Tb,           &
       reHeIIITb, bremTb

  !--------------
  ! locals
  INTG_PREC :: Tidx, Tidxp
  REAL*8 :: afac, c, hp, kb, mp, lTempS, lTempE, dlTemp
  REAL*8 :: nu0_HI, nu0_HeI, nu0_HeII, HIconst, HeIconst, HeIIconst
  REAL*8 :: min_temp, min_ni, min_rad
  REAL*8 :: HIval, HeIval, HeIIval, ecval, Erval, rhoval
  REAL*8 :: nH, nHI, nHII, nHe, nHeI, nHeII, nHeIII, ne
  REAL*8 :: HIanal, HeIanal, HeIIanal, Eranal, ecanal
  REAL*8 :: grey, T, lamT, lTemp, Tl, Tr, Tfac, kE, P1, Q1, cond, dx_sc
  REAL*8 :: wts(7), taus(7), vals(7), fvals(7), ival
  REAL*8 :: ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeII, ciHeIS
  REAL*8 :: reHII, reHeII1, reHeII2, reHeIII, brem, G, Lambda
  REAL*8 :: P2, Q2, k1, k2, k3, k4, k5, k6, G_HI, G_HeI, G_HeII
  REAL*8 :: aHI, bHI, cHI, dHI, kHI, expArg, expVHI, sqD, rt1, rt2
  REAL*8 :: aHeI, bHeI, cHeI, dHeI, kHeI, expVHeI
  REAL*8 :: aHeII, bHeII, cHeII, dHeII, kHeII, expVHeII

  !=======================================================================

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  ecres = 0._RKIND
  if (Nchem > 0)  HIres = 0._RKIND
  if (Nchem > 1) then
     HeIres = 0._RKIND
     HeIIres = 0._RKIND
  end if

  ! we only have this enabled for Model 1 (case B HII recomb rate), which
  ! requires chemistry (at least Hydrogen, if not also Helium)
  if ((Model /= 1) .or. (Nchem == 0)) then
     write(*,*) 'AnalyticLocResid1 ERROR: only implemented for Model 1 w/ chem', &
          'called with Model =',Model,' and Nchem =',Nchem
     ier = 0
     return
  endif


  ! initialize constants
  afac = adot/a
  dx_sc = dx*LenUnits/a
  c  = c_light                 ! speed of light [cm/s]
  hp = hplanck                 ! Planck's constant [ergs*s]
  kb = kboltz                  ! boltzmann constant [erg/K]
  mp = mass_h                  ! Mass of a proton [g]
  nu0_HI = 13.6_RKIND*ev2erg/hp    ! ionization frequency of HI   [hz]
  nu0_HeI = 24.6_RKIND*ev2erg/hp   ! ionization frequency of HeI  [hz]
  nu0_HeII = 54.4_RKIND*ev2erg/hp  ! ionization frequency of HeII [hz]
  min_temp = 1._RKIND              ! minimum temperature [K]
  min_ni   = 0._RKIND              ! minimum chem number density [cm^{-3}]
  min_rad  = 0._RKIND              ! minimum radiation density [ergs/cm^3]
  HIconst   = c*(IsEsHI - nu0_HI*IsEsHInu)/IsE
  HeIconst  = c*(IsEsHeI - nu0_HeI*IsEsHeInu)/IsE
  HeIIconst = c*(IsEsHeII - nu0_HeII*IsEsHeIInu)/IsE
  grey = 1._RKIND                  ! grey vs monochromatic coeff for eqns
  if (ESpectrum == -1)  grey = 0._RKIND

  ! radiation time integral weights and integration nodes
!!$  wts  = (/ 1._RKIND, 4._RKIND, 2._RKIND, 4._RKIND, 2._RKIND, 4._RKIND, 1._RKIND   /)*dt/18._RKIND
!!$  taus = (/ 6._RKIND, 5._RKIND, 4._RKIND, 3._RKIND, 2._RKIND, 1._RKIND, 1.e-11_RKIND /)*dt/6._RKIND

  ! lookup table constants
  lTempS = log(TempStart)
  lTempE = log(TempEnd)
  dlTemp = (lTempE - lTempS)/(1._RKIND*NTempBins - 1._RKIND)


  ! decide between none vs Hydrogen vs Hydrogen+Helium chemistry
  !    Hydrogen chemistry
  if (Nchem == 1) then

     ! get shortcut values time-centered variables
     Erval = (Er0 + Er)*0.5_RKIND
     ecval = (ec0 + ec)*0.5_RKIND
     HIval = (HI0 + HI)*0.5_RKIND
     HeIval = (HeI0 + HeI)*0.5_RKIND
     HeIIval = (HeII0 + HeII)*0.5_RKIND
     rhoval = rho*DenUnits
     nH = Hfrac*rhoval/mp
     nHI = HIval*NiUnits
     nHII = max(nH - nHI, 0.d0)
     ne = nHII

     ! compute temperature and ODE terms
     T = (gamma-1._RKIND)*rhoval/(nHI+nHII+ne)*(eint+ecval*ecUnits)/kb
     T = max(T,min_temp)
     lamT = 3.15614e5_RKIND/T
     lTemp = min(max(log(T), lTempS), lTempE)
     Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
     Tidxp = Tidx+1
     Tl = lTempS + (Tidx-1)*dlTemp
     Tr = lTempS +  Tidx*dlTemp
     Tfac = (lTemp - Tl)/(Tr - Tl)

     ! compute radiation ODE rates
     kE = HIval*NiUnits*IsEsHI/IsE
     P1 = grey*afac + c*kE
     Q1 = src_Er/ErUnits
     cond = 16._RKIND*c/kE/3._RKIND

     ! compute gas ODE rates
     ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
     ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
     reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
     brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac
     G = nHI/rhoval*Erval*ErUnits*HIconst
     Lambda = ne/rhoval*(ceHI*nHI + ciHI*nHI + reHII*nHII        &
          + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp) + brem*nHII)
     P2 = 2._RKIND*afac
     Q2 = (G - Lambda + src_ec)/ecUnits

     ! compute HI chemistry ODE rates
     k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
     k2 = 2.753e-14_RKIND*lamT**(1.5_RKIND) *                 &
          (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-2.242_RKIND)
     G_HI   = c*Erval*ErUnits/hp*IsEsHInu/IsE
     aHI = (k1+k2)*NiUnits
     bHI = -(k1 + 2._RKIND*k2)*nH - G_HI
     cHI = (k2*nH*nH + src_HI)/NiUnits



    ! compute quasi-steady-state solution for Er, place in Eranal
    if (abs(P1) < 1.e-14_RKIND) then
       Eranal = Er0 + Q1*dt
    else
       if (P1*dt > 7.e2_RKIND) then
          Eranal = Q1/P1
       else
          Eranal = (Er0 - Q1/P1)*exp(-P1*dt) + Q1/P1
       end if
    end if
!!$     ! The following analytical solution is derived using the Green's function
!!$     ! solution for the forced, constant-coefficient heat equation, where we 
!!$     ! only consider the spatial domain of dependence to be this cell (and 
!!$     ! ignore flow from farther distances).  While this is only a crude 
!!$     ! approximation to the solution, it is nonetheless very fast and provides
!!$     ! the exact analytical solution for the cell containing an emissivity 
!!$     ! source; as this is where the local nonlinear problem is hardest to solve, 
!!$     ! it should prove quite beneficial.
!!$     !    perform the time integral only if source is nonzero
!!$     Eranal = 0._RKIND
!!$     if (abs(Q1) > 1.d-13) then
!!$        ! evaluate the Green's function integrand at tau values, 
!!$        ! and combine to form the time-integral
!!$        vals = dx_sc/sqrt(cond*taus)
!!$        call erf_vec(fvals, vals, 7_IKIND)
!!$        ival = sum(wts * exp(-P1*taus) * fvals**3)
!!$
!!$        ! compute contribution due to sources
!!$        Eranal = ival*Q1
!!$!!$     else
!!$!!$        vals(1) = dx_sc/sqrt(cond*dt)
!!$!!$        call erf_vec(fvals, vals, 1_IKIND)
!!$     end if
!!$     !    compute contribution due to previous time solution
!!$!!$     !    (re-use fvals since it was last called with tau=0)
!!$!!$     Eranal = Eranal + Er0*exp(-P1*dt)*fvals(1)**3
!!$     Eranal = Eranal + Er0*exp(-P1*dt)
     

     ! compute quasi-steady-state solution for ec, place in ecanal
     if (abs(P2) < 1.e-14_RKIND) then
        ecanal = ec0 + Q2*dt
     else
        if (P2*dt > 7.e2_RKIND) then
           ecanal = Q2/P2
        else
           ecanal = (ec0 - Q2/P2)*exp(-P2*dt) + Q2/P2
        end if
     end if


     ! compute quasi-steady-state solution for HI, place in HIres
     dHI = bHI*bHI - 4._RKIND*aHI*cHI

     !    if the second-order terms are very small, treat it as a linear ODE
     if (abs(aHI*HIval**2/(cHI-bHI*HIval)) < 1.e-8_RKIND) then
        if (bHI*dt < -7.e2_RKIND) then
           HIanal = -cHI/bHI
        else
           HIanal = (HI0 + cHI/bHI)*exp(bHI*dt) - cHI/bHI
        end if

     !    otherwise, go ahead with the quadratic ODE solution
     else

        ! no roots to quadratic
        if (dHI/bHI/bHI < -1.e-8_RKIND) then
           sqD = sqrt(-dHI)
           kHI = 2._RKIND/sqD*atan((2._RKIND*aHI*HI0+bHI)/sqD)
           HIanal = (sqD*tan((dt+kHI)*0.5_RKIND*sqD)-bHI)*0.5_RKIND/aHI

        ! two roots to quadratic
        elseif (dHI/bHI/bHI > 1.e-8_RKIND) then
           rt1 = (-bHI+sqrt(dHI))*0.5_RKIND/aHI
           rt2 = (-bHI-sqrt(dHI))*0.5_RKIND/aHI 
           kHI = log(abs((HI0-rt2)/(HI0-rt1)))/aHI/(rt2-rt1)
           expArg = aHI*(dt+kHI)*(rt2-rt1)
           if (expArg < -7.e2_RKIND) then
              HIanal = rt2
           elseif (expArg > 7.e2_RKIND) then
              HIanal = rt1
           else
              expVHI = exp(expArg)
              if ((HI0-rt1)*(HI0-rt2) > 0._RKIND) then
                 HIanal = (rt2-rt1*expVHI)/(1._RKIND-expVHI)
              else
                 HIanal = (rt2+rt1*expVHI)/(1._RKIND+expVHI)
              end if
           end if

        ! one double root to quadratic
        else
           rt1 = -bHI*0.5_RKIND/aHI
           kHI = 1._RKIND/aHI/(rt1-HI0)
           HIanal = rt1 - 1._RKIND/aHI/(dt+kHI)

        end if  ! roots
     end if  ! quadratic


     ! enforce bounds
     Eranal = max(Eranal, min_rad/ErUnits)
     T = (gamma-1.d0)*rhoval/(nHI+nHII+ne)*(eint+ecanal*ecUnits)/kb
     if (T < min_temp)  &
          ecanal = (min_temp/(gamma-1._RKIND)/rhoval*(nHI+nHII+ne)*kb-eint)/ecUnits
     HIanal = min(HIanal, nH/NiUnits)
     HIanal = max(HIanal, min_ni/NiUnits)


     ! compute residuals
     Erres = Er - Eranal
     ecres = ec - ecanal
     HIres = HI - HIanal

     ! check some things
     !    this statement checks if Eranal = NaN
     if (Eranal /= Eranal) then
        print *,'NaN encountered in AnalyticLocResid1 (Er)!!'
        print '(2(A,es12.5))', '   P =',P1,', Q =',Q1
        print '(2(A,es12.5))', '   dt =',dt,', Ersrc =',src_Er
        print '(2(A,es12.5))', '   grey =',grey,', kE =',kE
        print '(2(A,es12.5))', '   Er0 =',Er0,', ErUnits =',ErUnits
        ier = 0
        return
     end if

     !    this statement checks if ecanal = NaN
     if (ecanal /= ecanal) then
        print *,'NaN encountered in AnalyticLocResid1 (ec)!!'
        print '(2(A,es12.5))', '   P =',P2,', Q =',Q2
        print '(2(A,es12.5))', '   dt =',dt,', ecsrc =',src_ec
        print '(2(A,es12.5))', '   G =',G,', Lambda =',Lambda
        print '(2(A,es12.5))', '   ec0 =',ec0,', ecUnits =',ecUnits
        ier = 0
        return
     end if

     !    this statement checks if HIanal = NaN
     if (HIanal /= HIanal) then
        print *,'NaN encountered in AnalyticLocResid1 (HI)!!'
        print '(2(A,es12.5))', '   A =',aHI,', B =',bHI
        print '(2(A,es12.5))', '   C =',cHI,', D =',dHI
        if (dHI/bHI/bHI < -1.0d-8) then
           print '(1(A,es12.5))', '   sqD =',sqrt(-dHI)
        elseif (dHI/bHI/bHI > 1.0d-8) then
           print '(3(A,es12.5))', '   rt1 =',(-bHI+sqrt(dHI))*0.5_RKIND/aHI, &
                ', rt2 =',(-bHI-sqrt(dHI))*0.5_RKIND/aHI,', expVHI =',expVHI
        else
           print '(1(A,es12.5))', '   rt1 =',-bHI*0.5_RKIND/aHI
        end if
        print '(2(A,es12.5))', '   K =',kHI,', G_HI =',G_HI
        print '(2(A,es12.5))', '  HIval =',HIval,', HI0 =',HI0
        print '(2(A,es12.5))', '  nH =',nH,', dt =',dt
        print '(2(A,es12.5))', '  NiUnits =',NiUnits,',  T =',T
        print '(2(A,es12.5))', '  k1 =',k1,', k2 =',k2
        ier = 0
        return
     end if


  !    Hydrogen + Helium chemistry
  else

     ! get shortcut values time-centered variables
     Erval = (Er0 + Er)*0.5_RKIND
     ecval = (ec0 + ec)*0.5_RKIND
     HIval = (HI0 + HI)*0.5_RKIND
     HeIval = (HeI0 + HeI)*0.5_RKIND
     HeIIval = (HeII0 + HeII)*0.5_RKIND
     rhoval = rho*DenUnits
     nH = Hfrac*rhoval/mp
     nHI = HIval*NiUnits
     nHII = max(nH - nHI, 0.d0)
     nHe = (1.d0-Hfrac)*rhoval/mp
     nHeI = HeIval*NiUnits
     nHeII = HeIIval*NiUnits
     nHeIII = max(nHe - nHeI - nHeII, 0.d0)
     ne = nHII + 0.25d0*nHeII + 0.5d0*nHeIII

     ! compute temperature and ODE terms
     T = (gamma-1._RKIND)*rhoval/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne) &
          * (eint+ecval*ecUnits)/kb
     T = max(T,min_temp)
     lamT = 3.15614d5/T
     lTemp = min(max(log(T), lTempS), lTempE)
     Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
     Tidxp = Tidx+1
     Tl = lTempS + (Tidx-1)*dlTemp
     Tr = lTempS +  Tidx*dlTemp
     Tfac = (lTemp - Tl)/(Tr - Tl)

     ! compute radiation ODE rates
     kE = (HIval*IsEsHI/IsE + HeIval*IsEsHeI/IsE + HeIIval*IsEsHeII/IsE)*NiUnits
     P1  = grey*afac + c*kE
     Q1  = src_Er/ErUnits
     cond = 16._RKIND*c/kE/3._RKIND

     ! compute gas ODE rates
     ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
     ceHeI = ceHeITb(Tidx) + (ceHeITb(Tidxp) - ceHeITb(Tidx))*Tfac
     ceHeII = ceHeIITb(Tidx) + (ceHeIITb(Tidxp) - ceHeIITb(Tidx))*Tfac
     ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
     ciHeI = ciHeITb(Tidx) + (ciHeITb(Tidxp) - ciHeITb(Tidx))*Tfac
     ciHeII = ciHeIITb(Tidx) + (ciHeIITb(Tidxp) - ciHeIITb(Tidx))*Tfac
     ciHeIS = ciHeISTb(Tidx) + (ciHeISTb(Tidxp) - ciHeISTb(Tidx))*Tfac
     reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
     reHeII1 = reHeII1Tb(Tidx) + (reHeII1Tb(Tidxp) - reHeII1Tb(Tidx))*Tfac
     reHeII2 = reHeII2Tb(Tidx) + (reHeII2Tb(Tidxp) - reHeII2Tb(Tidx))*Tfac
     reHeIII = reHeIIITb(Tidx) + (reHeIIITb(Tidxp) - reHeIIITb(Tidx))*Tfac
     brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac
     G = Erval*ErUnits/rhoval*(nHI*HIconst + nHeI*HeIconst + nHeII*HeIIconst)
     Lambda = ne/rhoval*(ceHI*nHI + ciHI*nHI + reHII*nHII        &
          + Comp1*(T-Comp2) + Comp_xray*(T-Comp_temp)            &
          + brem*(nHII+nHeII/4._RKIND+nHeIII)                        &
          + 0.25_RKIND*(ceHeI*nHeII*ne + ceHeII*nHeII + ciHeI*nHeI   &
          + ciHeII*nHeII + ciHeIS*nHeII*ne + reHeII1*nHeII       &
          + reHeII2*nHeII + reHeIII*nHeIII))
     P2 = 2._RKIND*afac
     Q2 = (G - Lambda + src_ec)/ecUnits

     ! compute chemistry ODE rates
     k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
     k2 = 2.753e-14_RKIND*lamT**(1.5_RKIND) *                 &
          (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-2.242_RKIND)
     k3 = k3Tb(Tidx) + (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac
     k4 = k4Tb(Tidx) + (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac
     k5 = k5Tb(Tidx) + (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac
     k6 = k6Tb(Tidx) + (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac
     G_HI   = c*Erval*ErUnits/hp*IsEsHInu/IsE
     G_HeI  = c*Erval*ErUnits/hp*IsEsHeInu/IsE
     G_HeII = c*Erval*ErUnits/hp*IsEsHeIInu/IsE
     aHI = (k1+k2)*NiUnits
     bHI = -(k1 + 2._RKIND*k2)*nH - (k1+k2)*(nHeII + 2._RKIND*nHeIII)/4._RKIND - G_HI
     cHI = (k2*nH*nH + k2*nH*(nHeII + 2._RKIND*nHeIII)/4._RKIND + src_HI)/NiUnits
     aHeI = k3*0.5_RKIND*NiUnits
     bHeI = -G_HeI - k3*(nHII + 0.5_RKIND*nHe - 0.25_RKIND*nHeII) - k4*0.5_RKIND*nHeII
     cHeI = (k4*nHeII*(nHII + 0.5_RKIND*nHe - 0.25_RKIND*nHeII) + src_HeI)/NiUnits
     aHeII = 0.25_RKIND*(k4+k5+k6)*NiUnits
     bHeII = -G_HeII - (k4+k5+k6)*(nHII + 0.5_RKIND*(nHe-nHeI)) &
          - 0.25_RKIND*(k3*nHeI + k6*(nHe-nHeI))
     cHeII = ((k3*nHeI + k6*(nHe-nHeI))*(nHII + 0.5_RKIND*(nHe-nHeI)) &
          + nHeI*G_HeI + src_HeII)/NiUnits



     ! compute quasi-steady-state solution for Er, place in Eranal
     if (abs(P1) < 1.e-14_RKIND) then
        Eranal = Er0 + Q1*dt
     else
        if (P1*dt > 7.e2_RKIND) then
           Eranal = Q1/P1
        else
           Eranal = (Er0 - Q1/P1)*exp(-P1*dt) + Q1/P1
        end if
     end if
!!$     !    perform time integral only if sources are nonzero
!!$     Eranal = 0._RKIND
!!$     if (abs(Q1) > 1.d-13) then
!!$        ! evaluate the Green's function integrand at tau values, 
!!$        ! and combine to form the time-integral
!!$        vals = dx_sc/sqrt(cond*taus)
!!$        call erf_vec(fvals, vals, 7_IKIND)
!!$        ival = sum(wts * exp(-P1*taus) * fvals**3)
!!$
!!$        ! compute contribution due to sources
!!$        Eranal = ival*Q1
!!$!!$     else
!!$!!$        vals(1) = dx_sc/sqrt(cond*dt)
!!$!!$        call erf_vec(fvals, vals, 1_IKIND)
!!$     end if
!!$     !    compute contribution due to previous time solution
!!$!!$     !    (re-use fvals since it was last called with tau=0)
!!$!!$     Eranal = Eranal + Er0*exp(-P1*dt)*fvals(1)**3
!!$     Eranal = Eranal + Er0*exp(-P1*dt)


     ! compute quasi-steady-state solution for ec, place in ecanal
     if (abs(P2) < 1.e-14_RKIND) then
        ecanal = ec0 + Q2*dt
     else
        if (P2*dt > 7.e2_RKIND) then
           ecanal = Q2/P2
        else
           ecanal = (ec0 - Q2/P2)*exp(-P2*dt) + Q2/P2
        end if
     end if


     ! compute quasi-steady-state solution for HI, place in HIres
     dHI = bHI*bHI - 4._RKIND*aHI*cHI

     !    if the second-order terms are very small, treat it as a linear ODE
     if (abs(aHI*HIval**2/(cHI-bHI*HIval)) < 1.e-8_RKIND) then
        if (bHI*dt < -7.e2_RKIND) then
           HIanal = -cHI/bHI
        else
           HIanal = (HI0 + cHI/bHI)*exp(bHI*dt) - cHI/bHI
        end if

     !    otherwise, go ahead with the quadratic ODE solution
     else

        ! no roots to quadratic
        if (dHI/bHI/bHI < -1.e-8_RKIND) then
           sqD = sqrt(-dHI)
           kHI = 2._RKIND/sqD*atan((2._RKIND*aHI*HI0+bHI)/sqD)
           HIanal = (sqD*tan((dt+kHI)*0.5_RKIND*sqD)-bHI)*0.5_RKIND/aHI

        ! two roots to quadratic
        elseif (dHI/bHI/bHI > 1.e-8_RKIND) then
           rt1 = (-bHI+sqrt(dHI))*0.5_RKIND/aHI
           rt2 = (-bHI-sqrt(dHI))*0.5_RKIND/aHI 
           kHI = log(abs((HI0-rt2)/(HI0-rt1)))/aHI/(rt2-rt1)
           expArg = aHI*(dt+kHI)*(rt2-rt1)
           if (expArg < -7.e2_RKIND) then
              HIanal = rt2
           elseif (expArg > 7.e2_RKIND) then
              HIanal = rt1
           else
              expVHI = exp(expArg)
              if ((HI0-rt1)*(HI0-rt2) > 0._RKIND) then
                 HIanal = (rt2-rt1*expVHI)/(1._RKIND-expVHI)
              else
                 HIanal = (rt2+rt1*expVHI)/(1._RKIND+expVHI)
              end if
           end if

        ! one double root to quadratic
        else
           rt1 = -bHI*0.5_RKIND/aHI
           kHI = 1._RKIND/aHI/(rt1-HI0)
           HIanal = rt1 - 1._RKIND/aHI/(dt+kHI)

        end if  ! roots
     end if  ! quadratic


     ! compute quasi-steady-state solution for HeI, place in HeIres
     dHeI = bHeI*bHeI - 4._RKIND*aHeI*cHeI

     !    if the second-order terms are very small, treat it as a linear ODE
     if (abs(aHeI*HeIval**2/(cHeI-bHeI*HeIval)) < 1.e-8_RKIND) then
        if (bHeI*dt < -7.e2_RKIND) then
           HeIanal = -cHeI/bHeI
        else
           HeIanal = (HeI0 + cHeI/bHeI)*exp(bHeI*dt) - cHeI/bHeI
        end if

     !    otherwise, go ahead with the quadratic ODE solution
     else

        ! no roots to quadratic
        if (dHeI/bHeI/bHeI < -1.e-8_RKIND) then
           sqD = sqrt(-dHeI)
           kHeI = 2._RKIND/sqD*atan((2._RKIND*aHeI*HeI0+bHeI)/sqD)
           HeIanal = (sqD*tan((dt+kHeI)*0.5_RKIND*sqD)-bHeI)*0.5_RKIND/aHeI

        ! two roots to quadratic
        elseif (dHeI/bHeI/bHeI > 1.e-8_RKIND) then
           rt1 = (-bHeI+sqrt(dHeI))*0.5_RKIND/aHeI
           rt2 = (-bHeI-sqrt(dHeI))*0.5_RKIND/aHeI 
           kHeI = log(abs((HeI0-rt2)/(HeI0-rt1)))/aHeI/(rt2-rt1)
           expArg = aHeI*(dt+kHeI)*(rt2-rt1)
           if (expArg < -7.0d2) then
              HeIanal = rt2
           elseif (expArg > 7.0d2) then
              HeIanal = rt1
           else
              expVHeI = exp(expArg)
              if ((HeI0-rt1)*(HeI0-rt2) > 0._RKIND) then
                 HeIanal = (rt2-rt1*expVHeI)/(1._RKIND-expVHeI)
              else
                 HeIanal = (rt2+rt1*expVHeI)/(1._RKIND+expVHeI)
              end if
           end if

        ! one double root to quadratic
        else
           rt1 = -bHeI*0.5_RKIND/aHeI
           kHeI = 1._RKIND/aHeI/(rt1-HeI0)
           HeIanal = rt1 - 1._RKIND/aHeI/(dt+kHeI)

        end if  ! roots
     end if  ! quadratic


     ! compute quasi-steady-state solution for HeII, place in HeIIres
     dHeII = bHeII*bHeII - 4._RKIND*aHeII*cHeII

     !    if the second-order terms are very small, treat it as a linear ODE
     if (abs(aHeII*HeIIval**2/(cHeII-bHeII*HeIIval)) < 1.e-8_RKIND) then
        if (bHeII*dt < -7.e2_RKIND) then
           HeIIanal = -cHeII/bHeII
        else
           HeIIanal = (HeII0 + cHeII/bHeII)*exp(bHeII*dt) - cHeII/bHeII
        end if

     !    otherwise, go ahead with the quadratic ODE solution
     else

        ! no roots to quadratic
        if (dHeII/bHeII/bHeII < -1.e-8_RKIND) then
           sqD = sqrt(-dHeII)
           kHeII = 2._RKIND/sqD*atan((2._RKIND*aHeII*HeII0+bHeII)/sqD)
           HeIIanal = (sqD*tan((dt+kHeII)*0.5_RKIND*sqD)-bHeII)*0.5_RKIND/aHeII

        ! two roots to quadratic
        elseif (dHeII/bHeII/bHeII > 1.e-8_RKIND) then
           rt1 = (-bHeII+sqrt(dHeII))*0.5_RKIND/aHeII
           rt2 = (-bHeII-sqrt(dHeII))*0.5_RKIND/aHeII 
           kHeII = log(abs((HeII0-rt2)/(HeII0-rt1)))/aHeII/(rt2-rt1)
           expArg = aHeII*(dt+kHeII)*(rt2-rt1)
           if (expArg < -7.e2_RKIND) then
              HeIIanal = rt2
           elseif (expArg > 7.e2_RKIND) then
              HeIIanal = rt1
           else
              expVHeII = exp(expArg)
              if ((HeII0-rt1)*(HeII0-rt2) > 0._RKIND) then
                 HeIIanal = (rt2-rt1*expVHeII)/(1._RKIND-expVHeII)
              else
                 HeIIanal = (rt2+rt1*expVHeII)/(1._RKIND+expVHeII)
              end if
           end if

        ! one double root to quadratic
        else
           rt1 = -bHeII*0.5_RKIND/aHeII
           kHeII = 1._RKIND/aHeII/(rt1-HeII0)
           HeIIanal = rt1 - 1._RKIND/aHeII/(dt+kHeII)

        end if  ! roots
     end if  ! quadratic


     ! enforce bounds
     Eranal = max(Eranal, min_rad/ErUnits)
     T = (gamma-1._RKIND)*rhoval/(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne) &
          *(eint+ecanal*ecUnits)/kb
     if (T < min_temp)  &
          ecanal = (min_temp/(gamma-1._RKIND)/rhoval &
          *(0.25_RKIND*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)*kb-eint)/ecUnits
     HIanal = min(HIanal, nH/NiUnits)
     HIanal = max(HIanal, min_ni/NiUnits)
     HeIanal = min(HeIanal, nHe/NiUnits)
     HeIanal = max(HeIanal, min_ni/NiUnits)
     HeIIanal = min(HeIIanal, nHe/NiUnits)
     HeIIanal = max(HeIIanal, min_ni/NiUnits)


     ! compute residuals
     Erres = Er - Eranal
     ecres = ec - ecanal
     HIres = HI - HIanal
     HeIres = HeI - HeIanal
     HeIIres = HeII - HeIIanal

     ! check some things
     !    this statement checks if Eranal = NaN
     if (Eranal /= Eranal) then
        print *,'NaN encountered in AnalyticLocResid1 (Er)!!'
        print '(2(A,es12.5))', '   P =',P1,', Q =',Q1
        print '(2(A,es12.5))', '   dt =',dt,', Ersrc =',src_Er
        print '(2(A,es12.5))', '   grey =',grey,', kE =',kE
        print '(2(A,es12.5))', '   Er0 =',Er0,', ErUnits =',ErUnits
        ier = 0
        return
     end if

     !    this statement checks if ecanal = NaN
     if (ecanal /= ecanal) then
        print *,'NaN encountered in AnalyticLocResid1 (ec)!!'
        print '(2(A,es12.5))', '   P =',P2,', Q =',Q2
        print '(2(A,es12.5))', '   dt =',dt,', ecsrc =',src_ec
        print '(2(A,es12.5))', '   G =',G,', Lambda =',Lambda
        print '(2(A,es12.5))', '   ec0 =',ec0,', ecUnits =',ecUnits
        ier = 0
        return
     end if

     !    this statement checks if HIanal = NaN
     if (HIanal /= HIanal) then
        print *,'NaN encountered in AnalyticLocResid1 (HI)!!'
        print '(2(A,es12.5))', '   A =',aHI,', B =',bHI
        print '(2(A,es12.5))', '   C =',cHI,', D =',dHI
        if (dHI/bHI/bHI < -1.e-8_RKIND) then
           print '(1(A,es12.5))', '   sqD =',sqrt(-dHI)
        elseif (dHI/bHI/bHI > 1.e-8_RKIND) then
           print '(3(A,es12.5))', '   rt1 =',(-bHI+sqrt(dHI))*0.5_RKIND/aHI, &
                ', rt2 =',(-bHI-sqrt(dHI))*0.5_RKIND/aHI,', expVHI =',expVHI
        else
           print '(1(A,es12.5))', '   rt1 =',-bHI*0.5_RKIND/aHI
        end if
        print '(2(A,es12.5))', '   K =',kHI,', G_HI =',G_HI
        print '(2(A,es12.5))', '  HIval =',HIval,', HI0 =',HI0
        print '(2(A,es12.5))', '  nH =',nH,', dt =',dt
        print '(2(A,es12.5))', '  NiUnits =',NiUnits,',  T =',T
        print '(2(A,es12.5))', '  k1 =',k1,', k2 =',k2
        ier = 0
        return
     end if

     !    this statement checks if HeIanal = NaN
     if (HeIanal /= HeIanal) then
        print *,'NaN encountered in AnalyticLocResid1 (HeI)!!'
        print '(2(A,es12.5))', '   A =',aHeI,', B =',bHeI
        print '(2(A,es12.5))', '   C =',cHeI,', D =',dHeI
        if (dHeI/bHeI/bHeI < -1.e-8_RKIND) then
           print '(1(A,es12.5))', '   sqD =',sqrt(-dHeI)
        elseif (dHeI/bHeI/bHeI > 1.e-8_RKIND) then
           print '(3(A,es12.5))', '   rt1 =',(-bHeI+sqrt(dHeI))*0.5_RKIND/aHeI, &
                ', rt2 =',(-bHeI-sqrt(dHeI))*0.5_RKIND/aHeI,', expVHeI =',expVHeI
        else
           print '(1(A,es12.5))', '   rt1 =',-bHeI*0.5_RKIND/aHeI
        end if
        print '(2(A,es12.5))', '   K =',kHeI,', G_HeI =',G_HeI
        print '(2(A,es12.5))', '  HeIval =',HeIval,', HeI0 =',HeI0
        print '(2(A,es12.5))', '  HeIIval =',HeIIval,', HeII0 =',HeII0
        print '(2(A,es12.5))', '  nHe =',nHe,', dt =',dt
        print '(2(A,es12.5))', '  NiUnits =',NiUnits,',  T =',T
        print '(2(A,es12.5))', '  k3 =',k3,', k4 =',k4
        ier = 0
        return
     end if

     !    this statement checks if HeIIanal = NaN
     if (HeIIanal /= HeIIanal) then
        print *,'NaN encountered in AnalyticLocResid1 (HeII)!!'
        print '(2(A,es12.5))', '   A =',aHeII,', B =',bHeII
        print '(2(A,es12.5))', '   C =',cHeII,', D =',dHeII
        if (dHeII/bHeII/bHeII < -1.e-8_RKIND) then
           print '(1(A,es12.5))', '   sqD =',sqrt(-dHeII)
        elseif (dHeII/bHeII/bHeII > 1.e-8_RKIND) then
           print '(3(A,es12.5))', '   rt1 =',(-bHeII+sqrt(dHeII))*0.5_RKIND/aHeII, &
                ', rt2 =',(-bHeII-sqrt(dHeII))*0.5_RKIND/aHeII,', expVHeII =',expVHeII
        else
           print '(1(A,es12.5))', '   rt1 =',-bHeII*0.5_RKIND/aHeII
        end if
        print '(2(A,es12.5))', '   K =',kHeII,', G_HeII =',G_HeII
        print '(2(A,es12.5))', '  HeIval =',HeIval,', HeI0 =',HeI0
        print '(2(A,es12.5))', '  HeIIval =',HeIIval,', HeII0 =',HeII0
        print '(2(A,es12.5))', '  nHe =',nHe,', dt =',dt
        print '(2(A,es12.5))', '  NiUnits =',NiUnits,',  T =',T
        print '(3(A,es12.5))', '  k4 =',k4,', k5 =',k5,', k6 =',k6
        ier = 0
        return
     end if

  end if  ! Nchem

  ! exit subroutine
  return

end subroutine gFLDProblem_AnalyticLocResid1
!=======================================================================





subroutine gFLDProblem_AnalyticLocResid4(Erres, HIres, Er, HI, Er0, HI0, &
     dt, rho, src_Er, src_HI, HFrac, ESpectrum, a, adot, IsE, IsEsHI,    &
     IsEsHInu, k1, k2, DenUnits, ErUnits, NiUnits, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       July 2009
  !
  !  PURPOSE: analytical chemistry/radiation residuals (model 4, nchem 1)
  !
  !  INPUTS:
  !     Er         - radiation energy array
  !     HI         - Hydrogen I number density array
  !     Er0        - old radiation energy array
  !     HI0        - old Hydrogen I number density array
  !     dt         - time step size
  !     rho        - density array
  !     src_Er     - source function values for radiation equation
  !     src_HI     - source function values for HI eq.
  !     HFrac      - percentage of mass composed of Hydrogen
  !     ESpectrum  - flag denoting choice of radiation energy spectrum
  !     a          - cosmological expansion parameter
  !     adot       - da/dt
  !     IsE        - int_{nu0}^{inf} sigE dnu
  !     IsEsHI     - int_{nu0}^{inf} sigE*sigHI dnu
  !     IsEsHInu   - int_{nu0}^{inf} sigE*sigHI/nu dnu
  !     k1,k2      - Hydrogen chemistry rates
  !     *Units     - variable scaling constants
  !
  !  OUTPUT ARGUMENTS: 
  !     Erres   - radiation energy density residual
  !     HIres   - Hydrogen I number density residual
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
  INTG_PREC, intent(in)  :: ESpectrum
  INTG_PREC, intent(out) :: ier
  REAL*8, intent(in)   :: a, adot, dt, HFrac, DenUnits, ErUnits, NiUnits
  REAL*8, intent(in)   :: IsE, IsEsHI, IsEsHInu, k1, k2
  REAL*8, intent(in)   :: Er, HI, Er0, HI0, rho, src_Er, src_HI
  REAL*8, intent(out)  :: Erres, HIres

  !--------------
  ! locals
  REAL*8 :: afac, c, hp, mp, nu0_HI, nu0_HeI, nu0_HeII, HIconst
  REAL*8 :: min_ni, min_rad, HIval, Erval, rhoval, nH, nHI, nHII, ne
  REAL*8 :: HIanal, Eranal, grey, kE, P1, Q1, G_HI
  REAL*8 :: aHI, bHI, cHI, dHI, kHI, expArg, expVHI, sqD, rt1, rt2

  !=======================================================================

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  HIres = 0._RKIND

  ! initialize constants
  afac = adot/a
  c  = c_light                 ! speed of light [cm/s]
  hp = hplanck                 ! Planck's constant [ergs*s]
  mp = mass_h                  ! Mass of a proton [g]
  nu0_HI = 13.6_RKIND*ev2erg/hp    ! ionization frequency of HI   [hz]
  min_ni   = 0._RKIND              ! minimum chem number density [cm^{-3}]
  min_rad  = 0._RKIND              ! minimum radiation density [ergs/cm^3]
  HIconst  = c*(IsEsHI - nu0_HI*IsEsHInu)/IsE
  grey = 1._RKIND                  ! grey vs monochromatic coeff for eqns
  if (ESpectrum == -1)  grey = 0._RKIND


  ! get shortcut values time-centered variables
  Erval = (Er0 + Er)*0.5d0
  HIval = (HI0 + HI)*0.5d0
  rhoval = rho*DenUnits
  nH = Hfrac*rhoval/mp
  nHI = HIval*NiUnits
  nHII = max(nH - nHI, 0.d0)
  ne = nHII

  ! compute radiation ODE rates
  kE = HIval*NiUnits*IsEsHI/IsE
  P1 = grey*afac + c*kE
  Q1 = src_Er/ErUnits

  ! compute HI chemistry ODE rates
  G_HI = c*Erval*ErUnits/hp*IsEsHInu/IsE
  aHI = (k1+k2)*NiUnits
  bHI = -(k1 + 2._RKIND*k2)*nH - G_HI
  cHI = (k2*nH*nH + src_HI)/NiUnits


  ! compute quasi-steady-state solution for Er, place in Eranal
  if (abs(P1) < 1.e-14_RKIND) then
     Eranal = Er0 + Q1*dt
  else
     if (P1*dt > 7.e2_RKIND) then
        Eranal = Q1/P1
     else
        Eranal = (Er0 - Q1/P1)*exp(-P1*dt) + Q1/P1
     end if
  end if


  ! compute quasi-steady-state solution for HI, place in HIres
  dHI = bHI*bHI - 4._RKIND*aHI*cHI

  !    if the second-order terms are very small, treat it as a linear ODE
  if (abs(aHI*HIval**2/(cHI-bHI*HIval)) < 1.0d-8) then
     if (bHI*dt < -7.e2_RKIND) then
        HIanal = -cHI/bHI
     else
        HIanal = (HI0 + cHI/bHI)*exp(bHI*dt) - cHI/bHI
     end if

     !    otherwise, go ahead with the quadratic ODE solution
  else

     ! no roots to quadratic
     if (dHI/bHI/bHI < -1.e-8_RKIND) then
        sqD = sqrt(-dHI)
        kHI = 2._RKIND/sqD*atan((2._RKIND*aHI*HI0+bHI)/sqD)
        HIanal = (sqD*tan((dt+kHI)*0.5_RKIND*sqD)-bHI)*0.5_RKIND/aHI

        ! two roots to quadratic
     elseif (dHI/bHI/bHI > 1.e-8_RKIND) then
        rt1 = (-bHI+sqrt(dHI))*0.5_RKIND/aHI
        rt2 = (-bHI-sqrt(dHI))*0.5_RKIND/aHI 
        kHI = log(abs((HI0-rt2)/(HI0-rt1)))/aHI/(rt2-rt1)
        expArg = aHI*(dt+kHI)*(rt2-rt1)
        if (expArg < -7.e2_RKIND) then
           HIanal = rt2
        elseif (expArg > 7.e2_RKIND) then
           HIanal = rt1
        else
           expVHI = exp(expArg)
           if ((HI0-rt1)*(HI0-rt2) > 0._RKIND) then
              HIanal = (rt2-rt1*expVHI)/(1._RKIND-expVHI)
           else
              HIanal = (rt2+rt1*expVHI)/(1._RKIND+expVHI)
           end if
        end if

        ! one double root to quadratic
     else
        rt1 = -bHI*0.5_RKIND/aHI
        kHI = 1._RKIND/aHI/(rt1-HI0)
        HIanal = rt1 - 1._RKIND/aHI/(dt+kHI)

     end if  ! roots
  end if  ! quadratic


  ! enforce bounds
  Eranal = max(Eranal, min_rad/ErUnits)
  HIanal = min(HIanal, nH/NiUnits)
  HIanal = max(HIanal, min_ni/NiUnits)


  ! compute residuals
  Erres = Er - Eranal
  HIres = HI - HIanal

  ! check some things
  !    this statement checks if Eranal = NaN
  if (Eranal /= Eranal) then
     print *,'NaN encountered in AnalyticLocResid4 (Er)!!'
     print '(2(A,es12.5))', '   P =',P1,', Q =',Q1
     print '(2(A,es12.5))', '   dt =',dt,', Ersrc =',src_Er
     print '(2(A,es12.5))', '   grey =',grey,', kE =',kE
     print '(2(A,es12.5))', '   Er0 =',Er0,', ErUnits =',ErUnits
     ier = 0
     return
  end if

  !    this statement checks if HIanal = NaN
  if (HIanal /= HIanal) then
     print *,'NaN encountered in AnalyticLocResid4 (HI)!!'
     print '(2(A,es12.5))', '   A =',aHI,', B =',bHI
     print '(2(A,es12.5))', '   C =',cHI,', D =',dHI
     if (dHI/bHI/bHI < -1.e-8_RKIND) then
        print '(1(A,es12.5))', '   sqD =',sqrt(-dHI)
     elseif (dHI/bHI/bHI > 1.e-8_RKIND) then
        print '(3(A,es12.5))', '   rt1 =',(-bHI+sqrt(dHI))*0.5_RKIND/aHI, &
             ', rt2 =',(-bHI-sqrt(dHI))*0.5_RKIND/aHI,', expVHI =',expVHI
     else
        print '(1(A,es12.5))', '   rt1 =',-bHI*0.5_RKIND/aHI
     end if
     print '(2(A,es12.5))', '   K =',kHI,', G_HI =',G_HI
     print '(2(A,es12.5))', '  HIval =',HIval,', HI0 =',HI0
     print '(2(A,es12.5))', '  nH =',nH,', dt =',dt
     print '(1(A,es12.5))', '  NiUnits =',NiUnits
     print '(2(A,es12.5))', '  k1 =',k1,', k2 =',k2
     ier = 0
     return
  end if

  ! exit subroutine
  return

end subroutine gFLDProblem_AnalyticLocResid4
!=======================================================================

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
subroutine gFLDProblem_LocalRHS(rhs_Er, rhs_ec, rhs_HI, rhs_HeI,       &
     rhs_HeII, src_Er, src_ec, src_HI, src_HeI, src_HeII, time, vx,    &
     vy, vz, rhoa, eca, Era, n_HIa, n_HeIa, n_HeIIa, Tempa, eha,       &
     kappaPa, kappaEa, a, adot, gamma, HFrac, Model, AnalyticChem,     &
     ESpectrum, CompA, Comp_xray, Comp_temp, IsE, IsEsHI, IsEsHInu,    &
     IsEsHeI, IsEsHeInu, IsEsHeII, IsEsHeIInu, NTempBins, TempStart,   &
     TempEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb,     &
     ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb, reHIITb,           &
     reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, piHI, piHeI, piHeII,     &
     aUnits, DenUnits, VelUnits, LenUnits, ErUnits, ecUnits, NiUnits,  &
     ecScale, Nchem, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,   &
     NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       August, 2006
!
!  PURPOSE: Computes the spatially-local components of the rad-hydro 
!           system of equations.
!
!  INPUTS:
!     src_Er     - source function values for radiation equation
!     src_ec     - source function values for gas energy correction eq.
!     src_HI     - source function values for HI eq.
!     src_HeI    - source function values for HeI eq.
!     src_HeII   - source function values for HeII eq.
!     time       - simulation time for evaluation
!     vx,vy,vz   - velocity arrays in each direction
!     rhoa       - density array
!     eca        - fluid energy correction array
!     Era        - radiation energy array
!     n_HIa      - Hydrogen I number density array
!     n_HeIa     - Helium I number density array
!     n_HeIIa    - Helium II number density array
!     Tempa      - Fluid temperature array
!     eha        - total fluid energy array
!     kappaPa    - Planck mean opacity array
!     kappaEa    - Energy mean opacity array
!     a          - cosmological expansion parameter
!     adot       - da/dt
!     gamma      - constant in ideal gas law
!     HFrac      - percentage of mass composed of Hydrogen
!     ESpectrum  - radiation spectrum choice
!                       1 -> 1e5 black body spectrum
!                       0 -> power law spectrum
!                      -1 -> monochromatic 
!     Model      - flag denoting physical model to use
!     CompA      - Compton cooling coefficient 1 (multiplier)
!     Comp_xray  - X-ray Compton heating coefficient
!     Comp_temp  - X-ray Compton heating temperature 
!     IsE        - int_{nu0}^{inf} sigE dnu
!     IsEsHI     - int_{nu0}^{inf} sigE*sigHI dnu
!     IsEsHInu   - int_{nu0}^{inf} sigE*sigHI/nu dnu
!     IsEsHeI    - int_{nu0}^{inf} sigE*sigHeI dnu
!     IsEsHeInu  - int_{nu0}^{inf} sigE*sigHeI/nu dnu
!     IsEsHeII   - int_{nu0}^{inf} sigE*sigHeII dnu
!     IsEsHeIInu - int_{nu0}^{inf} sigE*sigHeII/nu dnu
!     *Units     - variable scaling constants
!     Nchem      - number of chemical species in simulation
!     dx,dy,dz   - mesh spacing (comoving) in each direction
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     rhs_Er      - local nonlinear rhs for the Gray radiation energy
!                  density equation
!     rhs_ec      - local nonlinear rhs for the fluid energy 
!                  correction equation
!     rhs_HI     - local nonlinear rhs for the HI chemistry equation
!     rhs_HeI    - local nonlinear rhs for the HeI chemistry eq.
!     rhs_HeII   - local nonlinear rhs for the HeII chemistry eq.
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
  INTG_PREC, intent(in)  :: Model, Nchem, ESpectrum, NTempBins, AnalyticChem
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a, adot
  R_PREC,    intent(in)  :: time, dx, dy, dz, gamma, HFrac
  R_PREC,    intent(in)  :: CompA, Comp_xray, Comp_temp
  R_PREC,    intent(in)  :: TempStart, TempEnd, piHI, piHeI, piHeII
  R_PREC,    intent(in)  :: aUnits, DenUnits, VelUnits, LenUnits, ErUnits, &
       ecUnits, NiUnits, ecScale
  R_PREC,    intent(in)  :: IsE, IsEsHI, IsEsHInu, IsEsHeI
  R_PREC,    intent(in)  :: IsEsHeInu, IsEsHeII, IsEsHeIInu
  R_PREC,    intent(in),                                              &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) ::   &
       src_Er, src_ec, src_HI, src_HeI, src_HeII, vx, vy, vz, rhoa, &
       eca, Era, n_HIa, n_HeIa, n_HeIIa, Tempa, eha, kappaEa, kappaPa
  R_PREC,    intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb, k4Tb,   &
       k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, &
       ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb
  R_PREC,    intent(out),                                           &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       rhs_Er, rhs_ec, rhs_HI, rhs_HeI, rhs_HeII

!--------------
! locals
  INTG_PREC :: i, j, k, Tidx, Tidxp
  R_PREC :: lTempS, lTempE, dlTemp, lTemp, Tl, Tr, Tfac
  R_PREC :: k1, k2, k3, k4, k5, k6
  R_PREC :: aval, afac, grey, c, hp, mp, zr, gam_1, StBz, pi
  R_PREC :: HIconst, HeIconst, HeIIconst, kappaP, kappaE
  R_PREC :: dxi2, dyi2, dzi2, DivV, GradRhoDotV
  R_PREC :: rho, ec, eh, Er, nH, nHI, nHII, nHe, nHeI, nHeII, nHeIII, ne
  R_PREC :: T, lamT, G, Lambda
  R_PREC :: alpha, beta, eta, nu0_HI, nu0_HeI, nu0_HeII
  R_PREC :: ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII
  R_PREC :: reHII, reHeII1, reHeII2, reHeIII, brem, Comp1, Comp2
  R_PREC :: G_HI, G_HeI, G_HeII

!=======================================================================

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  rhs_Er = 0._RKIND
  rhs_ec = 0._RKIND
  if (Nchem > 0) then
     rhs_HI = 0._RKIND
  endif
  if (Nchem == 3) then
     rhs_HeI  = 0._RKIND
     rhs_HeII = 0._RKIND
  endif
  
  ! check that chemistry constants make sense
  if ((Nchem /= 0) .and. (Nchem /= 1) .and. (Nchem /= 3)) then
     write(*,*) 'Chemistry ERROR: illegal value, Nchem = ',Nchem, &
          ',  Nchem must be one of {0, 1, 3}.  Returning!'
     ier = 0
     return
  endif
!!$  if ((Nchem == 1) .and. (HFrac /= 1._RKIND)) then
!!$     write(*,*) 'Chemistry ERROR: illegal value, HFrac = ',HFrac, &
!!$          ',  value must equal 1 for Hydrogen only case.  Returning!'
!!$     ier = 0
!!$     return     
!!$  endif
  if ((HFrac < 0._RKIND) .or. (HFrac > 1._RKIND)) then
     write(*,*) 'Chemistry ERROR: illegal value, HFrac = ',HFrac, &
          ',  value must be in the interval [0,1].  Returning!'
     ier = 0
     return     
  endif
  if ((AnalyticChem == 1) .and. (Nchem /= 1)) then
     write(*,*) 'Chemistry ERROR: AnaltyicChem requires Nchem = 1 (for now)'
     ier = 0
     return
  end if
  if ((AnalyticChem == 1) .and. (Model /= 1) .and. (Model /= 4)) then
     write(*,*) 'Chemistry ERROR: AnaltyicChem requires Model = 1 or 4 (for now)'
     ier = 0
     return
  end if


  ! initialize constants
  aval = a*aunits
  gam_1 = gamma-1._RKIND
  pi = pi_val
  c  = c_light                 ! speed of light [cm/s]
  hp = hplanck                 ! Planck's constant [ergs*s]
  mp = mass_h                  ! Mass of a proton [g]
  zr = 1._RKIND/(aval) - 1._RKIND      ! cosmological redshift
  StBz  = 5.6704e-5_RKIND      ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  afac  = adot/a               ! adot/a
  alpha = -0.4910454_RKIND         ! exponent in emissivity fitting
  beta  = 2.17748887e-24_RKIND     ! scaling in emissivity fitting
  nu0_HI = 13.6_RKIND*ev2erg/hp    ! ionization frequency of HI   [hz]
  nu0_HeI = 24.6_RKIND*ev2erg/hp   ! ionization frequency of HeI  [hz]
  nu0_HeII = 54.4_RKIND*ev2erg/hp  ! ionization frequency of HeII [hz]
  grey = 1._RKIND                  ! grey vs monochromatic coeff for eqns
  if (ESpectrum == -1)  grey = 0._RKIND

  !   lookup table constants
  lTempS = log(TempStart)
  lTempE = log(TempEnd)
  dlTemp = (lTempE - lTempS)/(1._RKIND*NTempBins - 1._RKIND)

  ! compute shortcuts
  dxi2 = 0.5_RKIND*a/(dx*LenUnits)  ! convert to proper units during division
  dyi2 = 0.5_RKIND*a/(dy*LenUnits)
  dzi2 = 0.5_RKIND*a/(dz*LenUnits)
  HIconst   = c*(IsEsHI   - nu0_HI*IsEsHInu)/IsE
  HeIconst  = c*(IsEsHeI  - nu0_HeI*IsEsHeInu)/IsE
  HeIIconst = c*(IsEsHeII - nu0_HeII*IsEsHeIInu)/IsE
  if (aval .ne. 1._RKIND) then        ! Compton cooling coefficients
     Comp1 = CompA*(1._RKIND + zr)**4
     Comp2 = 2.73_RKIND*(1._RKIND + zr)
  else
     Comp1 = 0._RKIND
     Comp2 = 0._RKIND
  endif


  ! compute right-hand sides depending on physical model

  !  First compute gas energy correction-based adjustments alone
  !  (as these involve derivatives, use loop for appropriate dimension)
!!$  ! 1D model
!!$  if (Ny == 1) then
!!$     do k=1,Nz,1
!!$        do j=1,Ny,1
!!$           do i=1,Nx,1
!!$
!!$              ! shortcuts
!!$              ec = eca(i,j,k) * ecUnits
!!$              
!!$              ! velocity divergence
!!$              DivV = dxi2*(vx(i+1,j,k)-vx(i-1,j,k))*VelUnits
!!$              
!!$              ! (grad density).dot.(velocity)/density
!!$              GradRhoDotV = dxi2*(rhoa(i+1,j,k)-rhoa(i-1,j,k)) &
!!$                   * vx(i,j,k) / rhoa(i,j,k) * VelUnits
!!$              
!!$              ! put it together
!!$              rhs_ec(i,j,k) = ec/aval*(DivV - gam_1*GradRhoDotV)
!!$              
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$  ! 2D model
!!$  else if (Nz == 1) then
!!$     do k=1,Nz,1
!!$        do j=1,Ny,1
!!$           do i=1,Nx,1
!!$              
!!$              ! shortcuts
!!$              ec = eca(i,j,k) * ecUnits
!!$              
!!$              ! velocity divergence
!!$              DivV = (dxi2*(vx(i+1,j,k)-vx(i-1,j,k))  &
!!$                    + dyi2*(vy(i,j+1,k)-vy(i,j-1,k))) * VelUnits
!!$              
!!$              ! (grad density).dot.(velocity)/density
!!$              GradRhoDotV = (dxi2*(rhoa(i+1,j,k)-rhoa(i-1,j,k))*vx(i,j,k)  &
!!$                           + dyi2*(rhoa(i,j+1,k)-rhoa(i,j-1,k))*vy(i,j,k)) &
!!$                           / rhoa(i,j,k)*VelUnits
!!$              
!!$              ! put it together
!!$              rhs_ec(i,j,k) = ec/aval*(DivV - gam_1*GradRhoDotV)
!!$              
!!$           enddo
!!$        enddo
!!$     enddo
!!$     
!!$  ! 3D model
!!$  else
!!$     do k=1,Nz,1
!!$        do j=1,Ny,1
!!$           do i=1,Nx,1
!!$              
!!$              ! shortcuts
!!$              ec = eca(i,j,k) * ecUnits
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
!!$              rhs_ec(i,j,k) = ec/aval*(DivV - gam_1*GradRhoDotV)
!!$              
!!$           enddo
!!$        enddo
!!$     enddo
!!$  endif



  !==================================================
  !  Now incorporate Model-specific equations

  !    decoupled ODE test case
  if (Model == 0) then

     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1
              
              ! set shortcut values for this spatial location
              ec = eca(i,j,k)*ecUnits
              eh = eha(i,j,k)*VelUnits*VelUnits

              ! put ODE rhs together
              rhs_ec(i,j,k) = (eh + ec + 2._RKIND*time*exp(2._RKIND*time))/ecUnits
              rhs_Er(i,j,k) = 10._RKIND*cos(20._RKIND*time)/ErUnits
              rhs_HI(i,j,k) = (-sin(2._RKIND*time) + 8._RKIND*time)/NiUnits


           enddo
        enddo
     enddo


  !==================================================
  !    coupled ODE test case
  else if (Model == 3) then

     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1
              
              ! set shortcut values for this spatial location
              ec = eca(i,j,k)*ecUnits
              eh = eha(i,j,k)*VelUnits*VelUnits
              Er = Era(i,j,k)*ErUnits
              nHI = n_HIa(i,j,k)*NiUnits

              ! put ODE rhs together
              rhs_ec(i,j,k) = 10._RKIND*((7._RKIND/3._RKIND)*(eh+ec)   &
                                 + (-2._RKIND/3._RKIND)*nHI       &
                                 + (-2._RKIND/3._RKIND)*Er)/ecUnits
              rhs_Er(i,j,k) = 10._RKIND*((-2._RKIND/3._RKIND)*(eh+ec)  &
                                 + (-1._RKIND/6._RKIND)*nHI       &
                                 + (11._RKIND/6._RKIND)*Er)/ErUnits
              rhs_HI(i,j,k) = 10._RKIND*((-2._RKIND/3._RKIND)*(eh+ec)  &
                                 + (11._RKIND/6._RKIND)*nHI       &
                                 - (1._RKIND/6._RKIND)*Er)/NiUnits

           enddo
        enddo
     enddo


  !==================================================
  !    isothermal case-B Hydrogen ionization test case
  else if (Model == 4) then

     ! replace theta method with analytical solver for ec, HI eqns
     if (AnalyticChem == 1) then
        
        ! iterate over the domain
        do k=1,Nz,1
           do j=1,Ny,1
              do i=1,Nx,1
                 
                 ! set shortcut values for this spatial location,
                 ! converting densities from comoving to proper, and 
                 ! put all shortcuts in ACTUAL units
                 Er  = Era(i,j,k)*ErUnits
                 nHI = n_HIa(i,j,k)*NiUnits
                 
                 ! opacity
                 kappaE = nHI*IsEsHI/IsE
              
                 ! put it all together
                 rhs_Er(i,j,k) = (src_Er(i,j,k)          &
                      - (afac*grey + c*kappaE)*Er)/ErUnits

              enddo
           enddo
        enddo

     else

     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1
              
              ! set shortcut values for this spatial location,
              ! converting densities from comoving to proper, and 
              ! put all shortcuts in CGS units
              Er = Era(i,j,k)*ErUnits
              rho = rhoa(i,j,k)*DenUnits
              nHI = n_HIa(i,j,k)*NiUnits
              nH = Hfrac*rho/mp
              nHII = max(nH - nHI, 0._RKIND)
              ne = nHII

              ! opacity
              kappaE = nHI*IsEsHI/IsE
              
              ! temperature shortcuts
              !*** For Model 4 with cosmology, the temperature is held in ecScale ***!
              if (adot == 0._RKIND) then
                 T = Tempa(i,j,k)
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

              ! compute case B Hydrogen recombination coefficient 
              ! [Hui & Gnedin, 1997: RI^B_{HII}]
              !   (still need this because table holds case A coefficient!)
              k2 = 2.753e-14_RKIND*lamT**(1.5_RKIND) *                 &
                   (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-2.242_RKIND)

              ! compute Hydrogen photoionization rate
              G_HI = c*Er/hp*IsEsHInu/IsE
              
              ! compute (comoving, scaled) rhs for HI species
              rhs_HI(i,j,k) = (src_HI(i,j,k) + k2*ne*nHII &
                   - nHI*(k1*ne + G_HI))/NiUnits
              
              ! compute (comoving, scaled) local radiation rhs
              rhs_Er(i,j,k) = (src_Er(i,j,k)          &
                   - (afac*grey + c*kappaE)*Er)/ErUnits

              ! decouple gas energy from all other variables
              rhs_ec(i,j,k) = 0._RKIND
              
           enddo
        enddo
     enddo

  endif  ! AnalyticChem

  !==================================================
  !    isothermal point source emissivity test case (decoupled chem, gas)
  else if (Model == 5) then

     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1
              
              ! set shortcuts
              Er = Era(i,j,k)*ErUnits
              
              ! opacity
              kappaE = nHI*IsEsHI/IsE
              
              ! compute local radiation right-hand side
              rhs_Er(i,j,k) = (src_Er(i,j,k)          &
                   - (afac*grey + c*kappaE)*Er)/ErUnits

              ! decouple gas energy from all other variables
              rhs_ec(i,j,k) = 0._RKIND
              
              ! decouple HI evolution from all other variables
              rhs_HI(i,j,k) = 0._RKIND
              
           enddo
        enddo
     enddo


  !==================================================
  !    Zeus test cases (no chemistry)
  else if ((Model == 10) .or. ((Model >= 20) .and. (Model <= 29))) then

     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1

              ! set shortcut values for this spatial location
              ! converting densities from comoving to proper
              rho = rhoa(i,j,k)*DenUnits
              Er = Era(i,j,k)*ErUnits
              ec = eca(i,j,k)*ecUnits
              kappaP = kappaPa(i,j,k)
              kappaE = kappaEa(i,j,k)
              T = Tempa(i,j,k)
              eta = kappaP*StBz/pi*T**4

              ! put it all together
              rhs_ec(i,j,k) = (src_ec(i,j,k) + rhs_ec(i,j,k)  &
                   - 2._RKIND*afac*ec + (c*kappaE*Er - 4._RKIND*pi*eta)/rho)/ecUnits
              rhs_Er(i,j,k) = (src_Er(i,j,k) + 4._RKIND*pi*eta    &
                   - (afac*grey + c*kappaE)*Er)/ErUnits

           enddo
        enddo
     enddo


  !==================================================
  !     case A HII recombination rate, with emissivity
  else if (Model == 2) then

     ! Hydrogen only case
     if (Nchem == 1) then
        ! iterate over the domain
        do k=1,Nz,1
           do j=1,Ny,1
              do i=1,Nx,1
                 
                 ! set shortcut values for this spatial location,
                 ! converting densities from comoving to proper, and 
                 ! put all shortcuts in ACTUAL units
                 ec   = eca(i,j,k)*ecUnits
                 Er  = Era(i,j,k)*ErUnits
                 rho = rhoa(i,j,k)*DenUnits
                 nHI = n_HIa(i,j,k)*NiUnits
                 nH = Hfrac*rho/mp
                 nHII = max(nH - nHI, 0._RKIND)
                 ne  = nHII
                 
                 ! opacity
                 kappaE = nHI*IsEsHI/IsE
              
                 ! shortcuts for temperature, ln(temp) and powers
                 T = Tempa(i,j,k)
                 
                 ! look up rates
                 lTemp = min(max(log(T), lTempS), lTempE)
                 Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
                 Tidxp = Tidx+1
                 Tl = lTempS + (Tidx-1)*dlTemp
                 Tr = lTempS +  Tidx*dlTemp
                 Tfac = (lTemp - Tl)/(Tr - Tl)
                 k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
                 k2 = k2Tb(Tidx) + (k2Tb(Tidxp) - k2Tb(Tidx))*Tfac
                 ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
                 ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
                 reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
                 brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac
                 
                 ! compute Hydrogen photoionization rate 
                 G_HI = c*Er/hp*IsEsHInu/IsE
                 
                 ! compute fluid cooling rate.  Terms (in order):
                 !    Collisional Excitations
                 !    Collisional Ionizations
                 !    Recombinations
                 !    Compton cooling or heating 
                 !    X-ray Compton heating
                 !    Bremsstrahlung
                 Lambda = ne/rho                 &
                      *(ceHI*nHI                 &
                      + ciHI*nHI                 &
                      + reHII*nHII               &
                      + Comp1*(T-Comp2)          &
                      + Comp_xray*(T-Comp_temp)  &
                      + brem*nHII                &
                      )
                 
                 ! compute fluid heating rate
                 G = nHI/rho*Er*HIconst
                 
                 ! compute emissivity using power law approximation
                 eta = nHII*ne*beta*T**alpha

                 ! put it all together
                 rhs_ec(i,j,k) = (src_ec(i,j,k) + rhs_ec(i,j,k) &
                      + G - Lambda - 2._RKIND*afac*ec)/ecUnits
                 rhs_Er(i,j,k) = (src_Er(i,j,k) + 4._RKIND*pi*eta   &
                      - (afac*grey + c*kappaE)*Er)/ErUnits
                 rhs_HI(i,j,k) = (src_HI(i,j,k) + k2*ne*nHII    &
                      - nHI*(k1*ne + G_HI))/NiUnits
                 
              enddo
           enddo
        enddo

     ! Hydrogen + Helium case        
     else if (Nchem == 3) then
        ! iterate over the domain
        do k=1,Nz,1
           do j=1,Ny,1
              do i=1,Nx,1
                 
                 ! set shortcut values for this spatial location,
                 ! converting densities from comoving to proper, and 
                 ! put all shortcuts in ACTUAL units
                 ec     = eca(i,j,k)*ecUnits
                 Er     = Era(i,j,k)*ErUnits
                 rho    = rhoa(i,j,k)*DenUnits
                 nHI    = n_HIa(i,j,k)*NiUnits
                 nH     = Hfrac*rho/mp
                 nHII   = max(nH - nHI, 0._RKIND)
                 nHe    = (1._RKIND-HFrac)*rho/4._RKIND/mp
                 nHeI   = n_HeIa(i,j,k)*NiUnits/4._RKIND
                 nHeII  = n_HeIIa(i,j,k)*NiUnits/4._RKIND
                 nHeIII = max(nHe - nHeI - nHeII, 0._RKIND)
                 ne     = nHII + 0.25_RKIND*nHeII + 0.5_RKIND*nHeIII
                 
                 ! opacity
                 kappaE = nHI*IsEsHI/IsE + nHeI*IsEsHeI/IsE + nHeII*IsEsHeII/IsE
              
                 ! shortcuts for temperature, ln(temp) and powers
                 T = Tempa(i,j,k)
                 
                 ! look up rates
                 lTemp = min(max(log(T), lTempS), lTempE)
                 Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
                 Tidxp = Tidx+1
                 Tl = lTempS + (Tidx-1)*dlTemp
                 Tr = lTempS +  Tidx*dlTemp
                 Tfac = (lTemp - Tl)/(Tr - Tl)
                 k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
                 k2 = k2Tb(Tidx) + (k2Tb(Tidxp) - k2Tb(Tidx))*Tfac
                 k3 = k3Tb(Tidx) + (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac
                 k4 = k4Tb(Tidx) + (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac
                 k5 = k5Tb(Tidx) + (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac
                 k6 = k6Tb(Tidx) + (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac
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

                 ! compute Hydrogen photoionization rate 
                 G_HI = c*Er/hp*IsEsHInu/IsE
                 
                 ! compute Helium photoionization rates 
                 G_HeI  = c*Er/hp*IsEsHeInu/IsE
                 G_HeII = c*Er/hp*IsEsHeIInu/IsE

                 ! compute fluid cooling rate.  Terms (in order):
                 !    Collisional Excitations (3)
                 !    Collisional Ionizations (4)
                 !    Recombinations (4)
                 !    Compton cooling or heating (1)
                 !    X-ray Compton heating (1)
                 !    Bremsstrahlung (1)
                 Lambda = ne/rho                       &
                      *(ceHI*nHI                       &
                      + ceHeI*nHeII*ne/4._RKIND            &
                      + ceHeII*nHeII/4._RKIND              &
                      + ciHI*nHI                       &
                      + ciHeI*nHeI/4._RKIND                &
                      + ciHeII*nHeII/4._RKIND              &
                      + ciHeIS*nHeII*ne/4._RKIND           &
                      + reHII*nHII                     &
                      + reHeII1*nHeII/4._RKIND             &
                      + reHeII2*nHeII/4._RKIND             &
                      + reHeIII*nHeIII/4._RKIND            &
                      + Comp1*(T-Comp2)                &
                      + Comp_xray*(T-Comp_temp)        &
                      + brem*(nHII+nHeII/4._RKIND+nHeIII)  &
                      )
                 
                 ! compute fluid heating rate
                 G = Er/rho*(nHI*HIconst + nHeI*HeIconst + nHeII*HeIIconst)
                 
                 ! compute emissivity 
                 ! [WE HAVE NO FORMULA FOR ETA IN MULTI-SPECIES CASE!!]
                 eta = 0._RKIND

                 ! put it all together
                 rhs_ec(i,j,k) = (src_ec(i,j,k) + rhs_ec(i,j,k)      &
                      + G - Lambda - 2._RKIND*afac*ec)/ecUnits
                 rhs_Er(i,j,k) = (src_Er(i,j,k) + 4._RKIND*pi*eta        &
                      - (afac*grey + c*kappaE)*Er)/ErUnits
                 rhs_HI(i,j,k) = (src_HI(i,j,k) + k2*ne*nHII         &
                      - nHI*(k1*ne + G_HI))/NiUnits
                 rhs_HeI(i,j,k) = (src_HeI(i,j,k) + k4*ne*nHeII      &
                      - k3*ne*nHeI - nHeI*G_HeI)*4._RKIND/NiUnits
                 rhs_HeII(i,j,k) = (src_HeII(i,j,k) + nHeI*G_HeI     &
                      - nHeII*G_HeII + ne*k3*nHeI + ne*k6*nHeIII     &
                      - ne*(k4+k5)*nHeII)*4._RKIND/NiUnits

              enddo
           enddo
        enddo
     else
        write(0,*) 'gFLDProblem_LocRHS: illegal Nchem =',Nchem,' for model =',Model
     endif  ! Nchem

  !==================================================
  !     case B HII recombination rate
  else if (Model == 1) then

     ! Hydrogen only case
     if (Nchem == 1) then
        
        ! replace theta method with analytical solver for ec, HI eqns
        if (AnalyticChem == 1) then
        
        ! iterate over the domain
        do k=1,Nz,1
           do j=1,Ny,1
              do i=1,Nx,1
                 
                 ! set shortcut values for this spatial location,
                 ! converting densities from comoving to proper, and 
                 ! put all shortcuts in ACTUAL units
                 Er  = Era(i,j,k)*ErUnits
                 nHI = n_HIa(i,j,k)*NiUnits
                 
                 ! opacity
                 kappaE = nHI*IsEsHI/IsE
              
                 ! put it all together
                 rhs_Er(i,j,k) = (src_Er(i,j,k)          &
                      - (afac*grey + c*kappaE)*Er)/ErUnits

              enddo
           enddo
        enddo

        else

        ! iterate over the domain
        do k=1,Nz,1
           do j=1,Ny,1
              do i=1,Nx,1
                 
                 ! set shortcut values for this spatial location,
                 ! converting densities from comoving to proper, and 
                 ! put all shortcuts in ACTUAL units
                 ec   = eca(i,j,k)*ecUnits
                 Er  = Era(i,j,k)*ErUnits
                 rho = rhoa(i,j,k)*DenUnits
                 nHI = n_HIa(i,j,k)*NiUnits
                 nH = Hfrac*rho/mp
                 nHII = max(nH - nHI, 0._RKIND)
                 ne  = nHII
                 
                 ! opacity
                 kappaE = nHI*IsEsHI/IsE
              
                 ! shortcuts for temperature, ln(temp) and powers
                 T = Tempa(i,j,k)
                 lamT = 3.15614d5/T
                 
                 ! look up rates
                 lTemp = min(max(log(T), lTempS), lTempE)
                 Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
                 Tidxp = Tidx+1
                 Tl = lTempS + (Tidx-1)*dlTemp
                 Tr = lTempS +  Tidx*dlTemp
                 Tfac = (lTemp - Tl)/(Tr - Tl)
                 k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
                 ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
                 ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
                 reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
                 brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac

                 ! compute case B Hydrogen recombination coefficient 
                 ! [Hui & Gnedin, 1997: RI^B_{HII}]
                 !   (still need this because table holds case A coefficient!)
                 k2 = 2.753e-14_RKIND*lamT**(1.5_RKIND) *                 &
                      (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(-2.242_RKIND)
                 
                 ! compute Hydrogen photoionization rate 
                 G_HI = c*Er/hp*IsEsHInu/IsE

                 ! compute fluid cooling rate.  Terms (in order):
                 !    Collisional Excitations
                 !    Collisional Ionizations
                 !    Recombinations
                 !    Compton cooling or heating 
                 !    X-ray Compton heating
                 !    Bremsstrahlung
                 Lambda = ne/rho                  &
                      *(ceHI*nHI                  &
                      + ciHI*nHI                  &
                      + reHII*nHII                &
                      + Comp1*(T-Comp2)           &
                      + Comp_xray*(T-Comp_temp)   &
                      + brem*nHII                 &
                      )
                 
                 ! compute fluid heating rate
                 G = nHI/rho*Er*HIconst
                 
                 ! compute emissivity 
                 ! [WE HAVE NO FORMULA FOR ETA IN CASE-B CASE!!]
                 eta = 0._RKIND

                 ! put it all together
                 rhs_ec(i,j,k) = (src_ec(i,j,k) + rhs_ec(i,j,k) &
                      + G - Lambda - 2._RKIND*afac*ec)/ecUnits
                 rhs_Er(i,j,k) = (src_Er(i,j,k) + 4._RKIND*pi*eta   &
                      - (afac*grey + c*kappaE)*Er)/ErUnits
                 rhs_HI(i,j,k) = (src_HI(i,j,k) + k2*ne*nHII    &
                      - nHI*(k1*ne + G_HI))/NiUnits

              enddo
           enddo
        enddo
        
        endif  ! AnalyticChem

     ! Hydrogen + Helium case        
     else if (Nchem == 3) then

        ! replace theta method with analytical solver for ec, HI eqns
        if (AnalyticChem == 1) then
        
        ! iterate over the domain
        do k=1,Nz,1
           do j=1,Ny,1
              do i=1,Nx,1
                 
                 ! set shortcut values for this spatial location,
                 ! converting densities from comoving to proper, and 
                 ! put all shortcuts in ACTUAL units
                 Er     = Era(i,j,k)*ErUnits
                 rho    = rhoa(i,j,k)*DenUnits
                 nHI    = n_HIa(i,j,k)*NiUnits
                 nHeI   = n_HeIa(i,j,k)*NiUnits/4._RKIND
                 nHeII  = n_HeIIa(i,j,k)*NiUnits/4._RKIND
                 
                 ! opacity
                 kappaE = nHI*IsEsHI/IsE + nHeI*IsEsHeI/IsE + nHeII*IsEsHeII/IsE
              
                 ! put it all together
                 rhs_Er(i,j,k) = (src_Er(i,j,k)          &
                      - (afac*grey + c*kappaE)*Er)/ErUnits

              enddo
           enddo
        enddo

        else

        ! iterate over the domain
        do k=1,Nz,1
           do j=1,Ny,1
              do i=1,Nx,1
                 
                 ! set shortcut values for this spatial location,
                 ! converting densities from comoving to proper, and 
                 ! put all shortcuts in ACTUAL units
                 ec     = eca(i,j,k)*ecUnits
                 Er     = Era(i,j,k)*ErUnits
                 rho    = rhoa(i,j,k)*DenUnits
                 nHI    = n_HIa(i,j,k)*NiUnits
                 nH     = Hfrac*rho/mp
                 nHII   = max(nH - nHI, 0._RKIND)
                 nHe    = (1._RKIND-HFrac)*rho/4._RKIND/mp
                 nHeI   = n_HeIa(i,j,k)*NiUnits/4._RKIND
                 nHeII  = n_HeIIa(i,j,k)*NiUnits/4._RKIND
                 nHeIII = max(nHe - nHeI - nHeII, 0._RKIND)
                 ne     = nHII + 0.25_RKIND*nHeII + 0.5_RKIND*nHeIII
                 
                 ! opacity
                 kappaE = nHI*IsEsHI/IsE + nHeI*IsEsHeI/IsE + nHeII*IsEsHeII/IsE
              
                 ! shortcuts for temperature, ln(temp) and powers
                 T = Tempa(i,j,k)
                 lamT = 3.15614e5_RKIND/T
                 
                 ! look up rates
                 lTemp = min(max(log(T), lTempS), lTempE)
                 Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp,IKIND)+1))
                 Tidxp = Tidx+1
                 Tl = lTempS + (Tidx-1)*dlTemp
                 Tr = lTempS +  Tidx*dlTemp
                 Tfac = (lTemp - Tl)/(Tr - Tl)
                 k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
                 k3 = k3Tb(Tidx) + (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac
                 k4 = k4Tb(Tidx) + (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac
                 k5 = k5Tb(Tidx) + (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac
                 k6 = k6Tb(Tidx) + (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac
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

                 ! compute case B Hydrogen recombination coefficient 
                 ! [Hui & Gnedin, 1997: RI^B_{HII}]
                 !   (still need this because table holds case A coefficient!)
                 k2 = 2.753e-14_RKIND*lamT**(1.5_RKIND) /                &
                      (1._RKIND+(lamT/2.74_RKIND)**(0.407_RKIND))**(2.242_RKIND)
                 
                 ! compute Hydrogen photoionization rate 
                 G_HI = c*Er/hp*IsEsHInu/IsE
                 
                 ! compute Helium photoionization rates 
                 G_HeI  = c*Er/hp*IsEsHeInu/IsE
                 G_HeII = c*Er/hp*IsEsHeIInu/IsE

                 ! compute fluid cooling rate.  Terms (in order):
                 !    Collisional Excitations (3)
                 !    Collisional Ionizations (4)
                 !    Recombinations (4)
                 !    Compton cooling or heating (1)
                 !    X-ray Compton heating (1)
                 !    Bremsstrahlung (1)
                 Lambda = ne/rho                       &
                      *(ceHI*nHI                       &
                      + ceHeI*nHeII*ne/4._RKIND            &
                      + ceHeII*nHeII/4._RKIND              &
                      + ciHI*nHI                       &
                      + ciHeI*nHeI/4._RKIND                &
                      + ciHeII*nHeII/4._RKIND              &
                      + ciHeIS*nHeII*ne/4._RKIND           &
                      + reHII*nHII                     &
                      + reHeII1*nHeII/4._RKIND             &
                      + reHeII2*nHeII/4._RKIND             &
                      + reHeIII*nHeIII/4._RKIND            &
                      + Comp1*(T-Comp2)                &
                      + Comp_xray*(T-Comp_temp)        &
                      + brem*(nHII+nHeII/4._RKIND+nHeIII)  &
                      )
                 
                 ! compute fluid heating rate
                 G = Er/rho*(nHI*HIconst + nHeI*HeIconst + nHeII*HeIIconst)
                 
                 ! compute emissivity 
                 ! [WE HAVE NO FORMULA FOR ETA IN CASE-B CASE!!]
                 eta = 0._RKIND

                 ! put it all together
                 rhs_ec(i,j,k) = (src_ec(i,j,k) + rhs_ec(i,j,k)      &
                      + G - Lambda - 2._RKIND*afac*ec)/ecUnits
                 rhs_Er(i,j,k) = (src_Er(i,j,k) + 4._RKIND*pi*eta        &
                      - (afac*grey + c*kappaE)*Er)/ErUnits
                 rhs_HI(i,j,k) = (src_HI(i,j,k) + k2*ne*nHII         &
                      - nHI*(k1*ne + G_HI))/NiUnits
                 rhs_HeI(i,j,k) = (src_HeI(i,j,k) + k4*ne*nHeII      &
                      - k3*ne*nHeI - nHeI*G_HeI)/NiUnits
                 rhs_HeII(i,j,k) = (src_HeII(i,j,k) + nHeI*G_HeI     &
                      - nHeII*G_HeII + ne*k3*nHeI + ne*k6*nHeIII     &
                      - ne*(k4+k5)*nHeII)/NiUnits

              enddo
           enddo
        enddo

        endif  ! AnalyticChem

     else
        write(0,*) 'gFLDProblem_LocRHS: illegal Nchem =',Nchem,' for model =',Model
     endif  ! Nchem


  !==================================================
  else

     write(0,*) 'gFLDProblem_LocRHS: Model =',Model,' undefined!'

  endif ! Model


  return
end subroutine gFLDProblem_LocalRHS
!=======================================================================

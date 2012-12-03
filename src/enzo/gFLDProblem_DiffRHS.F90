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
subroutine gFLDProblem_DiffRHS_3D(rhs, EgCur, EgOld, Temp, kappaE,     &
     LimType, a, aUnits, LenUnits, EgUnits, dx, dy, dz, Nx, Ny, Nz,    &
     NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, xlface, xrface, ylface,       &
     yrface, zlface, zrface, Model, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2006
!  modified1:  July 20, 2007, by John Hayes; applying transport correction
!              to diffusion coefficients on external boundaries.
!  modified2:  August 10, 2007, by John Hayes; appended "_3D" to routine name.
!  modified3:  December 21, 2007, by John Hayes; corrected formula for Marshak BC
!
!  PURPOSE: Computes the spatially-diffusive component of the 
!           nonlinear rhs for the Gray FLD problem,
!              -1/a^2*Div(D(Eg)*Grad(Eg))
!           where D(Eg) is a nonlinear flux-limiter depending on Eg.  
!           We define the values
!              R_i = |Grad(Eg)_i|/sigT/omega/Eg,
!              omega = B/E
!              sigT = total extinction coefficient.
!           The '_i' subscript implies the gradient in the ith 
!           direction; these quantities are all required at cell faces, 
!           as that is the location of the divergence calculations.
!           With these components, we allow any of the following three 
!           forms of the limiter, 
!             [Levermore-Pomraning, 1981],
!                 D_i(Eg) = c/sigT/R_i/omega*[coth(R_i)-1/R_i],
!             [rational approx. to above, Levermore-Pomraning, 1981],
!                 D_i(Eg) = c/sigT/omega*(2+R_i)/(6+3*R_i+R_i**2),
!             [Larsen n=2 limiter],
!                 D_i(Eg) = c/sqrt((3*sigT)**2 + R**2)
!             [Zeus form of rational approx. to LP],
!                 D_i(Eg) = 2/pi*c/sigT/R_i*atan(R_i*pi/6),
!                 R_i = |Grad(Eg)_i|/sigT/Eg,
!           Each of the above three forms has relative merits:
!             - The original LP formulation has been well-tested in 
!               the community; however it involves subtracting two
!               large numbers when R is close to 0, allowing for 
!               possibly large roundoff errors.  Additionally, the 
!               coth = sinh/cosh function may be costly to perform 
!               repeatedly.
!             - The rational approximation alleviates the cost of 
!               intrinsic functions and catastrophic floating-point 
!               cancellation errors, but introduces a 7.2% error 
!               away from the original limiter.
!             - The new approximation also alleviates any catastrophic 
!               floating-point cancellation errors, and introduces 
!               only a 4.8% error away from the original limiter; 
!               however it still involves the use of possibly-expensive 
!               intrinsic functions.
!
!  INPUTS:
!     EgCur      - Gray radiation energy density
!     EgOld      - Gray radiation energy density (old time step)
!     Temp       - fluid temperature
!     kappaE     - energy mean absorption coefficient in cell
!     LimType    - INTG_PREC flag denoting type of flux limiter:
!                       0 -> standard Levermore-Pomraning lim. (LP, 1981)
!                       1 -> rational approx. to LP lim. (LP, 1981)
!                       2 -> Larsen n=2 limiter
!                       3 -> turns off limiter (constant of 1/3)
!                       4 -> Zeus limiter
!     a          - cosmological expansion parameter
!     *Units     - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     *face      - flag denoting whether the {x,y,z}*{l,r} subdomain 
!                  face lies on the overall domain exterior  
!                     0->interior, 1->exterior
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     rhs        - diffusive nonlinear rhs for the radiation 
!                  energy equation (in internal Enzo units)
!     ier        - success/failure output flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: LimType
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr, xlface, xrface
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr, ylface, yrface
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr, zlface, zrface
  INTG_PREC, intent(in)  :: Model
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a
  R_PREC, intent(in) :: aUnits, LenUnits, EgUnits
  R_PREC, intent(in) :: dx, dy, dz
  R_PREC, intent(in), target, &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       EgCur, EgOld, Temp
  R_PREC, intent(in), dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       kappaE
  R_PREC, intent(out) :: rhs(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)

  !--------------
  ! locals
  INTG_PREC  :: i, j, k
  R_PREC :: c, pi, dxi, dyi, dzi, Dlim
  R_PREC :: Egf, omega, R, sigT, AGradEg, StBz, Tf, Rmin
  R_PREC, dimension(3) :: GradEgL, GradEgR, DEgL, DEgR

  !=======================================================================

  ! initialize output to zero, flag to success
  rhs = 0._RKIND
  ier = 1

  ! set shortcut values 
  dxi = a/dx/LenUnits
  dyi = a/dy/LenUnits
  dzi = a/dz/LenUnits
  c = c_light        ! speed of light [cm/s]
  pi = pi_val
  StBz = 5.6704e-5_RKIND   ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  Rmin = 1.e-20_RKIND

  ! compute radiation energy gradient over domain
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           !--------------
           ! x-directional limiter, lower face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j,k) - EgOld(i-1,j,k))*dxi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i-1,j,k))/2._RKIND

           !    absorption, scattering, total extinction coeffs on face
           sigT = (kappaE(i,j,k)+kappaE(i-1,j,k))/2._RKIND

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i-1,j,k))/2._RKIND
              omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3._RKIND
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif
           DEgL(1) = Dlim

           !--------------
           ! x-directional limiter, upper face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i+1,j,k) - EgOld(i,j,k))*dxi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i+1,j,k))/2._RKIND

           !    absorption, scattering, total extinction coeffs on face
           sigT = (kappaE(i,j,k)+kappaE(i+1,j,k))/2._RKIND

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i+1,j,k))/2._RKIND
              omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits

              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3._RKIND
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif
           DEgR(1) = Dlim
           
           ! correct diffusion coefficient for Marshak boundaries
           if ((Model>=20) .and. (Model<30) .and. (i==Nx) .and. (xrface==1)) then
              DEgR(1) = 2._RKIND*DEgR(1) / (1._RKIND + 4._RKIND/3._RKIND*dxi/sigT)
           endif

           !--------------
           ! y-directional limiter, lower face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j,k) - EgOld(i,j-1,k))*dyi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i,j-1,k))/2._RKIND

           !    absorption, scattering, total extinction coeffs on face
           sigT = (kappaE(i,j,k)+kappaE(i,j-1,k))/2._RKIND

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i,j-1,k))/2._RKIND
              omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3._RKIND
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif
           DEgL(2) = Dlim

           !--------------
           ! y-directional limiter, upper face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j+1,k) - EgOld(i,j,k))*dyi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i,j+1,k))/2._RKIND

           !    absorption, scattering, total extinction coeffs on face
           sigT = (kappaE(i,j,k)+kappaE(i,j+1,k))/2._RKIND

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i,j+1,k))/2._RKIND
              omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3._RKIND
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)              
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif
           DEgR(2) = Dlim
           
           ! correct diffusion coefficient for Marshak boundaries
           if ((Model>=20) .and. (Model<30) .and. (j==Ny) .and. (yrface==1)) then
              DEgR(2) = 2._RKIND*DEgR(2) / (1._RKIND + 4._RKIND/3._RKIND*dyi/sigT)
           endif

           !--------------
           ! z-directional limiter, lower face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j,k) - EgOld(i,j,k-1))*dzi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i,j,k-1))/2._RKIND

           !    absorption, scattering, total extinction coeffs on face
           sigT = (kappaE(i,j,k)+kappaE(i,j,k-1))/2._RKIND

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i,j,k-1))/2._RKIND
              omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3._RKIND
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif
           DEgL(3) = Dlim

           !--------------
           ! z-directional limiter, upper face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j,k+1) - EgOld(i,j,k))*dzi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i,j,k+1))/2._RKIND

           !    absorption, scattering, total extinction coeffs on face
           sigT = (kappaE(i,j,k)+kappaE(i,j,k+1))/2._RKIND

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i,j,k+1))/2._RKIND
              omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits

              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3._RKIND
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif
           DEgR(3) = Dlim
           
           ! correct diffusion coefficient for Marshak boundaries
           if ((Model>=20) .and. (Model<30) .and. (k==Nz) .and. (zrface==1)) then
              DEgR(3) = 2._RKIND*DEgR(3) / (1._RKIND + 4._RKIND/3._RKIND*dzi/sigT)
           endif

           !--------------
           ! compute gradients of current Eg
           GradEgL(1) = (EgCur(i,j,k) - EgCur(i-1,j,k))*dxi
           GradEgR(1) = (EgCur(i+1,j,k) - EgCur(i,j,k))*dxi
           GradEgL(2) = (EgCur(i,j,k) - EgCur(i,j-1,k))*dyi
           GradEgR(2) = (EgCur(i,j+1,k) - EgCur(i,j,k))*dyi
           GradEgL(3) = (EgCur(i,j,k) - EgCur(i,j,k-1))*dzi
           GradEgR(3) = (EgCur(i,j,k+1) - EgCur(i,j,k))*dzi

           ! put it all together
           rhs(i,j,k) = ((DEgR(1)*GradEgR(1) - DEgL(1)*GradEgL(1))*dxi &
                       + (DEgR(2)*GradEgR(2) - DEgL(2)*GradEgL(2))*dyi &
                       + (DEgR(3)*GradEgR(3) - DEgL(3)*GradEgL(3))*dzi)

        enddo
     enddo
  enddo

  return
end subroutine gFLDProblem_DiffRHS_3D
!=======================================================================




subroutine gFLDProblem_DiffRHS_2D(rhs, EgCur, EgOld, Temp, kappaE,      &
     LimType, a, aUnits, LenUnits, EgUnits, dx, dy, Nx, Ny, NGxl, NGxr, &
     NGyl, NGyr, xlface, xrface, ylface, yrface, Model, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2006
!  modified1:  July 20, 2007, by John Hayes; applying transport correction
!              to diffusion coefficients on external boundaries.
!  modified2:  August, 2007, by John Hayes; 2D version cloned from original
!              routine
!  modified3:  December 21, 2007, by John Hayes; corrected formula for Marshak BC
!
!  PURPOSE: Computes the spatially-diffusive component of the 
!           nonlinear rhs for the Gray FLD problem,
!              -1/a^2*Div(D(Eg)*Grad(Eg))
!           where D(Eg) is a nonlinear flux-limiter depending on Eg.  
!           We define the values
!              R_i = |Grad(Eg)_i|/sigT/omega/Eg,
!              omega = B/E
!              sigT = total extinction coefficient.
!           The '_i' subscript implies the gradient in the ith 
!           direction; these quantities are all required at cell faces, 
!           as that is the location of the divergence calculations.
!           With these components, we allow any of the following three 
!           forms of the limiter, 
!             [Levermore-Pomraning, 1981],
!                 D_i(Eg) = c/sigT/R_i/omega*[coth(R_i)-1/R_i],
!             [rational approx. to above, Levermore-Pomraning, 1981],
!                 D_i(Eg) = c/sigT/omega*(2+R_i)/(6+3*R_i+R_i**2),
!             [Larsen n=2 limiter],
!                 D_i(Eg) = c/sqrt((3*sigT)**2 + R**2)
!             [Zeus form of rational approx. to LP],
!                 D_i(Eg) = 2/pi*c/sigT/R_i*atan(R_i*pi/6),
!                 R_i = |Grad(Eg)_i|/sigT/Eg,
!           Each of the above three forms has relative merits:
!             - The original LP formulation has been well-tested in 
!               the community; however it involves subtracting two
!               large numbers when R is close to 0, allowing for 
!               possibly large roundoff errors.  Additionally, the 
!               coth = sinh/cosh function may be costly to perform 
!               repeatedly.
!             - The rational approximation alleviates the cost of 
!               intrinsic functions and catastrophic floating-point 
!               cancellation errors, but introduces a 7.2% error 
!               away from the original limiter.
!             - The new approximation also alleviates any catastrophic 
!               floating-point cancellation errors, and introduces 
!               only a 4.8% error away from the original limiter; 
!               however it still involves the use of possibly-expensive 
!               intrinsic functions.
!
!  INPUTS:
!     EgCur      - Gray radiation energy density
!     EgOld      - Gray radiation energy density (old time step)
!     Temp       - fluid temperature
!     kappaE     - energy mean absorption coefficient in cell
!     LimType    - INTG_PREC flag denoting type of flux limiter:
!                       0 -> standard Levermore-Pomraning lim. (LP, 1981)
!                       1 -> rational approx. to LP lim. (LP, 1981)
!                       2 -> Larsen n=2 limiter
!                       3 -> turns off limiter (constant of 1/3)
!                       4 -> Zeus limiter
!     a          - cosmological expansion parameter
!     *Units     - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     *face      - flag denoting whether the {x,y,z}*{l,r} subdomain 
!                  face lies on the overall domain exterior  
!                     0->interior, 1->exterior
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     rhs        - diffusive nonlinear rhs for the radiation 
!                  energy equation (in internal Enzo units)
!     ier        - success/failure output flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: LimType
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr, xlface, xrface
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr, ylface, yrface
  INTG_PREC, intent(in)  :: Model
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a
  R_PREC, intent(in) :: aUnits, LenUnits, EgUnits
  R_PREC, intent(in) :: dx, dy
  R_PREC, intent(in), target, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr) :: &
       EgCur, EgOld, Temp
  R_PREC, intent(in), dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr) :: kappaE
  R_PREC, intent(out) :: rhs(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr)

  !--------------
  ! locals
  INTG_PREC :: i, j, k
  R_PREC :: c, pi, dxi, dyi, Dlim
  R_PREC :: Egf, omega, R, sigT, AGradEg, StBz, Tf, Rmin
  R_PREC, dimension(2) :: GradEgL, GradEgR, DEgL, DEgR

  !=======================================================================

  ! initialize output to zero, flag to success
  rhs = 0._RKIND
  ier = 1

  ! set shortcut values
  dxi = a/dx/LenUnits
  dyi = a/dy/LenUnits
  c = c_light        ! speed of light [cm/s]
  pi = pi_val
  StBz = 5.6704e-5_RKIND   ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  Rmin = 1.e-20_RKIND

  ! compute radiation energy gradient over domain
  do j=1,Ny,1
     do i=1,Nx,1

        !--------------
        ! x-directional limiter, lower face
        !    radiation energy gradient on face
        AGradEg = abs(EgOld(i,j) - EgOld(i-1,j))*dxi

        !    face-centered radiation energy value
        Egf = (EgOld(i,j) + EgOld(i-1,j))/2._RKIND

        !    absorption, scattering, total extinction coeffs on face
        sigT = (kappaE(i,j)+kappaE(i-1,j))/2._RKIND

        !    compute R for limiter based on LimType
        if ((LimType == 4) .or. (LimType == 2)) then
           R = AGradEg/Egf
           R = max(R,Rmin)
        else                          ! all others
           !    scaling coefficient ('effective albedo' -- LP)
           Tf = (Temp(i,j)+Temp(i-1,j))/2._RKIND
           omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits

           !    face-centered R value
           R = AGradEg/Egf/omega
           R = max(R,Rmin)  ! force away from 0 to avoid NaN
        endif

        !    compute limiter
        if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
        else if (LimType == 2) then  ! Larsen n=2 lim.
           Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
        else if (LimType == 3) then  ! no limiter
           Dlim = c/sigT/3._RKIND
        else if (LimType == 4) then  ! Zeus limiter
           Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
        else                         ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
        endif
        DEgL(1) = Dlim

        !--------------
        ! x-directional limiter, upper face
        !    radiation energy gradient on face
        AGradEg = abs(EgOld(i+1,j) - EgOld(i,j))*dxi
        
        !    face-centered radiation energy value
        Egf = (EgOld(i,j) + EgOld(i+1,j))/2._RKIND
        
        !    absorption, scattering, total extinction coeffs on face
        sigT = (kappaE(i,j)+kappaE(i+1,j))/2._RKIND

        !    compute R for limiter based on LimType
        if ((LimType == 4) .or. (LimType == 2)) then
           R = AGradEg/Egf
           R = max(R,Rmin)
        else                          ! all others
           !    scaling coefficient ('effective albedo' -- LP)
           Tf = (Temp(i,j)+Temp(i+1,j))/2._RKIND
           omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits
           
           !    face-centered R value
           R = AGradEg/Egf/omega
           R = max(R,Rmin)  ! force away from 0 to avoid NaN
        endif

        !    compute limiter
        if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
        else if (LimType == 2) then  ! Larsen n=2 lim.
           Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
        else if (LimType == 3) then  ! no limiter
           Dlim = c/sigT/3._RKIND
        else if (LimType == 4) then  ! Zeus limiter
           Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
        else                         ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
        endif
        DEgR(1) = Dlim
           
        ! correct diffusion coefficient for Marshak boundaries
        if ((Model>=20) .and. (Model<30) .and. (i==Nx) .and. (xrface==1)) then
           DEgR(1) = 2._RKIND*DEgR(1) / (1._RKIND + 4._RKIND/3._RKIND*dxi/sigT)
        endif

        !--------------
        ! y-directional limiter, lower face
        !    radiation energy gradient on face
        AGradEg = abs(EgOld(i,j) - EgOld(i,j-1))*dyi

        !    face-centered radiation energy value
        Egf = (EgOld(i,j) + EgOld(i,j-1))/2._RKIND

        !    absorption, scattering, total extinction coeffs on face
        sigT = (kappaE(i,j)+kappaE(i,j-1))/2._RKIND

        !    compute R for limiter based on LimType
        if ((LimType == 4) .or. (LimType == 2)) then
           R = AGradEg/Egf
           R = max(R,Rmin)
        else                          ! all others
           !    scaling coefficient ('effective albedo' -- LP)
           Tf = (Temp(i,j)+Temp(i,j-1))/2._RKIND
           omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits
           
           !    face-centered R value
           R = AGradEg/Egf/omega
           R = max(R,Rmin)  ! force away from 0 to avoid NaN
        endif

        !    compute limiter
        if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
        else if (LimType == 2) then  ! Larsen n=2 lim.
           Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
        else if (LimType == 3) then  ! no limiter
           Dlim = c/sigT/3._RKIND
        else if (LimType == 4) then  ! Zeus limiter
           Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
        else                         ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
        endif
        DEgL(2) = Dlim

        !--------------
        ! y-directional limiter, upper face
        !    radiation energy gradient on face
        AGradEg = abs(EgOld(i,j+1) - EgOld(i,j))*dyi

        !    face-centered radiation energy value
        Egf = (EgOld(i,j) + EgOld(i,j+1))/2._RKIND

        !    absorption, scattering, total extinction coeffs on face
        sigT = (kappaE(i,j)+kappaE(i,j+1))/2._RKIND

        !    compute R for limiter based on LimType
        if ((LimType == 4) .or. (LimType == 2)) then
           R = AGradEg/Egf
           R = max(R,Rmin)
        else                          ! all others
           !    scaling coefficient ('effective albedo' -- LP)
           Tf = (Temp(i,j)+Temp(i,j+1))/2._RKIND
           omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits
           
           !    face-centered R value
           R = AGradEg/Egf/omega
           R = max(R,Rmin)  ! force away from 0 to avoid NaN
        endif

        !    compute limiter
        if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
        else if (LimType == 2) then  ! Larsen n=2 lim.
           Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
        else if (LimType == 3) then  ! no limiter
           Dlim = c/sigT/3._RKIND
        else if (LimType == 4) then  ! Zeus limiter
           Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
        else                         ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
        endif
        DEgR(2) = Dlim
        
        ! correct diffusion coefficient for Marshak boundaries
        if ((Model>=20) .and. (Model<30) .and. (j==Ny) .and. (yrface==1)) then
           DEgR(2) = 2._RKIND*DEgR(2) / (1._RKIND + 4._RKIND/3._RKIND*dyi/sigT)
        endif

        !--------------
        ! compute gradients of current Eg
        GradEgL(1) = (EgCur(i,j) - EgCur(i-1,j))*dxi
        GradEgR(1) = (EgCur(i+1,j) - EgCur(i,j))*dxi
        GradEgL(2) = (EgCur(i,j) - EgCur(i,j-1))*dyi
        GradEgR(2) = (EgCur(i,j+1) - EgCur(i,j))*dyi

        ! put it all together
        rhs(i,j) = ((DEgR(1)*GradEgR(1) - DEgL(1)*GradEgL(1))*dxi &
                  + (DEgR(2)*GradEgR(2) - DEgL(2)*GradEgL(2))*dyi)

     enddo
  enddo

  return
end subroutine gFLDProblem_DiffRHS_2D
!=======================================================================




subroutine gFLDProblem_DiffRHS_1D(rhs, EgCur, EgOld, Temp, kappaE, &
     LimType, a, aUnits, LenUnits, EgUnits, dx, Nx, NGxl, NGxr,    &
     xlface, xrface, Model, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2006
!  modified1:  July 20, 2007, by John Hayes; applying transport correction
!              to diffusion coefficients on external boundaries.
!  modified2:  August, 2007, by John Hayes; 1D version cloned from original
!              routine
!  modified3:  December 21, 2007, by John Hayes; corrected formula for Marshak BC
!
!  PURPOSE: Computes the spatially-diffusive component of the 
!           nonlinear rhs for the Gray FLD problem,
!              -1/a^2*Div(D(Eg)*Grad(Eg))
!           where D(Eg) is a nonlinear flux-limiter depending on Eg.  
!           We define the values
!              R_i = |Grad(Eg)_i|/sigT/omega/Eg,
!              omega = B/E
!              sigT = total extinction coefficient.
!           The '_i' subscript implies the gradient in the ith 
!           direction; these quantities are all required at cell faces, 
!           as that is the location of the divergence calculations.
!           With these components, we allow any of the following three 
!           forms of the limiter, 
!             [Levermore-Pomraning, 1981],
!                 D_i(Eg) = c/sigT/R_i/omega*[coth(R_i)-1/R_i],
!             [rational approx. to above, Levermore-Pomraning, 1981],
!                 D_i(Eg) = c/sigT/omega*(2+R_i)/(6+3*R_i+R_i**2),
!             [Larsen n=2 limiter],
!                 D_i(Eg) = c/sqrt((3*sigT)**2 + R**2)
!             [Zeus form of rational approx. to LP],
!                 D_i(Eg) = 2/pi*c/sigT/R_i*atan(R_i*pi/6),
!                 R_i = |Grad(Eg)_i|/sigT/Eg,
!           Each of the above three forms has relative merits:
!             - The original LP formulation has been well-tested in 
!               the community; however it involves subtracting two
!               large numbers when R is close to 0, allowing for 
!               possibly large roundoff errors.  Additionally, the 
!               coth = sinh/cosh function may be costly to perform 
!               repeatedly.
!             - The rational approximation alleviates the cost of 
!               intrinsic functions and catastrophic floating-point 
!               cancellation errors, but introduces a 7.2% error 
!               away from the original limiter.
!             - The new approximation also alleviates any catastrophic 
!               floating-point cancellation errors, and introduces 
!               only a 4.8% error away from the original limiter; 
!               however it still involves the use of possibly-expensive 
!               intrinsic functions.
!
!  INPUTS:
!     EgCur      - Gray radiation energy density
!     EgOld      - Gray radiation energy density (old time step)
!     Temp       - fluid temperature
!     kappaE     - energy mean absorption coefficient in cell
!     LimType    - INTG_PREC flag denoting type of flux limiter:
!                       0 -> standard Levermore-Pomraning lim. (LP, 1981)
!                       1 -> rational approx. to LP lim. (LP, 1981)
!                       2 -> Larsen n=2 limiter
!                       3 -> turns off limiter (constant of 1/3)
!                       4 -> Zeus limiter
!     a          - cosmological expansion parameter
!     *Units     - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     *face      - flag denoting whether the {x,y,z}*{l,r} subdomain 
!                  face lies on the overall domain exterior  
!                     0->interior, 1->exterior
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     rhs        - diffusive nonlinear rhs for the radiation 
!                  energy equation (in internal Enzo units)
!     ier        - success/failure output flag (0->failure, 1->success)
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
  INTG_PREC, intent(in)  :: LimType
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr, xlface, xrface
  INTG_PREC, intent(in)  :: Model
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a
  R_PREC, intent(in) :: aUnits, LenUnits, EgUnits
  R_PREC, intent(in) :: dx
  R_PREC, intent(in), target, dimension(1-NGxl:Nx+NGxr) :: EgCur, EgOld, Temp
  R_PREC, intent(in), dimension(1-NGxl:Nx+NGxr) :: kappaE
  R_PREC, intent(out) :: rhs(1-NGxl:Nx+NGxr)

  !--------------
  ! locals
  INTG_PREC  :: i
  R_PREC :: c, pi, dxi, Dlim
  R_PREC :: Egf, omega, R, sigT, AGradEg, StBz, Tf
  R_PREC :: GradEgL, GradEgR, DEgL, DEgR, Rmin

  !=======================================================================

  ! initialize output to zero, flag to success
  rhs = 0._RKIND
  ier = 1

  ! set shortcut values
  dxi = a/dx/LenUnits
  c = c_light        ! speed of light [cm/s]
  pi = pi_val
  StBz = 5.6704e-5_RKIND   ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  Rmin = 1.0e-20_RKIND

  ! compute radiation energy gradient over domain
  do i=1,Nx,1

     !--------------
     ! x-directional limiter, lower face
     !    radiation energy gradient on face
     AGradEg = abs(EgOld(i) - EgOld(i-1))*dxi

     !    face-centered radiation energy value
     Egf = (EgOld(i) + EgOld(i-1))/2._RKIND

     !    absorption, scattering, total extinction coeffs on face
     sigT = (kappaE(i)+kappaE(i-1))/2._RKIND

     !    compute R for limiter based on LimType
     if ((LimType == 4) .or. (LimType == 2)) then
        R = AGradEg/Egf
        R = max(R,Rmin)
     else                       ! all others
        !    scaling coefficient ('effective albedo' -- LP)
        Tf = (Temp(i)+Temp(i-1))/2._RKIND
        omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits
        
        !    face-centered R value
        R = AGradEg/Egf/omega
        R = max(R,Rmin)  ! force away from 0 to avoid NaN
     endif

     !    compute limiter
     if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
        Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
     else if (LimType == 2) then  ! Larsen n=2 lim.
        Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
     else if (LimType == 3) then  ! no limiter
        Dlim = c/sigT/3._RKIND
     else if (LimType == 4) then  ! Zeus limiter
        Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
     else                         ! standard Levermore-Pomraning (LP, 1981)
        Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
     endif
     DEgL = Dlim

     !--------------
     ! x-directional limiter, upper face
     !    radiation energy gradient on face
     AGradEg = abs(EgOld(i+1) - EgOld(i))*dxi

     !    face-centered radiation energy value
     Egf = (EgOld(i) + EgOld(i+1))/2._RKIND

     !    absorption, scattering, total extinction coeffs on face
     sigT = (kappaE(i)+kappaE(i+1))/2._RKIND

     !    compute R for limiter based on LimType
     if ((LimType == 4) .or. (LimType == 2)) then
        R = AGradEg/Egf
        R = max(R,Rmin)
     else                       ! all others
        !    scaling coefficient ('effective albedo' -- LP)
        Tf = (Temp(i)+Temp(i+1))/2._RKIND
        omega = (4._RKIND*StBz/c*Tf**4)/Egf/EgUnits
        
        !    face-centered R value
        R = AGradEg/Egf/omega
        R = max(R,Rmin)  ! force away from 0 to avoid NaN
     endif

     !    compute limiter
     if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
        Dlim = c/omega*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
     else if (LimType == 2) then  ! Larsen n=2 lim.
        Dlim = c/sqrt((3._RKIND*sigT)**2 + R**2)
     else if (LimType == 3) then  ! no limiter
        Dlim = c/sigT/3._RKIND
     else if (LimType == 4) then  ! Zeus limiter
        Dlim = c*(2._RKIND*sigT+R)/(6._RKIND*sigT*sigT+3._RKIND*sigT*R+R*R)
     else                         ! standard Levermore-Pomraning (LP, 1981)
        Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
     endif
     DEgR = Dlim
        
     ! correct diffusion coefficient for Marshak boundaries
     if ((Model<=20) .and. (Model<30) .and. (i==Nx) .and. (xrface==1)) then
        DEgR = 2._RKIND*DEgR / (1._RKIND + 4._RKIND/3._RKIND*dxi/sigT)
     endif

     !--------------
     ! compute gradients of current Eg
     GradEgL = (EgCur(i) - EgCur(i-1))*dxi
     GradEgR = (EgCur(i+1) - EgCur(i))*dxi

     ! put it all together
     rhs(i) = (DEgR*GradEgR - DEgL*GradEgL)*dxi

  enddo

  return
end subroutine gFLDProblem_DiffRHS_1D
!=======================================================================

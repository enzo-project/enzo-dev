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
subroutine gFLDProblem_MatrixEntries_3D(matentries, EgCur, EgOld, Temp,  &
     kappaE, adjvec, LimType, dt, a, theta, aUnits, LenUnits, EgUnits,   &
     dx, dy, dz, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr,   &
     NGyl, NGyr, NGzl, NGzr, xrface, yrface, zrface, Model, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2006
!  modified1:  July 20, 2007, by John Hayes; applying transport correction
!              to diffusion coefficients on external boundaries.
!  modified2:  August 10, 2007 by John Hayes; appended "_3D" to routine name.
!  modified3:  December 21, 2007 by John Hayes; corrected Marshak BC.
!
!  PURPOSE: Computes the array of matrix stencil elements for the 
!           Gray FLD problem,
!              -dt/a*Div(D(Eg)*Grad(Eg))
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
!
!           As the stencil has 7 non-zero elements per matrix row, we 
!           set these entries over the computational domain, with the 
!           proper adjustments due to the choice of limiter.
!
!  INPUTS:
!     EgCur      - Gray radiation energy density
!     EgOld      - Gray radiation energy density (old time step)
!     Temp       - fluid temperature
!     kappaE     - Energy mean absorption coefficient in cell
!     adjvec     - Schur complement adjustment vector (for diagonal)
!     LimType    - integer flag denoting type of flux limiter:
!                       0 -> standard Levermore-Pomraning lim. (LP, 1981)
!                       1 -> rational approx. to LP lim. (LP, 1981)
!                       2 -> Larsen n=2 limiter
!                       3 -> turns off the limiter (constant of 1/3)
!                       4 -> Zeus limiter
!     a          - cosmological expansion parameter
!     dt         - time step size
!     theta      - overall implicitness parameter
!     *Units     - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     x*{s,e}    - start/end indices of linear solver domain; 
!                  typically 1:Nx for standard dims, but Dirichlet 
!                  BCs may move these to 0:Nx, 1:Nx+1, etc.
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
!     matentries - array of stencil values over the active domain.  
!                  Since the stencil has 7 nonzero entries, and as 
!                  this array should not include ghost cells, it 
!                  has dimensions (7,Nx,Ny,Nz).
!     ier        - success/failure flag (0->failure, 1->success)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
#include "fortran.def"
  implicit none

!--------------
! argument declarations
  integer, intent(in)  :: LimType
  integer, intent(in)  :: Nx, NGxl, NGxr, xrface, x0s, x0e
  integer, intent(in)  :: Ny, NGyl, NGyr, yrface, x1s, x1e
  integer, intent(in)  :: Nz, NGzl, NGzr, zrface, x2s, x2e
  integer, intent(in)  :: Model
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a
  real,    intent(in) :: dx, dy, dz, dt, theta
  real,    intent(in) :: aUnits, LenUnits, EgUnits
  real, target, intent(in),                                &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) &
       :: EgCur, EgOld
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) &
       :: Temp, kappaE
  real*8, intent(out) :: matentries(7,x0s:x0e,x1s:x1e,x2s:x2e)
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: adjvec

!--------------
! locals
  integer :: i, j, k
  real :: c, pi, StBz, dxi, dyi, dzi, dtfac
  real :: Egf, omega, R, sigT, AGradEg, Tf
  real :: Dlim
  real :: Rmin, acoef

!=======================================================================
  
!!$  write(*,*) 'Entering gFLDProblem::MatrixEntries_3D routine'

  ! initialize output flag, and set matentries to have all zero values
  ier = 1
  matentries = 0.d0

  ! set shortcut values
  dtfac = dt*theta    ! time step conversion factor
  dxi   = a/dx/LenUnits
  dyi   = a/dy/LenUnits
  dzi   = a/dz/LenUnits
  c  = 2.99792458d10     ! speed of light [cm/s]
  StBz = 5.6704d-5       ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  pi = 4.d0*datan(1.d0)
  Rmin = 1.0e-20

  ! iterate over the active domain
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           ! initialize matrix entries
           matentries(:,i,j,k) = 0.d0
           matentries(4,i,j,k) = adjvec(i,j,k)

           !--------------
           ! x-directional limiter, lower face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j,k) - EgOld(i-1,j,k))*dxi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i-1,j,k))/2.d0

           !    total extinction coeff on face
           sigT = (kappaE(i,j,k)+kappaE(i-1,j,k))/2.d0

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i-1,j,k))/2.d0
              omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3.d0
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on x-left Eg
           matentries(3,i,j,k) = matentries(3,i,j,k) - dtfac*Dlim*dxi*dxi
           !       dep. on self Eg
           matentries(4,i,j,k) = matentries(4,i,j,k) + dtfac*Dlim*dxi*dxi


           !--------------
           ! x-directional limiter, upper face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i+1,j,k) - EgOld(i,j,k))*dxi
           
           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i+1,j,k))/2.d0
           
           !    total extinction coeff on face
           sigT = (kappaE(i,j,k)+kappaE(i+1,j,k))/2.d0

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i+1,j,k))/2.d0
              omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3.d0
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif

           ! correct diffusion coefficient for Marshak boundaries
           if ((Model>=20) .and. (Model<30) .and. (i==Nx) .and. (xrface==1)) then
              acoef = 4.d0/3.d0*dxi/kappaE(i,j,k)
              Dlim = 2.0*Dlim / (1.0 + acoef)
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on x-right Eg
           matentries(5,i,j,k) = matentries(5,i,j,k) - dtfac*Dlim*dxi*dxi
           !       dep. on self Eg
           matentries(4,i,j,k) = matentries(4,i,j,k) + dtfac*Dlim*dxi*dxi


           !--------------
           ! y-directional limiter, lower face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j,k) - EgOld(i,j-1,k))*dyi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i,j-1,k))/2.d0

           !    total extinction coeff on face
           sigT = (kappaE(i,j,k)+kappaE(i,j-1,k))/2.d0

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i,j-1,k))/2.d0
              omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3.d0
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on y-left Eg
           matentries(2,i,j,k) = matentries(2,i,j,k) - dtfac*Dlim*dyi*dyi
           !       dep. on self Eg
           matentries(4,i,j,k) = matentries(4,i,j,k) + dtfac*Dlim*dyi*dyi


           !--------------
           ! y-directional limiter, upper face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j+1,k) - EgOld(i,j,k))*dyi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i,j+1,k))/2.d0

           !    total extinction coeff on face
           sigT = (kappaE(i,j,k)+kappaE(i,j+1,k))/2.d0

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i,j+1,k))/2.d0
              omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3.d0
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif

           ! correct diffusion coefficient for Marshak boundaries
           if ((Model>=20) .and. (Model<30) .and. (j==Ny) .and. (yrface==1)) then
              acoef = 4.d0/3.d0*dyi/kappaE(i,j,k)
              Dlim = 2.0*Dlim / (1.0 + acoef)
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on y-right Eg
           matentries(6,i,j,k) = matentries(6,i,j,k) - dtfac*Dlim*dyi*dyi
           !       dep. on self Eg
           matentries(4,i,j,k) = matentries(4,i,j,k) + dtfac*Dlim*dyi*dyi


           !--------------
           ! z-directional limiter, lower face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j,k) - EgOld(i,j,k-1))*dzi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i,j,k-1))/2.d0

           !    total extinction coeff on face
           sigT = (kappaE(i,j,k)+kappaE(i,j,k-1))/2.d0

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i,j,k-1))/2.d0
              omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3.d0
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on z-left Eg
           matentries(1,i,j,k) = matentries(1,i,j,k) - dtfac*Dlim*dzi*dzi
           !       dep. on self Eg
           matentries(4,i,j,k) = matentries(4,i,j,k) + dtfac*Dlim*dzi*dzi


           !--------------
           ! z-directional limiter, upper face
           !    radiation energy gradient on face
           AGradEg = abs(EgOld(i,j,k+1) - EgOld(i,j,k))*dzi

           !    face-centered radiation energy value
           Egf = (EgOld(i,j,k) + EgOld(i,j,k+1))/2.d0

           !    total extinction coeff on face
           sigT = (kappaE(i,j,k)+kappaE(i,j,k+1))/2.d0

           !    compute R for limiter based on LimType
           if ((LimType == 4) .or. (LimType == 2)) then
              R = AGradEg/Egf
              R = max(R,Rmin)
           else                             ! all others
              !    scaling coefficient ('effective albedo' -- LP)
              Tf = (Temp(i,j,k)+Temp(i,j,k+1))/2.d0
              omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
              
              !    face-centered R value
              R = AGradEg/Egf/omega
              R = max(R,Rmin)  ! force away from 0 to avoid NaN
           endif

           !    compute limiter
           if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R) 
           else if (LimType == 2) then  ! Larsen n=2 lim.
              Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
           else if (LimType == 3) then  ! no limiter
              Dlim = c/sigT/3.d0
           else if (LimType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R) 
           else                         ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
           endif

           ! correct diffusion coefficient for Marshak boundaries
           if ((Model>=20) .and. (Model<30) .and. (k==Nz) .and. (zrface==1)) then
              acoef = 4.d0/3.d0*dzi/kappaE(i,j,k)
              Dlim = 2.0*Dlim / (1.0 + acoef)
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on z-right Eg
           matentries(7,i,j,k) = matentries(7,i,j,k) - dtfac*Dlim*dzi*dzi
           !       dep. on self Eg
           matentries(4,i,j,k) = matentries(4,i,j,k) + dtfac*Dlim*dzi*dzi

        enddo
     enddo
  enddo

  return
end subroutine gFLDProblem_MatrixEntries_3D
!=======================================================================





subroutine gFLDProblem_MatrixEntries_2D(matentries, EgCur, EgOld, Temp,  &
     kappaE, adjvec, LimType, dt, a, theta, aUnits, LenUnits, EgUnits,   &
     dx, dy, x0s, x0e, x1s, x1e, Nx, Ny, NGxl, NGxr, NGyl, NGyr, xrface, &
     yrface, Model, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2006
!  modified1:  July 20, 2007, by John Hayes; applying transport correction
!              to diffusion coefficients on external boundaries.
!  modified2:  August 10, 2007, by John Hayes; cloned 2D version from
!              original routine.
!  modified3:  December 21, 2007 by John Hayes; corrected Marshak BC.
!
!  PURPOSE: Computes the array of matrix stencil elements for the 
!           Gray FLD problem,
!              -dt/a*Div(D(Eg)*Grad(Eg))
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
!
!           As the stencil has 5 non-zero elements per matrix row, we 
!           set these entries over the computational domain, with the 
!           proper adjustments due to the choice of limiter.
!
!  INPUTS:
!     EgCur      - Gray radiation energy density
!     EgOld      - Gray radiation energy density (old time step)
!     Temp       - fluid temperature
!     kappaE     - Energy mean absorption coefficient in cell
!     adjvec     - Schur complement adjustment vector (for diagonal)
!     LimType    - integer flag denoting type of flux limiter:
!                       0 -> standard Levermore-Pomraning lim. (LP, 1981)
!                       1 -> rational approx. to LP lim. (LP, 1981)
!                       2 -> Larsen n=2 limiter
!                       3 -> turns off the limiter (constant of 1/3)
!                       4 -> Zeus limiter
!     a          - cosmological expansion parameter
!     dt         - time step size
!     theta      - overall implicitness parameter
!     *Units     - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     x*{s,e}    - start/end indices of linear solver domain; 
!                  typically 1:Nx for standard dims, but Dirichlet 
!                  BCs may move these to 0:Nx, 1:Nx+1, etc.
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
!     matentries - array of stencil values over the active domain.  
!                  Since the stencil has 7 nonzero entries, and as 
!                  this array should not include ghost cells, it 
!                  has dimensions (7,Nx,Ny,Nz).
!     ier        - success/failure flag (0->failure, 1->success)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
#include "fortran.def"
  implicit none

!--------------
! argument declarations
  integer, intent(in)  :: LimType
  integer, intent(in)  :: Nx, NGxl, NGxr, xrface, x0s, x0e
  integer, intent(in)  :: Ny, NGyl, NGyr, yrface, x1s, x1e
  integer, intent(in)  :: Model
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a
  real,    intent(in) :: dx, dy, dt, theta
  real,    intent(in) :: aUnits, LenUnits, EgUnits
  real, target, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), intent(in) &
       :: EgCur, EgOld
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), intent(in) &
       :: Temp, kappaE
  real*8, intent(out) :: matentries(5,x0s:x0e,x1s:x1e)
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr) :: adjvec

!--------------
! locals
  integer :: i, j
  real :: c, pi, StBz, dxi, dyi, dtfac
  real :: Egf, omega, R, sigT, AGradEg, Tf
  real :: Dlim
  real :: Rmin, acoef

!=======================================================================
  
!!$  write(*,*) 'Entering gFLDProblem::MatrixEntries_2D routine'

  ! initialize output flag, and set matentries to have all zero values
  ier = 1
  matentries = 0.d0

  ! set shortcut values
  dtfac = dt*theta       ! time step conversion factor
  dxi   = a/dx/LenUnits
  dyi   = a/dy/LenUnits
  c  = 2.99792458d10     ! speed of light [cm/s]
  StBz = 5.6704d-5       ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  pi = 4.d0*datan(1.d0)
  Rmin = 1.0e-20

  ! iterate over the active domain
  do j=1,Ny,1
     do i=1,Nx,1

        ! initialize matrix entries
        matentries(:,i,j) = 0.d0
        matentries(3,i,j) = adjvec(i,j)

        !--------------
        ! x-directional limiter, lower face
        !    radiation energy gradient on face
        AGradEg = abs(EgOld(i,j) - EgOld(i-1,j))*dxi

        !    face-centered radiation energy value
        Egf = (EgOld(i,j) + EgOld(i-1,j))/2.d0

        !    total extinction coeff on face
        sigT = (kappaE(i,j)+kappaE(i-1,j))/2.d0

        !    compute R for limiter based on LimType
        if ((LimType == 4) .or. (LimType == 2)) then
           R = AGradEg/Egf
           R = max(R,Rmin)
        else                             ! all others
           !    scaling coefficient ('effective albedo' -- LP)
           Tf = (Temp(i,j)+Temp(i-1,j))/2.d0
           omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
           
           !    face-centered R value
           R = AGradEg/Egf/omega
           R = max(R,Rmin)  ! force away from 0 to avoid NaN
        endif
        
        !    compute limiter
        if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
        else if (LimType == 2) then  ! Larsen n=2 lim.
           Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
        else if (LimType == 3) then  ! no limiter
           Dlim = c/sigT/3.d0
        else if (LimType == 4) then  ! Zeus limiter
           Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
        else                         ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
        endif

        !    set the relevant matrix entries. Note: the diffusive component 
        !    need not be rescaled, since scaling and chain rule cancel 
        !       dep. on x-left Eg
        matentries(2,i,j) = matentries(2,i,j) - dtfac*Dlim*dxi*dxi
        !       dep. on self Eg
        matentries(3,i,j) = matentries(3,i,j) + dtfac*Dlim*dxi*dxi


        !--------------
        ! x-directional limiter, upper face
        !    radiation energy gradient on face
        AGradEg = abs(EgOld(i+1,j) - EgOld(i,j))*dxi
        
        !    face-centered radiation energy value
        Egf = (EgOld(i,j) + EgOld(i+1,j))/2.d0
        
        !    total extinction coeff on face
        sigT = (kappaE(i,j)+kappaE(i+1,j))/2.d0

        !    compute R for limiter based on LimType
        if ((LimType == 4) .or. (LimType == 2)) then
           R = AGradEg/Egf
           R = max(R,Rmin)
        else                             ! all others
           !    scaling coefficient ('effective albedo' -- LP)
           Tf = (Temp(i,j)+Temp(i+1,j))/2.d0
           omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
           
           !    face-centered R value
           R = AGradEg/Egf/omega
           R = max(R,Rmin)  ! force away from 0 to avoid NaN
        endif
        
        !    compute limiter
        if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
        else if (LimType == 2) then  ! Larsen n=2 lim.
           Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
        else if (LimType == 3) then  ! no limiter
           Dlim = c/sigT/3.d0
        else if (LimType == 4) then  ! Zeus limiter
           Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
        else                         ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
        endif

        ! correct diffusion coefficient for Marshak boundaries
        if ((Model>=20) .and. (Model<30) .and. (i==Nx) .and. (xrface==1)) then
           acoef = 4.d0/3.d0*dxi/kappaE(i,j)
           Dlim = 2.0*Dlim / (1.0 + acoef)
        endif

        !    set the relevant matrix entries. Note: the diffusive component 
        !    need not be rescaled, since scaling and chain rule cancel 
        !       dep. on x-right Eg
        matentries(4,i,j) = matentries(4,i,j) - dtfac*Dlim*dxi*dxi
        !       dep. on self Eg
        matentries(3,i,j) = matentries(3,i,j) + dtfac*Dlim*dxi*dxi


        !--------------
        ! y-directional limiter, lower face
        !    radiation energy gradient on face
        AGradEg = abs(EgOld(i,j) - EgOld(i,j-1))*dyi

        !    face-centered radiation energy value
        Egf = (EgOld(i,j) + EgOld(i,j-1))/2.d0

        !    total extinction coeff on face
        sigT = (kappaE(i,j)+kappaE(i,j-1))/2.d0

        !    compute R for limiter based on LimType
        if ((LimType == 4) .or. (LimType == 2)) then
           R = AGradEg/Egf
           R = max(R,Rmin)
        else                             ! all others
           !    scaling coefficient ('effective albedo' -- LP)
           Tf = (Temp(i,j)+Temp(i,j-1))/2.d0
           omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
           
           !    face-centered R value
           R = AGradEg/Egf/omega
           R = max(R,Rmin)  ! force away from 0 to avoid NaN
        endif
        
        !    compute limiter
        if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
        else if (LimType == 2) then  ! Larsen n=2 lim.
           Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
        else if (LimType == 3) then  ! no limiter
           Dlim = c/sigT/3.d0
        else if (LimType == 4) then  ! Zeus limiter
           Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
        else                         ! standard Levermore-Pomraning lim. (LP, 1981)
           Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
        endif

        !    set the relevant matrix entries. Note: the diffusive component 
        !    need not be rescaled, since scaling and chain rule cancel 
        !       dep. on y-left Eg
        matentries(1,i,j) = matentries(1,i,j) - dtfac*Dlim*dyi*dyi
        !       dep. on self Eg
        matentries(3,i,j) = matentries(3,i,j) + dtfac*Dlim*dyi*dyi


        !--------------
        ! y-directional limiter, upper face
        !    radiation energy gradient on face
        AGradEg = abs(EgOld(i,j+1) - EgOld(i,j))*dyi

        !    face-centered radiation energy value
        Egf = (EgOld(i,j) + EgOld(i,j+1))/2.d0

        !    total extinction coeff on face
        sigT = (kappaE(i,j)+kappaE(i,j+1))/2.d0

        !    compute R for limiter based on LimType
        if ((LimType == 4) .or. (LimType == 2)) then
           R = AGradEg/Egf
           R = max(R,Rmin)
        else                             ! all others
           !    scaling coefficient ('effective albedo' -- LP)
           omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
           
           !    face-centered R value
           R = AGradEg/Egf/omega
           R = max(R,Rmin)  ! force away from 0 to avoid NaN
        endif
        
        !    compute limiter
        if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
        else if (LimType == 2) then  ! Larsen n=2 lim.
           Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
        else if (LimType == 3) then  ! no limiter
           Dlim = c/sigT/3.d0
        else if (LimType == 4) then  ! Zeus limiter
           Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
        else                         ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
        endif

        ! correct diffusion coefficient for Marshak boundaries
        if ((Model>=20) .and. (Model<30) .and. (j==Ny) .and. (yrface==1)) then
           acoef = 4.d0/3.d0*dyi/kappaE(i,j)
           Dlim = 2.0*Dlim / (1.0 + acoef)
        endif

        !    set the relevant matrix entries. Note: the diffusive component 
        !    need not be rescaled, since scaling and chain rule cancel 
        !       dep. on y-right Eg
        matentries(5,i,j) = matentries(5,i,j) - dtfac*Dlim*dyi*dyi
        !       dep. on self Eg
        matentries(3,i,j) = matentries(3,i,j) + dtfac*Dlim*dyi*dyi

      enddo
   enddo

  return
end subroutine gFLDProblem_MatrixEntries_2D
!=======================================================================




subroutine gFLDProblem_MatrixEntries_1D(matentries, EgCur, EgOld, Temp, &
     kappaE, adjvec, LimType, dt, a, theta, aUnits, LenUnits, EgUnits,  &
     dx, x0s, x0e, Nx, NGxl, NGxr, xrface, Model, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2006
!  modified1:  July 20, 2007, by John Hayes; applying transport correction
!              to diffusion coefficients on external boundaries.
!  modified2:  August 10, 2007, by John Hayes; cloned 1D version from
!              original routine.
!  modified3:  December 21, 2007 by John Hayes; corrected Marshak BC.
!
!  PURPOSE: Computes the array of matrix stencil elements for the 
!           Gray FLD problem,
!              -dt/a*Div(D(Eg)*Grad(Eg))
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
!
!           As the stencil has 5 non-zero elements per matrix row, we 
!           set these entries over the computational domain, with the 
!           proper adjustments due to the choice of limiter.
!
!  INPUTS:
!     EgCur      - Gray radiation energy density
!     EgOld      - Gray radiation energy density (old time step)
!     Temp       - fluid temperature
!     kappaE     - Energy mean absorption coefficient in cell
!     adjvec     - Schur complement adjustment vector (for diagonal)
!     LimType    - integer flag denoting type of flux limiter:
!                       0 -> standard Levermore-Pomraning lim. (LP, 1981)
!                       1 -> rational approx. to LP lim. (LP, 1981)
!                       2 -> Larsen n=2 limiter
!                       3 -> turns off the limiter (constant of 1/3)
!                       4 -> Zeus limiter
!     a          - cosmological expansion parameter
!     dt         - time step size
!     theta      - overall implicitness parameter
!     *Units     - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     x*{s,e}    - start/end indices of linear solver domain; 
!                  typically 1:Nx for standard dims, but Dirichlet 
!                  BCs may move these to 0:Nx, 1:Nx+1, etc.
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
!     matentries - array of stencil values over the active domain.  
!                  Since the stencil has 7 nonzero entries, and as 
!                  this array should not include ghost cells, it 
!                  has dimensions (7,Nx,Ny,Nz).
!     ier        - success/failure flag (0->failure, 1->success)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
#include "fortran.def"
  implicit none

!--------------
! argument declarations
  integer, intent(in)  :: LimType
  integer, intent(in)  :: Nx, NGxl, NGxr, xrface, x0s, x0e
  integer, intent(in)  :: Model
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a
  real,    intent(in)  :: dx, dt, theta
  real,    intent(in)  :: aUnits, LenUnits, EgUnits
  real,    intent(in), target, dimension(1-NGxl:Nx+NGxr) :: EgCur, EgOld
  real,    intent(in), dimension(1-NGxl:Nx+NGxr) :: Temp, kappaE
  real*8,  intent(out) :: matentries(3,x0s:x0e)
  real, dimension(1-NGxl:Nx+NGxr) :: adjvec

!--------------
! locals
  integer :: i
  real :: c, pi, StBz, dxi, dtfac
  real :: Egf, omega, R, sigT, AGradEg, Tf
  real :: Dlim
  real :: Rmin, acoef

!=======================================================================
  
!!$  write(*,*) 'Entering gFLDProblem::MatrixEntries_1D routine'

  ! initialize output flag, and set matentries to have all zero values
  ier = 1
  matentries = 0.d0

  ! set shortcut values
  dtfac = dt*theta       ! time step conversion factor
  dxi   = a/dx/LenUnits
  c  = 2.99792458d10     ! speed of light [cm/s]
  StBz = 5.6704d-5       ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]
  pi = 4.d0*datan(1.d0)
  Rmin = 1.0e-20

  ! iterate over the active domain
  do i=1,Nx,1

     ! initialize matrix entries
     matentries(:,i) = 0.d0
     matentries(2,i) = adjvec(i)

     !--------------
     ! x-directional limiter, lower face
     !    radiation energy gradient on face
     AGradEg = abs(EgOld(i) - EgOld(i-1))*dxi

     !    face-centered radiation energy value
     Egf = (EgOld(i) + EgOld(i-1))/2.d0

     !    total extinction coeff on face
     sigT = (kappaE(i)+kappaE(i-1))/2.d0

     !    compute R for limiter based on LimType
     if ((LimType == 4) .or. (LimType == 2)) then
        R = AGradEg/Egf
        R = max(R,Rmin)
     else                             ! all others
        !    scaling coefficient ('effective albedo' -- LP)
        Tf = (Temp(i)+Temp(i-1))/2.d0
        omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
        
        !    face-centered R value
        R = AGradEg/Egf/omega
        R = max(R,Rmin)  ! force away from 0 to avoid NaN
     endif
        
     !    compute limiter
     if (LimType == 1) then        ! rational approx. to LP lim. (LP, 1981)
        Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
     else if (LimType == 2) then  ! Larsen n=2 lim.
        Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
     else if (LimType == 3) then  ! no limiter
        Dlim = c/sigT/3.d0
     else if (LimType == 4) then  ! Zeus limiter
        Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
     else                         ! standard Levermore-Pomraning (LP, 1981)
        Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
     endif

     !    set the relevant matrix entries. Note: the diffusive component 
     !    need not be rescaled, since scaling and chain rule cancel 
     !       dep. on x-left Eg
     matentries(1,i) = matentries(1,i) - dtfac*Dlim*dxi*dxi
     !       dep. on self Eg
     matentries(2,i) = matentries(2,i) + dtfac*Dlim*dxi*dxi


     !--------------
     ! x-directional limiter, upper face
     !    radiation energy gradient on face
     AGradEg = abs(EgOld(i+1) - EgOld(i))*dxi
     
     !    face-centered radiation energy value
     Egf = (EgOld(i) + EgOld(i+1))/2.d0
     
     !    total extinction coeff on face
     sigT = (kappaE(i)+kappaE(i+1))/2.d0

     !    compute R for limiter based on LimType
     if ((LimType == 4) .or. (LimType == 2)) then
        R = AGradEg/Egf
        R = max(R,Rmin)
     else                             ! all others
        !    scaling coefficient ('effective albedo' -- LP)
        Tf = (Temp(i)+Temp(i+1))/2.d0
        omega = (4.d0*StBz/c*Tf**4)/Egf/EgUnits
        
        !    face-centered R value
        R = AGradEg/Egf/omega
        R = max(R,Rmin)  ! force away from 0 to avoid NaN
     endif
        
     !    compute limiter
     if (LimType == 1) then       ! rational approx. to LP lim. (LP, 1981)
        Dlim = c/omega*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
     else if (LimType == 2) then  ! Larsen n=2 lim.
        Dlim = c/sqrt((3.d0*sigT)**2 + R**2)
     else if (LimType == 3) then  ! no limiter
        Dlim = c/sigT/3.d0
     else if (LimType == 4) then  ! Zeus limiter
        Dlim = c*(2.d0*sigT+R)/(6.d0*sigT*sigT+3.d0*sigT*R+R*R)
     else                         ! standard Levermore-Pomraning (LP, 1981)
        Dlim = c/omega*(cosh(R/sigT)/sinh(R/sigT)-sigT/R)/R
     endif

     ! correct diffusion coefficient for Marshak boundaries
     if ((Model>=20) .and. (Model<30) .and. (i==Nx) .and. (xrface==1)) then
        acoef = 4.d0/3.d0*dxi/kappaE(i)
        Dlim = 2.0*Dlim / (1.0 + acoef)
     endif

     !    set the relevant matrix entries. Note: the diffusive component 
     !    need not be rescaled, since scaling and chain rule cancel 
     !       dep. on x-right Eg
     matentries(3,i) = matentries(3,i) - dtfac*Dlim*dxi*dxi
     !       dep. on self Eg
     matentries(2,i) = matentries(2,i) + dtfac*Dlim*dxi*dxi

   enddo

  return
end subroutine gFLDProblem_MatrixEntries_1D
!=======================================================================

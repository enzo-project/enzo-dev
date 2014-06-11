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


function gFLD_limiter(E1, E2, k1, k2, nUn, lUn, dxi)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       February 2013
!  modified:   
!
!  PURPOSE: Computes the flux limiter at a given face
!=======================================================================
  implicit none
#include "fortran_types.def"
  R_PREC, intent(in) :: E1, E2, k1, k2, nUn, lUn, dxi
  R_PREC :: gFLD_limiter, Eavg, kap, R, Emin, Rmin, Dmax
  
  ! set limiter bounds
  Rmin = 1.e-2_RKIND/lUn
  Rmin = min(Rmin, 1.e-20_RKIND)    ! 1st is astro/cosmo, 2nd is lab frame
  Emin = 1.e-30_RKIND
  Dmax = 0.002_RKIND * c_light * lUn
  Dmax = max(Dmax, 1.e20_RKIND)     ! 1st is astro/cosmo, 2nd is lab frame

  ! compute limiter
  Eavg = max((E1 + E2)*0.5_RKIND, Emin)
  kap = 2._RKIND*k1*k2/(k1+k2)*nUn        ! harmonic mean
  R = max(dxi*abs(E1 - E2)/Eavg, Rmin)
  gFLD_limiter = min(c_light/sqrt(9._RKIND*kap*kap + R*R), Dmax)

end function gFLD_limiter
!=======================================================================



subroutine gFLDSplit_SetupSystem(matentries, rhsentries, rhsnorm, E0,   &
     E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,     &
     theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, rank, dx, dy, dz,     &
     BCXl, BCXr, BCYl, BCYr, BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e,  &
     Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, xlface, xrface,    &
     ylface, yrface, zlface, zrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!  modified:   
!
!  PURPOSE: Computes the array of matrix stencil elements and vector of 
!           rhs entries for the Grey FLD radiation problem,
!              d_t E - Div(D(E)*Grad(E)) = -adot/a*E - c*kappa*E + eta + src
!           where D(E) is a nonlinear flux-limiter depending on E0
!           (time lagged).  As the stencil has {7,5,3} non-zero elements 
!           per matrix row (depending on whether the problem is 3D, 2D 
!           or 1D), we set these entries over the computational domain.
!
!           We in fact solve a scaled version of the equation.  Since 
!           the values of E are in fact in normalized units 
!           (E_true = E*rUn), we must scale src by rUn to achieve the
!           correct equation.  Moreover, we do not solve the equation 
!           directly, and instead solve for a correction to the current
!           state such that the corrected solution satisfies the above 
!           equation.  This helps with enforcement of boundary conditions,
!           since they may be directly placed onto the current state, 
!           and the correction need only refrain from interfering.
!
!  INPUTS:
!     E0         - Grey radiation energy density (prev time step)
!     E          - Grey radiation energy density (current guess)
!     Temp       - gas temperature for black-body radiation
!     Temp0      - gas temperature (prev time step)
!     kappa      - opacity array
!     src        - spatially-dependent radiation source
!     dt         - time step size
!     a,a0       - cosmological expansion factor (new and old time steps)
!     adot,adot0 - da/dt and da0/dt
!     ESpectrum  - flag denoting what type of radiation field we have:
!                  (-1=>monochromatic)
!     theta      - overall implicitness parameter
!     *Un,*Un0   - variable scaling constants (new and old time steps)
!     rank       - 1, 2 or 3; the dimensionality of the problem
!     dx,dy,dz   - mesh spacing in each direction
!     BC*        - boundary condition type in each direction, face
!                     0->periodic
!                     1->Dirichlet
!                     2->Neumann
!     x*{s,e}    - start/end indices of linear solver domain; 
!                  typically 1:Nx for standard dims, but Dirichlet 
!                  BCs may move these to 0:Nx, 1:Nx+1, etc.
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     *{l,r}face - integer flag denoting whether direction/face 
!                  is external to the domain (0->int, 1->ext)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     matentries - array of stencil values over the active domain.  
!                  Since the stencil has 7 nonzero entries, and as 
!                  this array should not include ghost cells, it 
!                  has dimensions (7,x0s:x0e,x1s:x1e,x2s:x2e).
!     rhsentries - array of rhs values over the active domain.  As 
!                  this array should not include ghost cells, it 
!                  has dimensions (x0s:x0e,x1s:x1e,x2s:x2e)
!     rhsnorm    - 2-norm of rhs array
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
  INTG_PREC, intent(in)  :: rank, ESpectrum
  INTG_PREC, intent(in)  :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  INTG_PREC, intent(in)  :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  INTG_PREC, intent(in)  :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  P_PREC,    intent(in)  :: a, a0, adot, adot0
  R_PREC,    intent(in)  :: dt, theta, dx, dy, dz
  R_PREC,    intent(in)  :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0
  R_PREC,    intent(in)  :: E0(*), E(*), Temp(*), Temp0(*), kappa(*), src(*)
  REAL*8,    intent(out) :: matentries(*)
  REAL*8,    intent(out) :: rhsentries(*)
  R_PREC,    intent(out) :: rhsnorm
  INTG_PREC, intent(out) :: ier

  !=======================================================================
  

  ! call the apprpriate routine based on rank
  if (rank == 3) then

     call gFLDSplit_SetupSystem3D(matentries, rhsentries, rhsnorm, E0,    &
          E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,  &
          theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, dx, dy, dz, BCXl,  &
          BCXr, BCYl, BCYr, BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, &
          Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, xlface, xrface,     &
          ylface, yrface, zlface, zrface, ier)

  elseif (rank == 2) then

     call gFLDSplit_SetupSystem2D(matentries, rhsentries, rhsnorm, E0,    &
          E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,  &
          theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, dx, dy, BCXl,      &
          BCXr, BCYl, BCYr, x0s, x0e, x1s, x1e, Nx, Ny, NGxl, NGxr, NGyl, &
          NGyr, xlface, xrface, ylface, yrface, ier)

  elseif (rank == 1) then

     call gFLDSplit_SetupSystem1D(matentries, rhsentries, rhsnorm, E0,    &
          E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,  &
          theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, dx, BCXl, BCXr,    &
          x0s, x0e, Nx, NGxl, NGxr, xlface, xrface, ier)

  else
     write(0,*) 'gFLDSplit_SetupSystem error: illegal rank =',rank
  end if

end subroutine gFLDSplit_SetupSystem
!=======================================================================






subroutine gFLDSplit_SetupSystem3D(matentries, rhsentries, rhsnorm, E0, &
     E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,     &
     theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, dx, dy, dz, BCXl,     &
     BCXr, BCYl, BCYr, BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx,    &
     Ny, Nz, NGxl, NGxr, NGyl, NGyr,  NGzl, NGzr, xlface, xrface,       &
     ylface, yrface, zlface, zrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June 2009
!  modified:   
!
!  PURPOSE: 3D version of the routine
!=======================================================================
  implicit none
#include "fortran_types.def"
  
  !--------------
  ! argument declarations
  INTG_PREC, intent(in) :: ESpectrum
  INTG_PREC, intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  INTG_PREC, intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  INTG_PREC, intent(in) :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  P_PREC,    intent(in) :: a, a0, adot, adot0
  R_PREC,    intent(in) :: dt, theta, dx, dy, dz
  R_PREC,    intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) &
       :: E0, E, kappa, src, Temp, Temp0
  REAL*8,    intent(out) :: matentries(7,x0s:x0e,x1s:x1e,x2s:x2e)
  REAL*8,    intent(out) :: rhsentries(x0s:x0e,x1s:x1e,x2s:x2e)
  R_PREC,    intent(out) :: rhsnorm
  INTG_PREC, intent(out) :: ier

  !--------------
  ! locals
  INTG_PREC :: i, j, k
  R_PREC  :: dxi, dxi0, dyi, dyi0, dzi, dzi0
  REAL*8  :: dtfac, dtfac0, kap, kap0, eta, eta0, StBz, afac, afac0
  REAL*8  :: dxfac, dyfac, dzfac, dxfac0, dyfac0, dzfac0
  REAL*8  :: D_xl, D0_xl, D_xr, D0_xr, E0d_xl, E0d_xr, Ed_xl, Ed_xr
  REAL*8  :: D_yl, D0_yl, D_yr, D0_yr, E0d_yl, E0d_yr, Ed_yl, Ed_yr
  REAL*8  :: D_zl, D0_zl, D_zr, D0_zr, E0d_zl, E0d_zr, Ed_zl, Ed_zr
  REAL*8, external :: gFLD_limiter

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  matentries = 0._RKIND
  rhsentries = 0._RKIND
  rhsnorm    = 0._RKIND
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta         ! time step conversion factor
  dtfac0 = dt*(1._RKIND-theta)  ! time step conversion factor
  if (ESpectrum == -1) then
     afac  = 0._RKIND
     afac0 = 0._RKIND
  else
     afac  = adot/a         ! expansion factor (new time)
     afac0 = adot0/a0       ! expansion factor (old time)
  endif
  dxi   = 1._RKIND/dx/lUn
  dyi   = 1._RKIND/dy/lUn
  dzi   = 1._RKIND/dz/lUn
  dxi0  = 1._RKIND/dx/lUn0
  dyi0  = 1._RKIND/dy/lUn0
  dzi0  = 1._RKIND/dz/lUn0
  dxfac = dtfac*dxi*dxi
  dyfac = dtfac*dyi*dyi
  dzfac = dtfac*dzi*dzi
  dxfac0 = dtfac0*dxi0*dxi0
  dyfac0 = dtfac0*dyi0*dyi0
  dzfac0 = dtfac0*dzi0*dzi0
  StBz  = 5.6704e-5_RKIND   ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]


  ! iterate over the active domain
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           !--------------
           ! z-directional limiter, lower face
           !    compute gradients of E0, Ediff
           E0d_zl = E0(i,j,k) - E0(i,j,k-1)
           Ed_zl  = E(i,j,k) - E(i,j,k-1)
           D0_zl  = gFLD_limiter(E0(i,j,k), E0(i,j,k-1), kappa(i,j,k), &
                            kappa(i,j,k-1), nUn0, lUn0, dzi)
           D_zl   = gFLD_limiter(E(i,j,k), E(i,j,k-1), kappa(i,j,k), &
                            kappa(i,j,k-1), nUn, lUn, dzi)

           !--------------
           ! y-directional limiter, lower face
           !    compute gradients of E0, Ediff
           E0d_yl = E0(i,j,k) - E0(i,j-1,k)
           Ed_yl  = E(i,j,k) - E(i,j-1,k)
           D0_yl  = gFLD_limiter(E0(i,j,k), E0(i,j-1,k), kappa(i,j,k), &
                            kappa(i,j-1,k), nUn0, lUn0, dyi)
           D_yl   = gFLD_limiter(E(i,j,k), E(i,j-1,k), kappa(i,j,k), &
                            kappa(i,j-1,k), nUn, lUn, dyi)

           !--------------
           ! x-directional limiter, lower face
           !    compute gradients of E0, Ediff
           E0d_xl = E0(i,j,k) - E0(i-1,j,k)
           Ed_xl  = E(i,j,k) - E(i-1,j,k)
           D0_xl  = gFLD_limiter(E0(i,j,k), E0(i-1,j,k), kappa(i,j,k), &
                            kappa(i-1,j,k), nUn0, lUn0, dxi)
           D_xl   = gFLD_limiter(E(i,j,k), E(i-1,j,k), kappa(i,j,k), &
                            kappa(i-1,j,k), nUn, lUn, dxi)

           !--------------
           ! x-directional limiter, upper face
           !    compute gradients of E0, Ediff
           E0d_xr = E0(i+1,j,k) - E0(i,j,k)
           Ed_xr  = E(i+1,j,k) - E(i,j,k)
           D0_xr  = gFLD_limiter(E0(i,j,k), E0(i+1,j,k), kappa(i,j,k), &
                            kappa(i+1,j,k), nUn0, lUn0, dxi)
           D_xr   = gFLD_limiter(E(i,j,k), E(i+1,j,k), kappa(i,j,k), &
                            kappa(i+1,j,k), nUn, lUn, dxi)

           !--------------
           ! y-directional limiter, upper face
           !    compute gradients of E0, Ediff
           E0d_yr = E0(i,j+1,k) - E0(i,j,k)
           Ed_yr  = E(i,j+1,k) - E(i,j,k)
           D0_yr  = gFLD_limiter(E0(i,j,k), E0(i,j+1,k), kappa(i,j,k), &
                            kappa(i,j+1,k), nUn0, lUn0, dyi)
           D_yr   = gFLD_limiter(E(i,j,k), E(i,j+1,k), kappa(i,j,k), &
                            kappa(i,j+1,k), nUn, lUn, dyi)

           !--------------
           ! z-directional limiter, upper face
           !    compute gradients of E0, Ediff
           E0d_zr = E0(i,j,k+1) - E0(i,j,k)
           Ed_zr  = E(i,j,k+1) - E(i,j,k)
           D0_zr  = gFLD_limiter(E0(i,j,k), E0(i,j,k+1), kappa(i,j,k), &
                            kappa(i,j,k+1), nUn0, lUn0, dzi)
           D_zr   = gFLD_limiter(E(i,j,k), E(i,j,k+1), kappa(i,j,k), &
                            kappa(i,j,k+1), nUn, lUn, dzi)

           ! opacity values in this cell
           kap  = kappa(i,j,k)*nUn
           kap0 = kappa(i,j,k)*nUn0

           ! black-body radiation in this cell (if applicable; otherwise Temp=0)
           eta  = 4._RKIND*kap*StBz*Temp(i,j,k)**4
           eta0 = 4._RKIND*kap0*StBz*Temp0(i,j,k)**4

           ! set the matrix entries.  Note: the diffusive component 
           ! need not be rescaled, since scaling and chain rule cancel 
           matentries(1,i,j,k) = -dzfac*D_zl         ! z-left
           matentries(2,i,j,k) = -dyfac*D_yl         ! y-left
           matentries(3,i,j,k) = -dxfac*D_xl         ! x-left
           matentries(4,i,j,k) = 1._RKIND + dtfac*(afac + c_light*kap)  &   ! self
                + dxfac*(D_xl+D_xr) + dyfac*(D_yl+D_yr) + dzfac*(D_zl+D_zr)
           matentries(5,i,j,k) = -dxfac*D_xr         ! x-right
           matentries(6,i,j,k) = -dyfac*D_yr         ! y-right
           matentries(7,i,j,k) = -dzfac*D_zr         ! z-right

           ! set the rhs entries
           rhsentries(i,j,k) =   dtfac/rUn*(eta + src(i,j,k))                     &
                               + dtfac0/rUn0*(eta0 + src(i,j,k))                  &
                               - (1._RKIND + dtfac*(afac+c_light*kap))*E(i,j,k)       &
                               + dxfac*(D_xr*Ed_xr-D_xl*Ed_xl)                    &
                               + dyfac*(D_yr*Ed_yr-D_yl*Ed_yl)                    &
                               + dzfac*(D_zr*Ed_zr-D_zl*Ed_zl)                    &
                               + (1._RKIND - dtfac0*(afac0+c_light*kap0))*E0(i,j,k)   &
                               + dxfac0*(D0_xr*E0d_xr-D0_xl*E0d_xl)               &
                               + dyfac0*(D0_yr*E0d_yr-D0_yl*E0d_yl)               &
                               + dzfac0*(D0_zr*E0d_zr-D0_zl*E0d_zl)
           rhsnorm = rhsnorm + rhsentries(i,j,k)**2

        enddo
     enddo
  enddo


  ! update matrix/rhs based on boundary conditions/location
  !    z-left face
  if (zlface == 1) then
     ! Dirichlet
     if (BCZl==1) then
        k = x2s
        do j=1,Ny
           do i=1,Nx
              matentries(1,i,j,k) = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCZl==2) then
        k = x2s
        do j=1,Ny
           do i=1,Nx
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(1,i,j,k)
              matentries(1,i,j,k) = 0._RKIND
           enddo
        enddo
     endif
  end if

  !    y-left face
  if (ylface == 1) then
     ! Dirichlet
     if (BCYl==1) then
        j = x1s
        do k=1,Nz
           do i=1,Nx
              matentries(2,i,j,k) = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCYl==2) then
        j = x1s
        do k=1,Nz
           do i=1,Nx
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(2,i,j,k)
              matentries(2,i,j,k) = 0._RKIND
           enddo
        enddo
     endif
  end if

  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        do k=1,Nz
           do j=1,Ny
              matentries(3,i,j,k) = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        do k=1,Nz
           do j=1,Ny
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(3,i,j,k)
              matentries(3,i,j,k) = 0._RKIND
           enddo
        enddo
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        do k=1,Nz
           do j=1,Ny
              matentries(5,i,j,k) = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        do k=1,Nz
           do j=1,Ny
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(5,i,j,k)
              matentries(5,i,j,k) = 0._RKIND
           enddo
        enddo
     endif
  endif

  !    y-right face
  if (yrface==1) then
     ! Dirichlet
     if (BCYr==1) then
        j = x1e
        do k=1,Nz
           do i=1,Nx
              matentries(6,i,j,k) = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCYr==2) then
        j = x1e
        do k=1,Nz
           do i=1,Nx
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(6,i,j,k)
              matentries(6,i,j,k) = 0._RKIND
           enddo
        enddo
     endif
  endif

  !    z-right face
  if (zrface==1) then
     ! Dirichlet
     if (BCZr==1) then
        k = x2e
        do j=1,Ny
           do i=1,Nx
              matentries(7,i,j,k) = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCZr==2) then
        k = x2e
        do j=1,Ny
           do i=1,Nx
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(7,i,j,k)
              matentries(7,i,j,k) = 0._RKIND
           enddo
        enddo
     endif
  endif

  return
end subroutine gFLDSplit_SetupSystem3D
!=======================================================================






subroutine gFLDSplit_SetupSystem2D(matentries, rhsentries, rhsnorm, E0,   &
     E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,       &
     theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, dx, dy, BCXl, BCXr,     &
     BCYl, BCYr, x0s, x0e, x1s, x1e, Nx, Ny, NGxl, NGxr, NGyl, NGyr,      &
     xlface, xrface, ylface, yrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!  modified:   
!
!  PURPOSE: 2D version of the routine
!=======================================================================
  implicit none
#include "fortran_types.def"
  
  !--------------
  ! argument declarations
  INTG_PREC, intent(in) :: ESpectrum
  INTG_PREC, intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  INTG_PREC, intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  P_PREC,    intent(in) :: a, a0, adot, adot0
  R_PREC,    intent(in) :: dt, theta, dx, dy
  R_PREC,    intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), intent(in) &
       :: E0, E, src, kappa, Temp, Temp0
  REAL*8,    intent(out) :: matentries(5,x0s:x0e,x1s:x1e)
  REAL*8,    intent(out) :: rhsentries(x0s:x0e,x1s:x1e)
  R_PREC,    intent(out) :: rhsnorm
  INTG_PREC, intent(out) :: ier

  !--------------
  ! locals
  INTG_PREC :: i, j
  R_PREC  :: dxi, dxi0, dyi, dyi0
  REAL*8  :: dtfac, dtfac0, kap, kap0, eta, eta0, StBz, afac, afac0
  REAL*8  :: dxfac, dyfac, dxfac0, dyfac0
  REAL*8  :: D_xl, D0_xl, D_xr, D0_xr, E0d_xl, E0d_xr, Ed_xl, Ed_xr
  REAL*8  :: D_yl, D0_yl, D_yr, D0_yr, E0d_yl, E0d_yr, Ed_yl, Ed_yr
  REAL*8, external :: gFLD_limiter

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  matentries = 0._RKIND
  rhsentries = 0._RKIND
  rhsnorm    = 0._RKIND
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta         ! time step conversion factor
  dtfac0 = dt*(1._RKIND-theta)  ! time step conversion factor
  if (ESpectrum == -1) then
     afac  = 0._RKIND
     afac0 = 0._RKIND
  else
     afac  = adot/a         ! expansion factor (new time)
     afac0 = adot0/a0       ! expansion factor (old time)
  endif
  dxi   = 1._RKIND/dx/lUn
  dyi   = 1._RKIND/dy/lUn
  dxi0  = 1._RKIND/dx/lUn0
  dyi0  = 1._RKIND/dy/lUn0
  dxfac = dtfac*dxi*dxi
  dyfac = dtfac*dyi*dyi
  dxfac0 = dtfac0*dxi0*dxi0
  dyfac0 = dtfac0*dyi0*dyi0
  StBz  = 5.6704e-5_RKIND   ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]

  ! iterate over the active domain
  do j=1,Ny,1
     do i=1,Nx,1

        !--------------
        ! y-directional limiter, lower face
        !    compute gradients of E0, Ediff
        E0d_yl = E0(i,j) - E0(i,j-1)
        Ed_yl  = E(i,j) - E(i,j-1)
        D0_yl  = gFLD_limiter(E0(i,j), E0(i,j-1), kappa(i,j), &
                         kappa(i,j-1), nUn0, lUn0, dyi)
        D_yl   = gFLD_limiter(E(i,j), E(i,j-1), kappa(i,j), &
                         kappa(i,j-1), nUn, lUn, dyi)

        !--------------
        ! x-directional limiter, lower face
        !    compute gradients of E0, Ediff
        E0d_xl = E0(i,j) - E0(i-1,j)
        Ed_xl  = E(i,j) - E(i-1,j)
        D0_xl  = gFLD_limiter(E0(i,j), E0(i-1,j), kappa(i,j), &
                         kappa(i-1,j), nUn0, lUn0, dxi)
        D_xl   = gFLD_limiter(E(i,j), E(i-1,j), kappa(i,j), &
                         kappa(i-1,j), nUn, lUn, dxi)

        !--------------
        ! x-directional limiter, upper face
        !    compute gradients of E0, Ediff
        E0d_xr = E0(i+1,j) - E0(i,j)
        Ed_xr  = E(i+1,j) - E(i,j)
        D0_xr  = gFLD_limiter(E0(i,j), E0(i+1,j), kappa(i,j), &
                         kappa(i+1,j), nUn0, lUn0, dxi)
        D_xr   = gFLD_limiter(E(i,j), E(i+1,j), kappa(i,j), &
                         kappa(i+1,j), nUn, lUn, dxi)

        !--------------
        ! y-directional limiter, upper face
        !    compute gradients of E0, Ediff
        E0d_yr = E0(i,j+1) - E0(i,j)
        Ed_yr  = E(i,j+1) - E(i,j)
        D0_yr  = gFLD_limiter(E0(i,j), E0(i,j+1), kappa(i,j), &
                         kappa(i,j+1), nUn0, lUn0, dyi)
        D_yr   = gFLD_limiter(E(i,j), E(i,j+1), kappa(i,j), &
                         kappa(i,j+1), nUn, lUn, dyi)

        ! opacity values in this cell
        kap  = kappa(i,j)*nUn
        kap0 = kappa(i,j)*nUn0

        ! black-body radiation in this cell (if applicable)
        eta  = 4._RKIND*kap*StBz*Temp(i,j)**4
        eta0 = 4._RKIND*kap0*StBz*Temp0(i,j)**4

        ! set the matrix entries.  Note: the diffusive component 
        ! need not be rescaled, since scaling and chain rule cancel 
        matentries(1,i,j) = -dyfac*D_yl         ! y-left
        matentries(2,i,j) = -dxfac*D_xl         ! x-left
        matentries(3,i,j) = 1._RKIND + dtfac*(afac + c_light*kap)   &   
             + dxfac*(D_xl+D_xr) + dyfac*(D_yl+D_yr)      ! self
        matentries(4,i,j) = -dxfac*D_xr         ! x-right
        matentries(5,i,j) = -dyfac*D_yr         ! y-right

        ! set the rhs entries
        rhsentries(i,j) =   dtfac/rUn*(eta + src(i,j))                   &
                          + dtfac0/rUn0*(eta0 + src(i,j))                &
                          - (1._RKIND + dtfac*(afac+c_light*kap))*E(i,j)     &
                          + dxfac*(D_xr*Ed_xr-D_xl*Ed_xl)                &
                          + dyfac*(D_yr*Ed_yr-D_yl*Ed_yl)                &
                          + (1._RKIND - dtfac0*(afac0+c_light*kap0))*E0(i,j) &
                          + dxfac0*(D0_xr*E0d_xr-D0_xl*E0d_xl)           &
                          + dyfac0*(D0_yr*E0d_yr-D0_yl*E0d_yl)  
        rhsnorm = rhsnorm + rhsentries(i,j)**2

     enddo
  enddo
  
  ! update matrix/rhs based on boundary conditions/location
  !    y-left face
  if (ylface == 1) then
     ! Dirichlet
     if (BCYl==1) then
        j = x1s
        do i=1,Nx
           matentries(1,i,j) = 0._RKIND
        enddo
        ! Neumann
     else if (BCYl==2) then
        j = x1s
        do i=1,Nx
           matentries(3,i,j) = matentries(3,i,j) + matentries(1,i,j)
           matentries(1,i,j) = 0._RKIND
        enddo
     endif
  end if

  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        do j=1,Ny
           matentries(2,i,j) = 0._RKIND
        enddo
        ! Neumann
     else if (BCXl==2) then
        i = x0s
        do j=1,Ny
           matentries(3,i,j) = matentries(3,i,j) + matentries(2,i,j)
           matentries(2,i,j) = 0._RKIND
        enddo
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        do j=1,Ny
           matentries(4,i,j) = 0._RKIND
        enddo
        ! Neumann
     else if (BCXr==2) then
        i = x0e
        do j=1,Ny
           matentries(3,i,j) = matentries(3,i,j) + matentries(4,i,j)
           matentries(4,i,j) = 0._RKIND
        enddo
     endif
  endif

  !    y-right face
  if (yrface==1) then
     ! Dirichlet
     if (BCYr==1) then
        j = x1e
        do i=1,Nx
           matentries(5,i,j) = 0._RKIND
        enddo
        ! Neumann
     else if (BCYr==2) then
        j = x1e
        do i=1,Nx
           matentries(3,i,j) = matentries(3,i,j) + matentries(5,i,j)
           matentries(5,i,j) = 0._RKIND
        enddo
     endif
  endif

  return
end subroutine gFLDSplit_SetupSystem2D
!=======================================================================






subroutine gFLDSplit_SetupSystem1D(matentries, rhsentries, rhsnorm, E0, &
     E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,     &
     theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, dx, BCXl, BCXr, x0s,  &
     x0e, Nx, NGxl, NGxr, xlface, xrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June 2009
!  modified:   
!
!  PURPOSE: 1D version of the routine
!=======================================================================
  implicit none
#include "fortran_types.def"
  
  !--------------
  ! argument declarations
  INTG_PREC, intent(in) :: ESpectrum
  INTG_PREC, intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  P_PREC,    intent(in) :: a, a0, adot, adot0
  R_PREC,    intent(in) :: dt, theta, dx
  R_PREC,    intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0
  R_PREC, dimension(1-NGxl:Nx+NGxr), intent(in) :: E0, E, src, kappa, Temp, Temp0
  REAL*8,    intent(out) :: matentries(3,x0s:x0e)
  REAL*8,    intent(out) :: rhsentries(x0s:x0e)
  R_PREC,    intent(out) :: rhsnorm
  INTG_PREC, intent(out) :: ier

  !--------------
  ! locals
  INTG_PREC :: i
  R_PREC  :: dxi, dxi0
  REAL*8  :: dtfac, dtfac0, kap, kap0, eta, eta0, StBz, afac, afac0
  REAL*8  :: dxfac, dxfac0
  REAL*8  :: D_xl, D0_xl, D_xr, D0_xr, E0d_xl, E0d_xr, Ed_xl, Ed_xr
  REAL*8, external :: gFLD_limiter

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  matentries = 0._RKIND
  rhsentries = 0._RKIND
  rhsnorm    = 0._RKIND
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta         ! time step conversion factor
  dtfac0 = dt*(1._RKIND-theta)  ! time step conversion factor
  if (ESpectrum == -1) then
     afac  = 0._RKIND
     afac0 = 0._RKIND
  else
     afac  = adot/a         ! expansion factor (new time)
     afac0 = adot0/a0       ! expansion factor (old time)
  endif
  dxi   = 1._RKIND/dx/lUn
  dxi0  = 1._RKIND/dx/lUn0
  dxfac = dtfac*dxi*dxi
  dxfac0 = dtfac0*dxi0*dxi0
  StBz  = 5.6704e-5_RKIND   ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]

  ! iterate over the active domain
  do i=1,Nx,1

     !--------------
     ! x-directional limiter, lower face
     !    compute gradients of E0, Ediff
     E0d_xl = E0(i) - E0(i-1)
     Ed_xl  = E(i) - E(i-1)
     D0_xl  = gFLD_limiter(E0(i), E0(i-1), kappa(i), kappa(i-1), nUn0, lUn0, dxi)
     D_xl   = gFLD_limiter(E(i), E(i-1), kappa(i), kappa(i-1), nUn, lUn, dxi)

     !--------------
     ! x-directional limiter, upper face
     !    compute gradients of E0, Ediff
     E0d_xr = E0(i+1) - E0(i)
     Ed_xr  = E(i+1) - E(i)
     D0_xr  = gFLD_limiter(E0(i), E0(i+1), kappa(i), kappa(i+1), nUn0, lUn0, dxi)
     D_xr   = gFLD_limiter(E(i), E(i+1), kappa(i), kappa(i+1), nUn, lUn, dxi)

     ! opacity values in this cell
     kap  = kappa(i)*nUn
     kap0 = kappa(i)*nUn0

     ! black-body radiation in this cell (if applicable)
     eta  = 4._RKIND*kap*StBz*Temp(i)**4
     eta0 = 4._RKIND*kap0*StBz*Temp0(i)**4

     ! set the matrix entries.  Note: the diffusive component 
     ! need not be rescaled, since scaling and chain rule cancel 
     matentries(1,i) = -dxfac*D_xl            ! x-left
     matentries(2,i) = 1._RKIND + dtfac*(afac + c_light*kap)  &  ! self
                       + dxfac*(D_xl+D_xr)
     matentries(3,i) = -dxfac*D_xr            ! x-right

     ! set the rhs entries
     rhsentries(i) =   dtfac/rUn*(eta + src(i))                   &
                     + dtfac0/rUn0*(eta0 + src(i))                &
                     - (1._RKIND + dtfac*(afac+c_light*kap))*E(i)     &
                     + dxfac*(D_xr*Ed_xr-D_xl*Ed_xl)              &
                     + (1._RKIND - dtfac0*(afac0+c_light*kap0))*E0(i) &
                     + dxfac0*(D0_xr*E0d_xr-D0_xl*E0d_xl)
     rhsnorm = rhsnorm + rhsentries(i)**2

  enddo

  ! update matrix/rhs based on boundary conditions/location
  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        matentries(1,i) = 0._RKIND
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        matentries(2,i) = matentries(2,i) + matentries(1,i)
        matentries(1,i) = 0._RKIND
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        matentries(3,i) = 0._RKIND
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        matentries(2,i) = matentries(2,i) + matentries(3,i)
        matentries(3,i) = 0._RKIND
     endif
  endif

  return
end subroutine gFLDSplit_SetupSystem1D
!=======================================================================

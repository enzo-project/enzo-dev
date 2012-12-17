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
subroutine FSProb_SetupSystem(mat, rhs, rhsnorm, E, E0, kappa_h2on,       &
     kappa_arr, kappa_c, eta, dt, a, a0, adot, adot0, theta, lUn, lUn0,   &
     rUn, rUn0, nUn, nUn0, rank, dx, dy, dz, BCXl, BCXr, BCYl, BCYr,      &
     BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr,    &
     NGyl, NGyr, NGzl, NGzr, xlface, xrface, ylface, yrface, zlface,      &
     zrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       January 2010
!  modified:   
!
!  PURPOSE: Computes the array of matrix stencil elements and vector of 
!           rhs entries for the free-streaming radiation problem,
!              dot(E) - Div(c*E*sign(grad(E))) + adot/a*E + kappa*E = eta
!           where we have the [small] parameter
!              kappa = absorption coefficient (can be array-valued).
!           As the stencil has {7,5,3} non-zero elements per matrix row 
!           (depending on whether the problem is 3D, 2D or 1D), we 
!           set these entries over the computational domain.
!
!           We in fact solve a scaled version of the equation.  Since 
!           the values of E are in fact in normalized units 
!           (E_true = E*rUn), we must scale eta by rUn to achieve the
!           correct equation.  Moreover, we do not solve the equation 
!           directly, and instead solve for a correction to the current 
!           state such that the corrected solution satisfies the above 
!           equations.  This helps with enforcement of boundary conditions, 
!           since they may be directly placed into the current state, and 
!           the corrections need only refrain from interfering.
!
!  INPUTS:
!     E,E0       - Free-streaming radiation energy density (new & old time steps)
!     kappa_h2on - 0 -> use constant opacity
!                  1 -> use spatially dependent opacity
!     kappa_c    - opacity (constant)
!     kappa_arr  - opacity (spatially dependent)
!     eta        - spatially-dependent emissivity source
!     dt         - time step size
!     a,a0       - cosmological expansion factor (new and old time steps)
!     adot,adot0 - da/dt and da0/dt
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
!     *{l,r}face - INTG_PREC flag denoting whether direction/face 
!                  is external to the domain (0->int, 1->ext)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     mat - array of stencil values over the active domain.  
!                  Since the stencil has 7 nonzero entries, and as 
!                  this array should not include ghost cells, it 
!                  has dimensions (7,x0s:x0e,x1s:x1e,x2s:x2e).
!     rhs - array of rhs values over the active domain.  As 
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
  INTG_PREC, intent(in) :: rank, kappa_h2on
  INTG_PREC, intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  INTG_PREC, intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  INTG_PREC, intent(in) :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  P_PREC, intent(in) :: a, a0, adot, adot0
  R_PREC,    intent(in) :: kappa_c, dt, theta, dx, dy, dz
  R_PREC,    intent(in) :: lUn, lUn0, rUn, rUn0, nUn, nUn0
  R_PREC,    intent(in) :: E(*), E0(*), eta(*), kappa_arr(*)
  REAL*8,  intent(out) :: mat(*)
  REAL*8,  intent(out) :: rhs(*)
  R_PREC,    intent(out) :: rhsnorm
  INTG_PREC, intent(out) :: ier

  !=======================================================================
  

  ! call the apprpriate routine based on rank
  if (rank == 3) then

     call FSProb_SetupSystem3D(mat, rhs, rhsnorm, E, E0, kappa_h2on,        &
          kappa_arr, kappa_c, eta, dt, a, a0, adot, adot0, theta, lUn,      &
          lUn0, rUn, rUn0, nUn, nUn0, dx, dy, dz, BCXl, BCXr, BCYl, BCYr,   &
          BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, &
          NGyl, NGyr, NGzl, NGzr, xlface, xrface, ylface, yrface, zlface,   &
          zrface, ier)

  elseif (rank == 2) then

     call FSProb_SetupSystem2D(mat, rhs, rhsnorm, E, E0, kappa_h2on,        &
          kappa_arr, kappa_c, eta, dt, a, a0, adot, adot0, theta, lUn,      &
          lUn0, rUn, rUn0, nUn, nUn0, dx, dy, BCXl, BCXr, BCYl, BCYr, x0s,  &
          x0e, x1s, x1e, Nx, Ny, NGxl, NGxr, NGyl, NGyr, xlface, xrface,    &
          ylface, yrface, ier)

  elseif (rank == 1) then

     call FSProb_SetupSystem1D(mat, rhs, rhsnorm, E, E0, kappa_h2on,        &
          kappa_arr, kappa_c, eta, dt, a, a0, adot, adot0, theta, lUn,      &
          lUn0, rUn, rUn0, nUn, nUn0, dx, BCXl, BCXr, x0s, x0e, Nx, NGxl,   &
          NGxr, xlface, xrface, ier)

  else
     write(0,*) 'FSProb_SetupSystem error: illegal rank =',rank
  end if

end subroutine FSProb_SetupSystem
!=======================================================================






subroutine FSProb_SetupSystem3D(mat, rhs, rhsnorm, E, E0, kappa_h2on,    &
     kappa, kappa_c, eta, dt, a, a0, adot, adot0, theta, lUn, lUn0, rUn, &
     rUn0, nUn, nUn0, dx, dy, dz, BCXl, BCXr, BCYl, BCYr, BCZl, BCZr,    &
     x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,   &
     NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       January 2010
!  modified:   
!
!  PURPOSE: 3D version of the routine
!=======================================================================
  implicit none
#include "fortran_types.def"
  
  !--------------
  ! argument declarations
  INTG_PREC, intent(in) :: kappa_h2on
  INTG_PREC, intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  INTG_PREC, intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  INTG_PREC, intent(in) :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  P_PREC, intent(in) :: a, a0, adot, adot0
  R_PREC,    intent(in) :: kappa_c, dt, theta, dx, dy, dz
  R_PREC,    intent(in) :: lUn, lUn0, rUn, rUn0, nUn, nUn0
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) &
                       :: E, E0, eta, kappa
  REAL*8,  intent(out) :: mat(7,x0s:x0e,x1s:x1e,x2s:x2e)
  REAL*8,  intent(out) :: rhs(x0s:x0e,x1s:x1e,x2s:x2e)
  R_PREC,    intent(out) :: rhsnorm
  INTG_PREC, intent(out) :: ier

  !--------------
  ! locals
  INTG_PREC :: i, j, k
  REAL*8  :: dtfac, dtfac0, kap, kap0, eps, mu
  REAL*8  :: c, dxi, dxi0, dyi, dyi0, dzi, dzi0
  REAL*8  :: afac, afac0, Edir(3), Emax, delta_nU
  REAL*8  :: E0d_x, Ed_x, E0d_y, Ed_y, E0d_z, Ed_z

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  mat = 0.d0
  rhs = 0.d0
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta          ! time step conversion factor
  dtfac0 = dt*(1.d0-theta)   ! time step conversion factor
  afac   = adot/a            ! expansion factor (new time)
  afac0  = adot0/a0          ! expansion factor (old time)
  c      = c_light           ! speed of light [cm/s]
  dxi    = a/dx/lUn
  dyi    = a/dy/lUn
  dzi    = a/dz/lUn
  dxi0   = a0/dx/lUn0
  dyi0   = a0/dy/lUn0
  dzi0   = a0/dz/lUn0
  eps    = 1.d-12
  mu     = 1.d-2/(dxi+dyi+dzi) ! artificial viscosity (multiplied by c)
  delta_nU = nUn0/nUn

  ! iterate over the active domain
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           ! compute direction of radiation propagation, 
           ! correcting for case of zero gradients (need a
           ! a nonzero in each direction for MG solver)
           Edir = (/ (E0(i+1,j,k) - E0(i-1,j,k))*0.5d0, &
                     (E0(i,j+1,k) - E0(i,j-1,k))*0.5d0, &
                     (E0(i,j,k+1) - E0(i,j,k-1))*0.5d0 /) / abs(E0(i,j,k))
           Emax = max(maxval(abs(Edir)), eps*1.0e-3)
           Edir(1) = sign(max(abs(Edir(1)), eps*Emax), Edir(1))
           Edir(2) = sign(max(abs(Edir(2)), eps*Emax), Edir(2))
           Edir(3) = sign(max(abs(Edir(3)), eps*Emax), Edir(3))
           Edir = Edir / sqrt(Edir(1)**2 + Edir(2)**2 + Edir(3)**2)

           ! compute the cell-centered opacity
           if (kappa_h2on == 1) then
              kap = kappa(i,j,k)
           else
              kap = kappa_c
           endif
           kap0 = kap * delta_nU

           ! initialize matrix/rhs at this point
           mat(1,i,j,k) = -dtfac*c*mu*dzi*dzi
           mat(2,i,j,k) = -dtfac*c*mu*dyi*dyi
           mat(3,i,j,k) = -dtfac*c*mu*dxi*dxi
           mat(4,i,j,k) = 1.d0 + dtfac*(afac + kap + c*mu*2.d0* &
                               (dxi*dxi + dyi*dyi + dzi*dzi))
           mat(5,i,j,k) = -dtfac*c*mu*dxi*dxi
           mat(6,i,j,k) = -dtfac*c*mu*dyi*dyi
           mat(7,i,j,k) = -dtfac*c*mu*dzi*dzi
           rhs(i,j,k) = (dtfac/rUn + dtfac0/rUn0)*eta(i,j,k)       &
                        + (1.d0 - dtfac0*(afac0+kap0))*E0(i,j,k)   &
                        - (1.d0 + dtfac*(afac+kap))*E(i,j,k)

           !--------------
           ! z-directional differencing
           !    compute direction of radiation propagation
           !    and gradients of old E and guess
           if (Edir(3) < 0.d0) then     ! left
              E0d_z = Edir(3)*(E0(i,j,k) - E0(i,j,k-1))*dzi0 
              Ed_z  = Edir(3)*(E(i,j,k)  - E(i,j,k-1))*dzi 
              mat(1,i,j,k) = mat(1,i,j,k) + Edir(3)*dtfac*c*dzi
              mat(4,i,j,k) = mat(4,i,j,k) - Edir(3)*dtfac*c*dzi
           else                         ! right
              E0d_z = Edir(3)*(E0(i,j,k+1) - E0(i,j,k))*dzi0
              Ed_z  = Edir(3)*(E(i,j,k+1)  - E(i,j,k))*dzi
              mat(7,i,j,k) = mat(7,i,j,k) - Edir(3)*dtfac*c*dzi
              mat(4,i,j,k) = mat(4,i,j,k) + Edir(3)*dtfac*c*dzi
           endif
           E0d_z = E0d_z + mu*(E0(i,j,k+1) - 2.d0*E0(i,j,k) + E0(i,j,k-1))*dzi0*dzi0
           Ed_z  = Ed_z  + mu*(E(i,j,k+1) - 2.d0*E(i,j,k) + E(i,j,k-1))*dzi*dzi
           rhs(i,j,k) = rhs(i,j,k) + c*dtfac0*E0d_z + c*dtfac*Ed_z

           !--------------
           ! y-directional differencing
           !    compute direction of radiation propagation
           !    and gradients of old E and guess
           if (Edir(2) < 0.d0) then     ! left
              E0d_y = Edir(2)*(E0(i,j,k) - E0(i,j-1,k))*dyi0
              Ed_y  = Edir(2)*(E(i,j,k)  - E(i,j-1,k))*dyi
              mat(2,i,j,k) = mat(2,i,j,k) + Edir(2)*dtfac*c*dyi
              mat(4,i,j,k) = mat(4,i,j,k) - Edir(2)*dtfac*c*dyi
           else                         ! right
              E0d_y = Edir(2)*(E0(i,j+1,k) - E0(i,j,k))*dyi0
              Ed_y  = Edir(2)*(E(i,j+1,k)  - E(i,j,k))*dyi
              mat(6,i,j,k) = mat(6,i,j,k) - Edir(2)*dtfac*c*dyi
              mat(4,i,j,k) = mat(4,i,j,k) + Edir(2)*dtfac*c*dyi
           endif
           E0d_y = E0d_y + mu*(E0(i,j+1,k) - 2.d0*E0(i,j,k) + E0(i,j-1,k))*dyi0*dyi0
           Ed_y  = Ed_y  + mu*(E(i,j+1,k) - 2.d0*E(i,j,k) + E(i,j-1,k))*dyi*dyi
           rhs(i,j,k) = rhs(i,j,k) + c*dtfac0*E0d_y + c*dtfac*Ed_y

           !--------------
           ! x-directional differencing
           !    compute direction of radiation propagation
           !    and gradients of old E and guess
           if (Edir(1) < 0.d0) then     ! left
              E0d_x = Edir(1)*(E0(i,j,k) - E0(i-1,j,k))*dxi0
              Ed_x  = Edir(1)*(E(i,j,k)  - E(i-1,j,k))*dxi
              mat(3,i,j,k) = mat(3,i,j,k) + Edir(1)*dtfac*c*dxi
              mat(4,i,j,k) = mat(4,i,j,k) - Edir(1)*dtfac*c*dxi
           else                         ! right
              E0d_x = Edir(1)*(E0(i+1,j,k) - E0(i,j,k))*dxi0
              Ed_x  = Edir(1)*(E(i+1,j,k)  - E(i,j,k))*dxi
              mat(5,i,j,k) = mat(5,i,j,k) - Edir(1)*dtfac*c*dxi
              mat(4,i,j,k) = mat(4,i,j,k) + Edir(1)*dtfac*c*dxi
           endif
           E0d_x = E0d_x + mu*(E0(i+1,j,k) - 2.d0*E0(i,j,k) + E0(i-1,j,k))*dxi0*dxi0
           Ed_x  = Ed_x  + mu*(E(i+1,j,k) - 2.d0*E(i,j,k) + E(i-1,j,k))*dxi*dxi
           rhs(i,j,k) = rhs(i,j,k) + c*dtfac0*E0d_x + c*dtfac*Ed_x

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
              mat(:,i,j,k) = 0.d0
              mat(4,i,j,k) = 1.d0
              rhs(i,j,k)   = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCZl==2) then
        k = x2s
        do j=1,Ny
           do i=1,Nx
              mat(4,i,j,k) = mat(4,i,j,k) + mat(1,i,j,k)
              mat(1,i,j,k) = 0.d0
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
              mat(:,i,j,k) = 0.d0
              mat(4,i,j,k) = 1.d0
              rhs(i,j,k)   = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCYl==2) then
        j = x1s
        do k=1,Nz
           do i=1,Nx
              mat(4,i,j,k) = mat(4,i,j,k) + mat(2,i,j,k)
              mat(2,i,j,k) = 0.d0
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
              mat(:,i,j,k) = 0.d0
              mat(4,i,j,k) = 1.d0
              rhs(i,j,k)   = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        do k=1,Nz
           do j=1,Ny
              mat(4,i,j,k) = mat(4,i,j,k) + mat(3,i,j,k)
              mat(3,i,j,k) = 0.d0
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
              mat(:,i,j,k) = 0.d0
              mat(4,i,j,k) = 1.d0
              rhs(i,j,k)   = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        do k=1,Nz
           do j=1,Ny
              mat(4,i,j,k) = mat(4,i,j,k) + mat(5,i,j,k)
              mat(5,i,j,k) = 0.d0
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
              mat(:,i,j,k) = 0.d0
              mat(4,i,j,k) = 1.d0
              rhs(i,j,k)   = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCYr==2) then
        j = x1e
        do k=1,Nz
           do i=1,Nx
              mat(4,i,j,k) = mat(4,i,j,k) + mat(6,i,j,k)
              mat(6,i,j,k) = 0.d0
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
              mat(:,i,j,k) = 0.d0
              mat(4,i,j,k) = 1.d0
              rhs(i,j,k)   = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCZr==2) then
        k = x2e
        do j=1,Ny
           do i=1,Nx
              mat(4,i,j,k) = mat(4,i,j,k) + mat(7,i,j,k)
              mat(7,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  rhsnorm = sum(rhs*rhs)

  return
end subroutine FSProb_SetupSystem3D
!=======================================================================






subroutine FSProb_SetupSystem2D(mat, rhs, rhsnorm, E, E0, kappa_h2on,    &
     kappa, kappa_c, eta, dt, a, a0, adot, adot0, theta, lUn, lUn0, rUn, &
     rUn0, nUn, nUn0, dx, dy, BCXl, BCXr, BCYl, BCYr, x0s, x0e, x1s,     &
     x1e, Nx, Ny, NGxl, NGxr, NGyl, NGyr, xlface, xrface, ylface, yrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       January 2010
!  modified:   
!
!  PURPOSE: 2D version of the routine
!=======================================================================
  implicit none
#include "fortran_types.def"
  
  !--------------
  ! argument declarations
  INTG_PREC, intent(in) :: kappa_h2on
  INTG_PREC, intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  INTG_PREC, intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  P_PREC, intent(in) :: a, a0, adot, adot0
  R_PREC,    intent(in) :: kappa_c, dt, theta, dx, dy
  R_PREC,    intent(in) :: lUn, lUn0, rUn, rUn0, nUn, nUn0
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), intent(in) :: E, E0, eta, kappa
  REAL*8,  intent(out) :: mat(5,x0s:x0e,x1s:x1e)
  REAL*8,  intent(out) :: rhs(x0s:x0e,x1s:x1e)
  R_PREC,    intent(out) :: rhsnorm
  INTG_PREC, intent(out) :: ier

  !--------------
  ! locals
  INTG_PREC :: i, j
  REAL*8  :: dtfac, dtfac0, kap, kap0, eps, mu
  REAL*8  :: c, dxi, dxi0, dyi, dyi0
  REAL*8  :: afac, afac0, Edir(2), Emax, delta_nU
  REAL*8  :: E0d_x, Ed_x, E0d_y, Ed_y

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  mat = 0.d0
  rhs = 0.d0
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta          ! time step conversion factor
  dtfac0 = dt*(1.d0-theta)   ! time step conversion factor
  afac   = adot/a            ! expansion factor (new time)
  afac0  = adot0/a0          ! expansion factor (old time)
  c      = c_light           ! speed of light [cm/s]
  dxi    = a/dx/lUn
  dyi    = a/dy/lUn
  dxi0   = a0/dx/lUn0
  dyi0   = a0/dy/lUn0
  eps    = 1.d-12
  mu     = 1.d-2/(dxi+dyi)   ! artificial viscosity (multiplied by c)
  delta_nU = nUn0/nUn

  ! iterate over the active domain
  do j=1,Ny,1
     do i=1,Nx,1
        
        ! compute direction of radiation propagation, 
        ! correcting for case of zero gradients (need at least 
        ! a bit of diffusion in each direction for MG solver)
        Edir = (/ (E0(i+1,j) - E0(i-1,j))*0.5d0*dxi, &
                  (E0(i,j+1) - E0(i,j-1))*0.5d0*dyi /) / abs(E0(i,j))
        Emax = max(maxval(abs(Edir)), eps*1.0e-3)
        Edir(1) = sign(max(abs(Edir(1)), eps*Emax), Edir(1))
        Edir(2) = sign(max(abs(Edir(2)), eps*Emax), Edir(2))
        Edir = Edir / sqrt(Edir(1)**2 + Edir(2)**2)

        ! compute the cell-centered opacity
        if (kappa_h2on == 1) then
           kap = kappa(i,j)
        else
           kap = kappa_c
        endif
        kap0 = kap * delta_nU
        
        ! initialize matrix/rhs at this point
        mat(1,i,j) = -dtfac*c*mu*dyi*dyi
        mat(2,i,j) = -dtfac*c*mu*dxi*dxi
        mat(3,i,j) = 1.d0 + dtfac*(afac + kap  &
             + c*mu*2.d0*(dxi*dxi + dyi*dyi))
        mat(4,i,j) = -dtfac*c*mu*dxi*dxi
        mat(5,i,j) = -dtfac*c*mu*dyi*dyi
        rhs(i,j) = (dtfac/rUn + dtfac0/rUn0)*eta(i,j)     &
                 + (1.d0 - dtfac0*(afac0+kap0))*E0(i,j)   &
                 - (1.d0 + dtfac*(afac+kap))*E(i,j)

        !--------------
        ! y-directional differencing
        !    compute direction of radiation propagation
        !    and gradients of old E and guess
        if (Edir(2) < 0.d0) then     ! left
           E0d_y = Edir(2)*(E0(i,j) - E0(i,j-1))*dyi0
           Ed_y  = Edir(2)*(E(i,j)  - E(i,j-1))*dyi
           mat(1,i,j) = mat(1,i,j) + Edir(2)*dtfac*c*dyi
           mat(3,i,j) = mat(3,i,j) - Edir(2)*dtfac*c*dyi
        else                         ! right
           E0d_y = Edir(2)*(E0(i,j+1) - E0(i,j))*dyi0
           Ed_y  = Edir(2)*(E(i,j+1)  - E(i,j))*dyi
           mat(5,i,j) = mat(5,i,j) - Edir(2)*dtfac*c*dyi
           mat(3,i,j) = mat(3,i,j) + Edir(2)*dtfac*c*dyi
        endif
        E0d_y = E0d_y + mu*(E0(i,j+1) - 2.d0*E0(i,j) + E0(i,j-1))*dyi0*dyi0
        Ed_y  = Ed_y  + mu*(E(i,j+1) - 2.d0*E(i,j) + E(i,j-1))*dyi*dyi
        rhs(i,j) = rhs(i,j) + c*dtfac0*E0d_y + c*dtfac*Ed_y
        
        !--------------
        ! x-directional differencing
        !    compute direction of radiation propagation
        !    and gradients of old E and guess
        if (Edir(1) < 0.d0) then     ! left
           E0d_x = Edir(1)*(E0(i,j) - E0(i-1,j))*dxi0
           Ed_x  = Edir(1)*(E(i,j)  - E(i-1,j))*dxi
           mat(2,i,j) = mat(2,i,j) + Edir(1)*dtfac*c*dxi
           mat(3,i,j) = mat(3,i,j) - Edir(1)*dtfac*c*dxi
        else                         ! right
           E0d_x = Edir(1)*(E0(i+1,j) - E0(i,j))*dxi0
           Ed_x  = Edir(1)*(E(i+1,j)  - E(i,j))*dxi
           mat(4,i,j) = mat(4,i,j) - Edir(1)*dtfac*c*dxi
           mat(3,i,j) = mat(3,i,j) + Edir(1)*dtfac*c*dxi
        endif
        E0d_x = E0d_x + mu*(E0(i+1,j) - 2.d0*E0(i,j) + E0(i-1,j))*dxi0*dxi0
        Ed_x  = Ed_x  + mu*(E(i+1,j) - 2.d0*E(i,j) + E(i-1,j))*dxi*dxi
        rhs(i,j) = rhs(i,j) + c*dtfac0*E0d_x + c*dtfac*Ed_x

     enddo
  enddo

  ! update matrix/rhs based on boundary conditions/location
  !    y-left face
  if (ylface == 1) then
     ! Dirichlet
     if (BCYl==1) then
        j = x1s
        do i=1,Nx
           mat(:,i,j) = 0.d0
           mat(3,i,j) = 1.d0
           rhs(i,j)   = 0.d0
        enddo
     ! Neumann
     else if (BCYl==2) then
        j = x1s
        do i=1,Nx
           mat(3,i,j) = mat(3,i,j) + mat(1,i,j)
           mat(1,i,j) = 0.d0
        enddo
     endif
  end if

  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        do j=1,Ny
           mat(:,i,j) = 0.d0
           mat(3,i,j) = 1.d0
           rhs(i,j)   = 0.d0
        enddo
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        do j=1,Ny
           mat(3,i,j) = mat(3,i,j) + mat(2,i,j)
           mat(2,i,j) = 0.d0
        enddo
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        do j=1,Ny
           mat(:,i,j) = 0.d0
           mat(3,i,j) = 1.d0
           rhs(i,j)   = 0.d0
        enddo
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        do j=1,Ny
           mat(3,i,j) = mat(3,i,j) + mat(4,i,j)
           mat(4,i,j) = 0.d0
        enddo
     endif
  endif

  !    y-right face
  if (yrface==1) then
     ! Dirichlet
     if (BCYr==1) then
        j = x1e
        do i=1,Nx
           mat(:,i,j) = 0.d0
           mat(3,i,j) = 1.d0
           rhs(i,j)   = 0.d0
        enddo
     ! Neumann
     else if (BCYr==2) then
        j = x1e
        do i=1,Nx
           mat(3,i,j) = mat(3,i,j) + mat(5,i,j)
           mat(5,i,j) = 0.d0
        enddo
     endif
  endif

  rhsnorm = sum(rhs*rhs)

  return
end subroutine FSProb_SetupSystem2D
!=======================================================================






subroutine FSProb_SetupSystem1D(mat, rhs, rhsnorm, E, E0, kappa_h2on,    &
     kappa, kappa_c, eta, dt, a, a0, adot, adot0, theta, lUn, lUn0, rUn, &
     rUn0, nUn, nUn0, dx, BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface,  &
     xrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       January 2010
!  modified:   
!
!  PURPOSE: 1D version of the routine
!=======================================================================
  implicit none
#include "fortran_types.def"
  
  !--------------
  ! argument declarations
  INTG_PREC, intent(in) :: kappa_h2on
  INTG_PREC, intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  P_PREC, intent(in) :: a, a0, adot, adot0
  R_PREC,    intent(in) :: kappa_c, dt, theta, dx
  R_PREC,    intent(in) :: lUn, lUn0, rUn, rUn0, nUn, nUn0
  R_PREC, dimension(1-NGxl:Nx+NGxr), intent(in) :: E, E0, eta, kappa
  REAL*8,  intent(out) :: mat(3,x0s:x0e)
  REAL*8,  intent(out) :: rhs(x0s:x0e)
  R_PREC,    intent(out) :: rhsnorm
  INTG_PREC, intent(out) :: ier

  !--------------
  ! locals
  INTG_PREC :: i
  REAL*8  :: dtfac, dtfac0, kap, kap0, eps, mu
  REAL*8  :: c, dxi, dxi0
  REAL*8  :: afac, afac0, Edir, delta_nU, ONE
  REAL*8  :: E0d_x, Ed_x

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  mat = 0.d0
  rhs = 0.d0
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta          ! time step conversion factor
  dtfac0 = dt*(1.d0-theta)   ! time step conversion factor
  afac   = adot/a            ! expansion factor (new time)
  afac0  = adot0/a0          ! expansion factor (old time)
  c      = c_light           ! speed of light [cm/s]
  dxi    = a/dx/lUn
  dxi0   = a0/dx/lUn0
  eps    = 1.d-12
  mu     = 1.d-2/dxi         ! artificial viscosity (multiplied by c)
  delta_nU = nUn0/nUn
  ONE    = 1.d0

  ! iterate over the active domain
  do i=1,Nx,1
        
     ! compute direction of radiation propagation, 
     ! correcting for case of zero gradients (need at least 
     ! a bit of diffusion in each direction for MG solver)
     Edir = (E0(i+1) - E0(i-1))*0.5d0*dxi / abs(E0(i))
     Edir = sign(ONE, Edir)

     ! compute the cell-centered opacity
     if (kappa_h2on == 1) then
        kap = kappa(i)
     else
        kap = kappa_c
     endif
     kap0 = kap * delta_nU
     
     ! initialize matrix/rhs at this point
     mat(1,i) = -dtfac*c*mu*dxi*dxi
     mat(2,i) = 1.d0 + dtfac*(afac + kap   &
          + c*mu*2.d0*dxi*dxi)
     mat(3,i) = -dtfac*c*mu*dxi*dxi
     rhs(i) = (dtfac/rUn + dtfac0/rUn0)*eta(i)     &
              + (1.d0 - dtfac0*(afac0+kap0))*E0(i) &
              - (1.d0 + dtfac*(afac+kap))*E(i)

     !--------------
     ! x-directional differencing
     !    compute direction of radiation propagation
     !    and gradients of old E and guess
     if (Edir < 0.d0) then     ! left
        E0d_x = Edir*(E0(i) - E0(i-1))*dxi0
        Ed_x  = Edir*(E(i)  - E(i-1))*dxi
        mat(1,i) = mat(1,i) + Edir*dtfac*c*dxi
        mat(2,i) = mat(2,i) - Edir*dtfac*c*dxi
     else                         ! right
        E0d_x = Edir*(E0(i+1) - E0(i))*dxi0
        Ed_x  = Edir*(E(i+1)  - E(i))*dxi
        mat(3,i) = mat(3,i) - Edir*dtfac*c*dxi
        mat(2,i) = mat(2,i) + Edir*dtfac*c*dxi
     endif
     E0d_x = E0d_x + mu*(E0(i+1) - 2.d0*E0(i) + E0(i-1))*dxi0*dxi0
     Ed_x  = Ed_x  + mu*(E(i+1) - 2.d0*E(i) + E(i-1))*dxi*dxi
     rhs(i) = rhs(i) + c*dtfac0*E0d_x + c*dtfac*Ed_x

  enddo

  ! update matrix/rhs based on boundary conditions/location
  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        mat(:,i) = 0.d0
        mat(2,i) = 1.d0
        rhs(i)   = 0.d0
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        mat(2,i) = mat(2,i) + mat(1,i)
        mat(1,i) = 0.d0
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        mat(:,i) = 0.d0
        mat(2,i) = 1.d0
        rhs(i)   = 0.d0
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        mat(2,i) = mat(2,i) + mat(3,i)
        mat(3,i) = 0.d0
     endif
  endif

  rhsnorm = sum(rhs*rhs)

  return
end subroutine FSProb_SetupSystem1D
!=======================================================================

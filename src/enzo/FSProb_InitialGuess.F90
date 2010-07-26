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
subroutine FSProb_InitialGuess(Ef, Ef0, eta, iguess, dt, kappa_h2on,   &
     kappa, kappa_c, a, adot, aUn, lUn, tUn, EUn, dUn, dx, dy, dz, Nx, &
     Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June, 2009
!
!  PURPOSE: Computes an initial guess at the time evolved solution 
!           for the free-streaming radiation equation.
!
!  INPUTS:
!     Ef0        - radiation energy density (old step)
!     eta        - emissivity sources over domain
!     iguess     - method to use in constructing initial guess
!     dt         - time step size
!     kappa_h2on - use spatially-dependent opacity (1=on, 0=off)
!     kappa      - spatially-dependent opacity
!     kappa0     - background opacity
!     a          - cosmological expansion parameter
!     adot       - d(a)/dt
!     *Un        - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     Ef         - guess at time-evolved solution
!     ier        - success/failure flag (1->success, 0->failure)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in) :: iguess, kappa_h2on
  integer, intent(in) :: Nx, NGxl, NGxr
  integer, intent(in) :: Ny, NGyl, NGyr
  integer, intent(in) :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in) :: a, adot
  real,    intent(in) :: dt, kappa_c, dx, dy, dz
  real,    intent(in) :: aUn, lUn, tUn, EUn, dUn
  real,    intent(in) :: eta(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  real,    intent(in) :: Ef0(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  real,    intent(in) :: kappa(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  real,    intent(out) :: Ef(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  
  !--------------
  ! locals
  integer :: i, j, k
  real :: P, c, pi

  !=======================================================================

  ! initialize output flag to success
  ier = 1

  ! set some constants
  c  = 2.99792458d10           ! speed of light [cm/s]
  pi = 4.D0*datan(1.D0)

  ! compute initial guess based on input parameter
  !    use analytical solution locally in space
  if (iguess == 1) then

     ! array-based opacity
     if (kappa_h2on == 1) then
        do k=1,Nz
           do j=1,Ny
              do i=1,Nx
                 ! if attenuation is negligible, use simpler analytical 
                 ! solution to avoid division by zero 
                 P = adot/a + c*kappa(i,j,k)
                 if (P < 1.d-14) then
                    Ef(i,j,k) = Ef0(i,j,k) + dt*eta(i,j,k)*4.d0*pi
                 ! solution to avoid flt pt overflow
                 elseif (P*dt > 7.0d2) then
                    Ef(i,j,k) = eta(i,j,k)*4.d0*pi/P
                 ! otherwise use full analytical solution
                 else
                    Ef(i,j,k) = (Ef0(i,j,k) - eta(i,j,k)*4.d0*pi/P)*exp(-P*dt) &
                              + eta(i,j,k)*4.d0*pi/P
                 endif
              enddo
           enddo
        enddo

     ! constant opacity
     else
        do k=1,Nz
           do j=1,Ny
              do i=1,Nx
                 ! if attenuation is negligible, use simpler analytical 
                 ! solution to avoid division by zero 
                 P = adot/a + c*kappa_c
                 if (P < 1.d-14) then
                    Ef(i,j,k) = Ef0(i,j,k) + dt*eta(i,j,k)*4.d0*pi
                 ! solution to avoid flt pt overflow
                 elseif (P*dt > 7.0d2) then
                    Ef(i,j,k) = eta(i,j,k)*4.d0*pi/P
                 ! otherwise use full analytical solution
                 else
                    Ef(i,j,k) = (Ef0(i,j,k) - eta(i,j,k)*4.d0*pi/P)*exp(-P*dt) &
                              + eta(i,j,k)*4.d0*pi/P
                 endif
              enddo
           enddo
        enddo
     endif

  !    otherwise use old-time solution as initial guess
  else 
     Ef = Ef0
  endif

  return
end subroutine FSProb_InitialGuess
!=======================================================================

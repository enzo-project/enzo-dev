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
#include "fortran_types.def"

  !--------------
  ! argument declarations
  INTG_PREC, intent(in) :: iguess, kappa_h2on
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in) :: a, adot
  R_PREC,    intent(in) :: dt, kappa_c, dx, dy, dz
  R_PREC,    intent(in) :: aUn, lUn, tUn, EUn, dUn
  R_PREC,    intent(in) :: eta(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  R_PREC,    intent(in) :: Ef0(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  R_PREC,    intent(in) :: kappa(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  R_PREC,    intent(out) :: Ef(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  
  !--------------
  ! locals
  INTG_PREC :: i, j, k
  R_PREC :: P, c, pi

  !=======================================================================

  ! initialize output flag to success
  ier = 1

  ! set some constants
  c  = c_light             ! speed of light [cm/s]
  pi = pi_val

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
                 if (P < 1.e-14_RKIND) then
                    Ef(i,j,k) = Ef0(i,j,k) + dt*eta(i,j,k)*4._RKIND*pi
                 ! solution to avoid flt pt overflow
                 elseif (P*dt > 7.e2_RKIND) then
                    Ef(i,j,k) = eta(i,j,k)*4._RKIND*pi/P
                 ! otherwise use full analytical solution
                 else
                    Ef(i,j,k) = (Ef0(i,j,k) - eta(i,j,k)*4._RKIND*pi/P)*exp(-P*dt) &
                              + eta(i,j,k)*4._RKIND*pi/P
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
                 if (P < 1.e-14_RKIND) then
                    Ef(i,j,k) = Ef0(i,j,k) + dt*eta(i,j,k)*4._RKIND*pi
                 ! solution to avoid flt pt overflow
                 elseif (P*dt > 7.e2_RKIND) then
                    Ef(i,j,k) = eta(i,j,k)*4._RKIND*pi/P
                 ! otherwise use full analytical solution
                 else
                    Ef(i,j,k) = (Ef0(i,j,k) - eta(i,j,k)*4._RKIND*pi/P)*exp(-P*dt) &
                              + eta(i,j,k)*4._RKIND*pi/P
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

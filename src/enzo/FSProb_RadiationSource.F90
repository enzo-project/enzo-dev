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
subroutine FSProb_RadiationSource(eta, time, a, ProbType, NGammaDot,   &
     EtaRadius, EtaCenter, aUn, lUn, tUn, rUn, Nx, Ny, Nz, NGxl, NGxr, &
     NGyl, NGyr, NGzl, NGzr, x0L, x0R, x1L, x1R, x2L, x2R, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June, 2009
!
!  PURPOSE: Computes the emissivity point sources for the 
!           free-streaming radiation equation,
!              Div(Df(Ef)*Grad(Ef)) - adot/a*Ef + eta
!           Note: we use eta alone in the equation (no 4*pi term), so 
!           it must be scaled appropriately.
!
!  INPUTS:
!     time       - simulation time for evaluation
!     a          - cosmological expansion parameter
!     ProbType   - flag denoting physical problem to run
!     NGammaDot  - ionization source
!     EtaRadius  - ionization source radius in cells
!     EtaCenter  - ionization source center (comoving, 3D coordinates in cm)
!     *Un        - variable scaling constants
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     x*L/x*R    - left/right subdomain boundaries (comoving, no ghosts)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     eta        - array of emissivity sources
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
  INTG_PREC, intent(in) :: ProbType
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  REAL*8,  intent(in) :: NGammaDot
  P_PREC, intent(in) :: a
  R_PREC,    intent(in) :: time, EtaRadius, EtaCenter(3)
  R_PREC,    intent(in) :: aUn, lUn, tUn, rUn
  R_PREC,    intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  R_PREC,    intent(out) :: eta(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  
!--------------
! locals
  INTEGER   :: seed(12)
  INTG_PREC :: i, j, k, l, nsrc
  R_PREC    :: h_nu0, etaconst, rnums(10)
  R_PREC    :: dx, dy, dz, cellXl, cellXr, cellYl, cellYr, cellZl, cellZr
  R_PREC    :: cellXc, cellYc, cellZc
  REAL*8    :: dV

!=======================================================================

  ! initialize output to have all zero values, flag to success
  eta = 0._RKIND
  ier = 1
  
  ! initialize constants
  dx    = (x0R-x0L)/Nx                ! mesh spacing (comoving), x0 direction
  dy    = (x1R-x1L)/Ny                ! mesh spacing (comoving), x1 direction
  dz    = (x2R-x2L)/Nz                ! mesh spacing (comoving), x2 direction
  dV    = dx*dy*dz*(dble(lUn))**3     ! cell volume (proper)
  h_nu0 = 13.6_RKIND*ev2erg               ! ionization energy of HI [ergs]

  ! compute point source emissivity for various problems

  !   multiple (random) emissivity sources (w/ mean NGammaDot photons/s/cm^2)
  if (ProbType == 450) then

     ! get the number of sources from EtaRadius
     nsrc = EtaRadius   ! cast to an INTG_PREC

     ! set the seed by casting the time to an INTG_PREC
     seed(1) = 5._RKIND*time/tUn + x0L/dx + x1L/dy + x2L/dz
     seed(1) = seed(1) + 13
!     print *,'random seed = ',seed(1),' time = ',time
     call random_seed(PUT=seed)

     ! loop over sources
     do l=1,nsrc

        ! get 4 random numbers for each source (3 location, 1 strength)
        call random_number(rnums)
        i = max(min(int(rnums(4)*Nx,IKIND), Nx-1), 2)
        j = max(min(int(rnums(6)*Ny,IKIND), Ny-1), 2)
        k = max(min(int(rnums(8)*Nz,IKIND), Nz-1), 2)
!        eta(i,j,k) = rnums(10)*h_nu0*REAL(NGammaDot/dV,RKIND)
        eta(i,j,k) = rnums(10)*h_nu0*NGammaDot/dV
!        print '(A,3(i2,1x),A,es9.2)', '   setting source at ',i,j,k,' with strength ',eta(i,j,k)

     enddo

  !   point-source emissivity at location EtaCenter(1:3)
  else if ((ProbType > 450) .and. (ProbType <= 460)) then

     ! one-cell source
     if (EtaRadius == 0._RKIND) then
        
        ! compute eta factor for given ionization source
        etaconst = h_nu0*NGammaDot/dV
        
        ! place ionization source in one cell
        do k=1,Nz,1
           
           ! z-boundaries (comoving) for this cell
           cellZl = x2L + (k-1)*dz
           cellZr = cellZl + dz
           
           do j=1,Ny,1
              
              ! y-boundaries (comoving) for this cell
              cellYl = x1L + (j-1)*dy
              cellYr = cellYl + dy
              
              do i=1,Nx,1
                 
                 ! x-boundaries (comoving) for this cell
                 cellXl = x0L + (i-1)*dx
                 cellXr = cellXl + dx
                 
                 ! see if domain center is in cell (or on left edge)
                 if ( (cellXl <= EtaCenter(1)) .and. (cellXr > EtaCenter(1)) .and. &
                      (cellYl <= EtaCenter(2)) .and. (cellYr > EtaCenter(2)) .and. &
                      (cellZl <= EtaCenter(3)) .and. (cellZr > EtaCenter(3)) ) then
                    eta(i,j,k) = etaconst
                 endif
                 
              enddo
           enddo
        enddo

     else

        ! compute eta factor for given ionization source
        etaconst = h_nu0*NGammaDot/dV/8._RKIND/(EtaRadius**3)
        
        ! place ionization source in center of domain
        do k=1,Nz,1
           
           ! z-center (comoving) for this cell
           cellZc = x2L + (k-0.5_RKIND)*dz
           
           do j=1,Ny,1
              
              ! y-center (comoving) for this cell
              cellYc = x1L + (j-0.5_RKIND)*dy
              
              do i=1,Nx,1
                 
                 ! x-center (comoving) for this cell
                 cellXc = x0L + (i-0.5_RKIND)*dx
                 
                 ! see if cell is within source region
                 if ( (abs(cellXc-EtaCenter(1)) < EtaRadius*dx) .and. &
                      (abs(cellYc-EtaCenter(2)) < EtaRadius*dy) .and. &
                      (abs(cellZc-EtaCenter(3)) < EtaRadius*dz) ) then
                    eta(i,j,k) = etaconst
                 endif
                 
              enddo
           enddo
        enddo

     endif ! EtaRadius == 0

  !   point-source emissivity at center of every processor
  elseif (ProbType == 452) then

     ! compute eta factor for given ionization source
     etaconst = h_nu0*NGammaDot/dV
        
     ! place ionization source in center of subdomain
     eta(int(Nx/2,IKIND),int(Ny/2,IKIND),int(Nz/2,IKIND)) = etaconst

  !   homogeneous emissivity field w/ strength h_nu0*NGammaDot/dV
  elseif (ProbType == 462) then

     ! compute eta factor for given ionization source
     eta = h_nu0*NGammaDot/dV
        
  endif ! ProbType

  return
end subroutine FSProb_RadiationSource
!=======================================================================

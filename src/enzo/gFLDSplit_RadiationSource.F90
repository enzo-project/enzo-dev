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
subroutine gFLDSplit_RadiationSource(Ersrc, time, a, ProbType, ESpectrum, &
     NGammaDot, EtaRadius, EtaCenter, aUnits, LenUnits, TimeUnits,        &
     ErUnits, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, x0L, x0R,   &
     x1L, x1R, x2L, x2R, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!
!  PURPOSE: Computes the emissivity point sources into the radiation 
!           energy equation
!
!  INPUTS:
!     time       - simulation time for evaluation
!     a          - cosmological expansion parameter
!     ProbType   - flag denoting physical problem to run
!     Nchem      - number of chemical species
!     HFrac      - percentage of mass composed of Hydrogen
!     ESpectrum  - radiation spectrum choice
!                       1 -> 1e5 black body spectrum
!                       0 -> power law spectrum
!                      -1 -> monochromatic 
!     NGammaDot  - ionization source
!     EtaRadius  - ionization source radius in cells
!     EtaCenter  - ionization source center (comoving, 3D coordinates in cm)
!     *Units     - variable scaling constants
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     x*L/x*R    - left/right subdomain boundaries (comoving, no ghosts)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     Ersrc      - array of emissivity sources
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
  INTG_PREC, intent(in) :: ESpectrum, ProbType
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in) :: a
  R_PREC,    intent(in) :: time, NGammaDot, EtaRadius
  R_PREC,    intent(in) :: EtaCenter(3)
  R_PREC,    intent(in) :: aUnits, LenUnits, TimeUnits, ErUnits
  R_PREC,    intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  R_PREC,    intent(out) :: Ersrc(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  
  !--------------
  ! locals
  INTG_PREC :: i, j, k
  R_PREC :: pi, h_nu0, etaconst, specconst
  R_PREC :: dx, dy, dz, dV, cellXl, cellXr, cellYl, cellYr, cellZl, cellZr
  R_PREC :: cellXc, cellYc, cellZc

!=======================================================================

!!$  write(*,*) 'Entering gFLDSplit::RadiationSource routine'
!!$  write(*,*) '   NGammaDot = ',NGammaDot
!!$  write(*,*) '   EtaRadius = ',EtaRadius
!!$  write(*,*) '   EtaCenter = ',EtaCenter

  ! initialize output to have all zero values, flag to success
  Ersrc = 0._RKIND
  ier = 1
  
  ! initialize constants
  pi    = pi_val
  dx    = (x0R-x0L)/Nx                ! mesh spacing (comoving), x0 direction
  dy    = (x1R-x1L)/Ny                ! mesh spacing (comoving), x1 direction
  dz    = (x2R-x2L)/Nz                ! mesh spacing (comoving), x2 direction
  dV    = dx*dy*dz*(LenUnits)**3      ! cell volume (proper)
  h_nu0 = 13.6_RKIND*ev2erg               ! ionization energy of HI [ergs]

  ! scaling factor for T=10^5 blackbody spectrum
  if (ESpectrum == 1) then
     specconst = 1.52877652583602_RKIND
  else
     specconst = 1._RKIND
  endif

  ! compute emissivity for various problems

  !   point-source emissivity at location (EtaCenter(1:3))
  if ((ProbType == 410) .or. (ProbType == 411) .or. &
      (ProbType == 413) .or. (ProbType == 415)) then

     ! one-cell source
     if (EtaRadius == 0._RKIND) then
        
        ! compute eta factor for given ionization source
        etaconst = h_nu0*NGammaDot*specconst/dV
        
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
                    Ersrc(i,j,k) = etaconst
                 endif
                 
              enddo
           enddo
        enddo

     else

        ! compute eta factor for given ionization source
        etaconst = h_nu0*NGammaDot*specconst/dV/8._RKIND/(EtaRadius**3)
        
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
                    Ersrc(i,j,k) = etaconst
                 endif
                 
              enddo
           enddo
        enddo

     endif ! EtaRadius == 0

  !   emissivity flux along x=0 wall (NGammaDot photons/s/cm^2)
  else if (ProbType == 412) then

     ! place ionization source along left wall (if on this subdomain)
     if (x0L == 0._RKIND) then

        ! compute eta factor for given ionization source, and put on wall
        etaconst = h_nu0*NGammaDot*specconst/dy
        do k=1,Nz,1
           do j=1,Ny,1
              Ersrc(1,j,k) = etaconst
           enddo
        enddo
     endif
     
  !   point-source emissivity at center of every processor
  elseif (ProbType == 414) then

     ! compute eta factor for given ionization source
     etaconst = h_nu0*NGammaDot*specconst/dV
        
     ! place ionization source in center of subdomain
     Ersrc(int(Nx/2,IKIND),int(Ny/2,IKIND),int(Nz/2,IKIND)) = etaconst
     
  endif ! ProbType

!!$  write(*,*) 'RadiationSource: individual source is ',etaconst
!!$  etaconst = sum(Ersrc)
!!$  write(*,*) 'RadiationSource: integrated source is ',etaconst
  
  return
end subroutine gFLDSplit_RadiationSource
!=======================================================================

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
subroutine gFLDProblem_RadiationSource(Ersrc, time, Era, eca, nHIa,     &
     nHeIa, nHeIIa, Temp, rhoa, eha, vxa, vya, vza, a, Model, ProbType, &
     Nchem, HFrac, ESpectrum, NGammaDot, EtaRadius, EtaCenter0,         &
     EtaCenter1, EtaCenter2, aUnits, DenUnits, VelUnits, LenUnits,      &
     TimeUnits, ErUnits, ecUnits, NiUnits, Nx, Ny, Nz, NGxl, NGxr,      &
     NGyl, NGyr, NGzl, NGzr, x0L, x0R, x1L, x1R, x2L, x2R, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July, 2007
!
!  PURPOSE: Computes the emissivity point sources into the radiation 
!           energy equation
!
!  INPUTS:
!     time       - simulation time for evaluation
!     Era        - radiation energy density array
!     eca        - specific gas energy correction array
!     nHIa       - Hydrogen I density array
!     nHeIa      - Helium I density array
!     nHeIIa     - Helium II density array
!     Temp       - fluid temperature array
!     rhoa       - fluid density array
!     eha        - specific gas energy array
!     vxa        - x-directional gas velocity array
!     vya        - y-directional gas velocity array
!     vza        - z-directional gas velocity array
!     a          - cosmological expansion parameter
!     Model      - flag denoting physical model to use
!     ProbType   - flag denoting physical problem to run
!     Nchem      - number of chemical species
!     HFrac      - percentage of mass composed of Hydrogen
!     ESpectrum  - radiation spectrum choice
!                       1 -> 1e5 black body spectrum
!                       0 -> power law spectrum
!                      -1 -> monochromatic 
!     NGammaDot  - ionization source
!     EtaRadius  - ionization source radius in cells
!     EtaCenter* - ionization source center (comoving, 3D coordinates in cm)
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
  INTG_PREC, intent(in) :: Nchem, Model, ESpectrum, ProbType
  INTG_PREC, intent(in) :: Nx, NGxl, NGxr
  INTG_PREC, intent(in) :: Ny, NGyl, NGyr
  INTG_PREC, intent(in) :: Nz, NGzl, NGzr
  P_PREC, intent(in) :: a
  R_PREC,    intent(in) :: time, HFrac, NGammaDot, EtaRadius
  R_PREC,    intent(in) :: EtaCenter0, EtaCenter1, EtaCenter2
  R_PREC,    intent(in) :: aUnits, DenUnits, VelUnits, LenUnits, TimeUnits, &
       ErUnits, ecUnits, NiUnits
  R_PREC,    intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  R_PREC, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) &
       :: Era, eca, nHIa, nHeIa, nHeIIa, Temp, rhoa, eha, vxa, vya, vza
  R_PREC, intent(out) :: Ersrc(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  INTG_PREC, intent(out) :: ier
  
!--------------
! locals
  INTG_PREC :: i, j, k
  R_PREC :: pi, h_nu0, etaconst
  R_PREC :: dx, dy, dz, dV, cellXl, cellXr, cellYl, cellYr, cellZl, cellZr
  R_PREC :: cellXc, cellYc, cellZc

!=======================================================================

!!$  write(*,*) 'Entering gFLDProblem::RadiationSource routine'
!!$  write(*,*) '   NGammaDot = ',NGammaDot
!!$  write(*,*) '   EtaRadius = ',EtaRadius

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

  ! compute point source emissivity for various problems

  !   point-source emissivity at location (EtaCenter0,EtaCenter1,EtaCenter2)
  if ((ProbType == 410) .or. (ProbType == 411) .or. &
      (ProbType == 413) .or. (ProbType == 415)) then

     ! one-cell source
     if (EtaRadius == 0._RKIND) then
        
        ! compute eta factor for given ionization source
        etaconst = h_nu0*NGammaDot/dV
        ! [the differing factors arise due to monochromatic vs grey 
        !  radiation energy densities, where the grey source has a 
        !  T=1e5 blackbody spectrum; this factor is essentially 
        !  (\int_{nu0}^{infty} chi(nu) dnu) / (\int_0^{infty} chi(nu)/nu dnu)]
        if (ESpectrum == 1) then
           etaconst = etaconst * 1.52877652583602_RKIND
        endif

        
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
                 if ( (cellXl <= EtaCenter0) .and. (cellXr > EtaCenter0) .and. &
                      (cellYl <= EtaCenter1) .and. (cellYr > EtaCenter1) .and. &
                      (cellZl <= EtaCenter2) .and. (cellZr > EtaCenter2) ) then
                    Ersrc(i,j,k) = etaconst
                 endif
                 
              enddo
           enddo
        enddo

     else

        ! compute eta factor for given ionization source
        etaconst = h_nu0*NGammaDot/dV/8._RKIND/(EtaRadius**3)
        ! (the differing factors arise due to monochromatic vs grey 
        !  radiation energy densities, where the grey source has a 
        !  T=1e5 blackbody spectrum)
        if (ESpectrum == 1) then
           etaconst = etaconst * 1.52877652583602_RKIND
        endif
        
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
                 if ( (abs(cellXc-EtaCenter0) < EtaRadius*dx) .and. &
                      (abs(cellYc-EtaCenter1) < EtaRadius*dy) .and. &
                      (abs(cellZc-EtaCenter2) < EtaRadius*dz) ) then
                    Ersrc(i,j,k) = etaconst
                 endif
                 
              enddo
           enddo
        enddo

     endif ! EtaRadius == 0

  !   emissivity flux along x=0 wall (NGammaDot photons/s/cm^2)
  else if (ProbType == 412) then

     ! place ionization source along left wall (if on this subdomain)
     if (x1L == 0._RKIND) then

        ! compute eta factor for given ionization source, and put on wall
        etaconst = h_nu0*NGammaDot/dy
        ! (the additional factor arises due to the grey source having a
        !  T=1e5 blackbody spectrum)
        etaconst = etaconst * 1.52877652583602_RKIND
        do k=1,Nz,1
           do j=1,Ny,1
              Ersrc(1,j,k) = etaconst
           enddo
        enddo
     endif
     
  !   point-source emissivity at center of every processor
  elseif (ProbType == 414) then

     ! compute eta factor for given ionization source
     etaconst = h_nu0*NGammaDot/dV
        
     ! place ionization source in center of subdomain
     Ersrc(int(Nx/2,IKIND),int(Ny/2,IKIND),int(Nz/2,IKIND)) = etaconst
     
  endif ! ProbType

!!$  write(*,*) 'RadiationSource: individual source is ',etaconst
!!$  etaconst = sum(Ersrc)
!!$  write(*,*) 'RadiationSource: integrated source is ',etaconst
  
  return
end subroutine gFLDProblem_RadiationSource
!=======================================================================

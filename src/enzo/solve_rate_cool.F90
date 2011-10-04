#include "fortran.def"
#include "phys_const.def"
#include "error.def"

!=======================================================================
!/////////////////////  SUBROUTINE SOLVE_RATE  \\\\\\\\\\\\\\\\\\\\\\\\\

  subroutine solve_rate_cool(  &
       d, e, ge, u, v, w, de, &
       HI, HII, HeI, HeII, HeIII, &
       in, jn, kn, nratec, iexpand, imethod, &
       idual, ispecies, imetal, imcool, idust, idim, &
       is, js, ks, ie, je, ke, imax, ih2co, ipiht, igammah, &
       dt, aye, temstart, temend, &
       utem, uxyz, uaye, urho, utim, &
       eta1, eta2, gamma, fh, dtoh, z_solar, &
       k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a, &
       k11a, k12a, k13a, k13dda, k14a, k15a, &
       k16a, k17a, k18a, k19a, k22a, &
       k24, k25, k26, k27, k28, k29, k30, k31, &
       k50a, k51a, k52a, k53a, k54a, k55a, k56a, &
       ndratec, dtemstart, dtemend, h2dusta, &
       ncrna, ncrd1a, ncrd2a, &
       ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa,  &
       ciHeISa, ciHeIIa, reHIIa, reHeII1a, &
       reHeII2a, reHeIIIa, brema, compa, gammaha, &
       comp_xraya, comp_temp, piHI, piHeI, piHeII, &
       HM, H2I, H2II, DI, DII, HDI, metal, &
       hyd01ka, h2k01a, vibha, rotha, rotla, & 
       gpldla, gphdla, hdltea, hdlowa, &
       gaHIa, gaH2a, gaHea, gaHpa, gaela, &
       gasgra, metala, n_xe, xe_start, xe_end, &
       inutot, iradtype, nfreq, imetalregen, &
       iradshield, avgsighp, avgsighep, avgsighe2p, &
       iradtrans, iradcoupled, iradstep, ierr, irt_honly, &
       kphHI, kphHeI, kphHeII, kdissH2I, photogamma, &
       ih2optical, iciecool, ithreebody, ciecoa,  &
       icmbTfloor, iClHeat, &
       clEleFra, clGridRank, clGridDim, &
       clPar1, clPar2, clPar3, clPar4, clPar5, &
       clDataSize, clCooling, clHeating)
!
!  SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1:  January, 1996 by Greg Bryan; converted to KRONOS
!  modified2:  October, 1996 by GB; adapted to AMR
!  modified3:  May,     1999 by GB and Tom Abel, 3bodyH2, solver, HD
!  modified4:  June,    2005 by GB to solve rate & cool at same time
!  modified5:  April,   2009 by JHW to include radiative transfer
!  modified6:  September, 2009 by BDS to include cloudy cooling
!  modified7:  June,    2011 by JHW conversion to F90
!
!  PURPOSE:
!    Solve the multi-species rate and cool equations.
!
!  INPUTS:
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!    de       - electron density field
!    HI,HII   - H density fields (neutral & ionized)
!    HeI/II/III - He density fields
!    DI/II    - D density fields (neutral & ionized)
!    HDI      - neutral HD molecule density field
!    HM       - H- density field
!    H2I      - H_2 (molecular H) density field
!    H2II     - H_2+ density field
!    kph*     - photoionization fields
!    gamma*   - photoheating fields
!
!    is,ie    - start and end indices of active region (zero based)
!    idual    - dual energy formalism flag (0 = off, 1 = on)
!    iexpand  - comoving coordinates flag (0 = off, 1 = on)
!    idim     - dimensionality (rank) of problem
!    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
!    iradshield - flag for crude radiative shielding correction 
!    iradtype - type of radiative field (only used if = 8)
!    imetal   - flag if metal field is active (0 = no, 1 = yes)
!    imcool   - flag if there is metal cooling
!    idust    - flag for H2 formation on dust grains
!    imethod  - Hydro method (0 = PPMDE, 2 = ZEUS-type)
!    ih2co    - flag to include H2 cooling (1 = on, 0 = off)
!    ipiht    - flag to include photoionization heating (1 = on, 0 = off)
!    iradtrans - flag to include radiative transfer (1 = on, 0 = off)
!    iradcoupled - flag to indicate coupled radiative transfer
!    iradstep - flag to indicate intermediate coupled radiative transfer timestep
!
!    fh       - Hydrogen mass fraction (typically 0.76)
!    dtoh     - Deuterium to H mass ratio
!    z_solar  - Solar metal mass fraction
!    dt       - timestep to integrate over
!    aye      - expansion factor (in code units)
!
!    utim     - time units (i.e. code units to CGS conversion factor)
!    uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
!    urho     - density units
!    uxyz     - length units
!    utem     - temperature(-like) units
!
!    temstart, temend - start and end of temperature range for rate table
!    nratec   - dimensions of chemical rate arrays (functions of temperature)
!    dtemstart, dtemend - start and end of dust temperature range
!    ndratec  - extra dimension for H2 formation on dust rate (dust temperature)
!
!    icmbTfloor - flag to include temperature floor from cmb
!    iClHeat    - flag to include cloudy heating
!    clEleFra   - parameter to account for additional electrons from metals 
!    clGridRank - rank of cloudy cooling data grid
!    clGridDim  - array containing dimensions of cloudy data
!    clPar1, clPar2, clPar3, clPar4, clPar5 - arrays containing cloudy grid parameter values
!    clDataSize - total size of flattened 1D cooling data array
!    clCooling  - cloudy cooling data
!    clHeating  - cloudy heating data
!
!  OUTPUTS:
!    update chemical rate densities (HI, HII, etc)
!
!  PARAMETERS:
!    itmax   - maximum allowed sub-cycle iterations
!    mh      - H mass in cgs units
!
!-----------------------------------------------------------------------

    implicit NONE

!  General Arguments

    integer, intent(in) :: in, jn, kn, is, js, ks, ie, je, ke, nratec, &
         imethod, idual, iexpand, ih2co, ipiht, ispecies, imetal, idim,&
         iradtype, nfreq, imetalregen, iradshield, iradtrans, &
         iradcoupled, iradstep, n_xe, imcool, idust, irt_honly, &
         igammah, ih2optical, iciecool, ithreebody, imax, &
         ndratec
    real, intent(in) :: dt, aye, temstart, temend, eta1, eta2, gamma, &
         utim, uxyz, uaye, urho, utem, fh, dtoh, xe_start, xe_end, &
         dtemstart, dtemend, z_solar
    integer, intent(out) :: ierr
    
!  Density, energy and velocity fields fields
    
    real, intent(inout), dimension(in,jn,kn) :: &
         de, HI, HII, HeI, HeII, HeIII, HM, H2I, H2II, DI, DII, HDI, &
         d, ge, e, u, v, w, metal

!  Radiation fields

    real, intent(in), dimension(in,jn,kn) :: &
         kphHI, kphHeI, kphHeII, kdissH2I, photogamma

!  Cooling tables (coolings rates as a function of temperature)

    real, intent(in), dimension(nratec) :: &
         hyd01ka, h2k01a, vibha, rotha, rotla, gpldla, gphdla, hdltea, &
         hdlowa, gaHIa, gaH2a, gaHea, gaHpa, gaela, ciecoa, ceHIa, &
         ceHeIa, ceHeIIa, ciHIa, ciHeIa, ciHeISa, ciHeIIa, reHIIa, &
         reHeII1a, reHeII2a, reHeIIIa, brema, gasgra
    real, intent(in), dimension(nratec, n_xe) :: metala
    real, intent(in) :: inutot(nfreq)
    real, intent(in) :: compa, piHI, piHeI, piHeII, comp_xraya, comp_temp, &
         avgsighp, avgsighep, avgsighe2p, gammaha

!  Chemistry tables (rates as a function of temperature)

    real, intent(in), dimension(nratec) :: &
         k1a , k2a , k3a , k4a , k5a , k6a , k7a , k8a , k9a , &
         k10a, k11a, k12a, k13a, k14a, k15a, k16a, k17a, k18a, &
         k19a, k22a, k50a, k51a, k52a, k53a, k54a, k55a, k56a, &
         ncrna, ncrd1a, ncrd2a
    real, intent(in) :: k13dda(nratec, 7)
    real, intent(in) :: k24, k25, k26, k27, k28, k29, k30, k31
    real, intent(in) :: h2dusta(nratec, ndratec)

!  Cloudy cooling data

    integer, intent(in) :: icmbTfloor, iClHeat, clGridRank, clDataSize
    integer, intent(in) :: clGridDim(clGridRank)
    real, intent(in) :: clEleFra
    real, intent(in) :: clPar1(clGridDim(1)), clPar2(clGridDim(2)), &
         clPar3(clGridDim(3)), clPar4(clGridDim(4)), &
         clPar5(clGridDim(5))
    real, intent(in), dimension(clDataSize) :: clCooling, clHeating
    
!  Parameters

    integer :: itmax
    parameter (itmax = 10000)

#ifdef CONFIG_BFLOAT_4
    real :: tolerance
    parameter (tolerance = 1.0e-05)
#endif

#ifdef CONFIG_BFLOAT_8
    real :: tolerance
    parameter (tolerance = 1.0e-10)
#endif

    double precision :: mh
    parameter (mh = mass_h)

!  Locals

    integer :: i, j, k, iter
    integer :: clGridDim1, clGridDim2, clGridDim3, clGridDim4, clGridDim5
    real :: ttmin, dom, energy, comp1, comp2, olddtit
    double precision :: coolunit, dbase1, tbase1, xbase1, chunit, uvel
    double precision :: heq1, heq2, eqk221, eqk222, eqk131, eqk132, &
         eqt1, eqt2, eqtdef, dheq, heq, dlogtem

    !  row temporaries

    integer, dimension(:), allocatable :: indixe
    real, dimension(:), allocatable :: &
         t1, t2, logtem, tdef, dtit, ttot, p2d, tgas, tgasold, &
         tdust, metallicity, rhoH

    !  Rate equation row temporaries

    real, dimension(:), allocatable :: &
         HIp, HIIp, HeIp, HeIIp, HeIIIp, HMp, H2Ip, H2IIp, dep, &
         dedot, HIdot, dedot_prev, DIp, DIIp, HDIp, HIdot_prev, &
         k24shield, k25shield, k26shield, &
         h2dust, ncrn, ncrd1, ncrd2
    real, dimension(:), allocatable :: &
         k1 , k2 , k3 , k4 , k5 , k6 , k7 , k8 , k9 , k10, k11, &
         k12, k13, k14, k15, k16, k17, k18, k19, k22, k50, k51, &
         k52, k53, k54, k55, k56
    real, allocatable :: k13dd(:, :)

    !  Cooling/heating row locals

    double precision, dimension(:), allocatable :: &
         ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII, reHII, &
         reHeII1, reHeII2, reHeIII, brem, edot, cieco
    real, dimension(:), allocatable :: &
         hyd01k, h2k01, vibh, roth, rotl, gpldl, gphdl, hdlte, &
         hdlow

    !  Iteration mask

    logical, allocatable :: itmask(:)
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================
! Allocations
!=======================================================================
    allocate(indixe(imax))
    allocate(t1(imax))
    allocate(t2(imax))
    allocate(logtem(imax))
    allocate(tdef(imax))
    allocate(dtit(imax))
    allocate(ttot(imax))
    allocate(p2d(imax))
    allocate(tgas(imax))
    allocate(tgasold(imax))
    allocate(HIp(imax))
    allocate(HIIp(imax))
    allocate(HeIp(imax))
    allocate(HeIIp(imax))
    allocate(HeIIIp(imax))
    allocate(HMp(imax))
    allocate(H2Ip(imax))
    allocate(H2IIp(imax))
    allocate(dep(imax))
    allocate(dedot(imax))
    allocate(HIdot(imax))
    allocate(dedot_prev(imax))
    allocate(DIp(imax))
    allocate(DIIp(imax))
    allocate(HDIp(imax))
    allocate(HIdot_prev(imax))
    allocate(k24shield(imax))
    allocate(k25shield(imax))
    allocate(k26shield(imax))
    allocate(k1(imax))
    allocate(k2(imax))
    allocate(k3(imax))
    allocate(k4(imax))
    allocate(k5(imax))
    allocate(k6(imax))
    allocate(k7(imax))
    allocate(k8(imax))
    allocate(k9(imax))
    allocate(k10(imax))
    allocate(k11(imax))
    allocate(k12(imax))
    allocate(k13(imax))
    allocate(k14(imax))
    allocate(k15(imax))
    allocate(k16(imax))
    allocate(k17(imax))
    allocate(k18(imax))
    allocate(k19(imax))
    allocate(k22(imax))
    allocate(k50(imax))
    allocate(k51(imax))
    allocate(k52(imax))
    allocate(k53(imax))
    allocate(k54(imax))
    allocate(k55(imax))
    allocate(k56(imax))
    allocate(k13dd(imax,7))
    allocate(ceHI(imax))
    allocate(ceHeI(imax))
    allocate(ceHeII(imax))
    allocate(ciHI(imax))
    allocate(ciHeI(imax))
    allocate(ciHeIS(imax))
    allocate(ciHeII(imax))
    allocate(reHII(imax))
    allocate(reHeII1(imax))
    allocate(reHeII2(imax))
    allocate(reHeIII(imax))
    allocate(brem(imax))
    allocate(edot(imax))
    allocate(hyd01k(imax))
    allocate(h2k01(imax))
    allocate(vibh(imax))
    allocate(roth(imax))
    allocate(rotl(imax))
    allocate(gpldl(imax))
    allocate(cieco(imax))
    allocate(gphdl(imax))
    allocate(hdlte(imax))
    allocate(hdlow(imax))
    allocate(itmask(imax))
    allocate(tdust(imax))
    allocate(metallicity(imax))
    allocate(rhoH(imax))
    allocate(h2dust(imax))
    allocate(ncrn(imax))
    allocate(ncrd1(imax))
    allocate(ncrd2(imax))
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set error indicator

      ierr = 0
      
!     Set units

      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      uvel     = uxyz / utim
!      chunit = (7.17775d-12)/(2.d0*uvel*uvel*mh)   ! 4.5 eV per H2 formed
      chunit = (1.60218d-12)/(2.d0*uvel*uvel*mh)   ! 1 eV per H2 formed

      dlogtem = (log(temend) - log(temstart))/real(nratec-1)

!  Convert densities from comoving to proper

      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII, &
           HM, H2I, H2II, DI, DII, HDI, metal, &
           is, ie, js, je, ks, ke, &
           in, jn, kn, ispecies, imetal, aye**(-3))

      call ceiling_species(d, de, HI, HII, HeI, HeII, HeIII,&
           HM, H2I, H2II, DI, DII, HDI, metal,&
           is, ie, js, je, ks, ke,&
           in, jn, kn, ispecies, imetal)

!  Loop over zones, and do an entire i-column in one go

      do k = ks+1, ke+1
      do j = js+1, je+1

!       tolerance = 1.0e-06 * dt

!       Set iteration mask to include only cells with radiation in the
!       intermediate coupled chemistry / energy step

         if (iradcoupled .eq. 1 .and. iradstep .eq. 1) then
            do i = is+1, ie+1
               if (kphHI(i,j,k) .gt. 0) then 
                  itmask(i) = .true. 
               else 
                  itmask(i) = .false.
               endif
            enddo
         endif

!     Normal rate solver, but don't double-count the cells with radiation
 
         if (iradcoupled .eq. 1 .and. iradstep .eq. 0) then
            do i = is+1, ie+1
               if (kphHI(i,j,k) .gt. 0) then 
                  itmask(i) = .false. 
               else 
                  itmask(i) = .true.
               endif
            enddo
         endif

!     No radiation timestep coupling

         if (iradcoupled .eq. 0 .or. iradtrans .eq. 0) then 
            do i = is+1, ie+1
               itmask(i) = .true.
            enddo
         endif

!        Set time elapsed to zero for each cell in 1D section

         do i = is+1, ie+1
            ttot(i) = 0.d0
         enddo

!        ------------------ Loop over subcycles ----------------

         do iter = 1, itmax

!           Compute the cooling rate, tgas, tdust, and metallicity for this row

            call cool1d_multi( &
                 d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII, &
                 in, jn, kn, nratec, idual, imethod,               &
                 iexpand, ispecies, imetal, imcool, idust, idim,   &
                 is, ie, j, k, ih2co, ipiht, iter, igammah,        &
                 aye, temstart, temend, z_solar,                   &
                 utem, uxyz, uaye, urho, utim,                     &
                 eta1, eta2, gamma,                                &
                 ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa,            &
                 ciHeISa, ciHeIIa, reHIIa, reHeII1a,               &
                 reHeII2a, reHeIIIa, brema, compa, gammaha,        &
                 comp_xraya, comp_temp,                            &
                 piHI, piHeI, piHeII, comp1, comp2,                &
                 HM, H2I, H2II, DI, DII, HDI, metal,               &
                 hyd01ka, h2k01a, vibha, rotha, rotla,             &
                 hyd01k, h2k01, vibh, roth, rotl,                  &
                 gpldla, gphdla, gpldl, gphdl,                     &
                 hdltea, hdlowa, hdlte, hdlow,                     &
                 gaHIa, gaH2a, gaHea, gaHpa, gaela,                &
                 metala, n_xe, xe_start, xe_end,                   &
                 ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII, &
                 reHII, reHeII1, reHeII2, reHeIII, brem,           &
                 indixe, t1, t2, logtem, tdef, edot,               &
                 tgas, tgasold, p2d, tdust, metallicity, rhoH,     &
                 inutot, iradtype, nfreq, imetalregen,             &
                 iradshield, avgsighp, avgsighep, avgsighe2p,      &
                 iradtrans, photogamma,                            &
                 ih2optical, iciecool, ciecoa, cieco,              &
                 icmbTfloor, iClHeat,                              &
                 clEleFra, clGridRank, clGridDim,                  &
                 clPar1, clPar2, clPar3, clPar4, clPar5,           &
                 clDataSize, clCooling, clHeating,                 &
                 itmask)

!        Look-up rates as a function of temperature for 1D set of zones
!         (maybe should add itmask to this call)

            call lookup_cool_rates1d(temstart, temend, nratec, j, k, &
                 is, ie, imax, iradtype, iradshield, ithreebody,     &
                 in, jn, kn, ispecies, idust,                        &
                 tgas, HI, HII, HeI, HeII, tdust, metallicity,       &
                 k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,  &
                 k11a, k12a, k13a, k13dda, k14a, k15a, k16a,         &
                 k17a, k18a, k19a, k22a,                             &
                 k50a, k51a, k52a, k53a, k54a, k55a, k56a,           &
                 ndratec, dtemstart, dtemend, h2dusta,               &
                 ncrna, ncrd1a, ncrd2a,                              &
                 avgsighp, avgsighep, avgsighe2p, piHI, piHeI,       &
                 k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,            &
                 k11, k12, k13, k14, k15, k16, k17, k18,             &
                 k19, k22, k24, k25, k26,                            &
                 k50, k51, k52, k53, k54, k55,                       &
                 k56, k13dd, k24shield, k25shield, k26shield,        &
                 h2dust, ncrn, ncrd1, ncrd2,                         &
                 t1, t2, tdef, logtem, indixe,                       &
                 dom, coolunit, tbase1, itmask)

!           Compute dedot and HIdot, the rates of change of de and HI
!             (should add itmask to this call)

            call rate_timestep(dedot, HIdot, ispecies, idust,         &
                 de, HI, HII, HeI, HeII, HeIII, d,             &
                 HM, H2I, H2II,                                &
                 in, jn, kn, is, ie, j, k,                     &
                 k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, &
                 k12, k13, k14, k15, k16, k17, k18, k19, k22,  &
                 k24, k25, k26, k27, k28, k29, k30, k31,       &
                 k50, k51, k52, k53, k54, k55, k56,                 &
                 h2dust, ncrn, ncrd1, ncrd2, rhoH, &
                 k56, k24shield, k25shield, k26shield,         &
                 iradtrans, irt_honly, kphHI, kphHeI, kphHeII, & 
                 kdissH2I, itmask, edot, chunit, dom)

!           Find timestep that keeps relative chemical changes below 10%

            do i = is+1, ie+1
               if (itmask(i)) then
!              Bound from below to prevent numerical errors
               
                  if (abs(dedot(i)) .lt. tiny) &
                       dedot(i) = min(tiny,de(i,j,k))
                  if (abs(HIdot(i)) .lt. tiny) &
                       HIdot(i) = min(tiny,HI(i,j,k))

!              If the net rate is almost perfectly balanced then set
!                  it to zero (since it is zero to available precision)

               if (min(abs(k1(i)* de(i,j,k)*HI(i,j,k)), &
                    abs(k2(i)*HII(i,j,k)*de(i,j,k)))/&
                    max(abs(dedot(i)),abs(HIdot(i))) .gt. 1.0e6) then
                  dedot(i) = tiny
                  HIdot(i) = tiny
               endif

!              If the iteration count is high then take the smaller of
!                the calculated dedot and last time step's actual dedot.
!                This is intended to get around the problem of a low
!                electron or HI fraction which is in equilibrium with high
!                individual terms (which all nearly cancel).

               if (iter .gt. 50) then
                  dedot(i) = min(abs(dedot(i)), abs(dedot_prev(i)))
                  HIdot(i) = min(abs(HIdot(i)), abs(HIdot_prev(i)))
               endif

!              compute minimum rate timestep

               olddtit = dtit(i)
               dtit(i) = min(abs(0.1d0*de(i,j,k)/dedot(i)), &
                    abs(0.1d0*HI(i,j,k)/HIdot(i)), &
                    dt-ttot(i), 0.5d0*dt)

              if (d(i,j,k)*dom.gt.1e8.and.edot(i).gt.0.0)then
                ! Equilibrium value for H is:
                !  H = (-1.0 / (4*k22)) * (k13 - sqrt(8 k13 k22 rho + k13^2))
                ! We now want this to change by 10% or less, but we're only
                ! differentiating by dT.  We have de/dt.  We need dT/de.
                ! T = (g-1)*p2d*utem/N; tgas == (g-1)(p2d*utem/N)
                ! dH_eq / dt = (dH_eq/dT) * (dT/de) * (de/dt)
                !   dH_eq / dT (see above; we can calculate the derivative here)
                !   dT / de = utem * (gamma - 1.0) / N == tgas / p2d
                !   de / dt = edot
                ! Now we use our estimate of dT/de to get the estimated
                ! difference in the equilibrium
                eqt2 = min(log(tgas(i)) + 0.1d0*dlogtem, t2(i))
                eqtdef = (eqt2 - t1(i))/(t2(i) - t1(i))
                eqk222 = k22a(indixe(i)) + &
                     (k22a(indixe(i)+1) -k22a(indixe(i)))*eqtdef
                eqk132 = k13a(indixe(i)) + &
                     (k13a(indixe(i)+1) -k13a(indixe(i)))*eqtdef
                heq2 = (-1.d0 / (4.d0*eqk222)) * (eqk132- &
                     sqrt(8.d0*eqk132*eqk222*fh*d(i,j,k)+eqk132**2.d0))

                eqt1 = max(log(tgas(i)) - 0.1d0*dlogtem, t1(i))
                eqtdef = (eqt1 - t1(i))/(t2(i) - t1(i))
                eqk221 = k22a(indixe(i)) + &
                 (k22a(indixe(i)+1) -k22a(indixe(i)))*eqtdef
                eqk131 = k13a(indixe(i)) + &
                 (k13a(indixe(i)+1) -k13a(indixe(i)))*eqtdef
                heq1 = (-1.d0 / (4.d0*eqk221)) * (eqk131- &
                    sqrt(8.d0*eqk131*eqk221*fh*d(i,j,k)+eqk131**2.d0))

                dheq = (abs(heq2-heq1)/(exp(eqt2) - exp(eqt1))) &
                    * (tgas(i)/p2d(i)) * edot(i)
                heq = (-1.d0 / (4.d0*k22(i))) * (k13(i)- &
                    sqrt(8.d0*k13(i)*k22(i)*fh*d(i,j,k)+k13(i)**2.d0))
                !write(0,*) heq2, heq1, eqt2, eqt1, tgas(i), p2d(i), 
!     &                     edot(i)
                if(d(i,j,k)*dom.gt.1e18.and.i.eq.4)write(0, *) &
                   HI(i,j,k)/heq, edot(i),tgas(i)
                dtit(i) = min(dtit(i), 0.1d0*heq/dheq)
              endif
              if (iter.gt.10) dtit(i) = min(olddtit*1.5d0, dtit(i))

#define DONT_WRITE_COOLING_DEBUG
#ifdef WRITE_COOLING_DEBUG
!              Output some debugging information if required
#ifndef _OPENMP
              if (dtit(i)/dt .lt. 1.0e-2 .and. iter .gt. 800 .and. &
                   abs((dt-ttot(i))/dt) .gt. 1.0d-3) then
                 write(4,1000) iter,i,j,k,dtit(i), &
                      ttot(i),dt,de(i,j,k),dedot(i),HI(i,j,k),HIdot(i), &
                   tgas(i), dedot_prev(i), HIdot_prev(i)
                 write(4,1100) HI(i,j,k),HII(i,j,k), &
                      HeI(i,j,k),HeII(i,j,k),HeIII(i,j,k), &
                      HM(i,j,k),H2I(i,j,k),H2II(i,j,k),de(i,j,k)
                 write(4,1100) &
                      -      k1(i) *de(i,j,k)    *HI(i,j,k),       &
                      -      k7(i) *de(i,j,k)    *HI(i,j,k),       &
                      -      k8(i) *HM(i,j,k)    *HI(i,j,k),       &
                      -      k9(i) *HII(i,j,k)   *HI(i,j,k),       &
                      -      k10(i)*H2II(i,j,k)  *HI(i,j,k)/2.d0,  &
                      - 2.d0*k22(i)*HI(i,j,k)**2 *HI(i,j,k),       &
                      +      k2(i) *HII(i,j,k)   *de(i,j,k) ,      &
                      + 2.d0*k13(i)*HI(i,j,k)    *H2I(i,j,k)/2.d0, &
                      +      k11(i)*HII(i,j,k)   *H2I(i,j,k)/2.d0, &
                      + 2.d0*k12(i)*de(i,j,k)    *H2I(i,j,k)/2.d0, &
                      +      k14(i)*HM(i,j,k)    *de(i,j,k),       &
                      +      k15(i)*HM(i,j,k)    *HI(i,j,k),       &
                      + 2.d0*k16(i)*HM(i,j,k)    *HII(i,j,k),      &
                      + 2.d0*k18(i)*H2II(i,j,k)  *de(i,j,k)/2.d0,  &
                      +      k19(i)*H2II(i,j,k)  *HM(i,j,k)/2.d0
              endif
#endif /* _OPENMP */
 1000          format(i5,3(i3,1x),1p,11(e11.3))
 1100          format(1p,20(e11.3))
#endif /* WRITE_COOLING_DEBUG */
            else               ! itmask
               dtit(i) = dt;
            endif
            enddo               ! end loop over i

!           Compute maximum timestep for cooling/heating

            do i = is+1, ie+1
               if (itmask(i)) then
!              Set energy per unit volume of this cell based in the pressure
!              (the gamma used here is the right one even for H2 since p2d 
!               is calculated with this gamma).

               energy = max(p2d(i)/(gamma-1.d0), tiny)

!              This is an alternate energy calculation, based directly on
!              the code's specific energy field, which differs from the above
!              only if using the dual energy formalism.

!              energy = max(ge(i,j,k)*d(i,j,k), p2d(i)/(gamma-1.0), 
!     &                    tiny)
!              if (energy .lt. tiny) energy = d(i,j,k)*(e(i,j,k) - 
!     &              0.5*(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))

!              If the temperature is at the bottom of the temperature look-up 
!              table and edot < 0, then shut off the cooling.

               if (tgas(i) .le. 1.01d0*temstart .and. edot(i) .lt. 0.d0)&
                    edot(i) = tiny
	       if (abs(edot(i)) .lt. tiny) edot(i) = tiny
!
!              Compute timestep for 10% change

!              if (iter .gt. 100) then
!                 dtit(i) = min(real(abs(0.1*energy/edot(i))), 
!     &                        dt-ttot(i), dtit(i))
!              else
               dtit(i) = min(real(abs(0.1d0*energy/edot(i))), &
                    dt-ttot(i), dtit(i))
!              endif

               if (dtit(i) .ne. dtit(i)) & !#####
                    write(6,*) 'HUGE dtit :: ', energy, edot(i), dtit(i), &
                    dt, ttot(i), abs(0.1d0*energy/edot(i)), &
                    real(abs(0.1d0*energy/edot(i)))

#define NO_FORTRAN_DEBUG
#ifdef FORTRAN_DEBUG
               if (ge(i,j,k) .le. 0.0 .and. idual .eq. 1) &
                    write(6,*) 'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter
               if (idual .eq. 1 .and. &
                    ge(i,j,k)+edot(i)/d(i,j,k)*dtit(i) .le. 0.d0) &
                    write(6,*) i,j,k,iter,ge(i,j,k),edot(i),tgas(i), &
                    energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k)
#endif /* FORTRAN_DEBUG */

#ifdef WRITE_COOLING_DEBUG
!              If the timestep is too small, then output some debugging info

#ifndef _OPENMP
               if (((dtit(i)/dt .lt. 1.0e-2 .and. iter .gt. 800) &
                    .or. iter .gt. itmax-100) .and. &
                    abs((dt-ttot(i))/dt) .gt. 1.0d-3) &
                    write(3,2000) i,j,k,iter,ge(i,j,k),edot(i),tgas(i), &
                    energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k),dtit(i)
#endif /* _OPENMP */
 2000          format(4(i4,1x),1p,10(e14.3))
#endif /* WRITE_COOLING_DEBUG */
            endif   ! itmask
            enddo   ! end loop over i

!           Update total and gas energy

            do i = is+1, ie+1
               if (itmask(i)) then
               e(i,j,k)  = e(i,j,k) + edot(i)/d(i,j,k)*dtit(i)
#ifdef WRITE_COOLING_DEBUG
               if (e(i,j,k) .ne. e(i,j,k)) &
                    write(3,*) edot(i),d(i,j,k),dtit(i)
#endif /* WRITE_COOLING_DEBUG */

!              If using the dual energy formalism, there are 2 energy fields

               if (idual .eq. 1) then
                  ge(i,j,k) = ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i)

!                 Alternate energy update schemes (not currently used)

!                 ge(i,j,k) = max(ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i),
!     &                      0.5d0*ge(i,j,k))
!                 if (ge(i,j,k) .le. tiny) ge(i,j,k) = (energy + 
!     &              edot(i)*dtit(i))/d(i,j,k)
#ifdef WRITE_COOLING_DEBUG
                  if (ge(i,j,k) .le. 0.d0) write(3,*) &
                       'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter
#endif WRITE_COOLING_DEBUG
               endif
            endif               ! itmask
            enddo

!           Solve rate equations with one linearly implicit Gauss-Seidel 
!           sweep of a backward Euler method ---

            call step_rate(de, HI, HII, HeI, HeII, HeIII, d,   &
                 HM, H2I, H2II, DI, DII, HDI, dtit,            &
                 in, jn, kn, is, ie, j, k, ispecies, idust,          &
                 k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, &
                 k12, k13, k14, k15, k16, k17, k18, k19, k22,  &
                 k24, k25, k26, k27, k28, k29, k30, k31,       &
                 k50, k51, k52, k53, k54, k55, k56,                 &
                 h2dust, rhoH, &
                 k24shield, k25shield, k26shield,         &
                 HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,          &
                 HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,            &
                 dedot_prev, HIdot_prev,                       &
                 iradtrans, irt_honly, kphHI, kphHeI, kphHeII, & 
                 kdissH2I, itmask)

!           Add the timestep to the elapsed time for each cell and find
!            minimum elapsed time step in this row

            ttmin = huge
            do i = is+1, ie+1
               ttot(i) = min(ttot(i) + dtit(i), dt)
               if (abs(dt-ttot(i)) .lt. 0.001d0*dt) itmask(i) = .false.
               if (ttot(i).lt.ttmin) ttmin = ttot(i)
            enddo

!           If all cells are done (on this slice), then exit

            if (abs(dt-ttmin) .lt. tolerance*dt) go to 9999

!           Next subcycle iteration

         enddo

 9999    continue

!       Abort if iteration count exceeds maximum

         if (iter .gt. itmax) then
	    write(0,*) 'inside if statement solve rate cool:',is,ie
            write(6,*) 'MULTI_COOL iter > ',itmax,' at j,k =',j,k
            write(0,*) 'FATAL error (2) in MULTI_COOL'
            write(0,'(" dt = ",1pe10.3," ttmin = ",1pe10.3)') dt, ttmin
            write(0,'((16(1pe8.1)))') (dtit(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (ttot(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (edot(i),i=is+1,ie+1)
            write(0,'((16(l3)))') (itmask(i),i=is+1,ie+1)
            WARNING_MESSAGE
         endif

         if (iter .gt. itmax/2) then
            write(6,*) 'MULTI_COOL iter,j,k =',iter,j,k
         end if
!     
!     Next j,k
!     
       enddo
      enddo

!     Convert densities back to comoving from proper

      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII, &
           HM, H2I, H2II, DI, DII, HDI, metal, &
           is, ie, js, je, ks, ke, &
           in, jn, kn, ispecies, imetal, aye**3)

!     Correct the species to ensure consistency (i.e. type conservation)

      call make_consistent(de, HI, HII, HeI, HeII, HeIII, &
           HM, H2I, H2II, DI, DII, HDI, metal, &
           d, is, ie, js, je, ks, ke, &
           in, jn, kn, imax, ispecies, imetal, fh, dtoh)

!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================
! Deallocations
!=======================================================================

      deallocate(indixe)
      deallocate(t1)
      deallocate(t2)
      deallocate(logtem)
      deallocate(tdef)
      deallocate(dtit)
      deallocate(ttot)
      deallocate(p2d)
      deallocate(tgas)
      deallocate(tgasold)
      deallocate(HIp)
      deallocate(HIIp)
      deallocate(HeIp)
      deallocate(HeIIp)
      deallocate(HeIIIp)
      deallocate(HMp)
      deallocate(H2Ip)
      deallocate(H2IIp)
      deallocate(dep)
      deallocate(dedot)
      deallocate(HIdot)
      deallocate(dedot_prev)
      deallocate(DIp)
      deallocate(DIIp)
      deallocate(HDIp)
      deallocate(HIdot_prev)
      deallocate(k24shield)
      deallocate(k25shield)
      deallocate(k26shield)
      deallocate(k1)
      deallocate(k2)
      deallocate(k3)
      deallocate(k4)
      deallocate(k5)
      deallocate(k6)
      deallocate(k7)
      deallocate(k8)
      deallocate(k9)
      deallocate(k10)
      deallocate(k11)
      deallocate(k12)
      deallocate(k13)
      deallocate(k14)
      deallocate(k15)
      deallocate(k16)
      deallocate(k17)
      deallocate(k18)
      deallocate(k19)
      deallocate(k22)
      deallocate(k50)
      deallocate(k51)
      deallocate(k52)
      deallocate(k53)
      deallocate(k54)
      deallocate(k55)
      deallocate(k56)
      deallocate(k13dd)
      deallocate(ceHI)
      deallocate(ceHeI)
      deallocate(ceHeII)
      deallocate(ciHI)
      deallocate(ciHeI)
      deallocate(ciHeIS)
      deallocate(ciHeII)
      deallocate(reHII)
      deallocate(reHeII1)
      deallocate(reHeII2)
      deallocate(reHeIII)
      deallocate(brem)
      deallocate(edot)
      deallocate(hyd01k)
      deallocate(h2k01)
      deallocate(vibh)
      deallocate(roth)
      deallocate(rotl)
      deallocate(cieco)
      deallocate(gpldl)
      deallocate(gphdl)
      deallocate(hdlte)
      deallocate(hdlow)
      deallocate(itmask)
      deallocate(tdust)
      deallocate(metallicity)
      deallocate(rhoH)
      deallocate(h2dust)
      deallocate(ncrn)
      deallocate(ncrd1)
      deallocate(ncrd2)


      return
    end subroutine solve_rate_cool

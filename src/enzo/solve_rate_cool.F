#include "fortran.def"
#include "phys_const.def"
#include "error.def"

!=======================================================================
!/////////////////////  SUBROUTINE SOLVE_RATE  \\\\\\\\\\\\\\\\\\\\\\\\\

      subroutine solve_rate_cool(d, e, ge, u, v, w, de,
     &                HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, ispecies, imetal, imcool, idust, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht, igammah,
     &                dt, aye, temstart, temend, 
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma, fh, dtoh, z_solar,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a,
     &                k16a, k17a, k18a, k19a, k22a,
     &                k24, k25, k26, k27, k28, k29, k30, k31,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                ndratec, dtemstart, dtemend, h2dusta,
     &                ncrna, ncrd1a, ncrd2a,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
     &                gpldla, gphdla, hdltea, hdlowa,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela, 
     &                gasgra, metala, n_xe, xe_start, xe_end,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, iradcoupled, iradstep, ierr,
     &                irt_honly,
     &                kphHI, kphHeI, kphHeII, kdissH2I, 
     &                photogamma, 
     &                ih2optical, iciecool, ithreebody, ciecoa, 
     &                icmbTfloor, iClHeat,
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating)

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
#include "fortran_types.def"

!  General Arguments

      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke, nratec, imethod,
     &        idual, iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        iradtype, nfreq, imetalregen, iradshield, iradtrans,
     &        iradcoupled, iradstep, n_xe, ierr, imcool, idust,
     &        irt_honly, igammah, ih2optical, iciecool, ithreebody,
     &        ndratec
      R_PREC    dt, aye, temstart, temend, eta1, eta2, gamma,
     &        utim, uxyz, uaye, urho, utem, fh, dtoh, z_solar, 
     &        xe_start, xe_end, dtemstart, dtemend

!  Density, energy and velocity fields fields

      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        metal(in,jn,kn)

!  Radiation fields

      R_PREC kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &     kdissH2I(in,jn,kn)
      R_PREC photogamma(in,jn,kn)

!  Cooling tables (coolings rates as a function of temperature)

      R_PREC    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec)
      R_PREC    gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec), gasgra(nratec), 
     &        ciecoa(nratec)
      R_PREC    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      R_PREC    metala(nratec, n_xe)
      R_PREC    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p
      R_PREC    gammaha 

!  Chemistry tables (rates as a function of temperature)

      R_PREC k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &     k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &     k54a(nratec), k55a(nratec), k56a(nratec),
     &     k13dda(nratec, 7), h2dusta(nratec, ndratec),
     &     ncrna(nratec), ncrd1a(nratec), ncrd2a(nratec),
     &     k24, k25, k26, k27, k28, k29, k30, k31

!  Cloudy cooling data

      INTG_PREC icmbTfloor, iClHeat, clGridRank, clDataSize
      INTG_PREC clGridDim(5)
      R_PREC clEleFra
      R_PREC clPar1(clGridDim(1)), clPar2(clGridDim(2)), 
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)), 
     &     clPar5(clGridDim(5))
      R_PREC clCooling(clDataSize), clHeating(clDataSize)

!  Parameters

      INTG_PREC itmax, ijk
      parameter (itmax = 10000, ijk = MAX_ANY_SINGLE_DIRECTION)

#ifdef CONFIG_BFLOAT_4
      R_PREC tolerance
      parameter (tolerance = 1.e-5_RKIND)
#endif

#ifdef CONFIG_BFLOAT_8
      R_PREC tolerance
      parameter (tolerance = 1.e-10_RKIND)
#endif

      real*8 mh
      parameter (mh = mass_h)

!  Locals

      INTG_PREC i, j, k, iter
      INTG_PREC clGridDim1, clGridDim2, clGridDim3, clGridDim4, 
     &     clGridDim5
      R_PREC ttmin, dom, energy, comp1, comp2
      real*8 coolunit, dbase1, tbase1, xbase1, chunit, uvel
      real*8 heq1, heq2, eqk221, eqk222, eqk131, eqk132,
     &                 eqt1, eqt2, eqtdef, dheq, heq, dlogtem

!  row temporaries

      INTG_PREC indixe(ijk)
      R_PREC t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), 
     &     dtit(ijk), ttot(ijk), p2d(ijk), tgas(ijk), tgasold(ijk),
     &     tdust(ijk), metallicity(ijk), rhoH(ijk), olddtit

!  Rate equation row temporaries

      R_PREC HIp(ijk), HIIp(ijk), HeIp(ijk), HeIIp(ijk), HeIIIp(ijk),
     &     HMp(ijk), H2Ip(ijk), H2IIp(ijk),
     &     dep(ijk), dedot(ijk),HIdot(ijk), dedot_prev(ijk),
     &     DIp(ijk), DIIp(ijk), HDIp(ijk), HIdot_prev(ijk),
     &     k24shield(ijk), k25shield(ijk), k26shield(ijk)
      R_PREC k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k22(ijk),
     &     k50(ijk), k51(ijk), k52(ijk), k53(ijk), k54(ijk),
     &     k55(ijk), k56(ijk), k13dd(ijk, 7), h2dust(ijk),
     &     ncrn(ijk), ncrd1(ijk), ncrd2(ijk)

!  Cooling/heating row locals

      real*8 ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &     ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &     reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk), 
     &     brem(ijk), edot(ijk)
      R_PREC hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), rotl(ijk), 
     &     gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk), cieco(ijk)

!  Iteration mask

      LOGIC_PREC itmask(ijk)
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
c      chunit = (7.17775e-12_RKIND)/(2._RKIND*uvel*uvel*mh)   ! 4.5 eV per H2 formed
      chunit = (1.60218e-12_RKIND)/(2._RKIND*uvel*uvel*mh)   ! 1 eV per H2 formed

      dlogtem = (log(temend) - log(temstart))/REAL(nratec-1,RKIND)

!  Convert densities from comoving to proper

      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**(-3))

      call ceiling_species(d, de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal,
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal)

!  Loop over zones, and do an entire i-column in one go

      do k = ks+1, ke+1
      do j = js+1, je+1

!       tolerance = 1.0e-06 * dt

!       Set iteration mask to include only cells with radiation in the
!       intermediate coupled chemistry / energy step

         if (iradcoupled == 1 .and. iradstep == 1) then
            do i = is+1, ie+1
               if (kphHI(i,j,k) > 0) then 
                  itmask(i) = .true. 
               else 
                  itmask(i) = .false.
               endif
            enddo
         endif

!     Normal rate solver, but don't double-count the cells with radiation
 
         if (iradcoupled == 1 .and. iradstep == 0) then
            do i = is+1, ie+1
               if (kphHI(i,j,k) > 0) then 
                  itmask(i) = .false. 
               else 
                  itmask(i) = .true.
               endif
            enddo
         endif

!     No radiation timestep coupling

         if (iradcoupled == 0 .or. iradtrans == 0) then 
            do i = is+1, ie+1
               itmask(i) = .true.
            enddo
         endif

!        Set time elapsed to zero for each cell in 1D section

         do i = is+1, ie+1
            ttot(i) = 0._RKIND
         enddo

!        ------------------ Loop over subcycles ----------------

         do iter = 1, itmax

!           Compute the cooling rate, tgas, tdust, and metallicity for this row

            call cool1d_multi(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, imetal, imcool, idust, idim,
     &                is, ie, j, k, ih2co, ipiht, iter, igammah,
     &                aye, temstart, temend, z_solar,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                gasgra, metala, n_xe, xe_start, xe_end,
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, p2d, tdust, metallicity, rhoH, 
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, photogamma,
     &                ih2optical, iciecool, ciecoa, cieco,
     &                icmbTfloor, iClHeat,
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating,
     &                itmask)

!        Look-up rates as a function of temperature for 1D set of zones
!         (maybe should add itmask to this call)

            call lookup_cool_rates1d(temstart, temend, nratec, j, k,
     &               is, ie, ijk, iradtype, iradshield, ithreebody,
     &               in, jn, kn, ispecies, idust, 
     &               tgas, HI, HII, HeI, HeII, tdust, metallicity, 
     &               k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &               k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &               k17a, k18a, k19a, k22a,
     &               k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &               ndratec, dtemstart, dtemend, h2dusta, 
     &               ncrna, ncrd1a, ncrd2a, 
     &               avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     &               k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &               k11, k12, k13, k14, k15, k16, k17, k18,
     &               k19, k22, k24, k25, k26, 
     &               k50, k51, k52, k53, k54, k55,
     &               k56, k13dd, k24shield, k25shield, k26shield,
     &               h2dust, ncrn, ncrd1, ncrd2, 
     &               t1, t2, tdef, logtem, indixe, 
     &               dom, coolunit, tbase1, itmask)

!           Compute dedot and HIdot, the rates of change of de and HI
!             (should add itmask to this call)

            call rate_timestep(dedot, HIdot, ispecies, idust,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55, k56, 
     &                     h2dust, ncrn, ncrd1, ncrd2, rhoH, 
     &                     k24shield, k25shield, k26shield,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, 
     &                     kdissH2I, itmask, edot, chunit, dom)

!           Find timestep that keeps relative chemical changes below 10%

            do i = is+1, ie+1
               if (itmask(i)) then
!              Bound from below to prevent numerical errors
               
	       if (abs(dedot(i)) < tiny) 
     &             dedot(i) = min(REAL(tiny,RKIND),de(i,j,k))
	       if (abs(HIdot(i)) < tiny)
     &             HIdot(i) = min(REAL(tiny,RKIND),HI(i,j,k))

!              If the net rate is almost perfectly balanced then set
!                  it to zero (since it is zero to available precision)

               if (min(abs(k1(i)* de(i,j,k)*HI(i,j,k)),
     &                 abs(k2(i)*HII(i,j,k)*de(i,j,k)))/
     &              max(abs(dedot(i)),abs(HIdot(i))) 
     &              > 1.e6) then
                  dedot(i) = tiny
                  HIdot(i) = tiny
               endif

!              If the iteration count is high then take the smaller of
!                the calculated dedot and last time step's actual dedot.
!                This is intended to get around the problem of a low
!                electron or HI fraction which is in equilibrium with high
!                individual terms (which all nearly cancel).

               if (iter > 50) then
                  dedot(i) = min(abs(dedot(i)), abs(dedot_prev(i)))
                  HIdot(i) = min(abs(HIdot(i)), abs(HIdot_prev(i)))
               endif

!              compute minimum rate timestep

               olddtit = dtit(i)
               dtit(i) = min(abs(0.1_RKIND*de(i,j,k)/dedot(i)), 
     &                       abs(0.1_RKIND*HI(i,j,k)/HIdot(i)),
     &                       dt-ttot(i), 0.5_RKIND*dt)

              if ( d(i,j,k)*dom > 1.e8 .and. 
     &             edot(i) > 0.0 .and. 
     &             ispecies > 1 ) then
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
                eqt2 = min(log(tgas(i)) + 0.1_RKIND*dlogtem, t2(i))
                eqtdef = (eqt2 - t1(i))/(t2(i) - t1(i))
                eqk222 = k22a(indixe(i)) +
     &            (k22a(indixe(i)+1) -k22a(indixe(i)))*eqtdef
                eqk132 = k13a(indixe(i)) +
     &            (k13a(indixe(i)+1) -k13a(indixe(i)))*eqtdef
                heq2 = (-1._RKIND / (4._RKIND*eqk222)) * (eqk132-
     &               sqrt(8._RKIND*eqk132*eqk222*fh*d(i,j,k)
     &               +eqk132**2._RKIND))

                eqt1 = max(log(tgas(i)) - 0.1_RKIND*dlogtem, t1(i))
                eqtdef = (eqt1 - t1(i))/(t2(i) - t1(i))
                eqk221 = k22a(indixe(i)) +
     &            (k22a(indixe(i)+1) -k22a(indixe(i)))*eqtdef
                eqk131 = k13a(indixe(i)) +
     &            (k13a(indixe(i)+1) -k13a(indixe(i)))*eqtdef
                heq1 = (-1._RKIND / (4._RKIND*eqk221)) * (eqk131-
     &               sqrt(8._RKIND*eqk131*eqk221*fh*d(i,j,k)
     &               +eqk131**2._RKIND))

                dheq = (abs(heq2-heq1)/(exp(eqt2) - exp(eqt1)))
     &               * (tgas(i)/p2d(i)) * edot(i)
                heq = (-1._RKIND / (4._RKIND*k22(i))) * (k13(i)-
     &               sqrt(8._RKIND*k13(i)*k22(i)*fh*d(i,j,k)
     &               +k13(i)**2._RKIND))
                !write(0,*) heq2, heq1, eqt2, eqt1, tgas(i), p2d(i), 
!     &                     edot(i)
                if(d(i,j,k)*dom>1e18 .and. i==4)write(0,*)
     &              HI(i,j,k)/heq, edot(i),tgas(i)
                dtit(i) = min(1.d0*dtit(i), 0.1d0*heq/dheq)
              endif
              if (iter>10) 
     &             dtit(i) = min(olddtit*1.5_RKIND, dtit(i))

#define DONT_WRITE_COOLING_DEBUG
#ifdef WRITE_COOLING_DEBUG
!              Output some debugging information if required

               if (dtit(i)/dt > 1.e-2 .and. iter > 800
     &             .and. abs((dt-ttot(i))/dt) > 1.e-3) then
                  write(4,1000) iter,i,j,k,dtit(i),
     &              ttot(i),dt,de(i,j,k),dedot(i),HI(i,j,k),HIdot(i),
     &              tgas(i), dedot_prev(i), HIdot_prev(i)
                  write(4,1100) HI(i,j,k),HII(i,j,k),
     &              HeI(i,j,k),HeII(i,j,k),HeIII(i,j,k),
     &              HM(i,j,k),H2I(i,j,k),H2II(i,j,k),de(i,j,k)
                  write(4,1100)
     &               -      k1(i) *de(i,j,k)    *HI(i,j,k)  ,
     &               -      k7(i) *de(i,j,k)    *HI(i,j,k),
     &               -      k8(i) *HM(i,j,k)    *HI(i,j,k),
     &               -      k9(i) *HII(i,j,k)   *HI(i,j,k),
     &               -      k10(i)*H2II(i,j,k)  *HI(i,j,k)/2._RKIND,
     &               - 2._RKIND*k22(i)*HI(i,j,k)**2 *HI(i,j,k),
     &               +      k2(i) *HII(i,j,k)   *de(i,j,k) ,
     &               + 2._RKIND*k13(i)*HI(i,j,k)   *H2I(i,j,k)/2._RKIND,
     &               +      k11(i)*HII(i,j,k)   *H2I(i,j,k)/2._RKIND,
     &               + 2._RKIND*k12(i)*de(i,j,k)   *H2I(i,j,k)/2._RKIND,
     &               +      k14(i)*HM(i,j,k)    *de(i,j,k),
     &               +      k15(i)*HM(i,j,k)    *HI(i,j,k),
     &               + 2._RKIND*k16(i)*HM(i,j,k)    *HII(i,j,k),
     &               + 2._RKIND*k18(i)*H2II(i,j,k)  *de(i,j,k)/2._RKIND,
     &               +      k19(i)*H2II(i,j,k)  *HM(i,j,k)/2._RKIND
               endif
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

               energy = max(p2d(i)/(gamma-1._RKIND), tiny)

!              This is an alternate energy calculation, based directly on
!              the code's specific energy field, which differs from the above
!              only if using the dual energy formalism.

!              energy = max(ge(i,j,k)*d(i,j,k), p2d(i)/(gamma-1.0), 
!     &                    tiny)
!              if (energy < tiny) energy = d(i,j,k)*(e(i,j,k) - 
!     &              0.5*(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))

!              If the temperature is at the bottom of the temperature look-up 
!              table and edot < 0, then shut off the cooling.

               if (tgas(i) <= 1.01_RKIND*temstart .and. 
     &              edot(i) < 0.0) edot(i) = tiny
	       if (abs(edot(i)) < tiny) edot(i) = tiny
!
!              Compute timestep for 10% change

!              if (iter > 100) then
!                 dtit(i) = min(REAL(abs(0.1*energy/edot(i)),RKIND), 
!     &                        dt-ttot(i), dtit(i))
!              else
               dtit(i) = min(REAL(abs(0.1_RKIND*energy/edot(i)),RKIND), 
     &              dt-ttot(i), dtit(i))
!              endif

               if (dtit(i) .ne. dtit(i)) !#####
     &           write(6,*) 'HUGE dtit :: ', energy, edot(i), dtit(i),
     &                      dt, ttot(i), abs(0.1_RKIND*energy/edot(i)), 
     &                      REAL(abs(0.1_RKIND*energy/edot(i)),RKIND)

#define NO_FORTRAN_DEBUG
#ifdef FORTRAN_DEBUG
               if (ge(i,j,k) <= 0.0 .and. idual == 1)
     &           write(6,*) 'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter
               if (idual == 1 .and.
     &           ge(i,j,k)+edot(i)/d(i,j,k)*dtit(i) <= 0.0)
     &           write(6,*) i,j,k,iter,ge(i,j,k),edot(i),tgas(i),
     &              energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k)
#endif /* FORTRAN_DEBUG */

#ifdef WRITE_COOLING_DEBUG
!              If the timestep is too small, then output some debugging info

               if (((dtit(i)/dt < 1.e-2 .and. iter > 800) 
     &               .or. iter > itmax-100) .and.
     &              abs((dt-ttot(i))/dt) > 1.e-3)
     &           write(3,2000) i,j,k,iter,ge(i,j,k),edot(i),tgas(i),
     &             energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k),dtit(i)
 2000          format(4(i4,1x),1p,10(e14.3))
#endif /* WRITE_COOLING_DEBUG */
            endif   ! itmask
            enddo   ! end loop over i

!           Update total and gas energy

            do i = is+1, ie+1
               if (itmask(i)) then
               e(i,j,k)  = e(i,j,k) + edot(i)/d(i,j,k)*dtit(i)
#ifdef WRITE_COOLING_DEBUG
               if (e(i,j,k) .ne. e(i,j,k)) 
     &              write(3,*) edot(i),d(i,j,k),dtit(i)
#endif /* WRITE_COOLING_DEBUG */

!              If using the dual energy formalism, there are 2 energy fields

               if (idual == 1) then
                  ge(i,j,k) = ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i)

!                 Alternate energy update schemes (not currently used)

!                 ge(i,j,k) = max(ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i),
!     &                      0.5_RKIND*ge(i,j,k))
!                 if (ge(i,j,k) <= tiny) ge(i,j,k) = (energy + 
!     &              edot(i)*dtit(i))/d(i,j,k)
#ifdef WRITE_COOLING_DEBUG
                  if (ge(i,j,k) <= 0.0) write(3,*)
     &                 'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter
#endif /* WRITE_COOLING_DEBUG */
               endif
            endif               ! itmask
            enddo

!           Solve rate equations with one linearly implicit Gauss-Seidel 
!           sweep of a backward Euler method ---

            call step_rate(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, dtit,
     &                     in, jn, kn, is, ie, j, k, ispecies, idust,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55, k56, 
     &                     h2dust, rhoH,
     &                     k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
     &                     dedot_prev, HIdot_prev,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, 
     &                     kdissH2I, itmask)

!           Add the timestep to the elapsed time for each cell and find
!            minimum elapsed time step in this row

            ttmin = huge
            do i = is+1, ie+1
               ttot(i) = min(ttot(i) + dtit(i), dt)
               if (abs(dt-ttot(i)) < 0.001_RKIND*dt) 
     &              itmask(i) = .false.
               if (ttot(i)<ttmin) ttmin = ttot(i)
            enddo

!           If all cells are done (on this slice), then exit

            if (abs(dt-ttmin) < tolerance*dt) go to 9999

!           Next subcycle iteration

         enddo

 9999    continue

!       Abort if iteration count exceeds maximum

         if (iter > itmax) then
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

         if (iter > itmax/2) then
            write(6,*) 'MULTI_COOL iter,j,k =',iter,j,k
         end if
!     
!     Next j,k
!     
       enddo
      enddo

!     Convert densities back to comoving from proper

      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**3)

!     Correct the species to ensure consistency (i.e. type conservation)

      call make_consistent(de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal, 
     &                     d, is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal, fh, dtoh)

      return
      end

c -----------------------------------------------------------
!   This routine scales the density fields from comoving to
!     proper densities (and back again).

      subroutine scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, factor)
c -------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!     Arguments

      INTG_PREC in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      R_PREC     d(in,jn,kn),metal(in,jn,kn)
      R_PREC    factor

!     locals

      INTG_PREC i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*factor
               de(i,j,k)    = de(i,j,k)*factor
               HI(i,j,k)    = HI(i,j,k)*factor
               HII(i,j,k)   = HII(i,j,k)*factor
               HeI(i,j,k)   = HeI(i,j,k)*factor
               HeII(i,j,k)  = HeII(i,j,k)*factor
               HeIII(i,j,k) = HeIII(i,j,k)*factor
            enddo
         enddo
      enddo
      if (ispecies > 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*factor
                  H2I(i,j,k)  = H2I(i,j,k)*factor
                  H2II(i,j,k) = H2II(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
      if (ispecies > 2) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*factor
                  DII(i,j,k) = DII(i,j,k)*factor
                  HDI(i,j,k) = HDI(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
      if (imetal == 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*factor
               enddo
            enddo
         enddo
      endif

      return
      end

c -----------------------------------------------------------
!   This routine ensures that the species aren't below tiny.

      subroutine ceiling_species(d, de, HI, HII, HeI, HeII, HeIII,
     &                           HM, H2I, H2II, DI, DII, HDI, metal,
     &                           is, ie, js, je, ks, ke,
     &                           in, jn, kn, ispecies, imetal)
c -------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!     Arguments

      INTG_PREC in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC     d(in,jn,kn)
      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      R_PREC metal(in,jn,kn)

!     locals

      INTG_PREC i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               de(i,j,k)    = max(de(i,j,k), tiny)
               HI(i,j,k)    = max(HI(i,j,k), tiny)
               HII(i,j,k)   = max(HII(i,j,k), tiny)
               HeI(i,j,k)   = max(HeI(i,j,k), tiny)
               HeII(i,j,k)  = max(HeII(i,j,k), tiny)
               HeIII(i,j,k) = max(HeIII(i,j,k), 1.e-5_RKIND*tiny)
            enddo
         enddo
      enddo
      if (ispecies > 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  HM(i,j,k)   = max(HM(i,j,k), tiny)
                  H2I(i,j,k)  = max(H2I(i,j,k), tiny)
                  H2II(i,j,k) = max(H2II(i,j,k), tiny)
               enddo
            enddo
         enddo
      endif
      if (ispecies > 2) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DI(i,j,k)  = max(DI(i,j,k), tiny)
                  DII(i,j,k) = max(DII(i,j,k), tiny)
                  HDI(i,j,k) = max(HDI(i,j,k), tiny)
               enddo
            enddo
         enddo
      endif
      if (imetal == 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = min(max(metal(i,j,k),
     &                               tiny), 0.9_RKIND*d(i,j,k))
               enddo
            enddo
         enddo
      endif

      return
      end



! -----------------------------------------------------------
! This routine uses the temperature to look up the chemical
!   rates which are tabulated in a log table as a function
!   of temperature.

      subroutine lookup_cool_rates1d(temstart, temend, nratec, j, k,
     &                is, ie, ijk, iradtype, iradshield, ithreebody,
     &                in, jn, kn, ispecies, idust,
     &                tgas1d, HI, HII, HeI, HeII, tdust, metallicity, 
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                ndratec, dtemstart, dtemend, h2dusta, 
     &                ncrna, ncrd1a, ncrd2a,
     &                avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k22, k24, k25, k26,
     &                k50, k51, k52, k53, k54, k55,
     &                k56, k13dd, k24shield, k25shield, k26shield,
     &                h2dust, ncrn, ncrd1, ncrd2,
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, itmask)
! -------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!     Arguments

      INTG_PREC is, ie, ijk, iradtype, iradshield, nratec, 
     &        in, jn, kn, ispecies, ithreebody, j, k,
     &        idust, ndratec
      R_PREC temstart, temend, tgas1d(in), dom,
     &        dtemstart, dtemend
      real*8 coolunit, tbase1
      LOGIC_PREC itmask(ijk)

!     Chemistry rates as a function of temperature

      R_PREC k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &     k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &     k54a(nratec), k55a(nratec), k56a(nratec),
     &     k13dda(nratec, 7), h2dusta(nratec, ndratec), 
     &     ncrna(nratec), ncrd1a(nratec), ncrd2a(nratec),
     &     k24, k25, k26,
     &     avgsighp, avgsighep, avgsighe2p, piHI, piHeI

!     Density fields

      R_PREC    HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn)

!     Returned rate values

      R_PREC k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k22(ijk),
     &     k50(ijk), k51(ijk), k52(ijk), k53(ijk), k54(ijk),
     &     k55(ijk), k56(ijk), k13dd(ijk, 7), h2dust(ijk), 
     &     ncrn(ijk), ncrd1(ijk), ncrd2(ijk),
     &     k24shield(ijk), k25shield(ijk), k26shield(ijk)

!     1D temporaries (passed in)

      INTG_PREC indixe(ijk)
      R_PREC t1(ijk), t2(ijk), logtem(ijk), tdef(ijk),
     &     tdust(ijk), metallicity(ijk)

!     1D temporaries (not passed in)

      INTG_PREC d_indixe(ijk)
      R_PREC d_t1(ijk), d_t2(ijk), d_logtem(ijk), d_tdef(ijk),
     &     dusti1(ijk), dusti2(ijk)

!     Parameters

      real*8 everg, e24, e26
      parameter(everg = ev2erg, e24 = 13.6_RKIND, 
     &          e26 = 24.6_RKIND)

!     locals

      INTG_PREC i, n1
      R_PREC factor, x, logtem0, logtem9, dlogtem, nh,
     &     d_logtem0, d_logtem9, d_dlogtem

!     Set log values of start and end of lookup tables

      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem = (log(temend) - log(temstart))/REAL(nratec-1,RKIND)

      do i = is+1, ie+1
         if (itmask(i)) then
!        Compute temp-centered temperature (and log)

!        logtem(i) = log(0.5*(tgas(i)+tgasold(i)))
         logtem(i) = log(tgas1d(i))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

!        Find index into tble and precompute interpolation values

         indixe(i) = min(nratec-1, max(1,
     &        int((logtem(i)-logtem0)/dlogtem,IKIND)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

!        Do linear table lookup (in log temperature)

         k1(i) = k1a(indixe(i)) +
     &           (k1a(indixe(i)+1) -k1a(indixe(i)))*tdef(i)
         k2(i) = k2a(indixe(i)) +
     &           (k2a(indixe(i)+1) -k2a(indixe(i)))*tdef(i)
         k3(i) = k3a(indixe(i)) +
     &           (k3a(indixe(i)+1) -k3a(indixe(i)))*tdef(i)
         k4(i) = k4a(indixe(i)) +
     &           (k4a(indixe(i)+1) -k4a(indixe(i)))*tdef(i)
         k5(i) = k5a(indixe(i)) +
     &           (k5a(indixe(i)+1) -k5a(indixe(i)))*tdef(i)
         k6(i) = k6a(indixe(i)) +
     &           (k6a(indixe(i)+1) -k6a(indixe(i)))*tdef(i)
      endif
      enddo

!     Look-up for 9-species model

      if (ispecies > 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k7(i) = k7a(indixe(i)) +
     &            (k7a(indixe(i)+1) -k7a(indixe(i)))*tdef(i)
            k8(i) = k8a(indixe(i)) +
     &            (k8a(indixe(i)+1) -k8a(indixe(i)))*tdef(i)
            k9(i) = k9a(indixe(i)) +
     &            (k9a(indixe(i)+1) -k9a(indixe(i)))*tdef(i)
            k10(i) = k10a(indixe(i)) +
     &            (k10a(indixe(i)+1) -k10a(indixe(i)))*tdef(i)
            k11(i) = k11a(indixe(i)) +
     &            (k11a(indixe(i)+1) -k11a(indixe(i)))*tdef(i)
            k12(i) = k12a(indixe(i)) +
     &            (k12a(indixe(i)+1) -k12a(indixe(i)))*tdef(i)
            k13(i) = k13a(indixe(i)) +
     &            (k13a(indixe(i)+1) -k13a(indixe(i)))*tdef(i)
            k14(i) = k14a(indixe(i)) +
     &            (k14a(indixe(i)+1) -k14a(indixe(i)))*tdef(i)
            k15(i) = k15a(indixe(i)) +
     &            (k15a(indixe(i)+1) -k15a(indixe(i)))*tdef(i)
            k16(i) = k16a(indixe(i)) +
     &            (k16a(indixe(i)+1) -k16a(indixe(i)))*tdef(i)
            k17(i) = k17a(indixe(i)) +
     &            (k17a(indixe(i)+1) -k17a(indixe(i)))*tdef(i)
            k18(i) = k18a(indixe(i)) +
     &            (k18a(indixe(i)+1) -k18a(indixe(i)))*tdef(i)
            k19(i) = k19a(indixe(i)) +
     &            (k19a(indixe(i)+1) -k19a(indixe(i)))*tdef(i)
            k22(i) = k22a(indixe(i)) +
     &            (k22a(indixe(i)+1) -k22a(indixe(i)))*tdef(i)

!     H2 formation heating terms.

            ncrn(i) = ncrna(indixe(i)) +
     &           (ncrna(indixe(i)+1) -ncrna(indixe(i)))*tdef(i)
            ncrd1(i) = ncrd1a(indixe(i)) +
     &           (ncrd1a(indixe(i)+1) -ncrd1a(indixe(i)))*tdef(i)
            ncrd2(i) = ncrd2a(indixe(i)) +
     &           (ncrd2a(indixe(i)+1) -ncrd2a(indixe(i)))*tdef(i)

         endif
         enddo

         do n1 = 1, 7
            do i = is+1, ie+1
               if (itmask(i)) then
               k13dd(i,n1) = k13dda(indixe(i),n1) +
     &             (k13dda(indixe(i)+1,n1) - 
     &               k13dda(indixe(i)  ,n1) )*tdef(i)
            endif
            enddo
         enddo         

      endif

!     Look-up for 12-species model

      if (ispecies > 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k50(i) = k50a(indixe(i)) +
     &            (k50a(indixe(i)+1) -k50a(indixe(i)))*tdef(i)
            k51(i) = k51a(indixe(i)) +
     &            (k51a(indixe(i)+1) -k51a(indixe(i)))*tdef(i)
            k52(i) = k52a(indixe(i)) +
     &            (k52a(indixe(i)+1) -k52a(indixe(i)))*tdef(i)
            k53(i) = k53a(indixe(i)) +
     &            (k53a(indixe(i)+1) -k53a(indixe(i)))*tdef(i)
            k54(i) = k54a(indixe(i)) +
     &            (k54a(indixe(i)+1) -k54a(indixe(i)))*tdef(i)
            k55(i) = k55a(indixe(i)) +
     &            (k55a(indixe(i)+1) -k55a(indixe(i)))*tdef(i)
            k56(i) = k56a(indixe(i)) +
     &            (k56a(indixe(i)+1) -k56a(indixe(i)))*tdef(i)
         endif
         enddo
      endif

!     Look-up for H2 formation on dust

      if (idust > 0) then

         d_logtem0 = log(dtemstart)
         d_logtem9 = log(dtemend)
         d_dlogtem = (log(dtemend) - log(dtemstart))/real(ndratec-1)

         do i = is+1, ie+1
            if (itmask(i)) then

!              Assume dust melting at T > 1500 K

               if (tdust(i) > dtemend) then
                  h2dust(i) = tiny
               else

!                 Get log dust temperature

                  d_logtem(i) = log(tdust(i))
                  d_logtem(i) = max(d_logtem(i), d_logtem0)
                  d_logtem(i) = min(d_logtem(i), d_logtem9)

!                 Find index into table and precompute interpolation values

                  d_indixe(i) = min(ndratec-1, max(1,
     &                 int((d_logtem(i)-d_logtem0)/d_dlogtem,IKIND)+1))
                  d_t1(i) = (d_logtem0 + (d_indixe(i) - 1)*d_dlogtem)
                  d_t2(i) = (d_logtem0 + (d_indixe(i)    )*d_dlogtem)
                  d_tdef(i) = (d_logtem(i) - d_t1(i)) / 
     &                 (d_t2(i) - d_t1(i))

!                 Get rate from 2D interpolation

                  dusti1(i) = h2dusta(indixe(i), d_indixe(i)) +
     &                 (h2dusta(indixe(i)+1, d_indixe(i)) - 
     &                 h2dusta(indixe(i),   d_indixe(i)))*tdef(i)
                  dusti2(i) = h2dusta(indixe(i), d_indixe(i)+1) +
     &                 (h2dusta(indixe(i)+1, d_indixe(i)+1) - 
     &                 h2dusta(indixe(i),   d_indixe(i)+1))*tdef(i)
                  h2dust(i) = dusti1(i) + 
     &                 (dusti2(i) - dusti1(i))*d_tdef(i)

!                 Multiply by metallicity

                  h2dust(i) = h2dust(i) * metallicity(i)

               endif
            endif
         enddo
      endif

!        Include approximate self-shielding factors if requested

#ifdef RADIATION
      do i = is+1, ie+1
         if (itmask(i)) then
         k24shield(i) = k24
         k25shield(i) = k25
         k26shield(i) = k26
      endif
      enddo
      if (iradshield == 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k24shield(i) = k24shield(i)*exp(-HI(i,j,k)*avgsighp*dom)
            k25shield(i) = k25shield(i)*exp(-HeII(i,j,k)*avgsighe2p*dom)
            k26shield(i) = k26shield(i)*exp(-HeI(i,j,k)*avgsighep*dom)
         endif
         enddo
      endif

!        If using a high-energy radiation field, then account for
!          effects of secondary elections (Shull * Steenberg 1985)
!          (see calc_rate.src)

      if (iradtype == 8) then
         do i = is+1, ie+1
            if (itmask(i)) then
            x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.e-4_RKIND)
            factor = 0.3908_RKIND*(1._RKIND - 
     &           x**0.4092_RKIND)**1.7592_RKIND
            k24shield(i) = k24shield(i) + factor*(piHI + 
     &           0.08_RKIND*piHeI)/(e24*everg) *coolunit*tbase1
            factor = 0.0554_RKIND*(1._RKIND - 
     &           x**0.4614_RKIND)**1.666_RKIND
            k26shield(i) = k26shield(i) + factor*(piHI/0.08_RKIND + 
     &           piHeI)/(e26*everg) *coolunit*tbase1
         endif
         enddo
      endif

#endif /* RADIATION */

!           If using H2, and using the density-dependent collisional
!             H2 dissociation rate, then replace the the density-independant
!                k13 rate with the new one.
!         May/00: there appears to be a problem with the density-dependent
!             collisional rates.  Currently turned off until further notice.

#define USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
#ifdef USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
            if (ispecies > 1 .and. ithreebody == 0) then
               do i = is+1, ie+1
                  if (itmask(i)) then
                  nh = min(HI(i,j,k)*dom, 1.e9_RKIND)
                  k13(i) = tiny
                  if (tgas1d(i) >= 500.0 .and. tgas1d(i) < 1.e6) then
                     k13(i) = k13dd(i,1)-k13dd(i,2)/
     &                          (1._RKIND+(nh/k13dd(i,5))**k13dd(i,7))
     &                     + k13dd(i,3)-k13dd(i,4)/
     &                          (1._RKIND+(nh/k13dd(i,6))**k13dd(i,7))
                     k13(i) = max(10._RKIND**k13(i), tiny)
                  endif
               endif
               enddo
            endif
#endif /*  USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE */

      return
      end

! -------------------------------------------------------------------
!  This routine calculates the electron and HI rates of change in
!    order to determine the maximum permitted timestep

      subroutine rate_timestep(dedot, HIdot, ispecies, idust,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55, k56, 
     &                     h2dust, ncrn, ncrd1, ncrd2, rhoH, 
     &                     k24shield, k25shield, k26shield,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, 
     &                     kdissH2I, itmask, edot, chunit, dom)

! -------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!     arguments

      INTG_PREC ispecies, idust, is, ie, j, k, in, jn, kn, 
     &     iradtrans, irt_honly
      R_PREC dedot(in), HIdot(in), dom
      real*8 edot(in)
      LOGIC_PREC itmask(in)

!     Density fields

      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)

!     Radiation fields

      R_PREC    kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn)

      real*8 chunit

!     Rate values

      R_PREC k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &     k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &     k11(in), k12(in), k13(in), k14(in), k15(in),
     &     k16(in), k17(in), k18(in), k19(in), k22(in),
     &     k50(in), k51(in), k52(in), k53(in), k54(in),
     &     k55(in), k56(in), h2dust(in), 
     &     ncrn(in), ncrd1(in), ncrd2(in), rhoH(in), 
     &     k24shield(in), k25shield(in), k26shield(in),
     &     k24, k25, k26, k27, k28, k29, k30, k31

!     locals

      INTG_PREC i
      real*8 h2heatfac(in), H2delta(in), H2dmag, atten, tau

      if (ispecies == 1) then

         do i = is+1, ie+1
            if (itmask(i)) then
!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k3(i)*HeI(i,j,k)*de(i,j,k)/4._RKIND
     &               + k5(i)*HeII(i,j,k)*de(i,j,k)/4._RKIND
     &               - k2(i)*HII(i,j,k)*de(i,j,k)
     &               - k4(i)*HeII(i,j,k)*de(i,j,k)/4._RKIND
     &               - k6(i)*HeIII(i,j,k)*de(i,j,k)/4._RKIND
#ifdef RADIATION
     &               +      ( k24shield(i)*HI(i,j,k)
     &               + k25shield(i)*HeII(i,j,k)/4._RKIND
     &               + k26shield(i)*HeI(i,j,k)/4._RKIND)
#endif /* RADIATION */

!     Compute the HI density rate-of-change

            HIdot(i) =
     &               - k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k2(i)*HII(i,j,k)*de(i,j,k)
#ifdef RADIATION
     &               -      k24shield(i)*HI(i,j,k)
#endif /* RADIATION */

         endif                  ! itmask
         enddo
      else

!         Include molecular hydrogen rates for HIdot

         do i = is+1, ie+1
            if (itmask(i)) then
               HIdot(i) = 
     &               -      k1(i) *de(i,j,k)    *HI(i,j,k)  
     &               -      k7(i) *de(i,j,k)    *HI(i,j,k)
     &               -      k8(i) *HM(i,j,k)    *HI(i,j,k)
     &               -      k9(i) *HII(i,j,k)   *HI(i,j,k)
     &               -      k10(i)*H2II(i,j,k)  *HI(i,j,k)/2._RKIND
     &               - 2._RKIND*k22(i)*HI(i,j,k)**2 *HI(i,j,k)
     &               +      k2(i) *HII(i,j,k)   *de(i,j,k) 
     &               + 2._RKIND*k13(i)*HI(i,j,k)    *H2I(i,j,k)/2._RKIND
     &               +      k11(i)*HII(i,j,k)   *H2I(i,j,k)/2._RKIND
     &               + 2._RKIND*k12(i)*de(i,j,k)    *H2I(i,j,k)/2._RKIND
     &               +      k14(i)*HM(i,j,k)    *de(i,j,k)
     &               +      k15(i)*HM(i,j,k)    *HI(i,j,k)
     &               + 2._RKIND*k16(i)*HM(i,j,k)    *HII(i,j,k)
     &               + 2._RKIND*k18(i)*H2II(i,j,k)  *de(i,j,k)/2._RKIND
     &               +      k19(i)*H2II(i,j,k)  *HM(i,j,k)/2._RKIND
#ifdef RADIATION
     &               -      k24shield(i)*HI(i,j,k)
#endif /* RADIATION */

!     Add H2 formation on dust grains

            if (idust > 0) then
               HIdot(i) = HIdot(i) 
     &              - 2.d0 * h2dust(i) * rhoH(i)
            endif

!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i) * HI(i,j,k)   * de(i,j,k)
     &               + k3(i) * HeI(i,j,k)  * de(i,j,k)/4._RKIND
     &               + k5(i) * HeII(i,j,k) * de(i,j,k)/4._RKIND
     &               + k8(i) * HM(i,j,k)   * HI(i,j,k)
     &               + k15(i)* HM(i,j,k)   * HI(i,j,k)
     &               + k17(i)* HM(i,j,k)   * HII(i,j,k)
     &               + k14(i)* HM(i,j,k)   * de(i,j,k)
     &               - k2(i) * HII(i,j,k)  * de(i,j,k)
     &               - k4(i) * HeII(i,j,k) * de(i,j,k)/4._RKIND
     &               - k6(i) * HeIII(i,j,k)* de(i,j,k)/4._RKIND
     &               - k7(i) * HI(i,j,k)   * de(i,j,k)
     &               - k18(i)* H2II(i,j,k) * de(i,j,k)/2._RKIND

#ifdef RADIATION
     &               + (k24shield(i)*HI(i,j,k)
     &               +  k25shield(i)*HeII(i,j,k)/4.e0
     &               +  k26shield(i)*HeI(i,j,k)/4.e0)
#endif /* RADIATION */

!     H2 formation heating

!     Equation 23 from Omukai (2000)
            h2heatfac(i) = (1._RKIND + (ncrn(i) / (dom *
     &           (HI(i,j,k) * ncrd1(i) +
     &           H2I(i,j,k) * 0.5_RKIND* ncrd2(i)))))**(-1._RKIND)

            H2delta(i) = 
     &          HI(i,j,k) *
     &           (  4.48_RKIND * k22(i)*HI(i,j,k)**2._RKIND
     &            - 4.48_RKIND * k13(i) * H2I(i,j,k)/2._RKIND)

            ! We only want to apply this if the formation dominates, but we
            ! need to apply it outside the delta calculation.
            if(H2delta(i) > 0.0) then
              H2delta(i) = H2delta(i) * h2heatfac(i)
            endif

            if (idust > 0) then
               H2delta(i) = H2delta(i) + 
     &              h2dust(i) * HI(i,j,k) * rhoH(i) * 
     &              (0.2_RKIND + 4.2_RKIND * h2heatfac(i))
            endif

!            H2dmag = abs(H2delta)/(
!     &          HI(i,j,k)*( k22(i) * HI(i,j,k)**2._RKIND
!     &                    + k13(i) * H2I(i,j,k)/2._RKIND))
!            tau = (H2dmag/1.e-5)**-1._RKIND
!            tau = max(tau, 1.e-5)
!            atten = min((1.-exp(-tau))/tau,1.0)
            atten = 1._RKIND
            edot(i) = edot(i) + chunit * H2delta(i) * atten
!     &       + H2I(i,j,k)*( k21(i) * HI(i,j,k)**2._RKIND
!     &                    - k23(i) * H2I(i,j,k))
!H * (k22 * H^2 - k13 * H_2) + H_2 * (k21 * H^2 - k23 * H_2) */
         endif                  ! itmask
         enddo
      endif

!     Add photo-ionization rates if needed

#ifdef RADIATION
      if (iradtrans == 1) then
         if (irt_honly == 0) then
            do i = is+1, ie+1
               if (itmask(i)) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
     &                 + kphHeI(i,j,k) * HeI(i,j,k) / 4._RKIND
     &                 + kphHeII(i,j,k) *HeII(i,j,k) / 4._RKIND
               endif
            enddo
         else
            do i = is+1, ie+1
               if (itmask(i)) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
               endif
            enddo
         endif
      endif
#endif /* RADIATION */

      return
      end


! -----------------------------------------------------------
!  This routine uses one linearly implicit Gauss-Seidel sweep of 
!   a backward-Euler time integrator to advance the rate equations 
!   by one (sub-)cycle (dtit).

      subroutine step_rate(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, dtit,
     &                     in, jn, kn, is, ie, j, k, ispecies, idust, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55, k56, 
     &                     h2dust, rhoH,
     &                     k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
     &                     dedot_prev, HIdot_prev,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII,
     &                     kdissH2I, itmask)
c -------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!     arguments

      INTG_PREC ispecies, idust, in, jn, kn, is, ie, j, k, 
     &     iradtrans, irt_honly
      R_PREC    dtit(in), dedot_prev(in), HIdot_prev(in)
      LOGIC_PREC itmask(in)

!     Density fields

      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)

!     Radiation fields

      R_PREC    kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn)

!     Rate values

      R_PREC k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &     k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &     k11(in), k12(in), k13(in), k14(in), k15(in),
     &     k16(in), k17(in), k18(in), k19(in), k22(in),
     &     k50(in), k51(in), k52(in), k53(in), k54(in),
     &     k55(in), k56(in), h2dust(in), rhoH(in),
     &     k24shield(in), k25shield(in), k26shield(in),
     &     k24, k25, k26, k27, k28, k29, k30, k31

!     temporaries (passed in)

      R_PREC HIp(in), HIIp(in), HeIp(in), HeIIp(in), HeIIIp(in),
     &     HMp(in), H2Ip(in), H2IIp(in), dep(in),
     &     DIp(in), DIIp(in), HDIp(in)

!     locals

      INTG_PREC i
      R_PREC scoef, acoef

!   A) the 6-species integrator
!      
      if (ispecies == 1) then

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI

            scoef  = k2(i)*HII(i,j,k)*de(i,j,k)
            acoef  = k1(i)*de(i,j,k)
#ifdef RADIATION
     &             + k24shield(i)
            if (iradtrans == 1) acoef = acoef + kphHI(i,j,k)
#endif /* RADIATION */
            HIp(i)  = (scoef*dtit(i) + HI(i,j,k)) / 
     $           (1._RKIND + acoef*dtit(i))
            if (HIp(i) .ne. HIp(i)) then
               write(*,*) '[0]HUGE HIp! :: ', i, j, k, HIp(i), HI(i,j,k),
     $              HII(i,j,k), de(i,j,k), kphHI(i,j,k), 
     $              scoef, acoef, dtit(i)
               ERROR_MESSAGE
            endif

!        2) HII
c 
            scoef  = k1(i)*HIp(i)*de(i,j,k)
#ifdef RADIATION
     &                   + k24shield(i)*HIp(i)
            if (iradtrans == 1) 
     &           scoef = scoef + kphHI(i,j,k)*HIp(i)
#endif /* RADIATION */
            acoef  = k2(i)*de (i,j,k)
            HIIp(i) = (scoef*dtit(i) + HII(i,j,k)) / 
     $           (1._RKIND +acoef*dtit(i))
!
            if (HIIp(i) <= 0.0) then  !#####
               write(*,*) 'negative HIIp! :: ', i, j, k, HIIp(i), 
     $              scoef, dtit(i), HII(i,j,k), acoef,
     $              k2(i), de(i,j,k),
     $              kphHI(i,j,k), HIp(i),
     $              k24shield(i)
            endif

!        3) Electron density

            scoef = 0._RKIND
#ifdef RADIATION
     &                 + k24shield(i)*HI(i,j,k)
     &                 + k25shield(i)*HeII(i,j,k)/4._RKIND
     &                 + k26shield(i)*HeI(i,j,k)/4._RKIND
            if (iradtrans == 1 .and. irt_honly == 0) 
     &           scoef = scoef + kphHI(i,j,k)   * HI(i,j,k)
     &           + kphHeI(i,j,k)  * HeI(i,j,k)/4._RKIND
     &           + kphHeII(i,j,k) * HeII(i,j,k)/4._RKIND
            if (iradtrans == 1 .and. irt_honly == 1) 
     &           scoef = scoef + kphHI(i,j,k)   * HI(i,j,k)
#endif /* RADIATION */
            acoef = -(k1(i)*HI(i,j,k)      - k2(i)*HII(i,j,k)
     &              + k3(i)*HeI(i,j,k)/4._RKIND  
     &              - k6(i)*HeIII(i,j,k)/4._RKIND
     &              + k5(i)*HeII(i,j,k)/4._RKIND 
     &              - k4(i)*HeII(i,j,k)/4._RKIND)
            dep(i)   = (scoef*dtit(i) + de(i,j,k))
     &                     / (1._RKIND + acoef*dtit(i))

         endif                  ! itmask
         enddo

      endif                     ! (ispecies == 1)

!  --- (B) Do helium chemistry in any case: (for all ispecies values) ---

      do i = is+1, ie+1
         if (itmask(i)) then

!        4) HeI

         scoef  = k4(i)*HeII(i,j,k)*de(i,j,k)
         acoef  = k3(i)*de(i,j,k)
#ifdef RADIATION
     &                + k26shield(i)
         if (iradtrans == 1 .and. irt_honly == 0) 
     $        acoef = acoef + kphHeI(i,j,k)
#endif /* RADIATION */
         HeIp(i)   = ( scoef*dtit(i) + HeI(i,j,k) ) 
     &              / ( 1._RKIND + acoef*dtit(i) )

!        5) HeII

         scoef  = k3(i)*HeIp(i)*de(i,j,k)
     &          + k6(i)*HeIII(i,j,k)*de(i,j,k)
#ifdef RADIATION
     &          + k26shield(i)*HeIp(i)
         if (iradtrans == 1 .and. irt_honly == 0) 
     $        scoef = scoef + kphHeI(i,j,k)*HeIp(i)
#endif /* RADIATION */
         acoef  = k4(i)*de(i,j,k) + k5(i)*de(i,j,k)
#ifdef RADIATION
     &          + k25shield(i)
         if (iradtrans == 1 .and. irt_honly == 0) 
     $        acoef = acoef + kphHeII(i,j,k)
#endif /* RADIATION */
         HeIIp(i)  = ( scoef*dtit(i) + HeII(i,j,k) )
     &              / ( 1._RKIND + acoef*dtit(i) )

!       6) HeIII

         scoef   = k5(i)*HeIIp(i)*de(i,j,k)
#ifdef RADIATION
     &           + k25shield(i)*HeIIp(i)
         if (iradtrans == 1 .and. irt_honly == 0) 
     $        scoef = scoef + kphHeII(i,j,k) * HeIIp(i)
#endif /* RADIATION */
         acoef   = k6(i)*de(i,j,k)
         HeIIIp(i)  = ( scoef*dtit(i) + HeIII(i,j,k) )
     &                / ( 1._RKIND + acoef*dtit(i) )

      endif                     ! itmask
      enddo

c --- (C) Now do extra 3-species for molecular hydrogen ---

      if (ispecies > 1) then

!        First, do HI/HII with molecular hydrogen terms

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI
!     
            scoef  =      k2(i) * HII(i,j,k) * de(i,j,k) 
     &             + 2._RKIND*k13(i)* HI(i,j,k)  * H2I(i,j,k)/2._RKIND
     &             +      k11(i)* HII(i,j,k) * H2I(i,j,k)/2._RKIND
     &             + 2._RKIND*k12(i)* de(i,j,k)  * H2I(i,j,k)/2._RKIND
     &             +      k14(i)* HM(i,j,k)  * de(i,j,k)
     &             +      k15(i)* HM(i,j,k)  * HI(i,j,k)
     &             + 2._RKIND*k16(i)* HM(i,j,k)  * HII(i,j,k)
     &             + 2._RKIND*k18(i)* H2II(i,j,k)* de(i,j,k)/2._RKIND
     &             +      k19(i)* H2II(i,j,k)* HM(i,j,k)/2._RKIND
#ifdef RADIATION
     &             + 2._RKIND*k31   * H2I(i,j,k)/2._RKIND
            if (iradtrans == 1)
     &           scoef = scoef + 
     &           2._RKIND*kdissH2I(i,j,k)* H2I(i,j,k)/2._RKIND
#endif /* RADIATION */
            acoef  =      k1(i) * de(i,j,k)
     &             +      k7(i) * de(i,j,k)  
     &             +      k8(i) * HM(i,j,k)
     &             +      k9(i) * HII(i,j,k)
     &             +      k10(i)* H2II(i,j,k)/2._RKIND
     &             + 2._RKIND*k22(i)* HI(i,j,k)**2
#ifdef RADIATION
     &             + k24shield(i)
            if (iradtrans == 1) acoef = acoef + kphHI(i,j,k)
#endif /* RADIATION */

            if (idust > 0) then
               acoef = acoef + 2.d0 * h2dust(i) * rhoH(i)
            endif

            HIp(i)  = ( scoef*dtit(i) + HI(i,j,k) ) / 
     &                      ( 1._RKIND + acoef*dtit(i) )
            if (HIp(i) .ne. HIp(i)) then
               write(*,*) '[1]HUGE HIp! :: ', i, j, k, HIp(i), HI(i,j,k),
     $              HII(i,j,k), de(i,j,k), H2I(i,j,k), kphHI(i,j,k)
               write(*,*) '[2]HUGE HIp! :: ', i, j, k, 
     $              scoef, acoef, dtit(i),
     $              k2(i), k11(i), k12(i), k13(i), k14(i), k15(i),
     $              k1(i), k7(i), k8(i), k9(i), k10(i), k22(i)
               ERROR_MESSAGE
            endif

!          2) HII

            scoef  =    k1(i)  * HI(i,j,k) * de(i,j,k)
     &             +    k10(i) * H2II(i,j,k)*HI(i,j,k)/2._RKIND
#ifdef RADIATION
     &             + k24shield(i)*HI(i,j,k)
            if (iradtrans == 1) 
     &           scoef = scoef + kphHI(i,j,k)*HI(i,j,k)
#endif /* RADIATION */
            acoef  =    k2(i)  * de(i,j,k)
     &             +    k9(i)  * HI(i,j,k)
     &             +    k11(i) * H2I(i,j,k)/2._RKIND
     &             +    k16(i) * HM(i,j,k)
     &             +    k17(i) * HM(i,j,k)
            HIIp(i)   = ( scoef*dtit(i) + HII(i,j,k) )
     &                      / ( 1._RKIND + acoef*dtit(i) )
!     
!          3) electrons:

            scoef =   k8(i) * HM(i,j,k) * HI(i,j,k)
     &             +  k15(i)* HM(i,j,k) * HI(i,j,k)
     &             +  k17(i)* HM(i,j,k) * HII(i,j,k)
!                  
#ifdef RADIATION
     &             + k24shield(i)*HIp(i)
     &             + k25shield(i)*HeIIp(i)/4._RKIND
     &             + k26shield(i)*HeIp(i)/4._RKIND
            if (iradtrans == 1 .and. irt_honly == 0) 
     &           scoef = scoef + kphHI(i,j,k)   * HIp(i)
     &           + kphHeI(i,j,k)  * HeIp(i)/4._RKIND
     &           + kphHeII(i,j,k) * HeIIp(i)/4._RKIND
            if (iradtrans == 1 .and. irt_honly == 1) 
     &           scoef = scoef + kphHI(i,j,k)   * HIp(i)
#endif /* RADIATION */
            acoef = - (k1(i) *HI(i,j,k)    - k2(i)*HII(i,j,k)
     &              +  k3(i) *HeI(i,j,k)/4._RKIND 
     &              - k6(i)*HeIII(i,j,k)/4._RKIND
     &              +  k5(i) *HeII(i,j,k)/4._RKIND
     &              - k4(i)*HeII(i,j,k)/4._RKIND
     &              +  k14(i)*HM(i,j,k)
     &              -  k7(i) *HI(i,j,k)
     &              -  k18(i)*H2II(i,j,k)/2._RKIND)
            dep(i)  = ( scoef*dtit(i) + de(i,j,k) )
     &                / ( 1._RKIND + acoef*dtit(i) )

!           7) H2

            scoef = 2._RKIND*(k8(i)  * HM(i,j,k)   * HI(i,j,k)
     &            +       k10(i) * H2II(i,j,k) * HI(i,j,k)/2._RKIND
     &            +       k19(i) * H2II(i,j,k) * HM(i,j,k)/2._RKIND
     &            +       k22(i) * HI(i,j,k) * (HI(i,j,k))**2._RKIND)
            acoef = ( k13(i)*HI(i,j,k) + k11(i)*HII(i,j,k)
     &              + k12(i)*de(i,j,k) )
#ifdef RADIATION
     &              + k29 + k31
            if (iradtrans == 1) acoef = acoef + kdissH2I(i,j,k)
#endif /* RADIATION */

            if (idust > 0) then
               scoef = scoef + 2.d0 * h2dust(i) * HI(i,j,k) * rhoH(i)
            endif

            H2Ip(i) = ( scoef*dtit(i) + H2I(i,j,k) )
     &                / ( 1._RKIND + acoef*dtit(i) )

!           8) H-

            HMp(i) = ( k7(i)*HIp(i)*dep(i) )
     &             / ( (k8(i)+k15(i))*HIp(i)
     &             + ( k16(i)+k17(i))*HIIp(i)+k14(i)*dep(i)
#ifdef RADIATION
     &             + k27
#endif /* RADIATION */
     &                 )

!           9) H2+

            H2IIp(i) = 2._RKIND*( k9 (i)*HIp(i)*HIIp(i)
     &                    +   k11(i)*H2Ip(i)/2._RKIND*HIIp(i)
     &                    +   k17(i)*HMp(i)*HIIp(i)
#ifdef RADIATION
     &                    + k29*H2Ip(i)
#endif /* RADIATION */
     &                    )
     &                 /  ( k10(i)*HIp(i) + k18(i)*dep(i)
     &                    + k19(i)*HMp(i)
#ifdef RADIATION
     &                    + (k28+k30)
#endif /* RADIATION */
     &                    )

         endif                  ! itmask
         enddo
!     
      endif                     ! H2

!  --- (D) Now do extra 3-species for molecular HD ---
!     
      if (ispecies > 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
!     
!         1) DI
!     
            scoef =   (       k2(i) * DII(i,j,k) * de(i,j,k)
     &                 +      k51(i)* DII(i,j,k) * HI(i,j,k)
     &                 + 2._RKIND*k55(i)* HDI(i,j,k) *HI(i,j,k)/3._RKIND
     &                 )
            acoef  =    k1(i) * de(i,j,k)
     &             +    k50(i) * HII(i,j,k)
     &             +    k54(i) * H2I(i,j,k)/2._RKIND
     &             +    k56(i) * HM(i,j,k)
#ifdef RADIATION
     &             + k24shield(i)
            if (iradtrans == 1) acoef = acoef + kphHI(i,j,k)
#endif /* RADIATION */
            DIp(i)    = ( scoef*dtit(i) + DI(i,j,k) ) / 
     &                  ( 1._RKIND + acoef*dtit(i) )

!         2) DII
c 
            scoef =   (   k1(i)  * DI(i,j,k) * de(i,j,k)
     &            +       k50(i) * HII(i,j,k)* DI(i,j,k)
     &            +  2._RKIND*k53(i) * HII(i,j,k)* HDI(i,j,k)/3._RKIND
     &            )
#ifdef RADIATION
     &            + k24shield(i)*DI(i,j,k)
            if (iradtrans == 1) 
     &           scoef = scoef + kphHI(i,j,k)*DI(i,j,k)
#endif /* RADIATION */
            acoef =    k2(i)  * de(i,j,k)
     &            +    k51(i) * HI(i,j,k)
     &            +    k52(i) * H2I(i,j,k)/2._RKIND

            DIIp(i)   = ( scoef*dtit(i) + DII(i,j,k) )
     &                 / ( 1._RKIND + acoef*dtit(i) )

!          3) HDI
c 
            scoef = 3._RKIND*(k52(i) * DII(i,j,k)* H2I(i,j,k)/4._RKIND
     &             +      k54(i) * DI(i,j,k) * H2I(i,j,k)/4._RKIND
     &             + 2._RKIND*k56(i) * DI(i,j,k) * HM(i,j,k)/2._RKIND
     &                 )
            acoef  =    k53(i) * HII(i,j,k)
     &             +    k55(i) * HI(i,j,k)

            HDIp(i)   = ( scoef*dtit(i) + HDI(i,j,k) )
     &                 / ( 1._RKIND + acoef*dtit(i) )

         endif                  ! itmask
         enddo
      endif

!   --- (E) Set densities from 1D temps to 3D fields ---

      do i = is+1, ie+1
         if (itmask(i)) then
         HIdot_prev(i) = abs(HI(i,j,k)-HIp(i))/max(dtit(i), tiny)
         HI(i,j,k)    = max(HIp(i), tiny)
         HII(i,j,k)   = max(HIIp(i), tiny)
         HeI(i,j,k)   = max(HeIp(i), tiny)
         HeII(i,j,k)  = max(HeIIp(i), tiny)
         HeIII(i,j,k) = max(HeIIIp(i), 1.e-5_RKIND*tiny)

!        de(i,j,k)    = dep(i)

!        Use charge conservation to determine electron fraction

         dedot_prev(i) = de(i,j,k)
         de(i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND + 
     &        HeIII(i,j,k)/2._RKIND
         if (ispecies > 1) 
     &      de(i,j,k) = de(i,j,k) - HM(i,j,k) + H2II(i,j,k)/2._RKIND
         dedot_prev(i) = abs(de(i,j,k)-dedot_prev(i))/
     &                         max(dtit(i),tiny)

         if (ispecies > 1) then
            HM(i,j,k)    = max(HMp(i), tiny)
            H2I(i,j,k)   = max(H2Ip(i), tiny)
            H2II(i,j,k)  = max(H2IIp(i), tiny)
         endif

         if (ispecies > 2) then
            DI(i,j,k)    = max(DIp(i), tiny)
            DII(i,j,k)   = max(DIIp(i), tiny)
            HDI(i,j,k)   = max(HDIp(i), tiny)
         endif
      endif                     ! itmask
!     

      if (HI(i,j,k) .ne. HI(i,j,k)) then
         write(*,*) 'HUGE HI! :: ', i, j, k, HI(i,j,k)
      endif

      enddo                     ! end loop over i

      return
      end

! ------------------------------------------------------------------
!   This routine correct the highest abundence species to
!     insure conservation of particle number and charge.

      subroutine make_consistent(de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal, d,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, fh, dtoh)
! -------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!     Arguments

      INTG_PREC in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn), metal(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)  
      R_PREC    fh, dtoh

!     Parameters

      INTG_PREC ijk
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)

!     locals

      INTG_PREC i, j, k
      R_PREC totalH(ijk), totalHe(ijk), 
     &     totalD, correctH, correctHe, correctD
      R_PREC metalfree(ijk)

!     Loop over all zones

      do k = ks+1, ke+1
      do j = js+1, je+1

!     Compute total densities of H and He
!         (ensure non-negativity)

      if (imetal == 1) then
         do i = is+1, ie+1
            metalfree(i) = d(i,j,k) - metal(i,j,k)
         enddo
      else
         do i = is+1, ie+1
            metalfree(i) = d(i,j,k)
         enddo
      endif

      do i = is+1, ie+1
         HI   (i,j,k) = abs(HI   (i,j,k))
         HII  (i,j,k) = abs(HII  (i,j,k))
         HeI  (i,j,k) = abs(HeI  (i,j,k))
         HeII (i,j,k) = abs(HeII (i,j,k))
         HeIII(i,j,k) = abs(HeIII(i,j,k))
         totalH(i) = HI(i,j,k) + HII(i,j,k)
         totalHe(i) = HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k)
      enddo

!     include molecular hydrogen

      if (ispecies > 1) then
         do i = is+1, ie+1
            HM   (i,j,k) = abs(HM   (i,j,k))
            H2II (i,j,k) = abs(H2II (i,j,k))
            H2I  (i,j,k) = abs(H2I  (i,j,k))
            totalH(i) = totalH(i) + HM(i,j,k) + H2I(i,j,k) + H2II(i,j,k)
         enddo
      endif

!     Correct densities by keeping fractions the same

      do i = is+1, ie+1
         correctH = fh*metalfree(i)/totalH(i)
         HI(i,j,k)  = HI(i,j,k)*correctH
         HII(i,j,k) = HII(i,j,k)*correctH

         correctHe = (1._RKIND - fh)*metalfree(i)/totalHe(i)
         HeI(i,j,k)   = HeI(i,j,k)*correctHe
         HeII(i,j,k)  = HeII(i,j,k)*correctHe
         HeIII(i,j,k) = HeIII(i,j,k)*correctHe

!     Correct molecular hydrogen-related fractions

         if (ispecies > 1) then
            HM   (i,j,k) = HM(i,j,k)*correctH
            H2II (i,j,k) = H2II(i,j,k)*correctH
            H2I  (i,j,k) = H2I(i,j,k)*correctH
         endif
      enddo

!     Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (ispecies > 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) + 
     &           2._RKIND/3._RKIND*HDI(i,j,k)
            correctD = fh*dtoh*metalfree(i)/totalD
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
         enddo
      endif

!       Set the electron density

      do i = is+1, ie+1
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND + 
     &        HeIII(i,j,k)/2._RKIND
         if (ispecies > 1) de(i,j,k) = de(i,j,k)
     &             - HM(i,j,k) + H2II(i,j,k)/2._RKIND
      enddo

      enddo  ! end loop over j
      enddo  ! end loop over k

      return
      end

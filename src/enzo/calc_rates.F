#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c//////////////////////  SUBROUTINE CALC_RATES  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine calc_rates(
     &             nratec, aye, temstart, temend, alpha0, f3, 
     &             iradtype, casebrates, threebody,
     &             utem, uxyz, uaye, urho, utim,
     &             ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &             ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &             reHeII2a, reHeIIIa, brema, compa, 
     &             gammahacgs, gammaha,
     &             piHI, piHeI, piHeII,
     &             hyd01ka, h2k01a, vibha, rotha, rotla, 
     &             gpldl, gphdl, hdltea, hdlowa, hdcool, ciecoa,
     &             gaHIa, gaH2a, gaHea, gaHpa, gaela, gasgra,
     &             k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &             k11a, k12a, k13a, k13dda, k14a, k15a, k16a, k17a,
     &             k18a, k19a, k20a, k21a, k22a, k23a,
     &             k24, k25, k26, k27, k28, k29, k30, k31, 
     &             k50a, k51a, k52a, k53a, k54a, k55a, k56a, 
     &             ndratec, dtemstart, dtemend, h2dusta, 
     &             ncrna, ncrd1a, ncrd2a, ioutput
     &                     )
c
c  COMPUTE MULTISPECIES RATE LOOKUP TABLE
c
c  written by: Yu Zhang, Peter Anninos and Tom Abel
c  date:       
c  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
c  modified2: October, 1996 by GB; moved to AMR
c  modified3: February, 2000 by GB; added Toms new collisional rates
c  modified4: July, 2010 by Dan Reynolds; added case-B recombination rates
c
c  PURPOSE:
c    Construct tables for rate coefficients and some cooling functions 
c    versus the logarithm of temperature.
c
c  UNITS:
c    Most rate coefficients are computed in cgs units and then converted
c    into more convenient units.  This would usually be code units
c    (see units.src) but in many cases the natural code units are not
c    near unity (e.g. particles/box volume), so we include units and
c    constants from the equations themselves.  This is explained for
c    each of two types used: rate coefficients, cooling coefficients.
c    Note also that the densities in the rate equations are true, not
c    comoving densities and so we must include a factor of a^3.
c    NOTE: if the spectrum changes, this routine must be called again.
c
c  PARAMETERS:
c    nnu      = number of frequency bins
c    everg    = ergs per eV (Rydberg constant) in cgs
c    evhz     = everg / h
c    mh       = mass of hydrogen atom
c    temstart = table start temperature (K)
c    temend   = table end temperature (K)
c    dtemstart = dust table start temperature (K)
c    dtemend  = dust table end temperature (K)
c    tevk     = conversion constant [eV] = [K] / tevk
c    ireco    = Recombination cooling flag (1 = on)
c    ibrco    = Brehmmstrahlung cooling flag (1 = on)
c    icico    = Collisional ionization flag (1 = on)
c    iceco    = Collisional excitation flag (1 = on)
c    iphoto_mol = Molecular photo-dissociation flag (1 = on)
c
c  INPUTS
c    nratec   = number of temperature bins in rate equations
c    ndratec  = number of dust temperature bins for H2+dust rate
c
c    tbase1   = code time units
c    xbase1   = code length units
c    kunit    = conversion factor xbase1**3/tbase1
c
c    alpha0   = power law slope of radiation flux (?)
c    f3       = amplitude of external flux at z = 3
c    iradtype = radiation field type (0 - none,
c                8 - hard spectrum w/ absorption,
c                otherwise ignored)
c    casebrates = use case-B recombination rates for k2, k4 and k6
c
c  RATES:
c
c  (New rates numbering as of Feb/2000, agrees with original
c   Abel etal 1997)
c    old   -> new              old   -> new
c    k1-10 -> k1-10            k16   -> k14
c    k11   -> k13              k17   -> k15
c    k12   -> k11              k18   -> k16
c    k13   -> removed          k19   -> k17
c    k14   -> k12              k20   -> k18
c    k15   -> removed          k21   -> k19
c
C ------- 1:    HI    + e   -> HII   + 2e
C ------- 2:    HII   + e   -> H     + p
C ------- 3:    HeI   + e   -> HeII  + 2e
C ------- 4:    HeII  + e   -> HeI   + p
C ------- 5:    HeII  + e   -> HeIII + 2e
C ------- 6:    HeIII + e   -> HeII  + p
C ------- 7:    HI    + e   -> HM    + p
C ------- 8:    HM    + HI  -> H2I*  + e
C ------- 9:    HI    + HII -> H2II  + p
C ------- 10:   H2II  + HI  -> H2I*  + HII
c
c
C --------------  old  ---------------------
C ------- 11:   H2I   + H   -> 3H
C ------- 12:   H2I   + HII -> H2II  + H
C ------- 13:   H2I   + e   -> HI    + HM
C ------- 14:   H2I   + e   -> 2HI   + e
C ------- 15:   H2I   + H2I -> H2I   + 2HI
C ------- 16:   HM    + e   -> HI    + 2e
C ------- 17:   HM    + HI  -> 2H    + e
C ------- 18:   HM    + HII -> 2HI
C ------- 19:   HM    + HII -> H2II  + e
C ------- 20:   H2II  + e   -> 2HI
C ------- 21:   H2II  + HM  -> HI    + H2I
C --------------  old  ---------------------
c
c --------------  new  ---------------------
C ---11--       H2I   + HII -> H2II  + H
C ---12--       H2I   + e   -> 2HI   + e
C ---13--       H2I   + H   -> 3H
C ---14--       HM    + e   -> HI    + 2e
C ---15--       HM    + HI  -> 2H    + e
C ---16--       HM    + HII -> 2HI
C ---17--       HM    + HII -> H2II  + e
C ---18--       H2II  + e   -> 2HI
C ---19--       H2II  + HM  -> HI    + H2I
c (20) - unused
c --------------  new  ---------------------
c
c
C ------- 21:   2H    + H2  -> H2I   + H2I
C ------- 22:   2H    + H   -> H2I   + H
C ------- 24:   HI    + p   -> HII   + e
C ------- 25:   HeII  + p   -> HeIII + e
C ------- 26:   HeI   + p   -> HeII  + e
C ------- 27:   HM    + p   -> HI    + e
C ------- 28:   H2II  + p   -> HI    + HII
C ------- 29:   H2I   + p   -> H2II  + e
C ------- 30:   H2II  + p   -> 2HII  + e
C ------- 31:   H2I   + p   -> 2HI
C
c ------- 50-56: deuterium rates (given below)
c
c ------- h2dust: 2H  + grain -> H2  + grain (eq 3.8 from Hollenbach & McKee 1989)
c
c ncrn, ncrd1, ncrd1: numerator and denominator terms in Eq 23 of Omukai ea 2000
c               for H2 formation heating.
c
c-----------------------------------------------------------------------
c
c  This define indicates if the Tom Abels new (as of Feb/2000) collisional
c    rates (k1-k19) should be used.
c
#define USE_NEW_COLLISIONAL_RATES
c
      implicit NONE
#include "fortran_types.def"
c
c  Arguments
c
      INTG_PREC nratec, iradtype, casebrates, threebody, ndratec
      R_PREC  aye, temstart, temend, alpha0, f3, dtemstart, dtemend,
     &        utem, uxyz, uaye, urho, utim
      R_PREC  ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      R_PREC  hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldl(nratec),
     &        gphdl(nratec), hdltea(nratec), hdlowa(nratec)
      R_PREC  gaHIa(nratec), gaH2a(nratec), gaHea(nratec), 
     &        gaHpa(nratec), gaela(nratec), gasgra(nratec),
     &        ciecoa(nratec)
      R_PREC  compa, gammaha, gammahacgs, piHI, piHeI, piHeII
      R_PREC  k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &        k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &        k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &        k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &        k17a(nratec), k18a(nratec), k19a(nratec), k20a(nratec), 
     &        k21a(nratec), k23a(nratec), k24, k25, k26, k27, k28, k29,
     &        k30, k31,
     &        k22a(nratec), k50a(nratec), k51a(nratec), k52a(nratec), 
     &        k53a(nratec), k54a(nratec), k55a(nratec), k56a(nratec)
      INTG_PREC ioutput
      R_PREC  k13dda(nratec, 7), hdcool(nratec, 5),
     &        h2dusta(nratec, ndratec),
     &        ncrna(nratec), ncrd1a(nratec), ncrd2a(nratec)
c
c  Parameters
c
      INTG_PREC nnu
      parameter (nnu = 400)
      INTG_PREC ireco, ibrco, icico, iceco, iphoto_mol
      parameter (ireco = 1, ibrco = 1, icico = 1, iceco = 1, 
     &           iphoto_mol = 0)

      real*8 everg, evhz, tevk, mh, pi, dhuge, kb
      parameter (everg = ev2erg, evhz  = 2.41838d14, 
     &           tevk = 1.1605d4, mh = mass_h, 
     &           pi=pi_val, dhuge = 1.0d30, 
     &           kb = kboltz)

c
c         Set various cutoff values in eV.
c
      real*8 e24,e25,e26,e27,e28a,e28b,e28c,e29a,e29b,e29c,
     &                 e30a,e30b
        PARAMETER(
     &    e24  = 13.6d0
     &   ,e25  = 54.4d0
     &   ,e26  = 24.6d0
     &   ,e27  = 0.755d0
     &   ,e28a = 2.65d0
     &   ,e28b = 11.27d0
     &   ,e28c = 21.d0
     &   ,e29a = 15.42d0
     &   ,e29b = 16.5d0
     &   ,e29c = 17.7d0
     &   ,e30a = 30.d0
     &   ,e30b = 70.d0)

c  Locals
c
      INTG_PREC i,j,n
      real*8 nu, delnu, hnu, logttt, ttt, tev, logtev, 
     &                 xx, dum, tbase1, xbase1, kunit, coolunit,
     &                 dbase1, dlogtem, kunit_3bdy, x, y, cierate,
     &                 grain_coef, d_ttt, d_logttt, d_dlogtem,
     &                 ttt2, d_ttt2
      real*8 sigma24(nnu), sigma25(nnu), sigma26(nnu), 
     &                 sigma27(nnu), sigma28(nnu), sigma29(nnu), 
     &                 sigma30(nnu), sigma31(nnu), fext(nnu)
      real*8 tm
      real*8 HDLR, HDLV, lt, t3, lt3
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c  Get conversion units
c
c     t/x/dbase1 is the number (z dependant) that converts from the
c       dimensionless code units to physical units.  Also, in the
c       code aye = 1 at z=zinit, so to convert the usual a (=1 at z=0)
c       to a~ (written in the code as aye), we use a = a~*[a] 
c
      tbase1 = utim
      xbase1 = uxyz/(aye*uaye)    ! uxyz is [x]*a     
      dbase1 = urho*(aye*uaye)**3 ! urho is [dens]/a^3
c
c  1) Set the dimensions of the (non-radiative) rate coefficients.  
c    Note that we have included the units that convert density to 
c    number density, so the rate equations should look like 
c    (in dimensionless units, hence the primes):
c
c       d(d0~)/dt~ = k~ * d1~ * d2~ / a~^3
c
c    where k~ is the dimenionless rate coefficients and d0-2~ are three
c     dimensionless densities (i.e. d = [dens]*d~) and a~ is the 
c     dimensionless expansion coefficient (see above).
c
c    rate eqn        : delta(n0)  = k  * n1        * n2        * dt     / a^3
c    rate eqn units  : [dens]/mh  = k  * [dens]/mh * [dens]/mh * [time] / [a]^3
c    rate eqn dimless: delta(n0~) = k~ * n1~       * n2~       * dt~    / a~^3
c    so: k = [k] * k~  where [k] = ( [a]^3 * mh ) / ( [dens] * [time] )  (~)
c    reminder: the number densities here are normalized with [dens] which
c              is not a constant (it has a factor a^3), so the number
c              densities must be converted from comoving to proper.
c
      kunit   = (uaye**3 * mh) / (dbase1 * tbase1)
      kunit_3bdy  = kunit * (uaye**3 * mh) / dbase1
c
c  2) Set the dimension of the cooling coefficients (including constants)
c     (this equation has a rho because e is the specific energy, not
c      energy/unit volume).
c       delta(e)  = L     * n1        * n2        * dt     / dens   / a^3
c       [e]       = L     * [dens]/mh * [dens]/mh * [time] / [dens] / [a]^3
c       delta(e~) = L~    * n1~       * n2~       * dt~    / dens~  / a~^3 [~]
c     so L = [L] * L~ where [L] = [e] * mh**2 * [a]^3 / ([dens] * [time]) [~]
c       but [e] = ([a]*[x])**2 / [time]**2 and ([a] = 1 / (1 + zri) )
c      [L] = ([a]**5 * [x]**2 * mh**2) / ([dens] * [time]**3)
c
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
c
c    Note: some of the coffiecients have only one power of n.  These
c          do not have the /a^3 factor, also they have units
c          [L1] = ([a]**2 * [x]**2 * mh) / [time]**3
c               = [L] * [dens] * [a]**3 / mh
c          This is done through the dom variable in cool.src
c         (some have three powers of n and they are different by the
c          reciprocal of the above factor multiplying [L]).
c
c  3) the units for the radiative rate coefficients is just 1/[time]
c
c  4) Energy transfer from gas to dust grains, following equation 2.15
c     of Hollenbach & McKee (1989).
      grain_coef = 1.2d-31 * 1.0d3**(-0.5d0)
c
c  Compute log spacing in temperature
c
      ttt    = temstart
      logttt = log(ttt)
      dlogtem= (log(temend) - log(temstart))/REAL(nratec-1,8)
c
c  Compute log spacing in dust temperature
c
      d_ttt     = dtemstart
      d_logttt  = log(d_ttt)
      d_dlogtem = (log(dtemend) - log(dtemstart))/real(ndratec-1)
c
c  Initialize constants to tiny
c
      do i = 1, nratec
c
        k1a (i) = tiny
        k2a (i) = tiny
        k3a (i) = tiny
        k4a (i) = tiny
        k5a (i) = tiny
        k6a (i) = tiny
        k7a (i) = tiny
        k8a (i) = tiny
        k9a (i) = tiny
        k10a(i) = tiny
        k11a(i) = tiny
        k12a(i) = tiny
        k13a(i) = tiny
        k14a(i) = tiny
        k15a(i) = tiny
        k16a(i) = tiny
        k17a(i) = tiny
        k18a(i) = tiny
        k19a(i) = tiny
        k20a(i) = tiny
        k21a(i) = tiny
        k22a(i) = tiny
        k50a(i) = tiny
        k51a(i) = tiny
        k52a(i) = tiny
        k53a(i) = tiny
        k54a(i) = tiny
        k55a(i) = tiny
        k56a(i) = tiny
c
        do j = 1, 7
           k13dda(i,j) = tiny
        enddo
c
        do j = 1, ndratec
           h2dusta(i,j) = tiny
        enddo
c
        ncrna(i)  = tiny
        ncrd1a(i) = tiny
        ncrd2a(i) = tiny
c
        ceHIa(i)    = tiny
        ceHeIa(i)   = tiny
        ceHeIIa(i)  = tiny
        ciHIa(i)    = tiny
        ciHeIa(i)   = tiny
        ciHeISa(i)  = tiny
        ciHeIIa(i)  = tiny
        reHIIa(i)   = tiny
        reHeII1a(i) = tiny
        reHeII2a(i) = tiny
        reHeIIIa(i) = tiny
        brema(i)    = tiny
        compa       = tiny
c
        hyd01ka(i)  = tiny
        h2k01a(i)   = tiny
        vibha(i)    = tiny
        rotha(i)    = tiny
        rotla(i)    = tiny
        hdltea(i)   = tiny
        hdlowa(i)   = tiny
        ciecoa(i)   = tiny
c
        gaHIa(i)    = tiny
        gaH2a(i)    = tiny
        gaHea(i)    = tiny
        gaHpa(i)    = tiny
        gaela(i)    = tiny
        gasgra(i)   = tiny
      enddo
c
c  Fill in tables over the range temstart to temend
c
c -------------------------------------------------
c  1) rate coefficients (excluding external radiation field)
c
      do i = 1, nratec
c
c       Compute temperature of this bin (in eV)
c
        logttt = log(temstart) + REAL(i-1,8)*dlogtem
        ttt = exp(logttt)
        ttt2 = ttt / 1.0d2
        tev = ttt/tevk
        logtev = log(tev)
c
#ifdef USE_NEW_COLLISIONAL_RATES
c
c       Call Tom Abels routine (from his web page, Feb/2000)
c
        call coll_rates(ttt, k1a(i), k2a(i), k3a(i), k4a(i), k5a(i),
     &                  k6a(i), k7a(i), k8a(i), k9a(i), k10a(i),
     &                  k11a(i), k12a(i), k13a(i), k14a(i), k15a(i), 
     &                  k16a(i), k17a(i), k18a(i), k19a(i), k23a(i),
     &                  kunit, casebrates)
c
#else /* USE_NEW_COLLISIONAL_RATES */
c
c       Compute rates in cgs: cm^3/s (note: tev, logtev is double)
c        (the k1-19 rates below are an older set)
c 
        if (tev .gt. 0.8d0) then
          k1a(i) = exp(-32.71396786375d0
     &           + 13.53655609057d0*logtev
     &           - 5.739328757388d0*logtev**2 
     &           + 1.563154982022d0*logtev**3
     &           - 0.2877056004391d0*logtev**4
     &           + 0.03482559773736999d0*logtev**5
     &           - 0.00263197617559d0*logtev**6
     &           + 0.0001119543953861d0*logtev**7
     &           - 2.039149852002d-6*logtev**8) / kunit
c
          k3a(i) = exp(-44.09864886561001d0
     &           + 23.91596563469d0*logtev
     &           - 10.75323019821d0*logtev**2
     &           + 3.058038757198d0*logtev**3
     &           - 0.5685118909884001d0*logtev**4
     &           + 0.06795391233790001d0*logtev**5
     &           - 0.005009056101857001d0*logtev**6
     &           + 0.0002067236157507d0*logtev**7
     &           - 3.649161410833d-6*logtev**8) / kunit
c
          k4a(i) = (1.54d-9*(1.d0+0.3d0 / 
     &         exp(8.099328789667d0/tev))
     &           / (exp(40.49664394833662d0/tev)*tev**1.5d0)
     &           + 3.92d-13/tev**0.6353d0) / kunit
c
          k5a(i) = exp(-68.71040990212001d0
     &           + 43.93347632635d0*logtev
     &           - 18.48066993568d0*logtev**2
     &           + 4.701626486759002d0*logtev**3
     &           - 0.7692466334492d0*logtev**4
     &           + 0.08113042097303d0*logtev**5
     &           - 0.005324020628287001d0*logtev**6
     &           + 0.0001975705312221d0*logtev**7
     &           - 3.165581065665d-6*logtev**8) / kunit
        else
          k1a(i) = tiny
          k3a(i) = tiny
          k4a(i) = 3.92d-13/tev**0.6353d0 / kunit
          k5a(i) = tiny
        endif
c
c       redefine k4a if case B recombination rates are requested
        if (casebrates) then
           k4a(i) = 1.26d-14 * (5.7067d5/ttt)**(0.75d0)
        endif
c
c       compute the HII rec. rate using either the case A or case B formula
        if (casebrates) then
           if (ttt < 1.0d9) then
              k2a(i) = 4.881357d-6*ttt**(-1.5d0) 
     &             * (1.d0+1.14813d2
     &             * ttt**(-0.407d0))**(-2.242d0) / kunit
           else
              k2a(i) = tiny
           endif
        else
           if ( ttt .gt. 5500.d0 ) then
              k2a(i) = exp(-28.61303380689232d0
     &             - 0.7241125657826851d0*logtev
     &             - 0.02026044731984691d0*logtev**2
     &             - 0.002380861877349834d0*logtev**3
     &             - 0.0003212605213188796d0*logtev**4
     &             - 0.00001421502914054107d0*logtev**5
     &             + 4.989108920299513d-6*logtev**6
     &             + 5.755614137575758d-7*logtev**7
     &             - 1.856767039775261d-8*logtev**8
     &             - 3.071135243196595d-9*logtev**9) / kunit
           else
              k2a(i) = k4a(i)
           endif
        endif
c
c       compute the HeIII rec. rate using either the case A or case B formula
        if (casebrates) then
           if (ttt < 1.0d9) then
              k6a(i) = 7.8155e-5d0*ttt**(-1.5d0) 
     &             * (1.d0+2.0189d2 
     &             * ttt**(-0.407d0))**(-2.242d0) / kunit
           else
              k6a(i) = tiny
           endif
        else
           k6a(i) = 3.36d-10/sqrt(ttt)/(ttt/1.d3)**0.2d0
     &          / (1.d0+(ttt/1.d6)**0.7d0) / kunit
        endif
c
c       Molecular hydrogen and associated rates
c
        k7a(i) = 6.77d-15*tev**0.8779d0 / kunit
c
        if (tev .gt. 0.1d0) then
          k8a(i) = exp(-20.06913897587003d0
     &           + 0.2289800603272916d0*logtev
     &           + 0.03599837721023835d0*logtev**2
     &           - 0.004555120027032095d0*logtev**3
     &           - 0.0003105115447124016d0*logtev**4
     &           + 0.0001073294010367247d0*logtev**5
     &           - 8.36671960467864d-6*logtev**6
     &           + 2.238306228891639d-7*logtev**7) / kunit
        else
          k8a(i) = 1.43d-9 / kunit
        endif
c
        k9a(i) = 1.85d-23*ttt**1.8d0 / kunit
        if(ttt .gt. 6.7d3) 
     &     k9a(i) = 5.81d-16*(ttt/56200.d0)**(-0.6657d0*
     &       log10(ttt/56200.d0)) / kunit
c
        k10a(i) = 6.0d-10 / kunit
c
        if (tev .gt. 0.3d0) then
          k13a(i) = 1.0670825d-10*tev**2.012d0
     &          / (exp(4.463d0/tev)
     &          * (1.d0+0.2472d0*tev)**3.512d0) / kunit
c
          k11a(i) = exp(-24.24914687731536d0
     &            + 3.400824447095291d0*logtev
     &            - 3.898003964650152d0*logtev**2
     &            + 2.045587822403071d0*logtev**3
     &            - 0.5416182856220388d0*logtev**4
     &            + 0.0841077503763412d0*logtev**5
     &            - 0.007879026154483455d0*logtev**6
     &            + 0.0004138398421504563d0*logtev**7
     &            - 9.36345888928611d-6*logtev**8) / kunit

          k12a(i) = 4.38d-10*exp(-102000.d0/ttt)
     &         * ttt**0.35d0 / kunit
        else
          k13a(i) = tiny
          k11a(i) = tiny
          k12a(i) = tiny
        endif
c
        if (tev .gt. 0.04d0) then
          k14a(i) = exp(-18.01849334273d0
     &            + 2.360852208681d0*logtev
     &            - 0.2827443061704d0*logtev**2
     &            + 0.01623316639567d0*logtev**3
     &            - 0.03365012031362999d0*logtev**4
     &            + 0.01178329782711d0*logtev**5
     &            - 0.001656194699504d0*logtev**6
     &            + 0.0001068275202678d0*logtev**7
     &            - 2.631285809207d-6*logtev**8) / kunit
        else
          k14a(i) =  tiny
        endif
c
        if (tev .gt. 0.1d0) then
          k15a(i) = exp(-20.37260896533324d0
     &            + 1.139449335841631d0*logtev
     &            - 0.1421013521554148d0*logtev**2
     &            + 0.00846445538663d0*logtev**3
     &            - 0.0014327641212992d0*logtev**4
     &            + 0.0002012250284791d0*logtev**5
     &            + 0.0000866396324309d0*logtev**6
     &            - 0.00002585009680264d0*logtev**7
     &            + 2.4555011970392d-6*logtev**8
     &            - 8.06838246118d-8*logtev**9) / kunit
        else
          k15a(i) = 2.56d-9*tev**1.78186d0 / kunit
        endif
c
        k16a(i) = 6.5d-9/sqrt(tev) / kunit
c
        k17a(i) = 1.0d-8*ttt**(-0.4d0) / kunit
        if(ttt.gt.1.0d4)
     &     k17a(i)=4.0d-4*ttt**(-1.4d0)*exp(-15100.d0/ttt) / kunit
c
        k18a(i) = 5.56396d-8/tev**0.6035d0 / kunit
        k19a(i) = 4.64d-8/sqrt(tev) / kunit
c
#endif /* USE_NEW_COLLISIONAL_RATES */
c
c       Compute the density-dependant collision H2 dissociation rates.
c          (this givens a 7-variable function that is used to compute
c           the log10 of the rate -- normalize by dividing by kunit).
c
        call colh2diss(ttt, k13dda(i,1), k13dda(i,2), k13dda(i,3), 
     &                 k13dda(i,4), k13dda(i,5), k13dda(i,6), 
     &                 k13dda(i,7))
        k13dda(i,1) = k13dda(i,1) - log10(kunit)
c
c       ------ 3-body H2 rate ----
c       The first bit is my fit to A.E. Orel 1987, J.Chem.Phys., 87, 
c       314. I then match it to the T^-1 of Palla etal (1983) 
c       Which is then 4 times smaller than the Palla rate ! Thats
c       molecule uncertainties ! :-)
c
c Our varying threebody rates from Simon
        if(threebody.eq.0)then
          if (ttt .le. 300.d0) then
             k22a(i) = 1.3d-32 * (ttt/300.d0)**(-0.38d0)
          else
             k22a(i) = 1.3d-32 * (ttt/300.d0)**(-1.d0)
          endif
        elseif(threebody.eq.1)then
c Inverse of PSS83 three-body rate
          k22a(i) = 5.5d-29 / ttt
          k13a(i) = (5.24d-7 / ttt**0.485d0) 
     &                   * exp(-5.2d4 / ttt)
        elseif(threebody.eq.2)then
c Inverse of CW83 rate
          k22a(i) = 8.8d-33
          k13a(i) = 8.4d-11 * ttt**0.515d0
     &         * exp(-5.2d4 / ttt)
        elseif(threebody.eq.3)then
c Rate from JGC67, used by FH07 to construct their three-body rate
          k22a(i) = 1.44d-26 / ttt**1.54d0
          k13a(i) = (1.38d-4 / ttt**1.025d0)
     &         * exp(-5.2d4 / ttt)
        elseif(threebody.eq.4)then
c High density rate from MSM96
          k22a(i) = 7.7d-31 / ttt**0.464d0
          k13a(i) = 10.d0**(-178.4239d0 - 68.42243d0 
     &         * log10(ttt) + 43.20243d0 * log10(ttt)**2
     &         - 4.633167d0 * log10(ttt)**3 
     &         + 69.70086d0 * log10(1d0 + 40870.38d0 / ttt)
     &         - (23705.7d0 / ttt))
        else
            write(0,*) "I didn't understand your threebody rate."
            stop
        endif
        k13a(i) = k13a(i) / kunit
        k22a(i) = k22a(i) / kunit_3bdy
        
        k21a(i) = 2.8d-31 * (ttt**(-0.6d0)) / kunit_3bdy
c        write(6,*),'calc_rates: check: ', ttt, k2a(i), k22a(i)
c
c       ------ Deuterium rates -----
c
c       50) H+  +  D  -> H  +  D+
c       51) H   +  D+ -> H+ +  D
c       52) H2  +  D+ -> HD +  H+
c       53) HD  +  H+ -> H2 +  D+
c       54) H2  +  D  -> HD +  H
c       55) HD  +  H  -> H2 +  D
c       56) D   +  H- -> HD +  e-
c      [57) D-  +  H  -> HD +  e-]  included by multiply 56 by 2
c
        k50a(i) = 1.0d-9 *exp(-4.1d1/ttt)  / kunit
        k51a(i) = 1.0d-9                   / kunit
        k52a(i) = 2.1d-9                   / kunit
        k53a(i) = 1.0d-9 *exp(-4.57d2/ttt) / kunit
        k54a(i) = 7.5d-11*exp(-3.82d3/ttt) / kunit
        k55a(i) = 7.5d-11*exp(-4.24d3/ttt) / kunit
        k56a(i) = 1.5d-9 *(ttt/300.d0)**(-0.1d0) / kunit
c
c     H2 formation on dust grains.
c     Loop over dust temperature.
c
        do j = 1, ndratec
c
c     Compute dust temperature of this bin
c
           d_logttt = log(dtemstart) + real(j-1)*d_dlogtem
           d_ttt = exp(d_logttt)
           d_ttt2 = d_ttt / 1.0d2
c
c     k23 from Omukai (2000).
c
           h2dusta(i,j) = 6.0d-17 * (ttt / 300.0d0)**(0.5d0) *
     &          ((1.0d0 + exp(7.5d2 * ((1.0d0 / 75.0d0) - 
     &                                (1.0d0 / d_ttt))))**(-1.0d0)) *
     &          ((1.0d0 + (4.0d-2 * (ttt + d_ttt)**(0.5d0)) + 
     &           (2.0d-3 * ttt) + 
     &           (8.0d-6 * ttt**(2.0d0)))**(-1.0d0)) / kunit
c
c     Equation 3.8 from Hollenbach & McKee (1979).
c
c           h2dusta(i,j) = 3.0d-17 * ttt2**0.5d0 /
c     &          (1.0d0 + 
c     &           (0.4d0 * (ttt2 + d_ttt2)**0.5d0) +
c     &           (0.2d0 * ttt2) +
c     &           (8.0d-2 * ttt2**2.0d0)) / kunit
c
        enddo
c
c       H2 formation heating terms from Equation 23 of Omukai (2000).
c
        ncrna(i)  = 1.0d6 * (ttt**(-0.5d0))
        ncrd1a(i) = 1.6d0 * exp(-(400.0d0 / ttt)**(2.d0))
        ncrd2a(i) = 1.4d0 * exp(-12000.0d0 / (ttt + 1200.0d0))
c
      enddo
c
c  Write out cgs rate coefficients
c
c      open(10, file='k1-6.out', status='unknown')
c      do i = 1, nratec
c         write(10,1010) exp(log(temstart) + REAL(i-1,RKIND)*dlogtem),
c     &               k1a(i),k2a(i),k3a(i),k4a(i),k5a(i),k6a(i)
c      enddo
c 1010 format(7(1pe11.4))
c      close(10)
c
c -------------------------------------------------
c  2) Cooling/heating rates (excluding external radiation field)
c
      do i = 1, nratec
c
        logttt = log(temstart) + REAL(i-1,RKIND)*dlogtem
        ttt = exp(logttt)
c
c    a) Collisional excitations (Black 1981; Cen 1992)
c
      if ( iceco .eq. 1 ) then
        ceHIa(i) = 7.5d-19*exp(-min(log(dhuge),
     &        118348.d0/ttt))/(1.d0+sqrt(ttt/1.d5)) 
     &        / coolunit
        ceHeIa(i) = 9.1d-27*exp(-min(log(dhuge),13179.d0/ttt))
     &       *ttt**(-0.1687d0)/(1.d0+sqrt(ttt/1.d5)) 
     &       / coolunit
        ceHeIIa(i) = 5.54d-17*exp(-min(log(dhuge),
     &       473638.d0/ttt))*ttt**(-0.397d0) 
     &       / (1.d0+sqrt(ttt/1.0d5)) / coolunit
      else
        ceHIa(i)    = tiny
        ceHeIa(i)   = tiny
        ceHeIIa(i)  = tiny
      endif
c
c    b) Collisional ionizations (Cen 1992 or Abel 1996)
c
      if ( icico .eq. 1 ) then
c
c        ciHIa(i)    = 1.27d-21*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-min(log(dhuge),157809.1/ttt)) / coolunit
c        ciHeIa(i)   = 9.38d-22*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-min(log(dhuge),285335.4/ttt)) / coolunit
c        ciHeIIa(i)  = 4.95d-22*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-min(log(dhuge),631515.0/ttt)) / coolunit
c
        ciHeISa(i)  = 5.01d-27*(ttt)**(-0.1687d0) 
     &        / (1.d0+sqrt(ttt/1.d5))
     &        *exp(-min(log(dhuge),55338.d0/ttt)) / coolunit
c
c       Collisional ionizations (polynomial fits from Tom Abel)
c
        ciHIa(i)    = 2.18d-11*k1a(i)*kunit / coolunit
        ciHeIa(i)   = 3.94d-11*k3a(i)*kunit / coolunit
        ciHeIIa(i)  = 8.72d-11*k5a(i)*kunit / coolunit
      else
        ciHIa(i)    = tiny
        ciHeIa(i)   = tiny
        ciHeIIa(i)  = tiny
        ciHeISa(i)  = tiny
      endif
c
c    c) Recombinations (Cen 1992)
c
      if ( ireco .eq. 1 ) then
        reHIIa(i)   = 8.70d-27*sqrt(ttt)
     &        * (ttt/1000.d0)**(-0.2d0)
     &        / (1.d0 + (ttt/1.0d6)**(0.7d0)) / coolunit
        reHeII1a(i) = 1.55d-26*ttt**0.3647d0 / coolunit
        reHeII2a(i) = 1.24d-13*ttt**(-1.5d0) *
     &       exp(-min(log(dhuge),470000.d0/ttt)) *
     &       (1.d0+0.3d0*exp(-min(log(dhuge),94000.d0/ttt))) 
     &       / coolunit
        reHeIIIa(i) = 3.48d-26*sqrt(ttt)
     &       * (ttt/1000.d0)**(-0.2d0)
     &       / (1.d0 + (ttt/1.d6)**(0.7d0)) / coolunit
      else
        reHIIa(i)   = tiny
        reHeII1a(i) = tiny
        reHeII2a(i) = tiny
        reHeIIIa(i) = tiny
      endif
c
c    d) Bremsstrahlung (Black 1981)(Spitzer & Hart 1979)
c
      if ( ibrco .eq. 1 ) then
        brema(i)    = 1.43d-27*sqrt(ttt)*(1.1d0+0.34d0
     &        *exp(-(5.5d0-log10(ttt))**2/3.d0)) / coolunit
      else
        brema(i)    = tiny
      endif
c
c         Bremsstrahlung (Shapiro & Kang 1987)(Spitzer 1978)
c
c        if(ttt .lt. 3.2d5) gaunt = 0.79464d0 + 0.1243d0*log10(ttt)
c        if(ttt .ge. 3.2d5) gaunt = 2.13164 - 0.1240d0*log10(ttt)
c        brem1a(i) = 1.426d-27*sqrt(ttt)*gaunt / coolunit
c        if(ttt .lt. 4.0*3.2d5) gaunt = 0.79464d0 + 0.1243d0*log10(ttt/4.d0)
c        if(ttt .ge. 4.0*3.2d5) gaunt = 2.13164d0 - 0.1240d0*log10(ttt/4.d0)
c        brem2a(i) = 1.426d-27*sqrt(ttt)*gaunt / coolunit
c
c    e) Molecular hydrogen cooling
c
c       (Which one is used is controlled by a flag in cool1d_mulit.src)
c
c    e part1) - Lepp & Shull rates
c
        xx = log10(ttt/1.d4)
        vibha  (i) = 1.1d-18*exp(-min(log(dhuge),6744.d0/ttt)) 
     &               / coolunit
c
        if( ttt .gt. 1635.d0) then
           dum = 1.0d-12*sqrt(ttt)*exp(-1000.d0/ttt)
        else
           dum = 1.4d-13*exp((ttt/125.d0)-(ttt/577.d0)**2)
        endif
        hyd01ka(i) = dum*exp(-min(log(dhuge),
     &                            8.152d-13/(kboltz*ttt) )) 
     &               / coolunit
c
        dum = 8.152d-13*(4.2d0 /
     &       (kboltz*(ttt+1190.d0)) +
     &       1.d0/(kboltz*ttt))
        h2k01a(i) = 1.45d-12*sqrt(ttt)*exp(-min(log(dhuge),dum)) 
     &              / coolunit
c
        if(ttt .gt. 4031.d0)then
          rotla(i) = 1.38d-22*exp(-9243.d0/ttt) / coolunit
        else
          rotla(i) = 10.d0**(-22.9d0 - 0.553d0*xx 
     &          - 1.148d0*xx*xx) / coolunit
        endif
c
        if(ttt .gt. 1087.d0)then
           rotha(i) = 3.9d-19*exp(-6118.d0/ttt) / coolunit
        else
           rotha(i) = 10.d0**(-19.24d0 + 0.474d0*xx
     &          - 1.247d0*xx*xx) / coolunit
        endif
c
c     e part2) - Galli and Palla (1999) rates as fit by Tom Abel
c
        tm = max(ttt, 13.d0)       ! no cooling below 13 Kelvin ...
        tm = min(tm, 1.d5)         ! fixes numerics
        lt = log10(tm)
c     low density limit from Galli and Palla
        gpldl(i)  = 10.d0**(-103.d0+97.59d0*lt
     &       - 48.05d0*lt**2+10.8d0*lt**3-0.9032d0*lt**4) 
     &       / coolunit
c     high density limit from HM79 (typo corrected Aug 30/2007)
        t3 = tm/1000.d0
        HDLR = ((9.5d-22*t3**3.76d0)
     &       / (1.d0+0.12d0*t3**2.1d0)*
     &       exp(-(0.13d0/t3)**3)+3.d-24*exp(-0.51d0/t3))
        HDLV = (6.7d-19*exp(-5.86d0/t3) + 1.6d-18 
     &       * exp(-11.7d0/t3))
        gphdl(i)  = (HDLR + HDLV) / coolunit
c
c     e part 3) - Low density rates from Glover & Abel (2008)
c
c  Excitation by HI
c
      tm  = max(ttt, 10.d0)
      tm  = min(tm, 1.d4)
      lt3 = log10(tm / 1.d3)  

      if (tm .lt. 1.d2) then
        gaHIa(i) = 10.d0**(-16.818342d0
     &         + 37.383713d0 * lt3
     &         + 58.145166d0 * lt3**2
     &         + 48.656103d0 * lt3**3
     &         + 20.159831d0 * lt3**4
     &         + 3.847961d0 * lt3**5) / coolunit
      elseif (tm .lt. 1.d3) then
        gaHIa(i) = 10.d0**(-24.311209d0
     &         + 3.5692468d0 * lt3
     &         - 11.33286d0 * lt3**2
     &         - 27.850082d0 * lt3**3
     &         - 21.328264d0 * lt3**4
     &         - 4.2519023d0 * lt3**5) / coolunit
      else
        gaHIa(i) = 10.d0**(-24.311209d0
     &         + 4.6450521d0 * lt3
     &         - 3.7209846d0 * lt3**2
     &         + 5.9369081d0 * lt3**3
     &         - 5.5108047d0 * lt3**4
     &         + 1.5538288d0 * lt3**5) / coolunit
      endif
c
c Excitation by H2
c
      gaH2a(i) = 10.d0**(-23.962112d0
     &        + 2.0943374d0  * lt3
     &        - 0.77151436d0  * lt3**2
     &        + 0.43693353d0  * lt3**3
     &        - 0.14913216d0  * lt3**4
     &        - 0.033638326d0 * lt3**5) / coolunit
c
c Excitation by He
c
      gaHea(i) = 10.d0**(-23.689237d0
     &        + 2.1892372d0  * lt3
     &        - 0.81520438d0 * lt3**2
     &        + 0.29036281d0 * lt3**3
     &        - 0.16596184d0 * lt3**4
     &        + 0.19191375d0 * lt3**5) / coolunit
c
c Excitation by H+
c
      gaHpa(i) = 10.d0**(-21.716699d0
     &        + 1.3865783d0   * lt3
     &        - 0.37915285d0  * lt3**2
     &        + 0.11453688d0  * lt3**3
     &        - 0.23214154d0  * lt3**4
     &        + 0.058538864d0 * lt3**5) / coolunit
c
c Excitation by electrons
c
      if (tm .lt. 200.d0) then
        gaela(i) = 10.d0**(-34.286155d0
     &          - 48.537163d0  * lt3
     &          - 77.121176d0  * lt3**2
     &          - 51.352459d0  * lt3**3
     &          - 15.16916d0  * lt3**4
     &          - 0.98120322d0 * lt3**5) / coolunit
      else
        gaela(i) = 10.d0**(-22.190316
     &          + 1.5728955d0  * lt3
     &          - 0.213351d0 * lt3**2
     &          + 0.96149759d0 * lt3**3
     &          - 0.91023195d0 * lt3**4
     &          + 0.13749749d0 * lt3**5) / coolunit
      endif
c
c     f) HD cooling
c
c        HD Cooling Function (ergs cm3 /s)
c
c Fit to Lepp and Shull 1984 LTE (ergs/s) -> hdlte (ergs cm3/s)
c
      hdltea(i) = -35.6998d0 + 15.35716d0*dlog10(ttt) -
     &            5.58513d0   * (dlog10(ttt))**2 +
     &            0.8561149d0 * (dlog10(ttt))**3 - 
     &            1.75538d-2  * (dlog10(ttt))**4
      hdltea(i) = (10.d0**min(hdltea(i), 0._RKIND)) / coolunit
c      hdlte=10**lhdlte/den
c
c Galli and Palla 1998 low density limit (erg cm3 /s)
c  uses HD-He collisional data, so reduced by 1.27
c
      hdlowa(i) = ((3.d0 * (4.4d-12 
     &     + 3.6d-13*ttt**0.77d0) *
     &     dexp(-128.d0/ttt) * 128.d0 +
     &     (5.d0/3.d0) * (4.1d-12 +
     &     2.1d-13*ttt**0.92d0) *
     &     dexp(-255.d0/ttt) * 255.d0) * kb/1.27d0 ) 
     &     / coolunit
c      hdcool=hdlte/(1.d0+hdlte/hdlow)
c
c     g) CIE cooling as discussed in Ripamonti & Abel (2003)
c        cie_thin_cooling_rate returns rate as erg s^-1 cm^3 g^2 so
c         divide my mh^2 to get usual erg s^-1 cm^3, which is then
c         converted to code units with the conversion factor coolunit.
c         To get cooling rate
c     
      call cie_thin_cooling_rate(ttt, cierate)
c      ciecoa(i) = cierate/mh/mh/coolunit
c
c     Matt Turk: cie_thin_cooling actually returns
c     erg per second times cm^3 per gram per (gram in H2 molecules)
c     To reproduce eq 5 in RA04, divide c_t_c_r by mh, so to get
c     erg/s cm^3 we multiply by mh.  We can either divide by 2 for 
c     the H2 mass->number or in cool1d_multi_bdf.src, but we might
c     as well do it here.
c
      ciecoa(i) = cierate * (mh/2.d0) / coolunit
c      write(6,*), "thin: ", ttt, cierate * 1e14*mh*2
c
c     h) Energy transfer from gas to dust grains, following 
c        equation 2.15 of Hollenbach & McKee (1989)
      gasgra(i) = grain_coef * (ttt**0.5d0) * 
     &     (1.0d0 - (0.8d0 * exp(-75.0d0 / ttt))) / coolunit
c
      enddo
c
c     i) Compton cooling (Peebles 1971)
c
      compa     = 5.65d-36 / coolunit
c
c     j) Photoelectric heating by UV-irradiated dust (Wolfire 1995)
c        Default is 8.5e-26 for epsilon=0.05, G_0=1.7 (rate in erg s^-1 cm^-3)
c
      gammaha   = gammahacgs / coolunit
c
c -------------------------------------------------
c  3) External radiative processes
c
c     Initialize to tiny
c
      k24 = tiny
      k25 = tiny
      k26 = tiny
      k27 = tiny
      k28 = tiny
      k29 = tiny
      k30 = tiny
      k31 = tiny
c      
      piHI   = tiny
      piHeI  = tiny
      piHeII = tiny
c
c     Loop over all frequency bins, compute cross-sections
c        nu is in eV (range is 0.74 eV to 7.2 keV if nnu=400)
c
      do n = 1, nnu
        nu = 10.d0**((n-14)*0.01d0)              ! 0.74 ev -- 7.24d9 ev
c
c       24) HI photo-ionization cross-section
c
        if (nu .gt. e24) then
           dum = sqrt(nu/e24-1.d0)
           sigma24(n) = 6.3d-18 * (e24/nu)**4
     &               *  exp(4.d0-4.d0*atan(dum)/dum)
     &               / (1.d0-exp(-2.d0*pi/dum))
        else
           sigma24(n) = 0.d0
        endif
c
c       25) HeII photo-ionization cross-section
c
        if (nu .gt. e25) then
           dum = sqrt(nu/e25-1.d0)
           sigma25(n) = 1.58d-18 * (e25/nu)**4
     &               * exp(4.d0-4.d0*atan(dum)/dum)
     &               / (1.d0-exp(-2.d0*pi/dum))
        else
           sigma25(n) = 0.d0
        endif
c
c       26) HeI photo-ionization cross-section
c
        if (nu .gt. e26) then
c           sigma26(n) = 7.42d-18*(1.66*(nu/e26)**(-2.05)
c     &                          - 0.66*(nu/e26)**(-3.05))
c
c            Now From Verland, et al (1996)
c
           x = nu/13.61d0 - 0.4434d0
           y = sqrt(x**2 + 2.136d0**2)
           sigma26(n) = 9.492d-16*((x-1)**2 + 2.039d0**2) *
     &                  y**(0.5d0*3.188d0-5.5d0) *
     &                  (1.d0 + sqrt(y/1.469d0))**(-3.188d0)        
        else
           sigma26(n) = 0.d0
        endif
c
c       27) HM    + p   -> HI    + e
c
        if (nu .gt. e27) then
          sigma27(n) = 2.11d-16*(nu-e27)**1.5d0/nu**3
        else
          sigma27(n) = 0.d0
        endif
c
c       28) H2II  + p   -> HI    + HII
c
        if (nu .gt. e28a .and. nu .le. e28b) then
          sigma28(n) = 10.d0**(-40.97d0+6.03d0*nu
     &          -0.504d0*nu**2+1.387d-2*nu**3)
        elseif (nu .gt. e28b .and. nu .lt. e28c) then
           sigma28(n) = 10.d0**(-30.26d0+2.79d0*nu
     &          -0.184d0*nu**2+3.535d-3*nu**3)
        else
          sigma28(n) = 0.d0
        endif
c
c       29) H2I   + p   -> H2II  + e
c
        if (nu .gt. e29a .and. nu .le. e29b) then
          sigma29(n) = 6.2d-18*nu - 9.4d-17
        elseif (nu .gt. e29b .and. nu .le. e29c) then
          sigma29(n) = 1.4d-18*nu - 1.48d-17
        elseif (nu .gt. e29c) then
          sigma29(n) = 2.5d-14*nu**(-2.71d0)
        else
          sigma29(n) = 0.d0
        endif
c
c       30) H2II  + p   -> 2HII  + e
c
        if (nu .ge. e30a .and. nu .lt. e30b) then
           sigma30(n) = 10.d0**(-16.926d0-4.528d-2*nu
     &          +2.238d-4*nu**2+4.245d-7*nu**3)
        else
          sigma30(n) = 0.d0
        endif
c
c       31) H2I   + p   -> 2HI
c
c        sigma31(n) = 0.d0

C	 BUG FIX: (n) subscripts below were (i)
C	 reported by Patrick McDonald, modified in revision 2310

        if (nu .ge. e28b .and. nu .lt. e24) THEN
           sigma31(n) = 3.71d-18
        else
           sigma31(n) = 0.d0
        endif
c
c       Set the radiation field (only used in iradtype = 5 or 8)

        if (iradtype .eq. 5) then
c
c          5) A featureless power-low spectrum (quasar-like)
c         (f3 is the amplitude of flux at the HI ionization edge)
c
           fext(n) = f3*(nu/e24)**(-alpha0)
c
        elseif (iradtype .eq. 8) then
c
c          8) power law with absorbing neutral medium with
c              10^22 cm^-2 column density
c         (f3 is the amplitude of flux at the HI ionization edge,
c          the min is to prevent under-runs).
c
           fext(n) = f3*(nu/12.86d0)**(-1.d0)*exp(-1.d0*
     &       min(1.d22*(sigma24(n) + 0.08d0*sigma26(n)), 
     &           100.d0))
        else
c
c          otherwise) no radiation field
c
           fext(n) = 0.d0
        endif
      enddo
c
c     Integrate over the frequency spectrum
c

!     call open_mpi_error_file( 'F89', 89, 'unknown' )

      do n = 1, nnu
c          delnu = (10.d0**((n-14)*0.01d0) - 10.d0**((n-15)*0.01d0))*evhz
c          hnu   = 0.5d0*(10.d0**((n-14)*0.01d0) + 10.d0**((n-15)*0.01d0))*everg

          delnu = (10.d0**((n-13.5d0)*0.01d0) - 
     &        10.d0**((n-14.5d0)*0.01d0))*evhz
          hnu   = 10.d0**((n-14)*0.01d0)*everg
c
c         Add to radiative rate coefficients
c
c            k24-k26) HI, HeII, HeI photo-ionization
c
          k24 = k24 + fext(n)*sigma24(n) * delnu/hnu * tbase1
          k25 = k25 + fext(n)*sigma25(n) * delnu/hnu * tbase1
          k26 = k26 + fext(n)*sigma26(n) * delnu/hnu * tbase1
c
c            k27-k31) HM, H2II, H2I, H2II, H2I photo-dissociation
c
          if (iphoto_mol .eq. 1) then
             k27 = k27 + fext(n)*sigma27(n) * delnu/hnu * tbase1
             k28 = k28 + fext(n)*sigma28(n) * delnu/hnu * tbase1
             k29 = k29 + fext(n)*sigma29(n) * delnu/hnu * tbase1
             k30 = k30 + fext(n)*sigma30(n) * delnu/hnu * tbase1
             k31 = k31 + fext(n)*sigma31(n) * delnu/hnu * tbase1
          endif
c
c            Add to HI,HeI,HeII photoionization heating coefficients
c
          dum = max(hnu - e24*everg, 0.d0)
          piHI = piHI + fext(n)*sigma24(n)*delnu/hnu*dum/coolunit

!         write(89,1089) hnu, fext(n),sigma24(n),k24, dum, piHI

 1089     format(1p,10(e14.4,1x))
c
          dum = hnu - e26*everg
          piHeI = piHeI + fext(n)*sigma26(n)*delnu/hnu*dum/coolunit
c
          dum = hnu - e25*everg
          piHeII = piHeII + fext(n)*sigma25(n)*delnu/hnu*dum/coolunit
      enddo

!     call close_mpi_error_file( 89 )

c
c     The secondary electrons produced by steep spectral sources
c       are accounted for here using data from Shull & Steenberg.
c       Note that this is only good for ionization fractions < 1e-3 or so
c          (now down in the solvers to cover all ionization fractions)
c
#ifdef UNUSED
      if (iradtype .eq. 8) then
c
c        Add two terms which account for the amount of energy that
c          goes into secondary ionizations of HI and HeI (ignore HeII)
c          (We use piHI and piHeI which is the available energy and
c           then divide by the ionization energy of HI).  We also need
c           to convert from heating units (erg/s) to rate units (/s)
c          (The 0.08 accounts for the ratio n(HI)/n(HeI) where we are
c           assuming full neutrality)
c
         k24 = k24 + 0.39d0*(piHI + 0.08d0*piHeI)/(e24*everg)
     &         * coolunit * tbase1
         k26 = k26 + 0.055d0*(piHI/0.08d0 + piHeI)/(e26*everg)
     &         * coolunit * tbase1
c
c        Now correct the photo-heating rates to reflect that much of
c          the energy does not go into heating. (ignore HeII)
c
         piHI  = piHI *0.13d0
         piHeI = piHeI*0.13d0
c         
      endif
#endif /* UNUSED */

#define OUTPUT_COOLING_RATE_DATA
#ifdef OUTPUT_COOLING_RATE_DATA
c
c     Write out cooling rate data
c

      if (ioutput.eq.1) then
        open(10, file='cool_rates.out', status='unknown')

!      call open_mpi_error_file( 'cool_rates.out', 10, 'unknown' )

        do i=1, nratec
           logttt = log(temstart) + REAL(i-1,RKIND)*dlogtem
           ttt = exp(logttt)
           write(10,1000) ttt, ceHIa(i), ceHeIa(i), ceHeIIa(i),
     &                    ciHIa(i), ciHeIa(i), ciHeISa(i),
     &                    ciHeIIa(i), reHIIa(i), reHeII1a(i),
     &                    reHeII2a(i), reHeIIIa(i), brema(i),
     &                    compa       
 1000      format(1p,30(e10.3,1x))
        enddo
        write(10,*) 'piHI=', piHI*coolunit
        write(10,*) 'piHeI=', piHeI*coolunit
        write(10,*) 'piHeII=', piHeII*coolunit

!      call close_mpi_error_file( 10 )

        close(10)

c
c     Write out species rate data
c

        open(10, file='rates.out', status='unknown')

!      call open_mpi_error_file( 'rates.out', 10, 'unknown' )

        do i=1, nratec
           logttt = log(temstart) + REAL(i-1,RKIND)*dlogtem
           ttt = exp(logttt)
           write(10,1000) ttt, k1a(i), k2a(i), k3a(i), k4a(i), k5a(i),
     &                    k6a(i), k7a(i), k8a(i), k9a(i), k10a(i),
     &                    k11a(i), k12a(i), k13a(i), k14a(i), k15a(i), 
     &                    k16a(i), k17a(i), k18a(i), k19a(i), k22a(i)
        enddo
        write(10,*) 'k24=', k24/tbase1
        write(10,*) 'k25=', k25/tbase1
        write(10,*) 'k26=', k26/tbase1
        write(10,*) 'k27=', k27/tbase1
        write(10,*) 'k28=', k28/tbase1
        write(10,*) 'k29=', k29/tbase1
        write(10,*) 'k30=', k30/tbase1
        write(10,*) 'k31=', k31/tbase1

!      call close_mpi_error_file( 10 )

       close(10)

      endif

#endif /* OUTPUT_COOLING_RATE_DATA */

c
      return
      end

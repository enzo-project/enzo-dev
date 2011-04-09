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
     &             reHeII2a, reHeIIIa, brema, compa, gammaha,
     &             piHI, piHeI, piHeII,
     &             hyd01ka, h2k01a, vibha, rotha, rotla, 
     &             gpldl, gphdl, hdltea, hdlowa, hdcool, ciecoa,
     &             gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &             k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &             k11a, k12a, k13a, k13dda, k14a, k15a, k16a, k17a,
     &             k18a, k19a, k20a, k21a, k22a, k23a,
     &             k24, k25, k26, k27, k28, k29, k30, k31, 
     &             k50a, k51a, k52a, k53a, k54a, k55a, k56a, ioutput
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
c    tevk     = conversion constant [eV] = [K] / tevk
c    ireco    = Recombination cooling flag (1 = on)
c    ibrco    = Brehmmstrahlung cooling flag (1 = on)
c    icico    = Collisional ionization flag (1 = on)
c    iceco    = Collisional excitation flag (1 = on)
c    iphoto_mol = Molecular photo-dissociation flag (1 = on)
c
c  INPUTS
c    nratec   = number of temperature bins in rate equations
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
      INTG_PREC nratec, iradtype, casebrates, threebody
      R_PREC    aye, temstart, temend, alpha0, f3,
     &        utem, uxyz, uaye, urho, utim
      R_PREC    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      R_PREC    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldl(nratec),
     &        gphdl(nratec), hdltea(nratec), hdlowa(nratec)
      R_PREC    gaHIa(nratec), gaH2a(nratec), gaHea(nratec), 
     &        gaHpa(nratec), gaela(nratec),
     &        ciecoa(nratec)
      R_PREC    compa, gammaha, piHI, piHeI, piHeII
      R_PREC    k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &        k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &        k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &        k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &        k17a(nratec), k18a(nratec), k19a(nratec), k20a(nratec), 
     &        k21a(nratec), k23a(nratec), k24, k25, k26, k27, k28, k29,
     &        k30, k31,
     &        k22a(nratec), k50a(nratec), k51a(nratec), k52a(nratec), 
     &        k53a(nratec), k54a(nratec), k55a(nratec), k56a(nratec)
      INTG_PREC ioutput
      R_PREC    k13dda(nratec, 7), hdcool(nratec, 5)
c
c  Parameters
c
      INTG_PREC nnu
      parameter (nnu = 400)
      INTG_PREC ireco, ibrco, icico, iceco, iphoto_mol
      parameter (ireco = 1, ibrco = 1, icico = 1, iceco = 1, 
     &           iphoto_mol = 0)

      double precision everg, evhz, tevk, mh, pi, tiny_cgs, dhuge, kb
      parameter (everg = ev2erg, evhz  = 2.41838e14_RKIND, 
     &           tevk = 1.1605e4_RKIND, mh = mass_h, 
     &           pi=pi_val, tiny_cgs = 1.0e-37_RKIND, 
     &           dhuge = 1.0e30_RKIND, kb = kboltz)

c
c         Set various cutoff values in eV.
c
      double precision e24,e25,e26,e27,e28a,e28b,e28c,e29a,e29b,e29c,
     &                 e30a,e30b
        PARAMETER(
     &    e24  = 13.6_RKIND
     &   ,e25  = 54.4_RKIND
     &   ,e26  = 24.6_RKIND
     &   ,e27  = 0.755_RKIND
     &   ,e28a = 2.65_RKIND
     &   ,e28b = 11.27_RKIND
     &   ,e28c = 21._RKIND
     &   ,e29a = 15.42_RKIND
     &   ,e29b = 16.5_RKIND
     &   ,e29c = 17.7_RKIND
     &   ,e30a = 30._RKIND
     &   ,e30b = 70._RKIND)

c  Locals
c
      INTG_PREC i,j,n
      double precision nu, delnu, hnu, logttt, ttt, tev, logtev, 
     &                 xx, dum, tbase1, xbase1, kunit, coolunit,
     &                 dbase1, dlogtem, kunit_3bdy, x, y, cierate
      double precision sigma24(nnu), sigma25(nnu), sigma26(nnu), 
     &                 sigma27(nnu), sigma28(nnu), sigma29(nnu), 
     &                 sigma30(nnu), sigma31(nnu), fext(nnu)
      R_PREC tm
      double precision HDLR, HDLV, lt, t3, lt3
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
c    where k~ is the dimenionless rate coefficients and _RKIND-2~ are three
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
c  Compute log spacing in temperature
c
      ttt    = temstart
      logttt = log(ttt)
      dlogtem= (log(temend) - log(temstart))/REAL(nratec-1,RKIND)
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
        logttt = log(temstart) + REAL(i-1,RKIND)*dlogtem
        ttt = exp(logttt)
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
        if (tev .gt. 0.8_RKIND) then
          k1a(i) = exp(-32.71396786375_RKIND
     &           + 13.53655609057_RKIND*logtev
     &           - 5.739328757388_RKIND*logtev**2 
     &           + 1.563154982022_RKIND*logtev**3
     &           - 0.2877056004391_RKIND*logtev**4
     &           + 0.03482559773736999_RKIND*logtev**5
     &           - 0.00263197617559_RKIND*logtev**6
     &           + 0.0001119543953861_RKIND*logtev**7
     &           - 2.039149852002e-6_RKIND*logtev**8) / kunit
c
          k3a(i) = exp(-44.09864886561001_RKIND
     &           + 23.91596563469_RKIND*logtev
     &           - 10.75323019821_RKIND*logtev**2
     &           + 3.058038757198_RKIND*logtev**3
     &           - 0.5685118909884001_RKIND*logtev**4
     &           + 0.06795391233790001_RKIND*logtev**5
     &           - 0.005009056101857001_RKIND*logtev**6
     &           + 0.0002067236157507_RKIND*logtev**7
     &           - 3.649161410833e-6_RKIND*logtev**8) / kunit
c
          k4a(i) = (1.54e-9_RKIND*(1._RKIND+0.3_RKIND / 
     &         exp(8.099328789667_RKIND/tev))
     &           / (exp(40.49664394833662_RKIND/tev)*tev**1.5_RKIND)
     &           + 3.92e-13_RKIND/tev**0.6353_RKIND) / kunit
c
          k5a(i) = exp(-68.71040990212001_RKIND
     &           + 43.93347632635_RKIND*logtev
     &           - 18.48066993568_RKIND*logtev**2
     &           + 4.701626486759002_RKIND*logtev**3
     &           - 0.7692466334492_RKIND*logtev**4
     &           + 0.08113042097303_RKIND*logtev**5
     &           - 0.005324020628287001_RKIND*logtev**6
     &           + 0.0001975705312221_RKIND*logtev**7
     &           - 3.165581065665e-6_RKIND*logtev**8) / kunit
        else
          k1a(i) = tiny
          k3a(i) = tiny
          k4a(i) = 3.92e-13_RKIND/tev**0.6353_RKIND / kunit
          k5a(i) = tiny
        endif
c
c       redefine k4a if case B recombination rates are requested
        if (casebrates) then
           k4a(i) = 1.26e-14_RKIND * (5.7067e5_RKIND/ttt)**(0.75_RKIND)
        endif
c
c       compute the HII rec. rate using either the case A or case B formula
        if (casebrates) then
           if (ttt < 1.0d9) then
              k2a(i) = 4.881357e-6_RKIND*ttt**(-1.5_RKIND) 
     &             * (1._RKIND+1.14813e2_RKIND
     &             * ttt**(-0.407_RKIND))**(-2.242_RKIND) / kunit
           else
              k2a(i) = tiny
           endif
        else
           if ( ttt .gt. 5500._RKIND ) then
              k2a(i) = exp(-28.61303380689232_RKIND
     &             - 0.7241125657826851_RKIND*logtev
     &             - 0.02026044731984691_RKIND*logtev**2
     &             - 0.002380861877349834_RKIND*logtev**3
     &             - 0.0003212605213188796_RKIND*logtev**4
     &             - 0.00001421502914054107_RKIND*logtev**5
     &             + 4.989108920299513e-6_RKIND*logtev**6
     &             + 5.755614137575758e-7_RKIND*logtev**7
     &             - 1.856767039775261e-8_RKIND*logtev**8
     &             - 3.071135243196595e-9_RKIND*logtev**9) / kunit
           else
              k2a(i) = k4a(i)
           endif
        endif
c
c       compute the HeIII rec. rate using either the case A or case B formula
        if (casebrates) then
           if (ttt < 1.0d9) then
              k6a(i) = 7.8155e-5_RKIND*ttt**(-1.5_RKIND) 
     &             * (1._RKIND+2.0189e2_RKIND 
     &             * ttt**(-0.407_RKIND))**(-2.242_RKIND) / kunit
           else
              k6a(i) = tiny
           endif
        else
           k6a(i) = 3.36e-10_RKIND/sqrt(ttt)/(ttt/1.e3_RKIND)**0.2_RKIND
     &          / (1._RKIND+(ttt/1.e6_RKIND)**0.7_RKIND) / kunit
        endif
c
c       Molecular hydrogen and associated rates
c
        k7a(i) = 6.77e-15_RKIND*tev**0.8779_RKIND / kunit
c
        if (tev .gt. 0.1_RKIND) then
          k8a(i) = exp(-20.06913897587003_RKIND
     &           + 0.2289800603272916_RKIND*logtev
     &           + 0.03599837721023835_RKIND*logtev**2
     &           - 0.004555120027032095_RKIND*logtev**3
     &           - 0.0003105115447124016_RKIND*logtev**4
     &           + 0.0001073294010367247_RKIND*logtev**5
     &           - 8.36671960467864e-6_RKIND*logtev**6
     &           + 2.238306228891639e-7_RKIND*logtev**7) / kunit
        else
          k8a(i) = 1.43e-9_RKIND / kunit
        endif
c
        k9a(i) = 1.85e-23_RKIND*ttt**1.8_RKIND / kunit
        if(ttt .gt. 6.7e3_RKIND) 
     &     k9a(i) = 5.81e-16_RKIND*(ttt/56200._RKIND)**(-0.6657_RKIND*
     &       log10(ttt/56200._RKIND)) / kunit
c
        k10a(i) = 6.0e-10_RKIND / kunit
c
        if (tev .gt. 0.3_RKIND) then
          k13a(i) = 1.0670825e-10_RKIND*tev**2.012_RKIND
     &          / (exp(4.463_RKIND/tev)
     &          * (1._RKIND+0.2472_RKIND*tev)**3.512_RKIND) / kunit
c
          k11a(i) = exp(-24.24914687731536_RKIND
     &            + 3.400824447095291_RKIND*logtev
     &            - 3.898003964650152_RKIND*logtev**2
     &            + 2.045587822403071_RKIND*logtev**3
     &            - 0.5416182856220388_RKIND*logtev**4
     &            + 0.0841077503763412_RKIND*logtev**5
     &            - 0.007879026154483455_RKIND*logtev**6
     &            + 0.0004138398421504563_RKIND*logtev**7
     &            - 9.36345888928611e-6_RKIND*logtev**8) / kunit

          k12a(i) = 4.38e-10_RKIND*exp(-102000._RKIND/ttt)
     &         * ttt**0.35_RKIND / kunit
        else
          k13a(i) = tiny
          k11a(i) = tiny
          k12a(i) = tiny
        endif
c
        if (tev .gt. 0.04_RKIND) then
          k14a(i) = exp(-18.01849334273_RKIND
     &            + 2.360852208681_RKIND*logtev
     &            - 0.2827443061704_RKIND*logtev**2
     &            + 0.01623316639567_RKIND*logtev**3
     &            - 0.03365012031362999_RKIND*logtev**4
     &            + 0.01178329782711_RKIND*logtev**5
     &            - 0.001656194699504_RKIND*logtev**6
     &            + 0.0001068275202678_RKIND*logtev**7
     &            - 2.631285809207e-6_RKIND*logtev**8) / kunit
        else
          k14a(i) =  tiny
        endif
c
        if (tev .gt. 0.1_RKIND) then
          k15a(i) = exp(-20.37260896533324_RKIND
     &            + 1.139449335841631_RKIND*logtev
     &            - 0.1421013521554148_RKIND*logtev**2
     &            + 0.00846445538663_RKIND*logtev**3
     &            - 0.0014327641212992_RKIND*logtev**4
     &            + 0.0002012250284791_RKIND*logtev**5
     &            + 0.0000866396324309_RKIND*logtev**6
     &            - 0.00002585009680264_RKIND*logtev**7
     &            + 2.4555011970392e-6_RKIND*logtev**8
     &            - 8.06838246118e-8_RKIND*logtev**9) / kunit
        else
          k15a(i) = 2.56e-9_RKIND*tev**1.78186_RKIND / kunit
        endif
c
        k16a(i) = 6.5e-9_RKIND/sqrt(tev) / kunit
c
        k17a(i) = 1.0e-8_RKIND*ttt**(-0.4_RKIND) / kunit
        if(ttt.gt.1.0d4)
     &     k17a(i)=4.0e-4_RKIND*ttt**(-1.4_RKIND)*exp(-15100._RKIND/ttt) / kunit
c
        k18a(i) = 5.56396e-8_RKIND/tev**0.6035_RKIND / kunit
        k19a(i) = 4.64e-8_RKIND/sqrt(tev) / kunit
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
          if (ttt .le. 300._RKIND) then
             k22a(i) = 1.3e-32_RKIND * (ttt/300._RKIND)**(-0.38_RKIND)
          else
             k22a(i) = 1.3e-32_RKIND * (ttt/300._RKIND)**(-1._RKIND)
          endif
        elseif(threebody.eq.1)then
c Inverse of PSS83 three-body rate
          k22a(i) = 5.5e-29_RKIND / ttt
          k13a(i) = (5.24e-7_RKIND / ttt**0.485_RKIND) 
     &                   * exp(-5.2e4_RKIND / ttt)
        elseif(threebody.eq.2)then
c Inverse of CW83 rate
          k22a(i) = 8.8e-33_RKIND
          k13a(i) = 8.4e-11_RKIND * ttt**0.515_RKIND
     &         * exp(-5.2e4_RKIND / ttt)
        elseif(threebody.eq.3)then
c Rate from JGC67, used by FH07 to construct their three-body rate
          k22a(i) = 1.44e-26_RKIND / ttt**1.54_RKIND
          k13a(i) = (1.38e-4_RKIND / ttt**1.025_RKIND)
     &         * exp(-5.2e4_RKIND / ttt)
        elseif(threebody.eq.4)then
c High density rate from MSM96
          k22a(i) = 7.7e-31_RKIND / ttt**0.464_RKIND
          k13a(i) = 10._RKIND**(-178.4239_RKIND - 68.42243_RKIND 
     &         * log10(ttt) + 43.20243_RKIND * log10(ttt)**2
     &         - 4.633167_RKIND * log10(ttt)**3 
     &         + 69.70086_RKIND * log10(1_RKIND + 40870.38_RKIND / ttt)
     &         - (23705.7_RKIND / ttt))
        else
            write(0,*) "I didn't understand your threebody rate."
            stop
        endif
        k13a(i) = k13a(i) / kunit
        k22a(i) = k22a(i) / kunit_3bdy
        
        k21a(i) = 2.8e-31_RKIND * (ttt**(-0.6_RKIND)) / kunit_3bdy
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
        k50a(i) = 1.0e-9_RKIND *exp(-4.1e1_RKIND/ttt)  / kunit
        k51a(i) = 1.0e-9_RKIND                   / kunit
        k52a(i) = 2.1e-9_RKIND                   / kunit
        k53a(i) = 1.0e-9_RKIND *exp(-4.57e2_RKIND/ttt) / kunit
        k54a(i) = 7.5e-11_RKIND*exp(-3.82e3_RKIND/ttt) / kunit
        k55a(i) = 7.5e-11_RKIND*exp(-4.24e3_RKIND/ttt) / kunit
        k56a(i) = 1.5e-9_RKIND *(ttt/300._RKIND)**(-0.1_RKIND) / kunit
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
        ceHIa(i) = 7.5e-19_RKIND*exp(-min(log(dhuge),
     &        118348._RKIND/ttt))/(1._RKIND+sqrt(ttt/1.e5_RKIND)) 
     &        / coolunit
        ceHeIa(i) = 9.1e-27_RKIND*exp(-min(log(dhuge),13179._RKIND/ttt))
     &       *ttt**(-0.1687_RKIND)/(1._RKIND+sqrt(ttt/1.e5_RKIND)) 
     &       / coolunit
        ceHeIIa(i) = 5.54e-17_RKIND*exp(-min(log(dhuge),
     &       473638._RKIND/ttt))*ttt**(-0.397_RKIND) 
     &       / (1._RKIND+sqrt(ttt/1.0d5)) / coolunit
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
c        ciHIa(i)    = 1.27e-21_RKIND*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-min(log(dhuge),157809.1/ttt)) / coolunit
c        ciHeIa(i)   = 9.38e-22_RKIND*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-min(log(dhuge),285335.4/ttt)) / coolunit
c        ciHeIIa(i)  = 4.95e-22_RKIND*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-min(log(dhuge),631515.0/ttt)) / coolunit
c
        ciHeISa(i)  = 5.01e-27_RKIND*(ttt)**(-0.1687_RKIND) 
     &        / (1._RKIND+sqrt(ttt/1.e5_RKIND))
     &        *exp(-min(log(dhuge),55338._RKIND/ttt)) / coolunit
c
c       Collisional ionizations (polynomial fits from Tom Abel)
c
        ciHIa(i)    = 2.18e-11_RKIND*k1a(i)*kunit / coolunit
        ciHeIa(i)   = 3.94e-11_RKIND*k3a(i)*kunit / coolunit
        ciHeIIa(i)  = 8.72e-11_RKIND*k5a(i)*kunit / coolunit
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
        reHIIa(i)   = 8.70e-27_RKIND*sqrt(ttt)
     &        * (ttt/1000._RKIND)**(-0.2_RKIND)
     &        / (1._RKIND + (ttt/1.0d6)**(0.7_RKIND)) / coolunit
        reHeII1a(i) = 1.55e-26_RKIND*ttt**0.3647_RKIND / coolunit
        reHeII2a(i) = 1.24e-13_RKIND*ttt**(-1.5_RKIND) *
     &       exp(-min(log(dhuge),470000._RKIND/ttt)) *
     &       (1._RKIND+0.3_RKIND*exp(-min(log(dhuge),94000._RKIND/ttt))) 
     &       / coolunit
        reHeIIIa(i) = 3.48e-26_RKIND*sqrt(ttt)
     &       * (ttt/1000._RKIND)**(-0.2_RKIND)
     &       / (1._RKIND + (ttt/1.d6)**(0.7_RKIND)) / coolunit
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
        brema(i)    = 1.43e-27_RKIND*sqrt(ttt)*(1.1_RKIND+0.34_RKIND
     &        *exp(-(5.5_RKIND-log10(ttt))**2/3._RKIND)) / coolunit
      else
        brema(i)    = tiny
      endif
c
c         Bremsstrahlung (Shapiro & Kang 1987)(Spitzer 1978)
c
c        if(ttt .lt. 3.2e5_RKIND) gaunt = 0.79464_RKIND + 0.1243_RKIND*log10(ttt)
c        if(ttt .ge. 3.2e5_RKIND) gaunt = 2.13164 - 0.1240_RKIND*log10(ttt)
c        brem1a(i) = 1.426e-27_RKIND*sqrt(ttt)*gaunt / coolunit
c        if(ttt .lt. 4.0*3.2e5_RKIND) gaunt = 0.79464_RKIND + 0.1243_RKIND*log10(ttt/4._RKIND)
c        if(ttt .ge. 4.0*3.2e5_RKIND) gaunt = 2.13164_RKIND - 0.1240_RKIND*log10(ttt/4._RKIND)
c        brem2a(i) = 1.426e-27_RKIND*sqrt(ttt)*gaunt / coolunit
c
c    e) Molecular hydrogen cooling
c
c       (Which one is used is controlled by a flag in cool1d_mulit.src)
c
c    e part1) - Lepp & Shull rates
c
        xx = log10(ttt/1.e4_RKIND)
        vibha  (i) = 1.1e-18_RKIND*exp(-min(log(dhuge),6744._RKIND/ttt)) 
     &               / coolunit
c
        if( ttt .gt. 1635._RKIND) then
           dum = 1.0e-12_RKIND*sqrt(ttt)*exp(-1000._RKIND/ttt)
        else
           dum = 1.4e-13_RKIND*exp((ttt/125._RKIND)-(ttt/577._RKIND)**2)
        endif
        hyd01ka(i) = dum*exp(-min(log(dhuge),
     &                            8.152e-13_RKIND/(kboltz*ttt) )) 
     &               / coolunit
c
        dum = 8.152e-13_RKIND*(4.2_RKIND /
     &       (kboltz*(ttt+1190._RKIND)) +
     &       1._RKIND/(kboltz*ttt))
        h2k01a(i) = 1.45e-12_RKIND*sqrt(ttt)*exp(-min(log(dhuge),dum)) 
     &              / coolunit
c
        if(ttt .gt. 4031._RKIND)then
          rotla(i) = 1.38e-22_RKIND*exp(-9243._RKIND/ttt) / coolunit
        else
          rotla(i) = 10._RKIND**(-22.9_RKIND - 0.553_RKIND*xx 
     &          - 1.148_RKIND*xx*xx) / coolunit
        endif
c
        if(ttt .gt. 1087._RKIND)then
           rotha(i) = 3.9e-19_RKIND*exp(-6118._RKIND/ttt) / coolunit
        else
           rotha(i) = 10._RKIND**(-19.24_RKIND + 0.474_RKIND*xx
     &          - 1.247_RKIND*xx*xx) / coolunit
        endif
c
c     e part2) - Galli and Palla (1999) rates as fit by Tom Abel
c
        tm = max(ttt, 13._RKIND)       ! no cooling below 13 Kelvin ...
        tm = min(tm, 1.e5_RKIND)       ! fixes numerics
        lt = log10(tm)
c     low density limit from Galli and Palla
        gpldl(i)  = 10._RKIND**(-103._RKIND+97.59_RKIND*lt
     &       - 48.05_RKIND*lt**2+10.8_RKIND*lt**3-0.9032_RKIND*lt**4) 
     &       / coolunit
c     high density limit from HM79 (typo corrected Aug 30/2007)
        t3 = tm/1000._RKIND
        HDLR = ((9.5e-22_RKIND*t3**3.76_RKIND)
     &       / (1._RKIND+0.12_RKIND*t3**2.1_RKIND)*
     &       exp(-(0.13_RKIND/t3)**3)+3.e-24_RKIND*exp(-0.51_RKIND/t3))
        HDLV = (6.7e-19_RKIND*exp(-5.86_RKIND/t3) + 1.6e-18_RKIND 
     &       * exp(-11.7_RKIND/t3))
        gphdl(i)  = (HDLR + HDLV) / coolunit
c
c     e part 3) - Low density rates from Glover & Abel (2008)
c
c  Excitation by HI
c
      tm  = max(ttt, 10._RKIND)
      tm  = min(tm, 1.e4_RKIND)
      lt3 = log10(tm / 1.e3_RKIND)  

      if (tm .lt. 1.e2_RKIND) then
        gaHIa(i) = 10._RKIND**(-16.818342_RKIND
     &         + 37.383713_RKIND * lt3
     &         + 58.145166_RKIND * lt3**2
     &         + 48.656103_RKIND * lt3**3
     &         + 20.159831_RKIND * lt3**4
     &         + 3.847961_RKIND * lt3**5) / coolunit
      elseif (tm .lt. 1.e3_RKIND) then
        gaHIa(i) = 10._RKIND**(-24.311209_RKIND
     &         + 3.5692468_RKIND * lt3
     &         - 11.33286_RKIND * lt3**2
     &         - 27.850082_RKIND * lt3**3
     &         - 21.328264_RKIND * lt3**4
     &         - 4.2519023_RKIND * lt3**5) / coolunit
      else
        gaHIa(i) = 10._RKIND**(-24.311209_RKIND
     &         + 4.6450521_RKIND * lt3
     &         - 3.7209846_RKIND * lt3**2
     &         + 5.9369081_RKIND * lt3**3
     &         - 5.5108047_RKIND * lt3**4
     &         + 1.5538288_RKIND * lt3**5) / coolunit
      endif
c
c Excitation by H2
c
      gaH2a(i) = 10._RKIND**(-23.962112_RKIND
     &        + 2.0943374_RKIND  * lt3
     &        - 0.77151436_RKIND  * lt3**2
     &        + 0.43693353_RKIND  * lt3**3
     &        - 0.14913216_RKIND  * lt3**4
     &        - 0.033638326_RKIND * lt3**5) / coolunit
c
c Excitation by He
c
      gaHea(i) = 10._RKIND**(-23.689237_RKIND
     &        + 2.1892372_RKIND  * lt3
     &        - 0.81520438_RKIND * lt3**2
     &        + 0.29036281_RKIND * lt3**3
     &        - 0.16596184_RKIND * lt3**4
     &        + 0.19191375_RKIND * lt3**5) / coolunit
c
c Excitation by H+
c
      gaHpa(i) = 10._RKIND**(-21.716699_RKIND
     &        + 1.3865783_RKIND   * lt3
     &        - 0.37915285_RKIND  * lt3**2
     &        + 0.11453688_RKIND  * lt3**3
     &        - 0.23214154_RKIND  * lt3**4
     &        + 0.058538864_RKIND * lt3**5) / coolunit
c
c Excitation by electrons
c
      if (tm .lt. 200._RKIND) then
        gaela(i) = 10._RKIND**(-34.286155_RKIND
     &          - 48.537163_RKIND  * lt3
     &          - 77.121176_RKIND  * lt3**2
     &          - 51.352459_RKIND  * lt3**3
     &          - 15.16916_RKIND  * lt3**4
     &          - 0.98120322_RKIND * lt3**5) / coolunit
      else
        gaela(i) = 10._RKIND**(-22.190316
     &          + 1.5728955_RKIND  * lt3
     &          - 0.213351_RKIND * lt3**2
     &          + 0.96149759_RKIND * lt3**3
     &          - 0.91023195_RKIND * lt3**4
     &          + 0.13749749_RKIND * lt3**5) / coolunit
      endif
c
c     f) HD cooling
c
c        HD Cooling Function (ergs cm3 /s)
c
c Fit to Lepp and Shull 1984 LTE (ergs/s) -> hdlte (ergs cm3/s)
c
      hdltea(i) = -35.6998_RKIND + 15.35716_RKIND*dlog10(ttt) -
     &            5.58513_RKIND   * (dlog10(ttt))**2 +
     &            0.8561149_RKIND * (dlog10(ttt))**3 - 
     &            1.75538e-2_RKIND  * (dlog10(ttt))**4
      hdltea(i) = (10._RKIND**min(hdltea(i), 0._RKIND)) / coolunit
c      hdlte=10**lhdlte/den
c
c Galli and Palla 1998 low density limit (erg cm3 /s)
c  uses HD-He collisional data, so reduced by 1.27
c
      hdlowa(i) = ((3._RKIND * (4.4e-12_RKIND 
     &     + 3.6e-13_RKIND*ttt**0.77_RKIND) *
     &     dexp(-128._RKIND/ttt) * 128._RKIND +
     &     (5._RKIND/3._RKIND) * (4.1e-12_RKIND +
     &     2.1e-13_RKIND*ttt**0.92_RKIND) *
     &     dexp(-255._RKIND/ttt) * 255._RKIND) * kb/1.27_RKIND ) 
     &     / coolunit
c      hdcool=hdlte/(1._RKIND+hdlte/hdlow)
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
      ciecoa(i) = cierate * (mh/2._RKIND) / coolunit
c      write(6,*), "thin: ", ttt, cierate * 1e14*mh*2
c

      enddo
c
c     h) Compton cooling (Peebles 1971)
c
      compa     = 5.65e-36_RKIND / coolunit
c
c     i) Photoelectric heating by UV-irradiated dust (Wolfire 1995)
c        with epsilon=0.05, G_0=1.7 (rate in erg s^-1 cm^-3)
c
      gammaha   = 8.5e-26_RKIND / coolunit
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
        nu = 10._RKIND**((n-14)*0.01_RKIND)              ! 0.74 ev -- 7.24d9 ev
c
c       24) HI photo-ionization cross-section
c
        if (nu .gt. e24) then
           dum = sqrt(nu/e24-1._RKIND)
           sigma24(n) = 6.3e-18_RKIND * (e24/nu)**4
     &               *  exp(4._RKIND-4._RKIND*atan(dum)/dum)
     &               / (1._RKIND-exp(-2._RKIND*pi/dum))
        else
           sigma24(n) = 0._RKIND
        endif
c
c       25) HeII photo-ionization cross-section
c
        if (nu .gt. e25) then
           dum = sqrt(nu/e25-1._RKIND)
           sigma25(n) = 1.58e-18_RKIND * (e25/nu)**4
     &               * exp(4._RKIND-4._RKIND*atan(dum)/dum)
     &               / (1._RKIND-exp(-2._RKIND*pi/dum))
        else
           sigma25(n) = 0._RKIND
        endif
c
c       26) HeI photo-ionization cross-section
c
        if (nu .gt. e26) then
c           sigma26(n) = 7.42e-18_RKIND*(1.66*(nu/e26)**(-2.05)
c     &                          - 0.66*(nu/e26)**(-3.05))
c
c            Now From Verland, et al (1996)
c
           x = nu/13.61_RKIND - 0.4434_RKIND
           y = sqrt(x**2 + 2.136_RKIND**2)
           sigma26(n) = 9.492e-16_RKIND*((x-1)**2 + 2.039_RKIND**2) *
     &                  y**(0.5_RKIND*3.188_RKIND-5.5_RKIND) *
     &                  (1._RKIND + sqrt(y/1.469_RKIND))**(-3.188_RKIND)        
        else
           sigma26(n) = 0._RKIND
        endif
c
c       27) HM    + p   -> HI    + e
c
        if (nu .gt. e27) then
          sigma27(n) = 2.11e-16_RKIND*(nu-e27)**1.5_RKIND/nu**3
        else
          sigma27(n) = 0._RKIND
        endif
c
c       28) H2II  + p   -> HI    + HII
c
        if (nu .gt. e28a .and. nu .le. e28b) then
          sigma28(n) = 10._RKIND**(-40.97_RKIND+6.03_RKIND*nu
     &          -0.504_RKIND*nu**2+1.387e-2_RKIND*nu**3)
        elseif (nu .gt. e28b .and. nu .lt. e28c) then
           sigma28(n) = 10._RKIND**(-30.26_RKIND+2.79_RKIND*nu
     &          -0.184_RKIND*nu**2+3.535e-3_RKIND*nu**3)
        else
          sigma28(n) = 0._RKIND
        endif
c
c       29) H2I   + p   -> H2II  + e
c
        if (nu .gt. e29a .and. nu .le. e29b) then
          sigma29(n) = 6.2e-18_RKIND*nu - 9.4e-17_RKIND
        elseif (nu .gt. e29b .and. nu .le. e29c) then
          sigma29(n) = 1.4e-18_RKIND*nu - 1.48e-17_RKIND
        elseif (nu .gt. e29c) then
          sigma29(n) = 2.5e-14_RKIND*nu**(-2.71_RKIND)
        else
          sigma29(n) = 0._RKIND
        endif
c
c       30) H2II  + p   -> 2HII  + e
c
        if (nu .ge. e30a .and. nu .lt. e30b) then
           sigma30(n) = 10._RKIND**(-16.926_RKIND-4.528e-2_RKIND*nu
     &          +2.238e-4_RKIND*nu**2+4.245e-7_RKIND*nu**3)
        else
          sigma30(n) = 0._RKIND
        endif
c
c       31) H2I   + p   -> 2HI
c
c        sigma31(n) = 0._RKIND

C	 BUG FIX: (n) subscripts below were (i)
C	 reported by Patrick McDonald, modified in revision 2310

        if (nu .ge. e28b .and. nu .lt. e24) THEN
           sigma31(n) = 3.71e-18_RKIND
        else
           sigma31(n) = 0._RKIND
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
           fext(n) = f3*(nu/12.86_RKIND)**(-1._RKIND)*exp(-1._RKIND*
     &       min(1.e22_RKIND*(sigma24(n) + 0.08_RKIND*sigma26(n)), 
     &           100._RKIND))
        else
c
c          otherwise) no radiation field
c
           fext(n) = 0._RKIND
        endif
      enddo
c
c     Integrate over the frequency spectrum
c

!     call open_mpi_error_file( 'F89', 89, 'unknown' )

      do n = 1, nnu
c          delnu = (10._RKIND**((n-14)*0.01_RKIND) - 10._RKIND**((n-15)*0.01_RKIND))*evhz
c          hnu   = 0.5_RKIND*(10._RKIND**((n-14)*0.01_RKIND) + 10._RKIND**((n-15)*0.01_RKIND))*everg

          delnu = (10._RKIND**((n-13.5_RKIND)*0.01_RKIND) - 
     &        10._RKIND**((n-14.5_RKIND)*0.01_RKIND))*evhz
          hnu   = 10._RKIND**((n-14)*0.01_RKIND)*everg
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
          dum = max(hnu - e24*everg, 0._RKIND)
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
         k24 = k24 + 0.39_RKIND*(piHI + 0.08_RKIND*piHeI)/(e24*everg)
     &         * coolunit * tbase1
         k26 = k26 + 0.055_RKIND*(piHI/0.08_RKIND + piHeI)/(e26*everg)
     &         * coolunit * tbase1
c
c        Now correct the photo-heating rates to reflect that much of
c          the energy does not go into heating. (ignore HeII)
c
         piHI  = piHI *0.13_RKIND
         piHeI = piHeI*0.13_RKIND
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

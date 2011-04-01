#include "fortran.def"
#define EQOFFSET 1.5
#define EQOVERLAP 3.0
#define TOLERANCE 1e-20
c=======================================================================
c/////////////////////  SUBROUTINE SOLVE_RATE  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine primordial_solver(d, e, ge, u, v, w, de,
     &                HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, ispecies, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht,
     &                dt, aye, temstart, temend, 
     &                utem, uxyz, uaye, urho, utim, uvel,
     &                eta1, eta2, gamma, fh, dtoh,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a,
     &                k16a, k17a, k18a, k19a, k21a, k22a, k23a,
     &                k24, k25, k26, k27, k28, k29, k30, k31,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, 
     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
     &                gpldla, gphdla, hdltea, hdlowa, hdcoola, ciecoa, 
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                inutot, iradtype, nfreq, 
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iciecool, ih2optical, errcode, omaskflag,
     &                threebody, subgridmask
#ifdef UNUSED_TABULATED_EQ
     &                ,HIeqtable, HIIeqtable, H2Ieqtable,
     &                nrhobins, nebins,
     &                rhostart, rhoend, estart, eend
#endif
#ifdef PYFORT
     &                ,tgasouti, tgasoutf, geout
#endif
     &                )
c
c  SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
c
c  written by: Matt Turk, advisor: Tom Abel
c  API and inspiration taken from module by
c   Greg Bryan, Yu Zhang, Peter Anninos and Tom Abel
c
c  PURPOSE:
c    Solve the multi-species rate and cool equations.
c
c  INPUTS:
c    in,jn,kn - dimensions of 3D fields
c
c    d        - total density field
c    de       - electron density field
c    HI,HII   - H density fields (neutral & ionized)
c    HeI/II/III - He density fields
c    DI/II    - D density fields (neutral & ionized)
c    HDI      - neutral HD molecule density field
c    HM       - H- density field
c    H2I      - H_2 (molecular H) density field
c    H2II     - H_2+ density field
c
c    is,ie    - start and end indices of active region (zero based)
c    idual    - dual energy formalism flag (0 = off, 1 = on)
c    iexpand  - comoving coordinates flag (0 = off, 1 = on)
c    idim     - dimensionality (rank) of problem
c    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
c    iradshield - flag for crude radiative shielding correction 
c    iradtype - type of radiative field (only used if = 8)
c    imethod  - Hydro method (0 = PPMDE, 2 = ZEUS-type)
c    ih2co    - flag to include H2 cooling (1 = on, 0 = off)
c    ipiht    - flag to include photoionization heating (1 = on, 0 = off)
c    iciecool - flag to include CIE cooling (1 = on, 0 = off)
c    ih2optical - flag to include H2 optical depth fudge (1 = on, 0 = off)
c
c    fh       - Hydrogen mass fraction (typically 0.76)
c    dtoh     - Deuterium to H mass ratio
c    dt       - timestep to integrate over
c    aye      - expansion factor (in code units)
c
c    utim     - time units (i.e. code units to CGS conversion factor)
c    uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
c    urho     - density units
c    uxyz     - length units
c    utem     - temperature(-like) units
c
c    temstart, temend - start and end of temperature range for rate table
c    nratec   - dimensions of chemical rate arrays (functions of temperature)
c
c  OUTPUTS:
c    update chemical rate densities (HI, HII, etc)
c
c  PARAMETERS:
c    itmax   - maximum allowed sub-cycle iterations
c    mh      - H mass in cgs units
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c  General Arguments
c
      integer in, jn, kn, is, js, ks, ie, je, ke, nratec, imethod,
     &        idual, iexpand, ih2co, ipiht, ispecies, idim,
     &        iradtype, nfreq, iradshield, iciecool, 
     &        ih2optical, errcode, omaskflag, threebody
      real    aye, temstart, temend, eta1, eta2, gamma,
     &        utim, uxyz, uaye, urho, utem, fh, dtoh, uvel
#ifdef PYFORT
      real    dt(in,jn,kn)
#else
      real    dt
#endif    
c
c  Density, energy and velocity fields fields
c
CF2PY INTENT(INOUT) ::  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn)
CF2PY INTENT(INOUT) :: HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
CF2PY INTENT(INOUT) ::  HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
CF2PY INTENT(INOUT) ::  DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
CF2PY INTENT(INOUT) ::  d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn)
CF2PY INTENT(INOUT) ::  u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn)
CF2PY INTENT(INOUT) ::  subgridmask(in,jn,kn)
#ifdef UNUSED_TABULATED_EQ
CF2PY INTENT(INOUT) ::  HIeqtable(nebins,nrhobins)
CF2PY INTENT(INOUT) ::  HIIeqtable(nebins,nrhobins)
CF2PY INTENT(INOUT) ::  H2Ieqtable(nebins,nrhobins)
#endif
CF2PY INTENT(INOUT) ::  tgasouti(in,jn,kn)
CF2PY INTENT(INOUT) ::  tgasoutf(in,jn,kn)
CF2PY INTENT(INOUT) ::  geout(in,jn,kn)
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        tgasouti(in,jn,kn)
#ifdef PYFORT
     &        ,tgasoutf(in,jn,kn),geout(in,jn,kn)
#endif
#ifdef UNUSED_TABULATED_EQ
      integer nebins,nrhobins
      real    HIeqtable(nrhobins,nebins), HIIeqtable(nrhobins,nebins),
     &        H2Ieqtable(nrhobins,nebins),
     &        rhostart, rhoend,
     &        estart, eend
#endif
#ifdef PYFORT
      integer subgridmask(in,jn,kn)
#else
      integer subgridmask
#endif
      integer subgridtempmask(in,jn,kn)
c
c  Cooling tables (coolings rates as a function of temperature)
c
      real    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec),
     $        hdcoola(nratec, 5), ciecoa(nratec)
      real    gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec)
      real    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      real    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p
c
c  Chemistry tables (rates as a function of temperature)
c
      real k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k21a(nratec),
     &     k22a(nratec), k23a(nratec), k50a(nratec), k51a(nratec), 
     &     k52a(nratec), k53a(nratec), k54a(nratec), k55a(nratec), 
     &     k56a(nratec), k13dda(nratec, 7),
     &     k24, k25, k26, k27, k28, k29, k30, k31
      
c
c  Parameters
c
      integer itmax, ijk
      parameter (itmax = 1e8, ijk = MAX_ANY_SINGLE_DIRECTION)
      double precision mh
      parameter (mh = 1.67d-24)
c
c  Locals
c
      integer i, j, k, n, m, iter, inneriter, imin, ifedtgas
      real ttmin, dom, energy(ijk), comp1, comp2, maxchange, fH2, dfH2,
     &     scaling, prevfH2, HIfrac(ijk), HIIfrac(ijk), H2Ifrac(ijk)
      double precision coolunit, dbase1, tbase1, xbase1, kunit,
     &                 kunit_3bdy,
     &                  chunit, ciHIunit, ciHeIunit, ciHeIIunit

c  Stuff for TWOSTEP
      real*8 s(ijk,NSPECIES), sp(ijk,NSPECIES), sp1(ijk,NSPECIES),
     &       dsp(ijk,NSPECIES), dsq(ijk,NSPECIES), fd(ijk),
     &       temp1, temp2, temp3, kwq, kws,
     &       tsm1, tsp, ts, ty, dsqge, dspge,
     &       temptsc2, temp2tsc, temptsc12
      real*8 qt, sm1(ijk,NSPECIES), 
     &       tsgamma(ijk), tsc(ijk), tscapy(ijk,NSPECIES)
      real*8 temp, dtot(ijk), rtot(ijk), h2heatp(ijk), h2heatm(ijk),
     &       hiheatp(ijk), heiheatp(ijk), heiiheatp(ijk),
     &       esterr(ijk), itererror
      integer iHI,iHII,ide,iHeI,iHeII,iHeIII,iH2I,iH2II,iHM,ige,
     &        iDI,iDII,iHDI,ncspecies,
     &        coriter, totaliter, nciter, oreas
      real*8 dtcoef(NSPECIES)
      logical capycare(NSPECIES), eqable(NSPECIES)
      parameter (ncspecies=11, nciter=1)
c
c  row temporaries
c 
      integer indixe(ijk)
      real t1(ijk), t2(ijk), logtem(ijk), tdef(ijk),
     &     dtit(ijk), ttot(ijk), p2d(ijk), tgas(ijk), tgasold(ijk),
     &     olddtit(ijk), correctH(ijk), 
     &     correctHe(ijk), correctD(ijk), geprime(ijk), dprime(ijk)
c
c  Rate equation row temporaries
c
      real gamma2(ijk), tto, scr1, scr2, ascr
      double precision tgastemp
      real*8 H2Idot_prev(ijk)
      real k24shield(ijk), k25shield(ijk), k26shield(ijk)
      real r1 (ijk), r2 (ijk), r3 (ijk), r4 (ijk), r5 (ijk),
     &     r6 (ijk), r7 (ijk), r8 (ijk), r9 (ijk), r10(ijk),
     &     r11(ijk), r12(ijk), r13(ijk), r14(ijk), r15(ijk),
     &     r16(ijk), r17(ijk), r18(ijk), r19(ijk), r21(ijk),
     &     r22(ijk), r23(ijk), r50(ijk), r51(ijk), r52(ijk), 
     &     r53(ijk), r54(ijk), r55(ijk), r56(ijk),
     &     r24, r25, r26, r27, r28, r29, r30, r31
      real k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k21(ijk),
     &     k22(ijk), k23(ijk), k50(ijk), k51(ijk), k52(ijk), 
     &     k53(ijk), k54(ijk), k55(ijk), k56(ijk), k13dd(ijk, 7)
c
c  Cooling/heating row locals
c
      double precision ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &     ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &     reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk),
     &     brem(ijk), edot(ijk), 
     &     edotplus(ijk), edotminus(ijk)
      real hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), rotl(ijk),
     &     gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk), 
     &     hdcool(ijk, 5), cieco(ijk)

#ifdef UNUSED_TABULATED_EQ
      real estartcode, eendcode, rhostartcode, rhoendcode
#endif

c
c     Iteration mask, output mask, and output flag
c     
      logical itmask(ijk), omask(ijk), eqmask(ijk), dmask(ijk),
     &         rejected(ijk), eulermask(ijk),
     &         eqflag, itflag, miscflag, doupdate

      integer totaleq, mins(ijk)

      double precision everg, e24, e26
      real dlogrho, dloge
      parameter(everg = 1.60184d-12, e24 = 13.6d0, e26 = 24.6d0)
      integer totalit
      parameter(iHI=1,iHII=2,iHeI=3,iHeII=4,iHeIII=5,
     &          iH2I=6,ige=7,iDI=8,iDII=9,iHDI=10,
     &          ide=11,iHM=12,iH2II=13 ) ! eq species
      totaliter = 0
c      iciecool = 0(
c      ih2optical = 0

c Initialize the dt coefficients and say whether or not
c to do this second-order BDF

      eqable(iHI   )    =  .true.
      eqable(iHII  )    =  .true.
      eqable(ide   )    =  .false.
      eqable(iHeI  )    =  .false.
      eqable(iHeII)     =  .false.
      eqable(iHeIII)    =  .false.
      eqable(iH2I)      =  .true.
      eqable(ige)       =  .false.
      eqable(iDI)       =  .false.
      eqable(iDII)      =  .false.
      eqable(iHDI)      =  .false.
      eqable(iH2II)     =  .false.
      eqable(iHM)       =  .false.

      dtcoef(iHI)      =  0.1
      capycare(iHI)    = .true.
      dtcoef(iHII)     =  0.1 
      capycare(iHII)   = .true.
      dtcoef(ide)      =  0.1
      capycare(ide)    = .false.
      dtcoef(iHeI)     =  0.1
      capycare(iHeI)   = .false.
      dtcoef(iHeII)    =  0.1
      capycare(iHeII)  = .false.
      dtcoef(iHeIII)   =  0.1
      capycare(iHeIII) = .false.
      dtcoef(iH2I)     =  0.01
      capycare(iH2I)   = .true.
      dtcoef(ige)      =  0.01
      capycare(ige)    = .true.
      dtcoef(iDI)      =  1.0
      capycare(iDI)    = .false.
      dtcoef(iDII)     =  1.0
      capycare(iDII)   = .false.
      dtcoef(iHDI)     =  1.0
      capycare(iHDI)   = .false. ! Off for now; should track H2 directly
      dtcoef(iH2II)    =  1.0
      capycare(iH2II)  = .false.
      dtcoef(iHM)      =  1.0
      capycare(iHM)    = .false.

c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Set units
c
      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      kunit    = (uaye**3 * mh) / (dbase1 * tbase1)
      kunit_3bdy  = kunit * (uaye**3 * mh) / dbase1
c     4.48 eV = 7.17775e-12 ergs
      chunit = (7.17775e-12)/(2.0*uvel*uvel*mh)
      ciHIunit    = (2.18e-11)/(uvel*uvel*mh)
      ciHeIunit   = (3.94e-11)/(4.0*uvel*uvel*mh)
      ciHeIIunit  = (8.72e-11)/(4.0*uvel*uvel*mh)

#ifdef UNUSED_TABULATED_EQ
c      rhostartcode = log10((10**rhostart)/(urho*aye**3.0))
c      rhoendcode = log10((10**rhoend)/(urho*aye**3.0))
      estartcode = log10((10**estart)/(uvel**2.0))
      eendcode = log10((10**eend)/(uvel**2.0))
      rhostartcode = log10(10**rhostart/(fh*dom))
      rhoendcode = log10(10**rhoend/(fh*dom))
#endif

#ifdef PYFORT_OFF
      write(0,*) 'kunit',kunit
      write(0,*) 'kunit_3bdy',kunit_3bdy
      write(0,*) 'utim',utim
      write(0,*) 'dom',dom
      write(0,*) 'urho',urho
      write(0,*) 'aye',aye
      write(0,*) 'uaye',uaye
      write(0,*) 'uxyz',uxyz
#endif

 1000    format(I3,3X,50(e14.7,1x))
 1001    format(5(I10,3X),80(e24.17,1x))
c
c  Correct for flaws in flavor conservation and
c  convert densities from comoving to proper
c
      call primordial_make_consistent(
     &                     de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, d,
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, fh, dtoh, iter,
     &                     correctH, correctHe )
      call primordial_scale_fields(
     &                  d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, 
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, aye**(-3))
c
c  Loop over zones, and do an entire i-column in one go
c
      ifedtgas = 0
      totalit = 0
      totaleq = 0
      eqflag = .false.
      if(omaskflag.eq.1)write(0,*)'omaskflag=1'
      do k = ks+1, ke+1
      do j = js+1, je+1
c
c        Initialize variables
c
         iter = 1
         do i = is+1, ie+1
            ttot(i) = 0.0
            omask(i) = .false.
            rejected(i) = .false.
            fd(i) = fh*d(i,j,k)
#ifdef PYFORT
            if(subgridmask(i,j,k).eq.1
     &        .and.dt(i,j,k).gt.0.0d0)then
              itmask(i)=.true.
            else
               itmask(i)=.false.
            endif
            dtit(i) = dt(i,j,k)/5.0
#else
            dtit(i) = dt
            subgridmask = 1
            itmask(i) = .true.
#endif
            eqmask(i) = .false.
#ifdef UNUSED_TABULATED_EQ
            dprime(i) = -35
            geprime(i) = -35
            geprime(i) = log10((ge(i,j,k) +
     &                 (fh - H2I(i,j,k)/d(i,j,k))*chunit)*(uvel**2.0))
            dprime(i) = log10(d(i,j,k)*dom*fh)
c            write(0,*) rhostartcode, rhostart, estartcode, estart
            if ((dprime(i).gt.rhostart+3.5).and.(geprime(i).gt.estart))then
                totaleq = totaleq + 1
c                eqmask(i) = .true.
c                eqflag = .true.
c                write(0,*) geprime(i), estartcode, dprime(i), rhostartcode
c                write(0,*) (geprime(i).gt.estartcode)
c                stop
            endif
#endif
            s(i,iHI) = max(HI(i,j,k),tiny)
            s(i,iHII) = max(HII(i,j,k),tiny)
            s(i,ide) = max(de(i,j,k),tiny)
            s(i,iHeI) = max(HeI(i,j,k),tiny)
            s(i,iHeII) = max(HeII(i,j,k),tiny)
            s(i,iHeIII) = max(HeIII(i,j,k),tiny)
            s(i,iH2I) = max(H2I(i,j,k),tiny)
            s(i,iH2II) = max(H2II(i,j,k),tiny)
            s(i,iHM) = max(HM(i,j,k),tiny)
            s(i,ige) = max(ge(i,j,k),tiny)
            s(i,iDI) = max(DI(i,j,k),tiny)
            s(i,iDII) = max(DII(i,j,k),tiny)
            s(i,iHDI) = max(HDI(i,j,k),tiny)
            dtot(i) = fh*dtoh*d(i,j,k)
            rtot(i) = d(i,j,k)
            h2heatm(i) = tiny
            h2heatp(i) = tiny
            dmask(i) = itmask(i)
            correctH(i) = 1.0d0
            tsc(i) = 1.0d0
            tsgamma(i) = 1.0d0
#ifdef PYFORT
            tgasouti(i,j,k) = tgas(i)
#endif
         enddo
         do n= 1,NSPECIES
           do i = is+1, ie+1
             tscapy(i,n) = s(i,n)
           enddo
           do i = is+1, ie+1
             sm1(i,n) = s(i,n)
           enddo
           do i = is+1, ie+1
             sp(i,n) = s(i,n)
           enddo
         enddo
         

         doupdate = .false.
#include "primordial_derivstep.inc"

         do n = 1, ncspecies
           if(dtcoef(n).gt.0.0) then
             do i = is+1, ie+1
             if((itmask(i).eqv..true.).and.(eqmask(i).eqv..false.))then
               ! Figure out our dtit here
               kws = min(abs(dsq(i,n)-dsp(i,n)),huge)
               qt = max(s(i,n),tiny) / max(tiny, kws)
               temp1 = 0.1*abs(dtcoef(n)*qt)
               if(temp1.lt.dtit(i).and.temp1.gt.tiny)then
                 dtit(i) = temp1
                 mins(i) = n
               endif
#ifdef PYFORT
               tgasouti(i,j,k) = tgas(i)
c               dtit(i) = dt(i,j,k)/10.0
#endif
             endif
             enddo   ! end loop over i
           endif
         enddo
         do i = is+1, ie+1
         if ((itmask(i).eqv..true.))then!.and.(eqmask(i).eqv..false.))then
#ifdef PYFORT
           dtit(i) = min(dt(i,j,k)/2.0, abs(dt(i,j,k)-ttot(i)), dtit(i))
#else
           dtit(i) = min(dt/2.0, abs(dt-ttot(i)), dtit(i))
#endif
           olddtit(i) = dtit(i) ! First iteration!
         endif
         enddo

  100       continue


         coriter = 0

            do n = 1, NSPECIES
              do i = is+1, ie+1
              if(itmask(i).eqv..true.)then
                if(rejected(i).eqv..false.)then
                  sm1(i,n) = sp(i,n)
                  sp(i,n) = s(i,n)
                else
                  s(i,n) = sp(i,n)
                endif
              endif
              enddo
            enddo

            do i = is+1, ie+1
            if ((itmask(i).eqv..true.))then!.and.(eqmask(i).eqv..false.))then
              if(iter.eq.1e6)then
#ifdef PYFORT
                write(0,*) 't',i,j,k,dtit(i)/dt(i,j,k),mins(i)
#else
                write(0,*) 't',i,j,k,dtit(i)/dt,mins(i),ttot(i)/dt
#endif
                if(mins(i).gt.0)then
#ifdef PYFORT
                  write(0,*)'s',s(i,mins(i)),d(i,j,k),dt(i,j,k),dtit(i)
#else
                  write(0,*)'s',s(i,mins(i)),d(i,j,k),dt,dtit(i)
#endif
                  write(0,*)'q',dsp(i,mins(i)),dsq(i,mins(i)),
     & (dsq(i,mins(i))-dsp(i,mins(i)))/(dsq(i,mins(i))+dsp(i,mins(i)))
                  write(0,*)'p',tgas(i),gamma2(i),eqmask(i),
     & s(i,mins(i))/d(i,j,k)
                  write(0,*)'r',s(i,iH2I)/d(i,j,k),rejected(i),
     & correctH(i),correctHe(i)
                  write(0,*)'x',esterr(i)
                  do n = 1, ncspecies
                    write(0,*) 'aa',n,s(i,n)/d(i,j,k),dsq(i,n),dsp(i,n),
     &                        (dsq(i,n)-dsp(i,n))*dtit(i)/sp(i,n)
                  enddo
                endif
              endif
              rejected(i) = .false.
              if(iter.gt.1)then
                tsc(i) = olddtit(i)/dtit(i) ! c
                tsgamma(i) = (tsc(i)+1.0d0)/(tsc(i)+2.0d0) ! gamma
                do n = 1, ncspecies
                  temp2tsc = tsc(i)+tsc(i)
                  temptsc2 = tsc(i)*tsc(i)
                  temptsc12 = (tsc(i) + 1.0d0)*(tsc(i) + 1.0d0)
                  tscapy(i,n) = (temptsc12*sp(i,n)-sm1(i,n))
     &                        / (temptsc2 + temp2tsc)
                  if(tscapy(i,n).lt.0.0d0.and.
     &               capycare(n).eqv..true.)then
                    if(eqmask(i).eqv..true.)then!.and.eqable(n).eqv..true.)then
                    else
                    write(0,*) 'tscapy problems',eqmask(i),rejected(i),
     &                      iter,i,n,tsc(i),sp(i,n),sm1(i,n), esterr(i),
     &                      mins(i),eqmask(i),tscapy(i,n)
                    write(0,*) 'tscapy2', s(i,n)/d(i,j,k),
     &                          sp(i,n)/d(i,j,k),
     &                          sm1(i,n)/d(i,j,k),
     &                          s(i,n) * dom,
     &                          tgas(i)
                    write(0,*) 'tscapy3', correctH(i),
     &                          s(i,ige)/ge(i,j,k),
     &                          dsq(i,ige)/dsp(i,ige),
     &                          (dsq(i,ige)-dsp(i,ige))*dtit(i)
                    rejected(i) = .true.
c                    stop
                  endif
                  endif
                enddo
                if(itmask(i).eqv..false.)then
                  oreas = 3
#include "primordial_output.inc"
                endif
                
              endif
            endif
            enddo   ! end loop over i

            doupdate = .true.
            do coriter = 1, nciter
#include "primordial_derivstep.inc"
            enddo ! coriter
c         if(tgas(1).gt.9800)write(0,*)'aa',tgas(1),s(1,iH2I)/d(1,1,1),
c     & (dsq(1,ige)-dsp(1,ige))/(dsq(1,ige)+dsp(1,ige)),
c     & (dsq(1,iH2I)-dsp(1,iH2I))/(dsq(1,iH2I)+dsp(1,iH2I)),
c     & (h2heatp(1)-h2heatm(1))/(h2heatp(1)+h2heatm(1)), mins(1)

                

c         omaskflag = 0
 1002    format(6(I10,3X),9(e24.17,1x))
 1003    format(7(e24.17,1x))
 1004    format(1(I10,3X))
 1005    format(32(e24.17,1x))
c         write(0,*)j,k
         if(omaskflag.eq.1)then!.and.j.eq.1.and.k.eq.1)then
           do i=is+1,ie+1
c            i = 1
             if(itmask(i).eqv..true.) then
c             write(0,*)'writing output'
             oreas = 1
c             if(i.eq.1)write(0,*)tgas(i)
#include "primordial_output.inc"
             endif
           enddo
         endif
c         omaskflag = 0

            iter = iter + 1
            itflag = .false.
            do i=is+1,ie+1
            if((itmask(i).eqv..true.).and.(rejected(i).eqv..false.))then
              ttot(i) = ttot(i) + dtit(i)
#ifdef PYFORT
              if(dt(i,j,k).le.1.001*ttot(i))then
#else
              if(dt.le.1.001*ttot(i))then
#endif
                itmask(i) = .false.
c                write(0,*) d(i,j,k),s(i,iH2I)/d(i,j,k)
c                write(0,*) 'finishing',i,j,k,iter
             oreas = 0
#ifdef PYFORT
             if(subgridmask(i,j,k).eq.1)then
#include "primordial_output.inc"
             endif
                tgasoutf(i,j,k) = tgas(i)
#endif
              endif
            endif
            itflag = itmask(i) .or. itflag
            enddo

            if(iter.gt.itmax)then
              do i=is+1,ie+1
                if(itmask(i).eqv..true.)then
                  temp1 = ttot(i)/(
     &                  (1.0d0/sqrt(1.11e-31*dom*d(i,j,k)))/utim)
#ifndef PYFORT
c                  write(0,*)'FAILING?',i,j,k,mins(i),dtit(i),ttot(i),dt,olddtit(i),temp1
c                  write(0,*)'Tried for:',dt,'but only got',ttot(i)
c                  stop
#else
                  write(0,*)'FAIL',i,j,k,mins(i),dtit(i),ttot(i),
     &                          dt(i,j,k),olddtit(i),temp1
                  write(0,*)'Tried for:',dt(i,j,k),'but got',ttot(i)
                  if(mins(i).ne.-1)write(0,*)'did',s(i,mins(i))/d(i,j,k)
                  itmask(i) = .false.
                  oreas = 2
                  stop
#include "primordial_output.inc"
#endif
                endif
              enddo
c Uncomment this if you want to fail at all?
c              goto 9999
            endif

            if(itflag.eqv..false.)goto 9999

c           Guess stepsize for next iteration.
            do i = is+1, ie+1
            if(itmask(i).eqv..true.)then
              if(rejected(i).eqv..false.) olddtit(i) = dtit(i)
#ifdef PYFORT
              dtit(i) = min(dt(i,j,k)-ttot(i),max(
     &                  0.5,min(2.0,0.8/sqrt(esterr(i))))*dtit(i))
#else
              dtit(i) = min(dt-ttot(i),max(
     &                  0.5,min(2.0,0.8/sqrt(esterr(i))))*dtit(i))
#endif
            endif
            enddo

            goto 100

c     
 9999    continue

         do i=is+1, ie+1
           if(itmask(i).eqv..false.)then
           if(s(i,ige).le.0.0d0)write(6,*)'ge raspberried',i,j,k
           HI(i,j,k) = s(i,iHI)
           HII(i,j,k) = s(i,iHII)
           de(i,j,k) = s(i,ide)
           HeI(i,j,k) = s(i,iHeI)
           HeII(i,j,k) = s(i,iHeII)
           HeIII(i,j,k) = s(i,iHeIII)
           e(i,j,k) = e(i,j,k) + (s(i,ige)-ge(i,j,k))
           ge(i,j,k) = s(i,ige)
           H2I(i,j,k) = s(i,iH2I)
           H2II(i,j,k) = s(i,iH2II)
           HM(i,j,k) = s(i,iHM)
c           ge(i,j,k) = geprime(i) - (fh-s(i,iH2I)/d(i,j,k))*chunit
c           e(i,j,k) = ge(i,j,k) + 0.5*(u(i,j,k)**2.0 +
c     &                                 v(i,j,k)**2.0 +
c     &                                 w(i,j,k)**2.0)
           DI(i,j,k) = s(i,iDI)
           DII(i,j,k) = s(i,iDII)
           HDI(i,j,k) = s(i,iHDI)
           endif
         enddo
c     
c     Next j,k
c     

        do i=is+1,ie+1
            if(omask(i).eqv..true.)stop
        enddo
      totaliter = totaliter + iter
      enddo ! end of j
c        if(d(i,j,k)*dom.gt.1e9)then
c          write(6,*) 'finishing',k,totaliter/(ie-is)
c        endif
c        totaliter = 0
      enddo
c
c     Convert densities back to comoving from proper
c
      call primordial_scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, 
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, aye**3)
c
c     Correct the species to ensure consistency (i.e. type conservation)
c
      call primordial_make_consistent(
     &                     de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, d,
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, fh, dtoh, iter,
     &                     correctH, correctHe )


c      if (d(is,js,ks)*dom.gt.1e9.and.totaliter.gt.1e7)then
c        write(6,*) "solve_chemeq2: itcounts",
c     &        totaliter,DBLE(totaliter)/(jn*kn)
c      endif
      if (eqflag.eqv..true.) then
        write(6,*) "solve_chemeq2: total eqmask",totaleq,in*jn*kn
c       Now we check for NaN's
        do i=is+1,ie+1
        do j=js+1,je+1
        do k=ks+1,ke+1
          if(ge(i,j,k).ne.ge(i,j,k))write(6,*)'NAN! ge',i,j,k
          if(e(i,j,k).ne.e(i,j,k))write(6,*)'NAN! e',i,j,k
          if(HI(i,j,k).ne.HI(i,j,k))write(6,*)'NAN! HI',i,j,k
          if(HII(i,j,k).ne.HII(i,j,k))write(6,*)'NAN! HII',i,j,k
          if(HeI(i,j,k).ne.HeI(i,j,k))write(6,*)'NAN! HeI',i,j,k
          if(HeII(i,j,k).ne.HeII(i,j,k))write(6,*)'NAN! HeII',i,j,k
          if(HeIII(i,j,k).ne.HeIII(i,j,k))write(6,*)'NAN! HeIII',i,j,k
          if(H2I(i,j,k).ne.H2I(i,j,k))write(6,*)'NAN! H2I',i,j,k
          if(H2II(i,j,k).ne.H2II(i,j,k))write(6,*)'NAN! H2II',i,j,k
          if(HM(i,j,k).ne.HM(i,j,k))write(6,*)'NAN! HM',i,j,k
          if(de(i,j,k).ne.de(i,j,k))write(6,*)'NAN! de',i,j,k
          if(DI(i,j,k).ne.DI(i,j,k))write(6,*)'NAN! DI',i,j,k
          if(DII(i,j,k).ne.DII(i,j,k))write(6,*)'NAN! DII',i,j,k
          if(HDI(i,j,k).ne.HDI(i,j,k))write(6,*)'NAN! HDI',i,j,k
        enddo
        enddo
        enddo
      endif
      goto 9998

 9998 continue
c      if(omaskflag.eq.1)stop
      return
      end

c
c -----------------------------------------------------------
c   This routine scales the density fields from comoving to
c     proper densities (and back again).
c
      subroutine primordial_scale_fields(
     &                         d, de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, 
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, factor)
c -------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real     d(in,jn,kn)
      real    factor
c
c     locals
c
      integer i, j, k
c
c     Multiply all fields by factor (1/a^3 or a^3)
c
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
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*factor
                  H2I(i,j,k)  = H2I(i,j,k)*factor
                  H2II(i,j,k) = H2II(i,j,k)*factor
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*factor
                  DII(i,j,k) = DII(i,j,k)*factor
                  HDI(i,j,k) = HDI(i,j,k)*factor
               enddo
            endif
         enddo
      enddo
c
      return
      end


c ------------------------------------------------------------------
c   This routine correct the highest abundence species to
c     ensure conservation of particle number and charge.
c
      subroutine primordial_make_consistent(
     &                        de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, d,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, fh, dtoh, iter,
     &                        correctH, correctHe )
c -------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)  
      real    fh, dtoh, correctH(in), correctHe(in)
c
c     Parameters
c
      integer ijk, iter
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)
c
c     locals
c
      integer i, j, k
      real totalDen(ijk), totalH(ijk), totalHe(ijk), 
     &     totalD, correctD
c
c     Loop over all zones
c
      do k = ks+1, ke+1
      do j = js+1, je+1
c
c     Compute total densities of H and He
c         (ensure non-negativity)
c
      do i = is+1, ie+1
         HI   (i,j,k) = abs(HI   (i,j,k))
         HII  (i,j,k) = abs(HII  (i,j,k))
         HeI  (i,j,k) = abs(HeI  (i,j,k))
         HeII (i,j,k) = abs(HeII (i,j,k))
         HeIII(i,j,k) = abs(HeIII(i,j,k))
         totalH(i) = HI(i,j,k) + HII(i,j,k)
         totalHe(i) = HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k)
         totalDen(i) = totalH(i) + totalHe(i)
      enddo
c
c     include molecular hydrogen
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            HM   (i,j,k) = abs(HM   (i,j,k))
            H2II (i,j,k) = abs(H2II (i,j,k))
            H2I  (i,j,k) = abs(H2I  (i,j,k))
            totalH(i) = totalH(i) + HM(i,j,k) + H2I(i,j,k) + H2II(i,j,k)
            totalDen(i) = totalDen(i) + HM(i,j,k) + H2I(i,j,k)
     &                  + H2II(i,j,k)
         enddo
      endif
c
c     Correct densities by keeping fractions the same
c
      do i = is+1, ie+1
c         if (totalDen(i).ne.d(i,j,k))then
c            write(6,*) 'solve_rate_cool: totalDen0a', i,j,k
c            write(6,*) 'solve_rate_cool: totalDen0b', d(i,j,k), 
c     &                 totalDen
c         endif
         correctH(i) = fh*(d(i,j,k)/totalH(i))
c         if (abs(correctH(i)-1.0).gt.1e-1.and.d(i,j,k).gt.1e6) then
c            write(6,*) "correctH:",i,j,k,correctH(i),d(i,j,k),totalDen(i)
c         endif
         HI(i,j,k)  = HI(i,j,k)*correctH(i)
         HII(i,j,k) = HII(i,j,k)*correctH(i)
c
         correctHe(i) = (1.0 - fh)*(d(i,j,k)/totalHe(i))
         HeI(i,j,k)   = HeI(i,j,k)*correctHe(i)
         HeII(i,j,k)  = HeII(i,j,k)*correctHe(i)
         HeIII(i,j,k) = HeIII(i,j,k)*correctHe(i)
c
c     Correct molecular hydrogen-related fractions
c
         if (ispecies .gt. 1) then
            HM   (i,j,k) = HM(i,j,k)*correctH(i)
            H2II (i,j,k) = H2II(i,j,k)*correctH(i)
            H2I  (i,j,k) = H2I(i,j,k)*correctH(i)
         endif
      enddo
c
c     Do the same thing for deuterium (ignore HD) Assumes dtoh is small
c
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) + (2.0/3.0)*HDI(i,j,k)
            correctD = fh*dtoh*(d(i,j,k)/totalD)
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
c            HDI (i,j,k) = HDI(i,j,k)*(3.0/2.0)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
         enddo
      endif
c
c       Set the electron density and check total density
c
      do i = is+1, ie+1
         totalDen(i) = HI(i,j,k) + HII(i,j,k)
         totalDen(i) = totalDen(i) + HeI(i,j,k) + HeII(i,j,k)
     &               + HeIII(i,j,k)
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4. + HeIII(i,j,k)/2.
         if (ispecies .gt. 1) then 
            de(i,j,k) = de(i,j,k) - HM(i,j,k) + H2II(i,j,k)/2.
            totalDen(i) = totalDen(i) + HM(i,j,k) + H2I(i,j,k) 
     &                  + H2II(i,j,k)
         endif
         if (abs(log10(totalDen(i))-log10(d(i,j,k))).gt.0.01) then
            write(6,*) 'solve_rate_cool: totalDen0a', i,j,k
            write(6,*) 'solve_rate_cool: totalDen0b', d(i,j,k),
     &                 totalDen(i), iter
            write(6,*) 'solve_rate_cool: totalDen0c', HI(i,j,k),
     &       HII(i,j,k), H2I(i,j,k), H2II(i,j,k)
         endif
      enddo
c
      enddo  ! end loop over j
      enddo  ! end loop over k
c
      return
      end


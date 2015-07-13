#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c//////////////////////  SUBROUTINE COOL1D_SEP  \\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cool1d_sep(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, imetal, imcool, idim,
     &                is, ie, j, k, ih2co, ipiht, iter,
     &                aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, 
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow, 
     &                metala, metalc, n_xe, xe_start, xe_end, 
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot, nedot,
     &                tgas, tgasold, p2d,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, photogamma, itmask)
c
c  CALCULATE COOLING LUMINOSTIIES
c
c  written by: Yu Zhang, Peter Anninos and Tom Abel
c  date:       
c  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
c  modified2: October, 1996 by GB; moved to AMR
c  modified3: May, 2008 by JHW; separate edot into components
c
c  PURPOSE:
c    Solve the energy cooling equations.
c
c  INPUTS:
c    is,ie   - start and end indicies of active region (zero-based!)
c
c  PARAMETERS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c  Arguments
c
      INTG_PREC in, jn, kn, is, ie, j, k, nratec, imethod, idim,
     &        idual, iexpand, ih2co, ipiht, ispecies, imcool, 
     &        nfreq, iradshield, iradtype, imetalregen, 
     &        iradtrans, n_xe, nedot, imetal
      R_PREC    aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        eta1, eta2, gamma
      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn)
      R_PREC    photogamma(in,jn,kn)
      R_PREC    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec),
     $        metala(nratec, n_xe)
      R_PREC    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      R_PREC    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p, 
     &        xe_start, xe_end
c
c  Parameters
c
      real*8 mh
      R_PREC    ZSOLAR
      parameter (mh = mass_h, ZSOLAR = 0.02041_RKIND)
c
c  Locals
c
      INTG_PREC i, j1, iter, iradfield, tcmb_idx
      R_PREC dom, qq, vibl, logtem0, logtem9, dlogtem, energy, zr
      R_PREC dt2, ttmin, comp1, comp2, scoef, acoef, HIdot,
     &     hdlte1, hdlow1, gamma2, x, nH2, nother, fudge, fH2,
     &     gphdl1, factor, ciefudge, tau, dlogTcmb, logtcmb, metalc_cmb
      real*8 coolunit, dbase1, tbase1, xbase1, rtunits
c
c  Slice locals
c 
      INTG_PREC indixe(in)
      R_PREC t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &     tgas(in), tgasold(in)
      real*8 edot(in,nedot)
c
c  Cooling/heating slice locals
c
      real*8 ceHI(in), ceHeI(in), ceHeII(in),
     &     ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &     reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &     brem(in)
      R_PREC hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &     gpldl(in), gphdl(in), hdlte(in), hdlow(in), metalc(in)
c
c  Metal cooling locals
c
      INTG_PREC xe_idx
      R_PREC logxe0, logxe1, dlogxea, log_xe, xe1, dlogT, dlogxe, xi
      R_PREC xe_slope, xe_logtem0, xe_min, xe_max
      parameter (xe_slope = 5._RKIND, xe_logtem0 = 9.7859_RKIND)  ! e^9.78 = 10^4.25 K
      parameter (xe_min = -9.21_RKIND)  ! e^-9.21 = 1e-4
#ifdef CEN_METALS
      INTG_PREC nti, ndi, cmgenerate, NIT, NID, NIB
      R_PREC    TEMMIN, DELT, DENMIN, DELD, FREQDEL, FREQMIN
      PARAMETER(NIT=200,TEMMIN=3._RKIND,DELT=0.03_RKIND)
      PARAMETER(NID=300,DENMIN=-12._RKIND,DELD=0.05_RKIND)
      PARAMETER(NIB=400,FREQDEL=0.02_RKIND,FREQMIN=1._RKIND)
      R_PREC    metal_cool, metal_heat
      REAL*8    cbovcool(NIT,NID), cbovheat(NIT,NID), denwk(NID), 
     &        radt(NIB), eb(NIB)
      common  /cen_metal_com/ cbovcool, cbovheat, denwk
#endif /* CEN_METALS */

c     Iteration mask
      
      LOGIC_PREC itmask(in)

c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Set log values of start and end of lookup tables
c
      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem= (log(temend) - log(temstart))/REAL(nratec-1,RKIND)
c
c     Set units
c
      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      zr       = 1._RKIND/(aye*uaye) - 1._RKIND
      fudge    = 1._RKIND
      ciefudge = 1._RKIND
c
c     Set compton cooling coefficients (and temperature)
c
      if (iexpand .eq. 1) then
         comp1 = compa * (1._RKIND + zr)**4
         comp2 = 2.73_RKIND* (1._RKIND + zr)
      else
         comp1 = tiny
         comp2 = tiny
      endif
c
c     Compute Pressure
c
      if (imethod .eq. 2) then
c
c        Zeus - e() is really gas energy
c
         do i = is+1, ie+1
            if (itmask(i)) then
               p2d(i) = (gamma - 1._RKIND)*d(i,j,k)*e(i,j,k)
            endif
         enddo
      else
         if (idual .eq. 1) then
c
c           PPM with dual energy -- use gas energy
c
            do i = is+1, ie+1
               if (itmask(i)) then
                  p2d(i) = (gamma - 1._RKIND)*d(i,j,k)*ge(i,j,k)
               endif
            enddo
         else
c
c           PPM without dual energy -- use total energy
c
            do i = is+1, ie+1
               if (itmask(i)) then
                  p2d(i) = e(i,j,k) - 0.5_RKIND*u(i,j,k)**2
                  if (idim .gt. 1) 
     &                 p2d(i) = p2d(i) - 0.5_RKIND*v(i,j,k)**2
                  if (idim .gt. 2) 
     &                 p2d(i) = p2d(i) - 0.5_RKIND*w(i,j,k)**2
                  p2d(i) = max((gamma - 1._RKIND)*d(i,j,k)*p2d(i), 
     &                         tiny)
               endif
            enddo
         endif
      endif
c
c     Compute temperature
c
      do i = is+1, ie+1
         if (itmask(i)) then
            tgas(i) = 
     &           (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4._RKIND +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
         endif
      enddo
c
c          (include molecular hydrogen, but ignore deuterium)
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
               tgas(i) = tgas(i) +
     &              HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2._RKIND
            endif
         enddo
      endif
c
      do i = is+1, ie+1
         if (itmask(i)) then
            tgas(i) = max(p2d(i)*utem/tgas(i), temstart)
         endif
      enddo
c
c     Correct temperature for gamma from H2
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
               nH2 = 0.5_RKIND*(H2I(i,j,k) + H2II(i,j,k))
               nother = (HeI(i,j,k)+HeII(i,j,k)+HeIII(i,j,k))/4._RKIND
     &              + HI(i,j,k) + HII(i,j,k) + de(i,j,k)
               if (nH2/nother .gt. 1.0e-3_RKIND) then
                  x = 6100/tgas(i) ! not quite self-consistent
                  if (x .gt. 10._RKIND) then
                     gamma2 = 0.5_RKIND*5._RKIND
                  else
                    gamma2 = 0.5_RKIND*(5._RKIND + 2._RKIND*x**2
     &                    * exp(x)/(exp(x)-1)**2)
                  endif
               else
                  gamma2 = 2.5_RKIND
               endif
               gamma2 = 1._RKIND + (nH2 + nother)/
     &              (nH2*gamma2 + nother/(gamma-1._RKIND))
               tgas(i) = tgas(i) * (gamma2-1._RKIND)/(gamma-1._RKIND)
            endif
         enddo
      endif
c
c     If this is the first time through, just set tgasold to tgas
c
      if (iter .eq. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then 
               tgasold(i) = tgas(i)
            endif
         enddo
      endif
c
c     --- 6 species cooling ---
c
      do i = is+1, ie+1

         if (itmask(i)) then

c
c        Compute log temperature and truncate if above/below table max/min
c
         logtem(i) = log(0.5_RKIND*(tgas(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)
c
c        Compute index into the table and precompute parts of linear interp
c
         indixe(i) = min(nratec-1, max(1,
     &        int((logtem(i)-logtem0)/dlogtem,IKIND)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = 1._RKIND/(t2(i) - t1(i))
c
c        Lookup cooling values and do a linear temperature in log(T)
c
         ceHI(i) = ceHIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ceHIa(indixe(i)+1) -ceHIa(indixe(i)))*tdef(i)
         ceHeI(i) = ceHeIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ceHeIa(indixe(i)+1) -ceHeIa(indixe(i)))*tdef(i)
         ceHeII(i) = ceHeIIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ceHeIIa(indixe(i)+1) -ceHeIIa(indixe(i)))*tdef(i)
         ciHI(i) = ciHIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHIa(indixe(i)+1) -ciHIa(indixe(i)))*tdef(i)
         ciHeI(i) = ciHeIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHeIa(indixe(i)+1) -ciHeIa(indixe(i)))*tdef(i)
         ciHeIS(i) = ciHeISa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHeISa(indixe(i)+1) -ciHeISa(indixe(i)))*tdef(i)
         ciHeII(i) = ciHeIIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHeIIa(indixe(i)+1) -ciHeIIa(indixe(i)))*tdef(i)
         reHII(i) = reHIIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(reHIIa(indixe(i)+1) -reHIIa(indixe(i)))*tdef(i)
         reHeII1(i)=reHeII1a(indixe(i)) + (logtem(i) - t1(i))
     .        *(reHeII1a(indixe(i)+1)-reHeII1a(indixe(i)))*tdef(i)
         reHeII2(i)=reHeII2a(indixe(i)) + (logtem(i) - t1(i))
     .        *(reHeII2a(indixe(i)+1)-reHeII2a(indixe(i)))*tdef(i)
         reHeIII(i)=reHeIIIa(indixe(i)) + (logtem(i) - t1(i))
     .        *(reHeIIIa(indixe(i)+1)-reHeIIIa(indixe(i)))*tdef(i)
         brem(i) = brema(indixe(i)) + (logtem(i) - t1(i))
     .         *(brema(indixe(i)+1) -brema(indixe(i)))*tdef(i)

      endif
      enddo
c
c     Compute the cooling function
c
      do i = is+1, ie+1
         if (itmask(i)) then
c
c                    Collisional excitations
c
            edot(i,1)=-ceHI(i)*HI  (i,j,k)*de(i,j,k) ! ce of HI
            edot(i,2)=-ceHeI(i)*HeII(i,j,k)*de(i,j,k)**2*dom*0.25_RKIND ! ce of HeI
            edot(i,3)=-ceHeII(i)*HeII(i,j,k)*de(i,j,k)*0.25_RKIND ! ce of HeII
c
c                    Collisional ionizations
c
            edot(i,4)=-ciHI  (i)*HI  (i,j,k)*de(i,j,k) ! ci of HI
            edot(i,5)=-ciHeI (i)*HeI (i,j,k)*de(i,j,k)*0.25_RKIND ! ci of HeI
            edot(i,6)=-ciHeII(i)*HeII(i,j,k)*de(i,j,k)*0.25_RKIND ! ci of HeII
            edot(i,7)=-ciHeIS(i)*HeII(i,j,k)*de(i,j,k)**2*dom*0.25_RKIND ! ci of HeIS
c
c                    Recombinations
c
            edot(i,8)  = -reHII  (i)*HII  (i,j,k)*de(i,j,k) ! re of HII
            edot(i,9)  = -reHeII1(i)*HeII (i,j,k)*de(i,j,k)*0.25_RKIND ! re of HeII
            edot(i,10) = -reHeII2(i)*HeII (i,j,k)*de(i,j,k)*0.25_RKIND ! re of HeII
            edot(i,11) = -reHeIII(i)*HeIII(i,j,k)*de(i,j,k)*0.25_RKIND ! re of HeIII
c
c                    Compton cooling or heating
c
            edot(i,12) = -comp1*(tgas(i)-comp2)*de(i,j,k)/dom
c
c                    X-ray compton heating
c
            edot(i,13) = -comp_xraya * (tgas(i)-comp_temp)*de(i,j,k)/dom
c
c                    Bremsstrahlung
c
            edot(i,14) = -brem(i)*(HII(i,j,k)+HeII(i,j,k)*0.25_RKIND+
     $           HeIII(i,j,k))*de(i,j,k)

         endif
      enddo
c     
c     --- H2 cooling ---
c
      if (ispecies .gt. 1) then
c
#define USE_GALLI_PALLA1999
#define OPTICAL_DEPTH_FUDGE
c
c        Use the Galli and Palla (1999) cooling rates for molecular H.
c
#ifdef USE_GALLI_PALLA1999
c
         do i = is+1, ie+1
            if (itmask(i)) then
            gpldl(i) = gpldla(indixe(i)) + (logtem(i) - t1(i))
     .         *(gpldla(indixe(i)+1) - gpldla(indixe(i)))*tdef(i)
            gphdl(i) = gphdla(indixe(i)) + (logtem(i) - t1(i))
     .         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))*tdef(i)
c            cieco(i) = ciecoa(indixe(i)) + (logtem(i) - t1(i))
c     .         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))*tdef(i)
            endif
         enddo
c
         do i = is+1, ie+1

            if (itmask(i)) then

c
#ifdef OPTICAL_DEPTH_FUDGE
c            nH2 = 0.5*H2I(i,j,k)
c            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
c     .                HI(i,j,k) + HII(i,j,k) + de(i,j,k)
c            fH2 = nH2/(nH2 + nother)
c            fudge = sqrt((40.0 * 10**(4.8 * 
c     .               sqrt(max(log10(tgas(i)),2.0)-2.0)) / fH2**2)/
c     .               ((nH2 + nother)*dom) )
            fudge = (0.76_RKIND*d(i,j,k)*dom/8.0e9_RKIND)**(-0.45_RKIND)
            fudge = min(fudge, 1._RKIND)
#endif /* OPTICAL_DEPTH_FUDGE */
c
            gphdl1 = gphdl(i)/(HI(i,j,k)*dom)
            edot(i,15) = -REAL(ih2co,RKIND)*fudge*H2I(i,j,k)*
     .           gphdl(i)/(1._RKIND + gphdl1/gpldl(i)) / (2._RKIND*dom)
            endif
         enddo
c
#else /* USE_GALLI_PALLA1999 */
c
         do i = is+1, ie+1
            if (itmask(i)) then
            hyd01k(i) = hyd01ka(indixe(i)) + (logtem(i) - t1(i))
     .         *(hyd01ka(indixe(i)+1)-hyd01ka(indixe(i)))*tdef(i)
            h2k01(i) = h2k01a(indixe(i)) + (logtem(i) - t1(i))
     .         *(h2k01a(indixe(i)+1) - h2k01a(indixe(i)))*tdef(i)
            vibh(i) = vibha(indixe(i)) + (logtem(i) - t1(i))
     .         *(vibha(indixe(i)+1) - vibha(indixe(i)))*tdef(i)
            roth(i) = rotha(indixe(i)) + (logtem(i) - t1(i))
     .         *(rotha(indixe(i)+1) - rotha(indixe(i)))*tdef(i)
            rotl(i) = rotla(indixe(i)) + (logtem(i) - t1(i))
     .         *(rotla(indixe(i)+1) - rotla(indixe(i)))*tdef(i)
c            cieco(i) = ciecoa(indixe(i)) + (logtem(i) - t1(i))
c     .         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))*tdef(i)
            endif
         enddo
c
         do i = is+1, ie+1
            if (itmask(i)) then
            qq   = 1.2_RKIND*(HI(i,j,k)*dom)**0.77_RKIND + 
     .                (H2I(i,j,k)*dom/2._RKIND)**0.77_RKIND
            vibl = (HI(i,j,k)*hyd01k(i) + 
     .             H2I(i,j,k)/2._RKIND*h2k01(i))
     .             *dom*8.18e-13_RKIND
c
#ifdef OPTICAL_DEPTH_FUDGE
c            nH2 = 0.5*H2I(i,j,k)
c            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
c     .                HI(i,j,k) + HII(i,j,k) + de(i,j,k)
c            fH2 = nH2/(nH2 + nother)
c            fudge = sqrt((40.0 * 10**(4.8 * 
c     .               sqrt(max(log10(tgas(i)),2.0)-2.0)) / fH2**2)/
c     .               ((nH2 + nother)*dom) )
            fudge = (0.76_RKIND*d(i,j,k)*dom/8.e9_RKIND)**(-0.45_RKIND)
            fudge = min(fudge, 1._RKIND)
#endif /* OPTICAL_DEPTH_FUDGE */
c
            edot(i,15) = -REAL(ih2co,RKIND)*fudge*H2I(i,j,k)*(
     .           vibh(i)/(1._RKIND+vibh(i)/max(   vibl,tiny)) +
     .           roth(i)/(1._RKIND+roth(i)/max(qq*rotl(i),tiny))     
     .           )/2._RKIND/dom
            endif
         enddo
c
#endif /* USE_GALLI_PALLA1999 */
c

         do i = is+1, ie+1
#define NO_CIE_COOLING
#ifdef CIE_COOLING
c     cooling from H2-H2 and He-H2 collisional induced emission comes
c     with its own radiative transfer correction as discussed in
C     Ripamonti & Abel 2003
c            ciefudge = 1.
c            tau = (d(i,j,k)*dom/2e16)**2.8 ! 2e16 is in units of cm^-3
c            tau = max(tau, 1.e-5)
c            ciefudge = (1.-exp(-tau))/tau
c            edot(i,16) = -H2I(i,k,k)*d(i,j,k)*cieco(i)*ciefudge
#else /* CIE_COOLING */
            edot(i,16) = 0._RKIND
#endif
         enddo

      endif                     !  ispecies = 2
c
c     --- Cooling from HD ---
c
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
            hdlte(i) = hdltea(indixe(i)) + (logtem(i) - t1(i))
     .         *(hdltea(indixe(i)+1) - hdltea(indixe(i)))*tdef(i)
            hdlow(i) = hdlowa(indixe(i)) + (logtem(i) - t1(i))
     .         *(hdlowa(indixe(i)+1) - hdlowa(indixe(i)))*tdef(i)
            endif
         enddo
c
         do i = is+1, ie+1
            if (itmask(i)) then
               hdlte1 = hdlte(i)/(HDI(i,j,k)*dom/2._RKIND)
               hdlow1 = max(hdlow(i), tiny)
               edot(i,17) = -HDI(i,j,k)*
     .              (hdlte1/(1._RKIND + hdlte1/hdlow1)/(2._RKIND*dom))
            endif
         enddo
      else
         do i = is+1, ie+1
            edot(i,17) = 0._RKIND
         enddo
      endif
c
c     --- Cooling/heating due to metals ---
c

      if (imcool .eq. 1) then
         logxe0 = log(xe_start)
         logxe1 = log(xe_end)
         dlogxea = (logxe1 - logxe0) / REAL(n_xe-1,RKIND)
c
c  Determine temperature for CMB temperature so we can subtract metal
c  cooling at T=T_cmb
c
         logtcmb = log(comp2)
         tcmb_idx = max(1, int((logtcmb-logtem0)/dlogtem,IKIND)+1)
         dlogTcmb = (logtcmb - (logtem0 + (tcmb_idx-1)*dlogtem)) 
     $	/ dlogtem
         do i = is+1, ie+1
            if (itmask(i)) then
               log_xe = max(min(log(de(i,j,k) / d(i,j,k)), logxe1), 
     $              logxe0)
               if (logtem(i) .lt. xe_logtem0) then
                  xe_max = max(xe_slope * (logtem(i) - xe_logtem0), 
     $                         xe_min)
                  log_xe = min(log_xe, xe_max)
               endif
               xe_idx = min(n_xe-1, max(1, 
     $              int((log_xe - logxe0)/dlogxea,IKIND)+1))
               xe1 = logxe0 + (xe_idx - 1)*dlogxea
               dlogT = (logtem(i) - t1(i)) * tdef(i)
               dlogxe = (log_xe - xe1) / dlogxea
c
c  Interpolate metal cooling in temperature and electron fraction
c
               metalc(i) = 
     $              metala(indixe(i),   xe_idx  )
     $                 * (1._RKIND-dlogxe)*(1._RKIND-dlogT)+
     $              metala(indixe(i)+1, xe_idx  )
     $                 * (dlogxe  )*(1._RKIND-dlogT)   +
     $              metala(indixe(i),   xe_idx+1)
     $                 * (1._RKIND-dlogxe)*(dlogT)    +
     $              metala(indixe(i)+1, xe_idx+1)
     $                 * (dlogxe  )*(dlogT)
c
c  Subtract metal cooling at T_cmb
c
               metalc_cmb = 
     $              metala(tcmb_idx,   xe_idx  )
     $                 * (1._RKIND-dlogxe)*(1._RKIND-dlogTcmb)+
     $              metala(tcmb_idx+1, xe_idx  )
     $                 * (dlogxe  ) * (1._RKIND-dlogTcmb) +
     $              metala(tcmb_idx,   xe_idx+1)
     $                 * (1._RKIND-dlogxe) * (dlogTcmb)   +
     $              metala(tcmb_idx+1, xe_idx+1)
     $                 * (dlogxe  ) * (dlogTcmb)
               metalc(i) = metalc(i) - metalc_cmb

               xi = min(3._RKIND, metal(i,j,k)/(d(i,j,k)*ZSOLAR))
               edot(i,18) = -metalc(i) * d(i,j,k) * d(i,j,k) *
     $              (xi/JHW_METALS_NORMZ)
#ifdef UNUSED
               if (xi .gt. 1e-06) then
                  write(6,777) i,j,k, exp(log_xe), exp(logtem(i)), 
     $                 exp(logTcmb), indixe(i), tcmb_idx
                  write(6,778) d(i,j,k), xi, dlogxe, dlogT
                  write(6,779) xe_idx, metala(indixe(i), xe_idx),
     $                 metalc(i), metalc_cmb
                  write(6,780) edot(i,18),
     $                 metalc(i)*d(i,j,k)*(xi/JHW_METALS_NORMZ)/dom, 
     $                 coolunit
                  write(6,*) ''
               endif
 777           format('a0: ', 3(i5), 5x, 3(es14.3), 2(i5))
 778           format('a1: ', 4(es14.3))
 779           format('a2: ', i5, 3(es14.3))
 780           format('a3: ', 3(es14.3))
#endif

            endif
         enddo
      endif ! imcool == 1

#ifdef CEN_METALS
c
      if (imcool .eq. 2) then
c
c     Generate table if required
c
      if (imetalregen .eq. 1) then
         write(6,*) 'generating metallicity cooling table'
         if (iradtype .eq. 11 .or. iradtype .eq. 12) iradfield = 1
c
c        Clear table
c
         do j1=1, NID
            denwk(j1) = 10._RKIND**(DENMIN+(j1-0.5_RKIND)*DELD)
            do i=1, NIT
               cbovcool(i,j1) = 0.d0
               cbovheat(i,j1) = 0.d0
            enddo
         enddo
c
c        compute energy (in eV) of each bin
c
         if (iradfield .eq .1) then
            if (NIB .ne. nfreq) then
               write(6,*) 'cool1d_multi: NIB != nfreq'
               stop
            endif
c
            do i=1, NIB
               eb(i) = FREQMIN*10._RKIND**(FREQDEL*(i-1._RKIND))
            enddo
c
c        Convert to 10^{-21} erg/cm^2/s
c
            do i=1, NIB-1
               radt(i) = inutot(i)*4._RKIND*pi_val
     &                   *(eb(i+1)-eb(i))*(1.2e11_RKIND/6.625_RKIND)
            enddo
         else
            do i=1, NIB-1
               radt(i) = 0._RKIND
            enddo
         endif
c
         radt(NIB) = 0
         call mtlclht(denwk, radt, iradfield, cbovcool, cbovheat)
         imetalregen = 0
c
c        Convert from 10^-23 to our units
c
         do j1=1, NID
            do i=1, NIT
               cbovcool(i,j1) = cbovcool(i,j1) * (1.d-23/coolunit)
               cbovheat(i,j1) = cbovheat(i,j1) * (1.d-23/coolunit)
c               write(31,*) i,j1,cbovcool(i,j1),cbovheat(i,j1)
            enddo
         enddo
c
      endif
c
c     Look-up in table
c
      do i = is+1, ie+1
         if (itmask(i)) then
         ndi = (log10(de(i,j,k)*dom)-DENMIN)/DELD + 1
         ndi = min(max(1, ndi), NID)
         nti = (log10(0.5_RKIND*(tgas(i)+tgasold(i)))-TEMMIN)/DELT + 1
         nti = min(max(1, nti), NIT)
         xi = min(3._RKIND, metal(i,j,k)/(d(i,j,k)*ZSOLAR))
         metal_cool = cbovcool(nti,ndi) * xi / dom *
     &                (HI(i,j,k)+HII(i,j,k))
         metal_heat = cbovheat(nti,ndi) * xi / dom *
     &                (HI(i,j,k)+HII(i,j,k))
         if (tgas(i) .le. 1e4) then
            metal_cool = metal_cool*(tgas(i)/1e4)**3
            metal_heat = 0._RKIND
         endif
         edot(i,18) = metal_heat - metal_cool
         endif
      enddo
c
      endif
c
#endif /* CEN_METALS */
c
c     Set tgasold
c
      do i=is+1, ie+1
         if (itmask(i)) then
            tgasold(i) = tgas(i)
         endif
      enddo
c
      return
      end


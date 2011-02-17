#include "fortran.def"
#include "error.def"

!=======================================================================
!//////////////////////  SUBROUTINE COOL1D_MULTI  \\\\\\\\\\\\\\\\\\\\\\

      subroutine cool1d_multi(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, imetal, imcool, idim, 
     &                is, ie, j, k, ih2co, ipiht, iter, igammah,
     &                aye, temstart, temend,
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
     &                metala, n_xe, xe_start, xe_end, 
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, p2d,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, photogamma,
     &                ih2optical, iciecool, ciecoa, cieco,
     &                icmbTfloor, iClHeat, iClMMW,
     &                clMetNorm, clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating, clMMW,
     &                itmask)

!  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
!  modified2: October, 1996 by GB; moved to AMR
!  modified3: February, 2003 by Robert Harkness; iteration mask
!  modified6: September, 2009 by BDS to include cloudy cooling
!
!  PURPOSE:
!    Solve the energy cooling equations.
!
!  INPUTS:
!    is,ie   - start and end indicies of active region (zero-based!)
!
!  PARAMETERS:
!

!-----------------------------------------------------------------------

      implicit NONE

!  Arguments

      integer in, jn, kn, is, ie, j, k, nratec, imethod, idim,
     &        idual, iexpand, ih2co, ipiht, ispecies, imcool, 
     &        nfreq, iradshield, iradtype, imetalregen, iradtrans,
     &        n_xe, imetal, igammah, ih2optical, iciecool
      real    aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        eta1, eta2, gamma
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn)
      real    photogamma(in,jn,kn)
      real    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec),
     &        metala(nratec, n_xe)
      real    gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec), ciecoa(nratec)
      real    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      real    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p
      real    gammaha, xe_start, xe_end

!  Cloudy cooling data

      integer icmbTfloor, iClHeat, iClMMW, clGridRank, clDataSize
      integer clGridDim(clGridRank)
      real clMetNorm, clEleFra
      real clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)),
     &     clPar5(clGridDim(5))
      real clCooling(clDataSize), clHeating(clDataSize), 
     &     clMMW(clDataSize)

!  Parameters

      double precision mh
      real    ZSOLAR, mu_metal
      parameter (mh = 1.67d-24)      !DPC
      parameter (ZSOLAR = 0.02041)
      parameter (mu_metal = 16.0)    ! approx. mean molecular weight of metals

!  Locals

      integer i, j1, iter, iradfield
      real dom, qq, vibl, logtem0, logtem9, dlogtem, energy, zr
      real dt2, ttmin, comp1, comp2, scoef, acoef, HIdot,
     &     hdlte1, hdlow1, gamma2, x, nH2, nother, fudge, fH2,
     &     gphdl1, factor, dom_inv, tau, ciefudge
      double precision coolunit, dbase1, tbase1, xbase1, rtunits

!  Slice locals
 
      integer indixe(in)
      real t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &     tgas(in), tgasold(in)
      double precision edot(in)

!  Cooling/heating slice locals

      double precision ceHI(in), ceHeI(in), ceHeII(in),
     &     ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &     reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &     brem(in), cieco(in)
      real hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &     gpldl(in), gphdl(in), hdlte(in), hdlow(in), metalc(in)
      real gaHI(in), gaH2(in), gaHe(in), gaHp(in), gael(in),
     &     galdl(in)

!  Metal cooling locals

      integer xe_idx, tcmb_idx
      real logxe0, logxe1, dlogxea, xe, log_xe, xe1, dlogT, dlogxe, xi
      real logTcmb, dlogTcmb, metalc_cmb, xe_max, metalfree
      real xe_slope, xe_logtem0, xe_min
      parameter (xe_slope = 5.0, xe_logtem0 = 9.7859)  ! e^9.78 = 10^4.25 K
      parameter (xe_min = -9.21)  ! e^-9.21 = 1e-4
#ifdef CEN_METALS
      integer nti, ndi, cmgenerate, NIT, NID, NIB
      real    TEMMIN, DELT, DENMIN, DELD, FREQDEL, FREQMIN
      PARAMETER(NIT=200,TEMMIN=3.0,DELT=0.03)
      PARAMETER(NID=300,DENMIN=-12.0,DELD=0.05)
      PARAMETER(NIB=400,FREQDEL=0.02,FREQMIN=1.0)
      real    metal_cool, metal_heat, xi
      real    cbovcool(NIT,NID), cbovheat(NIT,NID), denwk(NID), 
     &        radt(NIB), eb(NIB)
      common  /cen_metal_com/ cbovcool, cbovheat, denwk
#endif /* CEN_METALS */

!  Iteration mask

      logical itmask(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set log values of start and end of lookup tables

      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem= (log(temend) - log(temstart))/real(nratec-1)

!     Set units

      dom      = urho*(aye**3)/mh
      dom_inv  = 1.0/dom
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      zr       = 1.0/(aye*uaye) - 1.0
      fudge    = 1.0
      iradfield = -1

!     Set compton cooling coefficients (and temperature)

      if (iexpand .eq. 1) then
         comp1 = compa * (1.0 + zr)**4
         comp2 = 2.73  * (1.0 + zr)
      else
         comp1 = tiny
         comp2 = tiny
      endif

!     Compute Pressure

      if (imethod .eq. 2) then

!        Zeus - e() is really gas energy

         do i = is+1, ie+1
            if ( itmask(i) ) then
            p2d(i) = (gamma - 1.0)*d(i,j,k)*e(i,j,k)
            end if
         enddo
      else
         if (idual .eq. 1) then

!           PPM with dual energy -- use gas energy

            do i = is+1, ie+1
               if ( itmask(i) ) then
               p2d(i) = (gamma - 1.0)*d(i,j,k)*ge(i,j,k)
               end if
            enddo
         else

!           PPM without dual energy -- use total energy

            do i = is+1, ie+1
               if ( itmask(i) ) then
               p2d(i) = e(i,j,k) - 0.5*u(i,j,k)**2
               if (idim .gt. 1) p2d(i) = p2d(i) - 0.5*v(i,j,k)**2
               if (idim .gt. 2) p2d(i) = p2d(i) - 0.5*w(i,j,k)**2
               p2d(i) = max((gamma - 1.0)*d(i,j,k)*p2d(i), tiny)
               end if
            enddo
         endif
      endif

!     Compute temperature

      do i = is+1, ie+1
         if ( itmask(i) ) then
         tgas(i) = 
     &           (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
     &            HI(i,j,k) + HII(i,j,k) + de(i,j,k)
         end if
      enddo

!          (include molecular hydrogen, but ignore deuterium)

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
            tgas(i) = tgas(i) +
     &            HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2.0
            end if
         enddo
      endif

!     Include metal species
      
      if (imetal .eq. 1) then
         do i = is+1, ie+1
            tgas(i) = tgas(i) + metal(i,j,k)/mu_metal
         enddo
      endif

      do i = is+1, ie+1
         if ( itmask(i) ) then
         tgas(i) = max(p2d(i)*utem/tgas(i), temstart)
         end if
      enddo

!     Correct temperature for gamma from H2

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
            nH2 = 0.5*(H2I(i,j,k) + H2II(i,j,k))
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            if (nH2/nother .gt. 1.0e-3) then
               x = 6100/tgas(i) ! not quite self-consistent
               if (x .gt. 10.0) then
                  gamma2 = 0.5*5.0
               else
                  gamma2 = 0.5*(5.0 + 2.0*x**2 * exp(x)/(exp(x)-1)**2)
               endif
            else
               gamma2 = 2.5
            endif
            gamma2 = 1.0 + (nH2 + nother)/
     &                     (nH2*gamma2 + nother/(gamma-1.0))
            tgas(i) = tgas(i) * (gamma2 - 1.0)/(gamma - 1.0)
            end if
         enddo
      endif

!     If this is the first time through, just set tgasold to tgas

      if (iter .eq. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
            tgasold(i) = tgas(i)
            end if
         enddo
      endif

!     --- 6 species cooling ---

      do i = is+1, ie+1
         if ( itmask(i) ) then

!        Compute log temperature and truncate if above/below table max/min

         logtem(i) = log(0.5*(tgas(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

!        Compute index into the table and precompute parts of linear interp

         indixe(i) = min(nratec-1,
     &                   max(1,int((logtem(i)-logtem0)/dlogtem)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

!        Lookup cooling values and do a linear temperature in log(T)

         ceHI(i) = ceHIa(indixe(i)) + tdef(i)
     &         *(ceHIa(indixe(i)+1) -ceHIa(indixe(i)))
         ceHeI(i) = ceHeIa(indixe(i)) + tdef(i)
     &         *(ceHeIa(indixe(i)+1) -ceHeIa(indixe(i)))
         ceHeII(i) = ceHeIIa(indixe(i)) + tdef(i)
     &         *(ceHeIIa(indixe(i)+1) -ceHeIIa(indixe(i)))
         ciHI(i) = ciHIa(indixe(i)) + tdef(i)
     &         *(ciHIa(indixe(i)+1) -ciHIa(indixe(i)))
         ciHeI(i) = ciHeIa(indixe(i)) + tdef(i)
     &         *(ciHeIa(indixe(i)+1) -ciHeIa(indixe(i)))
         ciHeIS(i) = ciHeISa(indixe(i)) + tdef(i)
     &         *(ciHeISa(indixe(i)+1) -ciHeISa(indixe(i)))
         ciHeII(i) = ciHeIIa(indixe(i)) + tdef(i)
     &         *(ciHeIIa(indixe(i)+1) -ciHeIIa(indixe(i)))
         reHII(i) = reHIIa(indixe(i)) + tdef(i)
     &         *(reHIIa(indixe(i)+1) -reHIIa(indixe(i)))
         reHeII1(i)=reHeII1a(indixe(i)) + tdef(i)
     &        *(reHeII1a(indixe(i)+1)-reHeII1a(indixe(i)))
         reHeII2(i)=reHeII2a(indixe(i)) + tdef(i)
     &        *(reHeII2a(indixe(i)+1)-reHeII2a(indixe(i)))
         reHeIII(i)=reHeIIIa(indixe(i)) + tdef(i)
     &        *(reHeIIIa(indixe(i)+1)-reHeIIIa(indixe(i)))
         brem(i) = brema(indixe(i)) + tdef(i)
     &         *(brema(indixe(i)+1) -brema(indixe(i)))

         end if
      enddo

!     Compute the cooling function

      do i = is+1, ie+1
         if ( itmask(i) ) then
         edot(i) = (

!                    Collisional excitations

     &             - ceHI  (i)*HI  (i,j,k)*de(i,j,k)             ! ce of HI
     &             - ceHeI (i)*HeII(i,j,k)*de(i,j,k)**2*dom/4.0  ! ce of HeI
     &             - ceHeII(i)*HeII(i,j,k)*de(i,j,k)/4.0         ! ce of HeII

!                    Collisional ionizations

     &             - ciHI  (i)*HI  (i,j,k)*de(i,j,k)             ! ci of HI
     &             - ciHeI (i)*HeI (i,j,k)*de(i,j,k)/4.0         ! ci of HeI
     &             - ciHeII(i)*HeII(i,j,k)*de(i,j,k)/4.0         ! ci of HeII
     &             - ciHeIS(i)*HeII(i,j,k)*de(i,j,k)**2*dom/4.0  ! ci of HeIS

!                    Recombinations

     &             - reHII  (i)*HII  (i,j,k)*de(i,j,k)          ! re of HII
     &             - reHeII1(i)*HeII (i,j,k)*de(i,j,k)/4.0      ! re of HeII
     &             - reHeII2(i)*HeII (i,j,k)*de(i,j,k)/4.0      ! re of HeII
     &             - reHeIII(i)*HeIII(i,j,k)*de(i,j,k)/4.0      ! re of HeIII

!                    Compton cooling or heating

     &             - comp1*(tgas(i)-comp2)*de(i,j,k)*dom_inv

!                    X-ray compton heating

     &             - comp_xraya*(tgas(i)-comp_temp)*de(i,j,k)*dom_inv

!                    Bremsstrahlung

     &             - brem(i)*(HII(i,j,k)+HeII(i,j,k)/4.0+HeIII(i,j,k))
     &                        *de(i,j,k)

!                    Photoelectric heating by UV-irradiated dust

     &             + float(igammah)*gammaha*(HI(i,j,k)+HII(i,j,k))
     &             *dom_inv)
        
         if (edot(i) .ne. edot(i)) then
            write(6,*) 'NaN in edot[1]: ', i, j, k, edot(i), 
     &           HI(i,j,k), HII(i,j,k), HeI(i,j,k), HeII(i,j,k), 
     &           HeIII(i,j,k), de(i,j,k), d(i,j,k), 
     &           tgas(i), p2d(i)
            ERROR_MESSAGE
         endif
         
         end if
      enddo
     
!     --- H2 cooling ---

      if (ispecies .gt. 1) then

#define USE_GLOVER_ABEL2008
#ifdef USE_GLOVER_ABEL2008
         do i = is+1, ie+1
            if ( itmask(i) ) then
            gaHI(i) = gaHIa(indixe(i)) + tdef(i)
     &         *(gaHIa(indixe(i)+1) - gaHIa(indixe(i)))
            gaH2(i) = gaH2a(indixe(i)) + tdef(i)
     &         *(gaH2a(indixe(i)+1) - gaH2a(indixe(i)))
            gaHe(i) = gaHea(indixe(i)) + tdef(i)
     &         *(gaHea(indixe(i)+1) - gaHea(indixe(i)))
            gaHp(i) = gaHpa(indixe(i)) + tdef(i)
     &         *(gaHpa(indixe(i)+1) - gaHpa(indixe(i)))
            gael(i) = gaela(indixe(i)) + tdef(i)
     &         *(gaela(indixe(i)+1) - gaHIa(indixe(i)))
            gphdl(i) = gphdla(indixe(i)) + tdef(i)
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then
#ifdef OPTICAL_DEPTH_FUDGE
            nH2 = 0.5*H2I(i,j,k)
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
     &                HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            fH2 = nH2/(nH2 + nother)
            fudge = sqrt((40.0 * 10**(4.8 * 
     &               sqrt(max(log10(tgas(i)),2.0)-2.0)) / fH2**2)/
     &               ((nH2 + nother)*dom) )
            fudge = min(fudge, 1.0)
#endif /* OPTICAL_DEPTH_FUDGE */
            ! Note that this optical depth approximation comes from
            ! RA04.
            if (ih2optical.eq.1) then
                fudge = (0.76*d(i,j,k)*dom/8.0e9)**(-0.45)
                fudge = min(fudge, 1.0)
            else
                fudge = 1.0
            endif
            galdl(i) = gaHI(i) * HI(i,j,k)  + gaH2(i) * H2I(i,j,k)
     &               + gaHe(i) * HeI(i,j,k) + gaHp(i) * HII(i,j,k)
     &               + gael(i) * de(i,j,k)
            gphdl1 = gphdl(i)/dom
            edot(i) = edot(i) - float(ih2co)*fudge*H2I(i,j,k)*
     &           gphdl(i)/(1.0 + gphdl1/galdl(i)) / (2.0*dom)

            end if
         enddo
#else

#define USE_GALLI_PALLA1999
#define NO_OPTICAL_DEPTH_FUDGE

!        Use the Galli and Palla (1999) cooling rates for molecular H.

#ifdef USE_GALLI_PALLA1999

         do i = is+1, ie+1
            if ( itmask(i) ) then
            gpldl(i) = gpldla(indixe(i)) + tdef(i)
     &         *(gpldla(indixe(i)+1) - gpldla(indixe(i)))
            gphdl(i) = gphdla(indixe(i)) + tdef(i)
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then

#ifdef OPTICAL_DEPTH_FUDGE
            nH2 = 0.5*H2I(i,j,k)
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
     &                HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            fH2 = nH2/(nH2 + nother)
            fudge = sqrt((40.0 * 10**(4.8 * 
     &               sqrt(max(log10(tgas(i)),2.0)-2.0)) / fH2**2)/
     &               ((nH2 + nother)*dom) )
            fudge = min(fudge, 1.0)
#endif /* OPTICAL_DEPTH_FUDGE */
            ! Note that this optical depth approximation comes from
            ! RA04.
            if (ih2optical.eq.1) then
                fudge = (0.76*d(i,j,k)*dom/8.0e9)**(-0.45)
                fudge = min(fudge, 1.0)
            else
                fudge = 1.0
            endif
            gphdl1 = gphdl(i)/(HI(i,j,k)*dom)
            edot(i) = edot(i) - float(ih2co)*fudge*H2I(i,j,k)*
     &           gphdl(i)/(1.0 + gphdl1/gpldl(i)) / (2.0*dom)

            end if
         enddo

#else /* USE_GALLI_PALLA1999 */

         do i = is+1, ie+1
            if ( itmask(i) ) then
            hyd01k(i) = hyd01ka(indixe(i)) + tdef(i)
     &         *(hyd01ka(indixe(i)+1)-hyd01ka(indixe(i)))
            h2k01(i) = h2k01a(indixe(i)) + tdef(i)
     &         *(h2k01a(indixe(i)+1) - h2k01a(indixe(i)))
            vibh(i) = vibha(indixe(i)) + tdef(i)
     &         *(vibha(indixe(i)+1) - vibha(indixe(i)))
            roth(i) = rotha(indixe(i)) + tdef(i)
     &         *(rotha(indixe(i)+1) - rotha(indixe(i)))
            rotl(i) = rotla(indixe(i)) + tdef(i)
     &         *(rotla(indixe(i)+1) - rotla(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then
            qq   = 1.2*(HI(i,j,k)*dom)**0.77 + 
     &                (H2I(i,j,k)*dom/2.)**0.77
            vibl = (HI(i,j,k)*hyd01k(i) + 
     &             H2I(i,j,k)/2.*h2k01(i))
     &             *dom*8.18e-13

#ifdef OPTICAL_DEPTH_FUDGE
            nH2 = 0.5*H2I(i,j,k)
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
     &                HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            fH2 = nH2/(nH2 + nother)
            fudge = sqrt((40.0 * 10**(4.8 * 
     &               sqrt(max(log10(tgas(i)),2.0)-2.0)) / fH2**2)/
     &               ((nH2 + nother)*dom) )
            fudge = min(fudge, 1.0)
#endif /* OPTICAL_DEPTH_FUDGE */

            edot(i) = edot(i) - float(ih2co)*fudge*H2I(i,j,k)*(
     &          vibh(i)/(1.0+vibh(i)/max(   vibl     ,tiny)) +
     &          roth(i)/(1.0+roth(i)/max(qq*rotl(i),tiny))     
     &                                                      )/2.0/dom
            end if
         enddo

#endif /* USE_GALLI_PALLA1999 */
#endif /* USE_GLOVER_ABEL2008 */

c     CIE
c     cooling from H2-H2 and He-H2 collisional induced emission comes
C     with its own radiative transfer correction as discussed in
C     Ripamonti & Abel 2003
         if (iciecool.eq.1) then
            do i = is+1, ie+1
            if (itmask(i)) then
c     Only calculate if H2I(i) is a substantial fraction
              if (d(i,j,k)*dom.gt.1e10) then
                ciefudge = 1.
                tau = ((d(i,j,k)/2e16)*dom)**2.8  ! 2e16 is in units of cm^-3
                tau = max(tau, 1.e-5)
                ciefudge = min((1.-exp(-tau))/tau,1.0)
c               Matt's attempt at a second exponentialier cutoff
                tau = ((d(i,j,k)/2e18)*dom)**8.0d0  ! 2e18 is in units of cm^-3
                tau = max(tau, 1.e-5)
                ciefudge = ciefudge*min((1.-exp(-tau))/tau,1.0)
c               ciefudge, which is applied to the continuum, is applied to edot
                edot(i) = ciefudge*(edot(i) - 
     &                  H2I(i,j,k)*(d(i,j,k)*cieco(i)))
              endif
            endif
            enddo
         endif

      endif

!     --- Cooling from HD ---

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
c CMB cooling floor
               if (tgas(i) .gt. comp2) then
                  hdlte(i) = hdltea(indixe(i)) + tdef(i)
     &            *(hdltea(indixe(i)+1) - hdltea(indixe(i)))
                  hdlow(i) = hdlowa(indixe(i)) + tdef(i)
     &            *(hdlowa(indixe(i)+1) - hdlowa(indixe(i)))
               else
                  hdlte(i) = tiny
                  hdlow(i) = tiny
               endif
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then
c  old (incorrect) way:
c               hdlte1 = hdlte(i)/(HDI(i,j,k)*dom/2.0)
c               hdlow1 = max(hdlow(i), tiny)
c               edot(i) = edot(i) - HDI(i,j,k)*
c     .                     (hdlte1/(1.0 + hdlte1/hdlow1)/(2.0*dom))
c  new (correct) way: (april 4, 2007)
               hdlte1 = hdlte(i)/(HI(i,j,k)*dom)
               hdlow1 = max(hdlow(i), tiny)
               edot(i) = edot(i) - HDI(i,j,k)*
     .                     (hdlte(i)/(1.0 + hdlte1/hdlow1)) / (3.0*dom)
            end if
         enddo
      endif

!     --- Compute (external) radiative heating terms ---

!                       Photoionization heating

#ifdef RADIATION
      if (iradshield .eq. 0) then

!        regular version

         if (iradtype .eq. 8) then

!           1) heating assuming high energy photons produces secondary
!              electrons which do the heating (Shull & Steenberg, 1985).

            do i = is+1, ie+1
               if ( itmask(i) ) then
               x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.0e-4)
               factor = 0.9971*(1.0-(1.0-x**0.2663)**1.3163)
               edot(i) = edot(i) + float(ipiht)*factor*(
     &                + piHI  *HI  (i,j,k)         ! pi of HI
     &                + piHeI *HeI (i,j,k)*0.25     ! pi of HeI
     &                + piHeII*HeII(i,j,k)*0.25     ! pi of HeII
     &              )/dom
               end if
            enddo

         else

!           2) standard heating

            do i = is+1, ie+1
               if ( itmask(i) ) then
               edot(i) = edot(i) + float(ipiht)*(
     &                + piHI  *HI  (i,j,k)         ! pi of HI
     &                + piHeI *HeI (i,j,k)*0.25     ! pi of HeI
     &                + piHeII*HeII(i,j,k)*0.25     ! pi of HeII
     &              )/dom
               end if
            enddo

         endif

      else

!        version with approximate self-shielding

         do i = is+1, ie+1
            if ( itmask(i) ) then

            edot(i) = edot(i) + float(ipiht)*(
     &                + piHI  *HI  (i,j,k)*
     &                   exp(-avgsighp*HI(i,j,k)*dom)
     &                + piHeI *HeI (i,j,k)*0.25*
     &                   exp(-avgsighep*0.25*HeI(i,j,k)*dom)
     &                + piHeII*HeII(i,j,k)*0.25*
     &                   exp(-avgsighe2p*0.25*HeII(i,j,k)*dom)
     &           )/dom

!            write(6,*) 'cool1d_multi: after PI heating', edot(i),
!     &           piHI, HI(i,j,k), avgsighp, dom

            end if
         enddo

      endif

!     Photoheating from radiative transfer

      if (iradtrans .eq. 1) then
c     for radiative transfer, convert from eV/s*TimeUnits to the coolunits use
         rtunits = 1.60217646e-12/utim/coolunit/dom
         do i = is+1, ie+1
            if (itmask(i)) then
               edot(i) = edot(i) + float(ipiht) * photogamma(i,j,k) * 
     &              rtunits * HI(i,j,k)
c               if (photogamma(i,j,k).gt.0) then
c                  print*, i,j,k,edot(i), photogamma(i,j,k),rtunits,dom,
c     $                 aye,utim,coolunit,d(i,j,k),hi(i,j,k),tgas(i)
c               endif
               if (edot(i) .ne. edot(i)) then
                  write(6,*) 'NaN in edot[2]: ', i, j, k, edot(i), 
     &                 photogamma(i,j,k),HI(i,j,k), de(i,j,k), d(i,j,k), 
     &                 ge(i,j,k), p2d(i), tgas(i), dom, urho, aye, mh
                  ERROR_MESSAGE
               endif
            endif
         enddo
      endif

#endif /* RADIATION */


!     --- Cooling/heating due to metals ---

      if (imcool .eq. 1) then

         logxe0 = log(xe_start)
         logxe1 = log(xe_end)
         dlogxea = (logxe1 - logxe0) / real(n_xe-1)
c
c  Determine temperature for CMB temperature so we can subtract metal
c  cooling at T=T_cmb
c
         logtcmb = log(comp2)
         tcmb_idx = max(1, int((logtcmb-logtem0)/dlogtem)+1)
         dlogTcmb = (logtcmb - (logtem0 + (tcmb_idx-1)*dlogtem)) 
     $        / dlogtem
         do i = is+1, ie+1
            if (itmask(i)) then
!
               if (de(i,j,k) .lt. 0.0) then !#####
                  write(6,*) 'cool1d_multi: negative d(i,j,k)', 
     $                 de(i,j,k), d(i,j,k), metal(i,j,k)
               endif
!
c               metalfree = d(i,j,k) - metal(i,j,k)
               xe = de(i,j,k) / d(i,j,k)
c               if (tgas(i) .lt. 1e4) xe = min(xe, 1e-2)
               log_xe = max(min(log(xe), logxe1), logxe0)
               if (logtem(i) .lt. xe_logtem0) then
                  xe_max = max(xe_slope * (logtem(i) - xe_logtem0), 
     $                         xe_min)
                  log_xe = min(log_xe, xe_max)
               endif
               xe_idx = min(n_xe-1,
     $              max(1, int((log_xe - logxe0)/dlogxea)+1))
               xe1 = logxe0 + (xe_idx - 1)*dlogxea
               dlogT = tdef(i)
               dlogxe = (log_xe - xe1) / dlogxea
c
c  Interpolate metal cooling in temperature and electron fraction
c
               metalc(i) = 
     $          metala(indixe(i),   xe_idx  ) * (1-dlogxe) * (1-dlogT) +
     $          metala(indixe(i)+1, xe_idx  ) * (dlogxe  ) * (1-dlogT) +
     $          metala(indixe(i),   xe_idx+1) * (1-dlogxe) * (dlogT)   +
     $          metala(indixe(i)+1, xe_idx+1) * (dlogxe  ) * (dlogT)
c
c  Subtract metal cooling at T_cmb
c
            metalc_cmb = 
     $        metala(tcmb_idx,   xe_idx  ) * (1-dlogxe) * (1-dlogTcmb) +
     $        metala(tcmb_idx+1, xe_idx  ) * (dlogxe  ) * (1-dlogTcmb) +
     $        metala(tcmb_idx,   xe_idx+1) * (1-dlogxe) * (dlogTcmb)   +
     $        metala(tcmb_idx+1, xe_idx+1) * (dlogxe  ) * (dlogTcmb)
               metalc(i) = metalc(i) - metalc_cmb

               xi = min(3.0, metal(i,j,k)/(d(i,j,k)*ZSOLAR))
               edot(i) = edot(i) - 
     $              metalc(i) * d(i,j,k) * (xi/JHW_METALS_NORMZ) / dom
#ifdef UNUSED
               if (xi .gt. 1e-4 .and. tgas(i) .lt. 1e4) then
                  write(6,'(a,3i5,2e12.4)') 'a0', i,j,k, exp(log_xe),
     $                 exp(logtem(i))
                  write(6,'(a,4e12.4)') 'a1', d(i,j,k), xi, 
     $                 dlogxe, dlogT
                  write(6,'(a,i5,4e12.4)') 'a2', xe_idx, 
     $                 metala(indixe(i), xe_idx),
     $                 metalc(i), metalc_cmb, coolunit
                  write(6,'(a,2e12.4)') 'a3', edot(i),
     $                 metalc(i)*d(i,j,k)*(xi/JHW_METALS_NORMZ)/dom

                  write(6,*) ''
               endif
#endif

               if (edot(i) .ne. edot(i)) then    !#####
                  write(6,*) 'cool1d_multi: edot = NaN!', 
     $                        metalc(i), d(i,j,k), xi, dom, edot(i), 
     $                        dlogxe, log_xe, xe1, dlogxea,
     $                        logxe1, logxe0, n_xe
                  write(6,*) 'more: ',
     $                        log_xe, log(xe), logxe1, logxe0, xe, 
     $                        de(i,j,k), metal(i,j,k), d(i,j,k),  
     $                        HII(i,j,k), HeII(i,j,k), HeIII(i,j,k),
     $                        xe_max, xe_min
                  edot(i) = 0.0
               endif

            endif
         enddo

      endif                     ! imcool == 1

#ifdef CEN_METALS

      if (imcool .eq. 2) then

!     Generate table if required

      if (imetalregen .eq. 1) then
         write(6,*) 'generating metallicity cooling table'
         if ( (iradtype .eq. 10) .or. (iradtype .eq. 11) ) iradfield = 1

!        Clear table

         do j1=1, NID
            denwk(j1) = 10.0**(DENMIN+(j1-0.5)*DELD)
            do i=1, NIT
               cbovcool(i,j1) = 0
               cbovheat(i,j1) = 0
            enddo
         enddo

!        compute energy (in eV) of each bin

!        write(6,*) 'radfield type ', iradfield  
!        write(6,*) 'radtype value ', iradtype

         if (iradfield .eq. 1) then

            if (NIB .ne. nfreq) then
               write(6,*) NIB, nfreq
               write(6,*) 'cool1d_multi: NIB != nfreq'
               ERROR_MESSAGE
            endif

            do i=1, NIB
               eb(i) = FREQMIN*10.0**(FREQDEL*(i-1.0))
            enddo

!        Convert to 10^{-21} erg/cm^2/s

            do i=1, NIB-1
               radt(i) = inutot(i)*4.0*3.14159*(eb(i+1)-eb(i))
     &                   *(1.2e11/6.625)
            enddo

         else

            do i=1, NIB-1
               radt(i) = 0.0
            enddo

         endif

         radt(NIB) = 0

!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!         call open_mpi_error_file('Metal', 57, 'unknown')
         call mtlclht(denwk, radt, iradfield, cbovcool, cbovheat)
!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!         call close_mpi_error_file( 57 )

         imetalregen = 0

!        Convert from 10^-23 to our units

         do j1=1, NID
            do i=1, NIT
               cbovcool(i,j1) = cbovcool(i,j1) * (1d-23/coolunit)  ! DPC
               cbovheat(i,j1) = cbovheat(i,j1) * (1d-23/coolunit)  ! DPC
!              write(31,*) i,j1,cbovcool(i,j1),cbovheat(i,j1)
            enddo
         enddo

      endif

!     Look-up in table

      do i = is+1, ie+1
         if ( itmask(i) ) then

         ndi = (log10(de(i,j,k)*dom)-DENMIN)/DELD + 1
         ndi = min(max(1, ndi), NID)
         nti = (log10(0.5*(tgas(i)+tgasold(i)))-TEMMIN)/DELT + 1
         nti = min(max(1, nti), NIT)
         xi = min(3.0, metal(i,j,k)/(d(i,j,k)*ZSOLAR))
         metal_cool = cbovcool(nti,ndi) * xi / dom *
     &                (HI(i,j,k)+HII(i,j,k))
         metal_heat = cbovheat(nti,ndi) * xi / dom *
     &                (HI(i,j,k)+HII(i,j,k))
         if (tgas(i) .le. 1e4) then
            metal_cool = metal_cool*(tgas(i)/1e4)**3
            metal_heat = 0.0
         endif

         edot(i) = edot(i) + metal_heat - metal_cool

         end if
      enddo

      endif

#endif /* CEN_METALS */

!     --- Cloudy metal cooling and heating ---

      if (imcool .eq. 3) then

         call cool1D_cloudy(d, de, HI, HII, HeI, HeII, HeIII,
     &        HM, H2I, H2II, DI, DII, HDI, metal,
     &        in, jn, kn, is, ie, j, k,
     &        logtem, edot, comp2, ispecies, dom, zr,
     &        icmbTfloor, iClHeat, iClMMW,
     &        clMetNorm, clEleFra, clGridRank, clGridDim,
     &        clPar1, clPar2, clPar3, clPar4, clPar5,
     &        clDataSize, clCooling, clHeating, clMMW,
     &        itmask)

      endif

!     Set tgasold

      do i=is+1, ie+1
         if ( itmask(i) ) then
         tgasold(i) = tgas(i)
         end if
      enddo

      return
      end

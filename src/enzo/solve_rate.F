#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c/////////////////////  SUBROUTINE SOLVE_RATE  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine solve_rate(
     &                de, HI, HII, HeI, HeII, HeIII, tgas, d,
     &                in, jn, kn, nratec, ispecies,
     &                is, js, ks, ie, je, ke,
     &                dt, aye, temstart, temend, 
     &                uaye, utim, urho, uxyz, fh, dtoh,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a,
     &                k24, k25, k26, k27, k28, k29, k30, k31,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                HM, H2I, H2II, DI, DII, HDI,
     &                iradshield, avgsigh, avgsighe, avgsighe2,
     &                iradtype, piHI, piHeI,
     &                iradtrans, irt_honly, kphHI, kphHeI, kphHeII, 
     &                kdissH2I)
c
c  SOLVE MULTI-SPECIES RATE EQUATIONS
c
c  written by: Yu Zhang, Peter Anninos and Tom Abel
c  date:       
c  modified1:  January, 1996 by Greg Bryan; converted to KRONOS
c  modified2:  October, 1996 by GB; adapted to AMR
c  modified3:  May,     1999 by GB and Tom Abel, 3bodyH2, solver, HD
c  modified4:  November 2003 by Robert Harkness; tighten tolerance
c
c  PURPOSE:
c    Solve the multi-species rate equations.
c
c  INPUTS:
c
c  PARAMETERS:
c    itmax   - maximum allowed sub-cycle iterations
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c  Arguments
c
      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke, nratec, ispecies,
     &        iradshield, iradtype, iradtrans, irt_honly
      R_PREC    dt, aye, temstart, temend, utim, uaye, urho, uxyz,
     &        fh, dtoh
      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &      tgas(in,jn,kn),    d(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      R_PREC kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &     kdissH2I(in,jn,kn)
      R_PREC k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &     k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &     k54a(nratec), k55a(nratec), k56a(nratec),
     &     k13dda(nratec, 7),
     &     k24, k25, k26, k27, k28, k29, k30, k31,
     &     avgsigh, avgsighe, avgsighe2, piHI, piHeI
c
c  Parameters
c
      INTG_PREC itmax, ijk
      parameter (itmax = 5000, ijk = MAX_ANY_SINGLE_DIRECTION)

#ifdef CONFIG_BFLOAT_4
      R_PREC tolerance
      parameter (tolerance = 1.e-5_RKIND)
#endif

#ifdef CONFIG_BFLOAT_8
      R_PREC tolerance
      parameter (tolerance = 1.e-10_RKIND)
#endif

      real*8 everg, mh, e24, e26
      parameter (everg = ev2erg, mh = mass_h,
     &           e24 = 13.6_RKIND, e26 = 24.6_RKIND)            ! DPC
c
c  Locals
c
      INTG_PREC i, j, k, n, iter, n1
      R_PREC qq, vibl, logtem0, logtem9, dlogtem
      R_PREC ttmin, comp1, comp2, c2, scoef, acoef, nh, dom, factor, x
      R_PREC HIp(ijk), HIIp(ijk), HeIp(ijk), HeIIp(ijk), HeIIIp(ijk),
     &     HMp(ijk), H2Ip(ijk), H2IIp(ijk),
     &     dep(ijk), dedot(ijk),HIdot(ijk),
     &     DIp(ijk), DIIp(ijk), HDIp(ijk),
     &     k24shield(ijk), k25shield(ijk), k26shield(ijk),
     &     totalH(ijk), totalHe(ijk), totalD, 
     &     correctH, correctHe, correctD, metalfree

      real*8 coolunit, dbase1, tbase1, xbase1
c
c  Slice temporaries (passed in)
c 
      INTG_PREC indixe(ijk)
      R_PREC t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), 
     &     dtit(ijk), ttot(ijk)
c
c  Rate equation slice temporaries (passed in)
c
      R_PREC k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k22(ijk),
     &     k50(ijk), k51(ijk), k52(ijk), k53(ijk), k54(ijk),
     &     k55(ijk), k56(ijk), k13dd(ijk, 7)
c
c      R_PREC tgas(in,jn,kn), tgasold(in,jn,kn)
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
c
c     Set log values of start and end of lookup tables
c
      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem = (log(temend) - log(temstart))/REAL(nratec-1,RKIND)
c
c  Convert densities from comoving to proper
c
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               de(i,j,k)    = de(i,j,k)/aye**3
               HI(i,j,k)    = HI(i,j,k)/aye**3
               HII(i,j,k)   = HII(i,j,k)/aye**3
               HeI(i,j,k)   = HeI(i,j,k)/aye**3
               HeII(i,j,k)  = HeII(i,j,k)/aye**3
               HeIII(i,j,k) = HeIII(i,j,k)/aye**3
            enddo
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)/aye**3
                  H2I(i,j,k)  = H2I(i,j,k)/aye**3
                  H2II(i,j,k) = H2II(i,j,k)/aye**3
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)/aye**3
                  DII(i,j,k) = DII(i,j,k)/aye**3
                  HDI(i,j,k) = HDI(i,j,k)/aye**3
               enddo
            endif
         enddo
      enddo
c
c
c  B) Solve rate equations using BDF
c
      do k = ks+1, ke+1
      do j = js+1, je+1
c
         do i = is+1, ie+1
c
            ttot(i) = 0._RKIND
c
c           Compute temp-centered temperature (and log)
c
c           logtem(i) = log(0.5_RKIND*(tgas(i,j,k)+tgasold(i,j,k)))
            logtem(i) = log(tgas(i,j,k))
            logtem(i) = max(logtem(i), logtem0)
            logtem(i) = min(logtem(i), logtem9)
c
c           Find index into tble and precompute interpolation values
c
            indixe(i) = min(nratec-1, max(1,
     &           int((logtem(i)-logtem0)/dlogtem,IKIND)+1))
            t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
            t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
            tdef(i) = t2(i) - t1(i)
c
c           Do linear table lookup (in log temperature)
c
            k1(i) = k1a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k1a(indixe(i)+1) -k1a(indixe(i)))/tdef(i)
            k2(i) = k2a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k2a(indixe(i)+1) -k2a(indixe(i)))/tdef(i)
            k3(i) = k3a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k3a(indixe(i)+1) -k3a(indixe(i)))/tdef(i)
            k4(i) = k4a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k4a(indixe(i)+1) -k4a(indixe(i)))/tdef(i)
            k5(i) = k5a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k5a(indixe(i)+1) -k5a(indixe(i)))/tdef(i)
            k6(i) = k6a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k6a(indixe(i)+1) -k6a(indixe(i)))/tdef(i)
         enddo
c
c        Look-up for 9-species model
c
         if (ispecies .gt. 1) then
            do i = is+1, ie+1
               k7(i) = k7a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k7a(indixe(i)+1) -k7a(indixe(i)))/tdef(i)
               k8(i) = k8a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k8a(indixe(i)+1) -k8a(indixe(i)))/tdef(i)
               k9(i) = k9a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k9a(indixe(i)+1) -k9a(indixe(i)))/tdef(i)
               k10(i) = k10a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k10a(indixe(i)+1) -k10a(indixe(i)))/tdef(i)
               k11(i) = k11a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k11a(indixe(i)+1) -k11a(indixe(i)))/tdef(i)
               k12(i) = k12a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k12a(indixe(i)+1) -k12a(indixe(i)))/tdef(i)
               k13(i) = k13a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k13a(indixe(i)+1) -k13a(indixe(i)))/tdef(i)
               k14(i) = k14a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k14a(indixe(i)+1) -k14a(indixe(i)))/tdef(i)
               k15(i) = k15a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k15a(indixe(i)+1) -k15a(indixe(i)))/tdef(i)
               k16(i) = k16a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k16a(indixe(i)+1) -k16a(indixe(i)))/tdef(i)
               k17(i) = k17a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k17a(indixe(i)+1) -k17a(indixe(i)))/tdef(i)
               k18(i) = k18a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k18a(indixe(i)+1) -k18a(indixe(i)))/tdef(i)
               k19(i) = k19a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k19a(indixe(i)+1) -k19a(indixe(i)))/tdef(i)
               k22(i) = k22a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k22a(indixe(i)+1) -k22a(indixe(i)))/tdef(i)
            enddo
            do n1 = 1, 7
               do i = is+1, ie+1
                  k13dd(i,n1) = k13dda(indixe(i),n1) + 
     &                          (logtem(i) - t1(i))
     &                           *(k13dda(indixe(i)+1,n1) -
     &                             k13dda(indixe(i)  ,n1))/tdef(i)
               enddo
            enddo
         endif
c
c        Look-up for 12-species model
c
         if (ispecies .gt. 2) then
            do i = is+1, ie+1
               k50(i) = k50a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k50a(indixe(i)+1) -k50a(indixe(i)))/tdef(i)
               k51(i) = k51a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k51a(indixe(i)+1) -k51a(indixe(i)))/tdef(i)
               k52(i) = k52a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k52a(indixe(i)+1) -k52a(indixe(i)))/tdef(i)
               k53(i) = k53a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k53a(indixe(i)+1) -k53a(indixe(i)))/tdef(i)
               k54(i) = k54a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k54a(indixe(i)+1) -k54a(indixe(i)))/tdef(i)
               k55(i) = k55a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k55a(indixe(i)+1) -k55a(indixe(i)))/tdef(i)
               k56(i) = k56a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k56a(indixe(i)+1) -k56a(indixe(i)))/tdef(i)
            enddo
         endif
c
c        Include approximate self-shielding factors if requested
c          (avgsigh already has factor to convert code density to
c           particles/cm^3)
c
#ifdef RADIATION
         do i = is+1, ie+1
            k24shield(i) = k24
            k25shield(i) = k25
            k26shield(i) = k26
         enddo
         if (iradshield .eq. 1) then
            do i = is+1, ie+1
               k24shield(i) = k24shield(i)*exp(-HI(i,j,k)*avgsigh)
               k25shield(i) = k25shield(i)*exp(-HeII(i,j,k)*avgsighe2)
               k26shield(i) = k26shield(i)*exp(-HeI(i,j,k)*avgsighe)
            enddo
         endif
c
c        If using a high-energy radiation field, then account for
c          effects of secondary elections (Shull * Steenberg 1985)
c          (see calc_rate.src)
c
         if (iradtype .eq. 8) then
            do i = is+1, ie+1
               x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.e-4_RKIND)
               factor = 0.3908_RKIND*(1._RKIND - 
     &              x**0.4092_RKIND)**1.7592_RKIND
               k24shield(i) = k24shield(i) + 
     &              factor*(piHI + 0.08_RKIND*piHeI) / 
     &              (e24*everg)*coolunit*tbase1
               factor = 0.0554_RKIND*(1._RKIND - 
     &              x**0.4614_RKIND)**1.666_RKIND
               k26shield(i) = k26shield(i) + 
     &              factor*(piHI/0.08_RKIND + piHeI) / 
     &              (e26*everg)*coolunit*tbase1
            enddo
         endif
c
#endif /* RADIATION */
c
c        ------------------ Loop over subcycles ----------------
c
         do iter = 1, itmax
c
c           If using H2, and using the density-dependent collisional
c             H2 dissociation rate, then replace the the density-independant
c                k13 rate with the new one.
c         May/00: there appears to be a problem with the density-dependent
c             collisional rates.  Currently turned off until further notice.
c
#define USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
#ifdef USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  nh = min(HI(i,j,k)*dom, 1.e9_RKIND)
                  k13(i) = tiny
                  if (tgas(i,j,k) .ge. 500._RKIND .and.
     &                tgas(i,j,k) .lt. 1.e6_RKIND) then
                     k13(i) = k13dd(i,1)-k13dd(i,2)/
     &                          (1._RKIND+(nh/k13dd(i,5))**k13dd(i,7))
     &                     + k13dd(i,3)-k13dd(i,4)/
     &                          (1._RKIND+(nh/k13dd(i,6))**k13dd(i,7))
                     k13(i) = max(10._RKIND**k13(i), tiny)
                  endif
               enddo
            endif
#endif /*  USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE */
c
c
c           --- (1) compute timestep ---
c     
            if (ispecies .eq. 1) then 
               do i = is+1, ie+1
c
c     Compute the electron density rate-of-change
c
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
c
c     Compute the HI density rate-of-change
c     
                  HIdot(i) =
     &               - k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k2(i)*HII(i,j,k)*de(i,j,k)
#ifdef RADIATION
     &               -      k24shield(i)*HI(i,j,k)
#endif /* RADIATION */
c     
               enddo
            else
c
c         Include molecular hydrogen rates for HIdot
c
               do i = is+1, ie+1
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
c
c     Compute the electron density rate-of-change
c
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
c
#ifdef RADIATION
     &               +      ( k24shield(i)*HI(i,j,k)
     &               + k25shield(i)*HeII(i,j,k)/4._RKIND
     &               + k26shield(i)*HeI(i,j,k)/4._RKIND)
#endif /* RADIATION */
               enddo
            endif
c
c     Add photo-ionization rates if necessary
c
#ifdef RADIATION
            if (iradtrans.eq.1) then
               if (irt_honly.eq.0) then
                  do i = is+1, ie+1
                     HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                     dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
     $                    + kphHeI(i,j,k) * HeI(i,j,k) / 4._RKIND
     $                    + kphHeII(i,j,k) *HeII(i,j,k) / 4._RKIND
                  enddo
               else
                  do i = is+1, ie+1
                     HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                     dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
                  enddo
               endif
            endif
#endif /* RADIATION */
c
c           Find timestep that keeps relative changes below 10%
c
            do i = is+1, ie+1
c
c              Bound from below to prevent numerical errors
c
               if (abs(dedot(i)) .lt. tiny) dedot(i) = tiny
               if (abs(HIdot(i)) .lt. tiny) HIdot(i) = tiny
c
c              If the net rate is almost perfectly balanced then set
c                  it to zero (since it is zero to available precision)
c
               if (min(abs(k1(i)* de(i,j,k)*HI(i,j,k)),
     &                 abs(k2(i)*HII(i,j,k)*de(i,j,k)))/
     &             max(abs(dedot(i)),abs(HIdot(i))) .gt. 1.0d6) then
                  dedot(i) = tiny
                  HIdot(i) = tiny
               endif
c
               dtit(i) = min(abs(0.1*de(i,j,k)/dedot(i)), 
     &                       abs(0.1*HI(i,j,k)/HIdot(i)),
     &                       dt-ttot(i), 0.5_RKIND*dt)
c
               if (dtit(i)/dt .lt. 1.e-2_RKIND .and. iter .gt. 200 .and.
     &             abs((dt-ttot(i))/dt) .gt. 1.e-3_RKIND) then

!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!               call open_mpi_error_file( 'F4', 4, 'unknown' )

                  write(4,1000) iter,i,j,k,dtit(i),
     &              ttot(i),dt,de(i,j,k),dedot(i),HI(i,j,k),HIdot(i),
     &              tgas(i,j,k)
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
     &               + 2._RKIND*k13(i)*HI(i,j,k)  *H2I(i,j,k)/2._RKIND,
     &               +      k11(i)*HII(i,j,k)   *H2I(i,j,k)/2._RKIND,
     &               + 2._RKIND*k12(i)*de(i,j,k)  *H2I(i,j,k)/2._RKIND,
     &               +      k14(i)*HM(i,j,k)    *de(i,j,k),
     &               +      k15(i)*HM(i,j,k)    *HI(i,j,k),
     &               + 2._RKIND*k16(i)*HM(i,j,k)    *HII(i,j,k),
     &               + 2._RKIND*k18(i)*H2II(i,j,k)  *de(i,j,k)/2._RKIND,
     &               +      k19(i)*H2II(i,j,k)  *HM(i,j,k)/2._RKIND

!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!               call close_mpi_error_file( 4 )

               endif
 1000          format(i5,3(i3,1x),1p,8(e11.3))
 1100          format(20(e11.3))
            enddo
c
c         --- (2) solve 6-species rate equation with one linearly implicit
c             Gauss Seidel sweep of a backward Euler time integrator ---
c
            if (ispecies .eq. 1) then
c
               do i = is+1, ie+1
c
c                 1) HI
c
                  scoef  = k2(i)*HII(i,j,k)*de(i,j,k)
                  acoef  = k1(i)*de(i,j,k)
#ifdef RADIATION
     &             + k24shield(i)
                  if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
#endif /* RADIATION */
                  HIp(i)    = ( scoef*dtit(i) + HI(i,j,k) ) / 
     &                        ( 1._RKIND + acoef*dtit(i) )
c
c                 2) HII
c 
                  scoef  = k1(i)*HIp(i)*de(i,j,k)
#ifdef RADIATION
     &                   + k24shield(i)*HIp(i)
                  if (iradtrans .eq. 1) 
     $                 scoef = scoef + kphHI(i,j,k)*HIp(i)
#endif /* RADIATION */
                  acoef  = k2(i)*de (i,j,k)
                  HIIp(i)  = ( scoef*dtit(i) + HII(i,j,k) )
     &                     / ( 1._RKIND + acoef*dtit(i) )
c
c                 3) Electron density
c
                  scoef = 0._RKIND
#ifdef RADIATION
     &                 + k24shield(i)*HI(i,j,k)
     &                 + k25shield(i)*HeII(i,j,k)/4._RKIND
     &                 + k26shield(i)*HeI(i,j,k)/4._RKIND
                  if (iradtrans .eq. 1 .and. irt_honly.eq.0) 
     $                 scoef = scoef + kphHI(i,j,k)   * HI(i,j,k)
     $                 + kphHeI(i,j,k)  * HeI(i,j,k)/4._RKIND
     $                 + kphHeII(i,j,k) * HeII(i,j,k)/4._RKIND
                  if (iradtrans .eq. 1 .and. irt_honly.eq.1) 
     $                 scoef = scoef + kphHI(i,j,k)   * HI(i,j,k)
#endif /* RADIATION */
                  acoef = -(k1(i)*HI(i,j,k)      
     &                    - k2(i)*HII(i,j,k)
     &                    + k3(i)*HeI(i,j,k)/4._RKIND
     &                    - k6(i)*HeIII(i,j,k)/4._RKIND
     &                    + k5(i)*HeII(i,j,k)/4._RKIND
     &                    - k4(i)*HeII(i,j,k)/4._RKIND
     &                )
                  dep(i)   = ( scoef*dtit(i) + de(i,j,k) )
     &                     / ( 1._RKIND + acoef*dtit(i) )
c
               enddo
c
            endif               ! (ispecies .eq. 1)
c 
c           do helium chemistry in any case:
c
            do i = is+1, ie+1
c
c              4) HeI
c 
               scoef  = k4(i)*HeII(i,j,k)*de(i,j,k)
               acoef  = k3(i)*de(i,j,k)
#ifdef RADIATION
     &                + k26shield(i)
               if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $              acoef = acoef + kphHeI(i,j,k)
#endif /* RADIATION */
               HeIp(i)   = ( scoef*dtit(i) + HeI(i,j,k) ) 
     &              / ( 1._RKIND + acoef*dtit(i) )
c
c              5) HeII
c
               scoef  = k3(i)*HeIp(i)*de(i,j,k)
     &                + k6(i)*HeIII(i,j,k)*de(i,j,k)
#ifdef RADIATION
     &                + k26shield(i)*HeIp(i)
               if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $              scoef = scoef + kphHeI(i,j,k)*HeIp(i)
#endif /* RADIATION */
               acoef  = k4(i)*de(i,j,k) + k5(i)*de(i,j,k)
#ifdef RADIATION
     &                + k25shield(i)
               if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $              acoef = acoef + kphHeII(i,j,k)
#endif /* RADIATION */
               HeIIp(i)  = ( scoef*dtit(i) + HeII(i,j,k) )
     &              / ( 1._RKIND + acoef*dtit(i) )
c
c              6) HeIII
c
               scoef   = k5(i)*HeIIp(i)*de(i,j,k)
#ifdef RADIATION
     &                 + k25shield(i)*HeIIp(i)
               if (iradtrans.eq.1 .and. irt_honly.eq.0)
     $              scoef = scoef + kphHeII(i,j,k) * HeIIp(i)
#endif /* RADIATION */
               acoef   = k6(i)*de(i,j,k)
               HeIIIp(i)  = ( scoef*dtit(i) + HeIII(i,j,k) )
     &              / ( 1._RKIND + acoef*dtit(i) )
            enddo
c
c         --- (3) Now do extra 3-species for molecular hydrogen ---
c
            if (ispecies .gt. 1) then
c
c              First, do HI/HII with molecular hydrogen terms
c
               do i = is+1, ie+1
c
c                 1) HI
c
                  scoef  =      k2(i) * HII(i,j,k) * de(i,j,k) 
     &                   + 2._RKIND*k13(i)* HI(i,j,k) 
     &                     * H2I(i,j,k)/2._RKIND
     &                   +      k11(i)* HII(i,j,k) * H2I(i,j,k)/2._RKIND
     &                   + 2._RKIND*k12(i)* de(i,j,k)  
     &                     * H2I(i,j,k)/2._RKIND
     &                   +      k14(i)* HM(i,j,k)  * de(i,j,k)
     &                   +      k15(i)* HM(i,j,k)  * HI(i,j,k)
     &                   + 2._RKIND*k16(i)* HM(i,j,k)  * HII(i,j,k)
     &                   + 2._RKIND*k18(i)* H2II(i,j,k)
     &                     * de(i,j,k)/2._RKIND
     &                   +      k19(i)* H2II(i,j,k)* HM(i,j,k)/2._RKIND
#ifdef RADIATION
     &                   + 2._RKIND*k31   * H2I(i,j,k)/2._RKIND
                  if (iradtrans.eq.1)
     $                 scoef = scoef + 2._RKIND*kdissH2I(i,j,k)* 
     $                         H2I(i,j,k)/2._RKIND
#endif /* RADIATION */
                  acoef  =      k1(i) * de(i,j,k)
     &                   +      k7(i) * de(i,j,k)  
     &                   +      k8(i) * HM(i,j,k)
     &                   +      k9(i) * HII(i,j,k)
     &                   +      k10(i)* H2II(i,j,k)/2._RKIND
     &                   + 2._RKIND*k22(i)* HI(i,j,k)**2
#ifdef RADIATION
     &                   + k24shield(i)
                  if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
#endif /* RADIATION */
                  HIp(i)  = ( scoef*dtit(i) + HI(i,j,k) ) / 
     &                      ( 1._RKIND + acoef*dtit(i) )
c
c                 2) HII
c 
                  scoef  =    k1(i)  * HI(i,j,k) * de(i,j,k)
     &                   +    k10(i) * H2II(i,j,k)*HI(i,j,k)/2._RKIND
#ifdef RADIATION
     &                   + k24shield(i)*HI(i,j,k)
                  if (iradtrans .eq. 1) 
     $                 scoef = scoef + kphHI(i,j,k)*HI(i,j,k)
#endif /* RADIATION */
                  acoef  =    k2(i)  * de(i,j,k)
     &                   +    k9(i)  * HI(i,j,k)
     &                   +    k11(i) * H2I(i,j,k)/2._RKIND
     &                   +    k16(i) * HM(i,j,k)
     &                   +    k17(i) * HM(i,j,k)
                  HIIp(i)   = ( scoef*dtit(i) + HII(i,j,k) )
     &                      / ( 1._RKIND + acoef*dtit(i) )
c
c                 3) electrons:
c
                  scoef =   k8(i) * HM(i,j,k) * HI(i,j,k)
     &                   +  k15(i)* HM(i,j,k) * HI(i,j,k)
     &                   +  k17(i)* HM(i,j,k) * HII(i,j,k)
c                  
#ifdef RADIATION
     &                   + k24shield(i)*HIp(i)
     &                   + k25shield(i)*HeIIp(i)/4._RKIND
     &                   + k26shield(i)*HeIp(i)/4._RKIND
                  if (iradtrans .eq. 1 .and. irt_honly.eq.0) 
     $                 scoef = scoef + kphHI(i,j,k)   * HIp(i)
     $                 + kphHeI(i,j,k)  * HeIp(i)/4._RKIND
     $                 + kphHeII(i,j,k) * HeIIp(i)/4._RKIND
                  if (iradtrans .eq. 1 .and. irt_honly.eq.1) 
     $                 scoef = scoef + kphHI(i,j,k)   * HIp(i)
#endif /* RADIATION */
                  acoef = - (k1(i) *HI(i,j,k)    
     &                   -  k2(i)*HII(i,j,k)
     &                   +  k3(i) *HeI(i,j,k)/4._RKIND
     &                   -  k6(i)*HeIII(i,j,k)/4._RKIND
     &                   +  k5(i) *HeII(i,j,k)/4._RKIND
     &                   -  k4(i)*HeII(i,j,k)/4._RKIND
     &                   +  k14(i)*HM(i,j,k)
     &                   -  k7(i) *HI(i,j,k)
     &                   -  k18(i)*H2II(i,j,k)/2._RKIND
     &                 )
                  dep(i)  = ( scoef*dtit(i) + de(i,j,k) )
     &                    / ( 1._RKIND + acoef*dtit(i) )
c
c                 7) H2
c
                  scoef = 2._RKIND*(k8(i)  * HM(i,j,k)   * HI(i,j,k)
     &                  +     k10(i) * H2II(i,j,k) * HI(i,j,k)/2._RKIND
     &                  +     k19(i) * H2II(i,j,k) * HM(i,j,k)/2._RKIND
     &                  +     k22(i) * HI(i,j,k)**3             
     &                 )
                  acoef = ( k13(i)*HI(i,j,k) + k11(i)*HII(i,j,k)
     &                  + k12(i)*de(i,j,k) )
#ifdef RADIATION
     &                  + k29 + k31
                  if (iradtrans.eq.1) 
     &                 acoef = acoef+2._RKIND*kdissH2I(i,j,k)
#endif /* RADIATION */
c
                  H2Ip(i) = ( scoef*dtit(i) + H2I(i,j,k) )
     &                    / ( 1._RKIND + acoef*dtit(i) )
c
c                 8) H-
c
                  HMp(i) = ( k7(i)*HIp(i)*dep(i) )
     &                   / ( (k8(i)+k15(i))*HIp(i)
     &                   + ( k16(i)+k17(i))*HIIp(i)+k14(i)*dep(i)
#ifdef RADIATION
     &                   + k27
#endif /* RADIATION */
     &                 )
c
c                 9) H2+
c
                  H2IIp(i) = 2._RKIND*( k9 (i)*HIp(i)*HIIp(i)
     &                          + k11(i)*H2Ip(i)/2._RKIND*HIIp(i)
     &                          + k17(i)*HMp(i)*HIIp(i)
#ifdef RADIATION
     &                          + k29*H2Ip(i)
#endif /* RADIATION */
     &                          )
     &                        / ( k10(i)*HIp(i) + k18(i)*dep(i)
     &                          + k19(i)*HMp(i)
#ifdef RADIATION
     &                          + (k28+k30)
#endif /* RADIATION */
     &                          )
               enddo
c     
            endif               ! H2
c
c         --- (4) Now do extra 3-species for molecular HD ---
c     
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
c     
c                 1) DI
c     
                  scoef =   (    k2(i)* DII(i,j,k) * de(i,j,k)
     &                      +    k51(i)* DII(i,j,k) * HI(i,j,k)
     &                      + 2._RKIND*k55(i)* HDI(i,j,k) 
     &                        * HI(i,j,k)/3._RKIND
     &                 )
                  acoef  =    k1(i) * de(i,j,k)
     &                 +     k50(i) * HII(i,j,k)
     &                 +     k54(i) * H2I(i,j,k)/2._RKIND
     &                 +     k56(i) * HM(i,j,k)
#ifdef RADIATION
     &                 + k24shield(i)
                  if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
#endif /* RADIATION */
                  DIp(i)    = ( scoef*dtit(i) + DI(i,j,k) ) / 
     &                 ( 1._RKIND + acoef*dtit(i) )
c
c                 2) DII
c 
                  scoef = (   k1(i) * DI(i,j,k) * de(i,j,k)
     &                  +     k50(i) * HII(i,j,k)* DI(i,j,k)
     &                  +  2._RKIND*k53(i) * HII(i,j,k)
     &                     * HDI(i,j,k)/3._RKIND
     &                  )
#ifdef RADIATION
     &                 + k24shield(i)*DI(i,j,k)
                  if (iradtrans .eq. 1)
     $                 scoef = scoef + kphHI(i,j,k)*DI(i,j,k)
#endif /* RADIATION */
                  acoef =    k2(i)  * de(i,j,k)
     &                  +    k51(i) * HI(i,j,k)
     &                  +    k52(i) * H2I(i,j,k)/2._RKIND
c
                  DIIp(i)   = ( scoef*dtit(i) + DII(i,j,k) )
     &                 / ( 1._RKIND + acoef*dtit(i) )
c
c                 3) HDI
c 
                  scoef = 3._RKIND*(k52(i)*DII(i,j,k)
     &                      * H2I(i,j,k)/2._RKIND/2._RKIND
     &                   +    k54(i) * DI(i,j,k) * H2I(i,j,k)/4._RKIND
     &                   + 2._RKIND*k56(i) *DI(i,j,k) 
     &                      * HM(i,j,k)/2._RKIND
     &                 )
                  acoef  =    k53(i) * HII(i,j,k)
     &                   +    k55(i) * HI(i,j,k)
c
                  HDIp(i)   = ( scoef*dtit(i) + HDI(i,j,k) )
     &                 / ( 1._RKIND + acoef*dtit(i) )

               enddo
            endif
c
c         --- (5) Set densities and update elapsed time ---
c
            do i = is+1, ie+1
c  
               HI(i,j,k)    = HIp(i)
               HII(i,j,k)   = HIIp(i)
               HeI(i,j,k)   = HeIp(i)
               HeII(i,j,k)  = HeIIp(i)
               HeIII(i,j,k) = max(HeIIIp(i), 1.e-25_RKIND)
c
c               de(i,j,k)    = dep(i)
c
c              Use charge conservation to determine electron fraction
c 
               de(i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND + 
     &              HeIII(i,j,k)/2._RKIND
               if (ispecies .gt. 1) de(i,j,k) = de(i,j,k)
     &                   - HM(i,j,k) + H2II(i,j,k)/2._RKIND
c
c
               if (ispecies .gt. 1) then
                  HM(i,j,k)    = HMp(i)
                  H2I(i,j,k)   = H2Ip(i)
                  H2II(i,j,k)  = H2IIp(i)
               endif
c
               if (ispecies .gt. 2) then
                  DI(i,j,k)    = DIp(i)
                  DII(i,j,k)   = DIIp(i)
                  HDI(i,j,k)   = HDIp(i)
               endif
c
c              Add the timestep to the elapsed time for each cell
c     
               ttot(i) = ttot(i) + dtit(i)
c     
            enddo               ! end loop over i
c
c           Find the lowest elapsed time step in this slice
c     
            ttmin = huge
            do i = is+1, ie+1
               if (ttot(i).lt.ttmin) ttmin = ttot(i)
            enddo
c     
c           If all cells are done (on this slice), then exit
c     
            if (abs(dt-ttmin).lt. tolerance*dt) go to 9999
c     
c           Next subcycle iteration
c     
         enddo
c     
 9999    continue
c     
c        Check to see if we exceed the maximum iterations
c     
         if (iter .gt. itmax ) then
            write(6,*) 'RATE iter exceeds ',itmax,' at j,k =',j,k
            ERROR_MESSAGE
         endif
c     
c     Next j,k
c     
      enddo
      enddo
c
c     Convert densities back to comoving from proper
c
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               de(i,j,k)    = de(i,j,k)    * aye**3
               HI(i,j,k)    = HI(i,j,k)    * aye**3
               HII(i,j,k)   = HII(i,j,k)   * aye**3
               HeI(i,j,k)   = HeI(i,j,k)   * aye**3
               HeII(i,j,k)  = HeII(i,j,k)  * aye**3
               HeIII(i,j,k) = HeIII(i,j,k) * aye**3
            enddo
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)  * aye**3
                  H2I(i,j,k)  = H2I(i,j,k) * aye**3
                  H2II(i,j,k) = H2II(i,j,k)* aye**3
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)   * aye**3
                  DII(i,j,k) = DII(i,j,k)  * aye**3
                  HDI(i,j,k) = HDI(i,j,k)  * aye**3
               enddo
            endif
         enddo
      enddo
c
c     Correct the species to ensure consistency (i.e. type conservation)
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
         enddo
      endif
c
c     Correct densities by keeping fractions the same
c
      do i = is+1, ie+1
         metalfree = totalH(i) + totalHe(i)
         correctH = fh*metalfree/totalH(i)
         HI(i,j,k)  = HI(i,j,k)*correctH
         HII(i,j,k) = HII(i,j,k)*correctH
c
         correctHe = (1._RKIND - fh)*metalfree/totalHe(i)
         HeI(i,j,k)   = HeI(i,j,k)*correctHe
         HeII(i,j,k)  = HeII(i,j,k)*correctHe
         HeIII(i,j,k) = HeIII(i,j,k)*correctHe
      enddo
c
c     Correct molecular hydrogen-related fractions
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            HM   (i,j,k) = HM(i,j,k)*correctH
            H2II (i,j,k) = H2II(i,j,k)*correctH
            H2I  (i,j,k) = H2I(i,j,k)*correctH
         enddo
      endif
c
c     Do the same thing for deuterium (ignore HD) Assumes dtoh is small
c
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) 
     &             + 2._RKIND/3._RKIND*HDI(i,j,k)
            metalfree = totalH(i) + totalHe(i) + totalD
            correctD = fh*dtoh*metalfree/totalD
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
         enddo
      endif
c
c       Set the electron density
c
      do i = is+1, ie+1
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND 
     &        + HeIII(i,j,k)/2._RKIND
         if (ispecies .gt. 1) de(i,j,k) = de(i,j,k)
     &             - HM(i,j,k) + H2II(i,j,k)/2._RKIND
      enddo
c
      enddo
      enddo
c
      return
      end


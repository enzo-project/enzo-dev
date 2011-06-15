#include "fortran.def"
#include "phys_const.def"
#include "error.def"
c -----------------------------------------------------------
!   This routine scales the density fields from comoving to
!     proper densities (and back again).

      subroutine scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, factor)
c -------------------------------------------------------------------

      implicit NONE

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real     d(in,jn,kn),metal(in,jn,kn)
      real    factor

!     locals

      integer i, j, k

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
      if (ispecies .gt. 1) then
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
      if (ispecies .gt. 2) then
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
      if (imetal .eq. 1) then
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

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      real     d(in,jn,kn)
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real metal(in,jn,kn)

!     locals

      integer i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               de(i,j,k)    = max(de(i,j,k), tiny)
               HI(i,j,k)    = max(HI(i,j,k), tiny)
               HII(i,j,k)   = max(HII(i,j,k), tiny)
               HeI(i,j,k)   = max(HeI(i,j,k), tiny)
               HeII(i,j,k)  = max(HeII(i,j,k), tiny)
               HeIII(i,j,k) = max(HeIII(i,j,k), 1d-5*tiny)
            enddo
         enddo
      enddo
      if (ispecies .gt. 1) then
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
      if (ispecies .gt. 2) then
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
      if (imetal .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = min(max(metal(i,j,k), tiny),
     $                               0.9d0*d(i,j,k))
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
     &                is, ie, imax, iradtype, iradshield, ithreebody,
     &                in, jn, kn,
     &                ispecies, tgas1d, HI, HII, HeI, HeII,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k22, k24, k25, k26,
     &                k50, k51, k52, k53, k54, k55,
     &                k56, k13dd, k24shield, k25shield, k26shield,
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, itmask)
! -------------------------------------------------------------------

      implicit NONE

!     Arguments

      integer is, ie, imax, iradtype, iradshield, nratec, 
     &        in, jn, kn, ispecies, ithreebody, j, k
      real temstart, temend, tgas1d(in), dom
      double precision coolunit, tbase1
      logical itmask(imax)

!     Chemistry rates as a function of temperature

      real k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &     k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &     k54a(nratec), k55a(nratec), k56a(nratec),
     &     k13dda(nratec, 7),
     &     k24, k25, k26,
     &     avgsighp, avgsighep, avgsighe2p, piHI, piHeI

!     Density fields

      real    HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn)

!     Returned rate values

      real k1 (imax), k2 (imax), k3 (imax), k4 (imax), k5 (imax),
     &     k6 (imax), k7 (imax), k8 (imax), k9 (imax), k10(imax),
     &     k11(imax), k12(imax), k13(imax), k14(imax), k15(imax),
     &     k16(imax), k17(imax), k18(imax), k19(imax), k22(imax),
     &     k50(imax), k51(imax), k52(imax), k53(imax), k54(imax),
     &     k55(imax), k56(imax), k13dd(imax, 7),
     &     k24shield(imax), k25shield(imax), k26shield(imax)

!     1D temporaries (passed in)

      integer indixe(imax)
      real t1(imax), t2(imax), logtem(imax), tdef(imax)

!     Parameters

      double precision everg, e24, e26
      parameter(everg = ev2erg, e24 = 13.6d0, e26 = 24.6d0)

!     locals

      integer i, n1
      real factor, x, logtem0, logtem9, dlogtem, nh

!     Set log values of start and end of lookup tables

      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem = (log(temend) - log(temstart))/real(nratec-1)

      do i = is+1, ie+1
         if (itmask(i)) then
!        Compute temp-centered temperature (and log)

!        logtem(i) = log(0.5*(tgas(i)+tgasold(i)))
         logtem(i) = log(tgas1d(i))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

!        Find index into tble and precompute interpolation values

         indixe(i) = min(nratec-1,
     &        max(1,int((logtem(i)-logtem0)/dlogtem)+1))
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

      if (ispecies .gt. 1) then
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

      if (ispecies .gt. 2) then
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

!        Include approximate self-shielding factors if requested

#ifdef RADIATION
      do i = is+1, ie+1
         if (itmask(i)) then
         k24shield(i) = k24
         k25shield(i) = k25
         k26shield(i) = k26
      endif
      enddo
      if (iradshield .eq. 1) then
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

      if (iradtype .eq. 8) then
         do i = is+1, ie+1
            if (itmask(i)) then
            x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.0d-4)
            factor = 0.3908d0*(1.d0 - x**0.4092d0)**1.7592d0
            k24shield(i) = k24shield(i) + 
     &         factor*(piHI + 0.08d0*piHeI)/(e24*everg) *coolunit*tbase1
            factor = 0.0554d0*(1.d0 - x**0.4614d0)**1.6660d0
            k26shield(i) = k26shield(i) + 
     &         factor*(piHI/0.08d0 + piHeI)/(e26*everg) *coolunit*tbase1
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
            if (ispecies .gt. 1 .and. ithreebody .eq. 0) then
               do i = is+1, ie+1
                  if (itmask(i)) then
                  nh = min(HI(i,j,k)*dom, 1.0d9)
                  k13(i) = tiny
                  if (tgas1d(i) .ge. 500.d0 .and.
     &                tgas1d(i) .lt. 1.0d6) then
                     k13(i) = k13dd(i,1)-k13dd(i,2)/
     &                          (1.d0+(nh/k13dd(i,5))**k13dd(i,7))
     &                     + k13dd(i,3)-k13dd(i,4)/
     &                          (1.d0+(nh/k13dd(i,6))**k13dd(i,7))
                     k13(i) = max(10.d0**k13(i), tiny)
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

      subroutine rate_timestep(dedot, HIdot, ispecies,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, 
     &                     kdissH2I, itmask, edot, chunit)

! -------------------------------------------------------------------

      implicit NONE

!     arguments

      integer ispecies, is, ie, j, k, in, jn, kn, iradtrans, irt_honly
      real dedot(in), HIdot(in), edot(in)
      logical itmask(in)

!     Density fields

      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)

!     Radiation fields

      real    kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn)

      double precision chunit

!     Rate values

      real k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &     k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &     k11(in), k12(in), k13(in), k14(in), k15(in),
     &     k16(in), k17(in), k18(in), k19(in), k22(in),
     &     k50(in), k51(in), k52(in), k53(in), k54(in),
     &     k55(in), k56(in), 
     &     k24shield(in), k25shield(in), k26shield(in),
     &     k24, k25, k26, k27, k28, k29, k30, k31

!     locals

      integer i
      double precision H2delta, H2dmag, atten, tau

      if (ispecies .eq. 1) then

         do i = is+1, ie+1
            if (itmask(i)) then
!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k3(i)*HeI(i,j,k)*de(i,j,k)/4.d0
     &               + k5(i)*HeII(i,j,k)*de(i,j,k)/4.d0
     &               - k2(i)*HII(i,j,k)*de(i,j,k)
     &               - k4(i)*HeII(i,j,k)*de(i,j,k)/4.d0
     &               - k6(i)*HeIII(i,j,k)*de(i,j,k)/4.d0
#ifdef RADIATION
     &               +      ( k24shield(i)*HI(i,j,k)
     &               + k25shield(i)*HeII(i,j,k)/4.d0
     &               + k26shield(i)*HeI(i,j,k)/4.d0)
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
     &               -      k10(i)*H2II(i,j,k)  *HI(i,j,k)/2.d0
     &               - 2.d0*k22(i)*HI(i,j,k)**2 *HI(i,j,k)
     &               +      k2(i) *HII(i,j,k)   *de(i,j,k) 
     &               + 2.d0*k13(i)*HI(i,j,k)    *H2I(i,j,k)/2.d0
     &               +      k11(i)*HII(i,j,k)   *H2I(i,j,k)/2.d0
     &               + 2.d0*k12(i)*de(i,j,k)    *H2I(i,j,k)/2.d0
     &               +      k14(i)*HM(i,j,k)    *de(i,j,k)
     &               +      k15(i)*HM(i,j,k)    *HI(i,j,k)
     &               + 2.d0*k16(i)*HM(i,j,k)    *HII(i,j,k)
     &               + 2.d0*k18(i)*H2II(i,j,k)  *de(i,j,k)/2.d0
     &               +      k19(i)*H2II(i,j,k)  *HM(i,j,k)/2.d0
#ifdef RADIATION
     &               -      k24shield(i)*HI(i,j,k)
#endif /* RADIATION */

!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i) * HI(i,j,k)   * de(i,j,k)
     &               + k3(i) * HeI(i,j,k)  * de(i,j,k)/4.d0
     &               + k5(i) * HeII(i,j,k) * de(i,j,k)/4.d0
     &               + k8(i) * HM(i,j,k)   * HI(i,j,k)
     &               + k15(i)* HM(i,j,k)   * HI(i,j,k)
     &               + k17(i)* HM(i,j,k)   * HII(i,j,k)
     &               + k14(i)* HM(i,j,k)   * de(i,j,k)
     &               - k2(i) * HII(i,j,k)  * de(i,j,k)
     &               - k4(i) * HeII(i,j,k) * de(i,j,k)/4.d0
     &               - k6(i) * HeIII(i,j,k)* de(i,j,k)/4.d0
     &               - k7(i) * HI(i,j,k)   * de(i,j,k)
     &               - k18(i)* H2II(i,j,k) * de(i,j,k)/2.d0

#ifdef RADIATION
     &               + (k24shield(i)*HI(i,j,k)
     &               +  k25shield(i)*HeII(i,j,k)/4.e0
     &               +  k26shield(i)*HeI(i,j,k)/4.e0)
#endif /* RADIATION */

            H2delta = 
     &          HI(i,j,k)*( k22(i) * HI(i,j,k)**2.d0
     &                    - k13(i) * H2I(i,j,k)/2.d0)
            H2dmag = abs(H2delta)/(
     &          HI(i,j,k)*( k22(i) * HI(i,j,k)**2.d0
     &                    + k13(i) * H2I(i,j,k)/2.d0))
!            tau = (H2dmag/1e-5)**-1.0d0
!            tau = max(tau, 1.e-5)
!            atten = min((1.-exp(-tau))/tau,1.0)
            atten = 1.d0
            edot(i) = edot(i) + chunit * H2delta * atten
!     &       + H2I(i,j,k)*( k21(i) * HI(i,j,k)**2.0d0
!     &                    - k23(i) * H2I(i,j,k))
!H * (k22 * H^2 - k13 * H_2) + H_2 * (k21 * H^2 - k23 * H_2) */
         endif                  ! itmask
         enddo
      endif

!     Add photo-ionization rates if needed

#ifdef RADIATION
      if (iradtrans .eq. 1) then
         if (irt_honly .eq. 0) then
            do i = is+1, ie+1
               if (itmask(i)) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
     &                 + kphHeI(i,j,k) * HeI(i,j,k) / 4.d0
     &                 + kphHeII(i,j,k) *HeII(i,j,k) / 4.d0
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
     &                     in, jn, kn, is, ie, j, k, ispecies,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
     &                     dedot_prev, HIdot_prev,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII,
     &                     kdissH2I, itmask)
c -------------------------------------------------------------------

      implicit NONE

!     arguments

      integer ispecies, in, jn, kn, is, ie, j, k, iradtrans, irt_honly
      real    dtit(in), dedot_prev(in), HIdot_prev(in)
      logical itmask(in)

!     Density fields

      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)

!     Radiation fields

      real    kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn)

!     Rate values

      real k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &     k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &     k11(in), k12(in), k13(in), k14(in), k15(in),
     &     k16(in), k17(in), k18(in), k19(in), k22(in),
     &     k50(in), k51(in), k52(in), k53(in), k54(in),
     &     k55(in), k56(in), 
     &     k24shield(in), k25shield(in), k26shield(in),
     &     k24, k25, k26, k27, k28, k29, k30, k31

!     temporaries (passed in)

      real HIp(in), HIIp(in), HeIp(in), HeIIp(in), HeIIIp(in),
     &     HMp(in), H2Ip(in), H2IIp(in), dep(in),
     &     DIp(in), DIIp(in), HDIp(in)

!     locals

      integer i
      real scoef, acoef

!   A) the 6-species integrator
!      
      if (ispecies .eq. 1) then

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI

            scoef  = k2(i)*HII(i,j,k)*de(i,j,k)
            acoef  = k1(i)*de(i,j,k)
#ifdef RADIATION
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
#endif /* RADIATION */
            HIp(i)  = (scoef*dtit(i) + HI(i,j,k))/(1.d0 + acoef*dtit(i))
            if (HIp(i) .ne. HIp(i)) then
               write(*,*) 'HUGE HIp! :: ', i, j, k, HIp(i), HI(i,j,k),
     $              HII(i,j,k), de(i,j,k), kphHI(i,j,k), 
     $              scoef, acoef, dtit(i)
               ERROR_MESSAGE
            endif

!        2) HII
c 
            scoef  = k1(i)*HIp(i)*de(i,j,k)
#ifdef RADIATION
     &                   + k24shield(i)*HIp(i)
            if (iradtrans .eq. 1) scoef = scoef + kphHI(i,j,k)*HIp(i)
#endif /* RADIATION */
            acoef  = k2(i)*de (i,j,k)
            HIIp(i) = (scoef*dtit(i) + HII(i,j,k))/(1.d0 +acoef*dtit(i))
!
            if (HIIp(i) .le. 0.0) then  !#####
               write(*,*) 'negative HIIp! :: ', i, j, k, HIIp(i), 
     $              scoef, dtit(i), HII(i,j,k), acoef,
     $              k2(i), de(i,j,k),
     $              kphHI(i,j,k), HIp(i),
     $              k24shield(i)
            endif

!        3) Electron density

            scoef = 0.d0
#ifdef RADIATION
     &                 + k24shield(i)*HI(i,j,k)
     &                 + k25shield(i)*HeII(i,j,k)/4.d0
     &                 + k26shield(i)*HeI(i,j,k)/4.d0
            if (iradtrans .eq. 1 .and. irt_honly .eq. 0) 
     &           scoef = scoef + kphHI(i,j,k)   * HI(i,j,k)
     &           + kphHeI(i,j,k)  * HeI(i,j,k)/4.d0
     &           + kphHeII(i,j,k) * HeII(i,j,k)/4.d0
            if (iradtrans .eq. 1 .and. irt_honly .eq. 1) 
     &           scoef = scoef + kphHI(i,j,k)   * HI(i,j,k)
#endif /* RADIATION */
            acoef = -(k1(i)*HI(i,j,k)      - k2(i)*HII(i,j,k)
     &              + k3(i)*HeI(i,j,k)/4.d0  - k6(i)*HeIII(i,j,k)/4.d0
     &              + k5(i)*HeII(i,j,k)/4.d0 - k4(i)*HeII(i,j,k)/4.d0)
            dep(i)   = (scoef*dtit(i) + de(i,j,k))
     &                     / (1.d0 + acoef*dtit(i))

         endif                  ! itmask
         enddo

      endif                     ! (ispecies .eq. 1)

!  --- (B) Do helium chemistry in any case: (for all ispecies values) ---

      do i = is+1, ie+1
         if (itmask(i)) then

!        4) HeI

         scoef  = k4(i)*HeII(i,j,k)*de(i,j,k)
         acoef  = k3(i)*de(i,j,k)
#ifdef RADIATION
     &                + k26shield(i)
         if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $        acoef = acoef + kphHeI(i,j,k)
#endif /* RADIATION */
         HeIp(i)   = ( scoef*dtit(i) + HeI(i,j,k) ) 
     &              / ( 1.d0 + acoef*dtit(i) )

!        5) HeII

         scoef  = k3(i)*HeIp(i)*de(i,j,k)
     &          + k6(i)*HeIII(i,j,k)*de(i,j,k)
#ifdef RADIATION
     &          + k26shield(i)*HeIp(i)
         if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $        scoef = scoef + kphHeI(i,j,k)*HeIp(i)
#endif /* RADIATION */
         acoef  = k4(i)*de(i,j,k) + k5(i)*de(i,j,k)
#ifdef RADIATION
     &          + k25shield(i)
         if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $        acoef = acoef + kphHeII(i,j,k)
#endif /* RADIATION */
         HeIIp(i)  = ( scoef*dtit(i) + HeII(i,j,k) )
     &              / ( 1.d0 + acoef*dtit(i) )

!       6) HeIII

         scoef   = k5(i)*HeIIp(i)*de(i,j,k)
#ifdef RADIATION
     &           + k25shield(i)*HeIIp(i)
         if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $        scoef = scoef + kphHeII(i,j,k) * HeIIp(i)
#endif /* RADIATION */
         acoef   = k6(i)*de(i,j,k)
         HeIIIp(i)  = ( scoef*dtit(i) + HeIII(i,j,k) )
     &                / ( 1.d0 + acoef*dtit(i) )

      endif                     ! itmask
      enddo

c --- (C) Now do extra 3-species for molecular hydrogen ---

      if (ispecies .gt. 1) then

!        First, do HI/HII with molecular hydrogen terms

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI
!     
            scoef  =      k2(i) * HII(i,j,k) * de(i,j,k) 
     &             + 2.d0*k13(i)* HI(i,j,k)  * H2I(i,j,k)/2.d0
     &             +      k11(i)* HII(i,j,k) * H2I(i,j,k)/2.d0
     &             + 2.d0*k12(i)* de(i,j,k)  * H2I(i,j,k)/2.d0
     &             +      k14(i)* HM(i,j,k)  * de(i,j,k)
     &             +      k15(i)* HM(i,j,k)  * HI(i,j,k)
     &             + 2.d0*k16(i)* HM(i,j,k)  * HII(i,j,k)
     &             + 2.d0*k18(i)* H2II(i,j,k)* de(i,j,k)/2.d0
     &             +      k19(i)* H2II(i,j,k)* HM(i,j,k)/2.d0
#ifdef RADIATION
     &             + 2.d0*k31   * H2I(i,j,k)/2.d0
            if (iradtrans.eq.1)
     &           scoef = scoef + 2.d0*kdissH2I(i,j,k)* H2I(i,j,k)/2.d0
#endif /* RADIATION */
            acoef  =      k1(i) * de(i,j,k)
     &             +      k7(i) * de(i,j,k)  
     &             +      k8(i) * HM(i,j,k)
     &             +      k9(i) * HII(i,j,k)
     &             +      k10(i)* H2II(i,j,k)/2.d0
     &             + 2.d0*k22(i)* HI(i,j,k)**2
#ifdef RADIATION
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
#endif /* RADIATION */
            HIp(i)  = ( scoef*dtit(i) + HI(i,j,k) ) / 
     &                      ( 1. + acoef*dtit(i) )
            if (HIp(i) .ne. HIp(i)) then
               write(*,*) 'HUGE HIp! :: ', i, j, k, HIp(i), HI(i,j,k),
     $              HII(i,j,k), de(i,j,k), H2I(i,j,k), kdissH2I(i,j,k),
     $              kphHI(i,j,k)
            endif

!          2) HII

            scoef  =    k1(i)  * HI(i,j,k) * de(i,j,k)
     &             +    k10(i) * H2II(i,j,k)*HI(i,j,k)/2.d0
#ifdef RADIATION
     &             + k24shield(i)*HI(i,j,k)
            if (iradtrans .eq. 1) scoef = scoef + kphHI(i,j,k)*HI(i,j,k)
#endif /* RADIATION */
            acoef  =    k2(i)  * de(i,j,k)
     &             +    k9(i)  * HI(i,j,k)
     &             +    k11(i) * H2I(i,j,k)/2.d0
     &             +    k16(i) * HM(i,j,k)
     &             +    k17(i) * HM(i,j,k)
            HIIp(i)   = ( scoef*dtit(i) + HII(i,j,k) )
     &                      / ( 1.d0 + acoef*dtit(i) )
!     
!          3) electrons:

            scoef =   k8(i) * HM(i,j,k) * HI(i,j,k)
     &             +  k15(i)* HM(i,j,k) * HI(i,j,k)
     &             +  k17(i)* HM(i,j,k) * HII(i,j,k)
!                  
#ifdef RADIATION
     &             + k24shield(i)*HIp(i)
     &             + k25shield(i)*HeIIp(i)/4.d0
     &             + k26shield(i)*HeIp(i)/4.d0
            if (iradtrans .eq. 1 .and. irt_honly.eq.0) 
     &           scoef = scoef + kphHI(i,j,k)   * HIp(i)
     &           + kphHeI(i,j,k)  * HeIp(i)/4.d0
     &           + kphHeII(i,j,k) * HeIIp(i)/4.d0
            if (iradtrans .eq. 1 .and. irt_honly.eq.1) 
     &           scoef = scoef + kphHI(i,j,k)   * HIp(i)
#endif /* RADIATION */
            acoef = - (k1(i) *HI(i,j,k)    - k2(i)*HII(i,j,k)
     &              +  k3(i) *HeI(i,j,k)/4.d0 - k6(i)*HeIII(i,j,k)/4.d0
     &              +  k5(i) *HeII(i,j,k)/4.d0- k4(i)*HeII(i,j,k)/4.d0
     &              +  k14(i)*HM(i,j,k)
     &              -  k7(i) *HI(i,j,k)
     &              -  k18(i)*H2II(i,j,k)/2.d0)
            dep(i)  = ( scoef*dtit(i) + de(i,j,k) )
     &                / ( 1.d0 + acoef*dtit(i) )

!           7) H2

            scoef = 2.d0*(k8(i)  * HM(i,j,k)   * HI(i,j,k)
     &            +       k10(i) * H2II(i,j,k) * HI(i,j,k)/2.d0
     &            +       k19(i) * H2II(i,j,k) * HM(i,j,k)/2.d0
     &            +       k22(i) * HI(i,j,k)**3)
            acoef = ( k13(i)*HI(i,j,k) + k11(i)*HII(i,j,k)
     &              + k12(i)*de(i,j,k) )
#ifdef RADIATION
     &              + k29 + k31
            if (iradtrans.eq.1) acoef = acoef + kdissH2I(i,j,k)
#endif /* RADIATION */

            H2Ip(i) = ( scoef*dtit(i) + H2I(i,j,k) )
     &                / ( 1.d0 + acoef*dtit(i) )

!           8) H-

            HMp(i) = ( k7(i)*HIp(i)*dep(i) )
     &             / ( (k8(i)+k15(i))*HIp(i)
     &             + ( k16(i)+k17(i))*HIIp(i)+k14(i)*dep(i)
#ifdef RADIATION
     &             + k27
#endif /* RADIATION */
     &                 )

!           9) H2+

            H2IIp(i) = 2.d0*( k9 (i)*HIp(i)*HIIp(i)
     &                    +   k11(i)*H2Ip(i)/2.d0*HIIp(i)
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
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
!     
!         1) DI
!     
            scoef =   (       k2(i) * DII(i,j,k) * de(i,j,k)
     &                 +      k51(i)* DII(i,j,k) * HI(i,j,k)
     &                 + 2.d0*k55(i)* HDI(i,j,k) * HI(i,j,k)/3.d0
     &                 )
            acoef  =    k1(i) * de(i,j,k)
     &             +    k50(i) * HII(i,j,k)
     &             +    k54(i) * H2I(i,j,k)/2.d0
     &             +    k56(i) * HM(i,j,k)
#ifdef RADIATION
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
#endif /* RADIATION */
            DIp(i)    = ( scoef*dtit(i) + DI(i,j,k) ) / 
     &                  ( 1.d0 + acoef*dtit(i) )

!         2) DII
c 
            scoef =   (   k1(i)  * DI(i,j,k) * de(i,j,k)
     &            +       k50(i) * HII(i,j,k)* DI(i,j,k)
     &            +  2.d0*k53(i) * HII(i,j,k)* HDI(i,j,k)/3.d0
     &            )
#ifdef RADIATION
     &            + k24shield(i)*DI(i,j,k)
            if (iradtrans .eq. 1) scoef = scoef + kphHI(i,j,k)*DI(i,j,k)
#endif /* RADIATION */
            acoef =    k2(i)  * de(i,j,k)
     &            +    k51(i) * HI(i,j,k)
     &            +    k52(i) * H2I(i,j,k)/2.d0

            DIIp(i)   = ( scoef*dtit(i) + DII(i,j,k) )
     &                 / ( 1.d0 + acoef*dtit(i) )

!          3) HDI
c 
            scoef = 3.d0*(k52(i) * DII(i,j,k)* H2I(i,j,k)/2.d0/2.d0
     &             +      k54(i) * DI(i,j,k) * H2I(i,j,k)/2.d0/2.d0
     &             + 2.d0*k56(i) * DI(i,j,k) * HM(i,j,k)/2.d0
     &                 )
            acoef  =    k53(i) * HII(i,j,k)
     &             +    k55(i) * HI(i,j,k)

            HDIp(i)   = ( scoef*dtit(i) + HDI(i,j,k) )
     &                 / ( 1.d0 + acoef*dtit(i) )

         endif                  ! itmask
         enddo
      endif

!   --- (E) Set densities from 1D temps to 3D fields ---

      do i = is+1, ie+1
         if (itmask(i)) then
         HIdot_prev(i) = abs(HI(i,j,k)-HIp(i))/max(dtit(i),tiny)
         HI(i,j,k)    = max(HIp(i), tiny)
         HII(i,j,k)   = max(HIIp(i), tiny)
         HeI(i,j,k)   = max(HeIp(i), tiny)
         HeII(i,j,k)  = max(HeIIp(i), tiny)
         HeIII(i,j,k) = max(HeIIIp(i), 1d-5*tiny)

!        de(i,j,k)    = dep(i)

!        Use charge conservation to determine electron fraction

         dedot_prev(i) = de(i,j,k)
         de(i,j,k) = HII(i,j,k) + HeII(i,j,k)/4.d0 + HeIII(i,j,k)/2.d0
         if (ispecies .gt. 1) 
     &      de(i,j,k) = de(i,j,k) - HM(i,j,k) + H2II(i,j,k)/2.d0
         dedot_prev(i) = abs(de(i,j,k)-dedot_prev(i))/
     &                         max(dtit(i),tiny)

         if (ispecies .gt. 1) then
            HM(i,j,k)    = max(HMp(i), tiny)
            H2I(i,j,k)   = max(H2Ip(i), tiny)
            H2II(i,j,k)  = max(H2IIp(i), tiny)
         endif

         if (ispecies .gt. 2) then
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
     &     HM, H2I, H2II, DI, DII, HDI, metal, d,
     &     is, ie, js, je, ks, ke, in, jn, kn, imax, 
     &     ispecies, imetal, fh, dtoh)
! -------------------------------------------------------------------

      implicit NONE

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn), metal(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)  
      real    fh, dtoh

!     Parameters

      integer imax

!     locals

      integer i, j, k
      real totalH(imax), totalHe(imax), 
     &     totalD, correctH, correctHe, correctD
      real metalfree(imax)

!     Loop over all zones

      do k = ks+1, ke+1
      do j = js+1, je+1

!     Compute total densities of H and He
!         (ensure non-negativity)

      if (imetal .eq. 1) then
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

      if (ispecies .gt. 1) then
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

         correctHe = (1.d0 - fh)*metalfree(i)/totalHe(i)
         HeI(i,j,k)   = HeI(i,j,k)*correctHe
         HeII(i,j,k)  = HeII(i,j,k)*correctHe
         HeIII(i,j,k) = HeIII(i,j,k)*correctHe

!     Correct molecular hydrogen-related fractions

         if (ispecies .gt. 1) then
            HM   (i,j,k) = HM(i,j,k)*correctH
            H2II (i,j,k) = H2II(i,j,k)*correctH
            H2I  (i,j,k) = H2I(i,j,k)*correctH
         endif
      enddo

!     Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) + 2.d0/3.d0*HDI(i,j,k)
            correctD = fh*dtoh*metalfree(i)/totalD
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
         enddo
      endif

!       Set the electron density

      do i = is+1, ie+1
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4.d0 + HeIII(i,j,k)/2.d0
         if (ispecies .gt. 1) de(i,j,k) = de(i,j,k)
     &             - HM(i,j,k) + H2II(i,j,k)/2.d0
      enddo

      enddo  ! end loop over j
      enddo  ! end loop over k

      return
      end

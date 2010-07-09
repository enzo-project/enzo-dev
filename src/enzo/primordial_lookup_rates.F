#include "fortran.def"
c
c -----------------------------------------------------------
c This routine uses the temperature to look up the chemical
c   rates which are tabulated in a log table as a function
c   of temperature.
c
      subroutine primordial_lookup_rates(temstart, temend, nratec, j, k,
     $                is, ie, ijk, iradtype, iradshield, in, jn, kn,
     $                ispecies, tgas1d, tgasold, d, HI, HII, HeI, HeII,
     $                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     $                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     $                k17a, k18a, k19a, k21a, k22a, k23a,
     $                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     $                avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     $                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     $                k11, k12, k13, k14, k15, k16, k17, k18,
     $                k19, k21, k22, k23, k24, k25, k26,
     $                k50, k51, k52, k53, k54, k55,
     $                k56, k13dd, k24shield, k25shield, k26shield,
     $                t1, t2, tdef, logtem, indixe, 
     $                dom, coolunit, tbase1, itmask,
     $                eqmask, threebody )
c -------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer is, ie, ijk, iradtype, iradshield, nratec, 
     &        in, jn, kn, ispecies, j, k, threebody
      real temstart, temend, tgas1d(in), tgasold(in), dom
      double precision coolunit, tbase1
      logical itmask(in), eqmask(in)
c
c     Chemistry rates as a function of temperature
c
      real k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k21a(nratec), 
     &     k22a(nratec), k23a(nratec), k50a(nratec), k51a(nratec),
     &     k52a(nratec), k53a(nratec), k54a(nratec), k55a(nratec), 
     &     k56a(nratec), k13dda(nratec, 7),
     &     k24, k25, k26,
     &     avgsighp, avgsighep, avgsighe2p, piHI, piHeI
c
c     Density fields
c
c      real    d(in,jn,kn), HI(in,jn,kn),   HII(in,jn,kn),
c     &        HeI(in,jn,kn), HeII(in,jn,kn)
       real d(in,jn,kn)
       real*8 HI(ijk), HII(ijk), HeI(ijk), HeII(ijk)
c
c     Returned rate values
c
      real k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k21(ijk),
     &     k22(ijk), k23(ijk), k50(ijk), k51(ijk), k52(ijk), 
     &     k53(ijk), k54(ijk), k55(ijk), k56(ijk), k13dd(ijk, 7),
     &     k24shield(ijk), k25shield(ijk), k26shield(ijk)
c
c     1D temporaries (passed in)
c 
      integer indixe(ijk)
      real t1(ijk), t2(ijk), logtem(ijk), tdef(ijk)
c     &     k13dd_dot1(ijk,7), k13dd_dot2(ijk,7)
c
c     Parameters
c
      double precision everg, e24, e26, temp
      parameter(everg = 1.60184d-12, e24 = 13.6d0, e26 = 24.6d0)
c
c     locals
c
      integer i, n1
      real factor, x, logtem0, logtem9, dlogtem, nh
c     &     k13dd_dot(ijk, 7), dd1, dd2, dd3, dd4, dd5, dd6, dd7
c
c     Set log values of start and end of lookup tables
c
      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem = (log(temend) - log(temstart))/real(nratec-1)
c
      do i = is+1, ie+1
      if (itmask(i).eqv..true.) then
c         if(tgas1d(i).gt.11000.and.d(i,j,k)*dom.gt.1e17)eqmask(i)=.true.
c
c        Compute time-centered temperature (and log)
c
         logtem(i) = log(tgas1d(i))
c         logtem(i) = log(0.5*(tgas1d(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)
c
c        Find index into table and precompute interpolation values
c
         indixe(i) = min(nratec-1,
     &                max(1,int((logtem(i)-logtem0)/dlogtem)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = t2(i) - t1(i)
c
c        Do linear table lookup (in log temperature)
c
#define LEFT +1
#define RIGHT -0
         k1(i) = k1a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k1a(indixe(i)LEFT) -k1a(indixe(i)RIGHT))/tdef(i)
         k2(i) = k2a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k2a(indixe(i)LEFT) -k2a(indixe(i)RIGHT))/tdef(i)
c         if(d(i,j,k)*dom.gt.1e12)k2(i)=tiny
         k3(i) = k3a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k3a(indixe(i)LEFT) -k3a(indixe(i)RIGHT))/tdef(i)
         k4(i) = k4a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k4a(indixe(i)LEFT) -k4a(indixe(i)RIGHT))/tdef(i)
         k5(i) = k5a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k5a(indixe(i)LEFT) -k5a(indixe(i)RIGHT))/tdef(i)
         k6(i) = k6a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k6a(indixe(i)LEFT) -k6a(indixe(i)RIGHT))/tdef(i)
c
      endif
      enddo
c
c     Look-up for 9-species model
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
      if (itmask(i).eqv..true.) then
            k7(i) = k7a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k7a(indixe(i)LEFT) -k7a(indixe(i)RIGHT))/tdef(i)
            k8(i) = k8a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k8a(indixe(i)LEFT) -k8a(indixe(i)RIGHT))/tdef(i)
            k9(i) = k9a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k9a(indixe(i)LEFT) -k9a(indixe(i)RIGHT))/tdef(i)
            k10(i) = k10a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k10a(indixe(i)LEFT) -k10a(indixe(i)RIGHT))/tdef(i)
            k11(i) = k11a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k11a(indixe(i)LEFT) -k11a(indixe(i)RIGHT))/tdef(i)
            k12(i) = k12a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k12a(indixe(i)LEFT) -k12a(indixe(i)RIGHT))/tdef(i)
            k13(i) = k13a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k13a(indixe(i)LEFT) -k13a(indixe(i)RIGHT))/tdef(i)
            k14(i) = k14a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k14a(indixe(i)LEFT) -k14a(indixe(i)RIGHT))/tdef(i)
            k15(i) = k15a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k15a(indixe(i)LEFT) -k15a(indixe(i)RIGHT))/tdef(i)
            k16(i) = k16a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k16a(indixe(i)LEFT) -k16a(indixe(i)RIGHT))/tdef(i)
            k17(i) = k17a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k17a(indixe(i)LEFT) -k17a(indixe(i)RIGHT))/tdef(i)
            k18(i) = k18a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k18a(indixe(i)LEFT) -k18a(indixe(i)RIGHT))/tdef(i)
            k19(i) = k19a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k19a(indixe(i)LEFT) -k19a(indixe(i)RIGHT))/tdef(i)
            k21(i) = k21a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k21a(indixe(i)LEFT) -k21a(indixe(i)RIGHT))/tdef(i)
            k22(i) = k22a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k22a(indixe(i)LEFT) -k22a(indixe(i)RIGHT))/tdef(i)
            k23(i) = k23a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k23a(indixe(i)LEFT) -k23a(indixe(i)RIGHT))/tdef(i)
         endif
         enddo
c     
         do n1 = 1, 7
            do i = is+1, ie+1
            if (itmask(i).eqv..true.) then
               k13dd(i,n1) = k13dda(indixe(i),n1) + (logtem(i) - t1(i))
     &             *(k13dda(indixe(i)LEFT,n1) - 
     &               k13dda(indixe(i)RIGHT,n1) )/tdef(i)
            endif
            enddo
         enddo
      endif
c
c     Look-up for 12-species model
c
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
            k50(i) = k50a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k50a(indixe(i)LEFT) -k50a(indixe(i)RIGHT))/tdef(i)
            k51(i) = k51a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k51a(indixe(i)LEFT) -k51a(indixe(i)RIGHT))/tdef(i)
            k52(i) = k52a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k52a(indixe(i)LEFT) -k52a(indixe(i)RIGHT))/tdef(i)
            k53(i) = k53a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k53a(indixe(i)LEFT) -k53a(indixe(i)RIGHT))/tdef(i)
            k54(i) = k54a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k54a(indixe(i)LEFT) -k54a(indixe(i)RIGHT))/tdef(i)
            k55(i) = k55a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k55a(indixe(i)LEFT) -k55a(indixe(i)RIGHT))/tdef(i)
            k56(i) = k56a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k56a(indixe(i)LEFT) -k56a(indixe(i)RIGHT))/tdef(i)
         endif
         enddo
      endif
c
c        Include approximate self-shielding factors if requested
c
#ifdef RADIATION
      do i = is+1, ie+1
      if (itmask(i).eqv..true.) then
         k24shield(i) = k24
         k25shield(i) = k25
         k26shield(i) = k26
      endif
      enddo
      if (iradshield .eq. 1) then
         stop
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
            k24shield(i) = k24shield(i)*exp(-HI(i)*avgsighp*dom)
            k25shield(i) = k25shield(i)*exp(-HeII(i)*avgsighe2p*dom)
            k26shield(i) = k26shield(i)*exp(-HeI(i)*avgsighep*dom)
         endif
         enddo
      endif
c
c        If using a high-energy radiation field, then account for
c          effects of secondary elections (Shull * Steenberg 1985)
c          (see calc_rate.src)
c
      if (iradtype .eq. 8) then
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
            x = max(HII(i)/(HI(i)+HII(i)), 1.0e-4)
            factor = 0.3908*(1.0 - x**0.4092)**1.7592
            k24shield(i) = k24shield(i) + 
     &         factor*(piHI + 0.08*piHeI)/(e24*everg) * coolunit*tbase1
            factor = 0.0554*(1.0 - x**0.4614)**1.6660
            k26shield(i) = k26shield(i) + 
     &         factor*(piHI/0.08 + piHeI)/(e26*everg) * coolunit*tbase1
         endif
         enddo
      endif
c
#endif /* RADIATION */
c
c           If using H2, and using the density-dependent collisional
c             H2 dissociation rate, then replace the the density-independant
c                k13 rate with the new one.
c         May/00: there appears to be a problem with the density-dependent
c             collisional rates.  Currently turned off until further notice.
c
#define USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
#ifdef USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
            if ((ispecies .gt. 1).and.(threebody.eq.0)) then
               do i = is+1, ie+1
               if (itmask(i).eqv..true.) then
                  nh = min(HI(i)*dom, 1.0e9)
                  k13(i) = tiny
                  if (tgas1d(i) .ge. 500.0 .and.
     &                tgas1d(i) .lt. 1.0e6) then
                     k13(i) = k13dd(i,1)-k13dd(i,2)/
     &                          (1.0+(nh/k13dd(i,5))**k13dd(i,7))
     &                      + k13dd(i,3)-k13dd(i,4)/
     &                          (1.0+(nh/k13dd(i,6))**k13dd(i,7))
                     k13(i) = max(10.0**k13(i), tiny)
                  endif
               endif
               enddo
            endif
#endif /*  USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE */
      return
      end



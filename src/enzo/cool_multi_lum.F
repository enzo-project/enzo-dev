#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c////////////////////  SUBROUTINE COOL_MULTI_LUM  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine cool_multi_lum(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                  lum, nlum, dx,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, ispecies, imetal, imcool, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht,
     &                dt, aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                gpldla, gphdla, hdltea, hdlowa, 
     &                metala, n_xe, xe_start, xe_end,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, photogamma)
c
c  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
c
c  written by: Yu Zhang, Peter Anninos and Tom Abel
c  date:       
c  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
c  modified2: October, 1996 by GB; moved to AMR
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
      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke, nratec, imethod,
     &        idual, iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        iradtype, nfreq, imetalregen, iradshield, 
     &        iradtrans, n_xe, nlum, imcool
      R_PREC    dt, aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        eta1, eta2, gamma, xe_start, xe_end, dx
      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        lum(in,jn,kn,nlum)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn)
      R_PREC    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec),
     $        metala(nratec, n_xe)
      R_PREC    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      R_PREC    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p
      R_PREC photogamma(in,jn,kn)
c
c  Parameters
c
      INTG_PREC ijk
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)
c
c  Locals
c
      INTG_PREC i, j, k, n
      R_PREC comp1, comp2, energy
      real*8 dom, coolunit, dbase1, tbase1, xbase1, mh, ulum
      parameter (mh = mass_h)
c
c  Slice locals
c 
      INTG_PREC indixe(ijk)
      R_PREC t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), p2d(ijk),
     &     dtit(ijk), ttot(ijk), tgas(ijk), tgasold(ijk)
      real*8 edot(in,nlum)
c
c  Cooling/heating slice locals
c
      real*8 ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &     ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &     reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk),
     &     brem(ijk)
      R_PREC hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), rotl(ijk),
     &     gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk), metalc(ijk)

!  Iteration mask for multi_cool
      LOGIC_PREC itmask(ijk)
                
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Convert densities from comoving to 'proper'
c
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)/aye**3
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
            if (imetal .ge. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)/aye**3
               enddo
            endif
         enddo
      enddo

      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      
c     wrong... Conversion from cooling rate to luminosity in units of 1e40 erg/s
c      ulum = 1e-40 * (dble(uxyz)*dble(dx))**3 * coolunit
      ulum = coolunit

c 
c     Loop over slices (in the k-direction)
c
      do k = ks+1, ke+1
       do j = js+1, je+1

          do i = is+1, ie+1
             itmask(i) = .true.
          enddo

c
c        Compute the cooling rate
c
         call cool1d_sep(     
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, imetal, imcool, idim,
     &                is, ie, j, k, ih2co, ipiht, 1_IKIND,
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
     &                indixe, t1, t2, logtem, tdef, edot, nlum,
     &                tgas, tgasold, p2d,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, photogamma, itmask)
c
c        Compute the cooling luminosity on the slice.  Multiply by
c        density because edot is specific energy (erg/s), not per unit
c        volume.
c
         do n = 1, nlum
            do i = is+1, ie+1
               lum(i,j,k,n) = -edot(i,n)*ulum!*d(i,j,k)
            enddo
         enddo
c
       enddo
      enddo
c
c     Convert densities back to comoving from 'proper'
c
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*aye**3
               de(i,j,k)    = de(i,j,k)*aye**3
               HI(i,j,k)    = HI(i,j,k)*aye**3
               HII(i,j,k)   = HII(i,j,k)*aye**3
               HeI(i,j,k)   = HeI(i,j,k)*aye**3
               HeII(i,j,k)  = HeII(i,j,k)*aye**3
               HeIII(i,j,k) = HeIII(i,j,k)*aye**3
            enddo
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*aye**3
                  H2I(i,j,k)  = H2I(i,j,k)*aye**3
                  H2II(i,j,k) = H2II(i,j,k)*aye**3
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*aye**3
                  DII(i,j,k) = DII(i,j,k)*aye**3
                  HDI(i,j,k) = HDI(i,j,k)*aye**3
               enddo
            endif
            if (imetal .ge. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*aye**3
               enddo
            endif
         enddo
      enddo
c
      return
      end


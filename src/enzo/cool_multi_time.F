#include "fortran.def"
!=======================================================================
!////////////////////  SUBROUTINE COOL_MULTI_TIME  \\\\\\\\\\\\\\\\\\\\\

      subroutine cool_multi_time(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                  cooltime, coolonly,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, ispecies, imetal, imcool, idust, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht, igammah,
     &                dt, aye, redshift, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma, z_solar,
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
     &                iradtrans, photogamma,
     &                ih2optical, iciecool, ciecoa, 
     &                icmbTfloor, iClHeat,
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating)

!  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
!  modified2: October, 1996 by GB; moved to AMR
!  modified3: February, 2003 by Robert Harkness; iteration mask
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
#include "fortran_types.def"

!  Arguments

      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke, nratec, imethod,
     &        idual, iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        iradtype, nfreq, imetalregen, iradshield, iradtrans,
     &        n_xe, imcool, idust, igammah, ih2optical, iciecool,
     &        coolonly
      R_PREC    dt, aye, temstart, temend, redshift,
     &        utem, uxyz, uaye, urho, utim,
     &        eta1, eta2, gamma, xe_start, xe_end, z_solar
      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        cooltime(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      R_PREC    metal(in,jn,kn)
      R_PREC    photogamma(in,jn,kn)
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

!  Cloudy cooling data

      INTG_PREC icmbTfloor, iClHeat, clGridRank, clDataSize
      INTG_PREC clGridDim(5)
      R_PREC clEleFra
      R_PREC clPar1(clGridDim(1)), clPar2(clGridDim(2)), 
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)), 
     &     clPar5(clGridDim(5))
      R_PREC clCooling(clDataSize), clHeating(clDataSize)

!  Parameters

      INTG_PREC ijk
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)

!  Locals

      INTG_PREC i, j, k, n
      R_PREC comp1, comp2, energy

!  Slice locals
 
      INTG_PREC indixe(ijk)
      R_PREC t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), p2d(ijk),
     &     dtit(ijk), ttot(ijk), tgas(ijk), tgasold(ijk),
     &     tdust(ijk), metallicity(ijk), rhoH(ijk)
      real*8 edot(ijk)

!  Cooling/heating slice locals

      real*8 ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &     ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &     reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk),
     &     brem(ijk), cieco(ijk)
      R_PREC hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), 
     &     rotl(ijk), gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk)
          
!  Iteration mask for multi_cool

      LOGIC_PREC itmask(ijk)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================


!     Convert densities from comoving to 'proper'

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
            if (ispecies > 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)/aye**3
                  H2I(i,j,k)  = H2I(i,j,k)/aye**3
                  H2II(i,j,k) = H2II(i,j,k)/aye**3
               enddo
            endif
            if (ispecies > 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)/aye**3
                  DII(i,j,k) = DII(i,j,k)/aye**3
                  HDI(i,j,k) = HDI(i,j,k)/aye**3
               enddo
            endif
            if (imetal == 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)/aye**3
               enddo
            endif
         enddo
      enddo
 
!     Loop over slices (in the k-direction)

      do k = ks+1, ke+1
       do j = js+1, je+1
         do i = is+1, ie+1
            itmask(i) = .true.
         end do

!        Compute the cooling rate
         call cool1d_multi(     
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, imetal, imcool, idust, idim,
     &                is, ie, j, k, ih2co, ipiht, 1_IKIND, igammah,
     &                aye, redshift, temstart, temend, z_solar,
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

!        Compute the cooling time on the slice
!          (the gamma used here is the same as used to calculate the pressure
!           in cool1d_multi)
         if (coolonly .eq. 1) then
            do i = is+1, ie+1
               energy = max(p2d(i)/(gamma-1._RKIND), tiny)   
               if (edot(i) .lt. 0._RKIND) then
                  cooltime(i,j,k) = -energy/edot(i)
               else
                  cooltime(i,j,k) = huge;
               endif
            enddo
         else
            do i = is+1, ie+1
               energy = max(p2d(i)/(gamma-1._RKIND), tiny)   
               cooltime(i,j,k) = abs(energy/edot(i))
            enddo
         endif

       enddo
      enddo

!     Convert densities back to comoving from 'proper'
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
            if (ispecies > 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*aye**3
                  H2I(i,j,k)  = H2I(i,j,k)*aye**3
                  H2II(i,j,k) = H2II(i,j,k)*aye**3
               enddo
            endif
            if (ispecies > 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*aye**3
                  DII(i,j,k) = DII(i,j,k)*aye**3
                  HDI(i,j,k) = HDI(i,j,k)*aye**3
               enddo
            endif
            if (imetal == 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*aye**3
               enddo
            endif
         enddo
      enddo

      return
      end


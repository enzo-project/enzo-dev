#include "fortran.def"
c=======================================================================
c//////////////////////  SUBROUTINE COOL1D_MULTI  \\\\\\\\\\\\\\\\\\\\\\
c
      subroutine primordial_cool(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, idim,
     &                is, ie, j, k, ih2co, ipiht, iter,
     &                aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, 
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, 
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow, 
     &                hdcoola, hdcool, ciecoa, cieco, 
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edotplus,
     &                edotminus, tgas, tgasold, p2d,
     &                inutot, iradtype, nfreq, 
     &                iradshield, avgsighp, avgsighep, 
     &                avgsighe2p, itmask, iciecool,ih2optical,ifedtgas,
     &                gamma2, omask)!, 
!     &                sobnstart, sobnend, sobtstart, sobtend, 
!     &                sobnn, sobnt, sobratio, neff )
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
c
c  Arguments
c
      integer in, jn, kn, is, ie, j, k, nratec, imethod, idim,
     &        idual, iexpand, ih2co, ipiht, ispecies,
     &        nfreq, iradshield, iradtype, iciecool,
     &        ih2optical, ifedtgas, ijk
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)

      real    aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        eta1, eta2, gamma
      real    d(in,jn,kn), e(in,jn,kn)
      real    u(ijk),    v(ijk),     w(ijk)
      real*8
     &       de(ijk),   HI(ijk),   HII(ijk),
     &      HeI(ijk), HeII(ijk), HeIII(ijk),
     &       ge(ijk),
     &        HM(ijk),  H2I(ijk), H2II(ijk),
     &        DI(ijk),  DII(ijk), HDI(ijk)
      logical omask(ijk)
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
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p,
     &        ciecont, mciecont, medot, gamma2(in)
c
c  Parameters
c
      double precision mh
      real    ZSOLAR
      parameter (mh = 1.67d-24, ZSOLAR = 0.02041)
c
c  Locals
c
      integer i, m, j1, iter, iradfield, mciei
      real dom, qq, vibl, logtem0, logtem9, dlogtem, energy, zr
      real dt2, ttmin, comp1, comp2, scoef, acoef, HIdot,
     &     hdlte1, hdlow1, x, nH2, nother, fudge, fH2,
     &     gphdl1, factor, ciefudge, tau, comp_edot, comp_xray_edot
      double precision coolunit, dbase1, tbase1, xbase1, 
     &                 tm, lt, t3, HDLR, HDLV, cierate, ttt, turnoff
c
c  Slice locals
c 
      integer indixe(in)
      real t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &     tgas(in), tgasold(in), tdef_inv(in)
      real*8 edotplus(ijk), edotminus(ijk)
      double precision ttgas
      real temp_table(5)
c
c  Cooling/heating slice locals
c
      double precision ceHI(in), ceHeI(in), ceHeII(in),
     &     ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &     reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &     brem(in)
      real hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &     gpldl(in), gphdl(in), hdlte(in), hdlow(in), 
     &     hdcool(in, 5), hdcoolval(in), cieco(in)
      real gaHI(in), gaH2(in), gaHe(in), gaHp(in), gael(in),
     &     galdl(in)
c
c     Iteration mask
      
      logical itmask(in)
      integer inneriter
      real oldgamma2(in)

      real g2p(10)
      data g2p/0.0493452, 2.28021, 0.115357,
     &         0.114188,  2.90778, 0.689708,
     &         64.2416, -9.36334, -0.377266, 69.8091/
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Set log values of start and end of lookup tables
c
      mciecont = 0
      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem= (log(temend) - log(temstart))/real(nratec-1)
c
c     Set units
c
      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
c      write(6,*) "coolunit: ", coolunit
c      write(6,*) "dom:      ", dom
      zr       = 1.0/(aye*uaye) - 1.0
      fudge    = 1.0
      ciefudge = 1.0

      do i = is+1, ie+1
        edotminus(i) = 0.
        edotplus(i) = 0.
      enddo

c
c     Set compton cooling coefficients (and temperature)
c
      if (iexpand .eq. 1) then
         comp1 = compa * (1.0 + zr)**4
         comp2 = 2.73  * (1.0 + zr)
      else
         comp1 = tiny
         comp2 = tiny
      endif
c
c     Compute Pressure
c
      do i = is+1, ie+1
      if (itmask(i)) then
        p2d(i) = (gamma - 1.0)*d(i,j,k)*ge(i)
      endif
      enddo
c
c     Compute temperature
c
      do i = is+1, ie+1
         if (itmask(i)) then
            tgas(i) = 
     &           (HeI(i) + HeII(i) + HeIII(i))/4.0 +
     &           HI(i) + HII(i) + de(i)
         endif
      enddo
c
c          (include molecular hydrogen, but ignore deuterium)
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
               tgas(i) = tgas(i) +
     &              HM(i) + (H2I(i) + H2II(i))/2.0
            endif
         enddo
      endif
c
      do i = is+1, ie+1
        if (itmask(i)) then
           tgasold(i) = tgas(i)
           tgas(i) = max(p2d(i)*(utem/tgas(i)), temstart)
c           write(6,*), 'tgas ', tgas(i)
           if(tgas(i).ne.tgas(i))then
              write(6,*) 'cool1d_multi tgas0a:',p2d(i),tgas(i)
              write(6,*) 'cool1d_multi tgas0b:',HI(i),HII(i),HM(i)
              write(6,*) 'cool1d_multi tgas0c:',H2I(i)
           endif
        endif
      enddo
c   Omukai
c
c     Correct temperature for gamma from H2
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
           nH2 = 0.5*(H2I(i) + H2II(i))
           nother = (HeI(i) + HeII(i) + HeIII(i))/4.0 +
     &          HI(i) + HII(i) + de(i) + HM(i)
           if(ispecies.gt.2)then
             nother = nother + (DI(i)+DII(i))/2.0
     &                       + (HDI(i)/3.0)
           endif
#define OLDGAMMA2
#ifdef OLDGAMMA2
           oldgamma2(i) = gamma
           do inneriter=1,1
             if (nH2/nother .gt. 1.0e-3) then
               x = 6100/tgas(i) ! not quite self-consistent
               if (x .gt. 10.0) then
                 gamma2(i) = 0.5*5.0
               else
                 gamma2(i) = 0.5*(5.0 + 2.0*x**2 * exp(x)
     &                     /(exp(x)-1)**2)
               endif
             else
               gamma2(i) = 2.5
             endif
             gamma2(i) = 1.0 + (nH2 + nother)/
     &            (nH2*gamma2(i) + nother/(gamma-1.0))
             tgas(i) = tgas(i) * (gamma2(i) - 1.0)/(oldgamma2(i) - 1.0)
             oldgamma2(i) = gamma2(i)
           enddo
#else
           oldgamma2(i) = gamma
           do inneriter=1,1
             x = log10(tgas(i))
             if(nH2/nother.lt.-1.0e-3.or.tgas(i).lt.-50)then
               gamma2(i) = 2.5 !gamma
             else
               gamma2(i) = g2p(1)*exp(-(x-g2p(2))**2.0/g2p(3))
     &                   + g2p(4)*exp(-(x-g2p(5))**2.0/g2p(6))
     &                   + exp(-g2p(7)*(x**g2p(8)))
     &                     *(g2p(9)+x**(-1.0d0*g2p(10)))
     &                   + 5.0/3.0
             endif
             gamma2(i) = 1.0 + (nH2 + nother)/
     &            (nH2/(gamma2(i)-1.0) + nother/(gamma-1.0))
             x = tgas(i)
             tgas(i) = tgas(i) * (gamma2(i) - 1.0)/(oldgamma2(i) - 1.0)
             oldgamma2(i) = gamma2(i)
             if(tgas(i).ne.tgas(i).or.tgas(i).lt.0.0d0)then
               write(0,*)'tgas issue', tgas(i), gamma2(i), x, nH2,
     &                      nother, gamma
               stop
             endif
           enddo
#endif
c         write(0,*) tgas(i)
         endif
         enddo
      endif
c      write(0,*) tgas(8)
c
c     If this is the first time through, just set tgasold to tgas
c
      if (iter .eq. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then 
c               write(6,*) i, tgas(i)
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
c         logtem(i) = log(0.5*(tgas(i)+tgasold(i)))
c Changed to be non-time-centered
         logtem(i) = log(tgas(i))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)
c
c        Compute index into the table and precompute parts of linear interp
c
         indixe(i) = min(nratec-1,
     .                  max(1,int((logtem(i)-logtem0)/dlogtem)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = t2(i) - t1(i)
         tdef_inv(i) = 1.0 / tdef(i)
c
c        Lookup cooling values and do a linear temperature in log(T)
c
         turnoff = 1.0
         if(d(i,j,k)*dom.gt.1e18)turnoff=0.0
         ceHI(i) = turnoff*(ceHIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ceHIa(indixe(i)+1) -ceHIa(indixe(i)))*tdef_inv(i))
         ceHeI(i) = turnoff*(ceHeIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ceHeIa(indixe(i)+1) -ceHeIa(indixe(i)))*tdef_inv(i))
         ceHeII(i) = turnoff*(ceHeIIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ceHeIIa(indixe(i)+1) -ceHeIIa(indixe(i)))*tdef_inv(i))
         ciHI(i) = (ciHIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHIa(indixe(i)+1) -ciHIa(indixe(i)))*tdef_inv(i))!*turnoff
         ciHeI(i) = (ciHeIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHeIa(indixe(i)+1) -ciHeIa(indixe(i)))*tdef_inv(i))!*turnoff
         ciHeIS(i) = (ciHeISa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHeISa(indixe(i)+1) -ciHeISa(indixe(i)))*tdef_inv(i))!*turnoff
         ciHeII(i) = (ciHeIIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHeIIa(indixe(i)+1) -ciHeIIa(indixe(i)))*tdef_inv(i))!*turnoff
         reHII(i) = (reHIIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(reHIIa(indixe(i)+1) -reHIIa(indixe(i)))*tdef_inv(i))!*turnoff
         reHeII1(i)=(reHeII1a(indixe(i)) + (logtem(i) - t1(i))
     .        *(reHeII1a(indixe(i)+1)-reHeII1a(indixe(i)))*tdef_inv(i))!*turnoff
         reHeII2(i)=(reHeII2a(indixe(i)) + (logtem(i) - t1(i))
     .        *(reHeII2a(indixe(i)+1)-reHeII2a(indixe(i)))*tdef_inv(i))!*turnoff
         reHeIII(i)=(reHeIIIa(indixe(i)) + (logtem(i) - t1(i))
     .        *(reHeIIIa(indixe(i)+1)-reHeIIIa(indixe(i)))*tdef_inv(i))!*turnoff
         brem(i) = turnoff*(brema(indixe(i)) + (logtem(i) - t1(i))
     .         *(brema(indixe(i)+1) -brema(indixe(i)))*tdef_inv(i))
         cieco(i) = ciecoa(indixe(i)) + (logtem(i) - t1(i))
     .      *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))*tdef_inv(i)
      endif
      enddo
c
c     Compute the cooling function
c
      do i = is+1, ie+1
         if (itmask(i)) then
         
!         if(omask(i))write(0,*)i,j,k,'a0',edotminus(i),edotplus(i),tgas(i)
         edotminus(i) = edotminus(i) + abs(
c
c                    Collisional excitations
c
     .             - (ceHI  (i)*HI  (i))*de(i)           ! ce of HI
     .             - ((ceHeI (i)*HeII(i))*de(i))
     .                  *de(i)*dom/4.0                       ! ce of HeI
     .             - (ceHeII(i)*HeII(i))*de(i)/4.0         ! ce of HeII
     .                 ) 
!         if(omask(i))write(0,*)i,j,k,'a1',edotminus(i),edotplus(i),tgas(i)
c
c                    Collisional ionizations
c
         edotminus(i) = edotminus(i) + abs(
     .             - (ciHI  (i)*HI  (i))*de(i)             ! ci of HI
     .             - (ciHeI (i)*HeI (i))*de(i)/4.0         ! ci of HeI
     .             - (ciHeII(i)*HeII(i))*de(i)/4.0         ! ci of HeII
     .             - ((ciHeIS(i)*HeII(i))*de(i))
     .                  *de(i)*dom/4.0  ! ci of HeIS
     .                 ) 
!         if(omask(i))write(0,*)i,j,k,'a2',edotminus(i),edotplus(i),tgas(i)
c
c                    Recombinations
c
         edotminus(i) = edotminus(i) + abs(
     .             - (reHII  (i)*HII  (i))*de(i)          ! re of HII
     .             - (reHeII1(i)*HeII (i))*de(i)/4.0      ! re of HeII
     .             - (reHeII2(i)*HeII (i))*de(i)/4.0      ! re of HeII
     .             - (reHeIII(i)*HeIII(i))*de(i)/4.0      ! re of HeIII
     .                 ) 
!         if(omask(i))write(0,*)i,j,k,'a3',edotminus(i),edotplus(i),tgas(i)
!         if(omask(i))write(0,*)i,j,k,'re1',reHII(i),HII(i),de(i)
!         if(omask(i))write(0,*)i,j,k,'re2',reHeII1(i),HeII(i),de(i)/4.0
!         if(omask(i))write(0,*)i,j,k,'re3',reHeII2(i),HeII(i),de(i)/4.0
!         if(omask(i))write(0,*)i,j,k,'re4',reHeIII(i),HeIII(i),de(i)/4.0
c
c                    Bremsstrahlung
c
         edotminus(i) = edotminus(i) + abs(
     .             - (brem(i)*(HII(i)+HeII(i)/4.0+HeIII(i)))
     .                        *de(i)
     .                 ) 
!         if(omask(i))write(0,*)i,j,k,'a4',edotminus(i),edotplus(i),tgas(i)
          
c
c                    Compton cooling or heating
c
           if(omask(i))write(0,*)'stats',i,j,k,d(i,j,k),
     &                 HI(i),HII(i),de(i),HeI(i),HeII(i),HeIII(i)
           comp_edot = - comp1*(tgas(i)-comp2)*de(i)/dom
           if (comp_edot.lt.0.0d0) then
!         if(omask(i))write(0,*)i,j,k,'b',edotminus(i),edotplus(i),tgas(i)
             edotminus(i) = edotminus(i)
     .          + abs(comp_edot)
           else
!         if(omask(i))write(0,*)i,j,k,'c',edotminus(i),edotplus(i),tgas(i)
             edotplus(i) = edotplus(i)
     .          + abs(comp_edot)
           endif
c
c                    X-ray compton heating
c
           comp_xray_edot=-comp_xraya*(tgas(i)-comp_temp)*de(i)/dom
           if (comp_xray_edot.lt.0.0d0) then
!         if(omask(i))write(0,*)i,j,k,'d1',edotminus(i),edotplus(i),tgas(i)
             edotminus(i) = edotminus(i) 
     .             + abs(comp_xray_edot)
           else
!         if(omask(i))write(0,*)i,j,k,'d2',edotminus(i),edotplus(i),tgas(i)
             edotplus(i) = edotplus(i)
     .            +  abs(comp_xray_edot)
           endif
         endif
      enddo
c     
c     --- H2 cooling ---
c
      if (ispecies .gt. 1) then
c

#define USE_GLOVER_ABEL2008
#ifdef USE_GLOVER_ABEL2008
         do i = is+1, ie+1
            if ( itmask(i) ) then
            gaHI(i) = gaHIa(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaHIa(indixe(i)+1) - gaHIa(indixe(i)))*tdef_inv(i)
            gaH2(i) = gaH2a(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaH2a(indixe(i)+1) - gaH2a(indixe(i)))*tdef_inv(i)
            gaHe(i) = gaHea(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaHea(indixe(i)+1) - gaHea(indixe(i)))*tdef_inv(i)
            gaHp(i) = gaHpa(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaHpa(indixe(i)+1) - gaHpa(indixe(i)))*tdef_inv(i)
            gael(i) = gaela(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaela(indixe(i)+1) - gaHIa(indixe(i)))*tdef_inv(i)
            gphdl(i) = gphdla(indixe(i)) + (logtem(i) - t1(i))
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))*tdef_inv(i)
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then
            if (ih2optical.eq.1) then
c
c            nH2 = 0.5*H2I(i)
c            nother = (HeI(i) + HeII(i) + HeIII(i))/4.0 +
c     .                HI(i) + HII(i) + de(i)
c            fH2 = nH2/(nH2 + nother)
c            fudge = sqrt((40.0 * 10**(4.8 * 
c     .               sqrt(max(log10(tgas(i)),2.0)-2.0)) / fH2**2)/
c     .               ((nH2 + nother)*dom) )
              fudge = (0.76*d(i,j,k)*dom/8.0e9)**(-0.45)
              fudge = min(fudge, 1.0)
            endif

            galdl(i) = gaHI(i) * HI(i)  + gaH2(i) * H2I(i)
     &               + gaHe(i) * HeI(i) + gaHp(i) * HII(i)
     &               + gael(i) * de(i)
            gphdl1 = gphdl(i)/dom
            edotminus(i) = edotminus(i) + float(ih2co)*fudge*H2I(i)*
     &           gphdl(i)/(1.0 + gphdl1/galdl(i)) / (2.0*dom)

            end if
         enddo
#else

#define USE_GALLI_PALLA1999

c        Use the Galli and Palla (1999) cooling rates for molecular H.
c
#ifdef USE_GALLI_PALLA1999
c
c
         do i = is+1, ie+1
            if ( itmask(i) ) then
            gpldl(i) = gpldla(indixe(i)) + (logtem(i) - t1(i))
     &         *(gpldla(indixe(i)+1) - gpldla(indixe(i)))*tdef_inv(i)
            gphdl(i) = gphdla(indixe(i)) + (logtem(i) - t1(i))
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))*tdef_inv(i)
            cieco(i) = ciecoa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))*tdef_inv(i)
            endif
         enddo

         do i = is+1, ie+1

            if (itmask(i)) then
            fudge = 1.0d0
            if (ih2optical.eq.1) then
c
c            nH2 = 0.5*H2I(i)
c            nother = (HeI(i) + HeII(i) + HeIII(i))/4.0 +
c     .                HI(i) + HII(i) + de(i)
c            fH2 = nH2/(nH2 + nother)
c            fudge = sqrt((40.0 * 10**(4.8 * 
c     .               sqrt(max(log10(tgas(i)),2.0)-2.0)) / fH2**2)/
c     .               ((nH2 + nother)*dom) )
              fudge = (0.76*d(i,j,k)*dom/8.0e9)**(-0.45)
              fudge = min(fudge, 1.0)
              ciefudge = 1.
              tau = ((d(i,j,k)/2e16)*dom)**2.8  ! 2e16 is in units of cm^-3
              tau = max(tau, 1.e-5)
              ciefudge = min((1.-exp(-tau))/tau,1.0)
              fudge = fudge * ciefudge
c              ciefudge = 1.
c              tau = ((d(i,j,k)/2e16)*dom)**2.8  ! 2e16 is in units of cm^-3
c              tau = max(tau, 1.e-5)
c              ciefudge = min((1.-exp(-tau))/tau,1.0)
c              fudge = fudge! * ciefudge
c              if (fudge.lt.1e-2)
c     .          write(6,*) 'fudge:', fudge
            endif
c
            gphdl1 = gphdl(i)/(HI(i)*dom)
!         if(omask(i))write(0,*)i,j,k,'e',edotminus(i),edotplus(i),tgas(i)
            edotminus(i) = edotminus(i) + abs(
     .           float(ih2co)*fudge*H2I(i)*
     .                 gphdl(i)/(1.0 + gphdl1/gpldl(i)) / (2.0*dom) )
!         if(omask(i))write(0,*)i,j,k,'e2',edotminus(i),edotplus(i),tgas(i),gphdl(i),gphdl1,gpldl(i),H2I(i)
            endif
         enddo
c
#else /* USE_GALLI_PALLA1999 */
cc
         do i = is+1, ie+1
            if (itmask(i)) then
            hyd01k(i) = hyd01ka(indixe(i)) + (logtem(i) - t1(i))
     .         *(hyd01ka(indixe(i)+1)-hyd01ka(indixe(i)))*tdef_inv(i)
            h2k01(i) = h2k01a(indixe(i)) + (logtem(i) - t1(i))
     .         *(h2k01a(indixe(i)+1) - h2k01a(indixe(i)))*tdef_inv(i)
            vibh(i) = vibha(indixe(i)) + (logtem(i) - t1(i))
     .         *(vibha(indixe(i)+1) - vibha(indixe(i)))*tdef_inv(i)
            roth(i) = rotha(indixe(i)) + (logtem(i) - t1(i))
     .         *(rotha(indixe(i)+1) - rotha(indixe(i)))*tdef_inv(i)
            rotl(i) = rotla(indixe(i)) + (logtem(i) - t1(i))
     .         *(rotla(indixe(i)+1) - rotla(indixe(i)))*tdef_inv(i)
            cieco(i) = ciecoa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))*tdef_inv(i)
            endif
         enddo
c
         do i = is+1, ie+1
            if (itmask(i)) then
            qq   = 1.2*(HI(i)*dom)**0.77 + 
     .                (H2I(i)*dom/2.)**0.77
            vibl = (HI(i)*hyd01k(i) + 
     .             H2I(i)/2.*h2k01(i))
     .             *dom*8.18e-13
c
            if (ih2optical .eq. 1) then
c            nH2 = 0.5*H2I(i)
c            nother = (HeI(i) + HeII(i) + HeIII(i))/4.0 +
c     .                HI(i) + HII(i) + de(i)
c            fH2 = nH2/(nH2 + nother)
c            fudge = sqrt((40.0 * 10**(4.8 * 
c     .               sqrt(max(log10(tgas(i)),2.0)-2.0)) / fH2**2)/
c     .               ((nH2 + nother)*dom) )
              fudge = (0.76*d(i,j,k)*dom/8.0e9)**(-0.45)
              fudge = min(fudge, 1.0)
            endif
c
!         if(omask(i))write(0,*)i,j,k,'f',edotminus(i),edotplus(i),tgas(i)
            edotminus(i)=edotminus(i) + float(ih2co)*fudge*H2I(i)*(
     .          vibh(i)/(1.0+vibh(i)/max(   vibl     ,tiny)) +
     .          roth(i)/(1.0+roth(i)/max(qq*rotl(i),tiny))     
     .                                                      )/2.0/dom
            endif
         enddo
c
#endif /* USE_GALLI_PALLA1999 */
#endif
c

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
                tau = ((d(i,j,k)/2e18)*dom)**8.0d0  ! 2e16 is in units of cm^-3
                tau = max(tau, 1.e-5)
                ciefudge = ciefudge*min((1.-exp(-tau))/tau,1.0)
c We apply the ciefudge factor at the end, because it should apply to the
c entire continuum
                edotminus(i) = ciefudge*edotminus(i) 
     &                 + ciefudge*H2I(i)*(d(i,j,k)*cieco(i))
              endif
            endif
            enddo
         endif

      endif                     !  ispecies = 2
c
c     --- Cooling from HD ---
c
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
c            hdlte(i) = hdltea(indixe(i)) + (logtem(i) - t1(i))
c     .         *(hdltea(indixe(i)+1) - hdltea(indixe(i)))*tdef_inv(i)
c            hdlow(i) = hdlowa(indixe(i)) + (logtem(i) - t1(i))
c     .         *(hdlowa(indixe(i)+1) - hdlowa(indixe(i)))*tdef_inv(i)
#ifdef PYFORTFIEJ
            ttgas = tgas(i)
            call calc_hdcool(ttgas, temp_table)
            do m = 1, 5
              hdcool(i,m) = temp_table(m)
            enddo
#else
            do m = 1, 5
              hdcool(i,m) = hdcoola(indixe(i),m) + (logtem(i) - t1(i))
     .        *(hdcoola(indixe(i)+1,m)-hdcoola(indixe(i),m))*tdef_inv(i)
            enddo
#endif
            endif
         enddo
c
         do i = is+1, ie+1
            if (itmask(i)) then
c               hdlte1 = hdlte(i)/(HDI(i)*dom/2.0)
c               hdlow1 = max(hdlow(i), tiny)
c               edotminus(i) = edotminus(i) + HDI(i)*
c     .                     (hdlte1/(1.0 + hdlte1/hdlow1)/(2.0*dom))
c   Do we want to consider HII as well?
c   Note that, as per Flower 2000, we do not include H2 collisions
                hdcoolval(i) = 0.0e0
                do m = 1,5
                  hdcoolval(i) = hdcoolval(i) + (hdcool(i,m)
     .                          *(log10(HI(i))*dom)**(m-1))
                enddo
                ciefudge = 1.
                tau = ((d(i,j,k)/2e16)*dom)**2.8  ! 2e16 is in units of cm^-3
                tau = max(tau, 1.e-5)
                ciefudge = min((1.-exp(-tau))/tau,1.0)
c                omask(i) = .true.
!         if(omask(i))write(0,*)i,'h',edotminus(i),edotplus(i),tgas(i),
!     .                               (HDI(i)*10**hdcoolval(i)/coolunit)*ciefudge
c                omask(i) = .false.
                edotminus(i) = edotminus(i)
     .                   + 0*HDI(i)*(10**hdcoolval(i)/coolunit)*ciefudge
            endif
            omask(i) = .false.
         enddo
      endif
c
c     --- Compute (external) radiative heating terms ---
c
c                       Photoionization heating
c
#ifdef RADIATION
      if (iradshield .eq. 0) then
c
c        regular version
c
         if (iradtype .eq. 8) then
c
c           1) heating assuming high energy photons produces secondary
c              electrons which do the heating (Shull & Steenberg, 1985).
c
            do i = is+1, ie+1
               if (itmask(i)) then
               x = max(HII(i)/(HI(i)+HII(i)), 1.0e-4)
               factor = 0.9971*(1.0-(1.0-x**0.2663)**1.3163)
               edotplus(i) = edotplus(i) + float(ipiht)*factor*(
     .                + piHI  *HI  (i)         ! pi of HI
     .                + piHeI *HeI (i)*0.25     ! pi of HeI
     .                + piHeII*HeII(i)*0.25     ! pi of HeII
     .              )/dom
               endif
            enddo
         else
c
c           2) standard heating
c
            do i = is+1, ie+1
               if (itmask(i)) then
               edotplus(i) = edotplus(i) + float(ipiht)*(
     .                + piHI  *HI  (i)         ! pi of HI
     .                + piHeI *HeI (i)*0.25     ! pi of HeI
     .                + piHeII*HeII(i)*0.25     ! pi of HeII
     .              )/dom
               endif
            enddo
         endif
      else
c
c        version with approximate self-shielding
c
         do i = is+1, ie+1
            if (itmask(i)) then
            edotplus(i) = edotplus(i) + float(ipiht)*(
     .                + piHI  *HI  (i)*
     .                   exp(-avgsighp*HI(i)*dom)
     .                + piHeI *HeI (i)*0.25*
     .                   exp(-avgsighep*0.25*HeI(i)*dom)
     .                + piHeII*HeII(i)*0.25*
     .                   exp(-avgsighe2p*0.25*HeII(i)*dom)
     .           )/dom
            endif
         enddo
      endif
c
#endif /* RADIATION */
      return
      end


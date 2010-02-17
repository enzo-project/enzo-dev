#include "fortran.def"
#define TOLERANCE 1e-20
#define RTOL 5e-1
#define ITERTOL 1e-4
c -----------------------------------------------------------
c  This routine is very simple.  It calculates the derivatives
c  of the species.
c
      subroutine primordial_calc_derivs(
     &                 in,  is,  ie,   s,  sp, sm1, dsq, dsp,  fd,
     &                k01, k02, k03, k04, k05, k06, k07, k08, k09, k10,
     &                ar01,ar02,ar03,ar04,ar05,ar06,ar07,ar08,ar09,ar10,
     &                k11, k12, k13, k14, k15, k16, k17, k18, k19, k21,
     &                ar11,ar12,ar13,ar14,ar15,ar16,ar17,ar18,ar19,ar21,
     &                k22, k23, k24, k25, k26, k27, k28, k29, k30, 
     &                ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30,
     &                k31, k50, k51, k52, k53, k54, k55, k56,
     &                ar31,ar50,ar51,ar52,ar53,ar54,ar55,ar56,
     &                k24shield, k25shield, k26shield,
     &                dmask, chunit, ncspecies, doupdate,
     &                tsgamma, dtit, tscapy, checker, mins, iter,
     &                dtcoef, dtot, rtot,
     &                h2heatp, h2heatm, hiheatp, heiheatp, heiiheatp,
     &                tsc, esterr, tgas,
     &                rejected, itererror, eqmask, dom, dt, ttot )
c -------------------------------------------------------------------
c
      implicit NONE
c
c     arguments
c
      integer in, is, ie, ijk, ncspecies, n
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)
      integer mins(ijk), iter
      real*8  chunit, checker, dtcoef(NSPECIES)
      real*8 s(ijk,NSPECIES), sp(ijk,NSPECIES), sm1(ijk,NSPECIES),
     &       dsq(ijk,NSPECIES), dsp(ijk,NSPECIES), tempde, peqHq, peqHp,
     &       temp1, temp2, temp3, temp4, tsgamma(ijk),
     &       tscapy(ijk,NSPECIES), dtot(ijk), rtot(ijk), h2heatp(ijk),
     &       h2heatm(ijk), hiheatp(ijk), heiheatp(ijk), heiiheatp(ijk),
     &       tsc(ijk), esterr(ijk), qbar, qdot, tbar, prevval,
     &       itererror, fd(ijk)
      logical dmask(ijk), doupdate, rejected(ijk), eqmask(ijk)
      real dtit(ijk), dom, dt, ttot(ijk), tgas(ijk)
c
c     Rate values
c
      real*8 r01, r02, r03, r04, r05,
     &       r06, r07, r08, r09, r10,
     &       r11, r12, r13, r14, r15,
     &       r16, r17, r18, r19, r21,
     &       r22, r23, r50, r51, r52,
     &       r53, r54, r55, r56, r61, r62
      real ar01(in), ar02(in), ar03(in), ar04(in), ar05(in),
     &     ar06(in), ar07(in), ar08(in), ar09(in), ar10(in),
     &     ar11(in), ar12(in), ar13(in), ar14(in), ar15(in),
     &     ar16(in), ar17(in), ar18(in), ar19(in), ar21(in),
     &     ar22(in), ar23(in), ar50(in), ar51(in), ar52(in),
     &     ar53(in), ar54(in), ar55(in), ar56(in),
     &     ar24, ar25, ar26, ar27, ar28, ar29, ar30, ar31
      real k01(in), k02(in), k03(in), k04(in), k05(in),
     &     k06(in), k07(in), k08(in), k09(in), k10(in),
     &     k11(in), k12(in), k13(in), k14(in), k15(in),
     &     k16(in), k17(in), k18(in), k19(in), k21(in),
     &     k22(in), k23(in), k50(in), k51(in), k52(in),
     &     k53(in), k54(in), k55(in), k56(in), 
     &     k24shield(in), k25shield(in), k26shield(in),
     &     k24, k25, k26, k27, k28, k29, k30, k31
c
c     locals
c
      integer i, retloc
      logical stopthepresses, flipdmask

      integer HI,HII,de,HeI,HeII,HeIII,H2I,H2II,HM,ge,DI,DII,HDI
      parameter(HI=1,HII=2,HeI=3,HeII=4,HeIII=5,
     &          H2I=6,ge=7,DI=8,DII=9,HDI=10,
     &          de=11,HM=12,H2II=13 ) ! eq species

        stopthepresses = .false.
        flipdmask=.false.
        if(itererror.gt.1)flipdmask=.true.
        itererror = -1.0

        do i = is+1, ie+1
        if(dmask(i).eqv..true.) then
          rejected(i) = .false.
          temp3 = tiny
          mins(i) = -1

#include "primordial_reactions.inc"
c
c-------------HM-----------------------------------------------------
c
        n = HM
        prevval = s(i,n)
        dsq(i,n) = 
     &       + r07
        dsp(i,n) = 
     &       + r08
     &       + r14
     &       + r15
     &       + r16
     &       + r17
     &       + r19
        if(doupdate.eqv..true.)then
          s(i,n) = s(i,n) * (dsq(i,n)/dsp(i,n))
          s(i,n) = max(s(i,n),tiny*rtot(i))
#include "primordial_reactions.inc"
        endif

c
c-------------H2II---------------------------------------------------
c
        n = H2II
        prevval = s(i,n)
        dsq(i,n) = (
     &       + r09
     &       + r11
     &       + r17 ) * 2.0d0
        dsp(i,n) = (
     &       + r10
     &       + r18
     &       + r19 ) * 2.0d0
        if(doupdate.eqv..true.)then
          s(i,n) = s(i,n) * (dsq(i,n)/dsp(i,n))
          s(i,n) = max(s(i,n),tiny*rtot(i))
#include "primordial_reactions.inc"
        endif

c Notes: r10 and r11 included in HI and HII cause problems with the H2 and
c recombination rates.
c
c-------------HII----------------------------------------------------
c
        n = HII
        prevval = s(i,n)
        dsq(i,n) = 
     &       + r01
     &       + r10 !
     &       + r51 !
        dsp(i,n) = 
     &       + r02
     &       + r09 !
     &       + r11 !
     &       + r16 !
     &       + r17 !
     &       + r50 !
        if(doupdate.eqv..true..and.eqmask(i).eqv..false.)then
          s(i,n) = (dsq(i,n)*tsgamma(i)*dtit(i)+tscapy(i,n))
     &            / (1.0d0+tsgamma(i)*dtit(i)*dsp(i,n)/s(i,n))
          temp2 = abs((2.0d0/(tsc(i)+1.0d0)) *
     &         (tsc(i)*s(i,n) - (1.0d0+tsc(i))*sp(i,n) + sm1(i,n)))
     &      /  (RTOL*dtcoef(n)*sp(i,n))
          itererror = max(itererror, abs(s(i,n)-prevval)/prevval)
          if(temp2.gt.temp3)then
            temp3 = temp2
            mins(i) = n
          endif
#include "primordial_reactions.inc"
        endif

c
c-------------electrons----------------------------------------------
c
        n = de
        prevval = s(i,n)
        dsq(i,n) = 
     &       + r01
     &       + r03
     &       + r05
     &       + r08 !
     &       + r14 !
     &       + r15 !
     &       + r17 !
     &       + r61 !
        dsp(i,n) = 
     &       + r02
     &       + r04
     &       + r06
     &       + r07 !
     &       + r18 !
     &       + r62 !
        if(doupdate.eqv..true.)then
c        s(i,n) = (dsq(i,n)*dtit(i)+s(i,n))
c     &         / (1.0d0+dtit(i)*dsp(i,n)/s(i,n))
        s(i,de) = s(i,HII) + s(i,HeII)/4.0d0 + s(i,HeIII)/2.0d0
     &           + s(i,H2II)/2.0d0 - s(i,HM) + s(i,DII)/2.0d0
        s(i,de) = max(tiny,s(i,de))
#include "primordial_reactions.inc"
        endif

c
c-------------HI-----------------------------------------------------
c
        n = HI
        prevval = s(i,n)
        dsq(i,n) = (
     &       + r02
     &       + r11 !
     &       + r12  * 2.0d0
     &       + r14 !
     &       + r15 !
     &       + r16  * 2.0d0 !
     &       + r18  * 2.0d0 !
     &       + r19 !
     &       + r13  * 2.0d0
     &       + r23  * 2.0d0
     &       + r50 !
     &       + r54 !
     &       )
        dsp(i,n) = (
     &       + r01
     &       + r07 !
     &       + r08 !
     &       + r09 !
     &       + r10 !
     &       + r21  * 2.0d0
     &       + r22  * 2.0d0
     &       + r51 !
     &       + r55 !
     &       )
        if(doupdate.eqv..true..and.eqmask(i).eqv..false.)then
          s(i,n) = (dsq(i,n)*tsgamma(i)*dtit(i)+tscapy(i,n))
     &            / (1.0d0+tsgamma(i)*dtit(i)*dsp(i,n)/s(i,n))
          temp2 = abs((2.0d0/(tsc(i)+1.0d0)) *
     &         (tsc(i)*s(i,n) - (1.0d0+tsc(i))*sp(i,n) + sm1(i,n)))
     &      /  (RTOL*dtcoef(n)*sp(i,n))
          itererror = max(itererror, abs(s(i,n)-prevval)/prevval)
          if(temp2.gt.temp3)then
            temp3 = temp2
            mins(i) = n
          endif
#include "primordial_reactions.inc"
        endif
        hiheatp(i) = 1.0d0 * (r01)

c
c-------------H2I----------------------------------------------------
c
        n = H2I
        prevval = s(i,n)
        dsq(i,n) = (
     &       + r08
     &       + r10
     &       + r19
     &       + r21
     &       + r22
     &       + r53 !
     &       + r55 !
     &       ) * 2.0d0
        dsp(i,n) = (
     &       + r11
     &       + r12
     &       + r13
     &       + r23
     &       + r52 !
     &       + r54 !
     &       ) * 2.0d0
! With quantitiyes in number density:
! H2_eq = (sqrt(k13)*sqrt(4*k22*rho+k13)-k13)/(2*k22)
        if(doupdate.eqv..true..and.eqmask(i).eqv..false.)then
          s(i,n) = (dsq(i,n)*tsgamma(i)*dtit(i)+tscapy(i,n))
     &            / (1.0d0+tsgamma(i)*dtit(i)*dsp(i,n)/s(i,n))
          temp2 = abs((2.0d0/(tsc(i)+1.0d0)) *
     &         (tsc(i)*s(i,n) - (1.0d0+tsc(i))*sp(i,n) + sm1(i,n)))
     &      /  (RTOL*dtcoef(n)*s(i,n))
          itererror = max(itererror, abs(s(i,n)-prevval)/prevval)
          if(temp2.gt.temp3)then
            temp3 = temp2
            mins(i) = n
          endif
        endif
#include "primordial_reactions.inc"
        if(s(i,H2I)/fd(i).gt.1e-7)then
          h2heatp(i) = 2.0d0 * (r21 + r22)
          h2heatm(i) = 2.0d0 * (r13 + r23)
        endif
c        if((tgas(i).gt.8000).and.(rtot(i)*dom.gt.1e8))then
c          h2heatp(i) = max((s(i,H2I)-sp(i,H2I))/dtit(i),0.0d0)
c          h2heatm(i) = max((sp(i,H2I)-s(i,H2I))/dtit(i),0.0d0)
c        endif

c
c-------------HeI----------------------------------------------------
c
        n = HeI
        prevval = s(i,n)
        dsq(i,n) = (
     &       + r04
     &       ) * 4.0d0
        dsp(i,n) = (
     &       + r03
     &       ) * 4.0d0
        if(doupdate.eqv..true.)then
        s(i,n) = (dsq(i,n)*tsgamma(i)*dtit(i)+tscapy(i,n))
     &          / (1.0d0+tsgamma(i)*dtit(i)*dsp(i,n)/s(i,n))
        s(i,n) = max(s(i,n), tiny)
#include "primordial_reactions.inc"
        temp2 = abs((2.0d0/(tsc(i)+1.0d0)) *
     &       (tsc(i)*s(i,n) - (1.0d0+tsc(i))*sp(i,n) + sm1(i,n)))
     &    /  (RTOL*dtcoef(n)*sp(i,n))
c        itererror = max(itererror, abs(s(i,n)-prevval)/prevval)
        if(temp2.gt.temp3)then
          temp3 = temp2
          mins(i) = n
        endif
        endif
        heiheatp(i)  = 4.0d0 * (r03)
c        write(0,*) heiheatp(i)

c
c-------------HeII---------------------------------------------------
c
        n = HeII
        prevval = s(i,n)
        dsq(i,n) = (
     &       + r03
     &       + r06
     &       ) * 4.0d0
        dsp(i,n) = (
     &       + r04
     &       + r05
     &       ) * 4.0d0
        if(doupdate.eqv..true.)then
        s(i,n) = (dsq(i,n)*tsgamma(i)*dtit(i)+tscapy(i,n))
     &          / (1.0d0+tsgamma(i)*dtit(i)*dsp(i,n)/s(i,n))
        s(i,n) = max(tiny,s(i,n))
#include "primordial_reactions.inc"
        temp2 = abs((2.0d0/(tsc(i)+1.0d0)) *
     &       (tsc(i)*s(i,n) - (1.0d0+tsc(i))*sp(i,n) + sm1(i,n)))
     &    /  (RTOL*dtcoef(n)*sp(i,n))
c        itererror = max(itererror, abs(s(i,n)-prevval)/prevval)
        if(temp2.gt.temp3)then
c          temp3 = temp2
c          mins(i) = n
        endif
        endif
        heiiheatp(i) = 4.0d0 * (r05)

c
c-------------HeIII--------------------------------------------------
c
        n = HeIII
        prevval = s(i,n)
        dsq(i,n) = (
     &       + r05
     &       ) * 4.0d0
        dsp(i,n) = (
     &       + r06
     &       ) * 4.0d0
        if(doupdate.eqv..true.)then
        s(i,n) = (dsq(i,n)*tsgamma(i)*dtit(i)+tscapy(i,n))
     &          / (1.0d0+tsgamma(i)*dtit(i)*dsp(i,n)/s(i,n))
        s(i,n) = max(tiny,s(i,n))
#include "primordial_reactions.inc"
        temp2 = abs((2.0d0/(tsc(i)+1.0d0)) *
     &       (tsc(i)*s(i,n) - (1.0d0+tsc(i))*sp(i,n) + sm1(i,n)))
     &    /  (RTOL*dtcoef(n)*sp(i,n))
c        itererror = max(itererror, abs(s(i,n)-prevval)/prevval)
        if(temp2.gt.temp3)then
c          temp3 = temp2
c          mins(i) = n
        endif
        endif

c NB:
c   We assume that the charge exchange with the H/H+ species will completely
c   regulate the abundance of D/D+
c Recall:
c   nD = (nDTot - nHD) * (nH * k51) / (nH*k51 + nH*k50)
        temp1 = max(tiny, dtot(i) - (2.0/3.0) * s(i,HDI))
c        s(i,DI)  = max(tiny*rtot(i), temp1 * 
c     &     (s(i,HI) * k51(i)) / (s(i,HI)*k51(i) + s(i,HII)*k50(i)))
c        s(i,DI) = s(i,DI) * (r62 + r51 + r55) / (r61 + r50 + r54)
c        s(i,DII) = s(i,DII) * (r61 + r50) / (r62 + r51)
c        s(i,DI) = max(temp1, tiny*rtot(i))
c        s(i,DII) = max(temp1 - s(i,DI),tiny*rtot(i))
c#include "primordial_reactions.inc"
        
c
c-------------DII----------------------------------------------------
c
        n = DII
        dsq(i,DII) = (  ! DII
     &       + r61
     &       + r50
     &       + r53
     & ) * 2.0d0
        dsp(i,DII) = (  ! DII
     &       + r62
     &       + r51
     &       + r52
     & ) * 2.0d0
        if(doupdate.eqv..true.)then
        s(i,n) = (dsq(i,n)*dtit(i)+sp(i,n))
     &         / (1.0d0+dtit(i)*dsp(i,n)/s(i,n))
        s(i,n) = max(tiny, s(i,n))
#include "primordial_reactions.inc"
        endif

c
c-------------DI-----------------------------------------------------
c
        n = DI
        dsq(i,DI) = (  ! DI
     &       + r62
     &       + r51
     &       + r55
     & ) * 2.0d0
        dsp(i,DI) = (  ! DI
     &       + r61
     &       + r50
     &       + r54
     & ) * 2.0d0
        if(doupdate.eqv..true.)then
        s(i,n) = (dsq(i,n)*dtit(i)+sp(i,n))
     &         / (1.0d0+dtit(i)*dsp(i,n)/s(i,n))
        s(i,n) = max(tiny, s(i,n))
#include "primordial_reactions.inc"
        endif

c
c-------------HDI----------------------------------------------------
c
        n = HDI
        dsq(i,HDI) = (  ! HDI
     &       + r52
     &       + r54
     & ) * 3.0d0
        dsp(i,HDI) = (  ! HDI
     &       + r53
     &       + r55
     & ) * 3.0d0
        if(doupdate.eqv..true.)then
        s(i,n) = (dsq(i,n)*dtit(i)+sp(i,n))
     &         / (1.0d0+dtit(i)*dsp(i,n)/s(i,n))
        s(i,n) = max(tiny, s(i,n))
#include "primordial_reactions.inc"
        endif

c Charge conservation constraint step only gets applied when we actually
c update, except, we probably only want to do it once both corrections are
c done.
        if(doupdate.eqv..true.)then
          s(i,de) = s(i,HII) + s(i,HeII)/4.0d0 + s(i,HeIII)/2.0d0
     &            + s(i,H2II)/2.0d0 - s(i,HM)
          s(i,de) = max(s(i,de),tiny)
        endif

        if(flipdmask.and.(itererror.lt.ITERTOL))dmask(i)=.false.
        esterr(i) = temp3
        stopthepresses = .false.
        do n=1,NSPECIES
          if((dsp(i,n).ne.dsp(i,n)).or.
     &       (dsq(i,n).ne.dsq(i,n)).or.
     &       (s(i,n).ne.s(i,n)).or.
     &       (s(i,n).lt.0.0))then
            write(0,*) "CALC_DERIV PROBLEM",i,n,iter,eqmask(i)
            write(0,*)'a',dsp(i,n),dsq(i,n),s(i,n)
            write(0,*)'b',sp(i,n),tsgamma(i),tscapy(i,n)
            write(0,*)'c',mins(i),doupdate,dtit(i)
            write(0,*)'d',sm1(i,n),esterr(i)
            write(0,*)'e',rejected(i)
            stopthepresses = .true.
          else
            s(i,n) = max(s(i,n),tiny*rtot(i))
          endif
        enddo
        if(stopthepresses.eqv..true.) then
          do n=1,NSPECIES
            write(0,*) n, s(i,n), dsp(i,n), dsq(i,n), sp(i,n)
          enddo
          write(0,*)'done.'
          stop
          dmask(i) = .false.
        endif

#ifdef NOTIMPLEMENTED
        ar01(i) = r01
        ar02(i) = r02
        ar03(i) = r03
        ar04(i) = r04
        ar05(i) = r05
        ar06(i) = r06
        ar07(i) = r07
        ar08(i) = r08
        ar09(i) = r09
        ar10(i) = r10
        ar11(i) = r11
        ar12(i) = r12
        ar13(i) = r13
        ar14(i) = r14
        ar15(i) = r15
        ar16(i) = r16
        ar17(i) = r17
        ar18(i) = r18
        ar19(i) = r19
        ar20(i) = r20
        ar21(i) = r21
        ar22(i) = r22
        ar23(i) = r23
        ar50(i) = r50
        ar51(i) = r51
        ar52(i) = r52
        ar53(i) = r53
        ar54(i) = r54
        ar55(i) = r55
        ar56(i) = r56
#endif

        endif
        enddo

        return
        end

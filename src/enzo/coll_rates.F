#include "fortran.def"
c=======================================================================
c/////////////////////  SUBROUTINE COLL_RATES  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      SUBROUTINE coll_rates(T, k1, k2, k3, k4, k5, k6, k7, k8, k9,
     $     k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k23,
     $     kunit, casebrates)
c
c  COMPUTE MULTISPECIES COLLISIONAL RATES
c
c  written by: Tom Abel
c  date:       
c  modified1: Feb, 2000 by Greg Bryan; adapted to AMR
c  modified2: July, 2010 by Dan Reynolds; added case-B recombination rates
c
c  PURPOSE:
c    Computes various collisional rates (from Tom Abels web page)
c
c  UNITS:
c    cgs / kunit (kunit is a normalizing factor)
c
c  PARAMETERS:
c
c  INPUTS:
C     T is the gas temperature in Kelvin
c     kunit is a normalizing factor that (i.e. outputted rates are
c           cgs / kunit).
c
c  OUTPUTS:
c     k1-k19: rates as given below
c
C     the coefficient numbering is as in Abel etal 1997, NewA, 2.
C ---1:--       HI    + e   -> HII   + 2e
C ---2:--       HII   + e   -> H     + p
C ---3:--       HeI   + e   -> HeII  + 2e
C ---4:--       HeII  + e   -> HeI   + p
C ---5:--       HeII  + e   -> HeIII + 2e
C ---6:--       HeIII + e   -> HeII  + p
C ---7:--       HI    + e   -> HM    + p
C ---8:--       HM    + HI  -> H2I*  + e
C ---9:--       HI    + HII -> H2II  + p
C ---10--       H2II  + HI  -> H2I*  + HII
C ---11--       H2I   + HII -> H2II  + H
C ---12--       H2I   + e   -> 2HI   + e
C ---13--       H2I   + H   -> 3H
C ---14--       HM    + e   -> HI    + 2e
C ---15--       HM    + HI  -> 2H    + e
C ---16--       HM    + HII -> 2HI
C ---17--       HM    + HII -> H2II  + e
C ---18--       H2II  + e   -> 2HI
C ---19--       H2II  + HM  -> HI    + H2I
c
c-----------------------------------------------------------------------
c
#define NO_USE_SAVIN2004
c
      implicit NONE
c
c  Arguments
c
      REAL  k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13,
     $     k14, k15, k16, k17, k18, k19, k23
      double precision kunit, T
      INTEGER casebrates
c
c  Parameters
c
c
c  Locals
c
      INTEGER i
      double precision log_T, log_T_eV, T_ev, xx, dum
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
      
C     ------- Compute various values of T.
      log_T = LOG(T)
      T_eV = T/11605.d0
      log_T_eV = log(T_eV)
      
      k1 = exp(-32.71396786375d0 
     &     + 13.53655609057d0*log_T_eV
     &     - 5.739328757388d0*log_T_eV**2 
     &     + 1.563154982022d0*log_T_eV**3
     &     - 0.2877056004391d0*log_T_eV**4
     &     + 0.03482559773736999d0*log_T_eV**5
     &     - 0.00263197617559d0*log_T_eV**6
     &     + 0.0001119543953861d0*log_T_eV**7
     &     - 2.039149852002d-6*log_T_eV**8) / kunit
      
      IF (T_eV .GT. 0.8) THEN
         k3 = exp(-44.09864886561001d0
     &        + 23.91596563469d0*log_T_eV
     &        - 10.75323019821d0*log_T_eV**2
     &        + 3.058038757198d0*log_T_eV**3
     &        - 0.5685118909884001d0*log_T_eV**4
     &        + 0.06795391233790001d0*log_T_eV**5
     &        - 0.005009056101857001d0*log_T_eV**6
     &        + 0.0002067236157507d0*log_T_eV**7
     &        - 3.649161410833d-6*log_T_eV**8) / kunit
         
         k4 = (1.54d-9*(1.d0+0.3d0/exp(8.099328789667d0/T_eV))
     &        / (exp(40.49664394833662d0/T_eV)*T_eV**1.5d0)
     &        + 3.92d-13/T_eV**0.6353d0) / kunit
         
         k5 = exp(-68.71040990212001d0
     &        + 43.93347632635d0*log_T_eV
     &        - 18.48066993568d0*log_T_eV**2
     &        + 4.701626486759002d0*log_T_eV**3
     &        - 0.7692466334492d0*log_T_eV**4
     &        + 0.08113042097303d0*log_T_eV**5
     &        - 0.005324020628287001d0*log_T_eV**6
     &        + 0.0001975705312221d0*log_T_eV**7
     &        - 3.165581065665d-6*log_T_eV**8) / kunit
      ELSE
         k1 = max(tiny, k1)
         k3 = tiny
         k4 = 3.92d-13/T_eV**0.6353d0 / kunit
         k5 = tiny
      ENDIF

c     redefine k4 if case B recombination rates are requested
      if (casebrates.eq.1) then
         k4 = 1.26d-14 * (5.7067d5/T)**(0.75d0)
      endif

c     set HII recombination rate as either case A or case B
      IF (casebrates.eq.1) THEN
         IF (T .lt. 1.0d9) then
            k2 = 4.881357d-6*T**(-1.5d0) 
     &           * (1.d0+1.14813d2*T**(-0.407d0))**(-2.242d0) / kunit
         ELSE
            k2 = tiny
         ENDIF
      ELSE
         IF ( T .GT. 5500.d0 ) THEN
            k2 = exp(-28.61303380689232d0
     &           - 0.7241125657826851d0*log_T_eV
     &           - 0.02026044731984691d0*log_T_eV**2
     &           - 0.002380861877349834d0*log_T_eV**3
     &           - 0.0003212605213188796d0*log_T_eV**4
     &           - 0.00001421502914054107d0*log_T_eV**5
     &           + 4.989108920299513d-6*log_T_eV**6
     &           + 5.755614137575758d-7*log_T_eV**7
     &           - 1.856767039775261d-8*log_T_eV**8
     &           - 3.071135243196595d-9*log_T_eV**9) / kunit
         ELSE
            k2 = k4
         ENDIF
      ENDIF
      
c     set HeIII recombination rate as either case A or case B
      IF (casebrates.eq.1) THEN
         IF (T .lt. 1.0d9) then
            k6 = 7.8155d-5*T**(-1.5d0) 
     &           * (1.d0+2.0189d2*T**(-0.407d0))**(-2.242d0) / kunit
         ELSE
            k6 = tiny
         ENDIF
      ELSE
         k6 = 3.36d-10/sqrt(T)/(T/1.d3)**0.2d0/(1+(T/1.d6)
     $        **0.7d0) / kunit
      ENDIF
      
      k7 = 6.77d-15*T_eV**0.8779d0 / kunit
      
      IF (T_eV .GT. 0.1d0) THEN
         k8 = exp(-20.06913897587003d0
     &        + 0.2289800603272916d0*log_T_eV
     &        + 0.03599837721023835d0*log_T_eV**2
     &        - 0.004555120027032095d0*log_T_eV**3
     &        - 0.0003105115447124016d0*log_T_eV**4
     &        + 0.0001073294010367247d0*log_T_eV**5
     &        - 8.36671960467864d-6*log_T_eV**6
     &        + 2.238306228891639d-7*log_T_eV**7) / kunit
      ELSE
         k8 = 1.43d-9 / kunit
      ENDIF
      
      k9 = 1.85d-23*T**1.8d0 / kunit
      IF (T .GT. 6.7d3) 
     &     k9 = 5.81d-16*(T/56200.d0)**(-0.6657d0*log10(T/56200.d0)) / 
     &          kunit
      
      k10 = 6.0d-10 / kunit
      
      IF (T_eV .GT. 0.3d0) THEN
         k13 = 1.0670825d-10*T_eV**2.012d0/
     &        (exp(4.463d0/T_eV)*(1.d0+0.2472d0*T_eV)**3.512d0) / kunit
         
#ifdef USE_SAVIN2004
         k11 = (exp(-21237.15d0/T) * 
     $       (- 3.3232183d-07
     $        + 3.3735382d-07 * log_T
     $        - 1.4491368d-07 * log_T**2
     $        + 3.4172805d-08 * log_T**3
     $        - 4.7813720d-09 * log_T**4
     $        + 3.9731542d-10 * log_T**5
     $        - 1.8171411d-11 * log_T**6
     $        + 3.5311932d-13 * log_T**7)) / kunit
#else /* Abel et al. (1997) */
         k11 = exp(-24.24914687731536d0
     &        + 3.400824447095291d0*log_T_eV
     &        - 3.898003964650152d0*log_T_eV**2
     &        + 2.045587822403071d0*log_T_eV**3
     &        - 0.5416182856220388d0*log_T_eV**4
     &        + 0.0841077503763412d0*log_T_eV**5
     &        - 0.007879026154483455d0*log_T_eV**6
     &        + 0.0004138398421504563d0*log_T_eV**7
     &        - 9.36345888928611d-6*log_T_eV**8) / kunit
#endif /* USE_SAVIN2004 */
         
C     k12 = 4.38e-10*exp(-102000.0/T)*T**0.35 / kunit
         k12 = 5.6d-11*exp(-102124.d0/T)*T**0.5d0 / kunit
      ELSE
         k13 = tiny 
         k11 = tiny
         k12 = tiny
      ENDIF
      
      IF (T_eV .GT. 0.04d0) THEN
         k14 = exp(-18.01849334273d0
     &        + 2.360852208681d0*log_T_eV
     &        - 0.2827443061704d0*log_T_eV**2
     &        + 0.01623316639567d0*log_T_eV**3
     &        - 0.03365012031362999d0*log_T_eV**4
     &        + 0.01178329782711d0*log_T_eV**5
     &        - 0.001656194699504d0*log_T_eV**6
     &        + 0.0001068275202678d0*log_T_eV**7
     &        - 2.631285809207d-6*log_T_eV**8) / kunit
      ELSE
         k14 =  tiny
      ENDIF
      
      IF (T_eV .GT. 0.1d0) THEN
         k15 = exp(-20.37260896533324d0
     &        + 1.139449335841631d0*log_T_eV
     &        - 0.1421013521554148d0*log_T_eV**2
     &        + 0.00846445538663d0*log_T_eV**3
     &        - 0.0014327641212992d0*log_T_eV**4
     &        + 0.0002012250284791d0*log_T_eV**5
     &        + 0.0000866396324309d0*log_T_eV**6
     &        - 0.00002585009680264d0*log_T_eV**7
     &        + 2.4555011970392d-6*log_T_eV**8
     &        - 8.06838246118d-8*log_T_eV**9) / kunit
      ELSE
         k15 = 2.56d-9*T_eV**1.78186d0 / kunit
      ENDIF
      
      k16 = 6.5d-9/sqrt(T_eV) / kunit
      
      k17 = 1.d-8*T**(-0.4d0) / kunit
      IF (T .GT. 1.0d4)
     &     k17=4.0d-4*T**(-1.4d0)*exp(-15100.d0/T) / kunit
      
C     k18 = 5.56396e-8/T_eV**0.6035 / kunit
      k18 = 1.d-8 / kunit
      If (T .GT. 617.d0)
     $     k18 = 1.32d-6 * T**(-0.76d0) / kunit
      k19 = 5.d-7*sqrt(100.d0/T) / kunit
      k23 = ((8.125d-8 / sqrt(T))
     $     * exp(-52000.d0/T)
     $     * (1.d0 - exp(-6000.d0/T))) / kunit
      k23 = max(k23, tiny)
      
      RETURN
      END

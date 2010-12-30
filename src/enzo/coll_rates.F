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
      T_eV = T/11605
      log_T_eV = log(T_eV)
      
      k1 = exp(-32.71396786375 
     &     + 13.53655609057*log_T_eV
     &     - 5.739328757388*log_T_eV**2 
     &     + 1.563154982022*log_T_eV**3
     &     - 0.2877056004391*log_T_eV**4
     &     + 0.03482559773736999*log_T_eV**5
     &     - 0.00263197617559*log_T_eV**6
     &     + 0.0001119543953861*log_T_eV**7
     &     - 2.039149852002e-6*log_T_eV**8) / kunit
      
      IF (T_eV .GT. 0.8) THEN
         k3 = exp(-44.09864886561001
     &        + 23.91596563469*log_T_eV
     &        - 10.75323019821*log_T_eV**2
     &        + 3.058038757198*log_T_eV**3
     &        - 0.5685118909884001*log_T_eV**4
     &        + 0.06795391233790001*log_T_eV**5
     &        - 0.005009056101857001*log_T_eV**6
     &        + 0.0002067236157507*log_T_eV**7
     &        - 3.649161410833e-6*log_T_eV**8) / kunit
         
         k4 = (1.54e-9*(1.+0.3/exp(8.099328789667/T_eV))
     &        / (exp(40.49664394833662/T_eV)*T_eV**1.5)
     &        + 3.92e-13/T_eV**0.6353) / kunit
         
         k5 = exp(-68.71040990212001
     &        + 43.93347632635*log_T_eV
     &        - 18.48066993568*log_T_eV**2
     &        + 4.701626486759002*log_T_eV**3
     &        - 0.7692466334492*log_T_eV**4
     &        + 0.08113042097303*log_T_eV**5
     &        - 0.005324020628287001*log_T_eV**6
     &        + 0.0001975705312221*log_T_eV**7
     &        - 3.165581065665e-6*log_T_eV**8) / kunit
      ELSE
         k1 = max(tiny, k1)
         k3 = tiny
         k4 = 3.92e-13/T_eV**0.6353 / kunit
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
         IF ( T .GT. 5500.0 ) THEN
            k2 = exp(-28.61303380689232
     &           - 0.7241125657826851*log_T_eV
     &           - 0.02026044731984691*log_T_eV**2
     &           - 0.002380861877349834*log_T_eV**3
     &           - 0.0003212605213188796*log_T_eV**4
     &           - 0.00001421502914054107*log_T_eV**5
     &           + 4.989108920299513e-6*log_T_eV**6
     &           + 5.755614137575758e-7*log_T_eV**7
     &           - 1.856767039775261e-8*log_T_eV**8
     &           - 3.071135243196595e-9*log_T_eV**9) / kunit
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
         k6 = 3.36e-10/sqrt(T)/(T/1.e3)**0.2/(1+(T/1.e6)
     $        **0.7) / kunit
      ENDIF
      
      k7 = 6.77e-15*T_eV**0.8779 / kunit
      
      IF (T_eV .GT. 0.1) THEN
         k8 = exp(-20.06913897587003
     &        + 0.2289800603272916*log_T_eV
     &        + 0.03599837721023835*log_T_eV**2
     &        - 0.004555120027032095*log_T_eV**3
     &        - 0.0003105115447124016*log_T_eV**4
     &        + 0.0001073294010367247*log_T_eV**5
     &        - 8.36671960467864e-6*log_T_eV**6
     &        + 2.238306228891639e-7*log_T_eV**7) / kunit
      ELSE
         k8 = 1.43e-9 / kunit
      ENDIF
      
      k9 = 1.85e-23*T**1.8 / kunit
      IF (T .GT. 6.7e3) 
     &     k9 = 5.81e-16*(T/56200)**(-0.6657*log10(T/56200)) / kunit
      
      k10 = 6.0e-10 / kunit
      
      IF (T_eV .GT. 0.3) THEN
         k13 = 1.0670825e-10*T_eV**2.012/
     &        (exp(4.463/T_eV)*(1+0.2472*T_eV)**3.512) / kunit
         
         k11 = exp(-24.24914687731536
     &        + 3.400824447095291*log_T_eV
     &        - 3.898003964650152*log_T_eV**2
     &        + 2.045587822403071*log_T_eV**3
     &        - 0.5416182856220388*log_T_eV**4
     &        + 0.0841077503763412*log_T_eV**5
     &        - 0.007879026154483455*log_T_eV**6
     &        + 0.0004138398421504563*log_T_eV**7
     &        - 9.36345888928611e-6*log_T_eV**8) / kunit
         
C     k12 = 4.38e-10*exp(-102000.0/T)*T**0.35 / kunit
         k12 = 5.6e-11*exp(-102124/T)*T**0.5 / kunit
      ELSE
         k13 = tiny 
         k11 = tiny
         k12 = tiny
      ENDIF
      
      IF (T_eV .GT. 0.04) THEN
         k14 = exp(-18.01849334273
     &        + 2.360852208681*log_T_eV
     &        - 0.2827443061704*log_T_eV**2
     &        + 0.01623316639567*log_T_eV**3
     &        - 0.03365012031362999*log_T_eV**4
     &        + 0.01178329782711*log_T_eV**5
     &        - 0.001656194699504*log_T_eV**6
     &        + 0.0001068275202678*log_T_eV**7
     &        - 2.631285809207e-6*log_T_eV**8) / kunit
      ELSE
         k14 =  tiny
      ENDIF
      
      IF (T_eV .GT. 0.1) THEN
         k15 = exp(-20.37260896533324
     &        + 1.139449335841631*log_T_eV
     &        - 0.1421013521554148*log_T_eV**2
     &        + 0.00846445538663*log_T_eV**3
     &        - 0.0014327641212992*log_T_eV**4
     &        + 0.0002012250284791*log_T_eV**5
     &        + 0.0000866396324309*log_T_eV**6
     &        - 0.00002585009680264*log_T_eV**7
     &        + 2.4555011970392e-6*log_T_eV**8
     &        - 8.06838246118e-8*log_T_eV**9) / kunit
      ELSE
         k15 = 2.56e-9*T_eV**1.78186 / kunit
      ENDIF
      
      k16 = 6.5e-9/sqrt(T_eV) / kunit
      
      k17 = 1.0e-8*T**(-0.4) / kunit
      IF (T .GT. 1.0e4)
     &     k17=4.0e-4*T**(-1.4)*exp(-15100.0/T) / kunit
      
C     k18 = 5.56396e-8/T_eV**0.6035 / kunit
      k18 = 1.e-8 / kunit
      If (T .GT. 617.)
     $     k18 = 1.32e-6 * T**(-0.76) / kunit
      k19 = 5.e-7*sqrt(100./T) / kunit
      k23 = ((8.125e-8 / sqrt(T))
     $     * exp(-52000/T)
     $     * (1.0d0 - exp(-6000/T))) / kunit
      k23 = max(k23, tiny)
      
      RETURN
      END

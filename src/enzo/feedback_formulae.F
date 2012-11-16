c=======================================================================
c//////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine SNII_feedback(mtot, age, Edot, Mdot)
c
c  Fitting formula of Type II SNe feedback
c 
c  written by: Shikui Tang
c  date:  Febuary 9 2009
c
c  dE/dt = 1.e52 * (Tot_mass / 1.e11 Msun) * A(t),
c   where A(t) equals:
c     (age[i]/3.16e6 yr) * exp(1-age[i]/3.16e6 yr) , for 0 < t < 3.16e6 yr
c     pow(age[i]/3.16e6 yr, -1.1) , for 3.16e6 yr <= t < 3.16e7 yr
c     0 , for t >= 3.16e7 yr. The integration of A(t/3.16e6) is
c  \intA(x) = e*(1-exp(-x)*(x+1)) (x=t/3.16E6<=1)
c           = e*(1-2/e) + 10(1-x^(-0.1)) (x>1)
c
c  dM/dt = B * dE/dt ,
c  where B is chosen such that 25% of the original star particle mass
c  gets ejected after 3.16e7 yr. Revised it according to A(t):
c  dM/dt = B' * Tot_mass * A(x) = B' * A(t/3.16E6)
c  Given \intA(x<=10)=2.775, B' = 2.85095E-8 /yr/Msun
c         ??? seems 15% is more reasonable, which gives B'=1.71E-8
c
c  Thus dE/dt = 1D41/B' * dM/dt = 3.5075E48 * dM/dt  (25% ejected at x=10)
c                               = 5.8458E48 * dM/dt  (15% ejected at x=10)
c  dM_ZII/dt = YII * dM/dt ,
c  where YII is the metal yield for SN II.
c  But metal is not handled in this subroutine.
c
c  INPUTS:
c
c    mtot  - total mass of the star particle, in the unit of solar mass
c          - Should it be the initial mass of the gas particle ???
c    age   - the age of the star particle, in the unit of year
c
c  OUTPUTS:
c
c    Edot - energy input rate of the star particle, ergs/yr
c    Mdot - mass input rate of the star particle, Msun/yr
c
c  REMARKS:
c    We may use the integral form to make result independant on time step
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments:
c
       R_PREC mtot, age, Edot, Mdot
c
c
c  Locals:
c
       R_PREC tt, At
c
       if (age .gt. 0._RKIND) then
          tt = age/3.16e6_RKIND
          if (age .lt. 3.16e6_RKIND) then
             At = tt * exp(1-tt)
          else if (age .lt. 3.16e7_RKIND) then
             At = tt**(-1.1_RKIND)
          else
             At = 0._RKIND
          endif
c
          Edot = 1.d41 * mtot * At !ergs/yr
c     Mdot = 2.851E-8 * mtot * At     !Msun/yr,  25% return at tt=10
          Mdot = 1.71e-8_RKIND * mtot * At !Msun/yr, 15% return at tt=10
       else
          Edot = 0._RKIND
          Mdot = 0._RKIND
       endif
       return
       end
c       
c
c=======================================================================
c//////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c
c      subroutine SNIa_feedback(mtot, age, Edot, Mdot, tau)
      subroutine SNIa_feedback(age, Mdot, tau) 
c 
c  Fitting formula of Type Ia SNe feedback
c
c  written by: Shikui Tang
c  date:  Febuary 7 2009
c  modified1: March 7 2009
c      Keep the form of Mdot, Edot=Mdot*eSNIa*c^2
c
c  dE/dt = 3.e41 * (Tot_mass/1.e12 Msun) * pow(age[i]/1.5e10 yr,-1.1)
c         * F(t) * Year_in_sec,
c  where F(t) = pow(age[i]/tau, 1.5) / (1 + pow(age[i]/tau, 1.5)), and
c  age[i] is the age in yrs from t=0, i.e., the birth epoch of a stellar
c  population with a stellar mass "Tot_mass".  'tau' is what's called tdyn
c  or tdp in the Fortran source code.
c  dM/dt  = ??
c   1 solar mass of metals per 10^51 erg ??? Iron
c  dM_ZIa = YIa * dM/dt
c
c  The Arguments have the same meaning as those in SNII_feedback
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments:
c
c       R_PREC mtot, age, Edot, Mdot, tau
       R_PREC age, Mdot, tau
c
c
c  Locals:
c
       R_PREC tt
       if (age .le. 3.16e7_RKIND) then
         !!Edot = 0._RKIND
         Mdot = 0._RKIND
       else
         tt = (age/tau)**1.5_RKIND
         !!Edot = mtot * 9.467e36_RKIND * 
     !!&         (age/1.5e10_RKIND)**(-1.1_RKIND) * tt/(1._RKIND+tt) !!ergs/yr
         !!Edot = Edot *2._RKIND !!In order to match Ken's plot
         !!9.467E36 = 3.e41_RKIND/1.e12_RKIND * Year_in_sec
         !!Mdot = 1.4_RKIND* Edot/1e51_RKIND  !!Msun/yr, 1.4Msun per 1D51 ergs
         !!To match Ken's plot of SNIa mass loss, in unit of Msun/yr
         Mdot = 2.65076e-14_RKIND*(age/1.5e10_RKIND)**(-1.1_RKIND)
     &        * tt/(1._RKIND+tt)
       endif
c
       return
       end
c
c
c=======================================================================
c//////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine PN_feedback(mtot, age, Edot, Mdot)
c
c  Fitting formula of PN  feedback
c
c  written by: Shikui Tang
c  date:  Febuary 7 2009
c  dE/dt = 0
c  dM/dt = (Tot_mass / 1.e11 Msun) * P(t), where P(t) is given by:
c    pow(10, -0.9 * (x-9) + 0.6) , for 7.5 < x < 9.2
c    pow(10, -1.4 * (x-9) + 0.7) , for x > 9.2
c  where x=log(t [yr]), and t=0 is the formation time of a given star
c  particle in our simulations.
c
c  dM_Z/dt = 0 (metal output negligible)
c
c  The Arguments have the same meaning as those in SNII_feedback
c
c  REMARKS:
c    We may use the integral form to make result independant on time step
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments:
c
       R_PREC mtot, age, Edot, Mdot
c
c  Locals:
c
      R_PREC tt
c      
      Edot = 0._RKIND
      tt = log10(age)
      if ( tt .le. 7.5_RKIND) then
         Mdot = 0._RKIND
      else if (tt .le. 9.2_RKIND) then
         Mdot = mtot/1.e11_RKIND * 
     &        10._RKIND**(-0.9_RKIND*(tt-9._RKIND)+0.6_RKIND)
      else
         Mdot = mtot/1.e11_RKIND * 
     &        10._RKIND**(-1.4_RKIND*(tt-9._RKIND)+0.7_RKIND)
      endif
c
      return
      end
c
c=======================================================================
c//////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c      function rmf(age, tau)
      subroutine cal_rmf(age, rmf, mlossSNII, mlossPN)
c
c  The residual mass fraction of a star particle due to mass ejection
c  by SNe II and PN, using the differential dM/dt listed above,
c  normalized to a one Msun stellar population
c
c  INPUTS:
c     age - the age of a stellar population
c  OUTPUT:
c     rmf - return the residual mass of a one Msun stellar population
c     mlossSNII - the cumulative mass ejection due to Type II SNe
c     mlossPN - the cumulative mass ejection due to PN.
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments:
c
       R_PREC rmf, age
c
c  Locals:
c
       R_PREC x, mlossSNII, mlossPN
       R_PREC e
c
c      cumulative mass loss due to SNe II, 15% mass loss when x=10 
c
       if (age .gt. 0._RKIND) then
          x = age/3.16e6_RKIND
          e = exp(1._RKIND)
          if ( x .le. 1._RKIND) then
             mlossSNII = 0.054_RKIND * e * (1._RKIND 
     &            - (x+1._RKIND)/exp(x))
          else if (x .le. 10._RKIND) then
             mlossSNII = 0.054_RKIND * (e*(1._RKIND-2._RKIND/e) + 
     &            10._RKIND*(1._RKIND-x**(-0.1_RKIND)))
          else
             mlossSNII = 0.1498_RKIND
          endif
c     
c     cumulative mass loss due to PN
c
          x = log10(age)
          if ( x .le. 7.5_RKIND) then
             mlossPN = 0._RKIND
          else if ( x .le. 9.2_RKIND) then
             mlossPN = 10._RKIND**(-1.3_RKIND) * 
     &            (age**0.1_RKIND - 10._RKIND**0.75_RKIND)
          else
             mlossPN = 0.239_RKIND - 498.82_RKIND * age**(-0.4_RKIND)
          endif
c
          rmf = 1._RKIND - mlossSNII - mlossPN
       else
          mlossSNII = 0._RKIND
          mlossPN = 0._RKIND
          rmf = 1._RKIND
       endif
       return
       end       

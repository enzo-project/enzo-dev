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
c-----------------------------------------------------------------------
c
c  Arguments:
c
       real mtot, age, Edot, Mdot
c
c
c  Locals:
c
       real tt, At
c
       if (age .gt. 0.0) then
          tt = age/3.16E6       
          if (age .lt. 3.16E6) then
             At = tt * exp(1-tt)
          else if (age .lt. 3.16E7) then
             At = tt**(-1.1)
          else
             At = 0.
          endif
c
          Edot = 1.D41 * mtot * At !ergs/yr
c     Mdot = 2.851E-8 * mtot * At     !Msun/yr,  25% return at tt=10
          Mdot = 1.71D-8 * mtot * At !Msun/yr, 15% return at tt=10
       else
          Edot = 0.0
          Mdot = 0.0
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
c-----------------------------------------------------------------------
c
c  Arguments:
c
c       real mtot, age, Edot, Mdot, tau
       real age, Mdot, tau
c
c
c  Locals:
c
       real tt
       if (age .le. 3.16E7) then
         !!Edot = 0.0
         Mdot = 0.0
       else
         tt = (age/tau)**1.5
         !!Edot = mtot * 9.467E36 * (age/1.5E10)**(-1.1) * tt/(1+tt) !!ergs/yr
         !!Edot = Edot *2 !!In order to match Ken's plot
         !!9.467E36 = 3.e41/1.e12 * Year_in_sec
         !!Mdot = 1.4* Edot/1D51  !!Msun/yr, 1.4Msun per 1D51 ergs
         !!To match Ken's plot of SNIa mass loss, in unit of Msun/yr
         Mdot = 2.65076E-14*(age/1.5E10)**(-1.1)*tt/(1+tt)
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
c-----------------------------------------------------------------------
c
c  Arguments:
c
       real mtot, age, Edot, Mdot
c
c  Locals:
c
      real tt
c      
      Edot = 0.0
      tt = log10(age)
      if ( tt .le. 7.5) then
         Mdot = 0.0
      else if (tt .le. 9.2) then
         Mdot = mtot/1.E11 * 10**(-0.9*(tt-9)+0.6)
      else
         Mdot = mtot/1.E11 * 10**(-1.4*(tt-9)+0.7)
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
c-----------------------------------------------------------------------
c
c  Arguments:
c
       real rmf, age
c
c  Locals:
c
       real x, mlossSNII, mlossPN
       real e
c
c      cumulative mass loss due to SNe II, 15% mass loss when x=10 
c
       if (age .gt. 0.0) then
          x = age/3.16E6
          e = exp(1.0)
          if ( x .le. 1.0) then
             mlossSNII = 0.054 * e * (1 - (x+1)/exp(x))
          else if (x .le. 10.0) then
             mlossSNII = 0.054 * (e*(1-2/e) + 10*(1-x**(-0.1)))
          else
             mlossSNII = 0.1498
          endif
c     
c     cumulative mass loss due to PN
c
          x = alog10(age)
          if ( x .le. 7.5) then
             mlossPN = 0
          else if ( x .le. 9.2) then
             mlossPN = 10**(-1.3) * (age**0.1 - 10**0.75)
          else
             mlossPN = 0.239 - 498.82 * age**(-0.4)
          endif
c
          rmf = 1.0 - mlossSNII - mlossPN
       else
          mlossSNII = 0.0
          mlossPN = 0.0
          rmf = 1.0
       endif
       return
       end       

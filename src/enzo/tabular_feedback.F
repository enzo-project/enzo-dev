#include "fortran.def"
#include "phys_const.def"
#include "error.def"

      module tabular_feedback

      contains
c ===============================================
c 
c     Return mass density in code units yielded by
c     both Type Ia and II supernovae
c 
      subroutine sne_mass(mass_yield,
     &                  mass_init, metal_frac, age, 
     &                  dt, t1,
     &                  mass_table, metal_array, age_array,
     &                  metal_size, age_size)
c 
      implicit none
#include "fortran_types.def"
c 
c     Arguments
c 
      R_PREC mass_yield
      R_PREC mass_init, metal_frac, age
      R_PREC dt, t1
      INTG_PREC metal_size, age_size
      R_PREC mass_table(4, age_size, metal_size)
      R_PREC metal_array(metal_size), age_array(age_size)
c 
c      Locals
c 
      R_PREC mass_Ia, mass_II
      R_PREC dt_yr
      INTG_PREC iSNII, iSNIa
c
c     I'm not clever enough to auto sync these with typedefs.h
      iSNII = 1_IKIND
      iSNIa = 2_IKIND
      dt_yr = dt*t1/yr_s
c 
c      Mass yields are rates (per year) per initial population mass
c      Note that the initial mass is expressed as a mass density
      call integrate_yields(mass_II, metal_frac, age, dt_yr,
     &            mass_table, metal_array, age_array, 
     &            metal_size, age_size, 4_IKIND, iSNII)
       call integrate_yields(mass_Ia, metal_frac, age, dt_yr,
     &            mass_table, metal_array, age_array,
     &            metal_size, age_size, 4_IKIND, iSNIa)
c 
      mass_Ia = mass_Ia * mass_init
      mass_II = mass_II * mass_init
      mass_yield = mass_Ia + mass_II
      end subroutine
c 
c  ===============================================
c  
c      Return metal fraction yielded by Type II supernovae
c 
      subroutine sne_II_metal(fracII_yield,
     &                    mass_init, metal_frac, age, 
     &                    dt, t1,
     &                    metal_table, metal_array, age_array,
     &                    metal_size, age_size)
c 
      implicit none
#include "fortran_types.def"
c 
c      Arguments
c 
      R_PREC fracII_yield
      R_PREC mass_init, metal_frac, age
      R_PREC dt, t1
      INTG_PREC metal_size, age_size
      R_PREC metal_table(4, age_size, metal_size)
      R_PREC metal_array(metal_size), age_array(age_size)
c 
c      Locals
c 
      R_PREC dt_yr
      INTG_PREC iSNII
c 
c     I'm not clever enough to auto sync these with typedefs.h
      iSNII = 1_IKIND
      dt_yr = dt*t1/yr_s
c 
c      SNe numbers are are rates (per year) per initial population mass
c      Note that the initial mass is expressed as a mass density
       call integrate_yields(fracII_yield, metal_frac, age, dt_yr,
     &            metal_table, metal_array, age_array,
     &            metal_size, age_size, 4_IKIND, iSNII)
c 
c            Convert erg per Msun to erg
      fracII_yield = fracII_yield * mass_init
      end subroutine
c 
c  ===============================================
c  
c      Return metal fraction yielded by Type Ia supernovae
c 
      subroutine sne_Ia_metal(fracIa_yield,
     &                    mass_init, metal_frac, age, 
     &                    dt, t1,
     &                    metal_table, metal_array, age_array,
     &                    metal_size, age_size)
c 
      implicit none
#include "fortran_types.def"
c 
c      Arguments
c 
      R_PREC fracIa_yield
      R_PREC mass_init, metal_frac, age
      R_PREC dt, t1
      INTG_PREC metal_size, age_size
      R_PREC metal_table(4, age_size, metal_size)
      R_PREC metal_array(metal_size), age_array(age_size)
c 
c      Locals
c 
      R_PREC dt_yr
      INTG_PREC iSNIa
c 
c     I'm not clever enough to auto sync these with typedefs.h
      iSNIa = 2_IKIND
      dt_yr = dt*t1/yr_s
c 
c      SNe numbers are are rates (per year) per initial population mass
c      Note that the initial mass is expressed as a mass density
       call integrate_yields(fracIa_yield, metal_frac, age, dt_yr,
     &            metal_table, metal_array, age_array,
     &            metal_size, age_size, 4_IKIND, iSNIa)
c 
c            Convert erg per Msun to erg
      fracIa_yield = fracIa_yield * mass_init
      end subroutine
c 
c  ===============================================
c  
c      Return metal fraction yielded by by
c      both Type Ia and II supernovae
c 
      subroutine sne_metal(frac_yield,
     &                    mass_init, metal_frac, age, 
     &                    dt, t1,
     &                    metal_table, metal_array, age_array,
     &                    metal_size, age_size)
c 
      implicit none
#include "fortran_types.def"
c 
c      Arguments
c 
      R_PREC frac_yield
      R_PREC mass_init, metal_frac, age
      R_PREC dt, t1
      INTG_PREC metal_size, age_size
      R_PREC metal_table(4, age_size, metal_size)
      R_PREC metal_array(metal_size), age_array(age_size)
c 
c      Locals
c 
      R_PREC fracIa_yield, fracII_yield
c 
      call sne_Ia_metal(fracIa_yield, mass_init, metal_frac, age,
     &                  dt, t1, metal_table, metal_array, age_array,
     &                  metal_size, age_size)
      call sne_II_metal(fracII_yield, mass_init, metal_frac, age,
     &                  dt, t1, metal_table, metal_array, age_array,
     &                  metal_size, age_size)
      frac_yield = fracIa_yield + fracII_yield
      end subroutine
c 
c  ===============================================
c  
c      Return energy density in code units yielded by
c      both Type Ia and II supernovae
c 
      subroutine sne_energy(energy_yield,
     &                    mass_init, metal_frac, age, 
     &                    erg_per_snIa, erg_per_snII,
     &                    dx, dt, 
     &                    x1, t1, d1, v1,
     &                    energy_table, metal_array, age_array,
     &                    metal_size, age_size)
c 
      implicit none
#include "fortran_types.def"
c 
c      Arguments
c 
      R_PREC energy_yield
      R_PREC mass_init, metal_frac, age
      R_PREC erg_per_snIa, erg_per_snII
      R_PREC dx, dt, x1, t1, d1, v1
      INTG_PREC metal_size, age_size
      R_PREC energy_table(2, age_size, metal_size)
      R_PREC metal_array(metal_size), age_array(age_size)
c 
c      Locals
c 
      R_PREC num_Ia, num_II
      R_PREC m1, dt_yr, minit_msun
      INTG_PREC iSNII, iSNIa
c 
c     I'm not clever enough to auto sync these with typedefs.h
      iSNII = 1_IKIND
      iSNIa = 2_IKIND
      m1 = d1 * x1**3
      dt_yr = dt*t1/yr_s
c 
c      SNe numbers are are rates (per year) per initial population mass
c      Note that the initial mass is expressed as a mass density
      call integrate_yields(num_II, metal_frac, age, dt_yr,
     &            energy_table, metal_array, age_array, 
     &            metal_size, age_size, 2_IKIND, iSNII)
       call integrate_yields(num_Ia, metal_frac, age, dt_yr,
     &            energy_table, metal_array, age_array,
     &            metal_size, age_size, 2_IKIND, iSNIa)
c 
c            Convert erg per Msun to erg
      minit_msun = mass_init * dx**3 * m1 / SolarMass   
      energy_yield = (num_II*erg_per_snII 
     &             + num_Ia*erg_per_snIa) 
     &             * minit_msun
      energy_yield = energy_yield / (m1*v1**2 * (dx)**3)
c      write(6,*) 'in sne_energy: num_II, num_Ia:',num_II*minit_msun,
c     &           num_Ia*minit_msun
      end subroutine
c ===============================================
c 
c     Return mass density in code units yielded by
c     stellar winds
c 
      subroutine agb_mass(mass_yield,
     &                   mass_init, metal_frac, age, 
     &                   dt, t1,
     &                   mass_table, metal_array, age_array,
     &                   metal_size, age_size)
c 
      implicit none
#include "fortran_types.def"
c 
c     Arguments
c 
      R_PREC mass_yield
      R_PREC mass_init, metal_frac, age
      R_PREC dt, t1
      INTG_PREC metal_size, age_size
      R_PREC mass_table(4, age_size, metal_size)
      R_PREC metal_array(metal_size), age_array(age_size)
c 
c      Locals
c 
      R_PREC dt_yr
      INTG_PREC iAGB
c 
      iAGB = 3_IKIND
      dt_yr = dt*t1/yr_s
c 
c      Mass yields are rates (per year) per initial population mass
c      Note that the initial mass is expressed as a mass density
      call integrate_yields(mass_yield, metal_frac, age, dt_yr,
     &            mass_table, metal_array, age_array,
     &            metal_size, age_size, 4_IKIND, iAGB)
c
      mass_yield = mass_yield * mass_init
      end subroutine
 
c  ===============================================
  
c      Return metal fraction yielded by Type II supernovae
 
      subroutine agb_metal(fracAGB_yield,
     &                    mass_init, metal_frac, age, 
     &                    dt, t1,
     &                    metal_table, metal_array, age_array,
     &                    metal_size, age_size)
c 
      implicit none
#include "fortran_types.def"
c 
c      Arguments
c 
      R_PREC fracAGB_yield
      R_PREC mass_init, metal_frac, age
      R_PREC dt, t1
      INTG_PREC metal_size, age_size
      R_PREC metal_table(4, age_size, metal_size)
      R_PREC metal_array(metal_size), age_array(age_size)
c 
c      Locals
c 
      R_PREC dt_yr
      INTG_PREC iAGB
c 
      iAGB = 3_IKIND
      dt_yr = dt*t1/yr_s
c 
c      SNe numbers are are rates (per year) per initial population mass
c      Note that the initial mass is expressed as a mass density
      call integrate_yields(fracAGB_yield, metal_frac, age, dt_yr,
     &            metal_table, metal_array, age_array,
     &            metal_size, age_size, 4_IKIND, iAGB)
c 
      fracAGB_yield = fracAGB_yield * mass_init
      end subroutine

c  ===============================================
  
c      Time-integrate the yield from particular source type
c      given a table of quantites and the stellar population
c      initial metal fraction and age in years.
c      Timestep dt must also be in years.
c      Time integration will account for table entries that
c      fall between the pop age and age+dt
c 
      subroutine integrate_yields(yield, metal, age, dt,
     &            tab_yield, tab_metal, tab_age,
     &            ntab_metal, ntab_age, num_sources, source_id)
c 
      implicit none
#include "fortran_types.def"
c 
c      Arguments
c 
      INTG_PREC ntab_metal, ntab_age, num_sources, source_id
      R_PREC yield, metal, age, dt
      R_PREC tab_yield(num_sources, ntab_age, ntab_metal)
      R_PREC tab_metal(ntab_metal), tab_age(ntab_age)
c 
c      Locals
c 
      INTG_PREC z_ub, a1_ub, a2_ub ! table index right bounds (for interp)
      R_PREC age1, age2            ! time integration bounds
      R_PREC rate1, rate2, t1, t2 
c 
c      Set time integration bounds
c 
      age1 = age
      if (age .lt. 0) age1 = 0
      age2 = age1 + dt
c 
c      Yield table will need to be interpolated in age & metal frac.
c      Find index to the right (above) of the target values
c 
      z_ub = 1_IKIND
      do while (tab_metal(z_ub) .lt. metal 
     &          .and. z_ub .lt. ntab_metal)
         z_ub = z_ub + 1_IKIND
      enddo
      a1_ub = 1_IKIND
      do while (tab_age(a1_ub) .lt. age1 
     &                  .and. a1_ub .lt. ntab_age)
         a1_ub = a1_ub + 1_IKIND
      enddo
      a2_ub = 1_IKIND
      do while (tab_age(a2_ub) .lt. age2
     &            .and. a2_ub .lt. ntab_age)
         a2_ub = a2_ub + 1_IKIND
      enddo
c 
c      If particle falls outside of age table, exit
c 
      if (a1_ub .eq. ntab_age+1) then
         goto 5
      endif
c 
c      Otherwise, find yield
c 
      yield = 0
c 
c      Case 1: both age (age1) & age+dt (age2) are in the same table "cell"
c 
      if (a2_ub - a1_ub .eq. 0) then
         call bilinear_interp(rate1, metal, age1,
     &         tab_yield, tab_metal, tab_age,
     &         ntab_metal, ntab_age, num_sources, source_id,
     &         z_ub-1, z_ub, a1_ub-1, a1_ub)
        yield = rate1 * dt
     
c 
c      Case 2: age1 and age2 are in different table "cells"
c      Therefore, integrate over all intervening cells
c 
      else
         t1 = age1
c         2D interpolation for left bound of first interval
         call bilinear_interp(rate1, metal, t1,
     &         tab_yield, tab_metal, tab_age,
     &         ntab_metal, ntab_age, num_sources, source_id,
     &         z_ub-1_IKIND, z_ub, a1_ub-1_IKIND, a1_ub)
c 
         do while(a1_ub .lt. a2_ub)
            t2 = tab_age(a1_ub)
c            1D interp (metal only) for right bound of current internval
            call bilinear_interp(rate2, metal, t2,
     &            tab_yield, tab_metal, tab_age,
     &            ntab_metal, ntab_age, num_sources, source_id,
     &            z_ub-1, z_ub, a1_ub, a1_ub)
c            Integrate this interval with trapezoid rule
            yield = yield + (rate1 + rate2) * (t2-t1)/2
c            Shift right to the next interval
            rate1 = rate2
            t1 = t2
            a1_ub = a1_ub + 1
         enddo
c 
c         2D interpolation for right bound of last interval
         call bilinear_interp(rate2, metal, age2,
     &         tab_yield, tab_metal, tab_age,
     &         ntab_metal, ntab_age, num_sources, source_id,
     &         z_ub-1, z_ub, a2_ub-1, a2_ub)
         yield = yield + (rate1 + rate2) * (age2-t1)/2
      endif
 5    end subroutine
c 
c  ===============================================
c  
c      Perform bilinear interpolation using a table
c      Do not extrapolate
c 
      subroutine bilinear_interp(interp, x, y,
     &                           tabinterp, tabx, taby, 
     &                           nx, ny, nsrcs, srcid,
     &                           ix1, ix2, iy1, iy2)
c 
      implicit none
#include "fortran_types.def"
c 
c      Arguments
c 
      INTG_PREC nx, ny, nsrcs, srcid, ix1, ix2, iy1, iy2
      R_PREC x, y
      R_PREC tabinterp(nsrcs,ny,nx), tabx(nx), taby(ny)
c 
c      Locals
c 
      R_PREC x1, x2, y1, y2
      R_PREC fx1y1, fx1y2, fx2y1, fx2y2
c 
c      Return
c 
      R_PREC interp
c 
c      Reset indicies to be in bounds
c 
      if (ix1 .lt. 1) ix1 = 1
      if (iy1 .lt. 1) iy1 = 1
      if (ix2 .gt. nx) ix2 = nx
      if (iy2 .gt. ny) iy2 = ny
c 
c      Find x & y associated with lower & upper bounds
c 
      x1 = tabx(ix1)
      x2 = tabx(ix2)
      y1 = taby(iy1)
      y2 = taby(iy2)
c 
      if (ix1 .eq. ix2 .and. iy1 .eq. iy2) then
c         both x & y fall outside the table
         interp = tabinterp(srcid, iy1, ix1)
      elseif (ix1 .eq. ix2) then
c         only x falls outside the table; interpolate y
         fx1y1 = tabinterp(srcid, iy1, ix1)
         fx1y2 = tabinterp(srcid, iy2, ix1)
         interp = fx1y1*(y2-y) + fx1y2*(y-y1)
         interp = interp / (y2-y1)
      elseif (iy1 .eq. iy2) then
c         only y falls outside the table; interpolate x
         fx1y1 = tabinterp(srcid, iy1, ix1)
         fx2y1 = tabinterp(srcid, iy1, ix2)
         interp = fx1y1*(x2-x) + fx2y1*(x-x1)
         interp = interp / (x2-x1)
      else
c         everything is inside the table; interpolate x & y
         fx1y1 = tabinterp(srcid, iy1, ix1)
         fx1y2 = tabinterp(srcid, iy2, ix1)
         fx2y1 = tabinterp(srcid, iy1, ix2)
         fx2y2 = tabinterp(srcid, iy2, ix2)
c 
         interp = fx1y1*(x2-x)*(y2-y) + fx2y1*(x-x1)*(y2-y) 
     &          + fx1y2*(x2-x)*(y-y1) + fx2y2*(x-x1)*(y-y1)
c 
         interp = interp / ( (x2-x1)*(y2-y1) )
      endif
      end subroutine

       end module tabular_feedback
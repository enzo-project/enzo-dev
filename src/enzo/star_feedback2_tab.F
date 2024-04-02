#include "fortran.def"
#include "phys_const.def"
#include "error.def"


c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback2_tab(nx, ny, nz,
     &                      d, dm, te, ge, u, v, w, metal,
     &                      idual, imetal, hmethod, dt, r, dx, t, z,
     &                      d1, x1, v1, t1,
     &                      distrad, diststep, distcells,
     &                      npart, xstart, ystart, zstart, ibuff,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tcp, metalf, minit, type,
     &                      ergSNII, ergSNIa,
     &                      itracksrc, metalSNII, metalSNIa,
     &                      ntabZ, ntabAge,
     &                      tabZ, tabAge, tabMass, tabMetal, tabEvents)

c
c  modification of Cen & Ostriker feedback (method 1; star_maker2.F)
c  to use tabulated yields
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    te,ge - total energy and gas energy fields
c    u,v,w - velocity fields
c    metal - metallicity density field
c    r     - refinement field (0 if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    hmethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - current particle mass
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    minit    - initial particle mass
c    type     - particle type
c    npart    - particle array size specified by calling routine
c    distrad  - feedback distribution radius in cells
c    diststep - distance in walking steps to deposit feedback
c    distcells - total number of cells over which to distribute feedback 
c
c    ergSNII   - energy in ergs per SNII event
c    ergSNIa   - energy in ergs per SNIa event
c    itracksrc - track SNII & SNIa metals separately?
c    metalSNII - field for SNII metals
c    metalSNIa - field for SNIa metals
c    ntabZ     - number of initial metallicities
c    ntabAge   - number of population ages
c    tabZ      - metal indicies for FB yield table
c    tabAge    - age indicies for FB yield table
c    tabMass   - table of SNe mass rates
c    tabMetal  - table of SNe metal rates
c    tabEvents - table of SNe event rates
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, hmethod,
     &      distrad, diststep, distcells
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, justburn
      P_PREC    xstart, ystart, zstart, t
      P_PREC    xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      R_PREC    minit(npart)
      INTG_PREC type(npart)
      R_PREC    ergSNII, ergSNIa
      INTG_PREC itracksrc, ntabZ, ntabAge
      R_PREC    metalSNII(nx,ny,nz), metalSNIa(nx,ny,nz)
      R_PREC    tabZ(ntabZ), tabAge(ntabAge)
      R_PREC    tabMass(ntabZ, ntabAge, 4), tabMetal(ntabZ, ntabAge, 4)
      R_PREC    tabEvents(ntabZ, ntabAge, 2)
c
c  Locals
c
      INTG_PREC i, j, k, n, ic, jc, kc, stepk, stepj, cellstep, 
     &          iter, z_ub, a1_ub, a2_ub
      R_PREC energy, events_ii, events_ia,
     &     meject, meject_ii, meject_ia,
     &     mzeject, mzeject_ii, mzeject_ia,
     &     dratio, distmass, age, minit_msun, meject_g,
     &     m1, dt_yr
      INTG_PREC iSNII, iSNIa
c
      iSNII = 1_IKIND
      iSNIa = 2_IKIND
      m1 = d1 * dx**3
      dt_yr = dt*t1/yr_s
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
c      write(6,*) 'star_feedback2: start'
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2) then
c
c         Compute index of the cell that the star particle
c           resides in.
c 
         i = int((xp(n) - xstart)/dx,IKIND) + 1
         j = int((yp(n) - ystart)/dx,IKIND) + 1
         k = int((zp(n) - zstart)/dx,IKIND) + 1
c
c         check bounds - if star particle is outside of this grid
c         then exit and give a warning.
c
         if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &          .or. k .lt. 1 .or. k .gt. nz) then
            write(6,*) 'warning: star particle out of grid',i,j,k
            goto 100
         endif
c           
            age = (t-tcp(n)) * t1/yr_s
c
c         Integrate yield tables to get mass yield & total SN events
c
            call integrate_yields(meject_ii, metalf(n), age, dt_yr,
     &            tabMass, tabZ, tabAge, 
     &            ntabZ, ntabAge, iSNII)
            call integrate_yields(meject_ia, metalf(n), age, dt_yr,
     &            tabMass, tabZ, tabAge,
     &            ntabZ, ntabAge, iSNIa)
            call integrate_yields(events_ii, metalf(n), age, dt_yr,
     &            tabEvents, tabZ, tabAge, 
     &            ntabZ, ntabAge, iSNII)
            call integrate_yields(events_ia, metalf(n), age, dt_yr,
     &            tabEvents, tabZ, tabAge,
     &            ntabZ, ntabAge, iSNIa)
c
c          Calculate mass & energy added to each cell
c          Yields of mass and energy are rates (per year) per initial pop mass
c          That is, meject_ii & meject_ia are mass frac per year
c          Note that minit is secretly a mass density
c
            meject = (meject_ii + meject_ia) * minit(n)
            distmass = meject / distcells
            minit_msun = minit(n) * dx**3 * m1 / SolarMass
            meject_g = meject * dx**3 * m1
c           erg per Msun per yr to erg
            energy = (events_ii*ergSNII 
     &              + events_ia*ergSNIa) 
     &              * minit_msun
c           to specific energy density (divide by density later)
            energy = energy / v1**2 / dx**6 / distcells
c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
            mp(n) = mp(n) - meject
c
c          if using distributed feedback, check if particle is
c          too close to the boundary
c
            if (distrad .gt. 0) then
               i = max((1 + ibuff + distrad), 
     &              min((nx - ibuff - distrad), i))
               j = max((1 + ibuff + distrad), 
     &              min((ny - ibuff - distrad), j))
               k = max((1 + ibuff + distrad), 
     &              min((nz - ibuff - distrad), k))
      endif
c
c           Deposit into neighboring cells
c
            do kc = k-distrad,k+distrad
               stepk = abs(kc-k)
               do jc = j-distrad,j+distrad
                  stepj = stepk + abs(jc-j)
                  do ic = i-distrad,i+distrad
                     cellstep = stepj + abs(ic-i)
                     if (cellstep .le. diststep) then
                        dratio = 1._RKIND/(d(ic,jc,kc) + distmass)
c
c           Add energy to energy field
c                        
                        te(ic,jc,kc) = ((te(ic,jc,kc)*d(ic,jc,kc)) +
     &                                  energy*dratio) * dratio
                        if (idual .eq. 1)
     &                   ge(ic,jc,kc) =
     &                       ((ge(ic,jc,kc)*d(ic,jc,kc)) +
     &                        energy*dratio) * dratio
c
c           Metal feedback. Similar to total ejected mass, 
c           table has metal mass yield per initial population mass.
c
                         if (imetal .eq. 1) then
                              call integrate_yields(mzeject_ii,
     &                              metalf(n), age, dt_yr,
     &                              tabMetal, tabZ, tabAge,
     &                              ntabZ, ntabAge, iSNII)
                              call integrate_yields(mzeject_ia,
     &                              metalf(n), age, dt_yr,
     &                              tabMetal, tabZ, tabAge,
     &                              ntabZ, ntabAge, iSNIa)
                           mzeject_ii = mzeject_ii * minit(n)
                           mzeject_ia = mzeject_ia * minit(n)
                           mzeject = mzeject_ii + mzeject_ia
c
c           Field is converted to metal fraction by handling routine.
c
                           metal(ic,jc,kc) = 
     &                          (metal(ic,jc,kc)*d(ic,jc,kc) + 
     &                          mzeject / distcells) * dratio
                           if (itracksrc .eq. 1) then
                              metalSNII(ic,jc,kc) =
     &                          (metalSNII(ic,jc,kc)*d(ic,jc,kc) + 
     &                          mzeject_ii / distcells) * dratio
                              metalSNIa(ic,jc,kc) =
     &                          (metalSNIa(ic,jc,kc)*d(ic,jc,kc) + 
     &                          mzeject_ia / distcells) * dratio
                           endif
                        endif
c
c           Mass and momentum feedback
c
                        u(ic,jc,kc) = u(ic,jc,kc)*d(ic,jc,kc) +
     &                       distmass * up(n)
                        v(ic,jc,kc) = v(ic,jc,kc)*d(ic,jc,kc) +
     &                       distmass * vp(n)
                        w(ic,jc,kc) = w(ic,jc,kc)*d(ic,jc,kc) +
     &                       distmass * wp(n)
                        d(ic,jc,kc) = d(ic,jc,kc) + distmass
                        u(ic,jc,kc) = u(ic,jc,kc)/d(ic,jc,kc)
                        v(ic,jc,kc) = v(ic,jc,kc)/d(ic,jc,kc)
                        w(ic,jc,kc) = w(ic,jc,kc)/d(ic,jc,kc)
c
c           If te is really total energy (and it is unless hmethod=2),
c             then just set this value
c
                        if (hmethod .ne. 2 .and. idual .eq. 1) then
                           te(ic,jc,kc) = 0.5_RKIND*(u(ic,jc,kc)**2 + 
     &                          v(ic,jc,kc)**2 + w(ic,jc,kc)**2) +
     &                          ge(ic,jc,kc)
                        endif
                     endif
                  enddo
               enddo
            enddo
c
 10         continue
         endif
c
 100     continue
c
      enddo
c
c      write(6,*) 'star_feedback2: end'
      return
      end
c
c ===============================================
c 
c     Time-integrate the yield from particular source type
c     given a table of quantites and the stellar population
c     initial metal fraction and age in years.
c     Timestep dt must also be in years.
c     Time integration will account for table entries that
c     fall between the pop age and age+dt
c
      subroutine integrate_yields(yield, metal, age, dt,
     &            tab_yield, tab_metal, tab_age,
     &            ntab_metal, ntab_age, source_id)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC ntab_metal, ntab_age, source_id
      R_PREC yield, metal, age, dt
      R_PREC tab_yield(ntab_metal, ntab_age, 4)
      R_PREC tab_metal(ntab_metal), tab_age(ntab_age)
c
c     Locals
c
      R_PREC age1, age2         ! time integration bounds
      R_PREC z_ub, a1_ub, a2_ub ! table index right bounds (for interp)
      R_PREC rate1, rate2       
c
c     Set time integration bounds
c
      age1 = age
      age2 = age + dt
c
c     Yield table will need to be interpolated in age & metal frac.
c     Find index to the right (above) of the target values
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
c     If particle falls outside of age table, exit
c
      if (a1_ub .eq. 1 .or. a1_ub .eq. ntab_age+1) then
         goto 5
      endif
c
c     Otherwise, find yield
c
      yield = 0
c
c     Case 1: both age (age1) & age+dt (age2) are in the same table "cell"
c
      if (a2_ub - a1_ub .eq. 0) then
         call bilinear_interp(rate1, metal, age1,
     &         tab_yield, tab_metal, tab_age,
     &         ntab_metal, ntab_age, source_id,
     &         z_ub-1, z_ub, a1_ub-1, a1_ub)
      yield = rate1 * dt
     
c
c     Case 2: age1 and age2 are in different table "cells"
c     Therefore, integrate over all interveneing cells
c
      else
c        2D interpolation for left bound of first interval
         call bilinear_interp(rate1, metal, age1,
     &         tab_yield, tab_metal, tab_age,
     &         ntab_metal, ntab_age, source_id,
     &         z_ub-1_IKIND, z_ub, a1_ub-1_IKIND, a1_ub)
         do while(a1_ub .lt. a2_ub)
c           1D interp (metal only) for right bound of current internval
            call bilinear_interp(rate2, metal, tab_age(a1_ub),
     &            tab_yield, tab_metal, tab_age,
     &            ntab_metal, ntab_age, source_id,
     &            z_ub-1, z_ub, a1_ub, a1_ub)
c           Integrate this interval with trapezoid rule
            yield = yield + (rate1 + rate2) * dt/2
c           Shift right to the next interval
            rate1 = rate2
            a1_ub = a1_ub + 1
         enddo
c        2D interpolation for right bound of last interval
         call bilinear_interp(rate2, metal, age2,
     &         tab_yield, tab_metal, tab_age,
     &         ntab_metal, ntab_age, source_id,
     &         z_ub-1, z_ub, a2_ub-1, a2_ub)
         yield = yield + (rate1 + rate2) * dt/2
      endif
 5    end
c
c ===============================================
c 
c     Perform bilinear interpolation using a table
c     Do not extrapolate
c
      subroutine bilinear_interp(interp, x, y,
     &                           tabinterp, tabx, taby, 
     &                           nx, ny, srcid,
     &                           ix1, ix2, iy1, iy2)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, srcid, ix1, ix2, iy1, iy2
      R_PREC x, y
      R_PREC tabinterp(nx,ny,4), tabx(nx), taby(ny)
c
c     Locals
c
      R_PREC x1, x2, y1, y2
      R_PREC fx1y1, fx1y2, fx2y1, fx2y2
c
c     Return
c
      R_PREC interp
c
c     Reset indicies to be in bounds
c
      if (ix1 .lt. 1) ix1 = 1
      if (iy1 .lt. 1) iy1 = 1
      if (ix2 .gt. nx) ix2 = nx
      if (iy2 .gt. ny) iy2 = ny
c
c     Find x & y associated with lower & upper bounds
c
      x1 = tabx(ix1)
      x2 = tabx(ix2)
      y1 = taby(iy1)
      y2 = taby(iy2)
c      write(6,*) 'x, x1, x2',x,x1,x2
c      write(6,*) 'y, y1, y2',y,y1,y2
c
      if (ix1 .eq. ix2 .and. iy1 .eq. iy2) then
c        both x & y fall outside the table
      !   write(6,*) 'case: no interpolation'
         interp = tabinterp(ix1, iy1, srcid)
      elseif (ix1 .eq. ix2) then
c        only x falls outside the table; interpolate y
      !    write(6,*) 'case: interpolate y'
         fx1y1 = tabinterp(ix1, iy1, srcid)
         fx1y2 = tabinterp(ix1, iy2, srcid)
         interp = fx1y1*(y2-y) + fx1y2*(y-y1)
         interp = interp / (y2-y1)
      elseif (iy1 .eq. iy2) then
c        only y falls outside the table; interpolate x
      !   write(6,*) 'case: interpolate x'
         fx1y1 = tabinterp(ix1, iy1, srcid)
         fx2y1 = tabinterp(ix2, iy1, srcid)
         interp = fx1y1*(x2-x) + fx2y1*(x-x1)
         interp = interp / (x2-x1)
      else
c        everything is inside the table; interpolate x & y
      !   write(6,*) 'case: interpolate both'
         fx1y1 = tabinterp(ix1, iy1, srcid)
         fx1y2 = tabinterp(ix1, iy2, srcid)
         fx2y1 = tabinterp(ix2, iy1, srcid)
         fx2y2 = tabinterp(ix2, iy2, srcid)
c
         interp = fx1y1*(x2-x)*(y2-y) + fx2y1*(x-x1)*(y2-y) 
     &          + fx1y2*(x2-x)*(y-y1) + fx2y2*(x-x1)*(y-y1)
c
         interp = interp / ( (x2-x1)*(y2-y1) )
      endif
      return
      end
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
     &          iter, z_index, age_index
      R_PREC energy, 
     &     meject, meject_ii, meject_ia, dratio,
     &     distmass, age
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
c      write(6,*) 'star_feedback2: start'
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2) then
c
c         Determine metallicity and age bounds for the stellar population
c
            age = t-tcp(n)
c
            z_index = ntabZ
            do iter = 1_IKIND, ntabZ
               if (metalf(n) < ntabZ(iter)) then
                  z_index = iter-1
               endif
            enddo
c
            age_index = ntabAge
            do iter = 1_IKIND, ntabAge
               if (age < ntabAge(iter)) then
                  age_index = iter-1
               endif
            enddo
c
c         Interpolate yield tables to get mass, metal frac, & event rate
c         If particle metal & age falls outside table, then what?
c         gotta check for that or I'll segfault!
c         Final index 3 & 4 are unused by this routines
c
            meject_ii = bilinear_interp(metalf(n), age,
     &            z_index, z_index+1, age_index, age_index+1,
     &            tabMass(z_index,   age_index,   1),
     &            tabMass(z_index,   age_index+1, 1),
     &            tabMass(z_index+1, age_index,   1),
     &            tabMass(z_index+1, age_index+1, 1))
c
            meject_ia = bilinear_interp(metalf(n), age,
     &            z_index, z_index+1, age_index, age_index+1,
     &            tabMass(z_index,   age_index,   2),
     &            tabMass(z_index,   age_index+1, 2),
     &            tabMass(z_index+1, age_index,   2),
     &            tabMass(z_index+1, age_index+1, 2))
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
c          calculate mass added to each cell
c
            distmass = meject / distcells
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
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
            mp(n) = mp(n) - mform * m_eject
c
c           Add energy to energy field
c
            do kc = k-distrad,k+distrad
               stepk = abs(kc-k)
               do jc = j-distrad,j+distrad
                  stepj = stepk + abs(jc-j)
                  do ic = i-distrad,i+distrad
                     cellstep = stepj + abs(ic-i)
                     if (cellstep .le. diststep) then
                        dratio = 1._RKIND/(d(ic,jc,kc) + distmass)
                        te(ic,jc,kc) = ((te(ic,jc,kc)*d(ic,jc,kc)) +
     &                                  energy) * dratio
                        if (idual .eq. 1)
     &                   ge(ic,jc,kc) =
     &                       ((ge(ic,jc,kc)*d(ic,jc,kc)) +
     &                        energy) * dratio

c
c           Metal feedback (note that in this function gas metal is
c             a fraction (rho_metal/rho_gas) rather than a density.
c             The conversion has been done in the handling routine)
c
                        if (imetal .eq. 1) then
c
c           "Cen method".  This takes into account gas recycling.
c
                           metal(ic,jc,kc) = 
     &                          (metal(ic,jc,kc)*d(ic,jc,kc) + 
     &                          (mform / distcells) * 
     &                          (yield * (1._RKIND-metalf(n)) + 
     &                          m_eject * metalf(n))) * dratio
            !metal is a fraction
c
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
c     Perform bilinear interpolation
c
      subroutine bilinear_interp(x, y, x1, x2, y1, y2,
     &                           fx1y1, fx1y2, fx2y1, fx2y2,
     &                           interp)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      R_PREC x, y, x1, x2, y1, y2, fx1y1, fx1y2, fx2y1, fx2y2
c
c     Return
c
      R_PREC interp
c
      interp = fx1y1*(x2-x)*(y2-y) + fx2y1*(x-x1)*(y2-y) 
     &       + fx1y2*(x2-x)*(y-y1) + fx2y2*(x-x1)*(y-y1)
c
      interp = interp / ( (x2-x1)*(y2-y1) )
      return,
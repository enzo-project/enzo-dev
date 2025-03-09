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
      use tabular_feedback
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
      R_PREC    tabMass(4, ntabAge, ntabZ), tabMetal(4, ntabAge, ntabZ)
      R_PREC    tabEvents(2, ntabAge, ntabZ)
c
c  Locals
c
      INTG_PREC i, j, k, n, ic, jc, kc, stepk, stepj, cellstep, 
     &          iter, z_ub, a1_ub, a2_ub
      R_PREC energy, events_ii, events_ia,
     &     meject, meject_ii, meject_ia,
     &     mzeject, mzeject_ii, mzeject_ia,
     &     dratio, distmass, age
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
c         Integrate yield tables to get mass & energy yields
c
            call sne_mass(meject, minit(n), metalf(n), age,
     &              dt, t1,
     &              tabMass, tabZ, tabAge, ntabZ, ntabAge)
            call sne_energy(energy, minit(n), metalf(n), age,
     &              ergSNIa, ergSNII,
     &              dx, dt, x1, t1, d1, v1,
     &              tabEvents, tabZ, tabAge, ntabZ, ntabAge)
c
c          Distribute mass & energy among feedback zone
            distmass = meject / distcells
            energy = energy / distcells
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
     &                                  energy) * dratio
                        if (idual .eq. 1)
     &                    ge(ic,jc,kc) =
     &                       ((ge(ic,jc,kc)*d(ic,jc,kc)) +
     &                        energy) * dratio
c
c           Metal feedback. Similar to total ejected mass, 
c           table has metal mass yield per initial population mass.
c
                         if (imetal .eq. 1) then
                              call sne_II_metal(mzeject_ii,
     &                         minit(n), metalf(n), age,
     &                         dt, t1,
     &                         tabMetal, tabZ, tabAge, ntabZ, ntabAge)
                              call sne_Ia_metal(mzeject_ia,
     &                         minit(n), metalf(n), age,
     &                         dt, t1,
     &                         tabMetal, tabZ, tabAge, ntabZ, ntabAge)
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
                        !write(6,*) 'updated dens', d(ic,kc,jc)
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
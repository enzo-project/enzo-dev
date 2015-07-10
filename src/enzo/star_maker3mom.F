#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker3(nx, ny, nz,
     &                      d, dm, temp, u, v, w, cooltime,
     &                      dt, r, metal, zfield1, zfield2,
     &                      dx, t, z, procnum, 
     &                      dunits, x1, vunits, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, mintdyn,
     &                      odthresh, masseff, smthresh, level, np, 
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      imetalSNIa, metalSNIa, metalfSNIa)

c
c  CREATES GALAXY PARTICLES
c
c  written by: Brian O'Shea
c  date:       13 November 2002
c	This file was originally a copy of star_maker2.src,
c	which was originally written by Chris Loken.  As of
c	today, this is intended to be the unigrid version of
c	star_maker2, so the jeans mass and stochastic star
c	formation have been completely removed.  See
c	star_maker2.src for changes previous to 13 Nov. 2002.
c
c  modified1:  20 Dec 2002 by BWO
c       Stochastic star formation is added again.
c       The particle masses are averaged over several cells to avoid
c       the "runaway star particle" phenomenon
c  modified2: 27 Sept 2011 by JHW
c    added metalSNIa & metalfSNIa; included feedback from SN Ia/PN
c    (original changes by M. Ryan Joung)
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    temp  - temperature field
c    u,v,w - velocity fields
c    cooltime - cooling time in code units
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    dunits,x1,vunits,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    odthresh - overdensity threshold (some number * avg. density)
c    masseff - gas-to-mass conversion efficiency ( 0<=masseff<=1 )
c    smthresh - star mass threshold (only creates stars with mass >
c        smthresh unless (random number) < starmass/smthresh )
c    mintdyn  - minimum dynamical time, in years
c    level - current level of refinement
c    procnum - processor number (for output)
c    imetalSNIa - SN Ia metallicity flag (0 - none, 1 - yes)
c
c  OUTPUTS:
c
c    np   - number of particles created
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle
c    metalf   - metallicity fraction of particle
c    nmax     - particle array size specified by calling routine
c    metalfSNIa - metallicity fraction of particle (from SN Ia) ! MKRJ
c
c
c-----------------------------------------------------------------------
       implicit none
c-----------------------------------------------------------------------
c
c  Arguments
c
      integer nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      integer procnum, imetalSNIa
      real    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      real    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      real    r(nx,ny,nz), cooltime(nx,ny,nz)
      real    metal(nx,ny,nz), zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      real    dt, dx, z
      real    dunits, x1, vunits, t1
      REALSUB xstart, ystart, zstart, t
      REALSUB xp(nmax), yp(nmax), zp(nmax)
      real    up(nmax), vp(nmax), wp(nmax)
      real    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      real    metalSNIa(nx,ny,nz), metalfSNIa(nmax)
      real    odthresh, masseff, smthresh, mintdyn
c
      real   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      integer  i, j, k, ii
      real   div, tdyn, dtot
      real   pi, G, sndspdC
      real   isosndsp2, starmass, starfraction, bmass, jeanmass
      double precision msolar
      parameter (pi=pi_val, G=GravConst, sndspdC=1.3095d8,
     &           msolar=SolarMass)
c
      ii = np

!     print*,'star_maker3: imetal is:',imetal

c
c  for each zone, : "star" particle is created if answers to all the
c  following questions are affirmative:
c
c    is this the finest level of refinement ?
c    is the density greater than a critical density ?
c    is the flow convergent ?
c    is the cooling time less than a dynamical time ? 
c    is the gas mass greater than the Jeans mass?
c
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
c
c              1) is this finest level of refinement?
c
               if (r(i,j,k) .ne. 0.d0) goto 10
c
c              2) is density greater than threshold?
c
               if (d(i,j,k) .lt. odthresh) goto 10
c
c              3) is divergence negative?
c                 (the first calculation is face centered for ZEUS, 
c                  the second is cell-centered for PPM)
c
               if (imethod .eq. 2) then
                  div = u(i+1,j  ,k  ) - u(i,j,k)
     &                + v(i  ,j+1,k  ) - v(i,j,k)
     &                + w(i  ,j  ,k+1) - w(i,j,k)
               else
                  div = u(i+1,j  ,k  ) - u(i-1,j  ,k  )
     &                + v(i  ,j+1,k  ) - v(i  ,j-1,k  )
     &                + w(i  ,j  ,k+1) - w(i  ,j  ,k-1)
               endif
               if (div .ge. 0.d0) goto 10
c
c              4) t_cool < t_free-fall (if T < 1.1e4 skip this check)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*dunits
               tdyn  = sqrt(3.d0*pi/32.d0/G/dtot)/t1

               if (tdyn .lt. cooltime(i,j,k) .and. 
     &             temp(i,j,k) .gt. 1.1d4) goto 10
c
c              5) is M > M_Jeans? (this definition involves only baryons under
c                 the assumption that the dark matter is stable, which
c                 implies that the dark matter velocity dispersion is >> 
c                 the sound speed.  This will be true for small perturbations
c                 within large halos).
c
               bmass = d(i,j,k)*dble(dunits)*dble(x1*dx)**3 / msolar
               isosndsp2 = sndspdC * temp(i,j,k)
               jeanmass = pi/(6.d0*sqrt(d(i,j,k)*dble(dunits))) *
     &                    dble(pi * isosndsp2 / G)**1.5d0 / msolar

c
c  THIS IS COMMENTED OUT - NO JEANS MASS CRITERION IN THIS ALGORITHM!!!
c  BWO, 13 NOV 02 (fix 3 dec 02)
c               if (bmass .lt. jeanmass) goto 10
c
c              6) Check to see if star is above threshold (given
c                 in units of M_solar)
c
               starfraction = min(masseff*dt/tdyn, 0.9d0)
               tdyn = max(tdyn, mintdyn*3.15d7/t1)

c
c  STOCHASTIC STAR FORMATION HAS BEEN ADDED AGAIN - BWO 20 Dec 2002
c
#define STOCHASTIC_STAR_FORMATION
c
#ifdef STOCHASTIC_STAR_FORMATION
c
c                 Keep global count of "unfullfilled" star formation
c                 and when total is larger than threshold, then create
c                 a star particle with the threshold mass or 1/2 the
c                 gas in the cell, whichever is smaller.
c
               if (starfraction*bmass .lt. smthresh) then
                  sformsum = sformsum + starfraction*bmass
                  if (sformsum .lt. smthresh) goto 10
                  starfraction = min(smthresh/bmass, 0.5d0)
                  sformsum = sformsum - starfraction*bmass
               endif
#else
c
c              is star mass greater than threshold, then make it.
c              if it's less than threshold, go to the next cell.
c
               if (starfraction*bmass .lt. smthresh) goto 10
#endif
c
c              Create a star particle
c
               ii = ii + 1
               mp(ii)  = starfraction * d(i,j,k)
c               tcp(ii) = t
c               tdp(ii) = tdyn
               tcp(ii) = t + tdyn
               tdp(ii) = dt
               xp(ii) = xstart + (float(i)-0.5d0)*dx
               yp(ii) = ystart + (float(j)-0.5d0)*dx
               zp(ii) = zstart + (float(k)-0.5d0)*dx
c
c              Star velocities averaged over multiple cells to
c              avoid "runaway star particle" phenomenon
c              imethod = 2 is zeus, otherwise PPM

               if (imethod .eq. 2) then
                  up(ii) = 0.5*(u(i,j,k)+u(i+1,j,k))
                  vp(ii) = 0.5*(v(i,j,k)+v(i,j+1,k))
                  wp(ii) = 0.5*(w(i,j,k)+w(i,j,k+1))
               else
                  up(ii) = u(i,j,k)
                  vp(ii) = v(i,j,k)
                  wp(ii) = w(i,j,k)
               endif
c
c              Set the particle metal fraction
c
               if (imetal .eq. 1) then
!                 write(*,'("Setting metal fraction")')
                  metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
               else
!                 write(*,'("Zero metal fraction")')
                  metalf(ii) = 0.d0
               endif
c
c              MKRJ 2/20/08 Do the same for particle metal fraction from SN Ia
c
               if (imetalSNIa .eq. 1) then
                  metalfSNIa(ii) = metalSNIa(i,j,k)    ! in here metal is a fraction
               endif
c
c              Remove mass from grid
c
               d(i,j,k) = (1.d0 - starfraction)*d(i,j,k)
c
c               write(7+procnum,1000) level,bmass*starfraction,tcp(ii),
c     &                           tdp(ii)*t1,d(i,j,k)*dunits,z,metalf(ii)
c
 1000          format(i5,1x,6(1pe10.3,1x))
c
c              Do not generate more star particles than available
c
               if (ii .eq. nmax) goto 20

10          continue

            enddo
         enddo
      enddo
 20   continue
c	
      if (ii .ge. nmax) then
         write(6,*) 'star_maker3: reached max new particle count'
         ERROR_MESSAGE
      endif
      np = ii
c
c      if (np .ne. 0) then
c         write(6,*) 'Stars created: number,time,level: ', np, t, level
c      endif
c
      return
      end
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback3(nx, ny, nz,
     &               d, dm, te, ge, u, v, w,
     &               metal, zfield1, zfield2,
     &               idual, imetal, imulti_metals, imethod, 
     &               dt, r, dx, t, z,
     &               dunits, x1, vunits, t1, sn_param, m_eject, yield,
     &               distrad, diststep, distcells,
     &               npart, xstart, ystart, zstart, ibuff,
     &               xp, yp, zp, up, vp, wp,
     &               mp, tdp, tcp, metalf, type, justburn,
     &               kinf)

c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Brian O,Shea
c  date:       13 November 2002
c  	this is a copy of star_maker2.src, and is designated
c	as of this writing to be for unigrid runs with star 
c	formation, and for experimentation.  See star_maker2.src
c	for the changelog prior to 13 November 2002.
c
c  modified1:  Britton Smith
c  date:       August 2010
c              Added distributed feedback where gas, metals, and 
c              energy are injected evenly over multiple cells.
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
c    dunits,x1,vunits,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    imulti_metals - flag to use multi metals zfield 1 and 2
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c    distrad  - feedback distribution radius in cells
c    diststep - distance in walking steps to deposit feedback
c    distcells - total number of cells over which to distribute feedback
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    npart    - particle array size specified by calling routine
c    sn_param - fraction of stellar rest mass that goes to feedback
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c    type     - particle type
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c    justburn     - time-weighted mass of star formation (code units)
c
c
c-----------------------------------------------------------------------
       implicit none
c-----------------------------------------------------------------------
c
c  Arguments
c
      integer nx, ny, nz, ibuff, npart, idual, imetal, 
     &      imulti_metals, imethod,
     &      distrad, diststep, distcells
      real    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      real    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      real    r(nx,ny,nz), ge(nx,ny,nz)
      real    metal(nx,ny,nz)
      real    zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      real    dt, dx, z
      real    dunits, x1, vunits, t1, justburn
      REALSUB xstart, ystart, zstart, t
      REALSUB xp(npart), yp(npart), zp(npart)
      real    up(npart), vp(npart), wp(npart)
      real    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      integer type(npart)
      real 	  kinf
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      integer n, ic, jc, kc,
     &     iface, jface, kface, i, j, k
      real mform, clight, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, 
     &     dist_mass_cells, dist_mom_cells, mass_per_cell, te_per_cell,
     &     thermal_energy_per_cell, kinetic_energy_per_cell,
     &	   energy_per_cell, mom_per_cell, face_shift,
     &     xface, yface, zface, xpos, ypos, zpos,
     &     energy_before, energy_after, mass_before, mass_after,
     &     kin_energy_before, kin_energy_after,
     &     asum, bsum, csum, delta_ke, ke_injected, ke_after
      real dxf, dyf, dzf, dxc, dyc, dzc
      real u1(4,4,4), v1(4,4,4), w1(4,4,4), d1(4,4,4),
     &     te1(4,4,4), ge1(4,4,4), metal1(4,4,4),
     &     ke_before(4,4,4)
      parameter (clight = c_light, msolar_e51 = 1800.d0)
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
c      write(6,*) 'star_feedback3: start',npart
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2) then
c


c        The star particle creation algorithm partnered with this 
c          feedback algorithm creates a star particle instantaneously.
c          However, we do feedback as if the star particles are created 
c          over a long period of time (in code units), so the particle
c          actually loses mass over time in an exponentially decaying
c          way.

c        Determine how much of a given star particle would have been 
c          turned into stars during this timestep.  Then calculate the mass
c          which should have formed during this timestel dt using the integral
c          form of the Cen & Ostriker formula.

            
c            write(6,*) t, dt, tcp(n), tdp(n),mp(n)
            if (t+dt < tcp(n)) goto 10
            xv1 = max((t      - tcp(n))/tdp(n), 0.0)
c            if (imethod .eq. 3) xv1 = max((t - dt - tcp(n))/tdp(n), 0.0)
c            write(6,*) 'xv1:',xv1
            if (xv1 .gt. 10.d0) goto 10     ! t-tcp >> tdp so ignore
            xv2 = (t + dt - tcp(n))/tdp(n)
            if (imethod .eq. 3) xv2 = (t - tcp(n))/tdp(n)
            write(6,*) 'xv1,xv2:',xv1, xv2
            mform = ((1.d0 + xv1)*exp(-xv1) - (1.d0 + xv2)*exp(-xv2)) 
c
c          skip if very little mass is formed.
c
            if (mform .lt. 1e-2) goto 10

c        First calculate the initial mass of the star particle 
c          in question.
c
            minitial = mp(n) / 
     &              (1.d0 - m_eject*(1.d0 - (1.d0 + xv1)*exp(-xv1)))
c
c         Then, calculate the amount of mass that would have formed in
c         this timestep.
c
            mform = minitial * mform
            mform = max(min(mform, mp(n)), 0.d0)
            write(6,*) 'mform, minitial = ',mform,minitial,dt,tcp(n)
c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
            mp(n) = mp(n) - mform * m_eject
c
c           Record amount of star formation in this grid.
c
            justburn = justburn + mform * dt * dx**3
c
c           Assuming 3x3x3 cube, calculate cell distribution
c
            write(6,*) 'distrad:',distrad
            distrad = 3
            dist_mass_cells = distrad**3
            dist_mom_cells = distrad**2 * 6
c
c           Calculate mass per cell ejected
c
            mass_per_cell = mform * m_eject / dist_mass_cells
            write(6,*) 'mass_per_cell:',mass_per_cell,m_eject,
     &           dist_mass_cells
c
c           Calculate how much of the star formation in this
c           timestep would have gone into supernova energy.
c            
            energy = sn_param * mform * (clight/vunits)**2
            energy_per_cell = energy / dist_mass_cells
            write(6,*) 'energy,energy_per_cell:',energy,energy_per_cell

            thermal_energy_per_cell = (1.0 - kinf) * energy_per_cell
c            kinetic_energy_per_cell = kinf * energy_per_cell 

c            write(6,*) 'kinetic_energy_per_cell:',
c     &           kinetic_energy_per_cell

c           Calculate momentum of this much energy spread into cells
c             mom = sqrt(2 * mform * m_eject * energy)
c

c            mom_per_cell = sqrt(2.0 * m_eject * kinetic_energy_per_cell)

c
c         check bounds - if star particle is outside of this grid
c         then exit and give a warning.
c
            if (xp(n) .lt. xstart .or. xp(n) .gt. xstart+dx*nx .or.
     &          yp(n) .lt. ystart .or. yp(n) .gt. ystart+dx*ny .or.
     &          zp(n) .lt. zstart .or. zp(n) .gt. zstart+dx*ny) then
               write(6,*) 'warning: star particle out of grid',
     &              xp(n),yp(n),zp(n), xstart, ystart, zstart
               goto 100
            endif
c
c           If using zeus, then velocities are face-centered so shift
c
            face_shift = 0.0
            if (imethod .eq. 2) face_shift = 0.5
c
c         Compute index of the first cell to add momentum,
c             accounting for possible face-centering
c
c            write(6,*) 'xp,yp,zp',xp(n), yp(n), zp(n)
c            write(6,*) 'xstart,ystart,zstart', xstart, ystart, zstart
            xface = (xp(n) - xstart)/dx - 0.5 - face_shift
            yface = (yp(n) - ystart)/dx - 0.5 - face_shift
            zface = (zp(n) - zstart)/dx - 0.5 - face_shift
c            write(6,*) 'xface,yface,zface',xface, yface, zface
c
            iface = int(xface + 0.5)
            jface = int(yface + 0.5)
            kface = int(zface + 0.5)
c            write(6,*) 'iface,jface,kface', iface, jface, kface
c
            dxf = real(iface) + 0.5 - xface
            dyf = real(jface) + 0.5 - yface
            dzf = real(kface) + 0.5 - zface
c            write(6,*) 'dxf, dyf, dzf', dxf, dyf, dzf
c
c         Compute index of the first cell to add mass, assuming cell-centering
c 
            xpos = (xp(n) - xstart)/dx - 0.5
            ypos = (yp(n) - ystart)/dx - 0.5
            zpos = (zp(n) - zstart)/dx - 0.5
c            write(6,*) 'xpos, ypos, zpos', xpos, ypos, zpos
c
            ic = int(xpos + 0.5)
            jc = int(ypos + 0.5)
            kc = int(zpos + 0.5)
c            write(6,*) 'ic,jc,kc', ic,jc,kc
c
            dxc = real(ic) + 0.5 - xpos
            dyc = real(jc) + 0.5 - ypos
            dzc = real(kc) + 0.5 - zpos
c            write(6,*) 'dxc, dyc, dzc', dxc, dyc, dzc
c
c           Zero local dummy field and kinetic energy field
c
            do k = 1, 4
               do j = 1, 4
                  do i = 1, 4
                     u1(i,j,k) = 0.0
                     v1(i,j,k) = 0.0
                     w1(i,j,k) = 0.0
                     d1(i,j,k) = 0.0
                     ge1(i,j,k) = 0.0
                     te1(i,j,k) = 0.0
                     metal1(i,j,k) = 0.0
                  enddo
               enddo
            enddo

c
c     Compute the kinetic energy in the affected region
c     before momentum is added (except for ZEUS)
c
            if (imethod .ne. 2) then
               do k = -1, +2
                  do j = -1, +2
                     do i = -1, +2

                        ke_before(i+2, j+2, k+2) = 
     &                       0.5*d(ic+i, jc+j, kc+k)*
     &                    (u(ic+i ,jc+j ,kc+k)**2 + 
     &                     v(ic+i ,jc+j ,kc+k)**2 + 
     &                     w(ic+i ,jc+j ,kc+k)**2)
                     
                     enddo
                  enddo
               enddo
            endif
c
c           First convert velocities to momenta and transform
c            into frame comoving with particle
c
            call momentum(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                    nx, ny, nz, ic, jc, kc, 
     &                    iface, jface, kface, imethod, imetal, +1)
c
c           Sum mass and energy before
c
c     change this to sum_mass_kin_energy
            call sum_mass_energy(u, v, w, d, ge, te, nx, ny, nz,
     &                           iface, jface, kface, ic, jc, kc,
     &                           mass_before, energy_before, 
     &                           kin_energy_before, idual, imethod)
c            write(6,*) "m,e before:", mass_before, energy_before
c
c           Now add mass and momentum terms (normalization 1.0) to
c              local dummy fields
c
            call add_feedback(u1, v1, w1, d1, ge1, te1, metal1, 4, 4, 4, 
     &                        2, 2, 2, 2, 2, 2,
     &                        dxf, dyf, dzf, dxc, dyc, dzc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf(n),
     &                        mass_per_cell, 1.0, 0.0)
c
c           The kinetic energy after the feedback event is a quadratic equation
c           with delta_p, where delta_p is the total amount of momentum added
c           to the feedback region:
c           E_k,a = asum + bsum*delta_p + csum*delta_p*delta_p
c
c           asum = sum(p_i,b*p_i,b/(m_i,b + delta_m_i))
c
c
c           Sum a, b and c terms to get momentum normalization
c
            if (kinf .gt. 0.0) then
               call sum_abc(u, v, w, d, ge, u1, v1, w1, d1, 
     &              nx, ny, nz, iface, jface, kface, ic, jc, kc,
     &              asum, bsum, csum)
               write(6,*) 'abc1',asum, bsum, csum
               asum = asum - (kin_energy_before + kinf*energy)
               write(6,*) 'abc2',asum, bsum, csum
c
c           Calculate momentum contribution
c
               mom_per_cell = (-bsum + sqrt(bsum**2 - 4.0*asum*csum))/
     &              (2.0*csum)
            else
c
c           If the kinetic fraction is set to 0, then 0 momentum is added
c           to each cell.
c
               mom_per_cell = 0.0
            endif

            write(6,*) 'mom_per_cell:',mom_per_cell
c
c           Now add mass and momentum, using three-point CIC
c
c            te_per_cell = 0.0
            call add_feedback(u, v, w, d, ge, te, metal, nx, ny, nz, 
     &                        ic, jc, kc, iface, jface, kface, 
     &                        dxf, dyf, dzf, dxc, dyc, dzc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf(n),
     &                        mass_per_cell, mom_per_cell, 
     &			      thermal_energy_per_cell)
c
c           Sum mass and energy before
c
            call sum_mass_energy(u, v, w, d, ge, te, nx, ny, nz,
     &                           iface, jface, kface, ic, jc, kc,
     &                           mass_after, energy_after, 
     &                           kin_energy_after, idual, imethod)
c
            write(6,*) 'energy (b,a,preda):',energy_before,energy_after, 
     &           energy_before + energy,energy
            write(6,*) 'mass (b,a,preda):', mass_before, mass_after, 
     &           mass_before + mform*m_eject,mform*m_eject
            write(6,*) 'kinetic energy (b,a,preda):',kin_energy_before,
     &           kin_energy_after, kin_energy_before + kinf*energy, 
     &           kinf*energy
            write(6,*) 'idual:',idual
c
c           Convert momenta back to velocities and transform back to lab frame
c
            call momentum(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                    nx, ny, nz, ic, jc, kc, 
     &                    iface, jface, kface, imethod, imetal, -1)
c
c           Add the increase in the kinetic energy to the total energy field
c             (unless we're using Zeus).  If using dual energy formalism, we
c             might want to enforce consistency.
c     
            if (imethod .ne. 2) then
               ke_injected = 0.0
               do k = -1, +2
                  do j = -1, +2
                     do i = -1, +2

                        ke_after = 0.5d0*
     &                        d(ic+i, jc+j, kc+k) *
     &                       (u(ic+i ,jc+j ,kc+k)**2 + 
     &                        v(ic+i ,jc+j ,kc+k)**2 + 
     &                        w(ic+i ,jc+j ,kc+k)**2) 

                        delta_ke = ke_after - ke_before(i+2,j+2,k+2)

c                        write(6,*) 'delta_ke:',ke_before(i+2,j+2,k+2),
c     &                       ke_after,delta_ke

                        te(ic+i ,jc+j ,kc+k) = 
     &                       te(iface+i ,jc+j ,kc+k) + 
     &                       delta_ke/d(ic+i, jc+j, kc+k)

                        ke_injected = ke_injected + delta_ke

                     enddo
                  enddo
               enddo
c

               write(6,*) 'ke_injected:',ke_injected
               if (abs((ke_injected - kinf*energy)/
     &                    ke_injected) .gt. 0.01 .and. 
     &              kinf .ne. 0.0) then
                  write(6,*) 'star_maker3: ke_injected does not match'
                  write(6,*) ke_injected,kinf*energy
                  write(6,*) kin_energy_before,kin_energy_after
                  write(6,*) kin_energy_after - kin_energy_before
                  ERROR_MESSAGE
               endif
            endif
c
 10         continue
         endif
c
 100     continue
c
      enddo
c
c      write(6,*) 'star_feedback3: end'
      return
      end
c
c ==========================================================
c
c     Convert velocities to momentum and back
c
      subroutine momentum(u, v, w, d, metal, up, vp, wp, 
     &                    nx, ny, nz, ic, jc, kc, 
     &                    iface, jface, kface, imethod, imetal, idir)
c
      implicit none
c
c     Arguments
c
      integer nx, ny, nz, ic, jc, kc
      integer iface, jface, kface, imethod, imetal, idir
      real    d(nx, ny, nz), metal(nx,ny,nz)
      real    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      real    up, vp, wp
c
c     Locals
c
      integer i, j, k
c
c     Error check
c
      if (idir .ne. -1 .and. idir .ne. 1) then
         write(6,*) 'incorrect idir value in momentum call'
         ERROR_MESSAGE
      endif
c
c     Loop over velocities, multiplying by densities (or dividing if
c       converting back)
c
      do k = -1, +2
         do j = -1, +2
            do i = -1, +2
c
c              idir = +1: convert vel -> mom
c               
               if (idir. eq. +1) then
                  if (imethod .eq. 2) then
                     u(iface+i ,jc+j ,kc+k) = 
     &                            (u(iface+i ,jc+j ,kc+k) - up)
     &                    * 0.5 * (d(iface+i  , jc+j, kc+k) +
     &                             d(iface+i+1, jc+j, kc+k))
                     v(ic+i ,jface+j ,kc+k) = 
     &                            (v(ic+i ,jface+j ,kc+k) - vp)
     &                    * 0.5 * (d(ic+i, jface+j  , kc+k) + 
     &                             d(ic+i, jface+j+1, kc+k))
                     w(ic+i ,jc+j ,kface+k) = 
     &                            (w(ic+i ,jc+j ,kface+k) - wp)
     &                    * 0.5 * (d(ic+i, jc+j, kface+k  ) + 
     &                             d(ic+i, jc+j, kface+k+1))
                  else
                     u(ic+i ,jc+j ,kc+k) = (u(ic+i ,jc+j ,kc+k)-up) *
     &                                      d(ic+i, jc+j, kc+k)
                     v(ic+i ,jc+j ,kc+k) = (v(ic+i ,jc+j ,kc+k)-vp) *
     &                                      d(ic+i, jc+j, kc+k)
                     w(ic+i ,jc+j ,kc+k) = (w(ic+i ,jc+j ,kc+k)-wp) *
     &                                      d(ic+i, jc+j, kc+k)
                  endif
                  if (imetal .eq. 1) then
                     metal(ic+i,jc+j,kc+k) = metal(ic+i,jc+j,kc+k) *
     &                    d(ic+i,jc+j,kc+k)
                  endif
c
c              if idir = -1: convert mom -> vel
c               
               else
                  if (imethod .eq. 2) then
                     u(iface+i ,jc+j ,kc+k) = u(iface+i ,jc+j ,kc+k)
     &                    /( 0.5 * (d(iface+i  , jc+j, kc+k) +
     &                              d(iface+i+1, jc+j, kc+k))) + up
                     v(ic+i ,jface+j ,kc+k) = v(ic+i ,jface+j ,kc+k) 
     &                    /( 0.5 * (d(ic+i, jface+j  , kc+k) + 
     &                              d(ic+i, jface+j+1, kc+k))) + vp
                     w(ic+i ,jc+j ,kface+k) = w(ic+i ,jc+j ,kface+k) 
     &                    /( 0.5 * (d(ic+i, jc+j, kface+k  ) + 
     &                              d(ic+i, jc+j, kface+k+1))) + wp
                  else
                     u(ic+i ,jc+j ,kc+k) = u(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + up
                     v(ic+i ,jc+j ,kc+k) = v(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + vp
                     w(ic+i ,jc+j ,kc+k) = w(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + wp
                  endif
                  if (imetal .eq. 1) then
                     metal(ic+i,jc+j,kc+k) = metal(ic+i,jc+j,kc+k) /
     &                    d(ic+i,jc+j,kc+k)
                  endif
               endif
c
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine sum_mass_energy(pu, pv, pw, d, ge, te, nx, ny, nz,
     &                           iface, jface, kface, ic, jc, kc,
     &                           mass_sum, energy_sum, kin_energy_sum,
     &                           idual, imethod)
c
      implicit none
c
c     Arguments
c
      integer nx, ny, nz, iface, jface, kface, ic, jc, kc, idual
      integer imethod
      real    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      real    d(nx,ny,nz), ge(nx,ny,nz), te(nx,ny,nz)
      real    mass_sum, energy_sum, kin_energy_sum, kin_energy
c
c     Locals
c
      integer i, j, k
      real    mass_term, mom_term, gas_energy
c
      mass_sum = 0.0
      energy_sum = 0.0
      kin_energy_sum = 0.0
c
c     Sum mass and energy      
c
      do k = -1, +2
         do j = -1, +2
            do i = -1, +2
c 1099          format(3(i10,2x),4(g12.3,2x))
c               write(10,1099) i,j,k,d(ic+i,jc+j,kc+k),
c     &              pu(iface+i ,jc   +j, kc+k),
c     &              pv(ic   +i ,jface+j, kc+k   ),
c     &              pw(ic   +i ,jc   +j, kface+k)
               mass_term = d(ic   +i ,jc   +j ,kc+k)
               mom_term = pu(iface+i ,jc   +j, kc+k   )**2 + 
     &                    pv(ic   +i ,jface+j, kc+k   )**2 +
     &                    pw(ic   +i ,jc   +j, kface+k)**2
c
c              Compute total mass and energy
c              
               kin_energy = mom_term / (2.0 * mass_term)

               mass_sum = mass_sum + mass_term

               kin_energy_sum = kin_energy_sum + kin_energy
c
c              Add thermal energy
c              (te is really gas energy if using Zeus)
c
               if (imethod .eq. 2) gas_energy = 
     &              te(ic+i, jc+j, kc+k)*d(ic+i, jc+j, kc+k)
               if (idual .eq. 1) gas_energy = 
     &              ge(ic+i, jc+j, kc+k)*d(ic+i, jc+j, kc+k)
               if (imethod .ne. 2 .and. idual .eq. 0) gas_energy =
     &              te(ic+i, jc+j, kc+k)*d(ic+i, jc+j, kc+k) 
     &              - kin_energy
c
               energy_sum = energy_sum + kin_energy + gas_energy
c
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine sum_abc(pu, pv, pw, d, ge, pu1, pv1, pw1, d1, 
     &                   nx, ny, nz, iface, jface, kface, ic, jc, kc,
     &                   asum, bsum, csum)
c
      implicit none
c
c     Arguments
c
      integer nx, ny, nz, iface, jface, kface, ic, jc, kc
      real    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      real    d(nx,ny,nz), ge(nx,ny,nz)
      real    pu1(4,4,4), pv1(4,4,4), pw1(4,4,4), d1(4,4,4)
      real    asum, bsum, csum
c
c     Locals
c
      integer i, j, k, istart, jstart, kstart
      real    mass_term, mom_term, aterm, bterm, cterm
      real    mass_sum, energy_sum
c
      asum = 0.0
      bsum = 0.0
      csum = 0.0
      mass_sum = 0.0
      energy_sum = 0.0
c
      istart = 2
      jstart = 2
      kstart = 2
c
c     Loop over all affected cells
c
      do k = -1, +2
         do j = -1, +2
            do i = -1, +2
c
               mass_term = d(ic   +i ,jc   +j ,kc+k)
c               write(6,*) "mass_term:",mass_term,
c     &              d1(istart+i,jstart+j,kstart+k)
               mom_term = pu(iface+i ,jc   +j, kc+k   )**2 + 
     &                    pv(ic   +i ,jface+j, kc+k   )**2 +
     &                    pw(ic   +i ,jc   +j, kface+k)**2
c               write(6,*) "mom_term:",mom_term,
c     &              pu(iface+i ,jc   +j, kc+k   ), 
c     &              pv(ic   +i ,jface+j, kc+k   ),
c     &              pw(ic   +i ,jc   +j, kface+k)
c     &              d1(istart+i,jstart+j,kstart+k)

c
c              Do these sums for the quadratic momentum calculation
c                 asum = sum(p_i^2/(2*(m_i+dm_i)))
c                 bsum = 2*sum(p_i*x_i*f_i/(2*(m_i+dm_i)))
c                 csum = sum(x_i*x_i*f_i/(2*(m_i+dm_i)))
c
               mass_term = mass_term + d1(istart+i,jstart+j,kstart+k)
               asum = asum + mom_term/(2.0 * mass_term)
               bterm = 
     &           pu(iface+i,jc+j,kc+k)*pu1(istart+i,jstart+j,kstart+k) +
     &           pv(ic+i,jface+j,kc+k)*pv1(istart+i,jstart+j,kstart+k) +
     &           pw(ic+i,jc+j,kface+k)*pw1(istart+i,jstart+j,kstart+k)
               bsum = bsum + bterm/mass_term ! factors of 2 cancel
               cterm = 
     &           pu1(istart+i,jstart+j,kstart+k)**2 +
     &           pv1(istart+i,jstart+j,kstart+k)**2 +
     &           pw1(istart+i,jstart+j,kstart+k)**2

               
               csum = csum + cterm/(2.0*mass_term)
c               write(6,*) "cterm,mterm:",cterm,mass_term,
c     &              cterm/(2.0*mass_term),csum
c
c               mass_sum = mass_sum + mass_term
c               write(6,*) "mass_sum:",mass_sum
c               energy_sum = energy_sum + mom_term/(2.0*mass_term)
c               energy_sum = energy_sum + ge(ic+1,jc+j,kc+k)
c
            enddo
         enddo
      enddo
c
c      write(6,*) "m,e in sum:", mass_sum, energy_sum
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine add_feedback(pu, pv, pw, d, ge, te, metal, nx, ny, nz, 
     &                        ic, jc, kc, iface, jface, kface,
     &                        dxf, dyf, dzf, dxc, dyc, dzc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf,
     &                        mass_per_cell, mom_per_cell, 
     &                        therm_per_cell)
c
      implicit none
c
c     Arguments
c
      integer nx, ny, nz, ic, jc, kc
      integer iface, jface, kface, imethod, imetal, imulti_metals, idual
      real    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz), te(nx,ny,nz)
      real    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      real    mass_per_cell, mom_per_cell, therm_per_cell
      real    m_eject, yield, metalf
      real    dxf, dyf, dzf, dxc, dyc, dzc,
     &        dxf1, dyf1, dzf1, dxc1, dyc1, dzc1
c
c     Locals
c
      integer i, j, k, i1, j1, k1
      real    delta_mass, delta_pu, delta_pv, delta_pw, delta_therm,
     &        dratio
c
c     Error check
c
      if (imulti_metals .eq. 1) then
         write(6,*) "momentum: not supported"
         ERROR_MESSAGE
      endif
c
c     Loop "cells" in particle-frame
c
      do k = -1, +1
       do j = -1, +1
        do i = -1, +1
c
c            For each particle "cell", do CIC-like deposit
c                  compute zone and face centered weight factors
c
          do i1 = i, i+1
            dxf1 = dxf
            dxc1 = dxc
            if (i1 .eq. i+1) dxf1 = 1.0 - dxf
            if (i1 .eq. i+1) dxc1 = 1.0 - dxc
            do j1 = j, j+1
              dyf1 = dyf
              dyc1 = dyc
              if (j1 .eq. j+1) dyf1 = 1.0 - dyf
              if (j1 .eq. j+1) dyc1 = 1.0 - dyc
              do k1 = k, k+1
                dzf1 = dzf
                dzc1 = dzc
                if (k1 .eq. k+1) dzf1 = 1.0 - dzf
                if (k1 .eq. k+1) dzc1 = 1.0 - dzc
c
                delta_mass = mass_per_cell*dxc1*dyc1*dzc1
c                delme = delme + delta_mass
c                write(6,*) 'delta_mass:',i,j,k,delta_mass
                delta_pu   = i*mom_per_cell*dxf1*dyc1*dzc1
                delta_pv   = j*mom_per_cell*dxc1*dyf1*dzc1
                delta_pw   = k*mom_per_cell*dxc1*dyc1*dzf1
                delta_therm   = therm_per_cell*dxc1*dyc1*dzc1


c                write(6,*) 'delta:',delta_pu,delta_pv,delta_pw,
c     &               delta_therm
c
c               Add mass, momentum
c                (add thermal energy here)
c
                dratio = d(ic+i1, jc+j1, kc+k1)/
     &               (d(ic+i1, jc+j1, kc+k1) + delta_mass)
                d(ic+i1 ,jc+j1 ,kc+k1) = d(ic+i1 ,jc+j1 ,kc+k1)
     &               + delta_mass
                pu(iface+i1 ,jc+j1 ,kc+k1) = 
     &               pu(iface+i1 ,jc+j1 ,kc+k1)+ delta_pu
                pv(ic+i1 ,jface+j1 ,kc+k1) = 
     &               pv(ic+i1 ,jface+j1 ,kc+k1) + delta_pv
                pw(ic+i1 ,jc+j1 ,kface+k1) = 
     &               pw(ic+i1 ,jc+j1 ,kface+k1)+ delta_pw


c
c               Add thermal energy
c
                te(ic+i1 ,jc+j1 ,kc+k1) = 
     &               te(ic+i1 ,jc+j1 ,kc+k1)*dratio + 
     &               delta_therm / d(ic+i1, jc+j1, kc+k1)

                if (idual .eq. 1)
     &               ge(ic+i1 ,jc+j1 ,kc+k1) = 
     &               ge(ic+i1, jc+j1, kc+k1)*dratio +
     &               delta_therm / d(ic+i1, jc+j1, kc+k1)

c
c           Metal feedback (note that in this function gas metal is
c             a fraction (rho_metal/rho_gas) rather than a density.
c             The conversion has been done in the handling routine)
c
                if (imetal .eq. 1) then
c
c           "Cen method".  This takes into account gas recycling.
c                 (metal is not a fraction here)
c
     
                   metal(ic+i1, jc+j1, kc+k1) = 
     &                  metal(ic+i1, jc+j1, kc+k1) + 
     &                  (delta_mass/m_eject)
     &                  * (yield * (1.d0-metalf) + 
     &                  m_eject * metalf)

                   write(6,*) 'metals injected:',(delta_mass/m_eject)
     &                  * (yield * (1.d0-metalf) +
     &                  m_eject * metalf),metal(ic+i1,jc+j1,kc+k1)
c 
                endif
c
c             End loop over CIC-deposit
c
              enddo
            enddo
          enddo
c
c       End loop over "cells" in particle-frame
c
        enddo
       enddo
      enddo
c
      return
      end

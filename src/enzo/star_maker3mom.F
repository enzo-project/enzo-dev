#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker3mom(nx, ny, nz,
     &                      d, dm, temp, u, v, w, cooltime,
     &                      dt, r, metal, zfield1, zfield2,
     &                      dx, t, z, procnum, 
     &                      dunits, x1, vunits, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, mintdyn,
     &                      odthresh, masseff, smthresh, level, np, 
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      imetalSNIa, metalSNIa, metalfSNIa,
     &                      exptime)

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
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      INTG_PREC procnum, imetalSNIa
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz)
      R_PREC    metal(nx,ny,nz), zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      R_PREC    dt, dx, z, exptime
      R_PREC    dunits, x1, vunits, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      R_PREC    metalSNIa(nx,ny,nz), metalfSNIa(nmax)
      R_PREC    odthresh, masseff, smthresh, mintdyn
c
      R_PREC   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      INTG_PREC  i, j, k, ii
      R_PREC   div, tdyn, dtot
      R_PREC   pi, G, sndspdC
      R_PREC   isosndsp2, starmass, starfraction, bmass, jeanmass
      real*8 msolar
      parameter (pi=pi_val, G=GravConst, 
     &           sndspdC=1.3095e8_RKIND,
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
               if (r(i,j,k) .ne. 0._RKIND) goto 10
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
               if (div .ge. 0._RKIND) goto 10
c
c              4) t_cool < t_free-fall (if T < 1.1e4 skip this check)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*dunits
               tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1

               if (tdyn .lt. cooltime(i,j,k) .and. 
     &             temp(i,j,k) .gt. 1.1e4_RKIND) goto 10
c
c              5) is M > M_Jeans? (this definition involves only baryons under
c                 the assumption that the dark matter is stable, which
c                 implies that the dark matter velocity dispersion is >> 
c                 the sound speed.  This will be true for small perturbations
c                 within large halos).
c
               bmass = d(i,j,k)*dble(dunits)*dble(x1*dx)**3 / msolar
               isosndsp2 = sndspdC * temp(i,j,k)
               jeanmass = pi/(6._RKIND*sqrt(d(i,j,k)*dble(dunits))) *
     &                    dble(pi * isosndsp2 / G)**1.5_RKIND / msolar

c
c  THIS IS COMMENTED OUT - NO JEANS MASS CRITERION IN THIS ALGORITHM!!!
c  BWO, 13 NOV 02 (fix 3 dec 02)
c               if (bmass .lt. jeanmass) goto 10
c
c              6) Check to see if star is above threshold (given
c                 in units of M_solar)
c
               starfraction = min(masseff*dt/tdyn, 0.9_RKIND)
               tdyn = max(tdyn, mintdyn*3.15e7_RKIND/t1)

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
                  starfraction = min(smthresh/bmass, 0.5_RKIND)
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
               tcp(ii) = t
               tdp(ii) = tdyn
c              If discrete explosions are used, then use tdp as
c              a flag indicating whether the particle has done
c              feedback rather than dynamical time field
               if (exptime .ge. 0._RKIND) then
                  tdp(ii) = 1._RKIND
               endif
               xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
c
c              Star velocities averaged over multiple cells to
c              avoid "runaway star particle" phenomenon
c              imethod = 2 is zeus, otherwise PPM

               if (imethod .eq. 2) then
                  up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
                  vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
                  wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))
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
                  metalf(ii) = 0._RKIND
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
               d(i,j,k) = (1._RKIND - starfraction)*d(i,j,k)
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
      subroutine star_feedback3mom(nx, ny, nz,
     &               d, mu, dm, te, ge, u, v, w,
     &               metal, zfield1, zfield2,
     &               idual, imetal, imulti_metals, imethod, 
     &               dt, r, dx, t, z,
     &               dunits, x1, vunits, t1, sn_param, m_eject, yield,
c     &               distrad, diststep, distcells,
     &               npart, xstart, ystart, zstart, ibuff,
     &               xp, yp, zp, up, vp, wp,
     &               mp, tdp, tcp, metalf, type, justburn,
     &               kinf_in,exptime_in)

c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Greg Bryan, Christine Simpson
c  date:      July 2015
c  	this is a copy of star_maker3.F that has been modifed to 
c       include the deposition of kinetic as well as thermal energy
c       and metals into a patch surrounding the star particle via a 
c       CIC approach.
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
c    tdp      - dynamical time of zone in which particle created or
c               if exptime >= 0, a flag for whether the discrete 
c               explosion has occurred
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    npart    - particle array size specified by calling routine
c    sn_param - fraction of stellar rest mass that goes to feedback
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c    type     - particle type
c    kinf_fixed - fixed kinetic energy fraction for use when var is
c                 false
c    var      - variable kinetic energy flag
c    exptime  - delay time for discrete explosions (if set to -1,
c               continuous energy injection used)
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c    justburn     - time-weighted mass of star formation (code units)
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, 
     &      imulti_metals, imethod,
     &      distrad, diststep, distcells
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    metal(nx,ny,nz)
      R_PREC    mu(nx,ny,nz)
      R_PREC    zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    dunits, x1, vunits, t1, justburn
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      INTG_PREC type(npart)
      R_PREC 	  kinf,kinf_in,exptime_in,exptime
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      INTG_PREC n, ic, jc, kc, ip, jp, kp,
     &     iface, jface, kface, i, j, k
      R_PREC mform, clight, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, 
     &     dist_mass_cells, dist_mom_cells, mass_per_cell, te_per_cell,
     &     thermal_energy_per_cell, kinetic_energy_per_cell,
     &	   energy_per_cell, mom_per_cell, face_shift,
     &     xface, yface, zface, xpos, ypos, zpos,
     &     energy_before, energy_after, mass_before, mass_after,
     &     kin_energy_before, kin_energy_after,
     &     asum, bsum, csum, delta_ke, ke_injected, ke_after
      R_PREC dxf, dyf, dzf, dxc, dyc, dzc, xfc, yfc, zfc,
     &     xfcshift,yfcshift,zfcshift
      R_PREC u1(4,4,4), v1(4,4,4), w1(4,4,4), d1(4,4,4),
     &     te1(4,4,4), ge1(4,4,4), metal1(4,4,4),
     &     ke_before(4,4,4)
      R_PREC Zsol, d_ave, num_d, energy_51, eunits, R_PDS, t_PDS, 
     &	   R_resolve, KE_rem, e_const, ergs_51, realmass, mu_cell, fbuff
      double precision msolar, mH
      parameter (clight = c_light, msolar_e51 = 1800._RKIND, 
     &           msolar = SolarMass,
     &	         mH = mass_h, e_const = 2.71828_RKIND,
     &	         ergs_51 = 52.57_RKIND)
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
c      write(6,*) 'star_feedback3: start',npart
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2) then
            write(6,*) 'star_feedback3mom'
c
c        The star particle creation algorithm partnered with this 
c        feedback algorithm creates a star particle instantaneously.
c        This feedback routine can either be used to also inject energy
c        instantaneously or to do feedback over a longer period of time
c        as done in star_maker3 modeled on Cen & Ostriker's method that
c        accounts for the true (unsimulated) longer formation time of
c        a stellar population.  The instananeous injection mode is 
c        intended for use with low mass star particles that produce a
c        handful of SN.
c
c     Do instantaneous injeciton
c
            if (exptime_in .ge. 0._RKIND) then
               write(6,*) 'exptime',exptime_in
c              convert explosion time to code units
               exptime = exptime_in * 1.e6_RKIND * 3.15e7_RKIND / t1  
c               write(6,*) 'times: ',abs(t-tcp(n)),exptime+dt,
c     &              t,tcp(n),exptime,dt
               if (tdp(n) .eq. 0._RKIND
     &              .or. abs(t-tcp(n)) .lt. exptime) goto 10
               
               tdp(n) = 0._RKIND
               minitial = mp(n)
               mform = mp(n)
               write(6,*) 'mform, minitial = ',mform,minitial,dx
               goto 5
            endif

c     Determine how much of a given star particle would have been 
c     turned into stars during this timestep.  Then calculate the mass
c     which should have formed during this timestel dt using the integral
c     form of the Cen & Ostriker formula.

            xv1 = (t      - tcp(n))/tdp(n)
c            write(6,*) t,tcp(n),tdp(n),xv1
            if (xv1 .lt. 0._RKIND) goto 10
            if (xv1 .gt. 12._RKIND) goto 10 ! t-tcp >> tdp so ignore
            xv2 = (t + dt - tcp(n))/tdp(n)

c     First calculate the initial mass of the star particle 
c     in question.
            minitial = mp(n) / 
     &           (1._RKIND - 
     &            m_eject*(1._RKIND - (1._RKIND + xv1)*exp(-xv1)))
c     
c     Then, calculate the amount of mass that would have formed in
c     this timestep.
c     
            mform = minitial * ((1._RKIND + xv1)*exp(-xv1) - 
     &           (1._RKIND + xv2)*exp(-xv2))
            mform = max(min(mform, mp(n)), 0._RKIND)
            
c     Compute index of the cell that the star particle
c     resides in.
c 
 5          ip = int((xp(n) - xstart)/dx) + 1
            jp = int((yp(n) - ystart)/dx) + 1
            kp = int((zp(n) - zstart)/dx) + 1
 
c
c          skip if very little mass is formed.
c
            if (mform/d(ip,jp,kp) .lt. 1e-10_RKIND) goto 10
c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
            mp(n) = mp(n) - mform * m_eject
c     
c     Record amount of star formation in this grid.
c
            justburn = justburn + mform * dt * dx**3
c
c     Assuming 3x3x3 cube, calculate cell distribution
c

            distrad = 3
            dist_mass_cells = distrad**3
            dist_mom_cells = distrad**2 * 6
            write(6,*) 'distrad:',distrad
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

            if (xp(n) .lt. xstart .or. xp(n) .gt. xstart+dx*nx .or.
     &          yp(n) .lt. ystart .or. yp(n) .gt. ystart+dx*ny .or.
     &          zp(n) .lt. zstart .or. zp(n) .gt. zstart+dx*ny) then
               write(6,*) 'warning: star particle out of grid',
     &              xp(n),yp(n),zp(n), xstart, ystart, zstart
               goto 100
            endif

c
c	Set center of feedback zone
c
            xfc = xp(n)
            yfc = yp(n)
            zfc = zp(n)
            fbuff = ibuff + 2._RKIND
c
c         check bounds - if star particle is near grid edge
c         then shift center of feedback region
c
            write(6,*) 'buffers:',ibuff,fbuff

            if (xfc .lt. xstart+fbuff*dx .or. 
     &           xfc .gt. xstart+dx*nx-fbuff*dx .or.
     &           yfc .lt. ystart+fbuff*dx .or. 
     &           yfc .gt. ystart+dx*ny-fbuff*dx .or.
     &           zfc .lt. zstart+fbuff*dx .or. 
     &           zfc .gt. zstart+dx*nz-fbuff*dx) then
               write(6,*) 'warning1: star feedback zone shifted',
     &              xfc,yfc,zfc, xstart, ystart, zstart,fbuff
               
	       xfcshift = xfc
	       yfcshift = yfc 
	       zfcshift	= zfc

	       xfc = max(xfc,xstart+fbuff*dx) 
	       yfc = max(yfc,ystart+fbuff*dx) 
	       zfc = max(zfc,zstart+fbuff*dx) 
               
	       xfc = min(xfc,xstart+dx*nx-fbuff*dx) 
	       yfc = min(yfc,ystart+dx*ny-fbuff*dx) 
	       zfc = min(zfc,zstart+dx*nz-fbuff*dx)

	       xfcshift = xfcshift - xfc
	       yfcshift = yfcshift - yfc 
	       zfcshift	= zfcshift - zfc

	       write(6,*) 'warning2: star feedback zone shifted',
     &              xfc,yfc,zfc,dx,xfcshift,yfcshift,zfcshift, 
     &		    nx, ny,nz
            endif
c
c           If using zeus, then velocities are face-centered so shift
c
            face_shift = 0._RKIND
            if (imethod .eq. 2) face_shift = 0.5_RKIND
c
c         Compute index of the first cell to add momentum,
c             accounting for possible face-centering
c
c            write(6,*) 'xp,yp,zp',xp(n), yp(n), zp(n)
c            write(6,*) 'xstart,ystart,zstart', xstart, ystart, zstart
            xface = (xfc - xstart)/dx - 0.5_RKIND - face_shift
            yface = (yfc - ystart)/dx - 0.5_RKIND - face_shift
            zface = (zfc - zstart)/dx - 0.5_RKIND - face_shift
c            write(6,*) 'xface,yface,zface',xface, yface, zface
c
            iface = int(xface + 0.5_RKIND)
            jface = int(yface + 0.5_RKIND)
            kface = int(zface + 0.5_RKIND)
c            write(6,*) 'iface,jface,kface', iface, jface, kface
c
            dxf = real(iface) + 0.5_RKIND - xface
            dyf = real(jface) + 0.5_RKIND - yface
            dzf = real(kface) + 0.5_RKIND - zface
c            write(6,*) 'dxf, dyf, dzf', dxf, dyf, dzf
c
c         Compute index of the first cell to add mass, assuming cell-centering
c 
            xpos = (xfc - xstart)/dx - 0.5_RKIND
            ypos = (yfc - ystart)/dx - 0.5_RKIND
            zpos = (zfc - zstart)/dx - 0.5_RKIND
c            write(6,*) 'xpos, ypos, zpos', xpos, ypos, zpos
c
            ic = int(xpos + 0.5_RKIND)
            jc = int(ypos + 0.5_RKIND)
            kc = int(zpos + 0.5_RKIND)
c            write(6,*) 'ic,jc,kc', ic,jc,kc
c
            dxc = real(ic) + 0.5_RKIND - xpos
            dyc = real(jc) + 0.5_RKIND - ypos
            dzc = real(kc) + 0.5_RKIND - zpos
c            write(6,*) 'dxc, dyc, dzc', dxc, dyc, dzc

c
c     Use fixed kinf value, unless kinf < 0 - then compute
c     variable kinf
            kinf = kinf_in

            if (kinf .lt. 0._RKIND) then
c     Calculate variable kinf based on R_PDS from
c     Cioffi et al. 1988
c
c               realmass = mform * dunits * pc**3 /msolar
c
c               realmass = realmass * (dx* x1 / pc)**3

               realmass = mform * dunits * dx * x1 *dx *x1 *dx *x1
               realmass = realmass /msolar

               energy_51 = sn_param * msolar_e51 * realmass

               Zsol = 0._RKIND
               num_d = 0._RKIND
               d_ave = 0._RKIND
c     
c     Compute average density and metallcity around cell containing
c     particle
c     
               write(6,*) 'averaging over patch'
               do k = -1,1
                  do j = -1,1
                     do i = -1,1
c                        write(6,*) 'Zsol,num_d:',
c     &                       metal(ic+i,jc+j,kc+k)/0.02_RKIND, 
c     &                       d(ic+i,jc+j,kc+k)*
c     &                       dunits/mH/mu(ic+i,jc+j,kc+k) 
c                        write(6,*) 'd,mu:', d(ic+i,jc+j,kc+k)*dunits,
c     &                       mu(ic+i,jc+j,kc+k)
c     mu hardcoded for now, but this should be changed:
c                        mu_cell = 0.82_RKIND 
                        mu_cell = mu(ic+i,jc+j,kc+k)
                        Zsol = Zsol + 
     &                       metal(ic+i,jc+j,kc+k)/0.02_RKIND
                        num_d = num_d + 
     &                       d(ic+i,jc+j,kc+k)*dunits/mu_cell/mH
                        d_ave = d_ave + d(ic+i,jc+j,kc+k)*dunits
                     enddo
                  enddo
               enddo
               
               Zsol = Zsol/27._RKIND
               num_d = num_d/27._RKIND
               d_ave = d_ave/27._RKIND
   
c     compute time and radius of transition to PDS phase for gas
c     with the computed properties
c     Note: t_PDS has units of 1e3yrs and R_PDS has units of pc
c          
c     For metal poor gas
c     
               if (Zsol .lt. 0.01_RKIND) then
                  t_PDS = 3.06e2_RKIND * 
     &                 (energy_51**(1._RKIND/8._RKIND)) *
     &                 (num_d**(-3._RKIND/4._RKIND))
                  R_PDS = 49.3_RKIND * (energy_51**0.25_RKIND) * 
     &                 (num_d**(-0.5_RKIND))
               endif
c
c     For metal rich gas
c
               if (Zsol .ge. 0.01_RKIND) then
                  t_PDS = 26.5_RKIND * 
     &                 (energy_51**(3._RKIND/14._RKIND)) * 
     &                 (Zsol**(-5._RKIND/14._RKIND)) * 
     &                 (num_d**(-4._RKIND/7._RKIND))
                  R_PDS = 18.5_RKIND * energy_51**(2._RKIND/7._RKIND) 
     &                    * num_d**(-3._RKIND/7._RKIND)
     &                    * Zsol**(-1._RKIND/7._RKIND)
               endif
               R_resolve = dx*x1/pc
               if (R_PDS .gt. 4.5*R_resolve) then
                  kinf = 0._RKIND
                  write(6,*) 'Sedov phase resolved:',Zsol,num_d,
     &                 energy_51,R_PDS,R_resolve, realmass
               else                  
                  kinf = 3.97133e-6_RKIND * (d_ave/mH) 
     &                 * (R_resolve**(-2._RKIND)) 
     &                 * (R_PDS**7._RKIND) * (t_PDS**(-2._RKIND)) 
     &                 * (energy_51**(-1._RKIND))
                  write(6,*) 'Sedov phase not resolved:',kinf,Zsol,
     &                 num_d, d_ave, 
     &                 energy_51,t_PDS,R_PDS,R_resolve
               endif
       
            endif

            thermal_energy_per_cell = (1._RKIND - kinf) * 
     &                                energy_per_cell


            write(6,*) 'thermal_energy_per_cell:',
     &           thermal_energy_per_cell,kinf,energy_per_cell
c
c -------- Inject energy, mass & metals----------
c

c
c           Zero local dummy field and kinetic energy field
c
            do k = 1, 4
               do j = 1, 4
                  do i = 1, 4
                     u1(i,j,k) = 0._RKIND
                     v1(i,j,k) = 0._RKIND
                     w1(i,j,k) = 0._RKIND
                     d1(i,j,k) = 0._RKIND
                     ge1(i,j,k) = 0._RKIND
                     te1(i,j,k) = 0._RKIND
                     metal1(i,j,k) = 0._RKIND
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
     &                    iface, jface, kface, imethod, imetal, 
     &                    +1_IKIND)
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
            call add_feedback(u1, v1, w1, d1, ge1, te1, metal1, 4_IKIND, 
     &                        4_IKIND, 4_IKIND, 
     &                        2_IKIND, 2_IKIND, 2_IKIND, 2_IKIND, 
     &                        2_IKIND, 2_IKIND,
     &                        dxf, dyf, dzf, dxc, dyc, dzc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf(n),
     &                        mass_per_cell, 1._RKIND, 0._RKIND)
c
c           The kinetic energy after the feedback event is a quadratic equation
c           with delta_p, where delta_p is the total amount of momentum added
c           to the feedback region:
c           E_k,a = asum + bsum*delta_p + csum*delta_p*delta_p
c
c           Sum a, b and c terms to get momentum normalization
c
            if (kinf .gt. 0._RKIND) then
               call sum_abc(u, v, w, d, ge, u1, v1, w1, d1, 
     &              nx, ny, nz, iface, jface, kface, ic, jc, kc,
     &              asum, bsum, csum)
               write(6,*) 'abc1',asum, bsum, csum
               asum = asum - (kin_energy_before + kinf*energy)
               write(6,*) 'kin_energy_before,kinf*energy,energy:',
     &              kin_energy_before,kinf*energy,energy
               write(6,*) 'abc2',asum, bsum, csum
c
c           Calculate momentum contribution
c
               mom_per_cell = (-bsum + 
     &              sqrt(bsum**2._RKIND - 4._RKIND*asum*csum))/
     &              (2._RKIND*csum)
            else
c
c           If the kinetic fraction is set to 0, then 0 momentum is added
c           to each cell.
c
               mom_per_cell = 0._RKIND
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
     &                    iface, jface, kface, imethod, imetal,
     &                    -1_IKIND)
c
c           Add the increase in the kinetic energy to the total energy field
c             (unless we're using Zeus).  If using dual energy formalism, we
c             might want to enforce consistency.
c     
            if (imethod .ne. 2) then
               ke_injected = 0._RKIND
               do k = -1, +2
                  do j = -1, +2
                     do i = -1, +2

                        ke_after = 0.5_RKIND*
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
     &                    ke_injected) .gt. 0.01_RKIND .and. 
     &              kinf .ne. 0._RKIND) then
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
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      INTG_PREC iface, jface, kface, imethod, imetal, idir
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    up, vp, wp
c
c     Locals
c
      INTG_PREC i, j, k
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
               if (idir. eq. +1_IKIND) then
                  if (imethod .eq. 2) then
                     u(iface+i ,jc+j ,kc+k) = 
     &                            (u(iface+i ,jc+j ,kc+k) - up)
     &                    * 0.5_RKIND * (d(iface+i  , jc+j, kc+k) +
     &                             d(iface+i+1, jc+j, kc+k))
                     v(ic+i ,jface+j ,kc+k) = 
     &                            (v(ic+i ,jface+j ,kc+k) - vp)
     &                    * 0.5_RKIND * (d(ic+i, jface+j  , kc+k) + 
     &                             d(ic+i, jface+j+1, kc+k))
                     w(ic+i ,jc+j ,kface+k) = 
     &                            (w(ic+i ,jc+j ,kface+k) - wp)
     &                    * 0.5_RKIND * (d(ic+i, jc+j, kface+k  ) + 
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
     &                    /( 0.5_RKIND * (d(iface+i  , jc+j, kc+k) +
     &                              d(iface+i+1, jc+j, kc+k))) + up
                     v(ic+i ,jface+j ,kc+k) = v(ic+i ,jface+j ,kc+k) 
     &                    /( 0.5_RKIND * (d(ic+i, jface+j  , kc+k) + 
     &                              d(ic+i, jface+j+1, kc+k))) + vp
                     w(ic+i ,jc+j ,kface+k) = w(ic+i ,jc+j ,kface+k) 
     &                    /( 0.5_RKIND * (d(ic+i, jc+j, kface+k  ) + 
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
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, iface, jface, kface, ic, jc, kc, idual
      INTG_PREC imethod
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      R_PREC    d(nx,ny,nz), ge(nx,ny,nz), te(nx,ny,nz)
      R_PREC    mass_sum, energy_sum, kin_energy_sum, kin_energy
c
c     Locals
c
      INTG_PREC i, j, k
      R_PREC    mass_term, mom_term, gas_energy
c
      mass_sum = 0._RKIND
      energy_sum = 0._RKIND
      kin_energy_sum = 0._RKIND
c
c     Sum mass and energy      
c
      do k = -1, +2
         do j = -1, +2
            do i = -1, +2
               mass_term = d(ic   +i ,jc   +j ,kc+k)
               mom_term = pu(iface+i ,jc   +j, kc+k   )**2 + 
     &                    pv(ic   +i ,jface+j, kc+k   )**2 +
     &                    pw(ic   +i ,jc   +j, kface+k)**2
c
c              Compute total mass and energy
c              
               kin_energy = mom_term / (2._RKIND * mass_term)

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
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, iface, jface, kface, ic, jc, kc
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      R_PREC    d(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu1(4,4,4), pv1(4,4,4), pw1(4,4,4), d1(4,4,4)
      R_PREC    asum, bsum, csum
c
c     Locals
c
      INTG_PREC i, j, k, istart, jstart, kstart
      R_PREC    mass_term, mom_term, aterm, bterm, cterm
      R_PREC    mass_sum, energy_sum
c
      asum = 0._RKIND
      bsum = 0._RKIND
      csum = 0._RKIND
      mass_sum = 0._RKIND
      energy_sum = 0._RKIND
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

               
               csum = csum + cterm/(2._RKIND*mass_term)
c               write(6,*) "cterm,mterm:",cterm,mass_term,
c     &              cterm/(2.0*mass_term),csum
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
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, idual
      INTG_PREC iface, jface, kface, imethod, imetal, imulti_metals
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      R_PREC    mass_per_cell, mom_per_cell, therm_per_cell
      R_PREC    m_eject, yield, metalf
      R_PREC    dxf, dyf, dzf, dxc, dyc, dzc,
     &        dxf1, dyf1, dzf1, dxc1, dyc1, dzc1
c
c     Locals
c
      INTG_PREC i, j, k, i1, j1, k1
      R_PREC    delta_mass, delta_pu, delta_pv, delta_pw, delta_therm,
     &        dratio, tot_mass

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
      
      tot_mass = 0._RKIND

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
c                write(6,*) i,j,k,dxc1*dyc1*dzc1
c                write(6,*) ic+i1,jc+j1,kc+k1
c                write(6,*) nx-ic-i1,ny-jc-j1,nz-kc-k1
c                write(6,*) '   '
                

                dratio = d(ic+i1, jc+j1, kc+k1)/(d(ic+i1, jc+j1, kc+k1) 
     &               + delta_mass)
                d(ic+i1 ,jc+j1 ,kc+k1) = d(ic+i1 ,jc+j1 ,kc+k1)
     &               + delta_mass
                pu(iface+i1 ,jc+j1 ,kc+k1) = 
     &               pu(iface+i1 ,jc+j1 ,kc+k1)+ delta_pu
                pv(ic+i1 ,jface+j1 ,kc+k1) = 
     &               pv(ic+i1 ,jface+j1 ,kc+k1) + delta_pv
                pw(ic+i1 ,jc+j1 ,kface+k1) = 
     &               pw(ic+i1 ,jc+j1 ,kface+k1)+ delta_pw

                tot_mass = tot_mass + delta_mass
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
     &                  * (yield * (1._RKIND-metalf) + 
     &                  m_eject * metalf)

c                   write(6,*) 'metals injected:',(delta_mass/m_eject)
c     &                  * (yield * (1._RKIND-metalf) +
c     &                  m_eject * metalf),metal(ic+i1,jc+j1,kc+k1)
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

c      write(6,*) 'tot_mass:',tot_mass
c
      return
      end

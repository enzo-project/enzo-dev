#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker3mech(nx, ny, nz,
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
c  CREATES STAR PARTICLES FOR KINETIC FEEDBACK
c
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
      R_PREC   sndspdC
      R_PREC   isosndsp2, starmass, starfraction, bmass, jeanmass
      parameter (sndspdC=1.3095e8_RKIND)
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
               tdyn  = sqrt(3._RKIND*pi_val/32._RKIND/
     &                      GravConst/dtot)/t1

               if (tdyn .lt. cooltime(i,j,k) .and. 
     &             temp(i,j,k) .gt. 1.1e4_RKIND) goto 10
c
c              5) is M > M_Jeans? (this definition involves only baryons under
c                 the assumption that the dark matter is stable, which
c                 implies that the dark matter velocity dispersion is >> 
c                 the sound speed.  This will be true for small perturbations
c                 within large halos).
c
               bmass = d(i,j,k)*dble(dunits)*dble(x1*dx)**3 / SolarMass
               isosndsp2 = sndspdC * temp(i,j,k)
               jeanmass = pi_val/(6._RKIND*
     &              sqrt(d(i,j,k)*dble(dunits))) *
     &              dble(pi_val * isosndsp2 /
     &              GravConst)**1.5_RKIND / SolarMass

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
      subroutine star_feedback3mech(nx, ny, nz,
     &               d, mu, dm, te, ge, u, v, w,
     &               metal, zfield1, zfield2,
     &               idual, imetal, imulti_metals, imethod, 
     &               dt, r, dx, t, z,
     &               dunits, x1, vunits, t1, thermal, m_eject, yield,
     &               npart, xstart, ystart, zstart, ibuff,
     &               xp, yp, zp, up, vp, wp,
     &               mp, tdp, tcp, metalf, type, justburn,
     &               exptime_in)

c
c  RELEASES "STAR" PARTICLE MOMENTUM, MASS AND METALS
c
c  written by: Cassi Lochhaas, Claire Kopenhafer
c  date:      December 2023
c  	this is a copy of star_maker3mom.F that has been modifed to 
c       calculate the amount of momentum to inject from properties
c       of the host and surrounding cells, following Kimm & Cen (2014)
c       and Kimm et al. (2015).
c       The momentum is then injected with the CIC approach.
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
c    thermal  - thermal energy to deposit in addition to momentum
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c    type     - particle type
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
#define NO_FORTRAN_DEBUG
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, 
     &      imulti_metals, imethod
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
      R_PREC 	 exptime_in,exptime
      R_PREC    thermal
c
c  Locals
c
      INTG_PREC n, ic, jc, kc, ip, jp, kp,
     &     i, j, k
      R_PREC mform, energy,
     &     m_eject, yield, minitial, xv1, xv2, distrad,
     &     dist_mass_cells, mass_per_cell,
     &     mass_before(3,3,3), 
     &     mass_after(3,3,3), kin_energy_before(3,3,3),
     &     kin_energy_after(3,3,3), kin_energy_inj(3,3,3),
     &     delta_ke, ke_injected, ke_after,
     &     nsn, mass_end(3,3,3), kin_energy_end(3,3,3), 
     &     mom_start(3,3,3), mom_inj(3,3,3), mom_end(3,3,3)
      R_PREC xfc, yfc, zfc,
     &     xfcshift,yfcshift,zfcshift
      R_PREC ke_before(3,3,3), mom_per_cell(3,3,3), d_init(nx,ny,nz),
     &     metal_init(nx,ny,nz)
      R_PREC Z_avg, n_avg,
     &	  realmass, mu_cell, fbuff, therm_per_cell
      R_PREC te_sum, ke_sum, ge_sum, mom_sum, mass_sum, mom_term,
     &     mass_term
c
c-----------------------------------------------------------------------
c
c     Make copies of density and metal fields for calculating momentum deposition
c     without the ejected mass of nearby particles in the same timestep
c
      do k = 1_IKIND, nz
         do j = 1_IKIND, ny
            do i = 1_IKIND, nx
               d_init(i,j,k) = d(i,j,k)
               metal_init(i,j,k) = metal(i,j,k)
            enddo
         enddo
      enddo
c
c     Loop over particles
c
c      write(6,*) 'star_feedback3: start',npart
      do n=1_IKIND, npart
         if (tcp(n) .gt. 0._RKIND .and. mp(n) .gt. 0._RKIND .and. 
     &        type(n) .eq. 2_IKIND) then
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
c              convert explosion time to code units
               exptime = exptime_in * 1.e6_RKIND * 3.15e7_RKIND / t1  
#ifdef FORTRAN_DEBUG
               write(6,*) 'mform, minitial = ',mform,minitial
               write(6,*) 'exptime,t,tcp,tdp = ',exptime,t,tcp(n),tdp(n)
#endif
               if (tdp(n) .eq. 0._RKIND
     &              .or. abs(t-tcp(n)) .lt. exptime) goto 10
               tdp(n) = 0._RKIND
               minitial = mp(n)
               mform = mp(n)
               goto 5
            endif

c     Determine how much of a given star particle would have been 
c     turned into stars during this timestep.  Then calculate the mass
c     which should have formed during this timestel dt using the integral
c     form of the Cen & Ostriker formula.

            xv1 = (t      - tcp(n))/tdp(n)
#ifdef FORTRAN_DEBUG
            write(6,*) t,tcp(n),tdp(n),xv1
#endif
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
 5          ip = int((xp(n) - xstart)/dx) + 1_IKIND
            jp = int((yp(n) - ystart)/dx) + 1_IKIND
            kp = int((zp(n) - zstart)/dx) + 1_IKIND
#ifdef FORTRAN_DEBUG
            write(6,*) 'x1,dx,dunits,vunits:',x1,dx,dunits,vunits
            write(6,*) 'mp(n):',mp(n)*dunits*x1**3._RKIND*dx**3_IKIND
     &        /SolarMass
            write(6,*) 'minitial:',minitial*dunits*x1**3._RKIND
     &        /SolarMass*dx**3_IKIND
            write(6,*) 'mform:',mform*dunits*x1**3._RKIND*dx**3_IKIND
     &        /SolarMass
#endif
c
c     Assuming 3x3x3 cube, calculate cell distribution
c
            distrad = 3_IKIND
            dist_mass_cells = distrad**3_IKIND
c
c          skip if very little mass is formed.
c
            if (mform/dist_mass_cells/d(ip,jp,kp) .lt. 1e-10_RKIND) 
     &           goto 10
c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
            mp(n) = mp(n) - mform * m_eject
c     
c     Record amount of star formation in this grid.
c
            justburn = justburn + mform * dt * dx**3_IKIND
c     
c           Calculate mass per cell ejected
c     
            mass_per_cell = mform * m_eject / dist_mass_cells
#ifdef FORTRAN_DEBUG
            write(6,*) 'mass_per_cell:',mass_per_cell*
     &         dunits*x1**3._RKIND/SolarMass*dx**3_IKIND
#endif
            if (xp(n) .lt. xstart .or. xp(n) .gt. xstart+dx*nx .or.
     &          yp(n) .lt. ystart .or. yp(n) .gt. ystart+dx*ny .or.
     &          zp(n) .lt. zstart .or. zp(n) .gt. zstart+dx*nz) then
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
            if (xfc .lt. xstart+fbuff*dx .or. 
     &          xfc .gt. xstart+dx*nx-fbuff*dx .or.
     &          yfc .lt. ystart+fbuff*dx .or. 
     &          yfc .gt. ystart+dx*ny-fbuff*dx .or.
     &          zfc .lt. zstart+fbuff*dx .or. 
     &          zfc .gt. zstart+dx*nz-fbuff*dx) then
#ifdef FORTRAN_DEBUG
               write(6,*) 'warning1: star feedback zone shifted',
     &              xfc,yfc,zfc, xstart, ystart, zstart,fbuff
#endif
c               
	       xfcshift = xfc
	       yfcshift = yfc 
	       zfcshift	= zfc
c
	       xfc = max(xfc,xstart+fbuff*dx) 
	       yfc = max(yfc,ystart+fbuff*dx) 
	       zfc = max(zfc,zstart+fbuff*dx) 
c              
	       xfc = min(xfc,xstart+dx*nx-fbuff*dx) 
	       yfc = min(yfc,ystart+dx*ny-fbuff*dx) 
	       zfc = min(zfc,zstart+dx*nz-fbuff*dx)
c
	       xfcshift = xfcshift - xfc
	       yfcshift = yfcshift - yfc 
	       zfcshift	= zfcshift - zfc
c
#ifdef FORTRAN_DEBUG
	       write(6,*) 'warning2: star feedback zone shifted',
     &              xfc,yfc,zfc,dx,xfcshift,yfcshift,zfcshift, 
     &		    nx, ny,nz
#endif
            endif
c
c         Compute index of the center of feedback zone
c
            ic = int((xfc - xstart)/dx) + 1_RKIND
            jc = int((yfc - ystart)/dx) + 1_RKIND
            kc = int((zfc - zstart)/dx) + 1_RKIND
c
c           Calculate the number of SNe in this timestep from the stars formed
c           in this timestep as 1 SN per 100 Msun formed.
c           
            realmass = mform * dunits * x1**3._RKIND*dx**3_IKIND
            realmass = realmass /SolarMass
            nsn = realmass / 100._RKIND
            energy = nsn * 10._RKIND**51._RKIND / (dunits*
     &                  x1**3._RKIND*vunits**2._RKIND*dx**3._RKIND)
            therm_per_cell = thermal / dist_mass_cells
#ifdef FORTRAN_DEBUG
            write(6,*) 'realmass,nsn:',realmass,nsn
            write(6,*) 'ic,jc,kc:',ic, jc, kc
#endif
c
            Z_avg = 0._RKIND
            n_avg = 0._RKIND
c     
c     Compute average density and metallcity around cell containing
c     particle
c     
c            write(6,*) 'averaging over patch'
            do k = -1_IKIND,1_IKIND
               do j = -1_IKIND,1_IKIND
                  do i = -1_IKIND,1_IKIND
                     mu_cell = mu(ic+i,jc+j,kc+k)
                     Z_avg = Z_avg + 
     &                 metal_init(ic+i,jc+j,kc+k)/0.02_RKIND
                     n_avg = n_avg + 
     &                 d_init(ic+i,jc+j,kc+k)*dunits*
     &                 mu_cell/mass_h
                  enddo
               enddo
            enddo
c               
            Z_avg = max(Z_avg/27._RKIND, 0.01)
            n_avg = n_avg/27._RKIND
#ifdef FORTRAN_DEBUG
            write(6,*) 'n_avg, Z_avg:', n_avg, Z_avg
#endif
c
c     Inject energy, mass & metals
c
c     Compute the momentum that will be added to each cell.
c
            call get_mom_per_cell(d_init, mass_per_cell,
     &                  dist_mass_cells, mom_per_cell,
     &                  nsn, d_init(ic,jc,kc), dx, x1, n_avg, Z_avg,
     &                  dunits, vunits, t1,
     &                  nx, ny, nz, ic, jc, kc)
#ifdef FORTRAN_DEBUG
            write(6,*) 'mom_per_cell:',mom_per_cell*
     &           (dunits*(dx*x1)**3._RKIND)*vunits/SolarMass/1e5_RKIND
#endif
c           mom_per_cell has units of mass*velocity/volume
c
c     Compute the kinetic energy in the affected region
c     before momentum is added (except for ZEUS).  This
c     is needed at the end of the calculation to update the
c     total energy (te) field.
c
            if (imethod .ne. 2_IKIND) then
               do k = -1_IKIND, +1_IKIND
                  do j = -1_IKIND, +1_IKIND
                     do i = -1_IKIND, +1_IKIND
                        ke_before(i+2, j+2, k+2) = 
     &                       0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                                (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                 v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                 w(ic+i ,jc+j ,kc+k)**2_IKIND)
                     enddo
                  enddo
               enddo
            endif
#ifdef FORTRAN_DEBUG
            write(6,*) 'ke_before:',ke_before*dunits*x1**3._RKIND*
     &                  vunits**2._RKIND*dx**3._RKIND
#endif
c           ke_before has units of mass*velocity^2/volume
c
c     First convert velocities to momenta and transform
c     into frame comoving with particle
c
            call momentum_convert(u, v, w, d, metal, up(n), vp(n),
     &                    wp(n), nx, ny, nz, ic, jc, kc, 
     &                    imethod, imetal, 
     &                    +1_IKIND, dx, dunits, vunits, x1)
c
c           u,v,w have units of mass*velocity/volume
c           Now u,v,w are the particle-relative *momentum* directions
c
#ifdef FORTRAN_DEBUG
c
c           Calculate mass and energy before in affected region
c
            call mass_kinetic_energy(u, v, w, d, nx, ny,
     &           nz, ic, jc, kc, mass_before, kin_energy_before,
     &           mom_start)
            write(6,*) "m before:", mass_before*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
            write(6,*) "e before:", kin_energy_before*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
            write(6,*) 'mom_start:',mom_start*dunits*vunits
     &                  *(dx*x1)**3._RKIND/SolarMass/1e5_RKIND
c           mass_before has units of mass/volume
c           kin_energy_before has units of mass*velocity^2/volume
c           mom_start has units of mass*velocity/volume
c
#endif
c
c           Now deposit the momentum onto the grid.
c
            call add_feedback_mech(u, v, w, d, ge, te, metal,
     &                         nx, ny, nz, ic, jc, kc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf(n),
     &                        mass_per_cell, mom_per_cell, 
     &                        therm_per_cell, dx, x1, dunits, vunits)
#ifdef FORTRAN_DEBUG
c
c           Calculate the mass and energy in the affected region again.
c
            call mass_kinetic_energy(u, v, w, d, nx, ny,
     &                          nz, ic, jc, kc,
     &                          mass_end, kin_energy_end, mom_end)
            write(6,*) "m end:", mass_end*
     &               (dunits*(dx*x1)**3._RKIND)/SolarMass
            write(6,*) 'kin_energy_end:',kin_energy_end*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
            write(6,*) 'mom_end:',mom_end*dunits*vunits
     &                  *(dx*x1)**3._RKIND/SolarMass/1e5_RKIND
#endif
c
c           Convert momenta back to velocities and transform back to lab frame
c
            call momentum_convert(u, v, w, d, metal, up(n), vp(n), 
     &                    wp(n), nx, ny, nz, ic, jc, kc, 
     &                    imethod, imetal,
     &                    -1_IKIND, dx, dunits, vunits, x1)
c
c           Add the increase in the kinetic energy to the total energy field
c             (unless we're using Zeus).  If using dual energy formalism, we
c             might want to enforce consistency.
c     
            if (imethod .ne. 2_IKIND) then
               ke_injected = 0._RKIND
               do k = -1_IKIND, +1_IKIND
                  do j = -1_IKIND, +1_IKIND
                     do i = -1_IKIND, +1_IKIND
                        ke_after = 0.5_RKIND *
     &                        d(ic+i, jc+j, kc+k) *
     &                       (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                        v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                        w(ic+i ,jc+j ,kc+k)**2_IKIND) 
                        delta_ke = ke_after - ke_before(i+2,j+2,k+2)
                        te(ic+i ,jc+j ,kc+k) = 
     &                       te(ic+i ,jc+j ,kc+k) + 
     &                       delta_ke/d(ic+i, jc+j, kc+k)
                        ke_injected = ke_injected + delta_ke
c                    ke_after has units of mass*velocity^2/volume
c                    te has units of mass*velocity^2/volume/mass*volume = velocity^2
c                    te is *specific* energy (energy/mass)
                     enddo
                  enddo
               enddo
            endif
#ifdef FORTRAN_DEBUG
c
c           Sum and print the momentum, kinetic energy, thermal energy, and
c           total energy in a large region surrounding the particle
c
            te_sum = 0._RKIND
            ke_sum = 0._RKIND
            ge_sum = 0._RKIND
            mom_sum = 0._RKIND
            mass_sum = 0._RKIND
            do k = -5_IKIND, +5_IKIND
               do j = -5_IKIND, +5_IKIND
                  do i = -5_IKIND, +5_IKIND
                     mass_term = d(ic+i, jc+j, kc+k)
                     mom_term = u(ic+i, jc+j, kc+k)**2._RKIND + 
     &                    v(ic+i, jc+j, kc+k)**2._RKIND + 
     &                    w(ic+i, jc+j, kc+k)**2._RKIND
                     ke_sum = ke_sum + 0.5_RKIND*(mom_term*mass_term)
                     mom_sum = mom_sum + mass_term*mom_term**0.5_RKIND
                     mass_sum = mass_sum + mass_term
                     te_sum = te_sum + te(ic+i, jc+j, kc+k) * mass_term
                     ge_sum = ge_sum + ge(ic+i, jc+j, kc+k) * mass_term
c           mass_sum has units mass/volume
c           mom_sum has units mass*velocity/volume
c           ke_sum has units mass*velocity^2/volume
c           te_sum has units mass*velocity^2/volume
c           ge_sum has units mass*velocity^2/volume
                  enddo
               enddo
            enddo
            write(6,*) 'mass, mom, KE, GE, TE:',
     &        mass_sum*dunits*(x1*dx)**3._RKIND/SolarMass,
     &        mom_sum*dunits*(x1*dx)**3._RKIND*vunits/
     &        SolarMass/1e5_RKIND,
     &        ke_sum*dunits*(x1*dx)**3._RKIND*vunits**2._RKIND,
     &        ge_sum*dunits*(x1*dx)**3._RKIND*vunits**2._RKIND,
     &        te_sum*dunits*(x1*dx)**3._RKIND*vunits**2._RKIND
#endif
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
c     Calculate the amount of momentum to add to each neighboring cell
c
      subroutine get_mom_per_cell(d_init, mass_per_cell, 
     &               dist_mass_cells, mom_per_cell,
     &               nsn, d_fbck_cell, dx, x1,
     &               n_avg, Z_avg, dunits, vunits, t1,
     &               nx, ny, nz, ic, jc, kc)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      R_PREC    d_init(nx,ny,nz), dx
      R_PREC    t1, Z_avg, n_avg
      R_PREC    mom_per_cell(3_IKIND, 3_IKIND, 3_IKIND)
      R_PREC    mass_per_cell, dist_mass_cells
      R_PREC    nsn, x1, dunits, d_fbck_cell, vunits
c
c     Locals
c
      INTG_PREC i, j, k, dist
      R_PREC    dMswept, chi, chi_th, fe, pSN, ergs_51
      R_PREC    den_cell, Mej, dMej, Nnbors, betaSN
c
        ergs_51 = 1.e51_RKIND
        Nnbors = dist_mass_cells - 1._RKIND
        betaSN = 1._RKIND/dist_mass_cells
        Mej = mass_per_cell*dist_mass_cells*dunits*(x1*dx)**3._RKIND
        chi_th = 69.58_RKIND*(nsn)**(-2._RKIND/17._RKIND) * 
     &               (n_avg)**(-4._RKIND/17._RKIND) * 
     &               (Z_avg)**(-0.28_RKIND)
#ifdef FORTRAN_DEBUG
        write(6,*) 'chi_th:',chi_th
#endif
        do k = -1_IKIND, +1_IKIND
            do j = -1_IKIND, +1_IKIND
               do i = -1_IKIND, +1_IKIND
                  dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
                  if (dist .eq. 0_IKIND) then
                     pSN = 0._RKIND
                  else
c
c       Compute the mass fraction of swept-up to ejected mass.
c
                     den_cell = d_init(ic+i,jc+j,kc+k)*dunits
                     dMej = (1._RKIND-betaSN)*Mej/Nnbors
                     dMswept = den_cell * 
     &                     (dx*x1/2._RKIND)**(3._RKIND) + 
     &                     ((1._RKIND-betaSN)*
     &                     (d_fbck_cell*dunits)*(dx*x1)**3._RKIND) /
     &                     Nnbors + dMej
c
                     chi = dMswept/dMej
#ifdef FORTRAN_DEBUG
                     write(6,*) 'chi:',chi
#endif
c
c       Compare chi to threshold and set momentum injection.
c
                     if (chi .lt. chi_th) then
                        fe = 1._RKIND - ((chi - 1._RKIND)/
     &                        (3._RKIND*(chi_th - 1._RKIND)))
                        pSN = (2._RKIND*chi * Mej * fe *
     &                        nsn * ergs_51)**(0.5_RKIND)
#ifdef FORTRAN_DEBUG
                        write(6,*) 'lt chi_th, pSN:',
     &                        pSN/SolarMass/1e5_RKIND
#endif
                     else
                        pSN = 3.e10_RKIND*SolarMass * 
     &                        (nsn)**(16._RKIND/17._RKIND) * 
     &                        (n_avg)**(-2._RKIND/17._RKIND) * 
     &                        (Z_avg)**(-0.14_RKIND)
#ifdef FORTRAN_DEBUG
                        write(6,*) 'gt chi_th, pSN:',
     &                        pSN/SolarMass/1e5_RKIND
#endif
                     endif
c              pSN has units of cgs-mass*cgs-velocity, want it to have units of code-mass*code-velocity/code-volume
c              conversion of mass from code units is (dunits*(dx*x1)**3._RKIND) (density times volume)
                     pSN = pSN/(dunits*(dx*x1)**3._RKIND)/vunits/Nnbors
c              dunits = cgs-density/code-density = cgs-mass*code-volume/(cgs-volume*code-mass)
c              x1^3 = cgs-volume/code-volume
c              vunits = cgs-velocity/code-velocity
c              dx^3 = code-volume
c              dunits*(dx*x1)**3._RKIND = cgs-mass*code-volume/(cgs-volume*code-mass)*code-volume*cgs-volume/code-volume
c              dunits*(dx*x1)**3._RKIND = cgs-mass*code-volume/code-mass
c              pSN = cgs-mass*cgs-velocity/(cgs-mass*code-volume/code-mass)/(cgs-velocity/code-velocity)
c              pSN = code-mass*code-velocity/code-volume
c              pSN now has units code-mass*code-velocity/code-volume
c              So pSN is actually in the right units after all!
c
c        If cell is on edge or corner of injection cube, reduce the
c        amount of momentum that gets added to each cartesian direction
c        to avoid double-counting.
c
                     if (dist .eq. 2_IKIND) then
                        pSN = pSN/(2._RKIND**0.5_RKIND)
                     else if (dist .eq. 3_IKIND) then
                        pSN = pSN/(3._RKIND**0.5_RKIND)
                     endif
                  endif
                  mom_per_cell(i+2,j+2,k+2) = pSN
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
c     Convert velocities to momentum and back
c
      subroutine momentum_convert(u, v, w, d, metal, up, vp, wp, 
     &                    nx, ny, nz, ic, jc, kc,
     &                    imethod, imetal, idir, dx, dunits, vunits, x1)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      INTG_PREC imethod, imetal, idir
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    up, vp, wp, dx
      R_PREC    dunits, vunits, x1
c
c     Locals
c
      INTG_PREC i, j, k
c
c     Error check
c
      if (idir .ne. -1_IKIND .and. idir .ne. 1_IKIND) then
         write(6,*) 'incorrect idir value in momentum call'
         ERROR_MESSAGE
      endif
c
c     Loop over velocities, multiplying by densities (or dividing if
c       converting back)
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
c
c              idir = +1: convert vel -> mom
c               
               if (idir. eq. +1_IKIND) then
                  if (imethod .eq. 2_IKIND) then
                     write(6,*) 'Zeus not supported!'
                     ERROR_MESSAGE
                  else
                     u(ic+i ,jc+j ,kc+k) = (u(ic+i ,jc+j ,kc+k)-up) *
     &                       d(ic+i, jc+j, kc+k)
                     v(ic+i ,jc+j ,kc+k) = (v(ic+i ,jc+j ,kc+k)-vp) *
     &                       d(ic+i, jc+j, kc+k)
                     w(ic+i ,jc+j ,kc+k) = (w(ic+i ,jc+j ,kc+k)-wp) *
     &                       d(ic+i, jc+j, kc+k)
                  endif
                  if (imetal .eq. 1_IKIND) then
                     metal(ic+i,jc+j,kc+k) = metal(ic+i,jc+j,kc+k) *
     &                    d(ic+i,jc+j,kc+k)
                  endif
c
c              if idir = -1: convert mom -> vel
c               
               else
                  if (imethod .eq. 2_IKIND) then
                     write(6,*) 'Zeus not supported!'
                     ERROR_MESSAGE
                  else
                     u(ic+i ,jc+j ,kc+k) = u(ic+i ,jc+j ,kc+k) /
     &                    (d(ic+i, jc+j, kc+k)) + up
                     v(ic+i ,jc+j ,kc+k) = v(ic+i ,jc+j ,kc+k) /
     &                    (d(ic+i, jc+j, kc+k)) + vp
                     w(ic+i ,jc+j ,kc+k) = w(ic+i ,jc+j ,kc+k) /
     &                    (d(ic+i, jc+j, kc+k)) + wp
                  endif
                  if (imetal .eq. 1_IKIND) then
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
      subroutine mass_kinetic_energy(pu, pv, pw, d,
     &                           nx, ny, nz,
     &                           ic, jc, kc,
     &                           mass_grid, kin_energy_grid, mom_grid)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      R_PREC    d(nx,ny,nz)
      R_PREC    mass_grid(3,3,3), kin_energy_grid(3,3,3)
      R_PREC    mom_grid(3,3,3)
c
c     Locals
c
      INTG_PREC i, j, k
      R_PREC    mass_term, mom_term, kin_energy
c
c
c     Sum mass and energy      
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
               mass_term = d(ic   +i ,jc   +j ,kc   +k)
               mom_term = pu(ic   +i ,jc   +j, kc   +k)**2_IKIND + 
     &                    pv(ic   +i ,jc   +j, kc   +k)**2_IKIND + 
     &                    pw(ic   +i ,jc   +j, kc   +k)**2_IKIND
c              mass_term has units of mass/volume
c              mom_term has units of mass^2*velocity^2/volume^2
c
c              Compute total mass and energy
c              
               kin_energy = mom_term / (2._RKIND * mass_term)
c              kin_energy has units of mass^2*velocity^2/volume^2 / (mass/volume)
c                          = mass*velocity^2/volume
               mass_grid(i+2,j+2,k+2) = mass_term
               kin_energy_grid(i+2,j+2,k+2) = kin_energy
               mom_grid(i+2,j+2,k+2) = mom_term**0.5_RKIND
c              mom_grid has units of mass*velocity/volume
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Add mass, momentum, and energy to affected cells.
c       Note that pu, pv, pw are momenta
c
      subroutine add_feedback_mech(pu, pv, pw, d, ge, te, metal, 
     &                        nx, ny, nz, ic, jc, kc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf,
     &                        mass_per_cell, mom_per_cell, 
     &                        therm_per_cell, dx, x1, dunits, vunits)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, idual
      INTG_PREC imethod, imetal, imulti_metals
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      R_PREC    mass_per_cell, mom_per_cell(3_IKIND,3_IKIND,3_IKIND)
      R_PREC    m_eject, yield, metalf, therm_per_cell
      R_PREC    dunits, x1, vunits, dx
c
c     Locals
c
      INTG_PREC i, j, k
      R_PREC    delta_mass, delta_pu, delta_pv, delta_pw, delta_therm,
     &          dratio
      R_PREC    pu_start, pu_deposit, pv_start, pv_deposit, pw_start,
     &          pw_deposit
c
c     Error check
c
      if (imulti_metals .eq. 1_IKIND) then
         write(6,*) "momentum: not supported"
         ERROR_MESSAGE
      endif
c
c     Loop cells in particle-frame
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
c
               delta_mass =   mass_per_cell
c              delta_mass has units of mass/volume
               delta_pu   =  i*mom_per_cell(i+2,j+2,k+2)
               delta_pv   =  j*mom_per_cell(i+2,j+2,k+2)
               delta_pw   =  k*mom_per_cell(i+2,j+2,k+2)
               delta_therm = therm_per_cell
c              delta_p's have units of mass*velocity/volume
c              delta_therm has units of mass*velocity^2/volume
#ifdef FORTRAN_DEBUG
               write(6,*) 'i,j,k', i, j, k
               write(6,*) 'p_before:',pu(ic+i,jc+j,kc+k)*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            pv(ic+i,jc+j,kc+k)*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            pw(ic+i,jc+j,kc+k)*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND
               write(6,*) 'delta_p:',delta_pu*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            delta_pv*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            delta_pw*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND
#endif
c
c     Add mass, momentum
c     (add thermal energy here)
c           
               pu_start = pu(ic+i,jc+j,kc+k)
               pv_start = pv(ic+i,jc+j,kc+k)
               pw_start = pw(ic+i,jc+j,kc+k)
               dratio = d(ic+i, jc+j, kc+k)/
     &                       (d(ic+i, jc+j, kc+k) + delta_mass)
               d(ic+i,jc+j,kc+k) = d(ic+i,jc+j,kc+k)
     &                       + delta_mass
               pu(ic+i,jc+j,kc+k) = 
     &                       pu_start + delta_pu
               pv(ic+i,jc+j,kc+k) = 
     &                       pv_start + delta_pv
               pw(ic+i,jc+j,kc+k) = 
     &                       pw_start+ delta_pw
c                        
c        If absolute value of momentum after deposition is less than it was
c        before deposition, add the canceled momentum as thermal energy
c
               pu_deposit = pu(ic+i,jc+j,kc+k)**2._RKIND - 
     &                      pu_start**2._RKIND
               if (pu_deposit .lt. 0._RKIND)
     &             delta_therm = delta_therm - 0.5_RKIND * 
     &                     pu_deposit/d(ic+i, jc+j, kc+k)
               pv_deposit = pv(ic+i,jc+j,kc+k)**2._RKIND - 
     &                      pv_start**2._RKIND
               if (pv_deposit .lt. 0._RKIND)
     &            delta_therm = delta_therm - 0.5_RKIND * 
     &                     pv_deposit/d(ic+i, jc+j, kc+k)
               pw_deposit = pw(ic+i,jc+j,kc+k)**2._RKIND - 
     &                      pw_start**2._RKIND
               if (pw_deposit .lt. 0._RKIND)
     &            delta_therm = delta_therm - 0.5_RKIND * 
     &                     pw_deposit/d(ic+i, jc+j, kc+k)
#ifdef FORTRAN_DEBUG 
               write(6,*) 'p_after,delta_therm:',
     &           pu(ic+i, jc+j, kc+k)*dx**3._RKIND*
     &           (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &           pv(ic+i, jc+j, kc+k)*dx**3._RKIND*
     &           (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &           pw(ic+i, jc+j, kc+k)*dx**3._RKIND*
     &           (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &           delta_therm*dx**3._RKIND*
     &           (dunits*x1**3._RKIND)*vunits**2._RKIND
#endif
c
c     Add thermal energy
c     
#ifdef FORTRAN_DEBUG
               write(6,*) 'te_before',te(ic+i,jc+j,kc+k)*
     &            (dunits*(dx*x1)**3._RKIND*vunits**2._RKIND)
#endif
               te(ic+i, jc+j, kc+k) = 
     &                   te(ic+i, jc+j, kc+k)*dratio + 
     &                   delta_therm / d(ic+i, jc+j, kc+k)
#ifdef FORTRAN_DEBUG
               write(6,*) 'te_after',te(ic+i,jc+j,kc+k)*
     &            (dunits*(dx*x1)**3._RKIND*vunits**2._RKIND)                        
               write(6,*) 'ge_before',ge(ic+i,jc+j,kc+k)*
     &            (dunits*(dx*x1)**3._RKIND*vunits**2._RKIND)
#endif
               if (idual .eq. 1_IKIND)
     &                   ge(ic+i, jc+j, kc+k) = 
     &                   ge(ic+i, jc+j, kc+k)*dratio +
     &                   delta_therm / d(ic+i, jc+j, kc+k)
#ifdef FORTRAN_DEBUG
               write(6,*) 'ge_after',ge(ic+i,jc+j,kc+k)*
     &            (dunits*(dx*x1)**3._RKIND*vunits**2._RKIND)
#endif
c     
c     Metal feedback (note that in this function gas metal is
c     a fraction (rho_metal/rho_gas) rather than a density.
c     The conversion has been done in the handling routine)
c     
               if (imetal .eq. 1_IKIND) then
c     
c     "Cen method".  This takes into account gas recycling.
c     (metal is not a fraction here)
c     
                  metal(ic+i, jc+j, kc+k) = 
     &                   metal(ic+i, jc+j, kc+k) + 
     &                   (delta_mass/m_eject)
     &                    * (yield * (1._RKIND-metalf) + 
     &                    m_eject * metalf)
                     
               endif
c     
c     End loop over cells in particle-frame
c     
            enddo
         enddo
      enddo
      return
      end

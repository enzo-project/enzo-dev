#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c#define FORTRAN_DEBUG
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK6 \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback6(nx, ny, nz,
     &               d, mu, dm, te, ge, u, v, w,
     &               metal, zfield1, zfield2,
     &               idual, imetal, imulti_metals, imethod, 
     &               dt, r, dx, t, z, procnum,
     &               dunits, x1, vunits, t1, thermal, m_eject, yield,
     &               npart, xstart, ystart, zstart, ibuff,
     &               xp, yp, zp, up, vp, wp,
     &               mp, tdp, tcp, metalf, type,
     &               particleID, justburn,
     &               exptime_in, mom_mult,
     &               mom_canc, feedback_log,
     &               use_tabfbk, minit, ergSNII, ergSNIa,
     &               itracksrc, metalSNII, metalSNIa,
     &               ntabZ, ntabAge,
     &               tabZ, tabAge, tabMass, tabMetal, tabEvents)

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
c    particleID - particle ID number given by ParticleNumber
c    exptime  - delay time for discrete explosions (if set to -1,
c               continuous energy injection used)
c
c    minit     - initial particle mass
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
c    justburn     - time-weighted mass of star formation (code units)
c
c
c-----------------------------------------------------------------------
      use tabular_feedback
      implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c#define FORTRAN_DEBUG
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, 
     &      imulti_metals, imethod, procnum
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
      INTG_PREC type(npart), particleID(npart)
      R_PREC 	 exptime_in,exptime
      R_PREC    thermal
      INTG_PREC mom_canc
      INTG_PREC feedback_log
      R_PREC    mom_mult
      R_PREC    minit(npart)
      R_PREC    ergSNII, ergSNIa, age
      INTG_PREC itracksrc, ntabZ, ntabAge, use_tabfbk
      R_PREC    metalSNII(nx,ny,nz), metalSNIa(nx,ny,nz)
      R_PREC    tabZ(ntabZ), tabAge(ntabAge)
      R_PREC    tabMass(4, ntabAge, ntabZ), tabMetal(4, ntabAge, ntabZ)
      R_PREC    tabEvents(2, ntabAge, ntabZ)
c
c  Locals
c
      INTG_PREC n, ic, jc, kc, ip, jp, kp,
     &     i, j, k, iface, jface, kface
      R_PREC mform, energy,
     &     m_eject, yield, minitial, xv1, xv2, distrad,
     &     distcells, mass_per_cell,
     &     mass_before, mass_dummy,
     &     kin_energy_before, mom_dummy,
     &     delta_ke, ke_injected, ke_after, kin_energy_inj,
     &     nsn, mass_end, kin_energy_end,
     &     mom_start, mom_end, kin_energy_dummy
      R_PREC xfcshift, yfcshift, zfcshift, xfc, yfc, zfc,
     &     xface, yface, zface, xpos, ypos, zpos, face_shift
      R_PREC ke_before(3,3,3), mom_per_cell(3,3,3)
      R_PREC Z_avg, n_avg, Z_floor,
     &	  realmass, mu_cell, fbuff, add_therm
      R_PREC te_sum, ke_sum, ge_sum, mom_sum, mass_sum, mom_term,
     &     mass_term, therm_sum, lost_p
      R_PREC mzeject, mzeject_ii, mzeject_ia
      R_PREC u1(3,3,3), v1(3,3,3), w1(3,3,3), d1(3,3,3),
     &     te1(3,3,3), ge1(3,3,3), metal1(3,3,3)

      character*32 feedbacklog_filename
      character*8   post

      if (feedback_log .eq. 1_IKIND) then
         write(post,'(i8)') procnum
         do i=1,8
            if(post(i:i).eq.' ') post(i:i)='0'
         enddo
         feedbacklog_filename = 'feedbacklog_' // post // '.txt'
         open(unit=4,file=feedbacklog_filename,status='unknown',
     &        position='append')
      endif
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
c      write(6,*) 'star_feedback6: start',npart
      do n=1_IKIND, npart
         if (tcp(n) .gt. 0._RKIND .and. mp(n) .gt. 0._RKIND .and. 
     &        type(n) .eq. 2_IKIND) then
c
c     Compute index of the cell that the star particle
c     resides in.
c 
            ip = int((xp(n) - xstart)/dx) + 1_IKIND
            jp = int((yp(n) - ystart)/dx) + 1_IKIND
            kp = int((zp(n) - zstart)/dx) + 1_IKIND
#ifdef FORTRAN_DEBUG
            write(6,*) 'n, mp(n):',n,mp(n)*dunits*(dx*x1)**3._RKIND
     &        /SolarMass
#endif
            if (xp(n) .lt. xstart .or. xp(n) .gt. xstart+dx*nx .or.
     &          yp(n) .lt. ystart .or. yp(n) .gt. ystart+dx*ny .or.
     &          zp(n) .lt. zstart .or. zp(n) .gt. zstart+dx*nz) then
               write(6,*) 'warning: star particle out of grid',
     &              xp(n),yp(n),zp(n), xstart, ystart, zstart
               goto 100
            endif
c
c     Assuming 3x3x3 cube, calculate cell distribution
c
            distrad = 3_IKIND
            distcells = distrad**3_IKIND
c
            age = (t-tcp(n)) * t1/yr_s
            if (use_tabfbk .eq. 1_IKIND) then
#ifdef FORTRAN_DEBUG
               write(6,*) 'creation time, age, dt:', tcp*t1/yr_s,
     &                        age, dt*t1/yr_s
               write(6,*) 'metalf:',metalf(n)
#endif
c
c         Integrate yield tables to get mass, energy, and metal yields
c
               call sne_mass(m_eject, minit(n), metalf(n), age,
     &                        dt, t1,
     &                        tabMass, tabZ, tabAge, ntabZ, ntabAge)
               call sne_energy(energy, minit(n), metalf(n), age,
     &                        ergSNIa, ergSNII,
     &                        dx, dt, x1, t1, dunits, vunits,
     &                        tabEvents, tabZ, tabAge, ntabZ, ntabAge)
c
c           Metal feedback. Similar to total ejected mass, 
c           table has metal mass yield per initial population mass.
c
               if (imetal .eq. 1) then
                  call sne_II_metal(mzeject_ii,
     &                 minit(n), metalf(n), age,
     &                 dt, t1,
     &                 tabMetal, tabZ, tabAge, ntabZ, ntabAge)
                  call sne_Ia_metal(mzeject_ia,
     &                 minit(n), metalf(n), age,
     &                 dt, t1,
     &                 tabMetal, tabZ, tabAge, ntabZ, ntabAge)
                  mzeject = mzeject_ii + mzeject_ia
                  mzeject = mzeject / distcells
                  mzeject_ii = mzeject_ii / distcells
                  mzeject_ia = mzeject_ia / distcells
               endif
c
c           Normalize energy by 10^51 ergs for momentum injection routine
               nsn = energy/(10._RKIND**51._RKIND / (dunits*
     &                  x1**3._RKIND*vunits**2._RKIND*dx**3._RKIND))
#ifdef FORTRAN_DEBUG
               write(6,*) 'nsn,m_eject,energy:',nsn,m_eject*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              energy*(dunits*x1**3._RKIND*vunits**2._RKIND*
     &               dx**3._RKIND)
#endif
c
c          Skip feedback if fewer than 0.01 SNe in this timestep
               if (nsn .lt. 0.01_RKIND) goto 10
c          Distribute mass among feedback zone
               mass_per_cell = m_eject / distcells
c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
               mp(n) = mp(n) - m_eject
               minitial = minit(n)
c
            else
c
c        The star_maker2 algorithm creates a star particle instantaneously.
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
                  exptime = exptime_in * Myr_s / t1  
                  if (tdp(n) .eq. 0._RKIND
     &                .or. abs(t-tcp(n)) .lt. exptime) goto 10
                  tdp(n) = 0._RKIND
                  minitial = mp(n)
                  mform = mp(n)
               else
c
c     Determine how much of a given star particle would have been 
c     turned into stars during this timestep.  Then calculate the mass
c     which should have formed during this timestep dt using the integral
c     form of the Cen & Ostriker formula.
c
                  xv1 = (t - tcp(n))/tdp(n)
#ifdef FORTRAN_DEBUG
                  write(6,*) t,tcp(n),tdp(n),xv1
#endif
                  if (xv1 .lt. 0._RKIND) goto 10
                  if (xv1 .gt. 12._RKIND) goto 10 ! t-tcp >> tdp so ignore
                  xv2 = (t + dt - tcp(n))/tdp(n)
c
c     First calculate the initial mass of the star particle 
c     in question.
                  minitial = mp(n) / 
     &               (1._RKIND - 
     &               m_eject*(1._RKIND - (1._RKIND + xv1)*exp(-xv1)))
c     
c     Then, calculate the amount of mass that would have formed in
c     this timestep.
c     
                  mform = minitial * ((1._RKIND + xv1)*exp(-xv1) - 
     &               (1._RKIND + xv2)*exp(-xv2))
                  mform = max(min(mform, mp(n)), 0._RKIND)
               endif

c          skip if very little mass is formed.
c
               if (mform/distcells/d(ip,jp,kp) .lt. 1e-10_RKIND) 
     &             goto 10
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
               mass_per_cell = mform * m_eject / distcells
c
c           Calculate the number of SNe in this timestep from the stars formed
c           in this timestep as 1 SN per 100 Msun formed. If a small fraction, skip.
c           
               realmass = mform * dunits * x1**3._RKIND*dx**3_IKIND
               realmass = realmass /SolarMass
               nsn = realmass / 100._RKIND
               if (nsn .lt. 1e-2_RKIND) goto 10
               energy = nsn * 10._RKIND**51._RKIND / (dunits*
     &                  x1**3._RKIND*vunits**2._RKIND*dx**3._RKIND)
c
c           Calculate the metals to deposit later in add_feedback function
c
               mzeject = (mform / distcells) * 
     &                   (yield * (1._RKIND-metalf(n)) + 
     &                   m_eject * metalf(n))
            endif
c
#ifdef FORTRAN_DEBUG
            write(6,*) 'mass_per_cell:',mass_per_cell*
     &         dunits*x1**3._RKIND/SolarMass*dx**3_IKIND
            write(6,*) 'energy:',energy * (dunits*
     &                  x1**3._RKIND*vunits**2._RKIND*dx**3._RKIND)
#endif
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
     &		        nx, ny,nz
#endif
            endif
c
c           If using zeus, then velocities are face-centered so shift
c
            face_shift = 0._RKIND
            if (imethod .eq. 2) face_shift = 0.5_RKIND
c
c         Compute index of center of feedback zone for momentum,
c         accounting for possible face-centering
c
            xface = (xfc - xstart)/dx - face_shift
            yface = (yfc - ystart)/dx - face_shift
            zface = (zfc - zstart)/dx - face_shift
c
            iface = int(xface) + 1_IKIND
            jface = int(yface) + 1_IKIND
            kface = int(zface) + 1_IKIND
c            write(6,*) 'iface, jface, kface:', iface, jface, kface
c
c         Compute index of the center of feedback zone for mass
c         Note that this is different than star_feedback3mom because
c         this is center-indexed and star_feedback3mom is corner-indexed
c 
            xpos = (xfc - xstart)/dx
            ypos = (yfc - ystart)/dx
            zpos = (zfc - zstart)/dx
c
            ic = int(xpos) + 1_IKIND
            jc = int(ypos) + 1_IKIND
            kc = int(zpos) + 1_IKIND
c            write(6,*) 'ic, jc, kc:', ic, jc, kc
c
            add_therm = thermal
#ifdef FORTRAN_DEBUG
            write(6,*) 'ic,jc,kc:',ic, jc, kc
#endif
c
            Z_avg = 0._RKIND
            n_avg = 0._RKIND
c     
c     Compute average density and metallicity around cell containing
c     particle
c     
c            write(6,*) 'averaging over patch'
            do k = -1_IKIND,1_IKIND
               do j = -1_IKIND,1_IKIND
                  do i = -1_IKIND,1_IKIND
                     mu_cell = mu(ic+i,jc+j,kc+k)
                     Z_avg = Z_avg + 
     &                 metal(ic+i,jc+j,kc+k)/0.02_RKIND
                     n_avg = n_avg + 
     &                 d(ic+i,jc+j,kc+k)*dunits*
     &                 mu_cell/mass_h
                  enddo
               enddo
            enddo
c               
c        The floor on Z_avg is required by the Kimm & Cen (2014) scheme (eq A5)
            Z_floor = max(Z_avg/27._RKIND, 0.01)
            n_avg = n_avg/27._RKIND
#ifdef FORTRAN_DEBUG
            write(6,*) 'n_avg, Z_avg:', n_avg, Z_avg/27._RKIND
#endif
c
c     Inject energy, mass & metals
c
c     Compute the momentum that will be added to each cell.
c
            call get_mom_per_cell(d, mass_per_cell,
     &                  distcells, mom_per_cell,
     &                  nsn, d(ic,jc,kc), dx, x1, n_avg, Z_floor,
     &                  dunits, vunits, t1,
     &                  nx, ny, nz, ic, jc, kc,
     &                  mom_mult)
c           mom_per_cell has units of mass*velocity/volume
c           If not injecting any momentum, skip
            mom_sum = 0._RKIND
            do k = 1_IKIND,3_IKIND
               do j = 1_IKIND,3_IKIND
                  do i = 1_IKIND,3_IKIND
                     mom_sum = mom_sum + mom_per_cell(i, j, k)
                  enddo
               enddo
            enddo
            if (mom_sum .eq. 0._RKIND) then
#ifdef FORTRAN_DEBUG
               write(6,*) 'no momentum injection for particle', n
#endif
               goto 10
            endif
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
            write(6,*) 'ke_before:',ke_before
     &                  *dunits*(dx*x1)**3._RKIND*vunits**2._RKIND
#endif
c           ke_before has units of mass*velocity^2/volume
c
c     First convert velocities to momenta and transform
c     into frame comoving with particle
c
            call momentum_convert(u, v, w, d, metal, up(n), vp(n),
     &                    wp(n), nx, ny, nz, ic, jc, kc, 
     &                    iface, jface, kface, imethod, imetal, 
     &                    +1_IKIND, dx, dunits, vunits, x1)
c
c           u,v,w have units of mass*velocity/volume
c           Now u,v,w are the particle-relative *momentum* directions
c
#ifdef FORTRAN_DEBUG
c           Sum mass, energy, and momentum before
            call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &            ic, jc, kc, iface, jface, kface,
     &            mass_before, kin_energy_before, mom_start)

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
c     Zero local dummy field and kinetic energy field
c
            do k = 1_IKIND, 3_IKIND
               do j = 1_IKIND, 3_IKIND
                  do i = 1_IKIND, 3_IKIND
                     u1(i,j,k) = 0._RKIND
                     v1(i,j,k) = 0._RKIND
                     w1(i,j,k) = 0._RKIND
                     d1(i,j,k) = d(i,j,k)
                     ge1(i,j,k) = 0._RKIND
                     te1(i,j,k) = 0._RKIND
                     metal1(i,j,k) = 0._RKIND
                  enddo
               enddo
            enddo
c           Now add mass and momentum terms to
c              local dummy fields
c
            call add_feedback_mech(u1, v1, w1, d1, ge1, te1, metal1, 
     &                        3_IKIND, 3_IKIND, 3_IKIND,
     &                        2_IKIND, 2_IKIND, 2_IKIND,
     &                        2_IKIND, 2_IKIND, 2_IKIND,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf(n), itracksrc,
     &                        mzeject, mzeject_ia, mzeject_ii,
     &                        metalSNII, metalSNIa,
     &                        mass_per_cell, mom_per_cell, therm_sum,
     &                        0._RKIND, dx, x1, dunits, vunits,
     &                        mom_canc, lost_p)
c
c           Calculate mass and energy in dummy fields
c
            call sum_energy_mom(u1, v1, w1, d1,
     &          3_IKIND, 3_IKIND, 3_IKIND,
     &          2_IKIND, 2_IKIND, 2_IKIND, 2_IKIND, 2_IKIND, 2_IKIND,
     &            mass_dummy, kin_energy_dummy, mom_dummy)
#ifdef FORTRAN_DEBUG
            write(6,*) "m dummy:", mass_dummy*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
            write(6,*) "e dummy:", kin_energy_dummy*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
            write(6,*) 'mom_dummy:',mom_dummy*dunits*vunits
     &                  *(dx*x1)**3._RKIND/SolarMass/1e5_RKIND
#endif
c
c           If kin_energy_dummy is less than tabular SNe energy, add
c           the remaining energy as thermal into the host cell.
            if (kin_energy_dummy .lt. energy) then
               add_therm = add_therm + (energy - kin_energy_dummy)
#ifdef FORTRAN_DEBUG
               write(6,*) "leftover energy injected as thermal:",
     &                  (energy-kin_energy_dummy)*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
#endif
            endif
c
c           Now deposit the momentum onto the grid.
c
            therm_sum = 0._RKIND
            lost_p = 0._RKIND
            add_therm = add_therm/distcells
            call add_feedback_mech(u, v, w, d, ge, te, metal,
     &                         nx, ny, nz, ic, jc, kc,
     &                        iface, jface, kface,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf(n), itracksrc,
     &                        mzeject, mzeject_ia, mzeject_ii,
     &                        metalSNII, metalSNIa,
     &                        mass_per_cell, mom_per_cell, therm_sum,
     &                        add_therm, dx, x1, dunits, vunits,
     &                        mom_canc, lost_p)
            if (feedback_log .eq. 1_IKIND) then
               mom_sum = 0._RKIND
               do k = 1_IKIND, 3_IKIND
                  do j = 1_IKIND, 3_IKIND
                     do i = 1_IKIND, 3_IKIND
                        mom_sum = mom_sum + mom_per_cell(i,j,k)
                     enddo
                  enddo
               enddo
               write(4,*) t*t1/yr_s/1e6_RKIND, particleID(n),
     &           age/1e6_RKIND,
     &           minitial*(dunits*(x1*dx)**3._RKIND)/SolarMass, nsn,
     &           mass_per_cell*distcells*(dunits*(x1*dx)**3._RKIND)/
     &           SolarMass,
     &           mzeject*distcells*(dunits*(x1*dx)**3._RKIND)/
     &           SolarMass,
     &           mom_sum*
     &           (dunits*(dx*x1)**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &           therm_sum*(dunits*(dx*x1)**3._RKIND)*vunits**2._RKIND
            endif
c
c           Calculate the mass and energy in the affected region again.
c
#ifdef FORTRAN_DEBUG
            call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &                          ic, jc, kc, iface, jface, kface,
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
     &                    iface, jface, kface, imethod, imetal,
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
      if (feedback_log .eq. 1) then
         close(4)
      endif
      return
      end
c
c ==========================================================
c
c     Calculate the amount of momentum to add to each neighboring cell
c     Formulae taken from Kimm & Cen (2014), Appendix A
c
      subroutine get_mom_per_cell(d, mass_per_cell, 
     &               distcells, mom_per_cell,
     &               nsn, d_fbck_cell, dx, x1,
     &               n_avg, Z_floor, dunits, vunits, t1,
     &               nx, ny, nz, ic, jc, kc,
     &               mom_mult)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      R_PREC    d(nx,ny,nz), dx
      R_PREC    t1, Z_floor, n_avg
      R_PREC    mom_per_cell(3_IKIND, 3_IKIND, 3_IKIND)
      R_PREC    mass_per_cell, distcells
      R_PREC    nsn, x1, dunits, d_fbck_cell, vunits
      R_PREC    mom_mult
c
c     Locals
c
      INTG_PREC i, j, k, dist
      R_PREC    dMswept, chi, chi_th, fe, pSN, ergs_51
      R_PREC    den_cell, Mej, dMej, Nnbors, betaSN
c
        ergs_51 = 1.e51_RKIND
        Nnbors = distcells - 1._RKIND
        betaSN = 1._RKIND/distcells
        Mej = mass_per_cell*distcells*dunits*(x1*dx)**3._RKIND
        chi_th = 69.58_RKIND*(nsn)**(-2._RKIND/17._RKIND) * 
     &               (n_avg)**(-4._RKIND/17._RKIND) * 
     &               (Z_floor)**(-0.28_RKIND)
#ifdef FORTRAN_DEBUG
        write(6,*) 'chi_th:',chi_th
#endif
        do k = -1_IKIND, +1_IKIND
            do j = -1_IKIND, +1_IKIND
               do i = -1_IKIND, +1_IKIND
c
c       Compute the mass fraction of swept-up to ejected mass.
c
                  den_cell = d(ic+i,jc+j,kc+k)*dunits
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
     &                        (Z_floor)**(-0.14_RKIND)
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
                     dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
                     if (dist .eq. 0_IKIND) then
                        pSN = 0._RKIND
c                     else if (dist .eq. 2_IKIND) then
c                        pSN = pSN/(2._RKIND**0.5_RKIND)
c                     else if (dist .eq. 3_IKIND) then
c                        pSN = pSN/(3._RKIND**0.5_RKIND)
                     endif
                  mom_per_cell(i+2,j+2,k+2) = pSN*mom_mult
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
     &                    nx, ny, nz, ic, jc, kc, iface, jface, kface,
     &                    imethod, imetal, idir, dx, dunits, vunits, x1)
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
                     u(iface+i, jc+j, kc+k) = 
     &                    (u(iface+i, jc+j, kc+k) - up)
     &                    * 0.5_RKIND * (d(iface+i, jc+j, kc+k) +
     &                    d(iface+i+1_IKIND, jc+j, kc+k))
                     v(ic+i, jface+j, kc+k) = 
     &                    (v(ic+i, jface+j, kc+k) - vp)
     &                    * 0.5_RKIND * (d(ic+i, jface+j, kc+k) + 
     &                    d(ic+i, jface+j+1_IKIND, kc+k))
                     w(ic+i, jc+j, kface+k) = 
     &                    (w(ic+i, jc+j, kface+k) - wp)
     &                    * 0.5_RKIND * (d(ic+i, jc+j, kface+k) + 
     &                    d(ic+i, jc+j, kface+k+1_IKIND))
                  else
                     u(ic+i ,jc+j ,kc+k) = (u(ic+i ,jc+j ,kc+k)-up) *
     &                       d(ic+i, jc+j, kc+k)
                     v(ic+i ,jc+j ,kc+k) = (v(ic+i ,jc+j ,kc+k)-vp) *
     &                       d(ic+i, jc+j, kc+k)
                     w(ic+i ,jc+j ,kc+k) = (w(ic+i ,jc+j ,kc+k)-wp) *
     &                       d(ic+i, jc+j, kc+k)
                  endif
c
c              if idir = -1: convert mom -> vel
c               
               else
                  if (imethod .eq. 2_IKIND) then
                     u(iface+i, jc+j, kc+k) = u(iface+i, jc+j, kc+k)
     &                    /(0.5_RKIND * (d(iface+i, jc+j, kc+k) +
     &                             d(iface+i+1_IKIND, jc+j, kc+k))) + up
                     v(ic+i, jface+j, kc+k) = v(ic+i, jface+j, kc+k) 
     &                    /(0.5_RKIND * (d(ic+i, jface+j, kc+k) + 
     &                             d(ic+i, jface+j+1_IKIND, kc+k))) + vp
                     w(ic+i, jc+j, kface+k) = w(ic+i, jc+j, kface+k) 
     &                    /(0.5_RKIND * (d(ic+i, jc+j, kface+k) + 
     &                             d(ic+i, jc+j, kface+k+1_IKIND))) + wp
                  else
                     u(ic+i ,jc+j ,kc+k) = u(ic+i ,jc+j ,kc+k) /
     &                    (d(ic+i, jc+j, kc+k)) + up
                     v(ic+i ,jc+j ,kc+k) = v(ic+i ,jc+j ,kc+k) /
     &                    (d(ic+i, jc+j, kc+k)) + vp
                     w(ic+i ,jc+j ,kc+k) = w(ic+i ,jc+j ,kc+k) /
     &                    (d(ic+i, jc+j, kc+k)) + wp
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
c     Sum mass, momentum, and energy
c       Note that pu, pv, pw are momenta
c
      subroutine sum_energy_mom(pu, pv, pw, d,
     &                           nx, ny, nz,
     &                           ic, jc, kc, iface, jface, kface,
     &                           mass_sum, kin_energy_sum, mom_sum)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, iface, jface, kface
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      R_PREC    d(nx,ny,nz)
      R_PREC    mass_sum, kin_energy_sum
      R_PREC    mom_sum
c
c     Locals
c
      INTG_PREC i, j, k
      R_PREC    mass_term, mom_term, kin_energy
c
      mass_sum = 0._RKIND
      kin_energy_sum = 0._RKIND
      mom_sum = 0._RKIND
c
c     Sum mass and energy      
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
               mass_term = d(ic   +i ,jc   +j ,kc   +k)
               mom_term = pu(iface+i ,jc   +j, kc   +k)**2_IKIND + 
     &                    pv(ic   +i ,jface+j, kc   +k)**2_IKIND + 
     &                    pw(ic   +i ,jc   +j, kface+k)**2_IKIND
c              mass_term has units of mass/volume
c              mom_term has units of mass^2*velocity^2/volume^2
c
c              Compute total mass and energy
c              
               kin_energy = mom_term / (2._RKIND * mass_term)
c              kin_energy has units of mass^2*velocity^2/volume^2 / (mass/volume)
c                          = mass*velocity^2/volume
               mass_sum = mass_sum + mass_term
               kin_energy_sum = kin_energy_sum + kin_energy
               mom_sum = mom_sum + mom_term**0.5_RKIND
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
     &                        iface, jface, kface,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf, itracksrc,
     &                        mzeject, mzeject_ia, mzeject_ii,
     &                        metalSNII, metalSNIa,
     &                        mass_per_cell, mom_per_cell, therm_sum,
     &                        add_therm, dx, x1, dunits, vunits,
     &                        mom_canc, lost_p)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, idual, dist
      INTG_PREC iface, jface, kface, imethod, imetal, imulti_metals
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      R_PREC    mass_per_cell, mom_per_cell(3_IKIND,3_IKIND,3_IKIND)
      R_PREC    m_eject, yield, metalf, add_therm, therm_sum
      R_PREC    mzeject, mzeject_ia, mzeject_ii
      R_PREC    metalSNII(nx,ny,nz), metalSNIa(nx,ny,nz)
      R_PREC    dunits, x1, vunits, dx, mult
      R_PREC    p_start, p_end, p_inj, extra_p, lost_p
      INTG_PREC    mom_canc, itracksrc
c
c     Locals
c
      INTG_PREC i, j, k, i1, j1, k1
      R_PREC    delta_mass, delta_pu, delta_pv, delta_pw, delta_therm,
     &          dratio, delta_p, vel_kick
      R_PREC    pu_start, pu_deposit, pv_start, pv_deposit, pw_start,
     &          pw_deposit
c
      therm_sum = 0._RKIND
      lost_p = 0._RKIND
c
c     Loop cells in particle-frame
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
c
               dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
               if (dist .eq. 0_IKIND) then
                  mult = 0._RKIND
               else if (dist .eq. 1_IKIND) then
                  mult = 1._RKIND
               else if (dist .eq. 2_IKIND) then
                  mult = 1._RKIND/(2._RKIND**0.5_RKIND)
               else if (dist .eq. 3_IKIND) then
                  mult = 1._RKIND/(3._RKIND**0.5_RKIND)
               endif
               delta_mass =   mass_per_cell
               delta_therm = add_therm
c              delta_mass has units of mass/volume
               delta_pu   =  i*mult*mom_per_cell(i+2,j+2,k+2)
               delta_pv   =  j*mult*mom_per_cell(i+2,j+2,k+2)
               delta_pw   =  k*mult*mom_per_cell(i+2,j+2,k+2)
               delta_p = (delta_pu**2._RKIND + delta_pv**2._RKIND + 
     &            delta_pw**2._RKIND)**(0.5_RKIND)
c              If delta_p leads to a velocity kick of more than
c              1000 km/s, limit to 1000
               vel_kick = delta_p/(d(ic+i,jc+j,kc+k) + delta_mass)*
     &                    vunits/1e5_RKIND
               if (vel_kick .gt. 1000._RKIND) then
                  delta_p = 1000._RKIND/(vunits/1e5_RKIND)*
     &                   (d(ic+i,jc+j,kc+k) + delta_mass)
               endif
c              
c              delta_p's have units of mass*velocity/volume
c              delta_therm has units of mass*velocity^2/volume
#ifdef FORTRAN_DEBUG
               write(6,*) 'i,j,k, dist', i, j, k, dist
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
               pu_start = pu(iface+i,jc+j,kc+k)
               pv_start = pv(ic+i,jface+j,kc+k)
               pw_start = pw(ic+i,jc+j,kface+k)
               p_start = (pu_start**2._RKIND + pv_start**2._RKIND +
     &           pw_start**2.)**(0.5_RKIND)
               dratio = d(ic+i, jc+j, kc+k)/
     &                       (d(ic+i, jc+j, kc+k) + delta_mass)
               d(ic+i,jc+j,kc+k) = d(ic+i,jc+j,kc+k)
     &                       + delta_mass
               pu(iface+i,jc+j,kc+k) = 
     &                       pu_start + delta_pu
               pv(ic+i,jface+j,kc+k) = 
     &                       pv_start + delta_pv
               pw(ic+i,jc+j,kface+k) = 
     &                       pw_start + delta_pw
               p_end = (pu(iface+i,jc+j,kc+k)**2._RKIND + 
     &                  pv(ic+i,jface+j,kc+k)**2._RKIND +
     &                  pw(ic+i,jc+j,kface+k)**2._RKIND)**(0.5_RKIND)
               if (mom_canc .eq. 1_IKIND) then
c
c     If the difference between end and start momentum is less than the
c     total amount that should've been deposited (due to momentum cancellation),
c     add the extra as thermal energy
c
                  p_inj = ABS(p_end - p_start)
#ifdef FORTRAN_DEBUG
                  write(6,*) 'p_start,p_end,p_inj',
     &            p_start*(x1*dx)**3._RKIND*dunits*vunits
     &            /SolarMass/1e5_RKIND,
     &            p_end*(x1*dx)**3._RKIND*dunits*vunits
     &            /SolarMass/1e5_RKIND,
     &            p_inj*(x1*dx)**3._RKIND*dunits*vunits
     &            /SolarMass/1e5_RKIND
#endif
                  if (p_inj .lt. delta_p) then
                     extra_p = delta_p - p_inj
                     lost_p = lost_p + extra_p**2._RKIND /
     &                 d(ic+i,jc+j,kc+k)
#ifdef FORTRAN_DEBUG
                     write(6,*) 'lost_p', extra_p*(x1*dx)**3._RKIND*
     &                 dunits*vunits/SolarMass/1e5_RKIND
                     write(6,*) 'lost_p energy', extra_p**2._RKIND /
     &                 d(ic+i,jc+j,kc+k)*
     &                 (dunits*(dx*x1)**3._RKIND*vunits**2._RKIND)
#endif
                     delta_therm = delta_therm + extra_p**2._RKIND /
     &                           d(ic+i,jc+j,kc+k)
                  endif
               endif
               therm_sum = therm_sum + delta_therm
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
c           Field is converted to metal fraction by handling routine.
c
                  metal(ic+i,jc+j,kc+k) = 
     &               (metal(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) + 
     &               mzeject) /(d(ic+i, jc+j, kc+k) + delta_mass)
                  if (itracksrc .eq. 1) then
                     metalSNII(ic+i,jc+j,kc+k) =
     &                 (metalSNII(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) + 
     &                 mzeject_ii) /(d(ic+i, jc+j, kc+k) + delta_mass)
                     metalSNIa(ic+i,jc+j,kc+k) =
     &                 (metalSNIa(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) + 
     &                 mzeject_ia) /(d(ic+i, jc+j, kc+k) + delta_mass)
                  endif
               endif
c     
c     End loop over cells in particle-frame
c     
            enddo
         enddo
      enddo
      return
      end

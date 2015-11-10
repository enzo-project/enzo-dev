#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c////////////////////////  SUBROUTINE POP3_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine pop3_maker(nx, ny, nz,
     &                      d, dm, h2d, h2diss, kph, temp, u, v, w, 
     &                      cooltime, dt, r, metal, dx, t, z, procnum, 
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, h2crit, metalcrit,
     &                      odthresh, starmass, level, np,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      type, ctype, justburn, iradtrans)
c
c  CREATES PRIMORDIAL STAR PARTICLES
c
c  written by: Chris Loken
c  date:       3 March 1997
c  modified1: 2 March 1999 by Brian O''Shea
c    2 inputs were added: odthresh and masseff, and the star creation
c    code was implemented
c  modified2: 18 May 1999 by Brian O''Shea
c    1 input added: smthresh, and star particle selection code added
c  modified3: 30 August 2005 by John Wise
c    2 inputs added: type and ctype to distinugish between particles now
c    that we have multiple star particle types, and added BWO's fix on
c    the runaway star particle phenomenon by averaging the velocities
c    over 5 cells
c  modified4: 30 August 2005 by John Wise
c    modified star_maker2 for single primordial star formation
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    h2d   - molecular hydrogen field
c    temp  - temperature field
c    u,v,w - velocity fields
c    cooltime - cooling time in code units
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    odthresh - overdensity threshold (some number * avg. density)
c    starmass - primordial star mass (solar masses)
c    h2crit   - critical value for primordial star formation
c    metalcrit- critical value for primordial star formation

c    level - current level of refinement
c    procnum - processor number (for output)
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
c    type     - particle types
c    ctype    - current type to set the particles to
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
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      INTG_PREC procnum, ctype, iradtrans
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), h2d(nx,ny,nz)
      R_PREC    h2diss(nx,ny,nz), kph(nx,ny,nz)
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, justburn
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      INTG_PREC type(nmax)
      R_PREC    odthresh, masseff, starmass, mintdyn, h2crit, lifetime
      R_PREC    metalcrit
c
c  Locals:
c
      INTG_PREC  i, j, k, ii, irad, ib, jb, kb, maxradius, radius2
      R_PREC   div, tdyn, dtot, wgt
      R_PREC   pi, G, sndspdC, m1
      R_PREC   isosndsp2, bmass, jeanmass, lum
      real*8 msolar, mh
      parameter (pi=pi_val, G=GravConst, 
     &           sndspdC=1.3095e8_RKIND,
     &           msolar=SolarMass, mh=mass_h)
c
      ii = np
c
c     calculate how many solar masses are in d1*x1^3
      m1 = dble(x1*dx)**3 * dble(d1) / msolar
c
c  for each zone, : a primordial "star" particle is created if answers
c  to all the following questions are affirmative:
c
c    is this the finest level of refinement ?
c    is the density greater than a critical density ?
c    is the molecular hydrogen fraction greater than critical value ?
c    is the metallicity less than a critical value ?
c    is the flow convergent ?
c    is the cooling time less than a dynamical time ? 
c
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
c
c              1) finest level of refinement?
c
               if (r(i,j,k) .ne. 0._RKIND) goto 10
c
c              2) density greater than threshold
c
               if (d(i,j,k) .lt. odthresh) goto 10
c
c              3) divergence negative
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
c               write(6,*) 'b:', d(i,j,k), div, temp(i,j,k), h2d(i,j,k), h2crit
               if (div .ge. 0._RKIND) goto 10
c
c              4) t_cool < t_free-fall (if T < 1.1e4 skip this check)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*d1
               tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1
c               if (tdyn .lt. cooltime(i,j,k) .and. 
c     &             temp(i,j,k) .gt. 1.1e4) 
c     &            write(40,*) '1',tdyn,cooltime(i,j,k),temp(i,j,k)

c               write(6,*) 'b:', i,j,k,d(i,j,k), div, temp(i,j,k), 
c     $              tdyn, cooltime(i,j,k), h2d(i,j,k), h2crit

               if (tdyn .lt. cooltime(i,j,k) .and. 
     &             temp(i,j,k) .gt. 1.0e3_RKIND) goto 10
c
c              JHW (30 AUG 2005) : REMOVED JEANS MASS CRITERION
c
c              5) M > M_Jeans (this definition involves only baryons under
c                 the assumption that the dark matter is stable, which
c                 implies that the dark matter velocity dispersion is >> 
c                 the sound speed.  This will be true for small perturbations
c                 within large halos).
c
               bmass = d(i,j,k) * m1
c               isosndsp2 = sndspdC * temp(i,j,k)
c               jeanmass = pi/(6._RKIND*sqrt(d(i,j,k)*dble(d1))) *
c     &                    dble(pi * isosndsp2 / G)**1.5_RKIND / msolar
c
c               if (bmass .lt. jeanmass) 
c     &            write(40,*) '2',bmass,jeanmass
c               if (bmass .lt. jeanmass) goto 10
c
c
c              6) If baryon mass is greater than threshold, then make it
c
c               if (bmass .lt. starmass) goto 10
c
c              7) Check if metallicity is less than critical value
c
               if (imetal .eq. 1 .and. metal(i,j,k) .gt. metalcrit)
     $              goto 10
c
c              8) Check if molecular hydrogren fraction is greater than
c              critical value and there is no ionizing radiation
c
               if (h2d(i,j,k) .lt. h2crit) goto 10
               if (kph(i,j,k) .gt. tiny) goto 10
c
c               write(6,7777) i, j, k
c               write(6,*) 'pop3:', i,j,k, d(i,j,k), div, tdyn, 
c     &              cooltime(i,j,k)
c               write(6,*) temp(i,j,k), h2d(i,j,k), bmass
c               write(6,*) h2diss(i,j,k), kph(i,j,k)
c
c              Create a star particle
c
               call pop3_properties(starmass, lum, lifetime)

               ii = ii + 1
               if (iradtrans .eq. 0) then
c                  type(ii) = -ctype
                  mp(ii)   = 0   !starmass / m1
                  tdp(ii) = 1.5_RKIND*lifetime*3.15e7_RKIND/t1
               else
c              because it's possible that we'll be using an IMF, adjust
c              these later
                  mp(ii)   = min(0.5_RKIND * bmass, starmass/m1)
                  tdp(ii) = lifetime*3.15e7_RKIND/t1
c                  mp(ii)   = min(0.5_RKIND * bmass, starmass/m1)
c                  tdp(ii) = lifetime*3.15e7_RKIND/t1
               endif

               type(ii) = -ctype
               tcp(ii) = t
               xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
c
c     Record amount of star formation
               justburn = justburn + mp(ii)
c
c              Star velocities averaged over multiple cells to
c              avoid "runaway star particle" phenomenon
c              imethod = 2 is zeus, otherwise PPM
c
               if (imethod .eq. 2) then
c                  up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
c                  vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
c                  wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))

              up(ii) = ( 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))*d(i,j,k) +
     &          0.5_RKIND*(u(i-1,j,k)+u(i,j,k))*d(i-1,j,k)        +
     &          0.5_RKIND*(u(i+1,j,k)+u(i+2,j,k))*d(i+1,j,k)      +
     &          0.5_RKIND*(u(i,j+1,k)+u(i+1,j+1,k))*d(i,j+1,k)    +        
     &          0.5_RKIND*(u(i,j-1,k)+u(i+1,j-1,k))*d(i,j-1,k)    + 
     &          0.5_RKIND*(u(i,j,k+1)+u(i+1,j,k+1))*d(i,j,k+1)    +        
     &          0.5_RKIND*(u(i,j,k-1)+u(i+1,j,k-1))*d(i,j,k-1) ) / 
     &          ( d(i,j,k)+d(i-1,j,k)+d(i+1,j,k) +
     &           d(i,j-1,k)+d(i,j+1,k)           +
     &           d(i,j,k-1)+d(i,j,k+1) ) 

              vp(ii) = ( 0.5_RKIND*(v(i,j,k)+v(i+1,j,k))*d(i,j,k) +
     &          0.5_RKIND*(v(i-1,j,k)+v(i,j,k))*d(i-1,j,k)        +
     &          0.5_RKIND*(v(i+1,j,k)+v(i+2,j,k))*d(i+1,j,k)      +
     &          0.5_RKIND*(v(i,j+1,k)+v(i+1,j+1,k))*d(i,j+1,k)    +
     &          0.5_RKIND*(v(i,j-1,k)+v(i+1,j-1,k))*d(i,j-1,k)    +
     &          0.5_RKIND*(v(i,j,k+1)+v(i+1,j,k+1))*d(i,j,k+1)    +
     &          0.5_RKIND*(v(i,j,k-1)+v(i+1,j,k-1))*d(i,j,k-1) ) /
     &          ( d(i,j,k)+d(i-1,j,k)+d(i+1,j,k) +
     &           d(i,j-1,k)+d(i,j+1,k)           +
     &           d(i,j,k-1)+d(i,j,k+1) )

              wp(ii) = ( 0.5_RKIND*(w(i,j,k)+w(i+1,j,k))*d(i,j,k) +
     &          0.5_RKIND*(w(i-1,j,k)+w(i,j,k))*d(i-1,j,k)        +
     &          0.5_RKIND*(w(i+1,j,k)+w(i+2,j,k))*d(i+1,j,k)      +
     &          0.5_RKIND*(w(i,j+1,k)+w(i+1,j+1,k))*d(i,j+1,k)    +
     &          0.5_RKIND*(w(i,j-1,k)+w(i+1,j-1,k))*d(i,j-1,k)    +
     &          0.5_RKIND*(w(i,j,k+1)+w(i+1,j,k+1))*d(i,j,k+1)    +
     &          0.5_RKIND*(w(i,j,k-1)+w(i+1,j,k-1))*d(i,j,k-1) ) /
     &          ( d(i,j,k)+d(i-1,j,k)+d(i+1,j,k) +
     &           d(i,j-1,k)+d(i,j+1,k)           +
     &           d(i,j,k-1)+d(i,j,k+1) )
 
               else
c                  up(ii) = u(i,j,k)
c                  vp(ii) = v(i,j,k)
c                  wp(ii) = w(i,j,k)

                  up(ii) = (u(i,j,k)*d(i,j,k) +
     &                      u(i-1,j,k)*d(i-1,j,k) +
     &                      u(i+1,j,k)*d(i+1,j,k) +
     &                      u(i,j-1,k)*d(i,j-1,k) +
     &                      u(i,j+1,k)*d(i,j+1,k) +
     &                      u(i,j,k-1)*d(i,j,k-1) +
     &                      u(i,j,k+1)*d(i,j,k+1) ) / 
     &                    ( d(i,j,k) + d(i-1,j,k) + d(i+1,j,k) +
     &                      d(i,j-1,k) + d(i,j+1,k) + d(i,j,k-1) +
     &                      d(i,j,k+1) )

                  vp(ii) = (v(i,j,k)*d(i,j,k) +
     &                      v(i-1,j,k)*d(i-1,j,k) +
     &                      v(i+1,j,k)*d(i+1,j,k) +
     &                      v(i,j-1,k)*d(i,j-1,k) +
     &                      v(i,j+1,k)*d(i,j+1,k) +
     &                      v(i,j,k-1)*d(i,j,k-1) +
     &                      v(i,j,k+1)*d(i,j,k+1) ) / 
     &                    ( d(i,j,k) + d(i-1,j,k) + d(i+1,j,k) +
     &                      d(i,j-1,k) + d(i,j+1,k) + d(i,j,k-1) +
     &                      d(i,j,k+1) )

                  wp(ii) = (w(i,j,k)*d(i,j,k) +
     &                      w(i-1,j,k)*d(i-1,j,k) +
     &                      w(i+1,j,k)*d(i+1,j,k) +
     &                      w(i,j-1,k)*d(i,j-1,k) +
     &                      w(i,j+1,k)*d(i,j+1,k) +
     &                      w(i,j,k-1)*d(i,j,k-1) +
     &                      w(i,j,k+1)*d(i,j,k+1) ) / 
     &                    ( d(i,j,k) + d(i-1,j,k) + d(i+1,j,k) +
     &                      d(i,j-1,k) + d(i,j+1,k) + d(i,j,k-1) +
     &                      d(i,j,k+1) )

               endif
c
c              Set the particle metal fraction
c
               if (imetal .eq. 1) then
                  metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
c                  metalf(ii) = metal(i,j,k)/d(i,j,k)
               else
                  metalf(ii) = 0._RKIND
               endif
c
c              Remove mass from grid

c               d(i,j,k) = 0.5_RKIND * d(i,j,k)

c               write(6,777) d(i,j,k), xp(ii), yp(ii), zp(ii)
c               write(6,7777) i, j, k
c               write(6,778) up(ii), vp(ii), wp(ii), mp(ii)
c               write(6,779) type(ii), tcp(ii), tdp(ii), metalf(ii)
 777           format('pop3_maker: ', e10.3, 3(f12.6))
 7777          format('            ', 3(i5))
 778           format('            ', 4(e15.3))
 779           format('            ', i5, 3(f15.6))
c
c               write(7+procnum,1000) bmass*starfraction,tdp(ii),tcp(ii),
c     &                       metalf(ii)
c               write(7+procnum,1000) level,bmass*starfraction,tcp(ii),
c     &                           tdp(ii)*t1,d(i,j,k)*d1,z,metalf(ii)

 1000          format(i5,1x,6(1pe10.3,1x))
c
c              Do not generate more star particles than available
c
               if (ii .eq. nmax) goto 20
c
10          continue
c
            enddo
         enddo
      enddo
 20   continue
c	
c      if (ii .ge. nmax) then
c         write(6,*) 'pop3_maker: reached max new particle count'
c         stop
c      endif
      np = ii

c      if (np .ne. 0) then
c         write(6,*) 'pop3_maker: number,time,level: ', np, t, level, 
c     $        irad
c      endif
c
      return
      end
c
c=======================================================================
c/////////////////////  SUBROUTINE POP3_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine pop3_feedback (nx, ny, nz,
     &     d, dm, te, ge, u, v, w, metal,
     &     idual, imetal, imethod, dt, r, dx, t, z,
     &     d1, x1, v1, t1, 
     &     npart, xstart, ystart, zstart, ibuff,
     &     xp, yp, zp, up, vp, wp,
     &     mp, life, tcp, metalf, type, ctype)
c
c  RELEASES PRIMORDIAL "STAR" PARTICLE ENERGY
c
c  written by: Chris Loken & Greg Bryan
c  date:       3 March 1997
c  modified1:  31 August 2005 by John Wise
c     Changed the feedback to mimic single primordial stars
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
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    npart    - particle array size specified by calling routine
c    type     - particle type
c    ctype    - particle that provides feedback
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
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, imethod, ctype
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), life(npart), tcp(npart), metalf(npart)
      INTG_PREC type(npart)
c
c  Locals
c
      INTG_PREC i, j, k, n
      R_PREC energy, lum, sn_nrg, m_eject, yield, dratio, m1
      R_PREC lower_pisn, upper_pisn, lsun, he_core, starmass, lifetime
      real*8 msolar
      parameter (lower_pisn = 140._RKIND, upper_pisn = 260._RKIND, 
     $           lsun = 2.e33_RKIND)
      parameter (msolar = SolarMass)
      LOGIC_PREC death, pisn
c
c-----------------------------------------------------------------------
c
c     Calculate mass conversion factor
c
      m1 = dble(x1*dx)**3 * dble(d1) / msolar
c
c     Loop over particles
c
c      write(6,*) 'star_feedback2: start'
      do n=1, npart
         if (type(n) .gt. ctype .and. tcp(n) .gt. 0 .and. 
     $       mp(n) .gt. 0 .and. life(n) .gt. 0) then
c
c     Based on the stellar mass, calculate its luminosity, SN energy,
c     and ejecta mass
            he_core = (13._RKIND/24._RKIND) * (mp(n)*m1 - 20._RKIND)
            sn_nrg = 5._RKIND + 1.304_RKIND * (he_core - 64._RKIND)
            starmass = mp(n)*m1
            call pop3_properties(starmass, lum, lifetime)
c
c     Determine whether the star burns for the entire timestep.  If not,
c     energy is (life - time_elapsed) * lum + sn_nrg.  If lifetime is
c     over, ignore this particle.
c
            if (t - tcp(n) .gt. life(n)) goto 10
            if (t + dt - tcp(n) .gt. life(n)) then
               energy = (life(n) - (t - tcp(n))) * t1 * lum * lsun + 
     $              sn_nrg * 1.d51
               death = .true.
            else
               energy = dt * t1 * lum * lsun
               death = .false.
            endif
            energy = energy / dble(m1*msolar) / dble(v1)**2
c
c           Determine whether the star results in a pair-instability SN
c
            if (starmass.gt.lower_pisn.and.starmass.lt.upper_pisn) then
               pisn = .true.
            else
               pisn = .false.
            endif
c
c           Compute index of nearest cell
c 
            i = int((xp(n) - xstart)/dx,IKIND) + 1
            j = int((yp(n) - ystart)/dx,IKIND) + 1
            k = int((zp(n) - zstart)/dx,IKIND) + 1
c
            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &          .or. k .lt. 1 .or. k .gt. nz) then
               write(6,*) 'warning: star particle out of grid',i,j,k
               goto 100
            endif
c
c     DURING THE STELLAR LIFETIME
c
            if (.not. death) then
               m_eject = 0._RKIND
               yield = 0._RKIND
            endif
c
c     STELLAR ENDPOINT: PISN
c
            if (death .and. pisn) then
               m_eject = mp(n)
               yield = m_eject
               mp(n) = 0._RKIND       ! totally disrupted
            endif
c
c     STELLAR ENDPOINT: BH FORMATION
c
            if (death .and. .not. pisn) then
               m_eject = 0._RKIND 
               yield = 0._RKIND
               type(n) = 2         ! BH formation
               life(n) = 0._RKIND      ! no dynamical time
            endif
c
c           Add energy to energy field
c
            dratio = d(i,j,k)/(d(i,j,k) + m_eject)
            te(i,j,k) = te(i,j,k)*dratio + energy / (d(i,j,k) + m_eject)
            if (idual .eq. 1) 
     &           ge(i,j,k) = ge(i,j,k)*dratio + 
     &                       energy / (d(i,j,k) + m_eject)
c            if (te(i,j,k) .le. 0.0) then
c               write(6,*) 'star',te(i,j,k),energy,sn_param,d(i,j,k)
c               stop
c            endif
c
c           Metal feedback (in here metal is a fraction)
c
            if (imetal .eq. 1) then
               metal(i,j,k) = (metal(i,j,k)*d(i,j,k) + yield) / 
     $              (d(i,j,k) + m_eject)
            endif
c
c           Mass and momentum feedback
c
            u(i,j,k) = u(i,j,k)*d(i,j,k) + m_eject * up(n)
            v(i,j,k) = v(i,j,k)*d(i,j,k) + m_eject * vp(n)
            w(i,j,k) = w(i,j,k)*d(i,j,k) + m_eject * wp(n)
            d(i,j,k) = d(i,j,k) + m_eject
            u(i,j,k) = u(i,j,k)/d(i,j,k)
            v(i,j,k) = v(i,j,k)/d(i,j,k)
            w(i,j,k) = w(i,j,k)/d(i,j,k)
c
c           If te is really total energy (and it is unless imethod=2),
c             then just set this value
c
            if (imethod .ne. 2 .and. idual .eq. 1) then
               te(i,j,k) = 0.5_RKIND*(u(i,j,k)**2 + v(i,j,k)**2 + 
     &                            w(i,j,k)**2) + ge(i,j,k)
            endif
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

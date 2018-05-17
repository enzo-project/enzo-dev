#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c

      subroutine star_maker_ssn(nx, ny, nz,
     &     d, dm, temp, u, v, w,
     &     dt, r, metal, dx, t, z, procnum,
     &     d1, x1, v1, t1,
     &     nmax, xstart, ystart, zstart, ibuff,
     &     imetal, imethod, tindsf,
     &     odthresh, useodthresh, masseff, smthresh, level, np,
     &     xp, yp, zp, up, vp, wp,
     &     mp, tdp, tcp, metalf,
     &     imetalSNIa, metalSNIa, metalfSNIa,
     &     imetalSNII, metalSNII, metalfSNII,
     &     mfcell)

c
c  CREATES STAR PARTICLES
c
c  Created 2014, adapted from star_maker2
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    temp  - temperature field
c    u,v,w - velocity fields
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    tindsf   - if true, remove dt/t_dyn term in starfraction
c    odthresh - overdensity threshold (some number * avg. density)
c    masseff - gas-to-mass conversion efficiency ( 0<=masseff<=1 )
c    smthresh - star mass threshold (only creates stars with mass >
c        smthresh unless (random number) < starmass/smthresh )
c    level - current level of refinement
c    procnum - processor number (for output)
c    imetalSNIa - SN Ia metallicity flag (0 - none, 1 - yes)
c    mfcell - the maximum baryonic mass of a cell we allow to become a star particle.
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
c    metalfSNIa - metallicity fraction of particle (from SN Ia)
c    nmax     - particle array size specified by calling routine
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
      INTG_PREC imetalSNIa, procnum, tindsf, useodthresh
      INTG_PREC imetalSNII
      R_PREC    metalSNII(nx,ny,nz)
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    metalSNIa(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax),
     &     metalfSNIa(nmax), metalfSNII(nmax)
      R_PREC    odthresh, masseff, smthresh, mfcell
c
      R_PREC   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      INTG_PREC  i, j, k, ii
      INTG_PREC, save :: ran1_init = 0
      R_PREC   div, tdyn, dtot
      R_PREC   pi, G, sndspdC
      R_PREC   isosndsp2, starmass, starfraction, bmass, jeanmass
      R_PREC   random, pstar
      real*8 msolar
      parameter (pi=pi_val, G=GravConst,
     &           sndspdC=1.3095e8_RKIND,
     &           msolar=SolarMass)
c
      ii = np

c     Initialize RNG the first time this routine is called
c     Use MPI rank to seed the RNG, ensuring different streams
c     on each processor.
      if(ran1_init .eq. 0) then
         call init_random_seed(procnum)
         ran1_init = 1
      endif

c
c  for each zone, : "star" particle is created if the density is
c  greater than a critical density
c
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
c
c              1) is density greater than threshold (if we are using it)?
c
               if (useodthresh .eq. 1) then
                  if (d(i,j,k) .lt. odthresh) then
                     goto 10
                  endif
               else
c              Make sure we never use more than half the mass in the cell to
c              create a particle.  One can think of this as an implicit
c              resolution-dependent density threshold.
                  if (mfcell*bmass .lt. smthresh) then
                     goto 10
                  endif
               endif

c
c              2) Calculate useful intermediate values.
c

               dtot = ( d(i,j,k) + dm(i,j,k) )*d1
               tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1
               bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3 / msolar
               isosndsp2 = sndspdC * temp(i,j,k)
               jeanmass = pi/(6._RKIND*sqrt(d(i,j,k)*dble(d1))) *
     &           dble(pi * isosndsp2 / G)**1.5_RKIND / msolar

               if (jeanmass > bmass) then
                  goto 10
               endif

c
c              3) Calculate fraction of cell mass to convert into
c                 star particles.
c

               if (tindsf .eq. 1) then
                  starfraction = masseff
               else
                  starfraction = min(masseff*dt/tdyn, 1.0_RKIND)
               endif

c
c              4) If the star particle mass is below the minimum mass,
c                 form stars stochastically such that the average SFR
c                 is correct.
c
               if ((starfraction*bmass .lt. smthresh) .or.
     &             (useodthresh .ne. 1)) then
                  call random_number(random)
                  pstar = starfraction*bmass/smthresh
                  if (random .gt. pstar) then
                     goto 10
                  else
                     starfraction = smthresh/bmass
                  endif
               endif

               starfraction = min(starfraction, mfcell)

c
c              5) Create a star particle
c
               ii = ii + 1
               mp(ii)  = starfraction*d(i,j,k)
               tcp(ii) = t
               tdp(ii) = tdyn
               xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
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
                  metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
               else
                  metalf(ii) = 0._RKIND
               endif

c              Metallicity from Type Ia SNe

               if (imetalSNIa .eq. 1) then
                  metalfSNIa(ii) = metalSNIa(i,j,k)    ! in here metal is a fraction
               endif
               if (imetalSNII .eq. 1) then
                   metalfSNII(ii) = metalSNII(i,j,k)
               endif
c
c              Remove mass from grid
c
               d(i,j,k) = (1.0_RKIND - starfraction)*d(i,j,k)

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
         write(6,*) 'star_maker_ssn: reached max new particle count'
         ERROR_MESSAGE
      endif
      np = ii
c
c      if (np .ne. 0) then
c         write(6,*) 'Stars created: number,time,level: ', np, t, level
c      endif
cc
      return
      end


c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback_ssn(nx, ny, nz,
     &                      d, dm, te, ge, u, v, w, metal,
     &                      idual, imetal, imethod, dt, r, dx, t, z,
     &                      d1, x1, v1, t1, yield,
     &                      npart, xstart, ystart, zstart, ibuff, level,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, type,
     &                      explosionFlag,smthresh,
     &                      willExplode, soonestExplosion, mu,
     &                      te1, metalSNII,
     &                      metalfSNII, imetalSNII,
     &                      s49_tot, maxlevel,
     &                      distrad, diststep, distcells)
c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Chris Loken & Greg Bryan
c  date:       3 March 1997
c  modified1:  BWO
c              13 Nov 2002
c              Many changes have been made between the date of
c              initial creation and today - all of them unlogged,
c              unfortunately.  Bugs were fixed in calculation of
c              total energy and conservation of metal and gas density.
c              The code has been cleaned up in general to enhance
c              readability.  This is the stable version - star_maker1
c              and star_maker3 should be used for experimentation.
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
c    level - current level of refinement
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    npart    - particle array size specified by calling routine
c    yield    - fraction of stellar mass that is converted to metals
c    type  - particle types (currently in the cell)
c    s49_tot - the ionizing luminosity of the particle, in units
c              of 10^49 photons per second.
c    maxlevel - the simulation's maximum refinement level
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
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, imethod, level
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    metalSNII(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, smthresh
      R_PREC    te1, mu
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      R_PREC    metalfSNII(npart)
      INTG_PREC type(npart)
      INTG_PREC explosionFlag(npart)
      INTG_PREC willExplode(npart)
      INTG_PREC imetalSNII
      R_PREC soonestExplosion(npart)
      R_PREC s49_tot(npart)
      INTG_PREC r(nx,ny,nz)
      INTG_PREC maxlevel
      INTG_PREC distrad
      INTG_PREC diststep
      INTG_PREC distcells
c
c  Locals
c
      INTG_PREC i, j, k, n, coldWind, stepk, stepj, cellstep, ic, jc, kc
      R_PREC mform, tfactor, clight, energy,
     &     m_eject, yield, minitial, xv1, xv2, dEject, ZSN,
     &     m_ejSN, dSN, td7, td7p2, td7p4, ionized, theDiff,
     &     logfac, alpha, num, stromgren_radius, stromgren_volume,
     &     cell_volume, pi, mh, dEjectSN, year, xc, yc, zc, mc,
     &     wind_energy, sn_energy, p_sn, radius, rx, ry, rz, msun, kms,
     $     ex, ey, ez, px, py, pz, kesum, mv, pcell
      parameter (clight = c_light, pi=pi_val,
     &     mh=mass_h, year=3.15569d7, msun=SolarMass,
     &     kms=1d5)
      LOGICAL typ
c
      do n=1, npart
         typ = type(n) .eq. 2 .or. type(n) .eq. 4
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. typ .and.
     &          (t - tcp(n)) .lt. 3.7e7*3.15e7/t1) then

           dEject = 0.0_RKIND


c
c
c
c         Compute index of the cell that the star particle
c           resides in.
c
            i = int((xp(n) - xstart)/dx,IKIND) + 1
            j = int((yp(n) - ystart)/dx,IKIND) + 1
            k = int((zp(n) - zstart)/dx,IKIND) + 1
c
c         check bounds - if star particle is outside of this grid
c         then something is very wrong and we should crash
c
            if (i .lt. 1 .or. i .gt. nx .or. 
     &          j .lt. 1 .or. j .gt. ny .or.
     &          k .lt. 1 .or. k .gt. nz) then
               write(6,*) 'star particle out of grid; i,j,k,level,',
     &                     i,j,k,level
               ERROR_MESSAGE
            endif

c
c          if using distributed feedback, correct if particle
c          is too close to the grid boundary
c

           if (distrad .gt. 0) then
              i = max((1 + ibuff + distrad),
     &             min((nx - ibuff - distrad), i))
              j = max((1 + ibuff + distrad),
     &              min((ny - ibuff - distrad), j))
              k = max((1 + ibuff + distrad),
     &             min((nz - ibuff - distrad), k))
           endif


c	   If the zone in which the particle sits is further refined
c           then do nothing.
           if( r(i,j,k) .eq. 0 .and. level .ne. maxlevel) then
               goto 100
           endif

c
c           Add internal energy to account for photoionization heating.
c           Energy is added in proportion to the ratio of the stromgren
c           radius to the cell volume if the stromgren sphere is smaller
c           than a cell.
c

c           Case B recombination, assuming T = 10^4 K.
            alpha = 2.60d-13 ! cm^3/s

            num = d(i,j,k) * d1 / mu / mh
            stromgren_radius = (
     &           (3._RKIND*s49_tot(n)*1d49)/(4._RKIND*pi*alpha*num**2)
     &           )**(1._RKIND/3._RKIND)
            stromgren_volume =
     &         (4._RKIND/3._RKIND)*pi*stromgren_radius**3
            cell_volume = dx*x1*dx*x1*dx*x1

            ionized = kboltz*1e4 / mu / mh / x1**2 * t1**2
            if (stromgren_volume .le. cell_volume) then
               ionized = ionized * stromgren_volume/cell_volume
            endif
            if (willExplode(n) .eq. 1 .and. ge(i,j,k) .lt. ionized) then
               theDiff = ionized-ge(i,j,k)
               ge(i,j,k) = ge(i,j,k) + theDiff
               te(i,j,k) = te(i,j,k) + theDiff
            endif


c           In preparation for using upcoming fitting formulae,
c           compute the current delay time in units of 10^7 yrs
            td7 = (t-tcp(n))*t1/3.16d14
            td7p2 = td7*td7
            td7p4 = td7p2*td7p2

c
c           Here m_eject is the mass ejection rate in log solar masses per year.
c
c           First calculate the contribution from stellar winds
c           (we'll add SNe momentarily)
c
            m_eject = 0.0_RKIND
            if (td7 .le. 0.29_RKIND) then
                m_eject = m_eject - 2.72379457924_RKIND
                m_eject = m_eject + 9.03549102928_RKIND*td7
                m_eject = m_eject - 349.073894935_RKIND*td7p2
                m_eject = m_eject + 6133.48804337_RKIND*td7p2*td7
                m_eject = m_eject - 45526.5891824_RKIND*td7p4
                m_eject = m_eject + 160159.422053_RKIND*td7p4*td7
                m_eject = m_eject - 254942.557778_RKIND*td7p4*td7p2
                m_eject = m_eject + 133227.581992_RKIND*td7p4*td7p2*td7
            else if (td7 .gt. .29_RKIND .and. td7 .le. .447_RKIND) then
                m_eject = m_eject + 1024395.67006_RKIND
                m_eject = m_eject - 22826983.8411_RKIND*td7
                m_eject = m_eject + 221695566.585_RKIND*td7p2
                m_eject = m_eject - 1225731636.28_RKIND*td7p2*td7
                m_eject = m_eject + 4219889881.74_RKIND*td7p4
                m_eject = m_eject - 9263931021.04_RKIND*td7p4*td7
                m_eject = m_eject + 12664879027.3_RKIND*td7p4*td7p2
                m_eject = m_eject - 9858823353.65_RKIND*td7p4*td7p2*td7
                m_eject = m_eject + 3345877714.47_RKIND*td7p4*td7p4
            else if (td7.ge. 0.447_RKIND .and. td7 .le. 3.24_RKIND) then
                m_eject = m_eject - 69.9540656568_RKIND
                m_eject = m_eject + 540.489990705_RKIND*td7
                m_eject = m_eject - 1785.45996729_RKIND*td7p2
                m_eject = m_eject + 3267.21333796_RKIND*td7p2*td7
                m_eject = m_eject - 3721.01017711_RKIND*td7p4
                m_eject = m_eject + 2776.66619261_RKIND*td7p4*td7
                m_eject = m_eject - 1381.69750895_RKIND*td7p4*td7p2
                m_eject = m_eject + 454.445793525_RKIND*td7p4*td7p2*td7
                m_eject = m_eject - 94.8530920783_RKIND*td7p4*td7p4
                m_eject = m_eject + 11.3775633348_RKIND*td7p4*td7p4*td7
                m_eject = m_eject -
     &               0.597112042033_RKIND*td7p4*td7p4*td7p2
            else if (td7 .ge. 3.24_RKIND .and. td7 .le. 7.28_RKIND) then
                m_eject = m_eject + 14.8304028947_RKIND
                m_eject = m_eject - 16.0726787694_RKIND*td7
                m_eject = m_eject + 4.55554172673_RKIND*td7p2
                m_eject = m_eject - 0.522521195039_RKIND*td7p2*td7
                m_eject = m_eject + 0.0214246750892_RKIND*td7p4
            else if (td7 .ge. 7.28_RKIND .and. td7 .le. 7.8_RKIND) then
                m_eject = m_eject + 3533.50268935_RKIND
                m_eject = m_eject - 1430.45985176_RKIND*td7
                m_eject = m_eject + 192.934935096_RKIND*td7p2
                m_eject = m_eject - 8.67530777079_RKIND*td7p2*td7
            else if (td7 .ge. 7.8_RKIND) then
                m_eject = m_eject - 3.5441266853_RKIND
                m_eject = m_eject - 0.0101229134154_RKIND*td7
            endif !! End of delay time cases

c           These are log specific energies with units of log((cm/s)^2)
            wind_energy = 0.0
            if ( td7 .lt. 0.29) then
                wind_energy = wind_energy + 16.7458830136_RKIND
                wind_energy = wind_energy + 2.87170884625_RKIND*td7
                wind_energy = wind_energy - 260.188160495_RKIND*td7p2
                wind_energy = wind_energy +
     &               7588.41970548_RKIND*td7p2*td7
                wind_energy = wind_energy - 109128.119673_RKIND*td7p4
                wind_energy = wind_energy +
     &               834565.297424_RKIND*td7p4*td7
                wind_energy = wind_energy -
     &               3488638.06781_RKIND*td7p4*td7p2
                wind_energy = wind_energy +
     &               7534432.3913_RKIND*td7p4*td7p2*td7
                wind_energy = wind_energy -
     &               6577304.157_RKIND*td7p4*td7p4
            else if( td7 .ge. 0.29 .and. td7 .le. 0.518) then
                wind_energy = wind_energy-8756108.90226_RKIND
                wind_energy = wind_energy+249183578.797_RKIND*td7
                wind_energy = wind_energy-3210919258.32_RKIND*td7p2
                wind_energy = wind_energy+24729596291.2_RKIND*td7p2*td7
                wind_energy = wind_energy-126485776877.0_RKIND*td7p4
                wind_energy = wind_energy+451123199767.0_RKIND*td7p4*td7
                wind_energy = wind_energy-1.14486556168d12*td7p4*td7p2
                wind_energy = wind_energy+2.067395301d12*td7p4*
     &               td7p2*td7
                wind_energy = wind_energy-2.60335249368d12*td7p4*td7p4
                wind_energy = wind_energy+2.17720380795d12*td7p4*
     &               td7p4*td7
                wind_energy = wind_energy-1.08835588273d12*td7p4*
     &               td7p4*td7p2
                wind_energy = wind_energy+246366864343.0_RKIND*td7p4*
     &               td7p4*td7p2*td7
            else if(td7 .ge. 0.518 .and. td7 .le. 2.0) then
                wind_energy = wind_energy+300.659606389_RKIND
                wind_energy = wind_energy-2175.28137376_RKIND*td7
                wind_energy = wind_energy+7038.17965731_RKIND*td7p2
                wind_energy = wind_energy-12640.7809456_RKIND*td7p2*td7
                wind_energy = wind_energy+13818.3936865_RKIND*td7p4
                wind_energy = wind_energy-9434.54106014_RKIND*td7p4*td7
                wind_energy = wind_energy+
     &               3935.34399667_RKIND*td7p4*td7p2
                wind_energy = wind_energy-
     &               918.140140181_RKIND*td7p4*td7p2*td7
                wind_energy = wind_energy+
     &               91.8268783049_RKIND*td7p4*td7p4
            else if(td7 .ge. 2.0 .and. td7 .le. 3.23) then
                wind_energy = wind_energy+1955.41904193_RKIND
                wind_energy = wind_energy-4288.5933438_RKIND*td7
                wind_energy = wind_energy+3935.44109106_RKIND*td7p2
                wind_energy = wind_energy-1921.4747372_RKIND*td7p2*td7
                wind_energy = wind_energy+526.42694795_RKIND*td7p4
                wind_energy = wind_energy-76.729393462_RKIND*td7p4*td7
                wind_energy = wind_energy+
     &               4.64818353202_RKIND*td7p4*td7p2
            else if(td7 .ge. 3.23) then
                wind_energy = wind_energy+12.6540102838_RKIND
            endif !! End of delay time cases.

c           If the most massive star in this particle has a delay
c           time greater than 22 Myr (i.e. it is low mass)
c           turn on 10 km/s wind outflows.
            coldWind = 0
            if( soonestExplosion(n) .lt. 0) then
                coldWind=1
            endif
            if( soonestExplosion(n)*t1/3.16d14 .gt. 2.20_RKIND) then
                coldWind=1
            endif

            if( coldWind .eq. 1) then
                wind_energy = 12.6540102838_RKIND
            endif

c           Convert wind_energy from cm^2/s^2 to code units
            wind_energy = 10.0_RKIND**(wind_energy -2.0_RKIND*log10(v1))

c           Find mass ejection rate in solar masses per year
            m_eject = 10.0_RKIND**m_eject

c           This ejection rate is correct for the 10^6 solar mass 'cluster'
c           used to compute these rates with starburst 99.  Reduce to account
c           for the size of the star particle

            m_eject = m_eject * smthresh/1d6

c           Convert mass ejection rate from solar masses per year to solar
c           masses per second

            m_eject = m_eject / year

c           Convert m_eject into a mass in solar masses rather than
c           a mass ejection rate by multiplying through by dt in seconds

            m_eject = m_eject * dt * t1

c           m_eject is now the mass ejected in winds in this timestep in
c           solar masses.

c           Now if there's a supernova, add the mass from that.
            m_ejSN = 0.0_RKIND
            if (explosionFlag(n) .gt. 0) then
             if (td7 .lt. 0.513_RKIND) then
                    m_ejSN = m_ejSN + 3.40965833751_RKIND
                    m_ejSN = m_ejSN - 16.0346449798_RKIND*td7
                    m_ejSN = m_ejSN + 31.5091825735_RKIND*td7p2
                    m_ejSN = m_ejSN - 21.3218283568_RKIND*td7p2*td7
             else if(.513_RKIND .le. td7 .and. td7 .le. .918_RKIND) then
                    m_ejSN = m_ejSN - 314538.854117_RKIND
                    m_ejSN = m_ejSN + 4453582.08399_RKIND*td7
                    m_ejSN = m_ejSN - 28218211.3741_RKIND*td7p2
                    m_ejSN = m_ejSN + 105370186.068_RKIND*td7p2*td7
                    m_ejSN = m_ejSN - 256824281.305_RKIND*td7p4
                    m_ejSN = m_ejSN + 426986197.681_RKIND*td7p4*td7
                    m_ejSN = m_ejSN - 490461521.485_RKIND*td7p4*td7p2
                    m_ejSN = m_ejSN +
     &                   384394390.035_RKIND*td7p4*td7p2*td7
                    m_ejSN = m_ejSN - 196752045.251_RKIND*td7p4*td7p4
                    m_ejSN = m_ejSN +
     &                   59399337.5861_RKIND*td7p4*td7p4*td7
                    m_ejSN = m_ejSN -
     &                   8033095.66643_RKIND*td7p4*td7p4*td7p2
             else if (0.918_RKIND .le. td7 .and. td7 .le. 3.23) then
                    m_ejSN = m_ejSN + 1.74261906723_RKIND
                    m_ejSN = m_ejSN - 0.92589554122_RKIND*td7
                    m_ejSN = m_ejSN + 0.551250718292_RKIND*td7p2
                    m_ejSN = m_ejSN - 0.220085806978_RKIND*td7p2*td7
                    m_ejSN = m_ejSN + 0.0510662546479_RKIND*td7p4
                    m_ejSN = m_ejSN - 0.00504400687495_RKIND*td7p4*td7
             else if (td7 .ge. 3.23_RKIND) then
                    m_ejSN = m_ejSN + 2.67991943387_RKIND
                    m_ejSN = m_ejSN - 0.461075452846_RKIND*td7
                    m_ejSN = m_ejSN - 0.0326899620754_RKIND*td7p2
             endif !! End of delay time cases
             m_ejSN = 10.0_RKIND**m_ejSN * explosionFlag(n)
c           Current time (code units), Delay time (10^7 yrs),
c           Mass ejected by SN (MSun), Mass ejected by winds (MSun)
c           Particle mass (code density), x,y,z of particle (code)
c           Density, Gas energy, Total energy in cell prior to SN (code)
c           and finally sp. energy of stellar wind (code units)
c             write(6,*) 'Supernova! ', t,td7,m_ejSN,m_eject,
c     &         mp(n), xp(n),yp(n),zp(n),d(i,j,k),ge(i,j,k),
c     &            te(i,j,k),wind_energy
            endif !! End of explosionFlag>0

            dEjectSN = m_ejSN * 1.989d33 / (dx*x1)**3 / d1
            dEject = m_eject * 1.989d33 / (dx*x1)**3 / d1

            if((dEject + dEjectSN) .gt. mp(n)) then
                write(6,*) 'WARNING: losing too much mass'
                write(6,*) 'Setting particle mass to zero'
                write(6, *) mp(n), dEject, dEjectSN
                write(6,*) 'm_eject, m_ejSN, td7 ',m_eject,m_ejSN,td7
                write(6,*)  willExplode(n), explosionFlag(n),i,j,k,n
                mp(n)=0.0_RKIND
            else
                mp(n) = mp(n) - (dEject + dEjectSN)
            endif

c           sn_energy is in code units of energy

            if (explosionFlag(n) .gt. 0) then
               sn_energy = explosionFlag(n) *
     &              (1.0d51/(d1*x1**3)/v1/v1)
            else
               sn_energy = 0
            endif  !! end of explosionFlag>0
c           Now we spread the energy out over the whole cell:
            energy = wind_energy*dEject*dx**3 + sn_energy
            energy = energy / ((d(i,j,k) + dEject + dEjectSN)*dx**3)

c           Add wind and SN energy to the cell the particle lives in.

            te(i,j,k) = te(i,j,k) + energy
            if (idual .eq. 1)
     &           ge(i,j,k) = ge(i,j,k) + energy

c
c            Metal feedback (note that in this function gas metal is
c            a fraction (rho_metal/rho_gas) rather than a density.
c            The conversion has been done in the handling routine)
c
            if (imetal .eq. 1) then

            ZSN = 0.0_RKIND
            if(explosionFlag(n) .gt. 0) then
             if (td7 .le. 0.365_RKIND) then
                    ZSN = ZSN - 2.12788803787_RKIND
                    ZSN = ZSN - 3.08562080661_RKIND * td7
                    ZSN = ZSN + 83.9782331967_RKIND * td7p2
                    ZSN = ZSN - 158.867576409_RKIND * td7p2*td7
             else if(.365_RKIND .le. td7 .and. td7 .le. .442_RKIND) then
                    ZSN = ZSN + 0.584246828515_RKIND
                    ZSN = ZSN - 10.2313474777_RKIND * td7
                    ZSN = ZSN + 42.3869874214_RKIND * td7p2
                    ZSN = ZSN - 46.7483110349_RKIND * td7p2*td7
            else if(.442_RKIND .le. td7 .and. td7 .le. .517_RKIND) then
                    ZSN = ZSN - 6913398.58591_RKIND
                    ZSN = ZSN + 99253357.5895_RKIND * td7
                    ZSN = ZSN - 610064398.162_RKIND * td7p2
                    ZSN = ZSN + 2081033583.17_RKIND * td7p2*td7
                    ZSN = ZSN - 4254712179.17_RKIND * td7p4
                    ZSN = ZSN + 5213617301.35_RKIND * td7p4*td7
                    ZSN = ZSN - 3545287887.67_RKIND * td7p4*td7p2
                    ZSN = ZSN + 1032027574.24_RKIND * td7p4*td7p2*td7
             else if(.517_RKIND .le. td7 .and. td7 .le. .718_RKIND) then
                    ZSN = ZSN - 0.856562917033_RKIND
                    ZSN = ZSN + 2.72919708256_RKIND *td7
                    ZSN = ZSN - 1.22818002617_RKIND *td7p2
                    ZSN = ZSN - 0.536735943166_RKIND *td7p2*td7
             else if(.718_RKIND .le. td7 .and. td7 .le. 1.7_RKIND) then
                    ZSN = ZSN + 6.69150908774_RKIND
                    ZSN = ZSN - 51.5159032945_RKIND * td7
                    ZSN = ZSN + 172.980768687_RKIND * td7p2
                    ZSN = ZSN - 307.64633271_RKIND  * td7p2*td7
                    ZSN = ZSN + 311.333072359_RKIND * td7p4
                    ZSN = ZSN - 180.28789728_RKIND  * td7p4*td7
                    ZSN = ZSN + 55.6987751625_RKIND * td7p4*td7p2
                    ZSN = ZSN - 7.12566268516_RKIND * td7p4*td7p2*td7
             else if(1.7_RKIND .le. td7 .and. td7 .le. 3.2_RKIND) then
                    ZSN = ZSN + 0.173226556725_RKIND
                    ZSN = ZSN - 0.142758066129_RKIND * td7
                    ZSN = ZSN + 0.0568361245745_RKIND * td7p2
                    ZSN = ZSN - 0.00740131298973_RKIND *td7p2*td7
             else
                    ZSN = ZSN - 0.459335515363_RKIND
                    ZSN = ZSN + 0.51236004354_RKIND * td7
                    ZSN = ZSN - 0.194211146544_RKIND* td7p2
                    ZSN = ZSN + 0.0264297153731_RKIND *td7p2*td7
             endif !! End of delaytime cases
            endif !! End of explosionFlag>0


            metal(i,j,k) = metal(i,j,k) * d(i,j,k)
            metal(i,j,k) = metal(i,j,k) + metalf(n) * dEject
            metal(i,j,k) = metal(i,j,k) + ZSN * dEjectSN
            metal(i,j,k) = metal(i,j,k) / (dEject + dEjectSN + d(i,j,k))
            if (imetalSNII .eq. 1) then
               metalSNII(i,j,k) = metalSNII(i,j,k) * d(i,j,k)
               metalSNII(i,j,k) = metalSNII(i,j,k) +
     &              (metalfSNII(n) * dEject)
               metalSNII(i,j,k) = metalSNII(i,j,k) + ZSN * dEjectSN
               metalSNII(i,j,k) = metalSNII(i,j,k)
     &             / (dEject + dEjectSN + d(i,j,k))
            endif !! End of imetalSNII==1

         endif                  !! End of imetal==1


c
c              Mass and momentum feedback from winds and SN ejecta (NOT from SN blastwaves)
c
           u(i,j,k) = u(i,j,k)*d(i,j,k) + (dEject + dEjectSN) * up(n)
           v(i,j,k) = v(i,j,k)*d(i,j,k) + (dEject + dEjectSN) * vp(n)
           w(i,j,k) = w(i,j,k)*d(i,j,k) + (dEject + dEjectSN) * wp(n)
           d(i,j,k) = d(i,j,k) + dEject + dEjectSN
           u(i,j,k) = u(i,j,k)/d(i,j,k)
           v(i,j,k) = v(i,j,k)/d(i,j,k)
           w(i,j,k) = w(i,j,k)/d(i,j,k)
c
c              If te is really total energy (and it is unless imethod=2),
c              then just set this value
c
           if (imethod .ne. 2 .and. idual .eq. 1) then
              te(i,j,k) = 0.5_RKIND*(u(i,j,k)**2 + v(i,j,k)**2 +
     &             w(i,j,k)**2) + ge(i,j,k)
           endif

c
c          do not do distributed feedback if there are no SNe
c          or if it is not turned on
c

           if (sn_energy .eq. 0) then
              goto 100
           endif

           if  (distrad .eq. 0) then
              goto 100
           endif

c
c          Find momentum to deposit per SNe in code units.
c

           p_sn = explosionFlag(n) * 3d5 * msun / (d1*x1**3) * 
     &      kms / v1 / (distcells - 1)

c          Approximate that the explosion happens in the center of the
c          feedback zone

           ex = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
           ey = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
           ez = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx

c
c          do momentum feedback over a distributed region
c          distribute equally among zones participating in feedback
c          Limit velocity after feedback to 1000 km/s to avoid accelerating
c          material in cells with very little mass
c

           kesum = 0
           mv = 1d3*kms/v1

           do kc = k-distrad,k+distrad
              stepk = abs(kc-k)
              do jc = j-distrad,j+distrad
                 stepj = stepk + abs(jc-j)
                 do ic = i-distrad,i+distrad
                    cellstep = stepj + abs(ic-i)
                    if (cellstep .le. diststep) then

c                      Don't do momentum feedback in the cell the SN goes off in

                       if ((ic .eq. i) .and. (jc .eq. j) .and.
     &                      (kc .eq. k)) then
                          continue
                       else

c                         Define some useful intermediate variables

                          xc = xstart + (REAL(ic,RKIND)-0.5_RKIND)*dx
                          yc = ystart + (REAL(jc,RKIND)-0.5_RKIND)*dx
                          zc = zstart + (REAL(kc,RKIND)-0.5_RKIND)*dx

                          radius = sqrt((xc - ex)**2 + (yc - ey)**2
     &                         + (zc - ez)**2)

                          rx = (xc - ex)/radius
                          ry = (yc - ey)/radius
                          rz = (zc - ez)/radius

                          mc = d(ic,jc,kc)*dx**3

c                         Limit momentum feedback such that each cell can do
c                         no more than its fair share of the thermal energy
c                         from the SNe

                          if (p_sn**2/(2.0_RKIND*mc) .gt. 
     &                         sn_energy/(distcells-1)) then
                             pcell=2.0_RKIND*mc*sn_energy/(distcells-1)
                             pcell=SQRT(pcell)
                          else
                             pcell = p_sn
                          endif

c                         Limit maximum velocity boost to 1000 km/s

                          if ((p_sn/mc) .gt. mv) then
                             pcell = mc*mv
                          endif

c                         Add momentum feedback

                          u(ic,jc,kc) = u(ic,jc,kc) + (rx*pcell)/mc
                          v(ic,jc,kc) = v(ic,jc,kc) + (ry*pcell)/mc
                          w(ic,jc,kc) = w(ic,jc,kc) + (rz*pcell)/mc

                          if (imethod .ne. 2 .and. idual .eq. 1) then
                             te(ic,jc,kc) = ge(ic,jc,kc) + 0.5_RKIND*(
     &                            u(ic,jc,kc)**2 +
     &                            v(ic,jc,kc)**2 +
     &                            w(ic,jc,kc)**2)
                          endif

c                         Record kinetic energy deposited in this cell

                          kesum = kesum + 0.5_RKIND/mc*pcell**2
                       endif
                    endif
                 enddo
              enddo
           enddo

           if (kesum .gt. sn_energy) then
              write(6,*) 'ssn1: Kinetic energy from momentum feedback
     & exceeds thermal energy', kesum, sn_energy, wind_energy, energy,
     &             ge(i,j,k), te(i,j,k), kesum/(d(i,j,k)*dx**3)
              kesum = sn_energy
           endif

           if (kesum/(d(i, j, k)*dx**3) .ge. ge(i,j,k)) then
              write(6,*) 'ssn2: Kinetic energy from momentum feedback
     & exceeds thermal energy', kesum, sn_energy, wind_energy, energy,
     &             ge(i,j,k), te(i,j,k), kesum/(d(i,j,k)*dx**3)
              ERROR_MESSAGE
           endif

c          Assume kinetic energy from momentum feedback is extracted from
c          thermal energy in the host zone

           ge(i, j, k) = ge(i, j, k) - kesum/(d(i, j, k)*dx**3)
           te(i, j, k) = te(i, j, k) - kesum/(d(i, j, k)*dx**3)

c
c
c
 10         continue
         endif !! End of check to see if this is a star particle
c
 100     continue
c
      enddo
c
c      write(6,*) 'star_feedback7: end'
 200  continue

      return
      end

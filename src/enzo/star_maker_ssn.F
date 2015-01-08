#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c

      subroutine star_maker_ssn(nx, ny, nz,
     &     d, dm, temp, u, v, w, cooltime,
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
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz), metal(nx,ny,nz)
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
      integer*4 n
      integer*4, allocatable :: seed(:)
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
c     Use (MPI rank + 1) to seed the RNG, ensuring different streams
c     on each processor.  Offset by one to avoid the poor RNG stream
c     generated when the seed is set to 0.
      if(ran1_init .eq. 0) then
         call random_seed(size=n)
         allocate(seed(n))
         do i = 1, n
            seed(i) = procnum+1
         end do
         call random_seed(put=seed)
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


c               if (useodthresh .ne. 1) then
c	          Make sure we never use more than half the mass in the cell to create a particle.
c                 One can think of this as an implicit resolution-dependent density threshold.
                  if (mfcell*bmass .lt. smthresh) then
                     goto 10
                  endif
c               endif

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


      function typical_mass(td7)
c
c   This is a fitting function to the results of a starburst99
c    calculation. The input, td7, is the delay time in units of
c    10^7 years. The output, typical_mass, is the mass in solar
c    masses of a star that explodes with the given delay time. 
c   The intermediate quantity td77 = (td7)^7
      R_PREC typical_mass,td7,td77
      typical_mass = 0
      td77 = td7*td7*td7*td7*td7*td7*td7
      typical_mass = typical_mass + 4.42035634891
      typical_mass = typical_mass - 13.2890466089*td7
      typical_mass = typical_mass + 26.1103296098*td7*td7
      typical_mass = typical_mass - 30.1876007562*td7*td7*td7
      typical_mass = typical_mass + 21.8976126631*td7*td7*td7*td7
      typical_mass = typical_mass -10.2544493943*td7*td7*td7*td7*td7
      typical_mass = typical_mass+3.09621304958*td7*td7*td7*td7*td7*td7
      typical_mass = typical_mass-0.581870413299*td77
      typical_mass = typical_mass+0.0618795946119*td77*td7
      typical_mass = typical_mass-0.00284352366041*td77*td7*td7

      typical_mass = 10**typical_mass
      return
      end

c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback_ssn(nx, ny, nz,
     &                      d, dm, te, ge, u, v, w, metal,
     &                      idual, imetal, imethod, dt, r, dx, t, z,
     &                      d1, x1, v1, t1, sn_param, retfr, yield,
     &                      npart, xstart, ystart, zstart, ibuff, level,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, type, 
     &                      explosionFlag,smthresh,
     &                      willExplode, soonestExplosion, gam, mu, 
     &                      te1, metalSNII,
     &                      metalfSNII, imetalSNII,
     &                      s49_tot, maxlevel)
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
c    sn_param - fraction of stellar rest mass that goes to feedback
c    retfr     - fraction of mass ejected in winds [not used!].
c               used to be m_eject
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
      R_PREC    dt, dx, z, retfr
      R_PREC    d1, x1, v1, t1, smthresh
      R_PREC    te1, gam, mu
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
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      INTG_PREC i, j, k, n, coldWind
      R_PREC mform, tfactor, clight, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, dratio, dEject, ZSN,
     &     m_ejSN, dSN, td7, td7p2, td7p4, ionized, theDiff, q2,
     &     logfac, alpha, num, stromgren_radius, stromgren_volume,
     &     cell_volume, pi, mh
      parameter (clight = c_light, msolar_e51 = 1800._RKIND,
     &           pi=pi_val, mh=mass_h)
      LOGICAL typ, cond
      R_PREC dbe, dmi, den, geb, teb, gem, tem
c
      do n=1, npart
         typ = type(n) .eq. 2 .or. type(n) .eq. 4
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. typ) then

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
c         then exit and give a warning.
c
            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &          .or. k .lt. 1 .or. k .gt. nz) then
               write(6,*) 'star particle out of grid; i,j,k,level,',
     &                     i,j,k,level
               goto 100
            endif

c	   If the zone in which the particle sits is further refined
c           then do nothing.
           if( r(i,j,k) .eq. 0 .and. level .ne. maxlevel) then
               goto 100
           endif

           geb = ge(i,j,k)
           teb = te(i,j,k)
           dbe = d(i,j,k)
           cond = te(i,j,k) .gt. 100_RKIND .or. ge(i,j,k) .gt. 100_RKIND
           if( cond ) then
	     write(6,*) 'ssn low density 395: ',d(i,j,k),dEject,level, 
     &           willExplode(n), explosionFlag(n),i,j,k,n,dbe, geb,teb
           endif

c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c

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
            ionized = 3.96e-5_RKIND
            if (stromgren_volume .le. cell_volume) then
               ionized = ionized * stromgren_volume/cell_volume
            endif
            if (willExplode(n) .eq. 1 .and. ge(i,j,k) .lt. ionized) then
               theDiff = ionized-ge(i,j,k)
               ge(i,j,k) = ge(i,j,k) + theDiff
               te(i,j,k) = te(i,j,k) + theDiff
            endif

           cond = te(i,j,k) .gt. 100_RKIND .or. ge(i,j,k) .gt. 100_RKIND
           if( cond ) then
	     write(6,*) 'ssn low density 458: ',d(i,j,k),dEject,level, 
     &           willExplode(n), explosionFlag(n),i,j,k,n,dbe, geb,teb
           endif


c           In preparation for using upcoming fitting formulae,
c           compute the current delay time in units of 10^7 yrs
            td7 = (t-tcp(n))*t1/3.16e14
            td7p2 = td7*td7
            td7p4 = td7p2*td7p2

c           The mass ejected in Msun. First calculate the contr from
c               stellar winds (we'll add SNe momentarily)            
            m_eject = 0.0_RKIND
            if (td7 .le. 0.29_RKIND) then
                m_eject = m_eject - 2.72379457924 
                m_eject = m_eject + 9.03549102928*td7
                m_eject = m_eject - 349.073894935*td7p2
                m_eject = m_eject + 6133.48804337*td7p2*td7
                m_eject = m_eject - 45526.5891824*td7p4
                m_eject = m_eject + 160159.422053*td7p4*td7
                m_eject = m_eject - 254942.557778*td7p4*td7p2
                m_eject = m_eject + 133227.581992*td7p4*td7p2*td7
            else if (td7 .gt. .29_RKIND .and. td7 .le. .447_RKIND) then
                m_eject = m_eject + 1024395.67006
                m_eject = m_eject - 22826983.8411*td7
                m_eject = m_eject + 221695566.585*td7p2
                m_eject = m_eject - 1225731636.28*td7p2*td7
                m_eject = m_eject + 4219889881.74*td7p4
                m_eject = m_eject - 9263931021.04*td7p4*td7
                m_eject = m_eject + 12664879027.3*td7p4*td7p2
                m_eject = m_eject - 9858823353.65*td7p4*td7p2*td7
                m_eject = m_eject + 3345877714.47*td7p4*td7p4
            else if (td7.ge. 0.447_RKIND .and. td7 .le. 3.24_RKIND) then
                m_eject = m_eject - 69.9540656568
                m_eject = m_eject + 540.489990705*td7
                m_eject = m_eject - 1785.45996729*td7p2
                m_eject = m_eject + 3267.21333796*td7p2*td7
                m_eject = m_eject - 3721.01017711*td7p4
                m_eject = m_eject + 2776.66619261*td7p4*td7
                m_eject = m_eject - 1381.69750895*td7p4*td7p2
                m_eject = m_eject + 454.445793525*td7p4*td7p2*td7
                m_eject = m_eject - 94.8530920783*td7p4*td7p4
                m_eject = m_eject + 11.3775633348*td7p4*td7p4*td7
                m_eject = m_eject - 0.597112042033*td7p4*td7p4*td7p2
            else if (td7 .ge. 3.24_RKIND .and. td7 .le. 7.28_RKIND) then
                m_eject = m_eject + 14.8304028947
                m_eject = m_eject - 16.0726787694*td7
                m_eject = m_eject + 4.55554172673*td7p2
                m_eject = m_eject - 0.522521195039*td7p2*td7
                m_eject = m_eject + 0.0214246750892*td7p4
            else if (td7 .ge. 7.28_RKIND .and. td7 .le. 7.8_RKIND) then
                m_eject = m_eject + 3533.50268935
                m_eject = m_eject - 1430.45985176*td7
                m_eject = m_eject + 192.934935096*td7p2
                m_eject = m_eject - 8.67530777079*td7p2*td7
            else if (td7 .ge. 7.8_RKIND) then
                m_eject = m_eject - 3.5441266853
                m_eject = m_eject - 0.0101229134154*td7
            endif !! End of delay time cases

           if( d(i,j,k) .lt. 25.0) then
	     write(6,*) 'ssn low density 509: ',d(i,j,k),dEject,level, 
     &           willExplode(n), explosionFlag(n),i,j,k,n,dbe
           endif

c           At this point, energy is specific energy, i.e. vel sq
            energy = 0.0
            if ( td7 .lt. 0.29) then
                energy = energy + 16.7458830136
                energy = energy + 2.87170884625*td7
                energy = energy - 260.188160495*td7p2
                energy = energy + 7588.41970548*td7p2*td7
                energy = energy - 109128.119673*td7p4
                energy = energy + 834565.297424*td7p4*td7
                energy = energy - 3488638.06781*td7p4*td7p2
                energy = energy + 7534432.3913*td7p4*td7p2*td7
                energy = energy - 6577304.157*td7p4*td7p4
            else if( td7 .ge. 0.29 .and. td7 .le. 0.518) then
                energy = energy-8756108.90226
                energy = energy+249183578.797*td7
                energy = energy-3210919258.32*td7p2
                energy = energy+24729596291.2*td7p2*td7
                energy = energy-126485776877.0*td7p4
                energy = energy+451123199767.0*td7p4*td7
                energy = energy-1.14486556168e12*td7p4*td7p2
                energy = energy+2.067395301e+12*td7p4*td7p2*td7
                energy = energy-2.60335249368e+12*td7p4*td7p4
                energy = energy+2.17720380795e+12*td7p4*td7p4*td7
                energy = energy-1.08835588273e+12*td7p4*td7p4*td7p2
                energy = energy+246366864343.0*td7p4*td7p4*td7p2*td7
            else if(td7 .ge. 0.518 .and. td7 .le. 2.0) then
                energy = energy+300.659606389
                energy = energy-2175.28137376*td7
                energy = energy+7038.17965731*td7p2
                energy = energy-12640.7809456*td7p2*td7
                energy = energy+13818.3936865*td7p4
                energy = energy-9434.54106014*td7p4*td7
                energy = energy+3935.34399667*td7p4*td7p2
                energy = energy-918.140140181*td7p4*td7p2*td7
                energy = energy+91.8268783049*td7p4*td7p4
            else if(td7 .ge. 2.0 .and. td7 .le. 3.23) then
                energy = energy+1955.41904193
                energy = energy-4288.5933438*td7
                energy = energy+3935.44109106*td7p2
                energy = energy-1921.4747372*td7p2*td7
                energy = energy+526.42694795*td7p4
                energy = energy-76.729393462*td7p4*td7
                energy = energy+4.64818353202*td7p4*td7p2
            else if(td7 .ge. 3.23) then
                energy = energy+12.6540102838
            endif !! End of delay time cases.

            coldWind = 0
            if( soonestExplosion(n) .lt. 0) then
                coldWind=1
            endif
            if( soonestExplosion(n)*t1/3.16e14 .gt. 2.20) then
                coldWind=1
            endif

            if( coldWind .eq. 1) then
                energy = 12.6540102838
            endif

           if( d(i,j,k) .lt. 25.0) then
	     write(6,*) 'ssn low density 574: ',d(i,j,k),dEject,level, 
     &           willExplode(n), explosionFlag(n),i,j,k,n,dbe
           endif

            logfac = 0.0
            if (energy .gt. 15.416685) then
                logfac = energy - 15.416685
                energy = energy - logfac
            endif
c           Convert energy from cm^2/s^2 to code units            
            energy = 10.0**(energy -2.0*log10(v1))



            m_eject = m_eject + logfac

c           Using the current dt, find the mass ejected this timestep            
c           m_eject is the log base 10 of the rate in Msun/yr
c           multiply by dt in seconds, divide by seconds per year
c           and multply by the ratio of the star particle mass to 
c           10^6 solar masses, the mass of the 'cluster' used to
c           compute these rates with starburst 99.
            m_eject = 10.0**m_eject * dt * t1 * smthresh/ 3.16e13
           
c           Now if there's a supernova, add the mass from that.
            m_ejSN = 0.0
            if (explosionFlag(n) .gt. 0) then
             if (td7 .lt. 0.513_RKIND) then
                    m_ejSN = m_ejSN + 3.40965833751
                    m_ejSN = m_ejSN - 16.0346449798*td7
                    m_ejSN = m_ejSN + 31.5091825735*td7p2
                    m_ejSN = m_ejSN - 21.3218283568*td7p2*td7
             else if(.513_RKIND .le. td7 .and. td7 .le. .918_RKIND) then
                    m_ejSN = m_ejSN - 314538.854117
                    m_ejSN = m_ejSN + 4453582.08399*td7
                    m_ejSN = m_ejSN - 28218211.3741*td7p2
                    m_ejSN = m_ejSN + 105370186.068*td7p2*td7
                    m_ejSN = m_ejSN - 256824281.305*td7p4
                    m_ejSN = m_ejSN + 426986197.681*td7p4*td7
                    m_ejSN = m_ejSN - 490461521.485*td7p4*td7p2
                    m_ejSN = m_ejSN + 384394390.035*td7p4*td7p2*td7
                    m_ejSN = m_ejSN - 196752045.251*td7p4*td7p4
                    m_ejSN = m_ejSN + 59399337.5861*td7p4*td7p4*td7
                    m_ejSN = m_ejSN - 8033095.66643*td7p4*td7p4*td7p2
             else if (0.918_RKIND .le. td7 .and. td7 .le. 3.23) then
                    m_ejSN = m_ejSN + 1.74261906723
                    m_ejSN = m_ejSN - 0.92589554122*td7
                    m_ejSN = m_ejSN + 0.551250718292*td7p2
                    m_ejSN = m_ejSN - 0.220085806978*td7p2*td7
                    m_ejSN = m_ejSN + 0.0510662546479*td7p4
                    m_ejSN = m_ejSN - 0.00504400687495*td7p4*td7
             else if (td7 .ge. 3.23_RKIND) then
                    m_ejSN = m_ejSN + 2.67991943387
                    m_ejSN = m_ejSN - 0.461075452846*td7
                    m_ejSN = m_ejSN - 0.0326899620754*td7p2
             endif !! End of delay time cases
             m_ejSN = 10.0**m_ejSN * explosionFlag(n)
c           Current time (code units), Delay time (10^7 yrs), 
c           Mass ejected by SN (MSun), Mass ejected by winds (MSun)
c           Particle mass (code density), x,y,z of particle (code)
c           Density, Gas energy, Total energy in cell prior to SN (code)
c           and finally sp. energy of stellar wind (code units)
             write(6,*) 'Supernova! ', t,td7,m_ejSN,m_eject,
     &         mp(n), xp(n),yp(n),zp(n), d(i,j,k),ge(i,j,k),
     &         te(i,j,k),energy
             m_eject = m_eject + m_ejSN
            endif !! End of explosionFlag>0


           if( d(i,j,k) .lt. 25.0) then
	     write(6,*) 'ssn low density 643: ',d(i,j,k),dEject,level, 
     &           willExplode(n), explosionFlag(n),i,j,k,n,dbe
           endif

            dSN = m_ejSN * 1.989e33 / (dx*x1)**3 / d1
            dEject = m_eject * 1.989e33 / (dx*x1)**3 / d1
            if(dEject .gt. mp(n)) then
                write(6,*) 'WARNING: losing too much mass',mp(n),dEject
                write(6,*) 'm_eject, m_ejSN, td7 ',m_eject,m_ejSN,td7
                write(6,*)  willExplode(n), explosionFlag(n),i,j,k,n
            endif
            mp(n) = mp(n) - dEject

            if (explosionFlag(n) .gt. 0) then
c               Add wind and SN energy. Need to weight by their
c               respective masses.
                energy = energy*(m_eject-m_ejSN*explosionFlag(n)) 
                energy = energy + 1.0e18 / 1.989 /v1/v1
                energy = energy / m_eject
            endif  !! end of explosionFlag>0
c           Now we spread the energy out over the whole cell:
            energy = energy * dEject/ (d(i,j,k) + dEject)
c            energy = sn_param * mform * (clight/v1)**2 / 
c     &                (d(i,j,k)+mform*m_eject)
c
c
c
c
c            write(6,*) 'star_feedback7: mp, minitial, mform, energy',
c     &                  mp(n), minitial, mform, energy
c
c
c              Add energy to energy field
c
               dratio = d(i,j,k)/(d(i,j,k) + dEject)
c               write(6,*) 'ssn: ',dratio,te(i,j,k),ge(i,j,k),energy,t,td7, 
c     &          dEject,m_ejSN,d(i,j,k),t1,m_eject,dEject,dx,dt
               te(i,j,k) = te(i,j,k)*dratio + energy
               if (idual .eq. 1) 
     &            ge(i,j,k) = ge(i,j,k)*dratio + energy


           gem = (ge(i,j,k) - energy)/dratio
           tem = (te(i,j,k) - energy)/dratio
           cond = te(i,j,k) .gt. 100_RKIND .or. ge(i,j,k) .gt. 100_RKIND
           if( cond ) then
	     write(6,*) 'ssn low density 686: ',d(i,j,k),dEject,level, 
     &           willExplode(n), explosionFlag(n),i,j,k,n,dbe,energy,
     &             te(i,j,k),ge(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),
     &             dRatio, td7, geb, teb, gem, tem
           endif


c
c              Metal feedback (note that in this function gas metal is
c              a fraction (rho_metal/rho_gas) rather than a density.
c              The conversion has been done in the handling routine)
c
               if (imetal .eq. 1) then
c
c              "Cen method".  This takes into account gas recycling.
c
c               To be more specific, it looks like we add up metals from
c               1) The gas in the cell, 2) metals instantaneously
c               produced this time step by star formation, and 3) the
c               metals which were originally in the gas when the star
c               formed and are now being released in proportion to 
c               the mass of stars formed this time step. Let's not do
c               this.
c                  metal(i,j,k) = (metal(i,j,k)*d(i,j,k) + mform * 
c     &                 (yield * (1._RKIND-metalf(n)) + m_eject * 
c     &                 metalf(n))) / (d(i,j,k)+dEject)  ! metal is a fraction
c




            ZSN = 0.0
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
                    ZSN = ZSN - 0.856562917033
                    ZSN = ZSN + 2.72919708256 *td7
                    ZSN = ZSN - 1.22818002617 *td7p2
                    ZSN = ZSN - 0.536735943166 *td7p2*td7
             else if(.718_RKIND .le. td7 .and. td7 .le. 1.7_RKIND) then
                    ZSN = ZSN + 6.69150908774
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


           if( d(i,j,k) .lt. 25.0) then
	     write(6,*) 'ssn low density 763: ',d(i,j,k),dEject,level, 
     &           willExplode(n), explosionFlag(n),i,j,k,n,dbe
           endif


            q2 = dEject - dSN
            metal(i,j,k) = metal(i,j,k) * d(i,j,k) 
            metal(i,j,k) = metal(i,j,k) +   metalf(n) * q2
            metal(i,j,k) = metal(i,j,k) +   ZSN * dSN 
            metal(i,j,k) = metal(i,j,k)  / (dEject + d(i,j,k))
            if (imetalSNII .eq. 1) then
             metalSNII(i,j,k) = metalSNII(i,j,k) * d(i,j,k)
             metalSNII(i,j,k) = metalSNII(i,j,k) + metalfSNII(n) *q2
             metalSNII(i,j,k) = metalSNII(i,j,k) + ZSN * dSN
             metalSNII(i,j,k) = metalSNII(i,j,k) / (dEject + d(i,j,k))
            endif  !! End of imetalSNII==1

           endif !! End of imetal==1


               dmi = d(i,j,k)
               if( d(i,j,k) .lt. 25.0) then
		 write(6,*) 'ssn low density 787: ',d(i,j,k),dEject,level,
     &             willExplode(n),explosionFlag(n),i,j,k,n,dbe,dmi,
     &             metal(i,j,k), metalSNII(i,j,k)
               endif


c
c              Mass and momentum feedback
c
               u(i,j,k) = u(i,j,k)*d(i,j,k) + dEject * up(n)
               v(i,j,k) = v(i,j,k)*d(i,j,k) + dEject * vp(n)
               w(i,j,k) = w(i,j,k)*d(i,j,k) + dEject * wp(n)
               d(i,j,k) = d(i,j,k) + dEject
               u(i,j,k) = u(i,j,k)/d(i,j,k)
               v(i,j,k) = v(i,j,k)/d(i,j,k)
               w(i,j,k) = w(i,j,k)/d(i,j,k)
c
c              If te is really total energy (and it is unless imethod=2),
c              then just set this value
c
               if (imethod .ne. 2 .and. idual .eq. 1) then
                  te(i,j,k) = 0.5_RKIND*(u(i,j,k)**2 + v(i,j,k)**2 + 
     &                               w(i,j,k)**2) + ge(i,j,k)
               endif


	       den = d(i,j,k)
               cond = te(i,j,k) .gt. 100_RKIND
               cond = cond  .or. ge(i,j,k) .gt. 100_RKIND
               if( cond ) then
		 write(6,*) 'ssn low density 761: ',d(i,j,k),dEject,level,
     &             willExplode(n), explosionFlag(n),i,j,k,n,dbe,dmi,den,
     &             te(i,j,k),ge(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),
     &             energy, dRatio, td7
               endif
		
c
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

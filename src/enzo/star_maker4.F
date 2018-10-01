#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker4(nx, ny, nz,
     &                      d, dm, temp, u, v, w,
     &                      dt, r, metal, dx, t, z, 
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, mintdyn,
     &                      odthresh, smthresh, level, np,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      imetalSNIa, metalSNIa, metalfSNIa)

c
c  CREATES STAR PARTICLES
c
c  written by:  Brian O'Shea
c  date:  29 March 2005
c    this version of star_maker creates star particles using the method
c    of Andrey Kravtsov as described in "On the origin of the global
c    Schmidt law of star formation", 2003 ApJ 590:L1-L4.  As in this paper,
c    I hold the star formation timescale constant and have a fixed 
c    overdensity threshold.
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
c    odthresh - overdensity threshold (some number * avg. density)
c    smthresh - star mass threshold (only creates stars with mass >
c        smthresh unless (random number) < starmass/smthresh )
c    mintdyn  - minimum dynamical time, in years
c    level - current level of refinement
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
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      INTG_PREC imetalSNIa
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      R_PREC    metalSNIa(nx,ny,nz), metalfSNIa(nmax)
      R_PREC    odthresh, smthresh, mintdyn
c
      R_PREC   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      INTG_PREC  i, j, k, ii
      R_PREC   div, tdyn, dtot
      R_PREC   pi, G, sndspdC,mproton
      R_PREC   isosndsp2, starmass, starfraction, bmass, jeanmass
      R_PREC   densthresh, timeconstant, gasfrac
      real*8 msolar
      parameter (pi=pi_val, G=GravConst, 
     &           sndspdC=1.3095e8_RKIND,
     &           msolar=SolarMass,
     &           mproton=mass_h)
c
      ii = np
c
c  calculate density threshold.  odthresh is in proper particles per cc
c  and (d1/mproton) gives the mean density of the universe in 
c  particles/cm^3 (assuming hydrogen is dominant)
c
c
      densthresh = odthresh / (d1 / mproton)

c
c  calculate time constant for star formation.  This assumes that the
c  user input is in units of years
c

      timeconstant = mintdyn * 3.156e7_RKIND / t1
c
c  for each zone, : "star" particle is created if the density exceeds
c  some threshold and this is the highest level of refinement.  That's
c  it.
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
               if (d(i,j,k) .lt. densthresh) goto 10

c            make sure that we never give put than 90% of the cell's mass
c            into a star particle

               gasfrac = min( 0.9_RKIND, dt / timeconstant )

c            calculate star mass in solar masses.  If this is less than the
c            user-defined threshold mass, do NOT make a star in this cell.
c            This is not exactly in keeping with the spirit of the Kravtsov
c            algorithm, and is somewhat degenerate with the density threshold,
c            but we really don't want millions and millions of star particles.

               starmass = gasfrac*d(i,j,k)*dble(d1)*dble(x1*dx)**3
     &                   / msolar

               if(starmass .lt. smthresh) goto 10

c
c              If both of these criteria are met, create a star particle
c
               ii = ii + 1


c            fill in star particle attributes

               mp(ii)  = gasfrac * d(i,j,k)
               tcp(ii) = t
               tdp(ii) = timeconstant
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
               if (imetalSNIa .eq. 1) then
                  metalfSNIa(ii) = metalSNIa(i,j,k)
               endif

c              Remove mass from grid
c
               d(i,j,k) = (1._RKIND - gasfrac)*d(i,j,k)
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
         write(6,*) 'star_maker4: reached max new particle count'
         ERROR_MESSAGE
      endif
      np = ii

cc
      return
      end
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback4(nx, ny, nz,
     &                      d, dm, te, ge, u, v, w, metal,
     &                      idual, imetal, imethod, dt, r, dx, t, z,
     &                      d1, x1, v1, t1, sn_param, m_eject, yield,
     &                      npart, xstart, ystart, zstart, ibuff,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, type)
c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Brian O'Shea
c  date:       29 March 2005
c
c    This is a simplified version of stellar feedback, using the method
c    of Andrey Kravtsov as described in "On the origin of the global
c    Schmidt law of star formation", 2003 ApJ 590:L1-L4.  As in this paper,
c    the metals and energy are deposited instantaneously into the gas.
c

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
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    npart    - particle array size specified by calling routine
c    sn_param - fraction of stellar rest mass that goes to feedback
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
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
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, imethod
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      INTG_PREC type(npart)
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      INTG_PREC i, j, k, n
      R_PREC mform, tfactor, clight, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, dratio
      parameter (clight = c_light, msolar_e51 = 1800.0_RKIND)
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
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
c         skip if we're past the first timestep (only feed back
c         on the first timestep)
c

            if( (t-tcp(n)) .gt. dt ) goto 10

c
c          calculate fraction of mass from star particle put back
c          into gas in this cell
c

            mform = m_eject * mp(n)

c
c           Calculate how much of the star's energy should have
c           gone into supernovae energy:  Do this instantaneously.
c           use the star mass before ejected mass is removed
c
            energy = sn_param * mp(n) * (clight/v1)**2 / 
     &                 (d(i,j,k)+mform)

c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae, etc.)
c
            mp(n) = mp(n) - mform

c
c           Add energy to energy field
c
            dratio = d(i,j,k)/(d(i,j,k) + mform)
            te(i,j,k) = te(i,j,k)*dratio + energy
            if (idual .eq. 1) 
     &          ge(i,j,k) = ge(i,j,k)*dratio + energy

c
c           Metal feedback (note that in this function gas metal is
c             a fraction (rho_metal/rho_gas) rather than a density.
c             The conversion has been done in the handling routine)
c
            if (imetal .eq. 1) then
c
c  "Cen method".  This takes into account gas recycling.
c   Basically what this equation does is say:
c
c   new_metal_fraction = ((old metal in gas) + (newly processed metal ejected
c                          from star) + (metal already in star that's ejected))
c                          / ((old metal density)  + 
c                             (gas density added from star))
c
               metal(i,j,k) = (metal(i,j,k)*d(i,j,k) + mform * 
     &              (yield * (1._RKIND-metalf(n)) + metalf(n)))
     &              / (d(i,j,k)+mform)  ! metal is a fraction
c
            endif
c
c           Mass and momentum feedback
c
            u(i,j,k) = u(i,j,k)*d(i,j,k) + mform * up(n)
            v(i,j,k) = v(i,j,k)*d(i,j,k) + mform * vp(n)
            w(i,j,k) = w(i,j,k)*d(i,j,k) + mform * wp(n)
            d(i,j,k) = d(i,j,k) + mform 
            u(i,j,k) = u(i,j,k)/d(i,j,k)
            v(i,j,k) = v(i,j,k)/d(i,j,k)
            w(i,j,k) = w(i,j,k)/d(i,j,k)
c
c           If te is really total energy (and it is unless imethod=2),
c             then just set this value
c
            if (imethod .ne. 2 .and. idual .eq. 1) then
               te(i,j,k) = 0.5_RKIND*(u(i,j,k)**2 + v(i,j,k)**2 + 
     &                          w(i,j,k)**2) + ge(i,j,k)
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

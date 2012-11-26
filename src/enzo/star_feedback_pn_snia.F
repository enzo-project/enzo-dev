#include "fortran.def"
#include "phys_const.def"
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback_pn_snia(nx, ny, nz,
     &                      d, dm, te, ge, u, v, w, metal,
     &                      idual, imetal, imethod, dt, r, dx, t, z,
     &                      d1, x1, v1, t1, sn_param, m_eject, yield,
     &                      distrad, diststep, distcells,
     &                      npart, xstart, ystart, zstart, ibuff,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, type, justburn,
     &                      iPN, imetalSNIa, metalSNIa)

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
c  modified2: 20 Feb 2008 by M. Ryan Joung
c    added metalSNIa & metalfSNIa; included feedback from SN Ia/PN
c  modified3: 3 Apr 2008 by M. Ryan Joung
c    added smthresh to use as proxy for initial SP mass
c  modified4: 10 Feb 2009 by Shikui Tang
c    use the KENTARO NAGAMINE's stellar feedback fitting formula
c    http://www.physics.unlv.edu/~kn/SNIa_2/
c  modified5: 21 Aug 2011 by John Wise
c    Removed Type II SNe, so this can be used in conjunction with any
c    other feedback routine
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
c    distrad  - feedback distribution radius in cells
c    diststep - distance in walking steps to deposit feedback
c    distcells - total number of cells over which to distribute feedback
c
c    iPN        - planetary nebulae feedback (0 - none, 1 - yes)
c    imetalSNIa - SN Ia metallicity flag (0 - none, 1 - yes)
c    metalSNIa  - SN Ia metallicity density field
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
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, imethod,
     &      distrad, diststep, distcells, imetalSNIa, iPN
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, justburn
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      INTG_PREC type(npart)
      R_PREC    metalSNIa(nx,ny,nz)

c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      INTG_PREC i, j, k, n, ic, jc, kc, stepk, stepj, cellstep
      R_PREC mform, tfactor, clight, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, dratio,
     &     distmass, delta_z
      parameter (clight = c_light, msolar_e51 = 1800._RKIND)
c     TSK variables for stellar feedback fitting formula      
      R_PREC EdotSNII, EdotSNIa, EdotPN
      R_PREC MdotSNII, MdotSNIa, MdotPN
      R_PREC Edota, Edotb, Mdota, Mdotb, dvol
      R_PREC tau, rmf, dtt, madd, madd_cell, delta_e, tmid, tnow
      real*8 msolar, esnIa, eSNII
      parameter (tau=1.0e9_RKIND,  msolar=1.989D33) ! tau in years
      parameter (eSNIa=5.D-6, eSNII=5.D-6) 
c-----------------------------------------------------------------------
c
c     The star particle creation algorithm partnered with this feedback
c     algorithm creates a star particle instantaneously.  However, we do
c     feedback as if the star particles are created over a long period
c     of time (in code units), so the particle actually loses mass over
c     time in an exponentially decaying way.
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
c      write(6,*) 'star_feedback2: start'
c
c     Conversion from mass units to solar masses
c
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and.
     &        (type(n) .eq. 2 .or. type(n) .eq. 7)) then
c
         tnow = (t-tcp(n))*t1/3.15e7_RKIND ! tnow is in yr
         dtt = dt*t1/3.15e7_RKIND

c        !!TSK 2/10/09 recover the initial mass of each star particle based on
c          its current mass (ignore the mass loss due to Type Ia SNe);
c          determine the mass and energy injection during this step using
c          the intergral form to let the feedback calculation independent of
c          time steps; but there is no intergral form for the SN Ia feedback.
c        !!Pay attention to the units required by subroutines in Ken_formula.f
c
         call cal_rmf(tnow, rmf, Mdota, Mdotb)
         minitial = mp(n) / rmf
         call cal_rmf(tnow-dtt, rmf, MdotSNII, MdotPN)
c         MdotSNII = minitial * (Mdota - MdotSNII)

         if (iPN .eq. 1) then
            MdotPN = minitial * (Mdotb - MdotPN)
         else
            MdotPN = 0._RKIND
         endif

         tmid = max(tnow-0.5_RKIND*dtt, 0.5_RKIND*dtt)
         call SNIa_feedback(tmid, MdotSNIa, tau)
         if (imetalSNIa .eq. 1) then
            MdotSNIa = minitial * MdotSNIa * dtt
         else
            MdotSNIa = 0._RKIND
         endif
c
c        We define E_X = e_X * mp * c^2, X=SNIa, SNII
c        In the code it is convenient using E_X = e'_X * mloss_X * c^2
c        thus e'_X = mp/mloss_X * e_X
c        For the given Ken's mass loss rate formula, for an 13.6 Gyr SP,
c        mp/mloss_SNIa = 1./1.2343D-3
c        mp/mloss_SNII = 1./0.15
c        PN wind has a velocity of 10 km/sec
c        The unit of Energy input is ergs
c
c         EdotSNII = eSNII/0.150081 * MdotSNII*(clight/v1)**2
         if (imetalSNIa .eq. 1) then
            EdotSNIa = eSNIa/1.234e-3_RKIND * MdotSNIa*(clight/v1)**2
         else
            EdotSNIa = 0._RKIND
         endif
         if (iPN .eq. 1) then
            EdotPN = 0.5_RKIND * MdotPN * (1.e6_RKIND/v1)**2
         else
            EdotPN = 0._RKIND
         endif
c
c        Now update the mass of star particles and  get the energy
c        and mass changed per volume in code unit,
c        madd (mass density increment), delta_e (energy density increment)
c
         mform = MdotSNIa + MdotPN
         delta_e = EdotSNIa + EdotPN
c
c        Compute index of the cell that the star particle resides in.
c 
         i = int((xp(n) - xstart)/dx,IKIND) + 1
         j = int((yp(n) - ystart)/dx,IKIND) + 1
         k = int((zp(n) - zstart)/dx,IKIND) + 1
c
c         check bounds - if star particle is outside of this grid
c         then exit and give a warning.
c
         if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &        .or. k .lt. 1 .or. k .gt. nz) then
            write(6,*) 'warning: star particle out of grid',i,j,k
            goto 100
         endif
c
c        skip if very little mass is formed.
c
         if (mform/d(i,j,k) .lt. 1.e-10_RKIND) goto 10
c
c        calculate mass added to each cell
c
         distmass = m_eject * mform / distcells
c
c        if using distributed feedback, check if particle is too close
c        to the boundary
c
         if (distrad .gt. 0) then
            i = max((1 + ibuff + distrad), 
     &           min((nx - ibuff - distrad), i))
            j = max((1 + ibuff + distrad), 
     &           min((ny - ibuff - distrad), j))
            k = max((1 + ibuff + distrad), 
     &           min((nz - ibuff - distrad), k))
         endif
c
c        subtract ejected mass from particle (ejection due
c        to winds, supernovae)
c
         mp(n) = mp(n) - mform * m_eject
c
c        Record amount of star formation in this grid.
c
         justburn = justburn + mform * dt * dx**3
c
c        Calculate how much of the star formation in this
c        timestep would have gone into supernova energy.
c
         energy = delta_e / distcells

c         write(6,*) 'feed:', t, n, mform, mp(n), tnow, dtt, MdotSNIa,
c     &        energy, d(i,j,k), te(i,j,k), ge(i,j,k)
c
c        Add energy to energy field
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
     &                    energy) * dratio
                     if (idual .eq. 1)
     &                    ge(ic,jc,kc) = 
     &                    ((ge(ic,jc,kc)*d(ic,jc,kc)) + energy) * dratio
c
c        Metal feedback (note that in this function gas metal is a
c        fraction (rho_metal/rho_gas) rather than a density.  The
c        conversion has been done in the handling routine).  This only
c        applies the planetary nebulae feedback because the Type II SNe
c        were done in the original star_feedback routines.
c        !metal is a fraction
c
c        "Cen method".  This takes into account gas recycling.
c
                     if (iPN .eq. 1) then
                        metal(ic,jc,kc) = 
     &                       (metal(ic,jc,kc)*d(ic,jc,kc) + 
     &                       (MdotPN / distcells) * 
     &                       (m_eject * metalf(n))) * dratio
                     endif
c
c         Type Ia SNe ejecta are nearly all metals, i.e. yield = 1.0
c
                     if (imetalSNIa .eq. 1) then
                        delta_z = (mdotSNIa / distcells) * 
     &                       (1._RKIND * (1._RKIND-metalf(n)) + 
     &                       m_eject * metalf(n))
                        metalSNIa(ic,jc,kc) = 
     &                       (metalSNIa(ic,jc,kc)*d(ic,jc,kc) + 
     &                       delta_z) * dratio
                        metal(ic,jc,kc) = 
     &                       (metal(ic,jc,kc)*d(ic,jc,kc) +
     &                       delta_z) * dratio
                     endif
c
c           Mass and momentum feedback
c
                     u(ic,jc,kc) = u(ic,jc,kc)*d(ic,jc,kc) +
     &                    distmass * up(n)
                     v(ic,jc,kc) = v(ic,jc,kc)*d(ic,jc,kc) +
     &                    distmass * vp(n)
                     w(ic,jc,kc) = w(ic,jc,kc)*d(ic,jc,kc) +
     &                    distmass * wp(n)
                     d(ic,jc,kc) = d(ic,jc,kc) + distmass
                     u(ic,jc,kc) = u(ic,jc,kc)/d(ic,jc,kc)
                     v(ic,jc,kc) = v(ic,jc,kc)/d(ic,jc,kc)
                     w(ic,jc,kc) = w(ic,jc,kc)/d(ic,jc,kc)
c
c           If te is really total energy (and it is unless imethod=2),
c             then just set this value
c
                     if (imethod .ne. 2 .and. idual .eq. 1) then
                        te(ic,jc,kc) = 0.5_RKIND*(u(ic,jc,kc)**2 + 
     &                       v(ic,jc,kc)**2 + w(ic,jc,kc)**2) +
     &                       ge(ic,jc,kc)
                     endif
                  endif
               enddo
            enddo
         enddo
c
 10      continue
      endif
c
 100  continue
c
      enddo
c
c      write(6,*) 'star_feedback2: end'
      return
      end

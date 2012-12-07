#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c//////////////////////  SUBROUTINE CLUSTER_MAKER \\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cluster_maker(nx, ny, nz, formleft, formright,
     &                      d, dm, temp, u, v, w, cooltime,
     &                      dt, r, metal, dx, t, z, procnum, 
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, imulti, efficiency, 
     &                      metalcrit, odthresh, lifetime, level, np, 
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      type, ctype, justburn, iradtrans,
     &                      imetalSNIa, metalSNIa, metalfSNIa)
c
c  CREATES JEANS-RESOLVED STAR PARTICLES
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
c  modified5: 26 July 2007 by John Wise
c    modified pop3_maker to model stellar clusters using radiative feedback
c  modified6: 20 Feb 2008 by JHW
c    added metalSNIa & metalfSNIa; included feedback from SN Ia/PN
c    (original changes by M. Ryan Joung)
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
c    lifetime - lifetime of a stellar particle
c    h2crit   - critical value for primordial star formation
c    metalcrit- critical value for primordial star formation
c    imetalSNIa - SN Ia metallicity flag (0 - none, 1 - yes)
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
      INTG_PREC procnum, ctype, iradtrans, imulti, imetalSNIa
      R_PREC    formleft(3), formright(3)
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    metalSNIa(nx,ny,nz), metalfSNIa(nmax)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, justburn
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      INTG_PREC type(nmax)
      R_PREC    odthresh, masseff, starmass, min_tdyn, efficiency
      R_PREC    lifetime
      R_PREC    metalcrit
c
c  Locals:
c
      INTG_PREC  i, j, k, ii, irad, ib, jb, kb, maxradius, radius2
      INTG_PREC  dim
      R_PREC   div, tdyn, dtot, wgt
      R_PREC   pi, G, sndspdC, m1
      R_PREC   isosndsp2, bmass, jeanmass, lum, min_temp
      P_PREC this_pos(3)
      real*8 msolar, mh
      parameter (pi=pi_val,
     &           G=GravConst,
     &           sndspdC=1.3095e8_RKIND, 
     &           msolar=SolarMass, 
     &           mh=mass_h)
c
c      if (iradtrans .eq. 0) then
c         write(6,*) 'Radiative transfer required for cluster_maker!'
c         write(6,*) 'Check parameter file.'
c         stop
c      endif

c      write(6,*) 'cluster maker:'
c      write(6,*) odthresh

      ii = np
c
c     calculate how many solar masses are in d1*x1^3
      m1 = dble(x1*dx)**3 * dble(d1) / msolar
c
c     Different minimum temperature when to ignore tcool < tdyn
      if (imulti.le.1) min_temp = 1.e4_RKIND
      if (imulti.ge.2) min_temp = 1.e3_RKIND
c
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

               this_pos(1) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               this_pos(2) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               this_pos(3) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
c
c              1) finest level of refinement?
c
               if (r(i,j,k) .ne. 0.0) goto 10
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
               if (div .ge. 0._RKIND) goto 10
c
c              4) t_cool < t_free-fall (if T < 8000.0 skip this check)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*d1
               tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1
c               if (tdyn .gt. cooltime(i,j,k) .or.
c     &            temp(i,j,k) .lt. 1.1e4) 
c     &            write(6,*) 'C',i,j,k,tdyn,cooltime(i,j,k),temp(i,j,k)
c               if (tdyn .lt. cooltime(i,j,k)) goto 10

c     Lifetime > 0 (resolved SF. Don't use tdyn)
               if (temp(i,j,k) .gt. min_temp .and. 
     &              lifetime.gt.0) goto 10
c     Lifetime < 0 means unresolved SF.  Use tcool<tdyn.
               if (lifetime.lt.0 .and.
     &              tdyn .lt. cooltime(i,j,k) .and. 
     &              temp(i,j,k) .gt. min_temp) goto 10

c               if (tdyn .lt. cooltime(i,j,k) .and. 
c     &             temp(i,j,k) .gt. min_temp) goto 10
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
c               jeanmass = pi/(6.0*sqrt(d(i,j,k)*dble(d1))) *
c     &                    dble(pi * isosndsp2 / G)**1.5 / msolar
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
c              7) Check if metallicity is greater than critical value
c
               if (imetal .eq. 1 .and. metal(i,j,k) .lt. metalcrit)
     $              goto 10
c
c              8) Check if molecular hydrogren fraction is greater than
c              critical value
c
c               if (h2d(i,j,k) .lt. h2crit) goto 10
c
c              9) Check if cell is within the star forming region
c
               do dim=1,3
                  if (this_pos(dim) .lt. formleft(dim) .or.
     $                 this_pos(dim) .gt. formright(dim))
     $                 goto 10
               enddo
c
c              Create a star particle
c
               ii = ii + 1
c               mp(ii)   = efficiency*bmass
               mp(ii)  = efficiency*d(i,j,k)
               if (lifetime.lt.0._RKIND) then
                  tdp(ii) = 4._RKIND*tdyn
               else
                  tdp(ii) = lifetime*3.15e7_RKIND/t1 ! yr -> code
               endif

               type(ii) = -ctype
               tcp(ii) = t
               xp(ii) = this_pos(1)
               yp(ii) = this_pos(2)
               zp(ii) = this_pos(3)
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
               if (imetalSNIa .eq. 1) then
                  metalfSNIa(ii) = metalSNIa(i,j,k)
               endif
c
c              Remove mass from grid

               d(i,j,k) = (1._RKIND - efficiency) * d(i,j,k)

c               write(6,777) d(i,j,k), xp(ii), yp(ii), zp(ii)
c               write(6,7777) i, j, k
c               write(6,778) up(ii), vp(ii), wp(ii), mp(ii)
c               write(6,779) type(ii), tcp(ii), tdp(ii), metalf(ii)
 777           format('clst_maker: ', e10.3, 3(f12.6))
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

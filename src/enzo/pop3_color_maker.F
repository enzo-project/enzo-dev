#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c/////////////////////  SUBROUTINE POP3_COLOR_MAKER \\\\\\\\\\\\\\\\\\\\
c
      subroutine pop3_color_maker(nx, ny, nz,
     &                      d, dm, u, v, w, 
     &                      dt, r, dx, t, z, procnum, 
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imethod, odthresh, level, np,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tcp, type, ctype )
#ifndef CONFIG_PFLOAT_16
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
c  modified5: June, 2009 by Matthew Turk
c    modified pop3_maker for much simpler, overdensity star formation
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
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imethod
      INTG_PREC procnum, ctype
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      INTG_PREC type(nmax)
      R_PREC    odthresh
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
c               write(0,*) d(i,j,k)
               if (r(i,j,k) .ne. 0._RKIND) goto 10
c
c              2) density greater than threshold
c
               if (d(i,j,k)*d1 .lt. odthresh) goto 10
c
c              3) divergence negative
c                 (the first calculation is face centered for ZEUS, 
c                  the second is cell-centered for PPM)
c
c
c              Create a star particle
c
               ii = ii + 1

               type(ii) = -ctype
               tcp(ii) = t
               xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
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
               mp(ii) = tiny

c                write(6,777) d(i,j,k)*d1, xp(ii), yp(ii), zp(ii)
C                write(6,7777) i, j, k
C                write(6,778) up(ii), vp(ii), wp(ii), mp(ii)
C                write(6,779) type(ii), tcp(ii), tdp(ii)
 777           format('pop3_maker: ', e10.3, 3(f12.6))
 7777          format('            ', 3(i5))
 778           format('            ', 4(e15.3))
 779           format('            ', i5, 2(f15.6))
c
c               write(7+procnum,1000) bmass*starfraction,tdp(ii),tcp(ii),
c     &                       metalf(ii)
c               write(7+procnum,1000) level,bmass*starfraction,tcp(ii),
c     &                           tdp(ii)*t1,d(i,j,k)*d1,z,metalf(ii)

 1000          format(i5,1x,6(1pe10.3,1x))
c
c              Do not generate more star particles than available
c
               if (ii .ge. nmax) goto 20
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
#endif
      end

#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker5(nx, ny, nz,
     &                      d, dm, temp, coolrate, u, v, w,
     &                      dt, r, metal, dx, t, z, procnum, 
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, mintdyn,
     &                      odthresh, shdens,
     &                      smthresh, level, np,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      rr_left0, rr_left1, rr_left2,
     &                      rr_right0, rr_right1, rr_right2,
     &                      imetalSNIa, metalSNIa, metalfSNIa)

c
c  CREATES STAR PARTICLES
c
c  written by:  Sam Skillman & Brian O'Shea
c  date:  23 June 2006
c    This algorithm is taken from Springel & Hernquist 2003,
c	  MNRAS 339,289-311.
c    All Equation numbers listed below are from their paper. 
c
c  modified by: Robert harkness, August 12th 2006
c               Use Fortran90 random numbers
c  modified by: Stephen Skory, June 30, 2008
c               Fix units and signs
c  modified by: Samuel Skillman, Stephen Skory, January 17th, 2010
c               More unit fixes, added the shdens parameter
c
c  
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    temp  - temperature field
c    coolrate - cooling rate of gas in cgs [note - positive is heading, neg. is cooling]
c    u,v,w - velocity fields
c    r     - refinement field (non-zero if zone is further refined)
c    metal - metallicity field
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    odthresh - Overdensity of the cell required to turn on star formation.
c    shdesn   - Critical density (rho_crit) for star formation, see Springel & Hernquist
c    smthresh - mass of all stars that are created ~= 7.6e4 (msolar)
c    mintdyn  - timescale for star formation ~= 2.1e9 (yrs)
c    level - current level of refinement
c    procnum - processor number (for output)
c    rr_left, rr_right - refinement region left and right edges
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
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz)
      R_PREC    temp(nx,ny,nz), coolrate(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    metalSNIa(nx,ny,nz), metalfSNIa(nmax)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, msolar1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      R_PREC    odthresh, shdens, smthresh, mintdyn
      P_PREC rr_left0, rr_left1, rr_left2
      P_PREC rr_right0, rr_right1, rr_right2
c
      R_PREC   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      INTG_PREC  i, j, k, ii, ran1_init
      R_PREC   div, tdyn, dtot, dummy1
      R_PREC   newd, newdt, newcool, newtemp
      R_PREC   tempx, tempy, tempz
      character*15 fname
      R_PREC   usn, tstar, pstar, sqrtepssn, timeconstant
      R_PREC   mproton, RandomSeed, beta,x ,y, kb
      R_PREC   starfraction, bmass, densthresh 
      real*8 msolar
      parameter (msolar=SolarMass, mproton=mass_h, 
     &           beta=0.1_RKIND, sqrtepssn=2.e24, 
     &           RandomSeed=-123456, kb=kboltz)

      R_PREC random

      ii = 0
c     Initialize Random Number Generator Once
      if(ran1_init .eq. 0)then
        call random_seed
        ran1_init = 1
      endif
c      Uncomment the stuff below to write star formation stuff to disk.
c      write(fname, "(a6, i2.2, a4)") "smout", procnum, ".txt"
c      open (unit=1701,file=fname,access="APPEND",status="UNKNOWN")

c     Define the energy output from supernovae

      usn=(1._RKIND-beta)*sqrtepssn/beta/msolar ! Below Eq (2)
      usn=usn*sqrtepssn   ! this trick done to avoid floating-point issues (BWO)

c     Loop over all grids
      do k=1+ibuff,nz-ibuff   
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
c	       Only if this is the most refined cell at this location.
	       if(r(i,j,k) .ne. 0._RKIND) goto 10

c             Only over-dense enough cells.
              if(d(i,j,k) .le. odthresh) goto 10

c         Don't consider this cell if it's not in the refinement region.
              tempx = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
              tempy = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
              tempz = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
              if (tempx .lt. rr_left0 .or. tempx .gt. rr_right0 .or. 
     &          tempy .lt. rr_left1 .or. tempy .gt. rr_right1 .or. 
     &          tempz .lt. rr_left2 .or. tempz .gt. rr_right2) then
                    goto 10
               endif

c	       Calculate star formation timescale               
               tstar = 31556926._RKIND*mintdyn*(d(i,j,k)*dble(d1)/  
     &              shdens)**(-0.5_RKIND)			!Eq(21)
            
            
c     note: minus sign is because 'coolrate' is actual backwards: when
c     coolrate is negative, gas is cooling (which is the opposite of how
c     it's defined in springel & Hernquist 2003, MNRAS, 339, 289, eqn. 16
               y = -1._RKIND*tstar*coolrate(i,j,k)/(d(i,j,k)*dble(d1))/
     &              (beta*usn-(1._RKIND-beta)*1.5_RKIND*kb*temp(i,j,k)/
     &              0.6_RKIND/mproton) !Eq(16)

c	       Calculate the fraction of mass in cold clouds
               if (y .le. 0._RKIND) then
                  x = 0._RKIND
               else
                  x = 1._RKIND + 1._RKIND/2._RKIND/y - 
     &                 sqrt(1._RKIND/y + 1._RKIND/4._RKIND/y/y)	!Eq(18)
               endif

c              Calculate total baryon mass in the cell               
               bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3/msolar	     

c              Calculate a parameter which is compared to a random number
               pstar=(bmass/smthresh)*
     &              (1._RKIND-exp(-(1._RKIND-beta)*x*dt*dble(t1)/tstar)) !Eq(39)

               newd = d(i,j,k) * dble(d1)
               newdt = dt*dble(t1)
               newcool = coolrate(i,j,k)
               newtemp = temp(i,j,k)
c               Uncomment the below to write out star formation stuff to disk.
c               write(1701, *) bmass, newd, newcool, newtemp, newdt,pstar

!               if (pstar .gt. 100.) then
!                  write(6,*) 'the value of pstar that you are comparing
!     &             is large, consider raising star mass or timescale'
!               endif

               call random_number(dummy1)

c ----------If a random number < pstar, create a star!----------------

               if(dummy1 .gt. pstar) goto 10
               starfraction = min(smthresh/bmass, 0.5_RKIND)

               ii = ii + 1
	       mp(ii)  = starfraction * d(i,j,k) 
               tcp(ii) = t
               tdp(ii) = tdyn
               xp(ii) = tempx
               yp(ii) = tempy
               zp(ii) = tempz
               if (imethod .eq. 2) then
                  up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
                  vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
                  wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))
               else
                  up(ii) = u(i,j,k)
                  vp(ii) = v(i,j,k)
                  wp(ii) = w(i,j,k)
               endif

c	       Set metal Fraction
               if (imetal .eq. 1) then
                  metalf(ii) = metal(i,j,k)
               else
                  metalf(ii) = 0._RKIND
               endif
               if (imetalSNIa .eq. 1) then
                  metalfSNIa(ii) = metalSNIa(i,j,k)
               endif

c	       Remove mass from grid
               d(i,j,k) = (1._RKIND - starfraction)*d(i,j,k)
c
c              Do not generate more star particles than available
c
               if (ii .eq. nmax) goto 20

 10            continue
            enddo
         enddo
      enddo
 20   continue
c
c      Uncomment to write stuff to disk.
c      close(1701)

      if (ii .ge. nmax) then
         write(6,*) 'star_maker5: reached max new particle count'
         ERROR_MESSAGE
      endif

c
c     Make sure to return the number of particles!
c
      np = ii
c
      return
      end
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback5(nx, ny, nz,
     &     d, dm, te, ge, u, v, w, metal,
     &     idual, imetal, imethod, dt, r, dx, t, z,
     &     d1, x1, v1, t1, yield,
     &     npart, xstart, ystart, zstart, ibuff,
     &     xp, yp, zp, up, vp, wp,
     &     mp, tdp, tcp, metalf, type)
c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Sam Skillman & Brian O'Shea
c  date:       23 June 2006
c
c  This is a stellar feedback version from Springel & Hernquist.
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
     &     yield, minitial, xv1, xv2, dratio, sqrtepssn,
     &     msolar, beta
      parameter (clight = c_light, msolar_e51 = 1800._RKIND, 
     &     msolar=SolarMass, sqrtepssn=2.e24_RKIND, beta=0.1_RKIND)

c    Basically calculate some of the same stuff from above.
c    Dump it all into thermal making sure to take care of
c    energy formalisms.

      do n=1, npart
         if(tcp(n) .eq. t .and. mp(n) .gt. 0 .and. type(n) .eq. 2) then
            
            i = int((xp(n) - xstart)/dx,IKIND) +1
            j = int((yp(n) - ystart)/dx,IKIND) +1
            k = int((zp(n) - zstart)/dx,IKIND) +1

c
c         check bounds - if star particle is outside of this grid
c         then exit and give a warning.
c

            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &          .or. k .lt. 1 .or. k .gt. nz) then
               write(6,*) 'warning: star particle out of grid',i,j,k
               write(6,*) 'xp yp zp',xp(n),yp(n),zp(n)
               write(6,*) 'xstart ystart zstart',xstart,ystart,zstart
               goto 100
            endif


c	    Metal Feedback
            metal(i,j,k) = metal(i,j,k) +
     &            (1._RKIND-beta)*yield*mp(n)/d(i,j,k)		!Eq(40)

            energy = sqrtepssn*mp(n)/msolar/v1/v1/d(i,j,k)  
            energy = energy*sqrtepssn    ! this is due to floating-point issues (BWO)

            te(i,j,k) = te(i,j,k) + energy
            
            if(idual .eq. 1)
     &           ge(i,j,k) = ge(i,j,k) + energy
            
            if (imethod .ne. 2 .and. idual .eq. 1) then
               te(i,j,k) = 0.5_RKIND*(u(i,j,k)**2 + v(i,j,k)**2 +
     &              w(i,j,k)**2) + ge(i,j,k)
            endif
            
 100        continue            
         endif
      enddo

      return
      end

#include "fortran.def"
#include "phys_const.def"
#include "error.def"
c
c=======================================================================
c/////////////////  SUBROUTINE PARTICLE_SPLITTER  \\\\\\\\\\\\\\\\\\\\\\
c
      subroutine particle_splitter(nx, ny, nz,
     &                      idual, imetal, imethod, dt, r, dx, t, z,
     &                      d1, x1, v1, t1,
     &                      xstart, ystart, zstart, ibuff,
     &                      npart,
     &                      xpold, ypold, zpold, upold, vpold, wpold,
     &                      mpold, tdpold, tcpold, metalfold, typeold,
     &                      nmax, npartnew, children, level,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, type, 
     &                      iterations, separation, ran1_init,
     &                      rr_leftedge, rr_rightedge)
c
c  SPLIT PARTICLES INTO CHILDREN PARTICLES WHEN REQUESTED
c
c  written by: Ji-hoon Kim
c  date:       15 October 2009
c 
c  description: Adopted from star_maker and star_feedback routines
c               with many more additional parameters.  Keep in mind 
c               that below routines are not consistent with 
c               children != 12.  See Grid_ParticleSplitter.C or 
c               ParticleSplitter.C for physical/technical details.
c 		
c  modified1:  
c
c
c  PARAMETERS:
c
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
c
c    npart    - particle array size specified by calling routine
c    nmax     - the maximum number of particles to be created
c    npartnew - particle array size that will be created here, initially 0
c    children - number of children particle per parent = 12
c    level    - current level of refinement
c    iterations - ParticleSplitterIterations (> 1 not recommended)
c    separation - ParticleSplitterChildrenParticleSeparation (= 1 by default)
c    rr_leftedge - RefineRegionLeftEdge[]
c    rr_rightedge - RefineRegionRightEdge[]
c
c    xp,yp,zp - positions of particles created
c    up,vp,wp - velocities of particles created
c    mp       - mass of particles created
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle created (-1 if not a star particle)
c    metalf   - metal fraction of particle created 
c    *old     - array of particle properties that already exist
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart
      INTG_PREC  idual, imetal, imethod, iterations
      R_PREC    dt, dx, z, separation
      R_PREC    d1, x1, v1, t1
      R_PREC    r(nx,ny,nz)
      P_PREC xstart, ystart, zstart, t
      P_PREC xpold(npart), ypold(npart), zpold(npart)
      R_PREC    upold(npart), vpold(npart), wpold(npart)
      R_PREC    mpold(npart), tdpold(npart), tcpold(npart)
      R_PREC    metalfold(npart)
      INTG_PREC typeold(npart)
      INTG_PREC npartnew, nmax, children, level, ran1_init
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      INTG_PREC type(nmax)
      P_PREC rr_leftedge(3), rr_rightedge(3)
c
c  Locals
c
      INTG_PREC n, m, ii, i, j, k
      R_PREC    xp0(children), yp0(children), zp0(children)
      R_PREC    rad, pi
      R_PREC    alpha(3)
      R_PREC    l11, l12, l13
      R_PREC    l21, l22, l23
      R_PREC    l31, l32, l33      
      parameter (pi=pi_val)
c
c-----------------------------------------------------------------------
c
c      write(6,*) 'particle_splitter: start'
c
c     Initialize random number generator once
c
      if (ran1_init .eq. 0) then
         call random_seed()
         ran1_init = 1
      endif
c
c     The distance to the children particles rad is the same for all. It is 
c     currently set to be dx (=CellWidth), but can be changed. 
c     (i.e. ParticleSplitterChildrenParticleSeparation = 1 by default)
c
      rad = dx * separation**REAL(iterations,RKIND)
c
c
c     Loop over existing (parent) particles; It implicitly assumes that 
c     only DM and conventional star particles  get splitted.  Other 
c     particles - which usually become Star class particles - doesn't 
c     seem to have any reason to be splitted.  (as of Oct.2009)
c
      ii = 0
c
      do n = 1, npart
         if (mpold(n) .gt. 0 .and. typeold(n) .le. 2) then    
c
c           Check bounds [1] - if the particle is outside of RefineRegion
c           then give a warning and exit.
c
            if (xpold(n) .lt. rr_leftedge(1) .or. 
     &          xpold(n) .gt. rr_rightedge(1) .or.
     &          ypold(n) .lt. rr_leftedge(2) .or. 
     &          ypold(n) .gt. rr_rightedge(2) .or.
     &          zpold(n) .lt. rr_leftedge(3) .or. 
     &          zpold(n) .gt. rr_rightedge(3)) then
c               write(6,*) 'part_spl: parent particle outside RR :',
c     &                     xpold(n), ypold(n), zpold(n), 
c     &                     rr_leftedge(1), rr_leftedge(2), level
               goto 100
            endif
c
c           Compute index of the cell that the parent particle resides in.
c 
            i = int((xpold(n) - xstart)/dx,IKIND) + 1
            j = int((ypold(n) - ystart)/dx,IKIND) + 1
            k = int((zpold(n) - zstart)/dx,IKIND) + 1
c
c           Check bounds [2] - if the particle is outside of this grid
c           then give a warning and exit.  (This is very possible when 
c           particles are not properly distributed to the right grids)
c
            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. 
     &          j .gt. ny .or. k .lt. 1 .or. k .gt. nz) then
c               write(6,*) 'part_spl: parent particle out of grid :',
c     &                     i,j,k,nx,ny,nz, 
c     &                     xpold(n),ypold(n),zpold(n),
c     &                     xstart, ystart, zstart, dx, level
               goto 100
            endif
c
c           Check bounds [3] - is this finest level of refinement?
c           if the particles are not yet inside the finest grids 
c           then give a warning and exit
c
c            if (r(i,j,k) .ne. 0.0) then
c               write(6,*) 'part_spl: parent particle not at finest :',
c     &                     i,j,k,nx,ny,nz,level
c               goto 100
c            endif
c
c
c
cc               Below was a quick test
cc
cc               mpold(n) = mpold(n) * 0.5
cc
cc               ii = ii + 1
cc               xp(ii) = xpold(n) + 0.1*dx
cc               yp(ii) = ypold(n) 
cc               zp(ii) = zpold(n) 
cc
cc               mp(ii)  = mpold(n)
cc               tcp(ii) = tcpold(n)
cc               tdp(ii) = tdpold(n)
cc               up(ii) = upold(n)
cc               vp(ii) = vpold(n)
cc               wp(ii) = wpold(n)
cc               type(ii) = typeold(n)
cc
cc               if (imetal .eq. 1) then
cc                  metalf(ii) = metalfold(n)    
cc               else
cc                  metalf(ii) = 0.0
cc               endif
c
c
c
c
c
c           ================================================================
c                              CREATE CHILDREN PARTICLES
c           ================================================================
c
c
c           First reduce the mass of the parent down to 1/(children+parent) = 1/13
c 
c            mpold(n) = mpold(n) / 13._RKIND
            mpold(n) = mpold(n) / (children + 1._RKIND)
c
c
c           ===================================
c                  Step I - Positioning
c           ===================================
c  
c           Now carefully set the position.  Populate 12 children particles 
c           around the parent at (0,0,0) assuming HCP (hexagonal closed-packed) 
c           structure in an xyz space.  Keep in mind that Any three neighboring 
c           children and the parent will form a tetrahedron of edge length rad.  
c
c
c           Step[I-1]: 6 children on the same x-y plane as the parent 
c                      (counter-clockwise; again, here we assume the parent is at the origin)
c
            xp0(1) = rad
            yp0(1) = 0._RKIND 
            zp0(1) = 0._RKIND 
c
            xp0(2) = 0.5_RKIND*rad
            yp0(2) = 0.866_RKIND*rad
            zp0(2) = 0._RKIND
c
            xp0(3) = -0.5_RKIND*rad
            yp0(3) = 0.866_RKIND*rad
            zp0(3) = 0._RKIND
c
            xp0(4) = -rad
            yp0(4) = 0._RKIND
            zp0(4) = 0._RKIND
c
            xp0(5) = -0.5_RKIND*rad
            yp0(5) = -0.866_RKIND*rad
            zp0(5) = 0._RKIND
c
            xp0(6) = 0.5_RKIND*rad
            yp0(6) = -0.866_RKIND*rad
            zp0(6) = 0._RKIND
c
c           Step[I-2]: 3 children above the parent plane (x-y)
c                      (the height of the tetrahedron of edge r = sqrt(2/3)*r)
c
            xp0(7) = 0.5_RKIND*rad
            yp0(7) = 0.2887_RKIND*rad
            zp0(7) = 0.8165_RKIND*rad
c
            xp0(8) = -0.5_RKIND*rad
            yp0(8) = 0.2887_RKIND*rad
            zp0(8) = 0.8165*rad
c
            xp0(9) = 0._RKIND
            yp0(9) = -0.5774_RKIND*rad
            zp0(9) = +0.8165_RKIND*rad
c
c           Step[I-3]: 3 children below the parent plane (x-y)
c
            xp0(10) = 0.5_RKIND*rad
            yp0(10) = 0.2887_RKIND*rad
            zp0(10) = -0.8165_RKIND*rad
c
            xp0(11) = -0.5_RKIND*rad
            yp0(11) = 0.2887_RKIND*rad
            zp0(11) = -0.8165_RKIND*rad
c
            xp0(12) = 0._RKIND
            yp0(12) = -0.5774_RKIND*rad
            zp0(12) = -0.8165_RKIND*rad
c
c
c           ===================================
c                    Step II - Rotation
c           ===================================
c
c           Provide the particle with Euler rotation (z->x->z) to an arbitrary 
c           orientation.  See Eq.(11.99) of Marion & Thornton (1995)
c           -----------------------------------------------------------------   
c                         lambda_11 =  cosCcosA - cosBsinAsinC
c                         lambda_21 = -sinCcosA - cosBsinAcosC
c                         lambda_31 =  sinBsinA
c                         lambda_12 =  cosCsinA + cosBcosAsinC
c                         lambda_22 = -sinCsinA + cosBcosAcosC
c                         lambda_32 = -sinBcosA
c                         lambda_13 =  sinCsinB
c                         lambda_23 =  cosCsinB
c                         lambda_33 =  cosB
c           -----------------------------------------------------------------
c
c
c           Step[II-1]: First generate random numbers which will be in the range of [0,1].
c
            call random_number(alpha)
            do m = 1, 3
               alpha(m) = alpha(m) * pi * 2._RKIND
            enddo
c
c            write(6,*) 'random number:', alpha(1), alpha(2), alpha(3)
c            write(6,*) 'cos:', cos(alpha(1)), cos(alpha(2)), cos(alpha(3))
c
c
c           Step[II-2]: Calculate three Euler angles (alpha's).
c
            l11 = cos(alpha(3))*cos(alpha(1))
     &          - cos(alpha(2))*sin(alpha(1))*sin(alpha(3))
            l21 = -sin(alpha(3))*cos(alpha(1))
     &          - cos(alpha(2))*sin(alpha(1))*cos(alpha(3))
            l31 = sin(alpha(2))*sin(alpha(1))
c
            l12 = cos(alpha(3))*sin(alpha(1))
     &          + cos(alpha(2))*cos(alpha(1))*sin(alpha(3))
            l22 = -sin(alpha(3))*sin(alpha(1))
     &          + cos(alpha(2))*cos(alpha(1))*cos(alpha(3))
            l32 = -sin(alpha(2))*cos(alpha(1))
c
            l13 = sin(alpha(3))*sin(alpha(2))
            l23 = cos(alpha(3))*sin(alpha(2))
            l33 = cos(alpha(2))
c
c
c           Step[II-3]: Use [xyz]p0 and Euler angles to generate [xyz]p.
c                       Here we add the initial position of the parent [xyz]pold.
c            
            do m = 1, children
               xp(ii+m) = xpold(n) 
     &                  + l11*xp0(m) + l12*yp0(m) + l13*zp0(m)
               yp(ii+m) = ypold(n) 
     &                  + l21*xp0(m) + l22*yp0(m) + l23*zp0(m)
               zp(ii+m) = zpold(n) 
     &                  + l31*xp0(m) + l32*yp0(m) + l33*zp0(m)
            enddo
c
c
c           ===================================
c               Step III - Other Properties
c           ===================================
c
c           Lastly, assign the velocities and all other properties of the parent 
c           to all children particles.
c
            do m = 1, children
               mp(ii+m)  = mpold(n)
               tcp(ii+m) = tcpold(n)
               tdp(ii+m) = tdpold(n)
               up(ii+m) = upold(n)
               vp(ii+m) = vpold(n)
               wp(ii+m) = wpold(n)
               type(ii+m) = typeold(n)
c
               if (imetal .eq. 1) then
                  metalf(ii+m) = metalfold(n)    !in here metal is a fraction
               else
                  metalf(ii+m) = 0._RKIND
               endif
            enddo
c
c
c           Increase the number of new particles
c
            ii = ii + children
c            
c           Do not generate more particles than available
c
            if (ii .gt. nmax) goto 20
c
         endif
c
 100     continue
c
      enddo 
c
c
 20   if (ii .gt. nmax) then
         write(6,*) 'particle_splitter reached max new particle count'
         ERROR_MESSAGE
      endif
c
c
      npartnew = ii
c
c
      if (npartnew .ne. 0) then
         write(6,*) 'particles created in splitter:', npartnew, level
      endif
c
c     write(6,*) 'particle_splitter: end'
      return
      end

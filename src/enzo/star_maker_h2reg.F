#include "fortran.def"
#include "error.def"

!     set this if you want to log star particle info in a file (one per
!     processor)
#define NO_STAR_LOG

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker_h2reg(nx, ny, nz,
     &     d, fHI, temp, u, v, w, cooltime,
     &     dt, r, metal, dx, t, z, procnum, 
     &     d1, x1, v1, t1,
     &     nmax, xstart, ystart, zstart, ibuff, 
     &     imetal, imethod, 
     &     StarFormationEfficiency,
     &     StarFormationNumberDensityThreshold, 
     &     StarFormationMinimumMass, 
     &     MinimumH2FractionForStarFormation,
     &     StochasticStarFormation,
     &     UseSobolevColumn,
     &     SigmaOverR,
     &     AssumeColdWarmPressureBalance,
     &     H2DissociationFlux_MW, 
     &     H2FloorInColdGas,
     &     ColdGasTemperature,
     &     ran1_init, 
     &     level, grid_id,
     &     oldnp, oldxp, oldyp, oldzp, oldmp, oldtcp, oldtype,
     &     np, xp, yp, zp, up, vp, wp,
     &     mp, tdp, tcp, metalf)
c     
c  CREATES STAR PARTICLES
c
c  written by:  Michael Kuhlen
c  date:  November 2010
c
c    Molecular hydrogen regulated star formation.
c
c    This SF recipe is modeled after Krumholz & Tan (2007) and McKee &
c    Krumholz (2010). The SF time scale is the gas free fall time, and
c    thus the SFR density is proportional to the gas density to the
c    power of 1.5. 
c
c    d(M_*)/dt = eps * M_H2(>n) / t_ff(n)
c
c    The SF is proportional to the molecular hydrogen density, not the
c    total gas density. The H_2 fraction is estimated using the
c    prescription given by McKee & Krumholz (2010), which is based on 1D
c    radiative transfer calculations and depends on the neutral hydrogen
c    number density, the metallicity, and the H2 dissociating flux. If
c    AssumeColdWarmPressureBalance==1, then the additional assumption of
c    pressure balance between the Cold Neutral Medium and the Warm
c    Neutral Medium removes the dependence on the H2 dissociating flux
c    (Krumholz, McKee, & Tumlinson 2009).
c
c    Optionally, a proper number density threshold and/or an H2 fraction
c    threshold is applied, below which no star formation occurs.
c
c    Typically I use this with StarFormationOncePerRootGridTimeStep, in
c    which case this routine is called only once per root grid step and
c    only for grids on MaximumRefinementLevel, but with dt =
c    TopGridTimeStep (i.e. the root grid time step).
c
c    [I used star_maker4.F as a template for this routine.]
c
c  INPUTS:
c
c    d     - density field
c    fHI   - neutral hydrogen mass fraction (rho_HI / rho_tot)
c    temp  - temperature field
c    u,v,w - velocity fields
c    cooltime - cooling time in code units
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep  (root grid time step, if StarFormationOncePerRootGridTimeStep)
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c
c     StarFormationEfficiency -
c          Specific star formation efficiency per free-fall time
c          (typically ~0.01).
c
c     StarFormationNumberDensityThreshold -
c          In proper particle / cm^3. 
c
c     StarFormationMinimumMass - 
c          Star particles with mass less than this are NOT formed.
c
c     MinimumH2FractionForStarFormation -
c          Only allow star formation if the f_H2 is greater than this
c          value. This can be quite low, perhaps even zero (default
c          1e-5).
c     
c     StochasticStarFormation -
c          If the stellar mass is less than StarFormationMinimumMass,
c          then stochastically form a star particle of mass equal to
c          StarFormationMinimumMass with a probability of (stellar
c          mass)/StarFormationMinimumMass.
c
c     UseSobolevColumn - 
c          Use a Sobolev-like estimate of the neutral hydrogen column density
c          instead of just (rho_H * dx). For details look at the function
c          definition of SobolevColumn() at the end of this file.
c
c     SigmaOverR - 
c          sigma_{d,-21} / R_{-16.5}: the ratio of the dust cross
c          section per H nucleus to 1000 Angstroem radiation normalized
c          to 10^{-21} cm^{-2} (sigma_{d,-21}) to the rate coefficient
c          for H2 formation on dust grains normalized to the Milky Way
c          value of 10^{-16.5} cm^3 s^{-1}. Both are linearly
c          proportional to the dust to gas ratio and hence the ratio is
c          independent of metallicity. Although SigmaOverR is probably
c          close to unity in nature (see discussion in Krumholz, McKee,
c          & Tumlinson 2009), Krumholz & Gnedin argue that in
c          simulations with spatial resolution of ~50 pc, the value of R
c          should be increased by a factor of 30 in order to account for
c          the subgrid clumping of the gas.
c
c     AssumeColdWarmPressureBalance -
c          With the assumption of pressure balance between the cold
c          neutral medium and warm neutral medium, the dependence on the
c          dissociating flux drops out. (Krumholz & McKee 2010).
c
c     H2DissociationFlux_MW - 
c          The Lyman-Werner intensity in units of the Milky Way's
c          value. Only used, if AssumeColdWarmPressureBalance == FALSE.
c
c     H2FloorInColdGas - 
c          When this is set to >0.0, then a minimum f_H2 is imposed for
c          all cells that have a temperature less than
c          ColdGasTemperature.
c
c     ColdGasTemperature - 
c          Only gas with temperature below this parameter is eligible
c          for H2FloorInColdGas (see above).
c
c    level - current level of refinement
c    procnum - processor number (for output)
c    oldnp - number of particles in previous timestep (= NumberOfParticles)
c    old{x,y,z}p - positions of particles in previous timestep (= ParticlePosition[])
c    oldmp - mass of particles in previous timestep (= ParticleMass)
c    oldtcp - creation time of particles in previous timestep (= ParticleAttribute[0])
c    oldtype - type of particles in previous timestep (= ParticleType)
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
c
c
c     ******* More details on this star formation recipe. *******
c     

c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
#include "phys_const.def"
#define XH 0.76_RKIND
#define MU 1.22_RKIND

c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      INTG_PREC grid_id
      INTG_PREC procnum
      INTG_PREC ran1_init
      R_PREC    d(nx,ny,nz), fHI(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t

      INTG_PREC oldnp
      P_PREC oldxp(oldnp), oldyp(oldnp), oldzp(oldnp)
      R_PREC oldmp(oldnp), oldtcp(oldnp)
      INTG_PREC oldtype(oldnp)

      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)

      INTG_PREC StochasticStarFormation
      INTG_PREC UseSobolevColumn
      INTG_PREC AssumeColdWarmPressureBalance
      R_PREC StarFormationEfficiency
      R_PREC StarFormationNumberDensityThreshold
      R_PREC StarFormationMinimumMass
      R_PREC MinimumH2FractionForStarFormation
      R_PREC SigmaOverR
      R_PREC H2DissociationFlux_MW
      R_PREC H2FloorInColdGas
      R_PREC ColdGasTemperature

c
c  Functions:
c
      R_PREC SobolevColumn

c
c  Locals:
c
      INTG_PREC  i, j, k, n, ii
      R_PREC   pi, G, mproton, m1
      R_PREC   starmass, starmass_in_Msun
      R_PREC   densthresh, timeconstant, gasfrac
      R_PREC   tau_ff, tau_max
      R_PREC   random, critval

      R_PREC H2_fraction, Z_MW, nH
      R_PREC phi_CNM, chi, Sigma, tau_c, s
      R_PREC G0prime, t_avg, age
      R_PREC Sigma_cell

      R_PREC U_MW, D_MW, Dstar, alpha, gg, Lambda, x

      real*8 msolar
      parameter (pi=pi_val, G=6.67e-8_RKIND,
     $     msolar=1.989e33_RKIND,mproton=1.67e-24_RKIND)
      

#ifdef STAR_LOG      
      character*32 starlog_filename

#ifdef LARGE_TASK_COUNT
      character*8   post
#else
      character*4   post
#endif

! taken from mpi_error_file.F
#ifdef LARGE_TASK_COUNT
      write(post,'(i8)') procnum
      do i=1,8
      if(post(i:i).eq.' ') post(i:i)='0'
      end do
!     write(*,'(i8,4x,a8)') id,post
#else
      write(post,'(i4)') procnum
      do i=1,4
      if(post(i:i).eq.' ') post(i:i)='0'
      end do
!     write(*,'(i4,4x,a4)') id,post
#endif

      starlog_filename = 'starlog_' // post // '.txt'
c      print *,'starlog_filename = ',starlog_filename
      open(unit=4,file=starlog_filename,status='unknown',
     $     position='append')
#endif


      if (StochasticStarFormation .eq. 1) then 
c     Initialize random number generator once using the MPI rank
         if (ran1_init .eq. 0) then
            call init_random_seed(procnum)
            ran1_init = 1
         endif
      endif
c      

c     Calculate mass conversion factor
      m1 = dble(x1*dx)**3 * dble(d1) / msolar

      ii = np

c     Calculate density threshold in code units.
c
c     StarFormationNumberDensityThreshold is in proper particles
c     per cc and (d1 / (MU*mproton)) gives the mean number density of
c     the universe in particles/cm^3 (assuming a mean molecular weight
c     of 1.22, corresponding to neutral monatomic gas with X=0.7,
c     Y=0.24)
c
      densthresh = StarFormationNumberDensityThreshold / 
     $     ( d1 / (MU * mproton) )
      




c     loop over all cells
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
               
c     Note: The following check for whether cells are really at the
c     finest refinement level is not really necessary when
c     StarFormationOncePerRootGridTimeStep is set (which is how this
c     star_maker was initially designed to be run), since in that case
c     star formation is only turned on for grids on
c     MaximumRefinementLevel. Nevertheless I'm keeping this condition
c     around, in case someone wants to run H2-regulated SF without
c     StarFormationOncePerRootGridTimeStep.
c     
c     Is this finest level of refinement?
c     
               if (r(i,j,k) .ne. 0._RKIND) goto 10


c     Some conditions for SF to occur.

c     Is density greater than threshold?     
               if (d(i,j,k) .lt. densthresh) goto 10
                  
c               write(4,*) 'Found grid cell that is a star ', 
c     $              'particle candidate: ', d(i,j,k), densthresh
                  

c     Calculate the free fall time, which is the SF time scale

c     t_ff = [ 3 * pi / (32 * G * rho) ]^(-1/2)
c          = 6.66e6 yr * [ rho / 10^-22 g/cm^3 ]^(-1/2)
c          = 2.10e14 seconds * [ rho / 10^-22 g/cm^3 ]^(-1/2)

               timeconstant = 2.10e14_RKIND / t1 /
     $              sqrt( d(i,j,k)*dble(d1) / 1.0e-22_RKIND )

c     Calculate molecular hydrogen frection (if H2StarFormation)
               H2_fraction = 1.0

c     Metallicity normalized to Milky Way (solar neighborhood)
               Z_MW = metal(i,j,k) / 0.0204_RKIND ! I got this value from yt (source?)

c     [Actually, 0.024 is the solar metallicity, which is not
c     necessarily equal to the solar neighborhood metallicity. For
c     example, the stellar metallicity in the solar neighborhood
c     (measured in G dwarfs) is 0.63 * Z_sun. However, we want the
c     *gas-phase* solar neighborhood metallicity, and that is measured
c     to be consistent with Z_sun (Rodriguez & Delgado-Inglada
c     2011). Thanks to Nick Gnedin & Andrey Kravtsov for clearing this
c     up for me.]

c     Hydrogen number density in proper cm^-3
               nH = d(i,j,k) * fHI(i,j,k) * d1 / mproton

               if (AssumeColdWarmPressureBalance .eq. 1) then

c     Eq.(7) Krumholz, McKee, & Tumlinson (2009)
                  phi_CNM = 3._RKIND
                  chi = 2.3_RKIND * SigmaOverR / phi_CNM * 
     $                 ( 1._RKIND + 3.1_RKIND * Z_MW**0.365_RKIND )
               else

                  G0prime = H2DissociationFlux_MW
                     
c     Eq.(9) McKee & Krumholz (2010)
                  chi = 71._RKIND * SigmaOverR * G0prime / nH
               
               endif

c     Surface density of cell in Msun pc^-2 (1 g/cm^2 = 4788.03 Msun/pc^2)

               if (UseSobolevColumn .ne. 1) then
c              *** Grid cell based ***
c                  [ Sigma = rho * dx ]

c     total density column (Sigma_gas) [deprecated]
ccc   Sigma = d(i,j,k) * d1 * (x1*dx) * 4788.03

c     neutral hydrogen column (Sigma_HI)
                  Sigma = d(i,j,k) * fHI(i,j,k) * 
     $                 d1 * (x1*dx) * 4788.03_RKIND
               else
c              *** Sobolev-like approximation ***
c                  [ Sigma ~ rho * ( rho / grad(rho) ) ]
                  
                  Sigma = SobolevColumn(d, fHI, nx, ny, nz,
     $                 i, j, k, dx, d1, x1)
                  
                  Sigma_cell = d(i,j,k) * fHI(i,j,k) * 
     $                 d1 * (x1*dx) * 4788.03_RKIND
               endif
                     
c     Eq.(22) Krumholz, McKee, & Tumlinson (2009)
               tau_c = 0.067_RKIND * Z_MW * Sigma
               
c     Eq.(91) McKee & Krumholz (2010)
               s = log( 1._RKIND + 0.6_RKIND * chi + 0.01_RKIND*chi*chi) 
     $              / (0.6_RKIND * tau_c)
               
c     Eq.(93) McKee & Krumholz (2010)
               H2_fraction = 0._RKIND
               if (s .lt. 2._RKIND) then
                  H2_fraction = 1._RKIND - 0.75_RKIND * 
     $                 ( s / (1._RKIND + 0.25_RKIND*s) )
               endif
               
c               if (H2_fraction .eq. 0._RKIND) then
c                  write(4,*) 'f_H2=0: ',s,chi,tau_c,Z_MW,Sigma,
c     $                 Sigma_cell,nH,d(i,j,k), fHI(i,j,k)
c               endif
c               if (H2_fraction .eq. 1._RKIND) then
c                  print *, 'f_H2=1: ',s,chi,tau_c,Z_MW,Sigma
c               endif


c     If H2FloorInColdGas is non-zero, then enforce a floor in f_H2, but
c     only for cold gas (ColdGasTemperature).
               if (H2FloorInColdGas .gt. 0._RKIND) then
                  
                  if(H2_fraction .lt. H2FloorInColdGas) then
                     
c#ifdef STAR_LOG
c                     write(4,*) 'Found cell with fH2 = ',
c     $                    H2_fraction, ' < ',
c     $                    H2FloorInColdGas, ' T = ',
c     $                    temp(i,j,k)
c#endif


                     if(temp(i,j,k) .lt. ColdGasTemperature) then
c#ifdef STAR_LOG
c                        write(4,*) 'Enforced fH2 = 0.01 for cell ',
c     $                       'with T = ',temp(i,j,k),
c     $                       'and n = ',
c     $                       d(i,j,k)*( d1 / (MU * mproton))
c#endif

                        H2_fraction = H2FloorInColdGas
                     endif
                     
                  endif  ! if(H2_fraction .lt. H2FloorInColdGas)
                  
               endif  ! if (H2FloorInColdGas .gt. 0._RKIND)

               
c     No star formation if H2_fraction < MinimumH2FractionForStarFormation.
               if (H2_fraction .lt.
     $              MinimumH2FractionForStarFormation) then
c                  write(4,*) 'Star formation prevented because ',
c     $                 'H2 fraction too low: f_H2 = ',
c     $                 MinimumH2FractionForStarFormation
                  goto 10
               endif
               
c     Ensure that no more than 90% of the cell's mass is converted into
c     a star particle.
               gasfrac = min( 0.9_RKIND, 
     $              StarFormationEfficiency * H2_fraction *
     $              dt / timeconstant )
               
               
c     Calculate star mass (this is actually a density, since Enzo treats
c     particle masses as densities)
               starmass =  gasfrac * d(i,j,k) 
               
c     Calculate star mass in solar masses, which means multiplying
c     'starmass' by the cell volume to get an actual mass, and then
c     dividing by Msun. Only make a star particle if (a) this is greater
c     than StarFormationMinimumMass or (b) StochasticStarFormation is
c     ste and a random number is less than
c     starmass_in_Msun/StarFormationMinimumMass.
               
               starmass_in_Msun = starmass * m1
               
               if(starmass_in_Msun .lt. 
     $              StarFormationMinimumMass) then

#ifdef STAR_LOG
                  write(4,*) 'Star particle does not meet minimum ',
     $                 'mass threshold: ', starmass_in_Msun
c     $                 dt, timeconstant, 
c     $                 gasfrac, starmass, StarFormationMinimumMass
#endif
                     
                  if (StochasticStarFormation .eq. 1) then
                     call random_number(random)
                     critval = starmass_in_Msun / 
     $                    StarFormationMinimumMass
                     if (random .gt. critval ) then
                        goto 10
                     endif
                        
                     starmass_in_Msun = StarFormationMinimumMass
                     starmass = starmass_in_Msun / m1
#ifdef STAR_LOG
                     write(4,200) 'Star particle does not meet ',
     $                    'minimum mass threshold: M =', 
     $                    starmass_in_Msun,
     $                    ', but stochastic creation occurred.',
     $                    random, critval
                     
 200                 format(a,a,1pe10.3,a,2(1pe10.3))
#endif
                  endif
               endif



#ifdef STAR_LOG
               write(4,*) 'Star particle formed with M/Msun = ',
     $              starmass_in_Msun, starmass
#endif

               ii = ii + 1


c            fill in star particle attributes

               if (starmass .eq. 0) then
                  print *, 'starmass == 0!'
                  ERROR_MESSAGE
               endif
               
               mp(ii)  = starmass
               tcp(ii) = t
               tdp(ii) = timeconstant
               xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx

#define AVERAGE_PARTICLE_VELOCITY
#ifdef AVERAGE_PARTICLE_VELOCITY
               if (imethod .eq. 2) then
                  
                  up(ii) = ( 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))*d(i,j,k) +
     &                 0.5_RKIND*(u(i-1,j,k)+u(i,j,k))*d(i-1,j,k)      +
     &                 0.5_RKIND*(u(i+1,j,k)+u(i+2,j,k))*d(i+1,j,k)    +
     &                 0.5_RKIND*(u(i,j+1,k)+u(i+1,j+1,k))*d(i,j+1,k)  +        
     &                 0.5_RKIND*(u(i,j-1,k)+u(i+1,j-1,k))*d(i,j-1,k)  + 
     &                 0.5_RKIND*(u(i,j,k+1)+u(i+1,j,k+1))*d(i,j,k+1)  +        
     &                 0.5_RKIND*(u(i,j,k-1)+u(i+1,j,k-1))*d(i,j,k-1) )/ 
     &                 ( d(i,j,k)+d(i-1,j,k)+d(i+1,j,k) +
     &                 d(i,j-1,k)+d(i,j+1,k)           +
     &                 d(i,j,k-1)+d(i,j,k+1) ) 
                  
                  vp(ii) = ( 0.5_RKIND*(v(i,j,k)+v(i+1,j,k))*d(i,j,k) +
     &                 0.5_RKIND*(v(i-1,j,k)+v(i,j,k))*d(i-1,j,k)      +
     &                 0.5_RKIND*(v(i+1,j,k)+v(i+2,j,k))*d(i+1,j,k)    +
     &                 0.5_RKIND*(v(i,j+1,k)+v(i+1,j+1,k))*d(i,j+1,k)  +
     &                 0.5_RKIND*(v(i,j-1,k)+v(i+1,j-1,k))*d(i,j-1,k)  +
     &                 0.5_RKIND*(v(i,j,k+1)+v(i+1,j,k+1))*d(i,j,k+1)  +
     &                 0.5_RKIND*(v(i,j,k-1)+v(i+1,j,k-1))*d(i,j,k-1) )/
     &                 ( d(i,j,k)+d(i-1,j,k)+d(i+1,j,k) +
     &                 d(i,j-1,k)+d(i,j+1,k)           +
     &                 d(i,j,k-1)+d(i,j,k+1) )
                  
                  wp(ii) = ( 0.5_RKIND*(w(i,j,k)+w(i+1,j,k))*d(i,j,k) +
     &                 0.5_RKIND*(w(i-1,j,k)+w(i,j,k))*d(i-1,j,k)      +
     &                 0.5_RKIND*(w(i+1,j,k)+w(i+2,j,k))*d(i+1,j,k)    +
     &                 0.5_RKIND*(w(i,j+1,k)+w(i+1,j+1,k))*d(i,j+1,k)  +
     &                 0.5_RKIND*(w(i,j-1,k)+w(i+1,j-1,k))*d(i,j-1,k)  +
     &                 0.5_RKIND*(w(i,j,k+1)+w(i+1,j,k+1))*d(i,j,k+1)  +
     &                 0.5_RKIND*(w(i,j,k-1)+w(i+1,j,k-1))*d(i,j,k-1) )/
     &                 ( d(i,j,k)+d(i-1,j,k)+d(i+1,j,k) +
     &                 d(i,j-1,k)+d(i,j+1,k)           +
     &                 d(i,j,k-1)+d(i,j,k+1) )                  
                  
               else
                  
                  up(ii) = (u(i,j,k)*d(i,j,k) +
     &                 u(i-1,j,k)*d(i-1,j,k) +
     &                 u(i+1,j,k)*d(i+1,j,k) +
     &                 u(i,j-1,k)*d(i,j-1,k) +
     &                 u(i,j+1,k)*d(i,j+1,k) +
     &                 u(i,j,k-1)*d(i,j,k-1) +
     &                 u(i,j,k+1)*d(i,j,k+1) ) / 
     &                 ( d(i,j,k) + d(i-1,j,k) + d(i+1,j,k) +
     &                 d(i,j-1,k) + d(i,j+1,k) + d(i,j,k-1) +
     &                 d(i,j,k+1) )
                  
                  vp(ii) = (v(i,j,k)*d(i,j,k) +
     &                 v(i-1,j,k)*d(i-1,j,k) +
     &                 v(i+1,j,k)*d(i+1,j,k) +
     &                 v(i,j-1,k)*d(i,j-1,k) +
     &                 v(i,j+1,k)*d(i,j+1,k) +
     &                 v(i,j,k-1)*d(i,j,k-1) +
     &                 v(i,j,k+1)*d(i,j,k+1) ) / 
     &                 ( d(i,j,k) + d(i-1,j,k) + d(i+1,j,k) +
     &                 d(i,j-1,k) + d(i,j+1,k) + d(i,j,k-1) +
     &                 d(i,j,k+1) )
                  
                  wp(ii) = (w(i,j,k)*d(i,j,k) +
     &                 w(i-1,j,k)*d(i-1,j,k) +
     &                 w(i+1,j,k)*d(i+1,j,k) +
     &                 w(i,j-1,k)*d(i,j-1,k) +
     &                 w(i,j+1,k)*d(i,j+1,k) +
     &                 w(i,j,k-1)*d(i,j,k-1) +
     &                 w(i,j,k+1)*d(i,j,k+1) ) / 
     &                 ( d(i,j,k) + d(i-1,j,k) + d(i+1,j,k) +
     &                 d(i,j-1,k) + d(i,j+1,k) + d(i,j,k-1) +
     &                 d(i,j,k+1) )
                  
               endif               
#else
               if (imethod .eq. 2) then
                  up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
                  vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
                  wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))
               else
                  up(ii) = u(i,j,k)
                  vp(ii) = v(i,j,k)
                  wp(ii) = w(i,j,k)
               endif
#endif

c
c              Set the particle metal fraction
c
               if (imetal .eq. 1) then
                  metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
               else
                  metalf(ii) = 0._RKIND
               endif
c
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
      
      if (ii .ge. nmax) then
         write(6,*) 'star_maker_mqk: reached max new particle count'
         ERROR_MESSAGE
      endif
      np = ii

#ifdef STAR_LOG
      close(4)
#endif

      return
      end


      R_PREC function SobolevColumn(d,fHI,nx,ny,nz,i,j,k,dx,d1,x1)
c     
c  CALCULATE SOBOLEV-LIKE APPROXIMATION OF THE HI SURFACE DENSITY
c
c  written by:  Michael Kuhlen
c  date:  May 2011
c
c  Methodology:
c
c    rho^+ = 0.5 * (rho(i,j,k) + rho(i+1,j,k))
c    h^+ = rho^+ / abs[(rho(i+1,j,k) - rho(i,j,k)) / dx]
c    Sigma_x^+(i,j,k) = rho^+ * h^+
c
c    rho^- = 0.5 * (rho(i,j,k) + rho(i-1,j,k))
c    h^- = rho^- / abs[(rho(i,j,k) - rho(i-1,j,k)) / dx]
c    Sigma_x^-(i,j,k) = rho^- * h^-
c
c    and similarly for y, z directions.
c
c    Calculate Sigma by taking the harmonic average over the cardinal
c    directions:
c
c    <Sigma> = 6 / [ 1/Sigma_x^+ + 1/Sigma_x^- + 
c                    1/Sigma_y^+ + 1/Sigma_y^- + 
c                    1/Sigma_z^+ + 1/Sigma_z^- ]
c
c
c  INPUTS:
c
c    d         - density field
c    fHI       - neutral hydrogen mass fraction (rho_HI / rho_tot)
c    nx,ny,nz  - dimensions of field arrays
c    i,j,k     - index at which column is to be calculated
c    dx        - zone size (code units)
c    d1,x1     - factors to convert d, dx to physical units
c
c  RETURN:
c    Sigma_HI - Sobolev-like estimate of the neutral hydrogen column density
c
c-----------------------------------------------------------------------
      implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c  Arguments
c
      INTG_PREC nx, ny, nz, i, j, k
      R_PREC d(nx,ny,nz), fHI(nx,ny,nz)
      R_PREC dx, d1, x1
c
c  Locals:
c
      R_PREC dHI
      R_PREC dHI_xm1, dHI_xp1
      R_PREC dHI_ym1, dHI_yp1
      R_PREC dHI_zm1, dHI_zp1

      R_PREC grad_xm, grad_xp
      R_PREC grad_ym, grad_yp
      R_PREC grad_zm, grad_zp

      R_PREC dHIavg_xm, dHIavg_xp
      R_PREC dHIavg_ym, dHIavg_yp
      R_PREC dHIavg_zm, dHIavg_zp

      R_PREC Sigma
      R_PREC Sigma_xm, Sigma_xp
      R_PREC Sigma_ym, Sigma_yp
      R_PREC Sigma_zm, Sigma_zp

      R_PREC unit_factor


c     HI density at grid cell of interest and six surrounding neighbor
c     cells.
      dHI = d(i,j,k) * fHI(i,j,k)

      dHI_xm1 = d(i-1,j,k) * fHI(i-1,j,k)
      dHI_xp1 = d(i+1,j,k) * fHI(i+1,j,k)

      dHI_ym1 = d(i,j-1,k) * fHI(i,j-1,k)
      dHI_yp1 = d(i,j+1,k) * fHI(i,j+1,k)

      dHI_zm1 = d(i,j,k-1) * fHI(i,j,k-1)
      dHI_zp1 = d(i,j,k+1) * fHI(i,j,k+1)


c     Calculate "gradients"
      grad_xm = abs( dHI_xm1 - dHI ) / dx
      grad_xp = abs( dHI_xp1 - dHI ) / dx

      grad_ym = abs( dHI_ym1 - dHI ) / dx
      grad_yp = abs( dHI_yp1 - dHI ) / dx

      grad_zm = abs( dHI_zm1 - dHI ) / dx
      grad_zp = abs( dHI_zp1 - dHI ) / dx


c     Check for zero gradients
      if( grad_xm .eq. 0._RKIND ) grad_xm = 1e-20_RKIND
      if( grad_xp .eq. 0._RKIND ) grad_xp = 1e-20_RKIND

      if( grad_ym .eq. 0._RKIND ) grad_ym = 1e-20_RKIND
      if( grad_yp .eq. 0._RKIND ) grad_yp = 1e-20_RKIND

      if( grad_zm .eq. 0._RKIND ) grad_zm = 1e-20_RKIND
      if( grad_zp .eq. 0._RKIND ) grad_zp = 1e-20_RKIND


c     Calculate averaged densities at (i,j,k)
      dHIavg_xm = 0.5_RKIND * ( dHI_xm1 + dHI )
      dHIavg_xp = 0.5_RKIND * ( dHI_xp1 + dHI )

      dHIavg_ym = 0.5_RKIND * ( dHI_ym1 + dHI )
      dHIavg_yp = 0.5_RKIND * ( dHI_yp1 + dHI )

      dHIavg_zm = 0.5_RKIND * ( dHI_zm1 + dHI )
      dHIavg_zp = 0.5_RKIND * ( dHI_zp1 + dHI )


c     Convert from code units to Msun pc^-2
c     (1 g/cm^2 = 4788.03 Msun/pc^2)
      unit_factor = d1 * x1 * 4788.03_RKIND


c     Calculate directional columns
      Sigma_xm = dHIavg_xm * dHIavg_xm / grad_xm * unit_factor
      Sigma_xp = dHIavg_xp * dHIavg_xp / grad_xp * unit_factor

      Sigma_ym = dHIavg_ym * dHIavg_ym / grad_ym * unit_factor
      Sigma_yp = dHIavg_yp * dHIavg_yp / grad_yp * unit_factor

      Sigma_zm = dHIavg_zm * dHIavg_zm / grad_zm * unit_factor
      Sigma_zp = dHIavg_zp * dHIavg_zp / grad_zp * unit_factor


c     Total column is the harmonic mean of directional columns
      SobolevColumn = 6._RKIND / ( 1._RKIND/Sigma_xm + 1._RKIND/Sigma_xp 
     $                      + 1._RKIND/Sigma_ym + 1._RKIND/Sigma_yp +
     $                        1._RKIND/Sigma_zm + 1._RKIND/Sigma_zp )
      

      return
      end

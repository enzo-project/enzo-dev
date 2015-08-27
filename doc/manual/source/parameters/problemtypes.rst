Problem Type Parameters
-----------------------

``ProblemType`` (external)
    This integer specifies the type of problem to be run. Its value
    causes the correct problem initializer to be called to set up the
    grid, and also may trigger certain boundary conditions or other
    problem-dependent routines to be called. The possible values are
    listed below. Default: none. 

For other problem-specific parameters follow the links below.  The problems
marked with "hydro_rk" originate from the MUSCL solver package in the enzo installation directory
``src/enzo/hydro_rk``.  For the 4xx radiation hydrodynamics problem types, see
the user guides in the installation directory ``doc/implicit_fld`` and ``doc/split_fld``.

============ ====================================
Problem Type Description and Parameter List
============ ====================================
1 	     :ref:`shocktube_param`
2	     :ref:`wavepool_param`
3 	     :ref:`shockpool_param`
4 	     :ref:`doublemach_param`
5 	     :ref:`shockinabox_param`
6 	     :ref:`implosion_param`
7 	     :ref:`sedovblast_param`
8 	     :ref:`khinstability_param`
9 	     :ref:`noh_param`
10 	     :ref:`rotatingcylinder_param`
11 	     :ref:`radiatingshock_param`
12 	     :ref:`freeexpansion_param`
14           :ref:`rotatingsphere_param`
20 	     :ref:`zeldovichpancake_param`
21 	     :ref:`pressurelesscollapse_param`
22 	     :ref:`adiabaticexpansion_param`
23 	     :ref:`testgravity_param`
24 	     :ref:`sphericalinfall_param`
25 	     :ref:`testgravitysphere_param`
26 	     :ref:`gravityequilibriumtest_param`
27 	     :ref:`collapsetest_param`
28 	     :ref:`testgravitymotion_param`
29 	     :ref:`testorbit_param`
30 	     :ref:`cosmologysimulation_param`
31 	     :ref:`galaxysimulation_param`
35 	     :ref:`shearingbox_param`
36	     Shearing Box 2D Simulation
37	     Stratifeid Shearing Box Simulation
40 	     :ref:`supernovarestart_param`
50 	     :ref:`photontest_param`
51	     Photon Test Restart
59       :ref:`stochastic_forcing_param`
60 	     :ref:`turbulence_param` 
61 	     :ref:`protostellar_param` 
62 	     :ref:`coolingtest_param`
63           One Zone Free Fall Test
70	     Conduction Test with Hydro Off
71	     Conduction Test with Hydro On
72	     Conduction Bubble Test
73	     Conduction Cloud Test
80           Explosion in a Stratified Medium Test
101          :ref:`3dcollapse_param`
102          :ref:`1dcollapse_param`
106          :ref:`mhdhydro_param`
107          :ref:`putsink_param`
108          :ref:`clustercoolingflow_param` 
200          :ref:`mhd1d_param`
201          :ref:`mhd2d_param`
202          :ref:`mhd3d_param`
203          :ref:`mhdtcollapse_param`
204          3D MHD Test
207          :ref:`galaxydisk_param`
208          :ref:`agndisk_param`
209	     MHD 1D Waves
210	     MHD Decaying Random Magnetic Fields
300          :ref:`poissonsolver_param`
400          :ref:`rhdtest1_param`
401          :ref:`rhdtest2_param`
402          :ref:`rhdtest3_param`
403          :ref:`rhdtest4_param`
404/405      :ref:`rhdtest5_param`
410/411	     :ref:`rhdtest10_param`
412 	     :ref:`rhdtest12_param`
413 	     :ref:`rhdtest13_param`
414/415	     :ref:`rhdtest14_param`
450-452	     Free-streaming radiation tests (to be removed)
============ ====================================

.. _shocktube_param:

Shock Tube (1: unigrid and AMR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Riemann problem or arbitrary discontinuity breakup problem. The
    discontinuity initially separates two arbitrary constant states:
    Left and Right. Default values correspond to the so called Sod
    Shock Tube setup (test 1.1). A table below contains a series of
    recommended 1D tests for hydrodynamic method, specifically designed
    to test the performance of the Riemann solver, the treatment of
    shock waves, contact discontinuities, and rarefaction waves in a
    variety of situations (Toro 1999, p. 129).

    It is also possible to set up a second discontinuity, creating three
    initial regions, rather than the two regions of the original Sod Shock
    Tube.

    ::

              Test  LeftDensity LeftVelocity LeftPressure RightDensity RightVelocity RightPressure
              1.1   1.0         0.0          1.0          0.125        0.0           0.1
              1.2   1.0         -2.0         0.4          1.0          2.0           0.4
              1.3   1.0         0.0          1000.0       1.0          0.0           0.01
              1.4   1.0         0.0          0.01         1.0          0.0           100.0
              1.5   5.99924     19.5975      460.894      5.99242      -6.19633      46.0950


``HydroShockTubesInitialDiscontinuity`` (external)
    The position of the initial discontinuity. Default: 0.5
``HydroShockTubesSecondDiscontinuity`` (external)
    The position of the second discontinuity, if a second discontinuity is 
    desired. Default: FLOAT_UNDEFINED, i.e. no second discontinuity.
``HydroShockTubesLeftDensity``, ``HydroShockTubesRightDensity``, ``HydroShockTubesCenterDensity`` (external)
    The initial gas density to the left and right of the discontinuity,
    and between the discontinuities if a second discontinuity has been 
    specified with HydroShockTubesSecondDiscontinuity.  Default: 1.0 for each
    value.
``HydroShockTubesLeftPressure``, ``HydroShockTubesRightPressure``, ``HydroShockTubesCenterPressure`` (external)
    The initial gas density to the left and right of the discontinuity,
    and between the discontinuities if a second discontinuity has been
    specified with HydroShockTubesSecondDiscontinuity.  Default: 1.0 for
    each of the left, right, and center regions.

``HydroShockTubesLeftVelocityX``, ``HydroShockTubesLeftVelocityY``, ``HydroShockTubesLeftVelocityZ`` (external)
    The initial gas velocity, in the x-, y-, and z-directions to the left of 
    the discontinuity.  Default: 0.0 for all directions.

``HydroShockTubesRightVelocityX``, ``HydroShockTubesRightVelocityY``, ``HydroShockTubesRightVelocityZ`` (external)
    The initial gas velocity, in the x-, y-, and z-directions to the right of 
    the discontinuity.  Default: 0.0 for all directions.

``HydroShockTubesCenterVelocityX``, ``HydroShockTubesCenterVelocityY``, ``HydroShockTubesCenterVelocityZ`` (external)
    The initial gas velocity, in the x-, y-, and z-directions between the 
    discontinuities, used if a second discontinuity has been specified with 
    HydroShockTubesSecondDiscontinuity. Default: 1.0 for all directions.

.. _wavepool_param:

Wave Pool (2)
~~~~~~~~~~~~~

    Wave Pool sets up a simulation with a 1D sinusoidal wave entering
    from the left boundary. The initial active region is uniform and
    the wave is entered via inflow boundary conditions.


``WavePoolAmplitude`` (external)
    The amplitude of the wave. Default: 0.01 - a linear wave.
``WavePoolAngle`` (external)
    Direction of wave propagation with respect to x-axis. Default: 0.0
``WavePoolDensity`` (external)
    Uniform gas density in the pool. Default: 1.0
``WavePoolNumberOfWaves`` (external)
    The test initialization will work for one wave only. Default: 1
``WavePoolPressure`` (external)
    Uniform gas pressure in the pool. Default: 1.0
``WavePoolSubgridLeft``, ``WavePoolSubgridRight`` (external)
    Start and end positions of the subgrid. Default: 0.0 and 0.0 (no
    subgrids)
``WavePoolVelocity1(2,3)`` (external)
    x-,y-, and z-velocities. Default: 0.0 (for all)
``WavePoolWavelength`` (external)
    The wavelength. Default: 0.1 (one-tenth of the box)

.. _shockpool_param:

Shock Pool (3: unigrid 2D, AMR 2D and unigrid 3D)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The Shock Pool test sets up a system which introduces a shock from
    the left boundary. The initial active region is uniform, and the
    shock wave enters via inflow boundary conditions. 2D and 3D
    versions available. (D. Mihalas & B.W. Mihalas, Foundations of
    Radiation Hydrodynamics, 1984, p. 236, eq. 56-40.)


``ShockPoolAngle`` (external)
    Direction of the shock wave propagation with respect to x-axis.
    Default: 0.0
``ShockPoolDensity`` (external)
    Uniform gas density in the preshock region. Default: 1.0
``ShockPoolPressure`` (external)
    Uniform gas pressure in the preshock region. Default: 1.0
``ShockPoolMachNumber`` (external)
    The ratio of the shock velocity and the preshock sound speed.
    Default: 2.0
``ShockPoolSubgridLeft``, ``ShockPoolSubgridRight`` (external)
    Start and end positions of the subgrid. Default: 0.0 and 0.0 (no
    subgrids)
``ShockPoolVelocity1(2,3)`` (external)
    Preshock gas velocity (the Mach number definition above assumes a
    zero velocity in the laboratory reference frame. Default: 0.0 (for
    all components)

.. _doublemach_param:

Double Mach Reflection (4)
~~~~~~~~~~~~~~~~~~~~~~~~~~

    A test for double Mach reflection of a strong shock (Woodward &
    Colella 1984). Most of the parameters are "hardwired": d0 = 8.0, e0
    = 291.25, u0 = 8.25\*sqrt(3.0)/2.0, v0 = -8.25\*0.5, w0 = 0.0


``DoubleMachSubgridLeft`` (external)
    Start position of the subgrid. Default: 0.0
``DoubleMachSubgridRight`` (external)
    End positions of the subgrid. Default: 0.0

.. _shockinabox_param:

Shock in a Box (5)
~~~~~~~~~~~~~~~~~~

    A stationary shock front in a static 3D subgrid (Anninos et al.
    1994). Initialization is done as in the Shock Tube test.


``ShockInABoxBoundary`` (external)
    Position of the shock. Default: 0.5
``ShockInABoxLeftDensity``, ``ShockInABoxRightDensity`` (external)
    Densities to the right and to the left of the shock front. Default:
    ``dL=1.0`` and ``dR = dL*((Gamma+1)*m^2)/((Gamma-1)*m^2 + 2)``, where
    ``m=2.0`` and ``speed=0.9*sqrt(Gamma*pL/dL)*m``.
``ShockInABoxLeftVelocity``, ``ShockInABoxRightVelocity`` (external)
    Velocities to the right and to the left of the shock front.
    Default: ``vL=shockspeed`` and
    ``vR=shockspeed-m*sqrt(Gamma*pL/dL)*(1-dL/dR)``, where ``m=2.0``,
    ``shockspeed=0.9*sqrt(Gamma*pL/dL)*m``.
``ShockInABoxLeftPressure``, ``ShockInABoxRightPressure`` (external)
    Pressures to the Right and to the Left of the shock
    front. Default: pL=1.0 and pR=pL*(2.0*Gamma*m^2 -
    (Gamma-1))/(Gamma+1), where m=2.0.
``ShockInABoxSubgridLeft``, ``ShockInABoxSubgridRight`` (external)
    Start and end positions of the subgrid. Default: 0.0 (for both)

.. _implosion_param:

Implosion (6)
~~~~~~~~~~~~~
 
    The implosion test sets up a converging shock problem in a square domain
    (x,y) \in (0, 0.3)x(0, 0.3) with gas initially at rest. Initial
    pressure and density is 1 everywhere except for a triangular region
    (0.15,0)(0.15,0) where d=0.125 and p=0.14. Reflecting boundary conditions
    at all boundaries. Adiabatic index gamma=1.4.
     
    If AMR is used, a hierarchy of subgrids (one per level) will be generated
    at start-up to properly resolve the initial discontinuity.
                      
    REFERENCE: Hui Li and Z. Li, JCP 153, 596, 1999.
               Chang et al. JCP 160, 89, 1999.



``ImplosionDensity`` (external)
   Initial density. Default: 1.0
``ImplosionPressure`` (external)
   Initial pressure. Default: 1.0
``ImplosionDimaondDensity`` (external)
   Initial density within diamond. Default: 0.125
``ImplosionDimaondPressure`` (external)
   Initial pressure within diamond. Default: 0.14
``ImplosionSubgridLeft``, ``ImplosionSubgridRight`` (external)
   Start and position of the subgrid. Default: 0.0 (for both)

.. _sedovblast_param:

Sedov Blast (7)
~~~~~~~~~~~~~~~

     Self-similar solution: L.I. Sedov (1946); 
     see also: Sedov (1959), Similarity and Dimensional Methods
     in Mechanics, pp. 210, 219, 228;
     see also: Landau & Lifshitz, Fluid Dynamics, Sect. 99 
     "The Propagation of Strong Shock Waves" (1959).
     Experiments, terrestrial/numerical: Taylor (1941, 1949).


``SedovBlastFullBox`` (external)
    Full box or one quadrant. Default: 0
``SedovBlastType`` (external)
    2D. Default: 0
``SedovBlastInitialTime`` (external)
    Initial time. Default: 0
``SedovBlastDensity`` (external)
    Initial density. Default: 1.0
``SedovBlastPressure`` (external)
    Initial pressure. Default: 1e-5
``SedovBlastInputEnergy`` (external)
    Energy input into system. Default: 1.0
``SedovBlastEnergyZones`` (external)
    Default: 3.5
``SedovBlastSubGridLeft``, ``SedovBlastSubGridRight`` (external)
    Start and end position of the subgrid. Default: 0.0 (for both)

.. _khinstability_param:

Kelvin-Helmholtz Instability (8)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This problem sets up a 2D box with periodic boundary conditions containing
    two fluids (inner fluid and outer fluid).  The inner fluid has a positive
    velocity and the outer fluid has a negative velocity with a difference of
    ``KHVelocityJump``.  The two fluids typically have different densities.
    The result is the build up of KH instabilities along the interface between
    the two fluids.

    Setting ``KHRamp`` to 0, creates the standard KH test problem
    where there is a discontinuous jump between the two fluids in
    x-velocity and density.  Random perturbations in y-velocity are the seeds 
    to the KH instability resulting in growth of multiple modes of the KHI.

    Setting ``KHRamp`` to 1 modifies the ICs so that there is a smooth
    ramp connecting the two fluids in x-velocity and density of width 
    ``KHRampWidth``.  A sinusoidal perturbation in y-velocity is the seed
    to the KH instability resulting in only growth of k=2 modes.  
    These results converge in behavior as resolution is increased, whereas 
    the standard ICs do not.  The ramped ICs are based on Robertson, Kravtsov, 
    Gnedin, Abel & Rudd 2010, but that work has a typo in the ramp equation, 
    and this implementation matches Robertson's actual ICs.  

``KHInnerDensity``, ``KHOuterDensity`` (external)
    Initial density. Default: 2.0 (inner) and 1.0 (outer)
``KHInnerPressure``, ``KHOuterPressure`` (external)
    Initial pressure. Default: 2.5 (for both)
``KHBulkVelocity`` (external)
    The bulk velocity of both fluids relative to the grid.  Default: 0.0
``KHVelocityJump`` (external)
    The difference in velocity between the outer fluid and the inner fluid.
    Inner fluid will have half this value and move to the right (positive),
    whereas outer fluid will have have this value and move to the left 
    (negative).  Total fluid velocities will combine this jump with 
    KHBulkVelocity.  Default: 1.0
``KHPerturbationAmplitude`` (external)
    Default: 0.1
``KHRamp`` (external)
    Whether to use ramped ICs or not.  Default: 1
``KHRampWidth`` (external)
    The width in y-space of the transition ramp.  Default: 0.05
``KHRandomSeed`` (external)
    The seed for the Mersennes random number generator.  This is only
    used in the case of the KHRamp=0 ICs.  By using the same seed
    from one run to the next, one can reproduce previous behavior with
    identical parameter files.  Default: 123456789


.. _noh_param:

2D/3D Noh Problem (9)
~~~~~~~~~~~~~~~~~~~~~
     
    Liska & Wendroff, 2003, SIAM J. Sci. Comput. 25, 995, 
    Section 4.5, Fig. 4.4.


``NohProblemFullBox`` (external)
    Default: 0
``NohSubgridLeft``, ``NohSubgridRight`` (external)
    Start and end positon of the subgrid. Default: 0.0 (for both)


.. _rotatingcylinder_param:

Rotating Cylinder (10)
~~~~~~~~~~~~~~~~~~~~~~

    A test for the angular momentum conservation of a collapsing
    cylinder of gas in an AMR simulation. Written by Brian O'Shea
    (`oshea@msu.edu <mailto:oshea@msu.edu>`_).


``RotatingCylinderOverdensity`` (external)
    Density of the rotating cylinder with respect to the
    background. Default: 20.0
``RotatingCylinderSubgridLeft``, ``RotatingCylinderSubgridRight`` (external)
    This pair of floating point numbers creates a subgrid region at the
    beginning of the simulation that will be refined to
    ``MaximumRefinementLevel``. It should probably encompass the whole
    cylinder. Positions are in units of the box, and it always creates
    a cube. No default value (meaning off).
``RotatingCylinderLambda`` (external)
    Angular momentum of the cylinder as a dimensionless quantity. This
    is identical to the angular momentum parameter lambda that is
    commonly used to describe cosmological halos. A value of 0.0 is
    non-rotating, and 1.0 means that the gas is already approximately
    rotating at the Keplerian value. Default: 0.05
``RotatingCylinderTotalEnergy`` (external)
    Sets the default gas energy of the ambient medium, in Enzo internal
    units. Default: 1.0
``RotatingCylinderRadius`` (external)
    Radius of the rotating cylinder in units of the box size. Note that
    the height of the cylinder is equal to the diameter. Default: 0.3
``RotatingCylinderCenterPosition`` (external)
    Position of the center of the cylinder as a vector of floats.
    Default: (0.5, 0.5, 0.5)

.. _radiatingshock_param:

Radiating Shock (11)
~~~~~~~~~~~~~~~~~~~~

    This is a test problem similar to the Sedov test problem documented
    elsewhere, but with radiative cooling turned on (and the ability to
    use ``MultiSpecies`` and all other forms of cooling). The main
    difference is that there are quite a few extras thrown in,
    including the ability to initialize with random density
    fluctuations outside of the explosion region, use a Sedov blast
    wave instead of just thermal energy, and some other goodies (as
    documented below).


``RadiatingShockInnerDensity`` (external)
    Density inside the energy deposition area (Enzo internal units).
    Default: 1.0
``RadiatingShockOuterDensity`` (external)
    Density outside the energy deposition area (Enzo internal units).
    Default: 1.0
``RadiatingShockPressure`` (external)
    Pressure outside the energy deposition area (Enzo internal units).
    Default: 1.0e-5
``RadiatingShockEnergy`` (external)
    Total energy deposited (in units of 1e51 ergs). Default: 1.0
``RadiatingShockSubgridLeft``, ``RadiatingShockSubgridRight`` (external)
    Pair of floats that defines the edges of the region where the
    initial conditions are refined to MaximumRefinementLevel. No
    default value.
``RadiatingShockUseDensityFluctuation`` (external)
    Initialize external medium with random density fluctuations.
    Default: 0
``RadiatingShockRandomSeed`` (external)
    Seed for random number geneator (currently using Mersenne Twister).
    Default: 123456789
``RadiatingShockDensityFluctuationLevel`` (external)
    Maximum fractional fluctuation in the density level. Default: 0.1
``RadiatingShockInitializeWithKE`` (external)
    Initializes the simulation with some initial kinetic energy if
    turned on (0 - off, 1 - on). Whether this is a simple sawtooth or a
    Sedov profile is controlled by the parameter
    ``RadiatingShockUseSedovProfile``. Default: 0
``RadiatingShockUseSedovProfile`` (external)
    If set to 1, initializes simulation with a Sedov blast wave profile
    (thermal and kinetic energy components). If this is set to 1, it
    overrides all other kinetic energy-related parameters. Default: 0
``RadiatingShockSedovBlastRadius`` (external)
    Maximum radius of the Sedov blast, in units of the box size.
    Default: 0.05
``RadiatingShockKineticEnergyFraction`` (external)
    Fraction of the total supernova energy that is deposited as kinetic
    energy. This only is used if ``RadiatingShockInitializeWithKE`` is set
    to 1. Default: 0.0
``RadiatingShockCenterPosition`` (external)
    Vector of floats that defines the center of the explosion. Default:
    (0.5, 0.5, 0.5)
``RadiatingShockSpreadOverNumZones`` (external)
    Number of cells that the shock is spread over. This corresponds to
    a radius of approximately N \* dx, where N is the number of cells
    and dx is the resolution of the highest level of refinement. This
    does not have to be an integer value. Default: 3.5

.. _freeexpansion_param:

Free Expansion (12)
~~~~~~~~~~~~~~~~~~~

This test sets up a blast wave in the free expansion stage. There
is only kinetic energy in the sphere with the radial velocity
proportional to radius. If let evolve for long enough, the problem
should turn into a Sedov-Taylor blast wave.

``FreeExpansionFullBox`` (external)
    Set to 0 to have the blast wave start at the origin with reflecting
    boundaries. Set to 1 to center the problem at the domain center
    with periodic boundaries. Default: 0
``FreeExpansionMass`` (external)
    Mass of the ejecta in the blast wave in solar masses. Default: 1
``FreeExpansionRadius`` (external)
    Initial radius of the blast wave. Default: 0.1
``FreeExpansionDensity`` (external)
    Ambient density of the problem. Default: 1
``FreeExpansionEnergy`` (external)
    Total energy of the blast wave in ergs. Default: 1e51
``FreeExpansionMaxVelocity`` (external)
    Maximum initial velocity of the blast wave (at the outer radius).
    If not set, a proper value is calculated using the formula in
    Draine & Woods (1991). Default: ``FLOAT_UNDEFINED``
``FreeExpansionTemperature`` (external)
    Ambient temperature of the problem in K. Default: 100
``FreeExapnsionBField`` (external)
    Initial uniform magnetic field. Default: 0 0 0
``FreeExpansionVelocity`` (external)
    Initial velocity of the ambient medium. Default: 0 0 0
``FreeExpansionSubgridLeft`` (external)
    Leftmost edge of the region to set the initial refinement. Default: 0
``FreeExpansionSubgridRight`` (external)
    Rightmost edge of the region to set the initial refinement.
    Default: 0

.. _rotatingsphere_param:

Rotating Sphere (14)
~~~~~~~~~~~~~~~~~~~~

A test originally created to study star formation. Sets up a rotating,
turbulent sphere of gas within an NFW halo. For details of the setup
process, see Meece (2014).


``RotatingSphereNFWMass`` (external)
    The mass of the NFW halo within R200 in solar masses.
    Default: 1.0e+7 M_sun
``RotatingSphereNFWConcentration`` (external)
    The NFW Concentration parameter, defined as virial radius over scale radius (R200/Rs).
    Default: 2.0
``RotatingSphereCoreRadius`` (external)
    Radius of the core region in code units. The core radius is used as the break in the
    density profile. Gas within the core is set up in HSE, while outside the core temperature
    increases adiabatically with density.
    Default: 16 pc
``RotatingSphereCentralDensity`` (external)
    This is the scaling density for the density profile in code units. The density profile is defined as
    rho(r) = rho_center * (r/Rc)^-alpha * (1+r/Rc)^(alpha-beta) where rho_center is this
    parameters, Rc is the core radius, alpha is the core exponent (below) and beta is the
    outer exponent (also below).
    Default: 1
``RotatingSphereCoreDensityExponent`` (external)
    The density scaling exponent in the core. Within the core, density approximately goes as
    (r/Rc)^-alpha, were alpha is this parameter.
    Default: 0.1
``RotatingSphereOuterDensityExponent`` (external)
    The density scaling exponent in the outer regions. Outside of the core, density
    approximately goes as (r/Rc)^-beta, were alpha is this parameter.
    Default: 2.5
``RotatingSphereExteriorTemperature`` (external)
    This is the temperature in K of gas outside the sphere, defined as the region where
    density would drop below the critical density.
    Default: 200.0
``RotatingSphereSpinParameter`` (external)
    The Baryonic spin parameter, defined as Lambda = (J * abs(E)^(1/2)) / (G M^(5/2)),
    where J is the total (gas) angular momentum, E is the binding energy of the gas due
    to the gas and dark matter, M is the gas mas, and G is the gravitational constant.
    All quantities are defined relative to the edge of the sphere defined above.
    Default: 0.05
``RotatingSphereAngularMomentumExponent`` (external)
    This is the power law index of the scaling relation for specific angular momentum
    as a function of mass enclosed. l scales as (M/M_T)^chi where chi is this parameter.
    Default: 0.9
``RotatingSphereUseTurbulence`` (external)
    0 = No Turbulence, 1 = Use Turbulence. If using turbulence, you need a file called
    turbulence.in, which can be generated using the file turbulence_generator.py in the
    RotatingSphere problem in the run directory.
    Default: 0
``RotatingSphereTurbulenceRMS`` (external)
    The RMS velocity of the turbulence is normalized to some fraction of the virial sound
    speed of the halo, as determined from the virial temperature of the halo. This parameter
    is that fraction. If RotatingSphereUseTurbulence == 0, this parameters is ignored.
    Default: 0.01
``RotatingSphereRedshift`` (external)
    The redshift is mainly used to determine the critical density of the universe. The problem
    generator assumes a cosmology with Omega_L=0.7, Omega_M = 0.3, and H0 = 70 km/s/mpc. Small
    variations in cosmology should not have a large effect on the properties of the sphere.
    Default: 20.0

.. _zeldovichpancake_param:

Zeldovich Pancake (20)
~~~~~~~~~~~~~~~~~~~~~~

    A test for gas dynamics, expansion terms and self-gravity in both
    linear and non-linear regimes [Bryan thesis (1996),
    Sect. 3.3.4-3.3.5; Norman & Bryan (1998), Sect. 4]


``ZeldovichPancakeCentralOffset`` (external)
    Offset of the pancake plane. Default: 0.0 (no offset)
``ZeldovichPancakeCollapseRedshift`` (external)
    A free parameter which determines the epoch of caustic formation.
    Default: 1.0
``ZeldovichPancakeDirection`` (external)
    Orientation of the pancake. Type: integer. Default: 0 (along the
    x-axis)
``ZeldovichPancakeInitialTemperature`` (external)
    Initial gas temperature. Units: degrees Kelvin. Default: 100
``ZeldovichPancakeOmegaBaryonNow`` (external)
    Omega Baryon at redshift z=0; standard setting. Default: 1.0
``ZeldovichPancakeOmegaCDMNow`` (external)
    Omega CDM at redshift z=0. Default: 0 (assumes no dark matter)

.. _pressurelesscollapse_param:

Pressureless Collapse (21)
~~~~~~~~~~~~~~~~~~~~~~~~~~

    An 1D AMR test for the gravity solver and advection routines: the
    two-sided one-dimensional collapse of a homogeneous plane parallel
    cloud in Cartesian coordinates. Isolated boundary conditions.
    Gravitational constant G=1; free fall time 0.399. The expansion
    terms are not used in this test. (Bryan thesis 1996, Sect. 3.3.1).


``PressurelessCollapseDirection`` (external)
    Coordinate direction. Default: 0 (along the x-axis).
``PressurelessCollapseInitialDensity`` (external)
    Initial density (the fluid starts at rest). Default: 1.0

.. _adiabaticexpansion_param:

Adiabatic Expansion (22)
~~~~~~~~~~~~~~~~~~~~~~~~

    A test for time-integration accuracy of the expansion terms (Bryan
    thesis 1996, Sect. 3.3.3).


``AdiabaticExpansionInitialTemperature`` (external)
    Initial temperature for Adiabatic Expansion test; test example
    assumes 1000 K. Default: 200. Units: degrees Kelvin
``AdiabaticExpansionInitialVelocity`` (external)
    Initial expansion velocity. Default: 100. Units: km/s
``AdiabaticExpansionOmegaBaryonNow`` (external)
    Omega Baryon at redshift z=0; standard value 1.0. Default: 1.0
``AdiabaticExpansionOmegaCDMNow`` (external)
    Omega CDM at redshift z=0; default setting assumes no dark matter.
    Default: 0.0

.. _testgravity_param:

Test Gravity (23)
~~~~~~~~~~~~~~~~~

    We set up a system in which there is one grid point with mass in
    order to see the resulting acceleration field. If finer grids are
    specified, the mass is one grid point on the subgrid as well.
    Periodic boundary conditions are imposed (gravity).


``TestGravityDensity`` (external)
    Density of the central peak. Default: 1.0
``TestGravityMotionParticleVelocity`` (external)
    Initial velocity of test particle(s) in x-direction. Default: 1.0
``TestGravityNumberOfParticles`` (external)
    The number of test particles of a unit mass. Default: 0
``TestGravitySubgridLeft``, ``TestGravitySubgridRight`` (external)
    Start and end positions of the subgrid. Default: 0.0 and 0.0 (no
    subgrids)
``TestGravityUseBaryons`` (external)
    Boolean switch. Type: integer. Default: 0 (FALSE)

.. _sphericalinfall_param:

Spherical Infall (24)
~~~~~~~~~~~~~~~~~~~~~

    A test based on Bertschinger's (1985) 3D self-similar spherical
    infall solution onto an initially overdense perturbation in an
    Einstein-de Sitter universe.


``SphericalInfallCenter`` (external)
    Coordinate(s) for the accretion center. Default: top grid center
``SphericalInfallFixedAcceleration`` (external)
    Boolean flag. Type: integer. Default: 0 (FALSE)
``SphericalInfallFixedMass`` (external)
    Mass used to calculate the acceleration from spherical infall
    (GM/(4*pi*r^3*a)). Default: If SphericalInfallFixedMass is
    undefined and ``SphericalInfallFixedAcceleration == TRUE``, then
    ``SphericalInfallFixedMass = SphericalInfallInitialPerturbation * TopGridVolume``
``SphericalInfallInitialPerturbation`` (external)
    The perturbation of initial mass density. Default: 0.1
``SphericalInfallOmegaBaryonNow`` (external)
    Omega Baryon at redshift z=0; standard setting. Default: 1.0
``SphericalInfallOmegaCDMNow`` (external)
    Omega CDM at redshift z=0. Default: 0.0 (assumes no dark matter)
    Default: 0.0
``SphericalInfallSubgridIsStatic`` (external)
    Boolean flag. Type: integer. Default: 0 (FALSE)
``SphericalInfallSubgridLeft``, ``SphericalInfallSubgridRight`` (external)
    Start and end positions of the subgrid. Default: 0.0 and 0.0 (no
    subgrids)
``SphericalInfallUseBaryons`` (external)
    Boolean flag. Type: integer. Default: 1 (TRUE)

.. _testgravitysphere_param:

Test Gravity: Sphere (25)
~~~~~~~~~~~~~~~~~~~~~~~~~

    Sets up a 3D spherical mass distribution and follows its evolution
    to test the gravity solver.


``TestGravitySphereCenter`` (external)
    The position of the sphere center. Default: at the center of the
    domain
``TestGravitySphereExteriorDensity`` (external)
    The mass density outside the sphere. Default: ``tiny_number``
``TestGravitySphereInteriorDensity`` (external)
    The mass density at the sphere center. Default: 1.0
``TestGravitySphereRadius`` (external)
    Radius of self-gravitating sphere. Default: 0.1
``TestGravitySphereRefineAtStart`` (external)
    Boolean flag. Type: integer. Default: 0 (FALSE)
``TestGravitySphereSubgridLeft``, ``TestGravitySphereSubgridRight`` (external)
    Start and end positions of the subgrid. Default: 0.0 and 0.0 (no
    subgrids)
``TestGravitySphereType`` (external)
    Type of mass density distribution within the sphere. Options
    include: (0) uniform density distrubution within the sphere radius;
    (1) a power law with an index -2.0; (2) a power law with an index
    -2.25 (the exact power law form is, e.g., r\ :sup:`-2.25`\ , where
    r is measured in units of ``TestGravitySphereRadius``). Default: 0
    (uniform density)
``TestGravitySphereUseBaryons`` (external)
    Boolean flag. Type: integer . Default: 1 (TRUE)

.. _gravityequilibriumtest_param:

Gravity Equilibrium Test (26)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Sets up a hydrostatic exponential atmosphere with the pressure=1.0
    and density=1.0 at the bottom. Assumes constant gravitational
    acceleration (uniform gravity field).


``GravityEquilibriumTestScaleHeight`` (external)
    The scale height for the exponential atmosphere . Default: 0.1

.. _collapsetest_param:

Collapse Test (27)
~~~~~~~~~~~~~~~~~~

    A self-gravity test.


``CollapseTestInitialTemperature`` (external)
    Initial gas temperature. Default: 1000 K. Units: degrees Kelvin
``CollapseTestNumberOfSpheres`` (external)
    Number of spheres to collapse; must be <= ``MAX_SPHERES=10`` (see
    ``Grid.h`` for definition). Default: 1
``CollapseTestRefineAtStart`` (external)
    Boolean flag. Type: integer. If TRUE, then initializing routine
    refines the grid to the desired level. Default: 1 (TRUE)
``CollapseTestUseColour`` (external)
    Boolean flag. Type: integer. Default: 0 (FALSE)
``CollapseTestUseParticles`` (external)
    Boolean flag. Type: integer. Default: 0 (FALSE)
``CollapseTestSphereCoreRadius`` (external)
    An array of core radii for collapsing spheres. Default: 0.1 (for
    all spheres)
``CollapseTestSphereDensity`` (external)
    An array of density values for collapsing spheres. Default: 1.0
    (for all spheres)
``CollapseTestSpherePosition`` (external)
    A two-dimensional array of coordinates for sphere centers. Type:
    float[``MAX_SPHERES``][``MAX_DIMENSION``]. Default for all spheres:
    0.5\*(``DomainLeftEdge[dim]`` + ``DomainRightEdge[dim]``)
``CollapseTestSphereRadius`` (external)
    An array of radii for collapsing spheres. Default: 1.0 (for all
    spheres)
``CollapseTestSphereTemperature`` (external)
    An array of temperatures for collapsing spheres. Default: 1.0.
    Units: degrees Kelvin
``CollapseTestSphereType`` (external)
    An integer array of sphere types. Default: 0
``CollapseTestSphereVelocity`` (external)
    A two-dimensional array of sphere velocities. Type:
    float[``MAX_SPHERES``][``MAX_DIMENSION``]. Default: 0.0
``CollapseTestUniformVelocity`` (external)
    Uniform velocity. Type: float[``MAX_DIMENSION``]. Default: 0 (for all
    dimensions)
``CollapseTestSphereMetallicity`` (external)
    Metallicity of the sphere in solar metallicity. Default: 0.
``CollapseTestFracKeplerianRot`` (external)
    Rotational velocity of the sphere in units of Keplerian velocity,
    i.e. 1 is rotationally supported. Default: 0.
``CollapseTestSphereTurbulence`` (external)
    Turbulent velocity field sampled from a Maxwellian distribution
    with the temperature specified in
    ``CollapseTestSphereTemperature``
    This parameter multiplies the turbulent velocities by its value.
    Default: 0.
``CollapseTestSphereDispersion`` (external)
    If using particles, this parameter multiplies the velocity
    dispersion of the particles by its value. Only valid in sphere type
    8 (cosmological collapsing sphere from a uniform density). Default:
    0.
``CollapseTestSphereCutOff`` (external)
    At what radius to terminate a Bonner-Ebert sphere. Units? Default:
    6.5
``CollapseTestSphereAng1`` (external)
    Controls the initial offset (at r=0) of the rotational axis. Units
    in radians. Default: 0.
``CollapseTestSphereAng2`` (external)
    Controls the outer offset (at ``r=SphereRadius`` of the rotational
    axis. In both ``CollapseTestSphereAng1`` and
    ``CollapseTestSphereAng2`` are set, the rotational axis linearly
    changes with radius between ``CollapseTestSphereAng1`` and
    ``CollapseTestSphereAng2``.  Units in radians. Default: 0.
``CollapseTestSphereConstantPressure`` (external)
    Constant pressure inside the sphere that is equal to the pressure
    at the outer radius.  Default: 0
``CollapseTestSphereSmoothSurface`` (external)
    The density interface between the ambient and sphere medium is
    smoothed with a hyperbolic tangent.  Default: 0
``CollapseTestSmoothRadius`` (external)
    The outer radius of the smoothed interface.  This parameter is in
    units of the sphere radius.  Default: 1.2
``CollapseTestSphereInitialLevel`` (external)
    Failed experiment to try to force refinement to a specified level.
    Not working. Default: 0.

.. _testgravitymotion_param:

Test Gravity Motion (28)
~~~~~~~~~~~~~~~~~~~~~~~~

``TestGravityMotionParticleVelocity`` (external)
    Initial velocity for particle. Default: 1.0

.. _testorbit_param:

Test Orbit (29)
~~~~~~~~~~~~~~~

``TestOrbitNumberOfParticles`` (external)
     Number of test particles. Default: 1
``TestOrbitRadius`` (external)
     Initial radius of orbit. Default: 0.2
``TestOrbitCentralMass`` (external)
     Central mass. Default: 1.0
``TestOrbitTestMass`` (external)
     Mass of the test particle. Default: 1.0e-6
``TestOrbitUseBaryons`` (external
     Boolean flag. (not implemented) Default: FALSE

.. _cosmologysimulation_param:

Cosmology Simulation (30)
~~~~~~~~~~~~~~~~~~~~~~~~~

    A sample cosmology simulation.


``CosmologySimulationDensityName`` (external)
    This is the name of the file which contains initial data for baryon
    density. Type: string. Example: ``GridDensity``. Default: none
``CosmologySimulationTotalEnergyName`` (external)
    This is the name of the file which contains initial data for total
    energy. Default: none
``CosmologySimulationGasEnergyName`` (external)
    This is the name of the file which contains initial data for gas
    energy. Default: none
``CosmologySimulationVelocity[123]Name`` (external)
    These are the names of the files which contain initial data for gas
    velocities. ``Velocity1`` - x-component; ``Velocity2`` - y-component;
    ``Velocity3`` - z-component. Default: none
``CosmologySimulationParticleMassName`` (external)
    This is the name of the file which contains initial data for
    particle masses. Default: none
``CosmologySimulationParticlePositionName`` (external)
    This is the name of the file which contains initial data for
    particle positions. Default: none
``CosmologySimulationParticleVelocityName`` (external)
    This is the name of the file which contains initial data for
    particle velocities. Default: none
``CosmologySimulationParticleVelocity[123]Name`` (external) This is
    the name of the file which contains initial data for particle
    velocities but only has one component per file. This is more
    useful with very large (>=2048\ :sup:`3`\ ) datasets. Currently
    one can only use this in conjunction with
    ``CosmologySimulationCalculatePositions``.  because it expects a
    3D grid structure instead of a 1D list of particles.  Default:
    None.
``CosmologySimulationCalculatePositions`` (external)
    If set to 1, Enzo will calculate the particle positions in one of
    two ways: 1) By using a linear Zeldo'vich approximation based on
    the particle velocities and a displacement factor [dln(growth
    factor) / dtau, where tau is the conformal time], which is stored
    as an attribute in the initial condition files, or 2) if the user
    has also defined either
    CosmologySimulationParticleDisplacementName or
    CosmologySimulationParticleDisplacement[123]Name, by reading in
    particle displacements from an external code and applying those
    directly.  The latter allows the use of non-linear displacements.
    Default: 0.
``CosmologySimulationParticleDisplacementName`` (external)
    This is the name of the file which contains initial data for
    particle displacements. Default: none
``CosmologySimulationParticleDisplacement[123]Name`` (external) This
    is the name of the file which contains initial data for particle
    displacements but only has one component per file. This is more
    useful with very large (>=2048\ :sup:`3`\ ) datasets. Currently
    one can only use this in conjunction with
    ``CosmologySimulationCalculatePositions``.  because it expects a
    3D grid structure instead of a 1D list of particles.  Default:
    None.
``CosmologySimulationNumberOfInitialGrids`` (external)
    The number of grids at startup. 1 means top grid only. If >1, then
    nested grids are to be defined by the following parameters.
    Default: 1
``CosmologySimulationSubgridsAreStatic`` (external)
    Boolean flag, defines whether the subgrids introduced at the
    startup are static or not. Type: integer. Default: 1 (TRUE)
``CosmologySimulationGridLevel`` (external)
    An array of integers setting the level(s) of nested subgrids. Max
    dimension ``MAX_INITIAL_GRIDS`` is defined in
    ``CosmologySimulationInitialize.C`` as 10. Default for all subgrids: 1,
    0 - for the top grid (grid #0)
``CosmologySimulationGridDimension[#]`` (external)
    An array (arrays) of 3 integers setting the dimensions of nested
    grids. Index starts from 1. Max number of subgrids
    ``MAX_INITIAL_GRIDS`` is defined in ``CosmologySimulationInitialize.C``
    as 10. Default: none
``CosmologySimulationGridLeftEdge[#]`` (external)
    An array (arrays) of 3 floats setting the left edge(s) of nested
    subgrids. Index starts from 1. Max number of subgrids
    ``MAX_INITIAL_GRIDS`` is defined in ``CosmologySimulationInitialize.C``
    as 10. Default: none
``CosmologySimulationGridRightEdge[#]`` (external)
    An array (arrays) of 3 floats setting the right edge(s) of nested
    subgrids. Index starts from 1. Max number of subgrids
    ``MAX_INITIAL_GRIDS`` is defined in ``CosmologySimulationInitialize.C``
    as 10. Default: none
``CosmologySimulationUseMetallicityField`` (external)
    Boolean flag. Type: integer. Default: 0 (FALSE)
``CosmologySimulationInitialFractionH2I`` (external)
    The fraction of molecular hydrogen (H_2) at ``InitialRedshift``. This
    and the following chemistry parameters are used if ``MultiSpecies`` is
    defined as 1 (TRUE). Default: 2.0e-20
``CosmologySimulationInitialFractionH2II`` (external)
    The fraction of singly ionized molecular hydrogen (H2+) at
    ``InitialRedshift``. Default: 3.0e-14
``CosmologySimulationInitialFractionHeII`` (external)
    The fraction of singly ionized helium at ``InitialRedshift``. Default:
    1.0e-14
``CosmologySimulationInitialFractionHeIII`` (external)
    The fraction of doubly ionized helium at ``InitialRedshift``. Default:
    1.0e-17
``CosmologySimulationInitialFractionHII`` (external)
    The fraction of ionized hydrogen at ``InitialRedshift``. Default:
    1.2e-5
``CosmologySimulationInitialFractionHM`` (external)
    The fraction of negatively charged hydrogen (H-) at
    ``InitialRedshift``. Default: 2.0e-9
``CosmologySimulationInitialFractionMetal`` (external)
    The fraction of metals at ``InitialRedshift``. Default: 1.0e-10
``CosmologySimulationInitialTemperature`` (external)
    A uniform temperature value at ``InitialRedshift`` (needed if the
    initial gas energy field is not supplied). Default: 550\*((1.0 +
    ``InitialRedshift``)/201)\ :sup:`2`\ 
``CosmologySimulationOmegaBaryonNow`` (external)
    This is the contribution of baryonic matter to the energy density
    at the current epoch (z=0), relative to the value required to
    marginally close the universe. Typical value 0.06. Default: 1.0
``CosmologySimulationOmegaCDMNow`` (external)
    This is the contribution of CDM to the energy density at the
    current epoch (z=0), relative to the value required to marginally
    close the universe. Typical value 0.24. Default: 0.0 (no dark
    matter)
``CosmologySimulationManuallySetParticleMassRatio`` (external)
    This binary flag (0 - off, 1 - on) allows the user to manually set
    the particle mass ratio in a cosmology simulation. Default: 0 (Enzo
    automatically sets its own particle mass)
``CosmologySimulationManualParticleMassRatio`` (external)
    This manually controls the particle mass in a cosmology simulation,
    when ``CosmologySimulationManuallySetParticleMassRatio`` is set to 1.
    In a standard Enzo simulation with equal numbers of particles and
    cells, the mass of a particle is set to
    ``CosmologySimulationOmegaCDMNow``/``CosmologySimulationOmegaMatterNow``,
    or somewhere around 0.85 in a WMAP-type cosmology. When a different
    number of particles and cells are used (128 particles along an edge
    and 256 cells along an edge, for example) Enzo attempts to
    calculate the appropriate particle mass. This breaks down when
    ``ParallelRootGridIO`` and/or ``ParallelParticleIO`` are turned on,
    however, so the user must set this by hand. If you have the ratio
    described above (2 cells per particle along each edge of a 3D
    simulation) the appropriate value would be 8.0 (in other words,
    this should be set to (number of cells along an edge) / (number of
    particles along an edge) cubed. Default: 1.0.

.. _galaxysimulation_param:

Isolated Galaxy Evolution (31)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Initializes an isolated galaxy, as per the Tasker & Bryan series of
    papers.


``GalaxySimulationRefineAtStart`` (external)
    Controls whether or not the simulation is refined beyond the root
    grid at initialization. (0 - off, 1 - on). Default: 1
``GalaxySimulationInitialRefinementLevel`` (external)
    Level to which the simulation is refined at initialization,
    assuming ``GalaxySimulationRefineAtStart`` is set to 1. Default: 0
``GalaxySimulationSubgridLeft``, ``GalaxySimulationSubgridRight`` (external)
    Vectors of floats defining the edges of the volume which is refined
    at start. No default value.
``GalaxySimulationUseMetallicityField`` (external)
    Turns on (1) or off (0) the metallicity field. Default: 0
``GalaxySimulationInitialTemperature`` (external)
    Initial temperature that the gas in the simulation is set to.
    Default: 1000.0
``GalaxySimulationUniformVelocity`` (external)
    Vector that gives the galaxy a uniform velocity in the ambient
    medium. Default: (0.0, 0.0, 0.0)
``GalaxySimulationDiskRadius`` (external)
    Radius (in Mpc) of the galax disk. Default: 0.2
``GalaxySimulationGalaxyMass`` (external)
    Dark matter mass of the galaxy, in Msun. Needed to initialize the
    NFW gravitational potential. Default: 1.0e+12
``GalaxySimulationGasMass`` (external)
    Amount of gas in the galaxy, in Msun. Used to initialize the
    density field in the galactic disk. Default: 4.0e+10
``GalaxySimulationDiskPosition`` (external)
    Vector of floats defining the center of the galaxy, in units of the
    box size. Default: (0.5, 0.5, 0.5)
``GalaxySimulationDiskScaleHeightz`` (external)
    Disk scale height, in Mpc. Default: 325e-6
``GalaxySimulationDiskScaleHeightR`` (external)
    Disk scale radius, in Mpc. Default: 3500e-6
``GalaxySimulationDarkMatterConcentrationParameter`` (external)
    NFW dark matter concentration parameter. Default: 12.0
``GalaxySimulationDiskTemperature`` (external)
    Temperature of the gas in the galactic disk. Default: 1.0e+4
``GalaxySimulationInflowTime`` (external)
    Controls inflow of gas into the box. It is strongly suggested that
    you leave this off. Default: -1 (off)
``GalaxySimulationInflowDensity`` (external)
    Controls inflow of gas into the box. It is strongly suggested that
    you leave this off. Default: 0.0
``GalaxySimulationAngularMomentum`` (external)
    Unit vector that defines the angular momentum vector of the galaxy
    (in other words, this and the center position define the plane of
    the galaxy). This _MUST_ be set! Default: (0.0, 0.0, 0.0)

.. _shearingbox_param:

Shearing Box Simulation (35)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``ShearingBoxProblemType`` (external)
    Value of 0 starts a sphere advection through the shearing box test.
    Value of 1 starts a standard Balbus & Hawley shearing box
    simulation. Default: 0
``ShearingBoxRefineAtStart`` (external)
    Refine the simulation at start. Default: 1.0
``ThermalMagneticRatio`` (external)
    Plasma beta (Pressure/Magnetic Field
    Energy) Default: 400.0
``FluctuationAmplitudeFraction`` (external)
    The magnitude of the sinusoidal velocity perturbations as a
    fraction of the angular velocity. Default: 0.1
``ShearingBoxGeometry`` (external)
    Defines the radius of the sphere for ``ShearingBoxProblemType`` =
    0, and the frequency of the velocity fluctuations (in units of
    2pi) for ``ShearingBoxProblemType`` = 1.  Default: 2.0

.. _supernovarestart_param:

Supernova Restart Simulation (40)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    All of the supernova parameters are to be put into a restart dump
    parameter file. Note that ProblemType must be reset to 40,
    otherwise these are ignored.

``SupernovaRestartEjectaCenter[#]`` (external)
    Input is a trio of coordinates in code units where the supernova's
    energy and mass ejecta will be centered. Default: ``FLOAT_UNDEFINED``
``SupernovaRestartEjectaEnergy`` (external)
    The amount of energy instantaneously output in the simulated
    supernova, in units of 1e51 ergs. Default: 1.0
``SupernovaRestartEjectaMass`` (external)
    The mass of ejecta in the supernova, in units of solar masses.
    Default: 1.0
``SupernovaRestartEjectaRadius`` (external)
    The radius over which the above two parameters are spread. This is
    important because if it's too small the timesteps basically go to
    zero and the simulation takes forever, but if it's too big then you
    loose information. Units are parsecs. Default: 1.0 pc
``SupernovaRestartName`` (external)
    This is the name of the restart data dump that the supernova
    problem is initializing from.
``SupernovaRestartColourField``
    Reserved for future use.

.. _photontest_param:

Photon Test (50)
~~~~~~~~~~~~~~~~

    This test problem is modeled after Collapse Test (27), and thus
    borrows all of its parameters that control the setup of spheres.
    Replace CollapseTest with PhotonTest in the sphere parameters, and
    it will be recognized. However there are parameters that control
    radiation sources, which makes this problem unique from collapse
    test. The radiation sources are fixed in space.


``PhotonTestNumberOfSources`` (external)
    Sets the number of radiation sources. Default: 1.
``PhotonTestSourceType`` (external)
    Sets the source type. No different types at the moment. Default: 0.
``PhotonTestSourcePosition`` (external)
    Sets the source position. Default: 0.5\*(``DomainLeftEdge`` + ``DomainRightEdge``)
``PhotonTestSourceLuminosity`` (external)
    Sets the source luminosity in units of photons per seconds.
    Default: 0.
``PhotonTestSourceLifeTime`` (external)
    Sets the lifetime of the source in units of code time. Default: 0.
``PhotonTestSourceRampTime`` (external)
    If non-zero, the source will exponentially increase its luminosity
    until it reaches the full luminosity when the age of the source
    equals this parameter. Default: 0.
``PhotonTestSourceEnergyBins`` (external)
    Sets the number of energy bins in which the photons are emitted
    from the source. Default: 4.
``PhotonTestSourceSED`` (external)
    An array with the fractional luminosity in each energy bin. The sum
    of this array must equal to one. Default: 1 0 0 0
``PhotonTestSourceEnergy`` (external)
    An array with the mean energy in each energy bin. Units are in eV.
    Default: 14.6 25.6 56.4 12.0 (i.e. HI ionizing, HeI ionizing, HeII
    ionizing, Lyman-Werner)
``PhotonTestSourceType`` (external)
    Indicates what radiation type (1 = isotropic, -2 = Beamed, -3 =
    Episodic). Default: 0
``PhotonTestSourceOrientation`` (external)
    Normal direction in Cartesian axes of beamed radiation (type =
    -2).  Default = 0 0 1
``PhotonTestInitialFractionHII`` (external)
    Sets the initial ionized fraction of hydrogen. Default: 1.2e-5
``PhotonTestInitialFractionHeII`` (external)
    Sets the initial singly-ionized fraction of helium. Default: 1e-14
``PhotonTestInitialFractionHeIII`` (external)
    Sets the initial doubly-ionized fraction of helium. Default: 1e-17
``PhotonTestInitialFractionHM`` (external)
    Sets the initial fraction of H\ :sup:`-`\ . Default: 2e-9
``PhotonTestInitialFractionH2I`` (external)
    Sets the initial neutral fraction of H2. Default: 2e-20
``PhotonTestInitialFractionH2II`` (external)
    Sets the initial ionized fraction of H2. Default: 3e-14
``PhotonTestOmegaBaryonNow`` (obsolete)
    Default: 0.05.

.. _stochastic_forcing_param:

Turbulence Simulation with Stochastic Forcing (59)
~~~~~~~~~~~~~~~~~~~~~~~~~~
    Typical quasi-isothermal "turbulence-in-a-box" problem with non-static driving field.
    For details on stochastic forcing, see Schmidt et al. 2009 A&A 494, 127-145 
    http://dx.doi.org/10.1051/0004-6361:200809967

    3D simulations with MUSCL hydro and MHD solver are tested.
    PPM, ZEUS and MHDCT unsupported at this time.

    Remember that in addition to the problem specific parameters below 
    UseDrivingField = 1 has to be turned on!


``DrivenFlowProfile`` (external)
    Shape of forcing power spectrum (1: delta peak, 2: band, 3: parabolic window).

``DrivenFlowAlpha`` (external)
    Ratio of domain length to integral length for each dimension (L = X/alpha).

``DrivenFlowBandWidth`` (external)
    Determines band width of the forcing spectrum relative to alpha (maximal value = 1).

``DrivenFlowMach`` (external)
    Characteristic velocity scale for each dimension (charcteristic force per unit mass F = V*V/L).

``DrivenFlowAutoCorrl`` (external)
    Determines autocorrelation time of the stochastic force in units of the integral time scale T = L/V.

``DrivenFlowWeight`` (external)
    Determines weight of solenoidal relative to dilatational modes (1 = purely solenoidal, 0 = purely dilatational).

``DrivenFlowSeed`` (external)
    Seed of random number generator.

``DrivenFlowDensity`` (external)
    Initial uniform density.

``DrivenFlowPressure`` (external)
    Initial uniform pressure.

``DrivenFlowMagField`` (external)
    Initial uniform magnetic field (x-direction)

.. _turbulence_param:

Turbulence Simulation (60)
~~~~~~~~~~~~~~~~~~~~~~~~~~

    Quasi-isothermal forced turbulence.

``TurbulenceSimulationsDensityName`` (external)

``TurbulenceSimulationTotalEnergyName`` (external)

``TurbulenceSimulationGasPressureName`` (external)

``TurbulenceSimulationGasEnergyName`` (external)

``TurbulenceSimulationVelocityName`` (external)

``TurbulenceSimulationRandomForcingName`` (external)

``TurbulenceSimulationMagneticName`` (external)

``TurbulenceSimulationInitialTemperature`` (external)    

``TurbulenceSimulationInitialDensity`` (external)

``TurbulenceSimulationSoundSpeed`` (external)

``TurbulenceSimulationInitialPressure`` (external)

``TurbulenceSimulationInitialDensityPerturbationAmplitude`` (external)

``TurbulenceSimulationNumberOfInitialGrids`` (external)
    Default: 1
``TurbulenceSimulationSubgridsAreStatic`` (external)
    Boolean flag. Default: 1
``TurbulenceSimulationGridLeftEdge[]`` (external)
    TBD
``TurbulenceSimulationGridRightEdge[]`` (external)
    TBD
``TurbulenceSimulationGridDimension[]`` (external)
    TBD
``TurbulenceSimulationGridLevel[]`` (external)
    TBD
``TurbulenceSimulationInitialMagneticField[i]`` (external)
    Initial magnetic field strength in the ith direction. Default: 5.0 (all)
``RandomForcing`` (external)
    This parameter is used to add random forcing field to create turbulence; see Mac Low 1999, ApJ 524, 169. Default: 0
``RandomForcingEdot`` (external)
    This parameter is used to define the value of such field; see TurbulenceSimulationInitialize.C and ComputeRandomForcingNormalization.C. Default: -1.0
``RandomForcingMachNumber`` (external)
    This parameter is used to define the value of such field; see Grid_TurbulenceSimulationInitialize.C and Grid_ComputeRandomForcingFields.C. Default: 0.0
``CycleSkipGlobalDataDump`` (external)
    Cycles to skip before global data (defined in ComputeRandomForcingNormalization.C) is dumped.

.. _protostellar_param:

Protostellar Collapse (61)
~~~~~~~~~~~~~~~~~~~~~~~~~~

     Bate 1998, ApJL 508, L95-L98

``ProtostellarCollapseCoreRadius`` (external)
     Radius of the core. Default: 0.005
``ProtostellarCollapseOuterDensity`` (external)
     Initial density. Default: 1.0
``ProtostellarCollapseAngularVelocity`` (external)
     Initial angular velocity. Default: 0
``ProtostellarCollapseSubgridLeft``, ``ProtostellarCollapseSubgridRight`` (external)
     Start and end position of subgrid. Default: 0 (for both)


.. _coolingtest_param:

Cooling Test (62)
~~~~~~~~~~~~~~~~~

    This test problem sets up a 3D grid varying smoothly in log-space in H
    number density (x dimension), metallicity (y-dimension), and temperature
    (z-dimension). The hydro solver is turned off. By varying the
    ``RadiativeCooling`` and ``CoolingTestResetEnergies`` parameters, two different
    cooling tests can be run. 1) Keep temperature constant, but iterate
    chemistry to allow species to converge. This will allow you to make plots
    of Cooling rate vs. T.  For this, set ``RadiativeCooling`` to 0 and
    ``CoolingTestResetEnergies`` to 1. 2) Allow gas to cool, allowing one to plot
    Temperature vs.  time. For this, set ``RadiativeCooling`` to 1 and
    ``CoolingTestResetEnergies`` to 0.


``CoolingTestMinimumHNumberDensity`` (external)
    The minimum density in code units at x=0. Default: 1
    [cm\ :sup:`-3`\ ].
``CoolingTestMaximumHNumberDensity`` (external)
    The maximum density in code units at
    x=``DomainRightEdge[0]``. Default: 1e6
    [cm\ :sup:`-3`\ ].
``CoolingTestMinimumMetallicity`` (external)
    The minimum metallicity at y=0. Default: 1e-6 [Z\ :sub:`sun`\ ].
``CoolingTestMaximumMetallicity`` (external)
    The maximum metallicity at
    y=``DomainRightEdge[1]``. Default: 1
    [Z\ :sub:`sun`\ ].
``CoolingTestMinimumTemperature`` (external)
    The minimum temperature in Kelvin at z=0. Default: 10.0 [K].
``CoolingTestMaximumTemperature`` (external)
    The maximum temperature in Kelvin at
    z=``DomainRightEdge[2]``. Default: 1e7 [K].
``CoolingTestResetEnergies`` (external)
    An integer flag (0 or 1) to determine whether the grid energies
    should be continually reset after every iteration of the chemistry
    solver such that the temperature remains constant as the mean
    molecular weight varies slightly. Default: 1.


.. _3dcollapse_param:

3D Collapse Test (101)
~~~~~~~~~~~~~~~~~~~~~~

``NumberOfSpheres`` (external)
``RefineAtStart``
``UseParticles``
``MediumDensity``
``MediumPressure``
``UniformVelocity``
``SphereType[]``
``SphereRadius[]``
``SphereCoreRadius[]``
``SphereDensity[]``
``SpherePressure[]``
``SphereSoundVelocity[]``
``SpherePosition[]``
``SphereVelocity[]``
``SphereAngVel[]``
``SphereTurbulence[]``
``SphereCutOff[]``
``SphereAng1[]``
``SphereAng2[]``
``SphereNumShells[]``


.. _1dcollapse_param:

1D Spherical Collapse Test (102)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``RefineAtStart`` (external)
    Boolean flag. Default: TRUE
``UseParticles`` (external)
    Boolean flag. Default: False
``MediumDensity`` (external)
    Initial density of the medium. Default: 1.0
``MediumPressure`` (external)
    Initial pressure of the medium. Default: 1.0
``SphereType`` (external)
    Default: 0
``SphereRadius`` (external)
    Radius of the sphere. Default: 1.0
``SphereCoreRadius`` (external)
    Radius of the core. Default: 0
``SphereDensity`` (external)
    Initial density of the sphere. Default: 1.0
``SpherePressure`` (external)
    Initial pressure of the sphere. Default: 1.0
``SphereSoundVelocity`` (external)
    Velocity of sound. Default: 1.0
``SphereAngVel`` (external)
    Angular velocity of the sphere. Default: 0.0

.. _mhdhydro_param:

Hydro and MHD Turbulence Simulation (106)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``RefineAtStart`` (external)
    Boolean flag. Default: TRUE
``PutSink`` (external)
    Boolean flag. Default: FALSE
``Density`` (external)
    Boolean flag. Default: TRUE
``SoundVelocity`` (external)
    Velocity of sound. Default: 1.0
``MachNumber`` (external)
    Default: 1.0
``AngularVelocity`` (external)
    Default: 0
``CloudRadius`` (external)
    Initial radius of the cloud. Default: 0.05
``SetTurbulence`` (external)
    Boolean flag. Default: TRUE
``InitialBfield`` (external)
    Initial magnetic field strength. Default: 0
``RandomSeed`` (external)
    Default: 52761
``CloudType`` (external)
    Default: 1


.. _putsink_param:

Put Sink from Restart (107)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``PutSinkRestartName`` (external)
     Filename to restart from. 


.. _clustercoolingflow_param:

Cluster Cooling Flow (108)
~~~~~~~~~~~~~~~~~~~~~~~~~~

``ClusterSMBHFeedback`` (external)
    Boolean flag. Default: FALSE
``ClusterSMBHJetMdot`` (external)
    Mdot of one Jet. Units: Solar mass per year. Default: 3.0
``ClusterSMBHJetVelocity`` (external)
    Units:km/s. Default: 10000.0
``ClusterSMBHJetRadius`` (external)
    The radius of the jet launching region. Units: cell width. Default: 6.0
``ClusterSMBHJetLaunchOffset`` (external)
    The distance of the jet launching plane to the center of the cluster. Units: cell width. Default: 10.0
``ClusterSMBHStartTime`` (external)
    The time to start feedback in code unit. Default: 1.0
``ClusterSMBHTramp`` (external)
    The ramp time in Myr. Default: 0.1
``ClusterSMBHJetOpenAngleRadius`` (external)
    Default: 0.0
``ClusterSMBHFastJetRadius`` (external)
    Default: 0.1
``ClusterSMBHFastJetVelocity`` (external)
    Unit: km/s. Default: 10000.0
``ClusterSMBHJetEdot`` (external)
    Unit: 10^44 ergs/s. Default: 1.0
``ClusterSMBHKineticFraction`` (external)
    The fraction of kinetic energy feedback; the rest is thermal feedback. Default: 1.0
``ClusterSMBHJetAngleTheta`` (external)
    The angle of the jet direction with respect to z-axis. Default: 0.0 (along the axis)
``ClusterSMBHJetAnglePhi`` (external)
    Default: 0.0
``ClusterSMBHJetPrecessionPeriod`` (external)
    Unit: Myr. Default: 0.0 (not precessing)
``ClusterSMBHCalculateGasMass`` (external)
    Type: integer. 1--Calculate the amount of cold gas around the SMBH and remove it at the rate of 2*Mdot; 2--Calculate Mdot based on the amount of cold gas around the SMBH; 0--off (do not remove cold gas). Default: 1.
``ClusterSMBHFeedbackSwitch`` (external)
    Boolean flag. When ClusterSMBHCalculateGasMass=1, ClusterSMBHFeedbackSwitch is turned on when there is enough cold gas (ClusterSMBHEnoughColdGas) around the SMBH. Default: FALSE
``ClusterSMBHEnoughColdGas`` (external)
    Unit: Solar mass. Default: 1.0e7
``ClusterSMBHAccretionTime`` (external)
    When ClusterSMBHCalculateGasMass = 2, Mdot = Mcold/ClusterSMBHAccretionTime. Default: 5.0 (Myr)
``ClusterSMBHJetDim`` (external)
    0--x; 1--y; 2--z. Default: 2
``ClusterSMBHAccretionEpsilon`` (external)
    Jet Edot = ClusterSMBHAccretionEpsilon * Mdot * c^2. Default: 0.001


.. _mhd1d_param:

1D MHD Test (200)
~~~~~~~~~~~~~~~~~

``RefineAtStart`` (external)
    Boolean flag. Default: TRUE
``LeftVelocityX``, ``RightVelocityX`` (external)
    Initial velocity x-direction. Default: 0 (for both)
``LeftVelocityY``, ``RightVelocityY`` (external)
    Initial velocity y-direction. Default: 0 (for both)
``LeftVelocityZ``, ``RightVelocityZ`` (external)
    Initial velocity z-direction. Default: 0 (for both)
``LeftPressure``, ``RightPressure`` (external)
    Initial pressure. Default: 1.0 (for both)
``LeftDensity``, ``RightDensity`` (external)
    Initial density. Default: 1.0 (for both)
``LeftBx``, ``RightBx`` (external)
    Initial magnetic field x-direction. Default: 0 (for both)
``LeftBy``, ``RightBy`` (external)
    Initial magnetic field y-direction. Default: 0 (for both)
``LeftBz``, ``RightBz``  (external)
    Initial magnetic field z-direction. Default: 0 (for both)

.. _mhd2d_param:

2D MHD Test (201)
~~~~~~~~~~~~~~~~~

``RefineAtStart`` (external)
    Boolean flag. Default: TRUE
``LowerVelocityX``, ``UpperVelocityX`` (external)
    Initial velocity x-direction. Default: 0 (for both)
``LowerVelocityY``, ``UpperVelocityY`` (external)
    Initial velocity y-direction. Default: 0 (for both)
``LowerPressure``, ``UpperPressure`` (external)
    Initial pressure. Default: 1.0 (for both)
``LowerDensity``, ``UpperDensity`` (external)
    Initial density. Default: 1.0 (for both)
``LowerBx``, ``UpperBx`` (external)
    Initial magnetic field x-direction. Default: 0 (for both)
``LowerBy``, ``UpperBy`` (external)
    Initial magnetic field y-direction. Default: 0 (for both)
``MHD2DProblemType`` (external)
    Default: 0
``RampWidth`` (external)
    Default: 0.05
``UserColour`` (external)
    Boolean flag. Default: FALSE

.. _mhd3d_param:

3D MHD Collapse Test (202)
~~~~~~~~~~~~~~~~~~~~~~~~~~


``RefineAtStart`` (external)
    Boolean flag. Default: FALSE
``LowerVelocityX``, ``UpperVelocityX`` (external)
    Initial velocity x-direction. Default: 0 (for both)
``LowerVelocityY``, ``UpperVelocityY`` (external)
    Initial velocity y-direction. Default: 0 (for both)
``LowerPressure``, ``UpperPressure`` (external)
    Initial pressure. Default: 1.0 (for both)
``LowerDensity``, ``UpperDensity`` (external)
    Initial density. Default: 1.0 (for both)
``LowerBx``, ``UpperBx`` (external)
    Initial magnetic field x-direction. Default: 0 (for both)
``LowerBy``, ``UpperBy`` (external)
    Initial magnetic field y-direction. Default: 0 (for both)
``MHD3DProblemType`` (external)
    Default: 0

.. _mhdtcollapse_param:

MHD Turbulent Collapse Test (203)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``RefineAtStart`` (external)
    Boolean flag. Default: TRUE
``Density`` (external)
    Initial density. Default: 1.0
``SoundVelocity`` (external)
    Speed of sound. Default: 1.0
``MachNumber`` (external)
    Default: 1.0
``InitialBfield`` (external)
    Initial magnetic field strength. Default: 0
``RandomSeed`` (external)
    Default: 0


.. _galaxydisk_param:

Galaxy Disk (207)
~~~~~~~~~~~~~~~~~

``NumberOfHalos`` (external)
    Number of Halos simulated. Default: 1
``RefineAtStart`` (external)
    Boolean flag. Default: TRUE
``UseParticles`` (external)
    Boolean flag. Default: FALSE
``UseGas`` (external)
    Boolean flag. Default: TRUE
``MediumTemperature`` (external)
    Temperature of the medium. Default: 1000
``MediumDensity`` (external)
    Density of the medium. Default: 1.0
``HaloMagneticField`` (external)
    Magnetic Field Strength. Default: 0
``UniformVelocity[i]`` (external)
    Velocity in all 3 dimensions. Default: 0 (all)
``GalaxyType[i]`` (external)
    Sppecifying galaxy type for the ith sphere. Default: 0 (all)
``HaloRadius[i]`` (external)
    Radius of the halo for the ith sphere. Default: 1 (all)
``HaloCoreRadius[i]`` (external)
    Core radius for the ith sphere. Default: 0.1 (all) 
``HaloDensity[i]`` (external)
    Density of the halo for the ith sphere. Default: 1 (all)
``HaloTemperature[i]`` (external)
    Temperature of the halo for the ith sphere. Default: 1 (all)
``HaloAngVel[i]`` (external)
    TBD
``HaloSpin[i]`` (external)
    TBD
``HaloPosition[i][j]`` (external)
    Position of the Halo. 
``HaloVelocity[i][j]`` (external)
    Velocity of the Halo.
``DiskRadius[i]`` (external)
    TBD
``DiskHeight[i]`` (external)
    TBD
``DiskDensity[i]`` (external)
    TBD
``DiskTemperature[i]`` (external)
    TBD
``DiskMassFraction[i]`` (external)
    Default: 0 (all)
``DiskFlaringParameter[i]`` (external)
    Default: 10 (all)

.. _agndisk_param:

AGN Disk (207)
~~~~~~~~~~~~~~

``DiskType`` (external)
    Default: 1
``RefineAtStart`` (external)
    Boolean flag. Default: 0
``BlackHoleMass`` (external)
    Initial mass of black hole. Default: 0
``UseGas`` (external)
    Boolean flag. Default: 1
``DiskDensity`` (external)
    Initial density of the disk. Default: 1
``DiskTemperature`` (external)
    Initial temperature of the disk. Default: 1
``DiskRadius`` (external)
    Initial radius of the disk. Default: 1
``DiskHeight`` (external)
    Initial height of the disk. Default: 1

.. _poissonsolver_param:

Poisson Solver Test (300)
~~~~~~~~~~~~~~~~~~~~~~~~~


``PoissonSolverTestType`` (external)
   Default: 0
``PoissonSolverTestGeometryControl`` (external)
   Default: 1
``PoissonSolverTestRefineAtStart`` (external)
   Boolean flag. Default: 0

.. _rhdtest1_param:

Radiation-Hydrodynamics Test 1 - Constant Fields (400)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Basic FLD radiation problem initializer, allowing setup of uniform
    fields throughout the computational domain, which are useful for
    testing radiation/material couplings. Test problem used for
    problem 4.2 in (Reynolds et al., "Self-consistent solution of
    cosmological radiation-hydrodynamics and chemical ionization,"
    JCP, 2009).

``RadHydroVelocity`` (external)
   Initialize velocity of ambient gas in the x,y,z directions. Default: 0 (all).
   Example RadHydroVelocity = 0.1 0.1 0.1
``RadHydroChemistry`` (external)
   Number of chemical species.  1 implies hydrogen only, 3 implies
   hydrogen and helium. Default: 1.
``RadHydroModel`` (external)
   Type of radiation/matter coupling: 1 implies a standard
   chemistry-dependent model, 4 implies an isothermal
   chemistry-dependent model, 10 implies a chemistry-independent model
   in thermodynamic equilibrium. Default: 1
``RadHydroDensity`` (external)
   Ambient density. Default: 10
``RadHydroTemperature`` (external)
   Ambient temperature. Default: 1
``RadHydroIEnergy`` (external)
   Ambient internal energy (replaces temperature, if specified).  
   Default: -1
``RadHydroRadiationEnergy`` (external)
   Ambient radiation energy. Default: 10
``RadHydroInitialFractionHII`` (external)
   Initial fraction of ionized hydrogen (in relation to all hydrogen). 
   Default: 0
``RadHydroHFraction`` (external)
   Initial fraction of hydrogen (in relation to the total density).
   Default: 1
``RadHydroInitialFractionHeII`` (external)
   Initial fraction of helium II (in relation to the total helium).
   Default: 0
``RadHydroInitialFractionHeIII`` (external)
   Initial fraction of helium III (in relation to the total helium).
   Default: 0

.. _rhdtest2_param:

Radiation-Hydrodynamics Test 2 - Streams (401)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Streaming radiation tests.  The problem utilizes a uniform density
    and a constant opacity, setting one face of the domain to have a
    radiation energy density of 1.  The radiation front propagates
    through the domain at the speed of light.  The sharpness of the
    radiation front is determined by the spatial resolution.  Test
    problem used for problem 4.1 in (Reynolds et al.,
    "Self-consistent solution of cosmological radiation-hydrodynamics
    and chemical ionization," JCP, 2009).

``RadHydroDensity`` (external)
   Ambient density. Default: 1.0
``RadHydroRadEnergy`` (external)
   Ambient radiation energy. Default 1.0e-10
``RadStreamDim`` (external)
   Dimension to test {0,1,2}. Default: 0
``RadStreamDir`` (external)
   Direction for streaming radiation. 0 for left to right. 1 for right to left.
   Default: 0

.. _rhdtest3_param:

Radiation-Hydrodynamics Test 3 - Pulse (402)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``RadHydroDensity`` (external)
   Ambient density. Default: 1.0
``RadHydroRadEnergy`` (external)
   Ambient radiation energy. Default 1.0e-10
``RadPulseDim`` (external)
   Dimension to test {0,1,2}. Default: 0

.. _rhdtest4_param:

Radiation-Hydrodynamics Test 4 - Grey Marshak Test (403)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Test problem used for problem 4.3 in (Reynolds et al.,
    "Self-consistent solution of cosmological radiation-hydrodynamics
    and chemical ionization," JCP, 2009).

``RadHydroDensity`` (external)
   Ambient density. Default: 1.0
``RadHydroRadEnergy`` (external)
   Ambient radiation energy. Default 1.0
``RadHydroGasEnergy`` (external)
   Ambient gas energy. Default: 1.0
``GreyMarshDir`` (external)
   Propagation coordinate for Marshak problem. {0,1,2}. Default: 0

.. _rhdtest5_param:

Radiation-Hydrodynamics Test 5 - Radiating Shock (404/405)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Test problem used for problem 4.4 in (Reynolds et al.,
    "Self-consistent solution of cosmological radiation-hydrodynamics
    and chemical ionization," JCP, 2009).

``DensityConstant`` (external)
   Ambient density. Default: 1.0
``GasTempConstant`` (external)
   Ambient gas temperature. Default: 1.0
``RadTempConstant`` (external)
   Ambient radiation temperature. Default: 1.0
``VelocityConstant`` (external)
   Imposed fluid velocity. Default: 1.0
``ShockDir`` (external)
   Propagation coordinate for shock. {0,1,2}. Default: 0
``CGSType`` (external)
   1 = Astrophysical Setup Parameters; 
   2 = "lab" setup parameters, after Lowrie; 
   Default: 1

.. _rhdtest10_param:

Radiation-Hydrodynamics Tests 10 and 11 - I-Front Tests (410/411)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Uniform density ionization front test problems.  These tests are
    used to replicate the isothermal and temperature-dependent I-front
    tests 1 and 2 from (Iliev et al., "Cosmological Radiative Transfer
    Codes Comparison Project I: The Static Density Field Tests,"
    MNRAS, 2006).  This test problem was used for problem 4.5 in
    (Reynolds et al., "Self-consistent solution of cosmological
    radiation-hydrodynamics and chemical ionization," JCP, 2009).

``RadHydroVelocity`` (external)
   Initial velocity of ambient gas in the x,y,z directions. Default: 0 (all). 
   Example RadHydroVelocity = 0.1 0.1 0.1
``RadHydroChemistry`` (external)
   Number of chemical species.  1 implies hydrogen only, 3 implies
   hydrogen and helium. Default: 1.
``RadHydroModel`` (external)
   Type of radiation/matter coupling: 1 implies a standard
   chemistry-dependent model, 4 implies an isothermal
   chemistry-dependent model. Default: 1
``RadHydroDensity`` (external)
   Ambient density. Default: 10
``RadHydroTemperature`` (external)
   Ambient temperature. Default: 1
``RadHydroIEnergy`` (external)
   Ambient internal energy (replaces temperature, if specified).  
   Default: -1
``RadHydroRadiationEnergy`` (external)
   Ambient radiation energy. Default: 10
``RadHydroInitialFractionHII`` (external)
   Initial fraction of ionized hydrogen (in relation to all hydrogen). 
   Default: 0
``RadHydroHFraction`` (external)
   Initial fraction of hydrogen (in relation to the total density).
   Default: 1 
``RadHydroInitialFractionHeII`` (external)
   Initial fraction of helium II (in relation to the total helium).
   Default: 0 
``RadHydroInitialFractionHeIII`` (external)
   Initial fraction of helium III (in relation to the total helium).
   Default: 0
``NGammaDot`` (external)
   Strength of ionization source, in number of photons per second.
   Default: 0
``EtaRadius`` (external)
   Radius of ionization source, in cells (0 implies a single-cell source).
   Default: 0
``EtaCenter`` (external)
   Location of ionization source, in scaled length units, in the x,y,z
   directions. Default: 0 (all). 
   Example EtaCenter = 0.5 0.5 0.5

.. _rhdtest12_param:

Radiation-Hydrodynamics Test 12 - HI ionization of a clump (412)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Ionization of a hydrogen clump, used to investigate I-front
    trapping in a dense clump, and the formation of a shadow.  This
    test replicates the test 3.4 from (Iliev et al., "Cosmological
    Radiative Transfer Codes Comparison Project I: The Static Density
    Field Tests," MNRAS, 2006).

``RadHydroVelocity`` (external)
   Initial velocity of ambient gas in the x,y,z directions. Default: 0 (all).
   Example RadHydroVelocity = 0.1 0.1 0.1
``RadHydroChemistry`` (external)
   Number of chemical species.  1 implies hydrogen only, 3 implies
   hydrogen and helium. Default: 1.
``RadHydroModel`` (external)
   Type of radiation/matter coupling: 1 implies a standard
   chemistry-dependent model, 4 implies an isothermal
   chemistry-dependent model. Default: 1
``RadHydroNumDensityIn`` (external)
   Number density inside the clump. Default: 0.04
``RadHydroNumDensityOut`` (external)
   Number density outside the clump. Default: 0.0002
``RadHydroTemperatureIn`` (external)
   Temperature inside the clump. Default: 40
``RadHydroTemperatureOut`` (external)
   Temperature outside the clump. Default: 8000
``RadHydroRadiationEnergy`` (external)
   Ambient radiation energy. Default: 10
``RadHydroInitialFractionHII`` (external)
   Initial fraction of ionized hydrogen (in relation to all hydrogen). 
   Default: 0
``ClumpCenter`` (external)
   Location of clump center, in cm, in the x,y,z directions. 
   Default: 1.54285e22 1.018281e22 1.018281e22
``ClumpRadius`` (external)
   Radius of clump, in cm.
   Default: 2.46856e21
``NGammaDot`` (external)
   Strength of ionization source along left wall, in number of photons
   per second.  Default: 0

.. _rhdtest13_param:

Radiation-Hydrodynamics Test 13 - HI ionization of a steep region (413)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Ionization of a steep density gradient, used to investigate HII
    region expansion along a 1/r^2 density profile.  This test
    replicates the test 3.2 from (Iliev et al., "Cosmological
    Radiative Transfer Comparison Project II: The
    Radiation-Hydrodynamic Tests," MNRAS, 2009).

``RadHydroVelocity`` (external)
   Initial velocity of ambient gas in the x,y,z directions. Default: 0 (all).
   Example RadHydroVelocity = 0.1 0.1 0.1
``RadHydroChemistry`` (external)
   Number of chemical species.  1 implies hydrogen only, 3 implies
   hydrogen and helium. Default: 1.
``RadHydroModel`` (external)
   Type of radiation/matter coupling: 1 implies a standard
   chemistry-dependent model, 4 implies an isothermal
   chemistry-dependent model. Default: 1
``RadHydroNumDensity`` (external)
   Number density inside the core of the dense region. Default: 3.2
``RadHydroDensityRadius`` (external)
   Radius of the dense region, in cm. Default: 2.8234155e+20
``RadHydroTemperature`` (external)
   Ambient temperature. Default: 100
``RadHydroRadiationEnergy`` (external)
   Ambient radiation energy. Default: 1e-20
``RadHydroInitialFractionHII`` (external)
   Initial fraction of ionized hydrogen (in relation to all hydrogen). 
   Default: 0
``EtaCenter`` (external)
   Center of the dense region (and ionization source), in cm, in the
   x,y,z directions.  Default: 0 0 0
``NGammaDot`` (external)
   Strength of ionization source, in number of photons per second.
   Default: 0

.. _rhdtest14_param:

Radiation-Hydrodynamics Tests 14/15 - Cosmological HI ionization (414/415)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HI ionization in a uniform density field.  This test problem was
    used for problems 4.6 and 4.8 in (Reynolds et al.,
    "Self-consistent solution of cosmological radiation-hydrodynamics
    and chemical ionization," JCP, 2009).  Test 4.6 utilized a single
    ionization source (test 415), whereas 4.8 replicated the test to
    the center of every processor for performing weak-scaling tests
    (test 414).

``RadHydroVelocity`` (external)
   Initial velocity of ambient gas in the x,y,z directions. Default: 0 (all).
   Example RadHydroVelocity = 0.1 0.1 0.1
``RadHydroChemistry`` (external)
   Number of chemical species.  1 implies hydrogen only, 3 implies
   hydrogen and helium. Default: 1.
``RadHydroModel`` (external)
   Type of radiation/matter coupling: 1 implies a standard
   chemistry-dependent model, 4 implies an isothermal
   chemistry-dependent model. Default: 1
``RadHydroTemperature`` (external)
   Ambient temperature in K. Default: 10000
``RadHydroRadiationEnergy`` (external)
   Ambient radiation energy in erg/cm^3. Default: 1.0e-32
``RadHydroInitialFractionHII`` (external)
   Initial fraction of ionized hydrogen (in relation to all hydrogen). 
   Default: 0
``RadHydroOmegaBaryonNow`` (external)
   Default: 0.2
``NGammaDot`` (external)
   Strength of ionization source, in number of photons per second.
   Default: 0
``EtaRadius`` (external)
   Radius of ionization source for test 415, in cells (0 implies a
   single-cell source). 
   Default: 0
``EtaCenter`` (external)
   Location of ionization source for test 415, in scaled length units,
   in the x,y,z directions. Default: 0 (all).
   Example EtaCenter = 0.5 0.5 0.5

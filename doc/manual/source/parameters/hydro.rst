Hydrodynamics Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

General
^^^^^^^

``UseHydro`` (external)
    This flag (1 - on, 0 - off) controls whether a hydro solver is used.  
    Default: 1
``HydroMethod`` (external)
    This integer specifies the hydrodynamics method that will be used.
    Currently implemented are

    ============ =============================
    Hydro method Description
    ============ =============================
    0            PPM DE (a direct-Eulerian version of PPM)
    1            [reserved]
    2            ZEUS (a Cartesian, 3D version of Stone & Norman). Note that if ZEUS is selected, it automatically turns off ``ConservativeInterpolation`` and the ``DualEnergyFormalism`` flags.
    3            Runge Kutta second-order based MUSCL solvers.
    4            Same as 3 but including Dedner MHD (Wang & Abel 2008). For 3 and 4 there are the additional parameters ``RiemannSolver`` and ``ReconstructionMethod`` you want to set.
    ============ =============================

    Default: 0

    More details on each of the above methods can be found at :ref:`hydro_methods`.
``RiemannSolver`` (external; only if ``HydroMethod`` is 3 or 4)
    This integer specifies the Riemann solver used by the MUSCL solver. Choice of

    ============== ===========================
    Riemann solver Description
    ============== ===========================
    0              [reserved]
    1              HLL (Harten-Lax-van Leer) a two-wave, three-state solver with no resolution of contact waves
    2              [reserved]
    3              LLF (Local Lax-Friedrichs)
    4              HLLC (Harten-Lax-van Leer with Contact) a three-wave, four-state solver with better resolution of contacts
    5              TwoShock
    ============== ===========================

    Default: 1 (HLL) for ``HydroMethod`` = 3; 5 (TwoShock) for
    ``HydroMethod`` = 0

``RiemannSolverFallback`` (external)
    If the euler update results in a negative density or energy, the
    solver will fallback to the HLL Riemann solver that is more
    diffusive only for the failing cell.  Only active when using the
    HLLC or TwoShock Riemann solver.  Default: OFF.
``ReconstructionMethod`` (external; only if ``HydroMethod`` is 3 or 4)
    This integer specifies the reconstruction method for the MUSCL solver. Choice of

    ===================== ====================
    Reconstruction Method Description
    ===================== ====================
    0                     PLM (piecewise linear)
    1                     PPM (piecwise parabolic)
    2                     [reserved]
    3                     [reserved]
    4                     [reserved]
    ===================== ====================

    Default: 0 (PLM) for ``HydroMethod`` = 3; 1 (PPM) for ``HydroMethod`` = 0

``Gamma`` (external)
    The ratio of specific heats for an ideal gas (used by all hydro
    methods). If using multiple species (i.e. ``MultiSpecies`` > 0), then
    this value is ignored in favor of a direct calculation (except for
    PPM LR) Default: 5/3.
``Mu`` (external)
    The molecular weight. Default: 0.6.
``ConservativeReconstruction`` (external)
    Experimental.  This option turns on the reconstruction of the
    left/right interfaces in the Riemann problem in the conserved
    variables (density, momentum, and energy) instead of the primitive
    variables (density, velocity, and pressure).  This generally gives
    better results in constant-mesh problems has been problematic in
    AMR simulations.  Default: OFF
``PositiveReconstruction`` (external)
    Experimental and not working.  This forces the Riemann solver to
    restrict the fluxes to always give positive pressure.  Attempts to
    use the Waagan (2009), JCP, 228, 8609 method.  Default: OFF
``CourantSafetyNumber`` (external)
    This is the maximum fraction of the CFL-implied timestep that will
    be used to advance any grid. A value greater than 1 is unstable
    (for all explicit methods). The recommended value is 0.4. Default:
    0.6.
``RootGridCourantSafetyNumber`` (external)
    This is the maximum fraction of the CFL-implied timestep that will
    be used to advance ONLY the root grid. When using simulations with
    star particle creation turned on, this should be set to a value of
    approximately 0.01-0.02 to keep star particles from flying all over
    the place. Otherwise, this does not need to be set, and in any case
    should never be set to a value greater than 1.0. Default: 1.0.
``DualEnergyFormalism`` (external)
    The dual energy formalism is needed to make total energy schemes
    such as PPM DE and PPM LR stable and accurate in the
    "hyper-Machian" regime (i.e. where the ratio of thermal energy to
    total energy < ~0.001). Turn on for cosmology runs with PPM DE and
    PPM LR. Automatically turned off when used with the hydro method
    ZEUS. Integer flag (0 - off, 1 - on). When turned on, there are two
    energy fields: total energy and thermal energy. Default: 0
``DualEnergyFormalismEta1``, ``DualEnergyFormalismEta2`` (external)
    These two parameters are part of the dual energy formalism and
    should probably not be changed. Defaults: 0.001 and 0.1
    respectively.
``PressureFree`` (external)
    A flag that is interpreted by the PPM DE hydro method as an
    indicator that it should try and mimic a pressure-free fluid. A
    flag: 1 is on, 0 is off. Default: 0
``PPMFlatteningParameter`` (external)
    This is a PPM parameter to control noise for slowly-moving shocks.
    It is either on (1) or off (0). Default: 0
``PPMDiffusionParameter`` (external)
    This is the PPM diffusion parameter (see the Colella and Woodward
    method paper for more details). It is either on (1) or off (0).
    Default: 1 [Currently disabled (set to 0)]
``PPMSteepeningParameter`` (external)
    A PPM modification designed to sharpen contact discontinuities. It
    is either on (1) or off (0). Default: 0
``ZEUSQuadraticArtificialViscosity`` (external)
    This is the quadratic artificial viscosity parameter C2 of Stone &
    Norman, and corresponds (roughly) to the number of zones over which
    a shock is spread. Default: 2.0
``ZEUSLinearArtificialViscosity`` (external)
    This is the linear artificial viscosity parameter C1 of Stone &
    Norman. Default: 0.0

Minimum Pressure Support Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``UseMinimumPressureSupport`` (external)
    When radiative cooling is turned on, and objects are allowed to
    collapse to very small sizes so that their Jeans length is no
    longer resolved, then they may undergo artificial fragmentation
    and angular momentum non-conservation.  To alleviate this problem,
    as discussed in more detail in Machacek, Bryan & Abel (2001), a
    very simple fudge was introduced: if this flag is turned on, then
    a minimum temperature is applied to grids with level ==
    ``MaximumRefinementLevel``. This minimum temperature is that
    required to make each cell Jeans stable multiplied by the
    parameter below.  More precisely, the temperature of a cell is set
    such that the resulting Jeans length is the square-root of the
    parameter ``MinimumPressureSupportParameter``.  So, for the
    default value of 100 (see below), this insures that the ratio of
    the Jeans length/cell size is at least 10.  Default: 0
``MinimumPressureSupportParameter`` (external)
    This is the numerical parameter discussed above. Default: 100

Magnetohydrodynamics (CT) Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Coming very soon...!

Magnetohydrodynamics (Dedner) Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``UseDivergenceCleaning`` (external)
    Method 1 and 2 are a failed experiment to do divergence cleaning
    using successive over relaxation. Method 3 uses conjugate gradient
    with a 2 cell stencil and Method 4 uses a 4 cell stencil. 4 is more
    accurate but can lead to aliasing effects. Default: 0
``DivergenceCleaningBoundaryBuffer`` (external)
    Choose to *not* correct in the active zone of a grid by a
    boundary of cells this thick. Default: 0
``DivergenceCleaningThreshold`` (external)
    Calls divergence cleaning on a grid when magnetic field divergence
    is above this threshold. Default: 0.001
``PoissonApproximateThreshold`` (external)
    Controls the accuracy of the resulting solution for divergence
    cleaning Poisson solver. Default: 0.001
``UseDrivingField`` (external)
    This parameter is used to add external driving force as a source term in some test problems; see hydro_rk/Grid_(MHD)SourceTerms.C. Default: 0
``DrivingEfficiency`` (external)
    This parameter is used to define the efficiency of such driving force; see hydro_rk/Grid_(MHD)SourceTerms.C. Default: 1.0
``UseConstantAcceleration`` (external)
    This parameter is used to add constant acceleration as a source term in some set-ups; see hydro_rk/Grid_(MHD)SourceTerms.C. Default: 0
``ConstantAcceleration[]`` (external)
    This parameter is used to define the value of such acceleration; see hydro_rk/Grid_(MHD)SourceTerms.C. 
``UseViscosity`` (external)
    This parameter is used to add viscosity and thereby update velocity in some set-ups (1 - constant viscosity, 2 - alpha viscosity); see ComputeViscosity in hydro_rk/Grid_AddViscosity.C.  Default: 0
``ViscosityCoefficient`` (external)
    This parameter is used to define the value of such viscosity for UseViscosity = 1; see ComputeViscosity in hydro_rk/Grid_AddViscosity.C. Default: 0.0
``UseGasDrag`` (external)
    This parameter is used to calculate velocity decrease caused by gas drag as a source term in some set-ups; see hydro_rk/Grid_(MHD)SourceTerms.C. Default: 0
``GasDragCoefficient`` (external)
    This parameter is used to define the value of such gas drag; see hydro_rk/Grid_(MHD)SourceTerms.C. Default: 0.0
``UseFloor`` (external)
    This parameter is used to impose the minimum energy based on MaximumAlvenSpeed in some set-ups; see hydro_rk/Grid_SetFloor.C. Default: 0
``MaximumAlvenSpeed`` (external)
    This parameter is used to define the value of such minimum; see hydro_rk/Grid_SetFloor.C. Default: 1e30
``UseAmbipolarDiffusion`` (external)
    This parameter is used to update magnetic fields by ambipolar diffusion in some set-ups; see hydro_rk/Grid_AddAmbipolarDiffusion.C. Default: 0
``UseResistivity`` (external)
    This parameter is used to add resistivity and thereby update magnetic fields in some set-ups; see ComputeResistivity in hydro_rk/Grid_AddResistivity.C.  Default: 0
``UsePhysicalUnit`` (external)
    For some test problems (mostly in hydro_rk), the relevant parameters could be defined in physical CGS units.  Default: 0
``SmallRho`` (external)
    Minimum value for density in hydro_rk/EvolveLevel_RK.C.  Default: 1e-30 (note that the default value assumes UsePhysicalUnit = 1)
``SmallT`` (external)
    Minimum value for temperature in hydro_rk/EvolveLevel_RK.C.  Default: 1e-10 (note that the default value assumes UsePhysicalUnit = 1)
``SmallP``
    [not used]
``RKOrder``
    [not used]
``Theta_Limiter`` (external)
    Flux limiter in the minmod Van Leer formulation.  Must be between 1 (most dissipative) and 2 (least dissipative). Default: 1.5
``Coordinate`` (external)
    Coordinate systems to be used in hydro_rk/EvolveLevel_RK.C.  Currently implemented are Cartesian and Spherical for HD_RK, and Cartesian and Cylindrical for MHD_RK.  See Grid_(MHD)SourceTerms.C.  Default: Cartesian
``EOSType`` (external)
    Types of Equation of State used in hydro_rk/EvolveLevel_RK.C (0 - ideal gas, 1 - polytropic EOS, 2 - another polytropic EOS, 3 - isothermal, 4 - pseudo cooling, 5 - another pseudo cooling, 6 - minimum pressure); see hydro_rk/EOS.h. Default: 0
``EOSSoundSpeed`` (external)
    Sound speed to be used in EOS.h for EOSType = 1, 2, 3, 4, 5.  Default: 2.65e4
``EOSCriticalDensity`` (external)
    Critical density to be used in EOS.h for EOSType = 1, 2, 4, 6. Default: 1e-13
``EOSGamma`` (external)
    Polytropic gamma to be used in EOS.h for EOSType = 1. Default: 1.667
``DivBDampingLength`` (external)
    From C_h (the Dedner wave speeds at which the div*B error is isotropically transferred; as defined in e.g. Matsumoto, PASJ, 2007, 59, 905) and this parameter, C_p (the decay rate of the wave) is calculated; see ComputeDednerWaveSpeeds.C  Default: 1.0
``UseCUDA`` (external)
    Set to 1 to use the CUDA-accelerated (M)HD solver.  Only works if compiled with cuda-yes. Default: 0
``ResetMagneticField`` (external)
    Set to 1 to reset the magnetic field in the regions that are denser
    than the critical matter density. Very handy when you want to
    re-simulate or restart the dumps with MHD. Default: 0
``ResetMagneticFieldAmplitude`` (external)
    The magnetic field values (in Gauss) that will be used for the
    above parameter. Default: 0.0 0.0 0.0
``CoolingCutOffDensity1``
    Reserved for future use
``CoolingCutOffDensity2``
    Reserved for future use
``CoolingCutOffTemperature``
    Reserved for future use
``CoolingPowerCutOffDensity1``
    Reserved for future use
``CoolingPowerCutOffDensity2``
    Reserved for future use


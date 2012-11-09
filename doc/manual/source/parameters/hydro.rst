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

Coming soon...

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
``ResetMagneticField`` (external)
    Set to 1 to reset the magnetic field in the regions that are denser
    than the critical matter density. Very handy when you want to
    re-simulate or restart the dumps with MHD. Default: 0
``ResetMagneticFieldAmplitude`` (external)
    The magnetic field values (in Gauss) that will be used for the
    above parameter. Default: 0.0 0.0 0.0


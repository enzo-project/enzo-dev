.. _hydrodynamics_parameters:

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

    ============== =============================
    Hydro method   Description
    ============== =============================
    0              PPM DE (a direct-Eulerian version of PPM)
    1              [reserved]
    2              ZEUS (a Cartesian, 3D version of Stone & Norman). Note that if ZEUS is selected, it automatically turns off ``ConservativeInterpolation`` and the ``DualEnergyFormalism`` flags.
    3              Runge Kutta second-order based MUSCL solvers.
    4              Same as 3 but including Dedner MHD (Wang & Abel 2008). For 3 and 4 there are the additional parameters ``RiemannSolver`` and ``ReconstructionMethod`` you want to set.
    5              No Hydro (Testing only)
    6              MHD with Constrained Transport.
    ============== =============================

    Default: 0

    More details on each of the above methods can be found at :ref:`hydro_methods`.
``FluxCorrection`` (external)
    This flag indicates if the flux fix-up step should be carried out
    around the boundaries of the sub-grid to preserve conservation (0 -
    off, 1 - on, 2 - direct correction for color fields). Strictly speaking
    this should always be used, but we have found it to lead to a less
    accurate solution for cosmological simulations because of the relatively
    sharp density gradients involved. However, it does appear to be
    important when radiative cooling is turned on and very dense structures
    are created. It does work with the ZEUS hydro method, but since velocity
    is face-centered, momentum flux is not corrected. If FluxCorrection = 1,
    species quantities are not flux corrected directly but are modified to
    keep the fraction constant based on the density change. If FluxCorrection
    = 2, species quantities are flux corrected directly in the same way as
    density and energy. Default: 1
``InterpolationMethod`` (external)
    There should be a whole section devoted to the interpolation
    method, which is used to generate new sub-grids and to fill in the
    boundary zones of old sub-grids, but a brief summary must suffice.
    The possible values of this integer flag are shown in the table
    below. The names specify (in at least a rough sense) the order of
    the leading error term for a spatial Taylor expansion, as well as a
    letter for possible variants within that order. The basic problem
    is that you would like your interpolation method to be:
    multi-dimensional, accurate, monotonic and conservative. There
    doesn't appear to be much literature on this, so I've had to
    experiment. The first one (ThirdOrderA) is time-consuming and
    probably not all that accurate. The second one (SecondOrderA) is
    the workhorse: it's only problem is that it is not always
    symmetric. The next one (SecondOrderB) is a failed experiment, and
    SecondOrderC is not conservative. FirstOrderA is everything except
    for accurate. If HydroMethod = 2 (ZEUS), this flag is ignored, and
    the code automatically uses SecondOrderC for velocities and
    FirstOrderA for cell-centered quantities. Default: 1
    ::

              0 - ThirdOrderA     3 - SecondOrderC
              1 - SecondOrderA    4 - FirstOrderA
              2 - SecondOrderB  

``ConservativeInterpolation`` (external)
    This flag (1 - on, 0 - off) indicates if the interpolation should
    be done in the conserved quantities (e.g. momentum rather than
    velocity). Ideally, this should be done, but it can cause problems
    when strong density gradients occur. This must(!) be set off for
    ZEUS hydro (the code does it automatically). Default: 1
``RiemannSolver`` (external)
    This integer specifies the Riemann solver. Solver options, and the relevant
    hydro method, are summarized as follows:

    ================ =========== ===========================
    Riemann solver   HydroMethod Description
    ================ =========== ===========================
    0                --          [reserved]
    1                0,3,4       HLL (Harten-Lax-van Leer) a two-wave, three-state solver with no resolution of contact waves
    2                            [reserved]
    3                3,4         LLF (Local Lax-Friedrichs)
    4                0,3         HLLC (Harten-Lax-van Leer with Contact) a three-wave, four-state solver with better resolution of contacts
    5                0           TwoShock 
    6                4,6         HLLD 
    ================ =========== ===========================

    Default: 1 (HLL) for ``HydroMethod`` = 3; 5 (TwoShock) for
    ``HydroMethod`` = 0; 6 (HLLD) for ``HydroMethod = 6``
``RiemannSolverFallback`` (external; only if ``HydroMethod`` is 0, 3 or 4)
    If the euler update results in a negative density or energy, the
    solver will fallback to the HLL Riemann solver that is more
    diffusive only for the failing cell.  Only active when using the
    HLLC or TwoShock Riemann solver.  Default: OFF.
``ReconstructionMethod`` (external; only if ``HydroMethod`` is 3 or 4)
    This integer specifies the reconstruction method for the MUSCL solver. Choice of

    ===================== ============ ===================
    Reconstruction Method HydroMethod  Description
    ===================== ============ ===================
    0                     0,3,4,6      PLM (piecewise linear) 
    1                     0            PPM (piecwise parabolic)
    2                                  [reserved]
    3                                  [reserved]
    4                                  [reserved]
    6                     6            MUSCL-Hancock (Non Runge-Kutta) 
    ===================== ============ ===================

    Default: 0 (PLM) for ``HydroMethod`` = 3; 1 (PPM) for ``HydroMethod`` = 0
``ConservativeReconstruction`` (external; only if ``HydroMethod`` is 3 or 4)
    Experimental.  This option turns on the reconstruction of the
    left/right interfaces in the Riemann problem in the conserved
    variables (density, momentum, and energy) instead of the primitive
    variables (density, velocity, and pressure).  This generally gives
    better results in constant-mesh problems has been problematic in
    AMR simulations.  Default: OFF
``PositiveReconstruction`` (external; only if ``HydroMethod`` is 3 or 4)
    Experimental and not working.  This forces the Riemann solver to
    restrict the fluxes to always give positive pressure.  Attempts to
    use the Waagan (2009), JCP, 228, 8609 method.  Default: OFF
``Gamma`` (external)
    The ratio of specific heats for an ideal gas (used by all hydro
    methods). If using multiple species (i.e. ``MultiSpecies`` > 0), then
    this value is ignored in favor of a direct calculation (except for
    PPM LR) Default: 5/3.
``Mu`` (external)
    The molecular weight. Default: 0.6.
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
``UseCoolingTimestep`` (external)
    This parameter will limit the timestep on each level by some fraction
    of the minimum cooling time on the level, where this fraction is
    set by ``CoolingTimestepSafetyFactor``.  In most cases, this will
    substantially decrease the timesteps, depending on the local
    cooling time, and thus increase the run time of any
    simulation. Default: OFF
``CoolingTimestepSafetyFactor`` (external)
    Described in ``UseCoolingTime``.  Default: 0.1
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
``SmallRho`` (external)
    Minimum value for density in code units. This is enforced in euler.F
    when using the PPM solver (``HydroMethod`` = 0) or in 
    hydro_rk/EvolveLevel_RK.C when ``HydroMethod`` is 3 or 4. Not enforced
    in other hydrodynamics methods. Default: 1e-30
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

``MHD_CT_Method`` (external) 
    Method for computing the electric field from the Riemann fluxes

    ========== ==========================================================================
    CT Method   Description  
    ========== ==========================================================================
    0           None (only for debugging)
    1           Balsara and Spicer 1999. First order average.
    2           Gardiner and Stone 2005. Second order Lax-Friedrichs type reconstruction.
                Uses ``CT_AthenaDissipation`` flag.
    3           Gardiner and Stone 2005.  Second order reconstruction using
                upwind switches
    ========== ==========================================================================

    Default: 3

``CT_AthenaDissipation``  (external) 
    For the Lax-Friedrichs CT method, this is the maximum wave speed.  (:math:`\alpha` in Gardiner & Stone 2005 eqn. 46). Default: 0.1

``EquationOfState`` (external, ct only) 
    0: standard adiabatic 1: Exactly isothermal
    equation of state.  This flag removes the total energy term completely, instead
    computing pressure as :math:`p = c^2 \rho`. This option only works with
    ``HydroMethod = 6`` and ``RiemannSolver = 6`` (HLLD) as this is the only purely
    isothermal Riemann solver in Enzo.  Default: 0

``IsothermalSoundSpeed`` (external, ct only) 
    When ``EquationOfState = 1``, this is the
    sound speed used for computation of pressure.  Default: 1

``MHDCTSlopeLimiter`` (external, ct only) 
    For computing derivatives for the reconstruction,
    this switches between zero slope (0), minmod (1), VanLeer (2), and
    characteristic  (3) characteristic with primitive limiting (4).  Default: 1

``ReconstructionMethod`` (external) 
    There are two reconstruction methods
    that work with MHDCT: Piecewise Linear Method (PLM) (0) and MUSCL-Hancock (6).  This
    formuation of MUSCL-Hancock is different from the 2nd order Runga Kutta used for
    ``HydroMethod = 3,4``.     

``RiemannSolver`` (external)  
    As with ``HydroMethod=4``, the prefered solver is
    HLLD (``RiemannSolver=6``).  Other solvers may be released if the DOE approves
    them.


``MHDCTUseSpecificEnergy`` (external) 
    Either specific energy is used internally
    (1) or conserved energy is used internally (0).  Minor difference in boundary
    condition update, included for comparison to old solutions.  Default: 1


``MHDCTDualEnergyMethod`` (external) 
    When ``DualEnergyFormalism = 1``, this switches
    between a method that solves an additional equation for the internal energy, as
    in the rest of Enzo, and method that updates the entropy.  


``MHD_WriteElectric`` (external)  
    Include the electric field in the output.
    Default: 0

``MHD_ProjectB`` (internal)  
    Project magnetic fields from fine to coarse.
    Should not be done in general, only used for initialization.  

``MHD_ProjectE`` (internal)  
    Project Electric fields from fine to coarse.
    Used for the time evolution of the fields.

.. _dedner_parameters:

Magnetohydrodynamics (Dedner) Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following parameters are considered only when ``HydroMethod`` is 3 or 4 (and occasionally only in some test problems).  
Because many of the following parameters are not actively being tested and maintained, users are encouraged to carefully examine the code before using it.

``UsePoissonDivergenceCleaning`` (external)
    Enables additional divergence cleaning by solving a Poisson equation.
    This works on top of the standard mixed hyperbolic/parabolic divergence cleaning
    and is in most cases not required.
    Works on indiviual grids, i.e., it's *not* a global divergence purge.
    Use with care as this feature is not extensively tested.
    Method 1 and 2 are a failed experiment to do divergence cleaning
    using successive over relaxation. Method 3 uses conjugate gradient
    with a 2 cell stencil and Method 4 uses a 4 cell stencil. 4 is more
    accurate but can lead to aliasing effects. Default: 0
``PoissonDivergenceCleaningBoundaryBuffer`` (external)
    Choose to *not* correct in the active zone of a grid by a
    boundary of cells this thick. Default: 0
``PoissonDivergenceCleaningThreshold`` (external)
    Calls divergence cleaning on a grid when magnetic field divergence
    is above this threshold. Default: 0.001
``PoissonApproximationThreshold`` (external)
    Controls the accuracy of the resulting solution for divergence
    cleaning Poisson solver. Default: 0.001
``PoissonBoundaryType`` (external)
    Controls the boundary conditions for divergence
    cleaning Poisson solver. 0 - Neumann (default). 1 - Dirichlet 
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
``MixSpeciesAndColors`` (external)
    This parameter enables color fields to be evolved as species in the MUSCL solvers. If ``PopIIISupernovaUseColour`` is on, this must also be turned on to trace the metal field. Default: 1

``SmallT`` (external)
    Minimum value for temperature in hydro_rk/EvolveLevel_RK.C.  Default: 1e-10 (note that the default value assumes UsePhysicalUnit = 1)
``SmallP``
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

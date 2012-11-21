Radiation Parameters
~~~~~~~~~~~~~~~~~~~~

.. _radiation_backgrounds:

Background Radiation Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``RadiationFieldType`` (external)
    This integer parameter specifies the type of radiation field that
    is to be used. Except for ``RadiationFieldType`` = 9, which should
    be used with ``MultiSpecies`` = 2, UV backgrounds can currently
    only be used with ``MultiSpecies`` = 1 (i.e. no molecular H
    support). The following values are used. Default: 0

   ::
  
     1  - Haardt & Madau spectrum with q_alpha = 1.5
     2  - Haardt & Madau spectrum with q_alpha = 1.8
     3  - Modified Haardt & Madau spectrum to match observations
     	  (Kirkman & Tytler 2005).
     4  - Haardt & Madau spectrum with q_alpha = 1.5 supplemented with an X-ray Compton heating
          background from Madau & Efstathiou (see astro-ph/9902080)
     9  - Constant molecular H2 photo-dissociation rate
     10 - Internally computed radiation field using the algorithm of Cen & Ostriker
     11 - Same as previous, but with very, very simple optical shielding fudge
     12 - Haardt & Madau spectrum with q_alpha = 1.57

``RadiationFieldLevelRecompute`` (external)
    This integer parameter is used only if the previous parameter is
    set to 10 or 11. It controls how often (i.e. the level at which)
    the internal radiation field is recomputed. Default: 0
``RadiationSpectrumNormalization`` (external)
    This parameter was initially used to normalize the photo-ionization
    and photo-heating rates computed in the function
    ``RadiationFieldCalculateRates()`` and then passed on to the
    ``calc_photo_rates()``, ``calc_rad()`` and ``calc_rates()`` routines.
    Later, the normalization as a separate input parameter was dropped
    for all cases by using the rates computed in
    ``RadiationFieldCalculateRates()`` with one exception: The molecular
    hydrogen (H2) dissociation rate. There a normalization is performed
    on the rate by multiplying it with ``RadiationSpectrumNormalization``.
    Default: 1e-21
``RadiationShield`` (external)
    This parameter specifies whether the user wants to employ
    approximate radiative-shielding. This parameter will be
    automatically turned on when RadiationFieldType is set to 11. See
    ``calc_photo_rates.src``. Default: 0
``RadiationFieldRedshift`` (external)
    This parameter specifies the redshift at which the radiation field
    is calculated.  Default: 0
``RadiationRedshiftOn`` (external) 
    The redshift at which the UV 
    background turns on. Default: 7.0.
``RadiationRedshiftFullOn`` (external) 
    The redshift at which the UV
    background is at full strength.  Between z =
    ``RadiationRedshiftOn`` and z = ``RadiationRedshiftFullOn``, the 
    background is gradually ramped up to full strength. Default: 6.0.
``RadiationRedshiftDropOff`` (external) 
    The redshift at which the 
    strength of the UV background is begins to gradually reduce,
    reaching zero by ``RadiationRedshiftOff``. Default: 0.0.
``RadiationRedshiftOff`` (external) 
    The redshift at which the UV 
    background is fully off. Default: 0.0.
``TabulatedLWBackground`` (external)
    When on, the amplitude of the Lyman-Werner background is read from the file LW_J21.in as a function of redshift.  Each line should have the redshift and LW background in units of 1e-21 erg/cm^3/s/Hz/sr.  Default: 0
``AdjustUVBackground`` (external)
    Add description. Default: 1.
``AdjustUVBackgroundHighRedshift`` (external)
    Add description. Default: 0.
``SetUVAmplitude`` (external)
    Add description. Default: 1.0.
``SetHeIIHeatingScale`` (external)
    Add description. Default: 1.8.
``RadiationSpectrumSlope`` (external)
    Add description. Default: 1.5.

.. _radiative_transfer_ray_tracing:

Radiative Transfer (Ray Tracing) Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``RadiativeTransfer`` (external)
    Set to 1 to turn on the adaptive ray tracing following Abel, Wise &
    Bryan 2007. Note that Enzo must be first recompiled after setting
    ``make photon-yes``. Default: 0.
``RadiativeTransferRadiationPressure`` (external)
    Set to 1 to turn on radiation pressure created from absorbed photon
    packages. Default: 0
``RadiativeTransferInitialHEALPixLevel`` (external)
    Chooses how many rays are emitted from radiation sources. The
    number of rays in Healpix are given through # =
    12x4\ :sup:`level`\ . Default: 3.
``RadiativeTransferRaysPerCell`` (external)
    Determines the accuracy of the scheme by giving the minimum number
    of rays to cross cells. The more the better (slower). Default: 5.1.
``RadiativeTransferSourceRadius`` (external)
    The radius at which the photons originate from the radiation
    source. A positive value results in a radiating sphere. Default: 0.
``RadiativeTransferPropagationRadius`` (external)
    The maximum distance a photon package can travel in one timestep.
    Currently unused. Default: 0.
``RadiativeTransferPropagationSpeed`` (external)
    The fraction of the speed of light at which the photons travel.
    Default: 1.
``RadiativeTransferCoupledRateSolver`` (external)
    Set to 1 to calculate the new ionization fractions and gas energies
    after every radiative transfer timestep. This option is highly
    recommended to be kept on. If not, ionization fronts will propagate too
    slowly. Default: 1.
``RadiativeTransferOpticallyThinH2`` (external)
    Set to 1 to include an optically-thin H_2 dissociating
    (Lyman-Werner) radiation field. Only used if ``MultiSpecies`` > 1. If
    ``MultiSpecies`` > 1 and this option is off, the Lyman-Werner radiation
    field will be calculated with ray tracing. Default: 1.
``RadiativeTransferSplitPhotonPackage`` (external)
    Once photons are past this radius, they can no longer split. In
    units of kpc. If this value is negative (by default), photons can
    always split. Default: ``FLOAT_UNDEFINED``.
``RadiativeTransferPhotonEscapeRadius`` (external)
    The number of photons that pass this distance from its source are
    summed into the global variable ``EscapedPhotonCount[]``. This variable
    also keeps track of the number of photons passing this radius
    multiplied by 0.5, 1, and 2. Units are in kpc. Not used if set to
    0. Default: 0.
``RadiativeTransferSourceClustering`` (external)
    Set to 1 to turn on ray merging from combined virtual sources on a
    binary tree. Default: 0.
``RadiativeTransferPhotonMergeRadius`` (external)
    The radius at which the rays will merge from their SuperSource,
    which is the luminosity weighted center of two sources. This radius
    is in units of the separation of two sources associated with one
    SuperSource. If set too small, there will be angular artifacts in
    the radiation field. Default: 2.5
``RadiativeTransferSourceBeamAngle`` (external)
    Rays will be emitted within this angle in degrees of the poles from sources with "Beamed" types.  Default: 30
``RadiativeTransferPeriodicBoundary`` (external)
    Set to 1 to turn on periodic boundary conditions for photon
    packages. Default: 0.
``RadiativeTransferTimestepVelocityLimit`` (external)
    Limits the radiative transfer timestep to a minimum value that is
    determined by the cell width at the finest level divided by this
    velocity. Units are in km/s. Default: 100.
``RadiativeTransferHIIRestrictedTimestep`` (external)
    Adaptive ray tracing timesteps will be restricted by a maximum change of 10% in neutral fraction if this parameter is set to 1.  If set to 2, then the incident flux can change by a maximum of 0.5 between cells.  See Wise & Abel (2011) in Sections 3.4.1 and 3.4.4 for more details.  Default: 0
``RadiativeTransferAdaptiveTimestep`` (external)
    Must be 1 when RadiativeTransferHIIRestrictedTimestep is non-zero.  When RadiativeTransferHIIRestrictedTimestep is 0, then the radiative transfer timestep is set to the timestep of the finest AMR level.  Default: 0
``RadiativeTransferLoadBalance`` (external)
    When turned on, the grids are load balanced based on the number of ray segments traced.  The grids are moved to different processors only for the radiative transfer solver.  Default: 0
``RadiativeTransferHydrogenOnly`` (external)
    When turned on, the photo-ionization fields are only created for hydrogen.  Default: 0
``RadiationXRaySecondaryIon`` (external)
    Set to 1 to turn on secondary ionizations and reduce heating from
    X-ray radiation (Shull & van Steenberg 1985). Currently only BH and
    MBH particles emit X-rays. Default: 0.
``RadiationXRayComptonHeating`` (external)
    Set to 1 to turn on Compton heating on electrons from X-ray
    radiation (Ciotti & Ostriker 2001). Currently only BH and MBH
    particles emit X-rays. Default: 0.
``RadiativeTransferInterpolateField`` (obsolete)
    A failed experiment in which we evaluate the density at the
    midpoint of the ray segment in each cell to calculate the optical
    depth. To interpolate, we need to calculate the vertex interpolated
    density fields. Default: 0.
``SimpleQ`` (external)
    Ionizing photon luminosity of a "simple radiating source" that is independent of mass.  In units of photons per second.  Default: 1e50
``SimpleRampTime`` (external)
    Time to exponential ramp up the luminosity of a simple radiating source.  In units of 1e6 years.  Default: 0.1
``RadiativeTransferTraceSpectrum`` (reserved)
    reserved for future experimentation. Default: 0.
``RadiativeTransferTraceSpectrumTable`` (reserved)
    reserved for future experimentation. Default: ``spectrum_table.dat``

.. _radiative_transfer_fld:

Radiative Transfer (FLD) Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``RadiativeTransferFLD`` (external)
    Set to 2 to turn on the fld-based radiation solvers following Reynolds,
    Hayes, Paschos & Norman, 2009. Note that you also have to compile
    the source using ``make photon-yes`` and a ``make
    hypre-yes``. Note that if FLD is turned on, it will force
    ``RadiativeCooling = 0``, ``GadgetEquilibriumCooling = 0``, and
    ``RadiationFieldType = 0`` to prevent conflicts. Default: 0.

    *IMPORTANT*: Set ``RadiativeTransfer = 0`` to avoid conflicts with the ray tracing solver above.
    Set ``RadiativeTransferOpticallyThinH2 = 0`` to avoid conflicts with the built-in optically-thin H_2 dissociating field from the ray-tracing solver. 
``ImplicitProblem`` (external)
    Set to 1 to turn on the implicit FLD solver, or 3 to turn on the
    split FLD solver. Default: 0.
``RadHydroParamfile`` (external)
    Names the (possibly-different) input parameter file containing
    solver options for the FLD-based solvers. These are described in
    the relevant User Guides, located in ``doc/implicit_fld`` and
    ``doc/split_fld``. Default: NULL.
``RadiativeTransferFLDCallOnLevel`` (reserved)
    The level in the static AMR hierarchy where the unigrid FLD solver
    should be called. Currently only works for 0 (the root grid).
    Default: 0.
``StarMakerEmissivityField`` (external)
    When compiled with the FLD radiation transfer >make emissivity-yes; make hypre-yes, setting this to 1 turns on the emissivity field to source the gray radiation. Default: 0
``uv_param`` (external)
    When using the FLD radiation transfer and StarMakerEmissivityFIeld = 1, this is the efficiency of mass to UV light ratio. Default: 0

Radiative Transfer (FLD) Implicit Solver Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    These parameters should be placed within the file named in
    ``RadHydroParamfile`` in the main parameter file. All are described in
    detail in the User Guide in ``doc/implicit_fld``.


``RadHydroESpectrum`` (external)
    Type of assumed radiation spectrum for radiation field, Default: 1.

   ::
 
    -1 - monochromatic spectrum at frequency h nu_{HI} = 13.6 eV
    0  - power law spectrum, (nu / nu_{HI} )^(-1.5) 
    1  - T = 1e5 blackbody spectrum

``RadHydroChemistry`` (external)
    Use of hydrogen chemistry in ionization model, set to 1 to turn on
    the hydrogen chemistry, 0 otherwise. Default: 1.
``RadHydroHFraction`` (external)
    Fraction of baryonic matter comprised of hydrogen. Default: 1.0.
``RadHydroModel`` (external)
    Determines which set of equations to use within the solver.
    Default: 1.

   ::
 
    1  - chemistry-dependent model, with case-B hydrogen II recombination coefficient.
    2  - chemistry-dependent model, with case-A hydrogen II recombination coefficient.
    4  - chemistry-dependent model, with case-A hydrogen II
       recombination coefficient, but assumes an isothermal gas energy.
    10 - no chemistry, instead uses a model of local thermodynamic
       equilibrium to couple radiation to gas energy.

``RadHydroMaxDt`` (external)
    maximum time step to use in the FLD solver. Default: 1e20 (no
    limit).
``RadHydroMinDt`` (external)
    minimum time step to use in the FLD solver. Default: 0.0 (no
    limit).
``RadHydroInitDt`` (external)
    initial time step to use in the FLD solver. Default: 1e20 (uses
    hydro time step).
``RadHydroDtNorm`` (external)
    type of p-norm to use in estimating time-accuracy for predicting
    next time step. Default: 2.0.    

   ::

     0 - use the max-norm.
    >0 - use the specified p-norm.
    <0 - illegal.

``RadHydroDtRadFac`` (external)
    Desired time accuracy tolerance for the radiation field. Default:
    1e20 (unused).
``RadHydroDtGasFac`` (external)
    Desired time accuracy tolerance for the gas energy field. Default:
    1e20 (unused).
``RadHydroDtChemFac`` (external)
    Desired time accuracy tolerance for the hydrogen I number density.
    Default: 1e20 (unused).
``RadiationScaling`` (external)
    Scaling factor for the radiation field, in case standard
    non-dimensionalization fails. Default: 1.0.
``EnergyCorrectionScaling`` (external)
    Scaling factor for the gas energy correction, in case standard
    non-dimensionalization fails. Default: 1.0.
``ChemistryScaling`` (external)
    Scaling factor for the hydrogen I number density, in case standard
    non-dimensionalization fails. Default: 1.0.
``RadiationBoundaryX0Faces`` (external)
    Boundary condition types to use on the x0 faces of the radiation
    field. Default: [0 0].

   ::
 
    0 - Periodic.
    1 - Dirichlet.
    2 - Neumann.

``RadiationBoundaryX1Faces`` (external)
    Boundary condition types to use on the x1 faces of the radiation
    field. Default: [0 0].
``RadiationBoundaryX2Faces`` (external)
    Boundary condition types to use on the x2 faces of the radiation
    field. Default: [0 0].
``RadHydroLimiterType`` (external)
    Type of flux limiter to use in the FLD approximation. Default: 4.

   ::

    0 - original Levermore-Pomraning limiter, à la Levermore & Pomraning, 1981 and Levermore, 1984.
    1 - rational approximation to LP limiter.
    2 - new approximation to LP limiter (to reduce floating-point cancellation error).
    3 - no limiter.
    4 - ZEUS limiter (limiter 2, but with no "effective albedo").

``RadHydroTheta`` (external)
    Time-discretization parameter to use, 0 gives explicit Euler, 1
    gives implicit Euler, 0.5 gives trapezoidal. Default: 1.0.
``RadHydroAnalyticChem`` (external)
    Type of time approximation to use on gas energy and chemistry
    equations. Default: 1 (if possible for model).

   ::

    0 - use a standard theta-method.
    1 - use an implicit quasi-steady state (IQSS) approximation.

``RadHydroInitialGuess`` (external)
    Type of algorithm to use in computing the initial guess for the
    time-evolved solution. Default: 0.

   ::
 
    0 - use the solution from the previous time step (safest).
    1 - use explicit Euler with only spatially-local physics (heating & cooling).
    2 - use explicit Euler with all physics.
    5 - use an analytic predictor based on IQSS approximation of
       spatially-local physics.

``RadHydroNewtTolerance`` (external)
    Desired accuracy for solution to satisfy nonlinear residual
    (measured in the RMS norm). Default: 1e-6.
``RadHydroNewtIters`` (external)
    Allowed number of Inexact Newton iterations to achieve tolerance
    before returning with FAIL. Default: 20.
``RadHydroINConst`` (external)
    Inexact Newton constant used in specifying tolerances for inner
    linear solver. Default: 1e-8.
``RadHydroMaxMGIters`` (external)
    Allowed number of iterations for the inner linear solver (geometric
    multigrid). Default: 50.
``RadHydroMGRelaxType`` (external)
    Relaxation method used by the multigrid solver. Default: 1.

    ::
    1 - Jacobi.
    2 - Weighted Jacobi.
    3 - Red/Black Gauss-Seidel (symmetric).
    4 - Red/Black Gauss-Seidel (non-symmetric).

``RadHydroMGPreRelax`` (external)
    Number of pre-relaxation sweeps used by the multigrid solver.
    Default: 1.
``RadHydroMGPostRelax`` (external)
    Number of post-relaxation sweeps used by the multigrid solver.
    Default: 1.
``EnergyOpacityC0``, ``EnergyOpacityC1``, ``EnergyOpacityC2``, ``EnergyOpacityC3``, ``EnergyOpacityC4`` (external)
    Parameters used in defining the energy-mean opacity used with
    ``RadHydroModel`` 10. Default: [1 1 0 1 0].
``PlanckOpacityC0``, ``PlanckOpacityC1``, ``PlanckOpacityC2``, ``PlanckOpacityC3``, ``PlanckOpacityC4`` (external)
    Parameters used in defining the Planck-mean opacity used with
    ``RadHydroModel`` 10. Default: [1 1 0 1 0].

Radiative Transfer (FLD) Split Solver Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    These parameters should be placed within the file named in
    ``RadHydroParamfile`` in the main parameter file. All are described in
    detail in the User Guide in ``doc/split_fld``.


``RadHydroESpectrum`` (external)
    Type of assumed radiation spectrum for radiation field, Default: 1.

   ::
 
    -1 - monochromatic spectrum at frequency h nu_{HI}= 13.6 eV
    0  - power law spectrum, (nu / nu_{HI})^(-1.5) 
    1  - T=1e5 blackbody spectrum

``RadHydroChemistry`` (external)
    Use of hydrogen chemistry in ionization model, set to 1 to turn on
    the hydrogen chemistry, 0 otherwise. Default: 1.
``RadHydroHFraction`` (external)
    Fraction of baryonic matter comprised of hydrogen. Default: 1.0.
``RadHydroModel`` (external)
    Determines which set of equations to use within the solver.
    Default: 1.

   ::
    1  - chemistry-dependent model, with case-B hydrogen II recombination
         coefficient.
    4  - chemistry-dependent model, with case-A hydrogen II recombination
         coefficient, but assumes an isothermal gas energy.
   10  - no chemistry, instead uses a model of local thermodynamic
       	 equilibrium to couple radiation to gas energy.


``RadHydroMaxDt`` (external)
    maximum time step to use in the FLD solver. Default: 1e20 (no
    limit).
``RadHydroMinDt`` (external)
    minimum time step to use in the FLD solver. Default: 0.0 (no
    limit).
``RadHydroInitDt`` (external)
    initial time step to use in the FLD solver. Default: 1e20 (uses
    hydro time step).
``RadHydroDtNorm`` (external)
    type of p-norm to use in estimating time-accuracy for predicting
    next time step. Default: 2.0.

   ::
    0  - use the max-norm.
    >0 - use the specified p-norm.
    <0 - illegal.

``RadHydroDtRadFac`` (external)
    Desired time accuracy tolerance for the radiation field. Default:
    1e20 (unused).
``RadHydroDtGasFac`` (external)
    Desired time accuracy tolerance for the gas energy field. Default:
    1e20 (unused).
``RadHydroDtChemFac`` (external)
    Desired time accuracy tolerance for the hydrogen I number density.
    Default: 1e20 (unused).
``RadiationScaling`` (external)
    Scaling factor for the radiation field, in case standard
    non-dimensionalization fails. Default: 1.0.
``EnergyCorrectionScaling`` (external)
    Scaling factor for the gas energy correction, in case standard
    non-dimensionalization fails. Default: 1.0.
``ChemistryScaling`` (external)
    Scaling factor for the hydrogen I number density, in case standard
    non-dimensionalization fails. Default: 1.0.
``RadiationBoundaryX0Faces`` (external)
    Boundary condition types to use on the x0 faces of the radiation
    field. Default: [0 0].

    ::

     0 - Periodic.
     1 - Dirichlet.
     2 - Neumann.

``RadiationBoundaryX1Faces`` (external)
    Boundary condition types to use on the x1 faces of the radiation
    field. Default: [0 0].
``RadiationBoundaryX2Faces`` (external)
    Boundary condition types to use on the x2 faces of the radiation
    field. Default: [0 0].
``RadHydroTheta`` (external)
    Time-discretization parameter to use, 0 gives explicit Euler, 1
    gives implicit Euler, 0.5 gives trapezoidal. Default: 1.0.
``RadHydroSolTolerance`` (external)
    Desired accuracy for solution to satisfy linear residual (measured
    in the 2-norm). Default: 1e-8.
``RadHydroMaxMGIters`` (external)
    Allowed number of iterations for the inner linear solver (geometric
    multigrid). Default: 50.
``RadHydroMGRelaxType`` (external)
    Relaxation method used by the multigrid solver. Default: 1.

    ::

     Jacobi.
     Weighted Jacobi.
     Red/Black Gauss-Seidel (symmetric).
     Red/Black Gauss-Seidel (non-symmetric).

``RadHydroMGPreRelax`` (external)
    Number of pre-relaxation sweeps used by the multigrid solver.
    Default: 1.
``RadHydroMGPostRelax`` (external)
    Number of post-relaxation sweeps used by the multigrid solver.
    Default: 1.
``EnergyOpacityC0``, ``EnergyOpacityC1``, ``EnergyOpacityC2`` (external)
    Parameters used in defining the energy-mean opacity used with
    RadHydroModel 10. Default: [1 1 0].


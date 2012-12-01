.. _initialization_parameters:

Initialization Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

``TopGridRank`` (external)
    This specifies the dimensionality of the root grid and by extension
    the entire hierarchy. It should be 1,2 or 3. Default: none
``TopGridDimensions`` (external)
    This is the dimension of the top or root grid. It should consist of
    1, 2 or 3 integers separated by spaces. For those familiar with the
    KRONOS or ZEUS method of specifying dimensions, these values do not
    include ghost or boundary zones. A dimension cannot be less than 3
    zones wide and more than ``MAX_ANY_SINGLE_DIRECTION`` -
    ``NumberOfGhostZones``\*2. ``MAX_ANY_SINGLE_DIRECTION`` is defined in
    ``fortran.def``. Default: none
``DomainLeftEdge``, ``DomainRightEdge`` (external)
    These float values specify the two corners of the problem domain
    (in code units). The defaults are: 0 0 0 for the left edge and 1 1
    1 for the right edge.
``LeftFaceBoundaryCondition``, ``RightFaceBoundaryCondition`` (external)
    These two parameters each consist of vectors of integers (of length
    ``TopGridRank``). They specify the boundary conditions for the top grid
    (and hence the entire hierarchy). The first integer corresponds to
    the x-direction, the second to the y-direction and the third, the
    z-direction. The possible values are: 0 - reflecting, 1 - outflow,
    2 - inflow, 3 - periodic, 4 - shearing. For inflow, the inflow
    values can be set through the next parameter, or more commonly are
    controlled by problem-specific code triggered by the ``ProblemType``.
    For shearing boundaries, the boundary pair in another direction
    must be periodic. Note that self gravity will not be consistent
    with shearing boundary conditions. Default: 0 0 0
``BoundaryConditionName`` (external)
    While the above parameters provide an easy way to set an entire
    side of grid to a given boundary value, the possibility exists to
    set the boundary conditions on an individual cell basis. This is
    most often done with problem specific code, but it can also be set
    by specifying a file which contains the information in the
    appropriate format. This is too involved to go into here. Default:
    none
``InitialTime`` (internal)
    The time, in code units, of the current step. For cosmology the
    units are in free-fall times at the initial epoch (see :ref:`EnzoOutputFormats`). Default: generally 0, depending on problem
``Initialdt`` (internal)
    The timestep, in code units, for the current step. For cosmology
    the units are in free-fall times at the initial epoch (see :ref:`EnzoOutputFormats`). Default: generally 0, depending on problem
``Unigrid`` (external)
    This parameter should be set to 1 (TRUE) for large cases--AMR as
    well as non-AMR--where the root grid is 512\ :sup:`3`\  or larger.
    This prevents initialization under subgrids at start up, which is
    unnecessary in cases with simple non-nested initial conditions.
    Unigrid must be set to 0 (FALSE) for cases with nested initial
    conditions. Default: 0 (FALSE). See also ``ParallelRootGridIO`` in :ref:`io_parameters`.
``UnigridTranspose`` (external)
    This parameter governs the fast FFT bookkeeping for Unigrid runs.
    Does not work with isolated gravity. Default: 0 (FALSE). See also
    ``Unigrid`` above.
``MaximumTopGridTimeStep`` (external)
    This parameter limits the maximum timestep on the root grid.  Default: huge_number.
``ShearingVelocityDirection`` (external)
    Select direction of shearing boundary. Default is x direction. Changing this is probably not a good idea.
``AngularVelocity`` (external)
    The value of the angular velocity in the shearing boundary.
    Default: 0.001
``VelocityGradient`` (external)
    The value of the per code length gradient in the angular velocity
    in the shearing boundary. Default: 1.0
``GridVelocity`` (external)
    The whole computational domain will have this velocity.  Experimental.  Default: 0 0 0
``StringKick`` (external)
    While this parameter was initially designed to describe the kick by cosmic strings in CosmologySimulation, it can be used to model the velocity (in km/s) that the baryons should move relative to dark matter at the initial redshift, in order to study the effect discussed by Tseliakhovich & Hirata (astro-ph:1005.2416). Default: 0
``StringKickDimension`` (external)
    This parameter is used to control the orthogonal direction of the flow.  Default: 0 (x-axis)
``MemoryLimit`` (external)
    If the memory usage on a single MPI process exceeds this number, then the simulation will halt after outputting.  Only used when the compile-time define MEM_TRACE is used. Default: 4e9
``HydrogenFractionByMass`` (external)
    This parameter is used to set up initial conditions in some test problems.  Default: 0.76
``DeuteriumToHydrogenRatio`` (external)
    This parameter is used to set up initial conditions in some test problems.  Default: 2.0*3.4e-5 (Burles & Tytler 1998, the parameter here is by mass, so multiply by 2)
``SolarMetalFractionByMass`` (external)
    This parameter is used to set up initial conditions in some test problems. Do NOT change this parameter unless you know exactly what you are doing. Default: 0.02041
``CoolDataIh2co`` (external)
    Whether to include molecular hydrogen cooling.  Do NOT change this parameter unless you know exactly what you are doing.  Default: 1
``CoolDataIpiht`` (external)
    Whether to include photoionization heating.  Do NOT change this parameter unless you know exactly what you are doing.  Default: 1
``CoolDataCompXray`` (external)
    Do NOT change this parameter unless you know exactly what you are doing.  Saved to CoolData.comp_xray. Default: 0
``CoolDataTempXray`` (external)
    Do NOT change this parameter unless you know exactly what you are doing.  Saved to CoolData.temp_xray. Default: 0
``NumberOfTemperatureBins`` (external)
    Do NOT change this parameter unless you know exactly what you are doing. Default: 600
``TemperatureStart`` (external)
    Do NOT change this parameter unless you know exactly what you are doing. Default: 10
``TemperatureEnd`` (external)
    Do NOT change this parameter unless you know exactly what you are doing. Default: 1e8
``ExternalBoundaryIO`` (external)
    not recommended for use at this point. Only works if compiled with ``ooc-boundary-yes``.  Default: 0
``ExternalBoundaryTypeIO`` (external)
    not recommended for use at this point. Default: 0
``ExternalBoundaryValueIO`` (external)
    not recommended for use at this point. Default: 0
``SimpleConstantBoundary`` (external)
    not recommended for use at this point. Default: 0
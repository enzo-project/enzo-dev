Hierarchy Control Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``StaticHierarchy`` (external)
    A flag which indicates if the hierarchy is static (1) or dynamic
    (0). In other words, a value of 1 takes the A out of AMR. Default:
    1
``RefineBy`` (external)
    This is the refinement factor between a grid and its subgrid. For
    cosmology simulations, we have found a ratio of 2 to be most useful.
    Default: 4
``MaximumRefinementLevel`` (external)
    This is the lowest (most refined) depth that the code will produce.
    It is zero based, so the total number of levels (including the root
    grid) is one more than this value. Default: 2
``CellFlaggingMethod`` (external)
    The method(s) used to specify when a cell should be refined. This
    is a list of integers, up to 9, as described by the following
    table. The methods combine in an "OR" fashion: if any of them
    indicate that a cell should be refined, then it is flagged. For
    cosmology simulations, methods 2 and 4 are probably most useful.
    Note that some methods have additional parameters which are
    described below. For more information about specific methods, see the
    method paper. Default: 1

    ================== ==========================================================
    CellFlaggingMethod Description
    ================== ==========================================================
    1                  Refine by slope
    2                  Refine by baryon mass
    3                  Refine by shocks
    4                  Refine by particle mass
    5                  Refine by baryon overdensity
    6                  Refine by Jeans length
    7                  Refine if (cooling time < cell width/sound speed)
    8                  Refine by must-refine particles
    9                  Refine by shear
    10                 Refine by optical depth (in RT calculation)
    11                 Refine by resistive length (in MHD calculation)
    12                 Refine by defined region "MustRefineRegion"
    13                 Refine by metallicity
    14                 Refine by shockwaves (found w/shock finder)
    15                 Refine by normalized second derivative
    16                 Refine by Jeans length from the inertial tensor
    100                Avoid refinement based on ForbiddenRefinement field
    101                Avoid refinement in regions defined in "AvoidRefineRegion"
    ================== ==========================================================

``RefineRegionLeftEdge``, ``RefineRegionRightEdge`` (external)
    These two parameters control the region in which refinement is
    permitted. Each is a vector of floats (of length given by the
    problem rank) and they specify the two corners of a volume.
    Default: set equal to ``DomainLeftEdge`` and ``DomainRightEdge``.
``RefineRegionAutoAdjust`` (external)
    This is useful for multiresolution simulations with particles in
    which the particles have varying mass. Set to 1 to automatically
    adjust the refine region at root grid timesteps to only contain
    high-resolution particles. This makes sure that the fine regions do
    not contain more massive particles which may lead to small
    particles orbiting them or other undesired outcomes. Setting to any
    integer (for example, 3) will make AdjustRefineRegion to work at
    (RefineRegionAutoAdjust-1)th level timesteps because sometimes the
    heavy particles are coming into the fine regions too fast that you
    need more frequent protection. Default: 0.
``RefineRegionTimeType`` (external)
    If set, this controls how the first column of a refinement region
    evolution file (see below) is interpreted, 0 for code time, 1 for
    redshift. Default: -1, which is equivalent to 'off'.
``RefineRegionFile`` (external)
    The name of a text file containing the corners of the time-evolving
    refinement region. The lines in the file change the values of
    ``RefineRegionLeft/RightEdge`` during the course of the simulation, and
    the lines are ordered in the file from early times to late times.
    The first column of data is the time index (in code units or
    redshift, see the parameter above) for the next six columns, which
    are the values of ``RefineRegionLeft/RightEdge``. For example, this
    might be two lines from the text file when time is indexed by
    redshift:
    ::

        0.60 0.530 0.612 0.185 0.591 0.667 0.208
        0.55 0.520 0.607 0.181 0.584 0.653 0.201

    In this case, the refinement region stays at the z=0.60 value
    until z=0.55, when the box moves slightly closer to the (0,0,0)
    corner. There is a maximum of 300 lines in the file and there is no
    comment header line. Default: None.
``MinimumOverDensityForRefinement`` (external)
    These float values (up to 9) are used if the
    ``CellFlaggingMethod`` is 2, 4 or 5. For method 2 and 4, the value is the density (baryon or particle), in code units, above which refinement occurs. When using method 5, it becomes rho [code] - 1. The elements in this array must match those in ``CellFlaggingMethod``. Therefore, if ``CellFlaggingMethod`` = 1 4 9 10, ``MinimumOverDensityForRefinement`` = 0 8.0 0 0.

    In practice, this value is converted into a mass by
    multiplying it by the volume of the top grid cell. The result is
    then stored in the next parameter (unless that is set directly in
    which case this parameter is ignored), and this defines the mass
    resolution of the simulation. Note that the volume is of a top grid
    cell, so if you are doing a multi-grid initialization, you must
    divide this number by r\ :sup:`(d\*l)`\  where r is the refinement
    factor, d is the dimensionality and l is the (zero-based) lowest
    level. For example, for a two grid cosmology setup where a cell should be
    refined whenever the mass exceeds 4 times the mean density of the
    subgrid, this value should be 4 / (2\ :sup:`(3\*1)`\ ) = 4 / 8 =
    0.5. Keep in mind that this parameter has no effect if it is
    changed in a restart output; if you want to change the refinement
    mid-run you will have to modify the next parameter. Up to 9
    numbers may be specified here, each corresponding to the respective
    ``CellFlaggingMethod``. Default: 1.5
``MinimumMassForRefinement`` (internal)
    This float is usually set by the parameter above and so is labeled
    internal, but it can be set by hand. For non-cosmological simulations, it can be the easier refinement criteria to specify. It is the mass above
    which a refinement occurs if the ``CellFlaggingMethod`` is
    appropriately set. For cosmological simulations, it is specified in units such
    that the entire mass in the computational volume is 1.0, otherwise it is in code units. There are 9 numbers here again, as per the
    above parameter. Default: none
``MinimumMassForRefinementLevelExponent`` (external).
    This parameter modifies the behaviour of the above parameter. As it
    stands, the refinement based on the ``MinimumMassForRefinement``
    (hereafter Mmin) parameter is complete Lagrangian. However, this
    can be modified. The actual mass used is
    Mmin\*r\ :sup:`(l\*alpha)`\  where r is the refinement factor, l is
    the level and alpha is the value of this parameter
    (``MinimumMassForRefinementLevelExponent``). Therefore a negative value
    makes the refinement super-Lagrangian, while positive values are
    sub-Lagrangian. There are up to 9 values specified here, as per
    the above two parameters. Default: 0.0
``SlopeFlaggingFields[#]`` (external)
    If ``CellFlaggingMethod`` is 1, and you only want to refine on the
    slopes of certain fields then you can enter the number IDs of the
    fields. Default: Refine on slopes of all fields.
``MinimumSlopeForRefinement`` (external)
    If ``CellFlaggingMethod`` is 1, then local gradients are used as the
    refinement criteria. All variables are examined and the relative
    slope is computed: abs(q(i+1)-q(i-1))/q(i). Where this value
    exceeds this parameter, the cell is marked for refinement. This
    causes problems if q(i) is near zero. This is a single integer (as
    opposed to the list of five for the above parameters). Entering
    multiple numbers here correspond to the fields listed in
    ``SlopeFlaggingFields``. Default: 0.3
``MinimumPressureJumpForRefinement`` (external)
    If refinement is done by shocks, then this is the minimum
    (relative) pressure jump in one-dimension to qualify for a shock.
    The definition is rather standard (see Colella and Woodward's PPM
    paper for example) Default: 0.33
``MinimumEnergyRatioForRefinement`` (external)
    For the dual energy formalism, and cell flagging by
    shock-detection, this is an extra filter which removes weak shocks
    (or noise in the dual energy fields) from triggering the shock
    detection. Default: 0.1
``MinimumShearForRefinement`` (external)
    It is the minimum shear above which a refinement occurs if the CellFlaggingMethod is appropriately set. Default: 0
``OldShearMethod`` (external)
    If using the shear refinement criterion, setting this variable to 1 enables 
    the old method for calculating the shear criterion, which actually 
    calculates it based on shear and vorticity and makes some assumptions
    about the simulations (c_s=1, etc.).  However, this is necessary
    if you want to reproduce some of the old enzo results 
    (e.g. Kritsuk et al. 2006).  Default: 0
``MetallicityRefinementMinMetallicity`` (external)
    This is the threshold metallicity (in units of solar metallicity)
    above which cells must be refined to a minimum level of
    ``MetallicityRefinementMinLevel``. Default: 1.0e-5
``MetallicityRefinementMinLevel`` (external)
    Sets the minimum level (maximum cell size) to which a cell enriched
    with metal above a level set by ``MetallicityRefinementMinMetallicity``
    will be refined. This can be set to any level up to and including
    ``MaximumRefinementLevel``. (No default setting)
``MetallicityRefinementMinDensity`` (external)
    It is the minimum density above which a refinement occurs when the cells are refined on metallicity.  Default: FLOAT_UNDEFINED
``ShockwaveRefinementMinMach`` (external)
    The minimum Mach number required to refine a level when using ShockwaveRefinement. Default: 1.3
``ShockwaveRefinementMinVelocity`` (external)
    The minimum shock velocity required to refine a level when using ShockwaveRefinement. Default: 1.0e7 (cm/s)
``ShockwaveRefinementMaxLevel`` (external)
    The maximum level to refine to using the ShockwaveRefinement criteria. Default: 0 (not used)
``SecondDerivativeFlaggingFields`` (external)
    The field indices (list of up to 7) that are used for the normalized second
    derivative refinement criteria. Default: INT_UNDEFINED
``MinimumSecondDerivativeForRefinement`` (external)
    The value of the second derivative above which a cell will be flagged for
    refinement. Each value in this list (of up to 7 values) should be between
    0.0 and 1.0.  Values between 0.3-0.8 are recommended.  Default: 0.3
``SecondDerivativeEpsilon`` (external)
    Used to avoid refining around oscillations/fluctuations in the normalized
    second derivative refinement method.  The higher the value, the more it
    will filter out.  For fluid instability simulations, a value of ~0.01 is
    good.  For full-physics simulations, values around ~0.2 are recommended. Be
    aware that fluctuations on this scale in initial conditions may cause
    immediate refinement to the maximum level.  Default: 1.0e-2
``RefineByJeansLengthSafetyFactor`` (external)
    If the Jeans length refinement criterion (see ``CellFlaggingMethod``)
    is being used, then this parameter specifies the number of cells
    which must cover one Jeans length. Default: 4
``JeansRefinementColdTemperature`` (external)
    If the Jeans length refinement criterion (see ``CellFlaggingMethod``)
    is being used, and this parameter is greater than zero, it will be
    used in place of the temperature in all cells. Default: -1.0
``RefineByResistiveLengthSafetyFactor`` (external)
    Resistive length is defined as the curl of the magnetic field over
    the magnitude of the magnetic field. We make sure this length is
    covered by this number of cells. i.w. The resistive length in a MHD simulation should not be smaller than CellWidth * RefineByResistiveLengthSafetyFactor.  Default: 2.0
``MustRefineParticlesRefineToLevel`` (external)
    The maximum level on which ``MustRefineParticles`` are required to
    refine to. Currently sink particles and MBH particles are required
    to be sitting at this level at all times. Default: 0
``MustRefineParticlesRefineToLevelAutoAdjust`` (external)
    The parameter above might not be handy in cosmological simulations
    if you want your ``MustRefineParticles`` to be refined to a certain
    physical length, not to a level whose cell size keeps changing.
    This parameter (positive integer in pc) allows you to do just that.
    For example, if you set ``MustRefineParticlesRefineToLevelAutoAdjust``
    = 128 (pc), then the code will automatically calculate
    ``MustRefineParticlesRefineToLevel`` using the boxsize and redshift
    information. Default: 0 (FALSE)
``MustRefineParticlesMinimumMass`` (external)
    This was an experimental parameter to set a minimum for ``MustRefineParticles``.  Default: 0.0
``MustRefineRegionMinRefinementLevel`` (external)
    Minimum level to which the rectangular solid volume defined by
    ``MustRefineRegionLeftEdge`` and ``MustRefineRegionRightEdge`` will be
    refined to at all times. (No default setting)
``MustRefineRegionLeftEdge`` (external)
    Bottom-left corner of refinement region. Must be within the overall
    refinement region. Default: 0.0 0.0 0.0
``MustRefineRegionRightEdge`` (external)
    Top-right corner of refinement region. Must be within the overall
    refinement region. Default: 1.0 1.0 1.0
``StaticRefineRegionLevel[#]`` (external)
    This parameter is used to specify regions of the problem that are
    to be statically refined, regardless of other parameters. This is mostly
    used as an internal mechanism to keep the initial grid hierarchy in
    place, but can be specified by the user. Up to 20 static regions
    may be defined (this number set in ``macros_and_parameters.h``), and
    each static region is labeled starting from zero. For each static
    refined region, two pieces of information are required: (1) the
    region (see the next two parameters), and (2) the level at which
    the refinement is to occurs (0 implies a level 1 region will always
    exist). Default: none
``StaticRefineRegionLeftEdge[#]``, ``StaticRefineRegionRightEdge[#]`` (external)
    These two parameters specify the two corners of a statically
    refined region (see the previous parameter). Default: none
``AvoidRefineRegionLevel[#]`` (external)
    This parameter is used to limit the refinement to this level in a
    rectangular region.  Up to MAX_STATIC_REGIONS regions can be used.
``AvoidRefineRegionLeftEdge[#]``, ``AvoidRefineRegionRightEdge[#]`` (external) 
    These two parameters specify the two corners of a region that
    limits refinement to a certain level (see the previous
    parameter). Default: none
``MultiRefineRegionGeometry[#]`` (external)
    This parameter (and the ones following) describe a physical region of the simulation box for which an 
    independent refinement maximum and minimum (separate from ``MaximumRefinementLevel``) can be specified.
``MultiRefineRegionGeometry[#]`` controls the geometry of the refined volume. Currently implemented 
    geometries are: (0) a rectangular region, (1) a ring of infinite height and (2) a cylinder of infinite 
    height. Up to 20 multi-refined regions may be defined (number the same as for ``StaticRefineRegion``)
    and each multi-refined region is labelled starting from zero. Default: -1 (no multi-regions)
``MultiRefineRegionLeftEdge[#]``, ``MultiRefineRegionRightEdge[#]`` (external)
    Used when ``MultiRefineRegionGeometry[#] = 0`` and specifies the two corners in code units of a 
    rectagular multi-region with a given maximum and minimum refinement level. Default: none.
``MultiRefineRegionCenter[#]`` (external)
    Used when ``MultiRefineRegionGeometry[#] = 1 or 2`` and specifies the center of the ring or cylinder 
    in code units. Default: none
``MultiRefineRegionRadius[#]`` (external)
    Used when ``MultiRefineRegionGeometry[#] = 1 or 2`` and specifies the radius of the ring or cylinder 
    in code units. In the case of the ring, this marks the distance to the middle of the ring's thickness. 
    The thickness is specified with ``MultiRefineRegionWidth``. Default: none
``MultiRefineRegionWidth[#]`` (external)
    Used when ``MultiRefineRegionGeometry[#] = 1`` and specifies the width (thickness) of the ring in 
    code units. Default: none
``MultiRefineRegionOrientation[#]`` (external)
    Used when ``MultiRefineRegionGeometry[#] = 1 or 2`` and is a unit vector pointing along the vertical
    direction of the ring or cylinder. Default: none.
``MultiRefineRegionStaggeredRefinement[#]`` (external)
    Used when ``MultiRefineRegionGeometry[#] = 1 or 2``. To avoid a sharp change in refinement at the edge of
    the ring or cylinder, the allowed refinement is staggered from the maximum allowed value outside the 
    region, ``MultiRefineRegionOuterMaximumLevel``, to the maximum allowed refinement inside the region, 
    ``MultiRefineRegionMaximumLevel``. This parameter is the length over which that staggering occurs in 
    code units. Default: 0.0 (no staggering)
``MultiRefineRegionMaximumLevel[#]``, ``MultiRefineRegionMinimumLevel[#]`` (external)
    Maximum and minimum allowed refinement inside the region. Default: ``MaximumRefinementLevel``, 0
``MultiRefineRegionMaximumOuterLevel``, ``MultiRefineRegionMinimumOuterLevel`` (external)
    Maximum and minimum allowed refinement outside all regions. Default: ``MaximumRefinementLevel``, 0
``MinimumEfficiency`` (external)
    When new grids are created during the rebuilding process, each grid
    is split up by a recursive bisection process that continues until a
    subgrid is either of a minimum size or has an efficiency higher
    than this value. The efficiency is the ratio of flagged zones
    (those requiring refinement) to the total number of zones in the
    grid. This is a number between 0 and 1 and should probably by
    around 0.4 for standard three-dimensional runs. Default: 0.2
``NumberOfBufferZones`` (external)
    Each flagged cell, during the regridding process, is surrounded by
    a number of zones to prevent the phenomenon of interest from
    leaving the refined region before the next regrid. This integer
    parameter controls the number required, which should almost always
    be one. Default: 1
``MinimumSubgridEdge`` (external)
    The minimum length of the edge of a subgrid.  See :ref:`running_large_simulations`. Default: 6
``MaximumSubgridSize`` (external)
    The maximum size (volume) of a subgrid.  See :ref:`running_large_simulations`. Default: 32768
``SubgridSizeAutoAdjust`` (external)
    See :ref:`running_large_simulations`.  Default: 1 (TRUE)
``OptimalSubgridsPerProcessor`` (external)
    See :ref:`running_large_simulations`.  Default: 16
``LoadBalancing`` (external)
    Set to 0 to keep child grids on the same processor as their
    parents. Set to 1 to balance the work on one level over all
    processors. Set to 2 or 3 to load balance the grids but keep them
    on the same node. Option 2 assumes grouped scheduling, i.e. proc #
    = (01234567) reside on node (00112233) if there are 4 nodes. Option
    3 assumes round-robin scheduling (proc = (01234567) -> node =
    (01230123)). Set to 4 for load balancing along a Hilbert
    space-filling curve on each level. See :ref:`running_large_simulations`. Default: 1
``LoadBalancingCycleSkip`` (external)
    This sets how many cycles pass before we load balance the root
    grids. Only works with LoadBalancing set to 2 or 3. NOT RECOMMENDED
    for nested grid calculations. Default: 10
``LoadBalancingMinLevel`` (external)
    Load balance the grids in levels greater than this parameter.  Default: 0
``LoadBalancingMaxLevel`` (external)
    Load balance the grids in levels less than this parameter.  Default: MAX_DEPTH_OF_HIERARCHY
``ResetLoadBalancing`` (external)
    When restarting a simulation, this parameter resets the processor number of each root grid to be sequential.  All child grids are assigned to the processor of their parent grid.  Only implemented for LoadBalancing = 1.  Default = 0
``NumberOfRootGridTilesPerDimensionPerProcessor`` (external)
    Splits the root grid into 2^(dimensions*this parameter) grids per MPI process.  Default: 1
``UserDefinedRootGridLayout`` (external)
   A three element array.  Splits the root grid into ``N`` subgrids where ``N``
   is the product of the supplied values.  The first entry corresponds to the
   number of root grid decompositions along the x axis of the simulation, the
   second element the number of decompositions along the y axis, and the third
   the number of decompositions along the z axis.

   This parameter is only used if all three elements of the array are set to a
   value different from the dummy default value.  If that is the case the root
   grid will be *manually* decomposed and the value supplied for
   ``NumberOfRootGridTilesPerDimensionPerProcessor`` will be ignored.  This is
   most useful when an automatic root grid decomposition is inefficient (for
   example, in a deeply nested isolated galaxy simulation).

   This parameter should be used with caution since it is possible to get into
   a situation where there are fewer grids than CPU cores.  Normally this can
   never happen since there will always be at least one root grid tile for every
   CPU.  Most simulations assume you will be running with as many root grid
   tiles as CPUs - if you instead opt to reduce the number of root grid tiles
   per CPU to a number less than one, Enzo might break in unpredictable ways.
   Default: -99999 -99999 -99999

``FastSiblingLocatorEntireDomain`` (external)
    In zoom-in calculations, the fast sibling locator doesn't need to search the entire domain.  Turning this parameter on restricts the finder to the inner nested grid.  Currently broken.  Default: 0
``MoveParticlesBetweenSiblings`` (external)
    During RebuildHierarchy, particles that have moved beyond the grid boundaries are moved to the correct grid.  Default: 1
``RebuildHierarchyCycleSkip`` (external)
    Set the number of cycles at a given level before rebuilding the hierarchy.  Example: RebuildHierarchyCycleSkip[1] = 4

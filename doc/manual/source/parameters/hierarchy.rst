Hierarchy Control Parameters
----------------------------

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
    described below. Default: 1

    :: 

       1 - refine by slope		       9 - refine by shear
       2 - refine by baryon mass	       10 - refine by optical depth (in radiation calculation)
       3 - refine by shocks 		       11 - refine by resistive length (in MHD calculation)
       4 - refine by particle mass	       12 - refine by defined region "MustRefineRegion"
       5 - refine by baryon overdensity	       13 - refine by metallicity
       	  (currently disabled)                 14 - refine around shockwaves
       6  - refine by Jeans length             16 - refine by Jeans length from the inertial tensor	       
       7  - refine if (cooling time < cell width/sound speed)
       8 - refine by must-refine particles
       100 - avoid refinement based on ForbiddenRefinement field
       101 - avoid refinement in regions defined in "AvoidRefineRegion" 

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
``MetallicityRefinementMinLevel`` (external)
    Sets the minimum level (maximum cell size) to which a cell enriched
    with metal above a level set by ``MetallicityRefinementMinMetallicity``
    will be refined. This can be set to any level up to and including
    ``MaximumRefinementLevel``. (No default setting)
``MetallicityRefinementMinMetallicity`` (external)
    This is the threshold metallicity (in units of solar metallicity)
    above which cells must be refined to a minimum level of
    ``MetallicityRefinementMinLevel``. Default: 1.0e-5
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
``FluxCorrection`` (external)
    This flag indicates if the flux fix-up step should be carried out
    around the boundaries of the sub-grid to preserve conservation (1 -
    on, 0 - off). Strictly speaking this should always be used, but we
    have found it to lead to a less accurate solution for cosmological
    simulations because of the relatively sharp density gradients
    involved. However, it does appear to be important when radiative
    cooling is turned on and very dense structures are created.
    It does work with the ZEUS
    hydro method, but since velocity is face-centered, momentum flux is
    not corrected. Species quantities are not flux corrected directly
    but are modified to keep the fraction constant based on the density
    change. Default: 1
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
``RefineByJeansLengthSafetyFactor`` (external)
    If the Jeans length refinement criterion (see ``CellFlaggingMethod``)
    is being used, then this parameter specifies the number of cells
    which must cover one Jeans length. Default: 4
``JeansRefinementColdTemperature`` (external)
    If the Jeans length refinement criterion (see ``CellFlaggingMethod``)
    is being used, and this parameter is greater than zero, it will be
    used in place of the temperature in all cells. Default: -1.0
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
``RefineByResistiveLength`` (external)
    Resistive length is defined as the curl of the magnetic field over
    the magnitude of the magnetic field. We make sure this length is
    covered by this number of cells. Default: 2
``LoadBalancing`` (external)
    Set to 0 to keep child grids on the same processor as their
    parents. Set to 1 to balance the work on one level over all
    processors. Set to 2 or 3 to load balance the grids but keep them
    on the same node. Option 2 assumes grouped scheduling, i.e. proc #
    = (01234567) reside on node (00112233) if there are 4 nodes. Option
    3 assumes round-robin scheduling (proc = (01234567) -> node =
    (01230123)). Set to 4 for load balancing along a Hilbert
    space-filling curve on each level. Default: 1
``LoadBalancingCycleSkip`` (external)
    This sets how many cycles pass before we load balance the root
    grids. Only works with LoadBalancing set to 2 or 3. NOT RECOMMENDED
    for nested grid calculations. Default: 10

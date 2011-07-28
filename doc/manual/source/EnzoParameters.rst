.. _parameters:

Enzo Parameter List
===================

The following is a largely complete list of the parameters that Enzo
understands, and a brief description of what they mean. They are grouped
roughly by meaning; an alphabetical list is also available. Parameters for
individual test problems are also listed here.

This parameter list has two purposes. The first is to describe and explain the
parameters that can be put into the initial parameter file that begins a run.
The second is to provide a comprehensive list of all parameters that the code
uses, including those that go into an output file (which contains a complete
list of all parameters), so that users can better understand these output
files.

The parameters fall into a number of categories:

**external**
    These are user parameters in the sense that they can be set in the
    parameter file, and provide the primary means of communication
    between Enzo and the user.
**internal**
    These are mostly not set in the parameter file (although strictly
    speaking they can be) and are generally used for program to
    communicate with itself (via the restart of output files).
**obsolete**
    No longer used.
**reserved**
    To be used later.

Generally the external parameters are the only ones that are modified or set,
but the internal parameters can provide useful information and can sometimes be
modified so I list them here as well. Some parameters are true/false or on/off
boolean flags.  Eventually, these may be parsed, but in the meantime, we use the
common convention of 0 meaning false or off and 1 for true or on.

This list includes parameters for the Enzo 2.0 release.

.. highlight:: none

Stopping Parameters
-------------------

``StopTime`` (external)
    This parameter specifies the time (in code units) when the
    calculation will halt. For cosmology simulations, this variable is
    automatically set by ``CosmologyFinalRedshift``. *No default.*
``StopCycle`` (external)
    The cycle (top grid timestep) at which the calculation stops. A
    value of zero indicates that this criterion is not be used.
    *Default: 100,000*
``StopFirstTimeAtLevel`` (external)
    Causes the simulation to immediately stop when a specified level is
    reached. Default value 0 (off), possible values are levels 1
    through maximum number of levels in a given simulation.
``NumberOfOutputsBeforeExit`` (external)
    After this many datadumps have been written, the code will exit.  If 
    set to 0 (default), this option will not be used.  Default: 0.
``StopCPUTime`` (external)
    Causes the simulation to stop if the wall time exceeds ``StopCPUTime``.
    The simulation will output if the wall time after the next
    top-level timestep will exceed ``StopCPUTime``, assuming that the wall
    time elapsed during a top-level timestep the same as the previous
    timestep. In units of seconds. Default: 2.592e6 (30 days)
``ResubmitOn`` (external)
    If set to 1, the simulation will stop if the wall time will exceed
    ``StopCPUTime`` within the next top-level timestep and run a shell
    script defined in ``ResubmitCommand`` that should resubmit the job
    for the user. Default: 0.
``ResubmitCommand`` (external)
    Filename of a shell script that creates a queuing (e.g. PBS)
    script from two arguments, the number of processors and parameter
    file.  This script is run by the root processor when stopping with
    ``ResubmitOn``. An example script can be found in
    input/resubmit.sh. Default: (null)

Initialization Parameters
-------------------------

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
6 	     Implosion
7 	     SedovBlast
8 	     KH Instability
9 	     2D/3D Noh Problem
10 	     :ref:`rotatingcylinder_param`
11 	     :ref:`radiatingshock_param`
12 	     :ref:`freeexpansion_param`
20 	     :ref:`zeldovichpancake_param`
21 	     :ref:`pressurelesscollapse_param`
22 	     :ref:`adiabaticexpansion_param`
23 	     :ref:`testgravity_param`
24 	     :ref:`sphericalinfall_param`
25 	     :ref:`testgravitysphere_param`
26 	     :ref:`gravityequilibriumtest_param`
27 	     :ref:`collapsetest_param`
28 	     TestGravityMotion
29 	     TestOrbit
30 	     :ref:`cosmologysimulation_param`
31 	     :ref:`galaxysimulation_param`
35 	     :ref:`shearingbox_param`
40 	     :ref:`supernovarestart_param`
50 	     :ref:`photontest_param`
60 	     Turbulence Simulation
61 	     Protostellar Collapse
62 	     :ref:`coolingtest_param`
101	     3D Collapse Test (hydro_rk)
102	     1D Spherical Collapse Test (hydro_rk)
106	     Hydro and MHD Turbulence Simulation (hydro_rk)
107 	     Put Sink from restart
200	     1D MHD Test
201	     2D MHD Test
202	     3D MHD Collapse Test
203	     MHD Turbulent Collapse Test
207	     Galaxy disk
208	     AGN disk
300	     Poisson solver test
400 	     Radiation-Hydrodynamics test 1 -- constant fields
401 	     Radiation-Hydrodynamics test 2 -- stream test
402 	     Radiation-Hydrodynamics test 3 -- pulse test
403 	     Radiation-Hydrodynamics test 4 -- grey Marshak test
404/405	     Radiation-Hydrodynamics test 5 -- radiating shock test
410/411	     Radiation-Hydrodynamics test 10/11 -- Static HI ionization
412 	     Radiation-Hydrodynamics test 12 -- HI ionization of a clump
413 	     Radiation-Hydrodynamics test 13 -- HI ionization of a steep region
414/415	     Radiation-Hydrodynamics test 14/15 -- Cosmological HI ionization
450-452	     Free-streaming radiation tests
============ ====================================

.. raw:: html

   <p></p>

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
``ShearingVelocityDirection`` (external)
    Select direction of shearing boundary. Default is x direction. Changing this is probably not a good idea.
``AngularVelocity`` (external)
    The value of the angular velocity in the shearing boundary.
    Default: 0.001
``VelocityGradient`` (external)
    The value of the per code length gradient in the angular velocity
    in the shearing boundary. Default: 1.0
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

Simulation Identifiers and UUIDs
--------------------------------

These parameters help to track, identify and group datasets. For reference,
`Universally Unique Identifiers
<http://en.wikipedia.org/wiki/Universally_Unique_Identifier>`_ (UUIDs) are
opaque identifiers using random 128-bit numbers, with an extremely low chance
of collision. (See :ref:`SimulationNamesAndIdentifiers` for a longer
description of these parameters.)

``MetaDataIdentifier`` (external)
    This is a character string without spaces (specifically, something
    that can be picked by "%s"), that can be defined in a parameter
    file, and will be written out in every following output, if it is
    found.
``MetaDataSimulationUUID`` (internal)
    A UUID that will be written out in all of the following outputs.
    Like ``MetaDataIdentifier``, an existing UUID will be kept, but if one
    is not found, and new one will be generated.
``MetaDataDatasetUUID`` (internal)
    A UUID created for each specific output.
``MetaDataRestartDatasetUUID`` (internal)
    If a ``MetaDataDatasetUUID`` UUID is found when the parameter file is
    read in, it will written to the following datasets. This is used to
    track simulations across restarts and parameter adjustments.
``MetaDataInitialConditionsUUID`` (internal)
    This is similar to ``MetaDataRestartDatasetUUID``, except it's used to
    track which initial conditions were used.

I/O Parameters
--------------

There are three ways to specify the frequency of outputs:
time-based, cycle-based (a cycle is a top-grid timestep), and, for
cosmology simulations, redshift-based. There is also a shortened
output format intended for visualization (movie format). Please
have a look at :ref:`controlling_data_output` for more information.

``dtDataDump`` (external)
    The time interval, in code units, between time-based outputs. A
    value of 0 turns off the time-based outputs. Default: 0
``CycleSkipDataDump`` (external)
    The number of cycles (top grid timesteps) between cycle-based
    outputs. Zero turns off the cycle-based outputs. Default: 0
``DataDumpName`` (external)
    The base file name used for both time and cycle based outputs.
    Default: data
``RedshiftDumpName`` (external)
    The base file name used for redshift-based outputs (this can be
    overridden by the ``CosmologyOutputRedshiftName`` parameter). Normally
    a four digit identification number is appended to the end of this
    name, starting from 0000 and incrementing by one for every output.
    This can be over-ridden by including four consecutive R's in the
    name (e.g. RedshiftRRRR) in which case the an identification number
    will not be appended but the four R's will be converted to a
    redshift with an implied decimal point in the middle (i.e. z=1.24
    becomes 0124). Default: RedshiftOutput
``CosmologyOutputRedshift[NNNN]`` (external)
    The time and cycle-based outputs occur regularly at constant
    intervals, but the redshift outputs are specified individually.
    This is done by the use of this statement, which sets the output
    redshift for a specific identification number (this integer is
    between 0000 and 9999 and is used in forming the name). So the
    statement ``CosmologyOutputRedshift[1] = 4.0`` will cause an output to
    be written out at z=4 with the name RedshiftOutput0001 (unless the
    base name is changed either with the previous parameter or the next
    one). This parameter can be repeated with different values for the
    number (NNNN) Default: none
``CosmologyOutputRedshiftName[NNNN]`` (external)
    This parameter overrides the parameter ``RedshiftOutputName`` for this
    (only only this) redshift output. Can be used repeatedly in the
    same manner as the previous parameter. Default: none
``OutputFirstTimeAtLevel`` (external)
    This forces Enzo to output when a given level is reached, and at
    every level thereafter. Default is 0 (off). User can usefully
    specify anything up to the maximum number of levels in a given
    simulation.
``FileDirectedOutput``
    If this parameter is set to 1, whenever the finest level has finished
    evolving Enzo will check for new signal files to output.  (See
    :ref:`force_output_now`.)  Default 1.
``XrayLowerCutoffkeV``, ``XrayUpperCutoffkeV``, ``XrayTableFileName`` (external)
    These parameters are used in 2D projections (``enzo -p ...``). The
    first two specify the X-ray band (observed at z=0) to be used, and
    the last gives the name of an ascii file that contains the X-ray
    spectral information. A gzipped version of this file good for
    bands within the 0.1 - 20 keV range is provided in the
    distribution in ``input/lookup_metal0.3.data``. If these
    parameters are specified, then the second field is replaced with
    integrated emissivity along the line of sight in units of 10\
    :sup:`-23` erg/cm\ :sup:`2`/s. Default: ``XrayLowerCutoffkeV =
    0.5``, ``XrayUpperCutoffkeV = 2.5``.
``ExtractFieldsOnly`` (external)
    Used for extractions (enzo -x ...) when only field data are needed
    instead of field + particle data. Default is 1 (TRUE).
``dtRestartDump``
    Reserved for future use.
``dtHistoryDump``
    Reserved for future use.
``CycleSkipRestartDump``
    Reserved for future use.
``CycleSkipHistoryDump``
    Reserved for future use.
``RestartDumpName``
    Reserved for future use.
``HistoryDumpName``
    Reserved for future use.
``ParallelRootGridIO`` (external)
    Normally for the mpi version, the root grid is read into the root
    processor and then partitioned to separate processors using communication.
    However, for
    very large root grids (e.g. 512\ :sup:`3`\ ), the root processor
    may not have enough memory. If this toggle switch is set on (i.e.
    to the value 1), then each processor reads its own section of the
    root grid. More I/O is required (to split up the grids and
    particles), but it is more balanced in terms of memory.
    ``ParallelRootGridIO`` and ``ParallelParticleIO`` MUST be set to 1 (TRUE)
    for runs involving > 64 cpus! Default: 0 (FALSE). See also ``Unigrid``
    below.
``Unigrid`` (external)
    This parameter should be set to 1 (TRUE) for large cases--AMR as
    well as non-AMR--where the root grid is 512\ :sup:`3`\  or larger.
    This prevents initialization under subgrids at start up, which is
    unnecessary in cases with simple non-nested initial conditions.
    Unigrid must be set to 0 (FALSE) for cases with nested initial
    conditions. Default: 0 (FALSE). See also ``ParallelRootGridIO`` above.
``UnigridTranspose`` (external)
    This parameter governs the fast FFT bookkeeping for Unigrid runs.
    Does not work with isolated gravity. Default: 0 (FALSE). See also
    ``Unigrid`` above.
``OutputTemperature`` (external)
    Set to 1 if you want to output a temperature field in the datasets.
    Always 1 for cosmology simulations. Default: 0.
``OutputCoolingTime`` (external)
    Set to 1 if you want to output the cooling time in the datasets.
    Default: 0.
``OutputSmoothedDarkMatter`` (external)
    Set to 1 if you want to output a dark matter density field,
    smoothed by an SPH kernel. Set to 2 to also output smoothed dark
    matter velocities and velocity dispersion. Set to 0 to turn off.
    Default: 0.
``OutputGriddedStarParticle`` (external)
    Set to 1 or 2 to write out star particle data gridded onto mesh.
    This will be useful e.g. if you have lots of star particles in a
    galactic scale simulation. 1 will output just
    ``star_particle_density``; and 2 will dump
    ``actively_forming_stellar_mass_density``, ``SFR_density``, etc.
    Default: 0.
``VelAnyl`` (external)
    Set to 1 if you want to output the divergence and vorticity of
    velocity. Works in 2D and 3D.
``BAnyl`` (external)
    Set to 1 if you want to output the divergence and vorticity of
    ``Bfield``. Works in 2D and 3D.
``SmoothedDarkMatterNeighbors`` (external)
    Number of nearest neighbors to smooth dark matter quantities over.
    Default: 32.

.. _streaming_param:

Streaming Data Format
~~~~~~~~~~~~~~~~~~~~~

``NewMovieLeftEdge``, ``NewMovieRightEdge`` (external)
    These two parameters control the region for which the streaming
    data are written. Default: ``DomainLeftEdge`` and ``DomainRightEdge``.
``MovieSkipTimestep`` (external)
    Controls how many timesteps on a level are skipped between outputs
    in the streaming data. Streaming format is off if this equals
    ``INT_UNDEFINED``. Default: ``INT_UNDEFINED``
``Movie3DVolume`` (external)
    Set to 1 to write streaming data as 3-D arrays. This should always
    be set to 1 if using the streaming format. A previous version had
    2D maximum intensity projections, which now defunct. Default: 0.
``MovieVertexCentered`` (external)
    Set to 1 to write the streaming data interpolated to vertices. Set
    to 0 for cell-centered data. Default: 0.
``NewMovieDumpNumber`` (internal)
    Counter for streaming data files. This should equal the cycle
    number.
``MovieTimestepCounter`` (internal)
    Timestep counter for the streaming data files.
``MovieDataField`` (external)
    A maximum of 6 data fields can be written in the streaming format.
    The data fields are specified by the array element of
    BaryonField, i.e. 0 = Density, 7 = HII
    Density. For writing temperature, a special value of 1000 is used.
    This should be improved to be more transparent in which fields will
    be written. Any element that equals ``INT_UNDEFINED`` indicates no
    field will be written. Default: ``INT_UNDEFINED`` x 6
``NewMovieParticleOn`` (external)
    Set to 1 to write all particles in the grids. Set to 2 to write
    ONLY particles that aren't dark matter, e.g. stars. Set to 3/4 to
    write ONLY particles that aren't dark matter into a file separate
    from the grid info. (For example, ``MoviePackParticle_P000.hdf5``,
    etc. will be the file name; this will be very helpful in speeding
    up the access to the star particle data, especially for the
    visualization or for the star particle. See ``AMRH5writer.C``) Set to 0
    for no particle output. Default: 0.

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

       1 - refine by slope		       6  - refine by Jeans length
       2 - refine by baryon mass	       7  - refine if (cooling time < cell width/sound speed)
       3 - refine by shocks 		       11 - refine by resistive length
       4 - refine by particle mass	       12 - refine by defined region "MustRefineRegion"
       5 - refine by baryon overdensity	       13 - refine by metallicity
       	  (currently disabled)
       101 - avoid refinement in regions
             defined in "AvoidRefineRegion"

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

Hydrodynamic Parameters
-----------------------

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

Magnetohydrodynamic Parameters
------------------------------

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

Cosmology Parameters
--------------------

``ComovingCoordinates`` (external)
    Flag (1 - on, 0 - off) that determines if comoving coordinates are
    used or not. In practice this turns on or off the entire cosmology
    machinery. Default: 0
``CosmologyFinalRedshift`` (external)
    This parameter specifies the redshift when the calculation will
    halt. Default: 0.0
``CosmologyOmegaMatterNow`` (external)
    This is the contribution of all non-relativistic matter (including
    HDM) to the energy density at the current epoch (z=0), relative to
    the value required to marginally close the universe. It includes
    dark and baryonic matter. Default: 0.279
``CosmologyOmegaLambdaNow`` (external)
    This is the contribution of the cosmological constant to the energy
    density at the current epoch, in the same units as above. Default:
    0.721
``CosmologyComovingBoxSize`` (external)
    The size of the volume to be simulated in Mpc/h (at z=0). Default:
    64.0
``CosmologyHubbleConstantNow`` (external)
    The Hubble constant at z=0, in units of 100 km/s/Mpc. Default:
    0.701
``CosmologyInitialRedshift`` (external)
    The redshift for which the initial conditions are to be generated.
    Default: 20.0
``CosmologyMaxExpansionRate`` (external)
    This float controls the timestep so that cosmological terms are
    accurate followed. The timestep is constrained so that the relative
    change in the expansion factor in a step is less than this value.
    Default: 0.01
``CosmologyCurrentRedshift`` (information only)
    This is not strictly speaking a parameter since it is never
    interpreted and is only meant to provide information to the user.
    Default: n/a

Gravity Parameters
------------------

``TopGridGravityBoundary`` (external)
    A single integer which specified the type of gravitational boundary
    conditions for the top grid. Possible values are 0 for periodic and
    1 for isolated (for all dimensions). The isolated boundary
    conditions have not been tested recently, so caveat emptor.
    Default: 0
``SelfGravity`` (external)
    This flag (1 - on, 0 - off) indicates if the baryons and particles
    undergo self-gravity.
``GravitationalConstant`` (external)
    This is the gravitational constant to be used in code units. For cgs units it
    should be 4\*pi\*G. For cosmology, this value must be 1 for the
    standard units to hold. A more detailed decription can be found at :ref:`EnzoInternalUnits`. Default: 4\*pi.
``GreensFunctionMaxNumber`` (external)
    The Green's functions for the gravitational potential depend on the
    grid size, so they are calculated on a as-needed basis. Since they
    are often re-used, they can be cached. This integer indicates the
    number that can be stored. They don't take much memory (only the
    real part is stored), so a reasonable number is 100. [Ignored in
    current version]. Default: 1
``GreensFunctionMaxSize``
    Reserved for future use.
``S2ParticleSize`` (external)
    This is the gravitational softening radius, in cell widths, in
    terms of the S2 particle described by Hockney and Eastwood in their
    book Computer Simulation Using Particles. A reasonable value is
    3.0. [Ignored in current version]. Default: 3.0
``GravityResolution`` (external)
    This was a mis-guided attempt to provide the capability to increase
    the resolution of the gravitational mesh. In theory it still works,
    but has not been recently tested. Besides, it's just not a good
    idea. The value (a float) indicates the ratio of the gravitational
    cell width to the baryon cell width. [Ignored in current version].
    Default: 1
``PotentialIterations`` (external)
    Number of iterations to solve the potential on the subgrids. Values
    less than 4 sometimes will result in slight overdensities on grid
    boundaries. Default: 4.
``BaryonSelfGravityApproximation`` (external)
    This flag indicates if baryon density is derived in a strange,
    expensive but self-consistent way (0 - off), or by a completely
    reasonable and much faster approximation (1 - on). This is an
    experiment gone wrong; leave on. Well, actually, it's important for
    very dense structures as when radiative cooling is turned on, so
    set to 0 if using many levels and radiative cooling is on [ignored
    in current version]. Default: 1
``MaximumGravityRefinementLevel`` (external)
    This is the lowest (most refined) depth that a gravitational
    acceleration field is computed. More refined levels interpolate
    from this level, provided a mechanism for instituting a minimum
    gravitational smoothing length. Default: ``MaximumRefinementLevel``
    (unless ``HydroMethod`` is ZEUS and radiative cooling is on, in which
    case it is ``MaximumRefinementLevel`` - 3).
``MaximumParticleRefinementLevel`` (external)
    This is the level at which the dark matter particle contribution to
    the gravity is smoothed. This works in an inefficient way (it
    actually smoothes the particle density onto the grid), and so is
    only intended for highly refined regions which are nearly
    completely baryon dominated. It is used to remove the discreteness
    effects of the few remaining dark matter particles. Not used if set
    to a value less than 0. Default: -1

External Gravity Source
~~~~~~~~~~~~~~~~~~~~~~~~

These parameters set-up an external static background gravity source that is
added to the acceleration field for the baryons and particles.

``PointSourceGravity`` (external)
    This parameter indicates that there is to be a
    (constant) gravitational field with a point source profile (``PointSourceGravity`` =
    1) or NFW profile (``PointSourceGravity`` = 2). Default: 0
``PointSourceGravityConstant`` (external)
    If ``PointSourceGravity`` = 1, this is the magnitude of the point
    source acceleration at a distance of 1
    length unit (i.e. GM in code units). If ``PointSourceGravity`` =
    2, then it takes the mass of the dark matter halo in CGS
    units. ``ProblemType`` = 31 (galaxy disk simulation) automatically calculates
    values for ``PointSourceGravityConstant`` and
    ``PointSourceGravityCoreRadius``. Default: 1
``PointSourceGravityCoreRadius`` (external)
    For ``PointSourceGravity`` = 1, this is the radius inside which
    the acceleration field is smoothed in code units. With ``PointSourceGravity`` =
    2, it is the scale radius, rs, in CGS units (see Navarro, Frank & White,
    1997). Default: 0
``PointSourceGravityPosition`` (external)
    If the ``PointSourceGravity`` flag is turned on, this parameter
    specifies the center of the point-source gravitational field in
    code units. Default: 0 0 0
``ExternalGravity`` (external)
   This fulfills the same purpose as ``PointSourceGravity`` but is
   more aptly named. Currently, it has only a single option
   ``ExternalGravity = 1`` which turns on an alternative
   implementation of the NFW profile. The profile properties are
   defined via the parameters ``HaloCentralDensity``, ``HaloConcentration`` and ``HaloVirialRadius``. Default: 0 
``ExternalGravityDensity`` 
   Reserved for future use.
``ExternalGravityRadius``
   Reserved for future use.
``UniformGravity`` (external)
    This flag (1 - on, 0 - off) indicates if there is to be a uniform
    gravitational field. Default: 0
``UniformGravityDirection`` (external)
    This integer is the direction of the uniform gravitational field: 0
    - along the x axis, 1 - y axis, 2 - z axis. Default: 0
``UniformGravityConstant`` (external)
    Magnitude (and sign) of the uniform gravitational acceleration.
    Default: 1

Particle Parameters
-------------------

``ParticleBoundaryType`` (external)
    The boundary condition imposed on particles. At the moment, this
    parameter is largely ceremonial as there is only one type
    implemented: periodic, indicated by a 0 value. Default: 0
``ParticleCourantSafetyNumber`` (external)
    This somewhat strangely named parameter is the maximum fraction of
    a cell width that a particle is allowed to travel per timestep
    (i.e. it is a constant on the timestep somewhat along the lines of
    it's hydrodynamic brother). Default: 0.5
``NumberOfParticles`` (obsolete)
    Currently ignored by all initializers, except for TestGravity and
    TestGravitySphere where it is the number of test points. Default: 0
``NumberOfParticleAttributes`` (internal)
    It is set to 3 if either ``StarParticleCreation`` or
    ``StarParticleFeedback`` is set to 1 (TRUE). Default: 0
``AddParticleAttributes`` (internal)
    If set to 1, additional particle attributes will be added and
    zeroed. This is handy when restarting a run, and the user wants to
    use star formation afterwards. Default: 0.
``ParallelParticleIO`` (external)
    Normally, for the mpi version, the particle data are read into the
    root processor and then distributed to separate processors.
    However, for very large number of particles, the root processor may
    not have enough memory. If this toggle switch is set on (i.e. to
    the value 1), then Ring i/o is turned on and each processor reads
    its own part of the particle data. More I/O is required, but it is
    more balanced in terms of memory. ``ParallelRootGridIO`` and
    ``ParallelParticleIO`` MUST be set for runs involving > 64 cpus!
    Default: 0 (FALSE).
``ParticleSplitterIterations`` (external)
    Set to 1 to split particles into 13 particles (= 12 children+1
    parent, Kitsionas & Whitworth (2002)). This should be ideal for
    setting up an low-resolution initial condition for a relatively low
    computational cost, running it for a while, and then restarting it
    for an extremely high-resolution simulation in a focused region.
    Currently it implicitly assumes that only DM (type=1) and
    conventional star particles (type=2) inside the ``RefineRegion`` get
    split. Other particles, which usually become Star class objects,
    seem to have no reason to be split. Default: 0
``ParticleSplitterChildrenParticleSeparation`` (external)
    This is the spacing between the child particles placed on a
    hexagonal close-packed (HCP) array. In the unit of a cell size
    which the parent particle resides in. Default: 1.0

Parameters for Additional Physics
---------------------------------

``RadiativeCooling`` (external)
    This flag (1 - on, 0 - off) controls whether or not a radiative
    cooling module is called for each grid. There are currently several
    possibilities, controlled by the value of another flag. See :ref:`cooling` 
    for more information on the various cooling methods.  Default: 0
    
    -  If the ``MultiSpecies`` flag is off, then equilibrium cooling is
       assumed and one of the following two will happen. If the parameter
       ``GadgetCooling`` is set to 1, the primordial equilibrium code is
       called (see below). If ``GadgetCooling`` is set to 0, a file called
       ``cool_rates.in`` is read to set a cooling curve. This file consists
       of a set of temperature and the associated cgs cooling rate; a
       sample compute with a metallicity Z=0.3 Raymond-Smith code is
       provided in ``input/cool_rates.in``. This has a cutoff at 10000 K
       (Sarazin & White 1987). Another choice will be
       ``input/cool_rates.in_300K`` which goes further down to 300 K (Rosen
       & Bregman 1995).
    -  If the ``MultiSpecies`` flag is on, then the cooling rate is
       computed directly by the species abundances. This routine (which
       uses a backward differenced multi-step algorithm) is borrowed
       from the Hercules code written by Peter Anninos and Yu Zhang,
       featuring rates from Tom Abel. Other varieties of cooling are
       controlled by the ``MetalCooling`` parameter, as discused below.

``GadgetCooling`` (external)
    This flag (1 - on, 0 - off) turns on (when set to 1) a set of
    routines that calculate cooling rates based on the assumption of a
    six-species primordial gas (H, He, no H2 or D) in equilibrium, and
    is valid for temperatures greater than 10,000 K. This requires the
    file ``TREECOOL`` to execute. Default: 0
``MetalCooling`` (external)
    This flag (0 - off, 1 - metal cooling from Glover & Jappsen 2007,
    2 - Cen et al (1995), 3 - Cloudy cooling from Smith, Sigurdsson, &
    Abel 2008) turns on metal cooling for runs that track
    metallicity. Option 1 is valid for temperatures between 100 K and
    10\ :sup:`8`\ K because it considers fine-structure line emission
    from carbon, oxygen, and silicon and includes the additional metal
    cooling rates from Sutherland & Dopita (1993). Option 2 is only
    valid for temperatures above 10\ :sup:`4`\ K. Option 3 uses
    multi-dimensional tables of heating/cooling values created with
    Cloudy and optionally coupled to the ``MultiSpecies``
    chemistry/cooling solver. This method is valid from 10 K to 10\
    :sup:`8`\ K. See the Cloudy Cooling parameters below.  Default: 0.
``MetalCoolingTable`` (internal)
    This field contains the metal cooling table required for
    ``MetalCooling`` option 1. In the top level directory input/, there are
    two files ``metal_cool.dat`` and ``metal_cool_pop3.dat`` that consider
    metal cooling for solar abundance and abundances from
    pair-instability supernovae, respectively. In the same directory,
    one can find an IDL routine (``make_Zcool_table.pro``) that generates
    these tables. Default: ``metal_cool.dat``
``MultiSpecies`` (external)
    If this flag (1, 2, 3- on, 0 - off) is on, then the code follows
    not just the total density, but also the ionization states of
    Hydrogen and Helium. If set to 2, then a nine-species model
    (including H2, H2+ and H-) will be computed, otherwise only six
    species are followed (H, H+, He, He+, He++, e-). If set to 3, then
    a 12 species model is followed, including D, D+ and HD. This
    routine, like the last one, is based on work done by Abel, Zhang
    and Anninos. Default: 0
``GadgetEquilibriumCooling`` (external)
    An implementation of the ionization equilibrium cooling code used
    in the GADGET code which includes both radiative cooling and a
    uniform metagalactic UV background specified by the ``TREECOOL`` file
    (in the ``amr_mpi/exe`` directory). When this parameter is turned on,
    ``MultiSpecies`` and ``RadiationFieldType`` are forced to 0 and
    ``RadiativeCooling`` is forced to 1.
    [Not in public release version]
``PhotoelectricHeating`` (external)
    If set to be 1, Gamma_pe = 5.1e-26 erg/s will be added uniformly
    to the gas without any shielding (Tasker & Bryan 2008). At the
    moment this is still experimental. Default: 0
``MultiMetals`` (external)
    This was added so that the user could turn on or off additional
    metal fields - currently there is the standard metallicity field
    (Metal_Density) and two additional metal fields (Z_Field1 and
    Z_Field2). Acceptable values are 1 or 0, Default: 0 (off).

.. _cloudy_cooling:

Cloudy Cooling
~~~~~~~~~~~~~~

Cloudy cooling from Smith, Sigurdsson, & Abel (2008) interpolates
over tables of precomputed cooling data. Cloudy cooling is turned
on by setting ``MetalCooling`` to 3. ``RadiativeCooling`` must also be set
to 1. Depending on the cooling data used, it can be coupled with
``MultiSpecies`` = 1, 2, or 3 so that the metal-free cooling comes from
the ``MultiSpecies`` machinery and the Cloudy tables provide only the
metal cooling. Datasets range in dimension from 1 to 5. Dim 1:
interpolate over temperature. Dim 2: density and temperature. Dim
3: density, metallicity, and temperature. Dim 4: density,
metallicity, electron fraction, and temperature. Dim 5: density,
metallicity, electron fraction, spectral strength, and temperature.
See Smith, Sigurdsson, & Abel (2008) for more information on
creating Cloudy datasets.

``CloudyCoolingGridFile`` (external)
    A string specifying the path to the Cloudy cooling dataset.

``IncludeCloudyHeating`` (external)
    An integer (0 or 1) specifying whether the heating rates are to be
    included in the calculation of the cooling. Some Cloudy datasets
    are made with the intention that only the cooling rates are to be
    used. Default: 0 (off).

``IncludeCloudyMMW`` (external)
    An integer (0 or 1) specifying whether the additional mean
    molecular weight contributed by the metals be used in the
    conversion from internal energy to temperature. These values will
    come from the Cloudy dataset. For metallicities less than solar,
    this addition will be negligible. Default: 0 (off).

``CMBTemperatureFloor`` (external)
    An integer (0 or 1) specifying whether a temperature floor is
    created at the temperature of the cosmic microwave background
    (T\ :sub:`CMB`\  = 2.72 (1 + z) K). This is accomplished in the
    code by subtracting the cooling rate at T\ :sub:`CMB`\  such that
    Cooling = Cooling(T) - Cooling(T\ :sub:`CMB`\ ). Default: 1 (on).

``CloudyMetallicityNormalization`` (external)
    A float value used in the conversion of metal density into
    metallicity. This value will change depending on the specific
    abundance patterns used to make the Cloudy dataset. The value of
    this factor is calculated as the sum of (A\ :sub:`i`\  \*
    m\ :sub:`i`\ ) over all elements i heavier than He, where
    A\ :sub:`i`\  is the solar number abundance relative to H and
    m\ :sub:`i`\  is the atomic mass. For the solar abundance pattern
    from the latest version of Cloudy, using all metals through Zn,
    this value is 0.018477. Default: 0.018477.

``CloudyElectronFractionFactor`` (external)
    A float value to account for additional electrons contributed by
    metals. This is only used with Cloudy datasets with dimension
    greater than or equal to 4. The value of this factor is calculated
    as the sum of (A\ :sub:`i`\  \* i) over all elements i heavier than
    He, where A\ :sub:`i`\  is the solar number abundance relative to
    H. For the solar abundance pattern from the latest version of
    Cloudy, using all metals through Zn, this value is 9.153959e-3.
    Default: 9.153959e-3.

Inline Halo Finding
~~~~~~~~~~~~~~~~~~~

Enzo can find dark matter (sub)halos on the fly with a
friends-of-friends (FOF) halo finder and a subfind method,
originally written by Volker Springel. All output files will be
written in the directory FOF/.

``InlineHaloFinder`` (external)
    Set to 1 to turn on the inline halo finder. Default: 0.
``HaloFinderSubfind`` (external)
    Set to 1 to find subhalos inside each dark matter halo found in the
    friends-of-friends method. Default: 0.
``HaloFinderOutputParticleList`` (external)
    Set to 1 to output a list of particle positions and IDs for each
    (sub)halo. Written in HDF5. Default: 0.
``HaloFinderMinimumSize`` (external)
    Minimum number of particles to be considered a halo. Default: 50.
``HaloFinderLinkingLength`` (external)
    Linking length of particles when finding FOF groups. In units of
    cell width of the finest static grid, e.g. unigrid -> root cell
    width. Default: 0.1.
``HaloFinderCycleSkip`` (external)
    Find halos every N\ :sup:`th`\  top-level timestep, where N is this
    parameter. Not used if set to 0. Default: 3.
``HaloFinderTimestep`` (external)
    Find halos every dt = (this parameter). Only evaluated at each
    top-level timestep. Not used if negative. Default: -99999.0
``HaloFinderLastTime`` (internal)
    Last time of a halo find. Default: 0.

Inline Python
~~~~~~~~~~~~~

``PythonSubcycleSkip`` (external)
    The number of times Enzo should reach the bottom of the hierarchy
    before exposing its data and calling Python. Only works with
    python-yes in compile settings.

.. _StarParticleParameters:

Star Formation and Feedback Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For details on each of the different star formation methods available in Enzo see :ref:`star_particles`.


``StarParticleCreation`` (external)
    This parameter is bitwise so that multiple types of star formation
    routines can be used in a single simulation. For example if methods
    1 and 3 are desired, the user would specify 10 (2\ :sup:`1`\  +
    2\ :sup:`3`\ ), or if methods 1, 4 and 7 are wanted, this would be
    146 (2\ :sup:`1`\  + 2\ :sup:`4`\  + 2\ :sup:`7`\ ). Default: 0
    
    ::

	0 - Cen & Ostriker (1992)
	1 - Cen & Ostriker (1992) with stocastic star formation
	2 - Global Schmidt Law / Kravstov et al. (2003)
	3 - Population III stars / Abel, Wise & Bryan (2007)
	4 - Sink particles: Pure sink particle or star particle with wind feedback depending on 
	    choice for HydroMethod / Wang et al. (2009)
	5 - Radiative star clusters  / Wise & Cen (2009)
	6 - [reserved]
	7 - Cen & Ostriker (1992) with no delay in formation
	8 - Springel & Hernquist (2003)
	9 - Massive Black Hole (MBH) particles insertion by hand / Kim et al. (2010)
	10 - Population III stellar tracers  

``StarParticleFeedback`` (external)
    This parameter works the same way as ``StarParticleCreation`` but only
    is valid for ``StarParticleCreation`` = 0, 1, 2, 7 and 8 because methods 3, 5 and 9
    use the radiation transport module and ``Star_*.C`` routines to
    calculate the feedback, 4 has explicit feedback and 10 does not use feedback. Default: 0.

``StarFeedbackDistRadius`` (external)
    If this parameter is greater than zero, stellar feedback will be
    deposited into the host cell and neighboring cells within this
    radius.  This results in feedback being distributed to a cube with
    a side of ``StarFeedbackDistRadius+1``. It is in units of cell
    widths of the finest grid which hosts the star particle.  Only
    implemented for ``StarFeedbackCreation`` = 0 or 1 with ``StarParticleFeedback`` =  1. (If ``StarParticleFeedback`` = 0, stellar feedback is only deposited into the cell in which the star particle lives).  Default: 0.

``StarFeedbackDistCellStep`` (external)
    In essence, this parameter controls the shape of the volume where
    the feedback is applied, cropping the original cube.  This volume
    that are within ``StarFeedbackDistCellSteps`` cells from the host
    cell, counted in steps in Cartesian directions, are injected with
    stellar feedback.  Its maximum value is ``StarFeedbackDistRadius``
    * ``TopGridRank``.  Only implemented for ``StarFeedbackCreation`` = 0
    or 1.  See :ref:`distributed_feedback` for an illustration.
    Default: 0.

Normal Star Formation
^^^^^^^^^^^^^^^^^^^^^

The parameters below are considered in ``StarParticleCreation`` method
0, 1, 2, 7 and 8.

``StarMakerOverDensityThreshold`` (external)
    The overdensity threshold in code units (for cosmological simulations, note that code units are relative to the total mean density, not
    just the dark matter mean density) before star formation will be
    considered. For ``StarParticleCreation`` = 7 in cosmological
    simulations, however, ``StarMakerOverDensityThreshold`` should be in
    particles/cc, so it is not the ratio with respect to the
    ``DensityUnits`` (unlike most other
    star_makers). This way one correctly represents the Jeans
    collapse and molecular cloud scale physics even in cosmological
    simulations. Default: 100
``StarMakerSHDensityThreshold`` (external)
    The critical density of gas used in Springel & Hernquist star
    formation ( \\rho_{th} in the paper) used to determine the star
    formation timescale in units of g cm\ :sup:`-3`\ . Only valid for ``StarParticleCreation`` = 8. Default: 7e-26.
``StarMakerMassEfficiency`` (external)
    The fraction of identified baryonic mass in a cell
    (Mass\*dt/t_dyn) that is converted into a star particle. Default:
    1
``StarMakerMinimumMass`` (external)
    The minimum mass of star particle, in solar masses. Note however,
    the star maker algorithm 2 has a (default off) "stochastic" star formation
    algorithm that will, in a pseudo-random fashion, allow star
    formation even for very low star formation rates. It attempts to do
    so (relatively successfully according to tests) in a fashion that
    conserves the global average star formation rate. Default: 1e9
``StarMakerMinimumDynamicalTime`` (external)
    When the star formation rate is computed, the rate is proportional
    to M_baryon \* dt/max(t_dyn, t_max) where t_max is this
    parameter. This effectively sets a limit on the rate of star
    formation based on the idea that stars have a non-negligible
    formation and life-time. The unit is years. Default: 1e6
``StarMassEjectionFraction`` (external)
    The mass fraction of created stars which is returned to the gas
    phase. Default: 0.25
``StarMetalYield`` (external)
    The mass fraction of metals produced by each unit mass of stars
    created (i.e. it is multiplied by mstar, not ejected). Default:
    0.02
``StarEnergyToThermalFeedback`` (external)
    The fraction of the rest-mass energy of the stars created which is
    returned to the gas phase as thermal energy. Default: 1e-5
``StarEnergyToStellarUV`` (external)
    The fraction of the rest-mass energy of the stars created which is
    returned as UV radiation with a young star spectrum. This is used when calculating the radiation background. Default: 3e-6
``StarEnergyToQuasarUV`` (external)
    The fraction of the rest-mass energy of the stars created which is returned as UV radiation with a quasar spectrum. This is used when calculating the radiation background. Default: 5e-6

Population III Star Formation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters below are considered in ``StarParticleCreation`` method 3.

``PopIIIStarMass`` (external)
    Stellar mass of Population III stars created in
    ``StarParticleCreation`` method 3. Units of solar masses. The
    luminosities and supernova energies are calculated from Schaerer
    (2002) and Heger & Woosley (2002), respectively.
``PopIIIBlackHoles`` (external)
    Set to 1 to create black hole particles that radiate in X-rays for
    stars that do not go supernova (< 140 solar masses and > 260 solar
    masses). Default: 0.
``PopIIIBHLuminosityEfficiency`` (internal)
    The radiative efficiency in which the black holes convert accretion
    to luminosity. Default: 0.1.
``PopIIIOverDensityThreshold`` (internal)
    The overdensity threshold (relative to the total mean density)
    before Pop III star formation will be considered. Default: 1e6.
``PopIIIH2CriticalFraction`` (internal)
    The H_2 fraction threshold before Pop III star formation will be
    considered. Default: 5e-4.
``PopIIIMetalCriticalFraction`` (internal)
    The metallicity threshold (relative to gas density, not solar)
    before Pop III star formation will be considered. Note: this should
    be changed to be relative to solar! Default: 1e-4.
``PopIIISupernovaRadius`` (internal)
    If the Population III star will go supernova (140<M<260 solar
    masses), this is the radius of the sphere to inject the supernova
    thermal energy at the end of the star's life. Units are in parsecs.
    Default: 1.
``PopIIISupernovaUseColour`` (internal)
    Set to 1 to trace the metals expelled from supernovae. Default: 0.

Radiative Star Cluster Star Formation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters below are considered in ``StarParticleCreation`` method 5.

``StarClusterUseMetalField`` (internal)
    Set to 1 to trace ejecta from supernovae. Default: 0.
``StarClusterMinDynamicalTime`` (internal)
    When determining the size of a star forming region, one method is
    to look for the sphere with an enclosed average density that
    corresponds to some minimum dynamical time. Observations hint that
    this value should be a few million years. Units are in years.
    Default: 1e7.
``StarClusterIonizingLuminosity`` (internal)
    The specific luminosity of the stellar clusters. In units of
    ionizing photons per solar mass. Default: 1e47.
``StarClusterSNEnergy`` (internal)
    The specific energy injected into the gas from supernovae in the
    stellar clusters. In units of ergs per solar mass. Default: 6.8e48
    (Woosley & Weaver 1986).
``StarClusterSNRadius`` (internal)
    This is the radius of the sphere to inject the supernova thermal
    energy in stellar clusters. Units are in parsecs. Default: 10.
``StarClusterFormEfficiency`` (internal)
    Fraction of gas in the sphere to transfer from the grid to the star
    particle. Recall that this sphere has a minimum dynamical time set
    by ``StarClusterMinDynamicalTime``. Default: 0.1.
``StarClusterMinimumMass`` (internal)
    The minimum mass of a star cluster particle before the formation is
    considered. Units in solar masses. Default: 1000.
``StarClusterCombineRadius`` (internal)
    It is possible to merge star cluster particles together within this
    specified radius. Units in parsecs. This is probably not necessary
    if ray merging is used. Originally this was developed to reduce the
    amount of ray tracing involved from galaxies with hundreds of these
    radiating particles. Default: 10.

Massive Black Hole Particle Formation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters below are considered in StarParticleCreation method 9.

``MBHInsertLocationFilename`` (external)
    The mass and location of the MBH particle that has to be inserted.
    For example, the content of the file should be in the following
    form. For details, see ``mbh_maker.src``. Default:
    ``mbh_insert_location.in``
    ::

        #order: MBH mass (in Ms), MBH location[3], MBH creation time
        100000.0      0.48530579      0.51455688      0.51467896      0.0

.. _radiation_backgrounds:

Background Radiation Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``RadiationFieldType`` (external)
    This integer parameter specifies the type of radiation field that
    is to be used. Except for ``RadiationFieldType`` = 9, which should
    be used with ``MultiSpecies`` = 2, UV backgrounds can currently
    only be used with ``MultiSpecies`` = 1 (i.e. no molecular H
    support). The following values are used. Default: 0

   ::
  
     1. Haardt & Madau spectrum with q_alpha=1.5
     2. Haardt & Madau spectrum with q_alpha = 1.8
     3. reserved for experimentation
     4. H&M spectrum (q_alpha=1.5. supplemented with an X-ray Compton heating
         background from Madau & Efstathiou (see astro-ph/9902080)
     9. a constant molecular H2 photo-dissociation rate
     10. internally computed radiation field using the algorithm of Cen & Ostriker
     11. same as previous, but with very, very simple optical shielding fudge
     12. Haardt & Madau spectrum with q_alpha=1.57

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
``RadiationFieldRedshift`` (external)
    This parameter specifies the redshift at which the radiation field
    is calculated.  Default: 0
``RadiationShield`` (external)
    This parameter specifies whether the user wants to employ
    approximate radiative-shielding. This parameter will be
    automatically turned on when RadiationFieldType is set to 11. See
    ``calc_photo_rates.src``. Default: 0
``RadiationRedshiftOn`` (external) The redshift at which the UV 
    background turns on. Default: 7.0.
``RadiationRedshiftFullOn`` (external) The redshift at which the UV
    background is at full strength.  Between z =
    ``RadiationRedshiftOn`` and z = ``RadiationRedshiftFullOn``, the 
    background is gradually ramped up to full strength. Default: 6.0.
``RadiationRedshiftDropOff`` (external) The redshift at which the 
    strength of the UV background is begins to gradually reduce,
    reaching zero by ``RadiationRedshiftOff``. Default: 0.0.
``RadiationRedshiftOff`` (external) The redshift at which the UV 
    background is fully off. Default: 0.0.
``AdjustUVBackground`` (external)
    Add description. Default: 1.
``SetUVAmplitude`` (external)
    Add description. Default: 1.0.
``SetHeIIHeatingScale`` (external)
    Add description. Default: 1.8.
``RadiationSpectrumSlope`` (external)
    Add description. Default: 1.5.

Minimum Pressure Support Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Radiative Transfer (Ray Tracing) Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``RadiativeTransfer`` (external)
    Set to 1 to turn on the adaptive ray tracing following Abel, Wise &
    Bryan 2007. Note that Enzo must be first recompiled after setting
    ``make photon-yes``. Default: 0.
``RadiativeTransferRadiationPressure`` (external)
    Set to 1 to turn on radiation pressure created from absorbed photon
    packages. Default: 0
``RadiativeTransferInitialHEALPixLevel`` (internal)
    Chooses how many rays are emitted from radiation sources. The
    number of rays in Healpix are given through # =
    12x4\ :sup:`level`\ . Default: 3.
``RadiativeTransferRaysPerCell`` (external)
    Determines the accuracy of the scheme by giving the minimum number
    of rays to cross cells. The more the better (slower). Default: 5.1.
``RadiativeTransferSourceRadius`` (external)
    The radius at which the photons originate from the radiation
    source. A positive value results in a radiating sphere. Default: 0.
``RadiativeTransferPropagationRadius`` (internal)
    The maximum distance a photon package can travel in one timestep.
    Currently unused. Default: 0.
``RadiativeTransferPropagationSpeed`` (internal)
    The fraction of the speed of light at which the photons travel.
    Default: 1.
``RadiativeTransferCoupledRateSolver`` (internal)
    Set to 1 to calculate the new ionization fractions and gas energies
    after every radiative transfer timestep. This option is highly
    recommended to be kept on. If not, ionization fronts will propagate too
    slowly. Default: 1.
``RadiativeTransferOpticallyThinH2`` (external)
    Set to 1 to include an optically-thin H_2 dissociating
    (Lyman-Werner) radiation field. Only used if ``MultiSpecies`` > 1. If
    ``MultiSpecies`` > 1 and this option is off, the Lyman-Werner radiation
    field will be calculated with ray tracing. Default: 1.
``RadiativeTransferSplitPhotonPackage`` (internal)
    Once photons are past this radius, they can no longer split. In
    units of kpc. If this value is negative (by default), photons can
    always split. Default: ``FLOAT_UNDEFINED``.
``RadiativeTransferPhotonEscapeRadius`` (internal)
    The number of photons that pass this distance from its source are
    summed into the global variable ``EscapedPhotonCount[]``. This variable
    also keeps track of the number of photons passing this radius
    multiplied by 0.5, 1, and 2. Units are in kpc. Not used if set to
    0. Default: 0.
``RadiativeTransferInterpolateField`` (obsolete)
    A failed experiment in which we evaluate the density at the
    midpoint of the ray segment in each cell to calculate the optical
    depth. To interpolate, we need to calculate the vertex interpolated
    density fields. Default: 0.
``RadiativeTransferSourceClustering`` (internal)
    Set to 1 to turn on ray merging. Not fully tested and may still be
    buggy. Default: 0.
``RadiativeTransferPhotonMergeRadius`` (internal)
    The radius at which the rays will merge from their SuperSource,
    which is the luminosity weighted center of two sources. This radius
    is in units of the separation of two sources associated with one
    SuperSource. If set too small, there will be angular artifacts in
    the radiation field. Default: 10.
``RadiativeTransferTimestepVelocityLimit`` (external)
    Limits the radiative transfer timestep to a minimum value that is
    determined by the cell width at the finest level divided by this
    velocity. Units are in km/s. Default: 100.
``RadiativeTransferPeriodicBoundary`` (external)
    Set to 1 to turn on periodic boundary conditions for photon
    packages. Default: 0.
``RadiativeTransferTraceSpectrum`` (external)
    reserved for experimentation. Default: 0.
``RadiativeTransferTraceSpectrumTable`` (external)
    reserved for experimentation. Default: ``spectrum_table.dat``
``RadiationXRaySecondaryIon`` (external)
    Set to 1 to turn on secondary ionizations and reduce heating from
    X-ray radiation (Shull & van Steenberg 1985). Currently only BH and
    MBH particles emit X-rays. Default: 0.
``RadiationXRayComptonHeating`` (external)
    Set to 1 to turn on Compton heating on electrons from X-ray
    radiation (Ciotti & Ostriker 2001). Currently only BH and MBH
    particles emit X-rays. Default: 0.

Radiative Transfer (FLD) Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``RadiativeTransferFLD`` (external)
    Set to 2 to turn on the fld-based radiation solvers following Reynolds,
    Hayes, Paschos & Norman, 2009. Note that you also have to compile
    the source using ``make photon-yes`` and a ``make
    hypre-yes``. Note that if FLD is turned on, it will force
    ``RadiativeCooling = 0``, ``GadgetEquilibriumCooling = 0``, and
    ``RadiationFieldType = 0`` to prevent conflicts. Default: 0.
``ImplicitProblem`` (external)
    Set to 1 to turn on the implicit FLD solver, or 3 to turn on the
    split FLD solver. Default: 0.
``RadHydroParamfile`` (external)
    Names the (possibly-different) input parameter file containing
    solver options for the FLD-based solvers. These are described in
    the relevant User Guides, located in ``doc/implicit_fld`` and
    ``doc/split_fld``. Default: NULL.
``RadiativeTransfer`` (external)
    Set to 0 to avoid conflicts with the ray tracing solver above.
    Default: 0.
``RadiativeTransferFLDCallOnLevel`` (reserved)
    The level in the static AMR hierarchy where the unigrid FLD solver
    should be called. Currently only works for 0 (the root grid).
    Default: 0.
``RadiativeTransferOpticallyThinH2`` (external)
    Set to 0 to avoid conflicts with the built-in optically-thin H_2
    dissociating field from the ray-tracing solver. Default: 1.

Radiative Transfer (FLD) Implicit Solver Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Massive Black Hole Physics Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following parameters are for the accretion and feedback from the
massive black hole particle (``PARTICLE_TYPE_MBH``). More details
will soon be described in Kim et al. (2010).

Accretion Physics
^^^^^^^^^^^^^^^^^

``MBHAccretion`` (external)
    Set to 1 to turn on accretion based on the Eddington-limited
    spherical Bondi-Hoyle formula (Bondi 1952). Set to 2 to turn on
    accretion based on the Bondi-Hoyle formula but with fixed
    temperature defined below. Set to 3 to turn on accretion with a
    fixed rate defined below. Set to 4 to to turn on accretion based on
    the Eddington-limited spherical Bondi-Hoyle formula, but without
    v_rel in the denominator. Set to 5 to turn on accretion based on
    Krumholz et al.(2006) which takes vorticity into account. Set to 6 
    to turn on alpha disk formalism based on DeBuhr et al.(2010).  
    7 and 8 are still failed experiment. Add 10 to each of these options 
    (i.e. 11, 12, 13, 14) to ignore the Eddington limit. See
    ``Star_CalculateMassAccretion.C``. Default: 0 (FALSE)
``MBHAccretionRadius`` (external)
    This is the radius (in pc) of a gas sphere from which the accreting
    mass is subtracted out at every timestep. Instead, you may want to
    try set this parameter to -1, in which case an approximate Bondi
    radius is calculated and used (from ``DEFAULT_MU`` and
    ``MBHAccretionFixedTemperature``). If set to -N, it will use N\*(Bondi
    radius). See ``CalculateSubtractionParameters.C``. Default: 50.0
``MBHAccretingMassRatio`` (external)
    There are three different scenarios you can utilize this parameter.
    (1) In principle this parameter is a nondimensional factor
    multiplied to the Bondi-Hoyle accretion rate; so 1.0 should give
    the plain Bondi rate. (2) However, if the Bondi radius is resolved
    around the MBH, the local density used to calculate Mdot can be
    higher than what was supposed to be used (density at the Bondi
    radius!), resulting in the overestimation of Mdot. 0.0 <
    ``MBHAccretingMassRatio`` < 1.0 can be used to fix this. (3) Or, one
    might try using the density profile of R\ :sup:`-1.5`\  to estimate
    the density at the Bondi radius, which is utilized when
    ``MBHAccretingMassRatio`` is set to -1. See
    ``Star_CalculateMassAccretion.C``. Default: 1.0
``MBHAccretionFixedTemperature`` (external)
    This parameter (in K) is used when ``MBHAccretion = 2``. A fixed gas
    temperature that goes into the Bondi-Hoyle accretion rate
    estimation formula. Default: 3e5
``MBHAccretionFixedRate`` (external)
    This parameter (in Msun/yr) is used when ``MBHAccretion = 3``. Default:
    1e-3
``MBHTurnOffStarFormation`` (external)
    Set to 1 to turn off star formation (only for ``StarParicleCreation``
    method 7) in the cells where MBH particles reside. Default: 0
    (FALSE)
``MBHCombineRadius`` (external)
    The distance (in pc) between two MBH particles in which two
    energetically-bound MBH particles merge to form one particle.
    Default: 50.0
``MBHMinDynamicalTime`` (external)
    Minimum dynamical time (in yr) for a MBH particle. Default: 1e7
``MBHMinimumMass`` (external)
    Minimum mass (in Msun) for a MBH particle. Default: 1e3

Feedback Physics
^^^^^^^^^^^^^^^^

``MBHFeedback`` (external)
    Set to 1 to turn on thermal feedback of MBH particles (``MBH_THERMAL``
    - not fully tested). Set to 2 to turn on mechanical feedback of MBH
    particles (``MBH_JETS``, bipolar jets along the total angular momentum
    of gas accreted onto the MBH particle so far). Set to 3 to turn on
    another version of mechanical feedback of MBH particles (``MBH_JETS``, 
    always directed along z-axis). Set to 4 to turn on experimental version of 
    mechanical feedback (`MBH_JETS`, bipolar jets along the total angular 
    momentum of gas accreted onto the MBH particle so far + 10 degree random 
    noise).  Set to 5 to turn on experimental version of mechanical feedback
    (``MBH_JETS``, launched at random direction). Note that, even when this
    parameter is set to 0, MBH particles still can be radiation sources
    if ``RadiativeTransfer`` is on. See ``Grid_AddFeedbackSphere.C``.
    Default: 0 (FALSE)

   ::
 
     ``RadiativeTransfer = 0`` & ``MBHFeedback = 0`` : no feedback at all
     ``RadiativeTransfer = 0`` & ``MBHFeedback = 1`` : purely thermal feedback
     ``RadiativeTransfer = 0`` & ``MBHFeedback = 2`` : purely mechanical feedback
     ``RadiativeTransfer = 1`` & ``MBHFeedback = 0`` : purely radiative feedback
     ``RadiativeTransfer = 1`` & ``MBHFeedback = 2`` : radiative and
       mechanical feedback combined (one has to change the following
       ``MBHFeedbackRadiativeEfficiency`` parameter accordingly, say from 0.1
       to 0.05, to keep the same total energy across different modes of
       feedback)

``MBHFeedbackRadiativeEfficiency`` (external)
    The radiative efficiency of a black hole. 10% is the widely
    accepted value for the conversion rate from the rest-mass energy of
    the accreting material to the feedback energy, at the innermost
    stable orbit of a non-spinning Schwarzschild black hole (Shakura &
    Sunyaev 1973, Booth & Schaye 2009). Default: 0.1
``MBHFeedbackEnergyCoupling`` (external)
    The fraction of feedback energy that is thermodynamically (for
    ``MBH_THERMAL``) or mechanically (for ``MBH_JETS``) coupled to the gas.
    0.05 is widely used for thermal feedback (Springel et al. 2005, Di
    Matteo et al. 2005), whereas 0.0001 or less is recommended for
    mechanical feedback depending on the resolution of the simulation
    (Ciotti et al. 2009). Default: 0.05
``MBHFeedbackMassEjectionFraction`` (external)
    The fraction of accreting mass that is returning to the gas phase.
    For either ``MBH_THERMAL`` or ``MBH_JETS``. Default: 0.1
``MBHFeedbackMetalYield`` (external)
    The mass fraction of metal in the ejected mass. Default: 0.02
``MBHFeedbackThermalRadius`` (external)
    The radius (in pc) of a sphere in which the energy from
    ``MBH_THERMAL`` feedback is deposited. If set to a negative value, the
    radius of a sphere gets bigger in a way that the sphere encloses
    the constant mass (=
    4/3\*pi\*(-``MBHFeedbackThermalRadius``)\ :sup:`3`\  Msun). The latter
    is at the moment very experimental; see ``Star_FindFeedbackSphere.C``.
    Default: 50.0
``MBHFeedbackJetsThresholdMass`` (external)
    The bipolar jets by ``MBH_JETS`` feedback are injected every time the
    accumulated ejecta mass surpasses ``MBHFeedbackJetsThresholdMass`` (in
    Msun). Although continuously injecting jets into the gas cells
    might sound great, unless the gas cells around the MBH are resolved
    down to Mdot, the jets make little or no dynamical impact on the
    surrounding gas. By imposing ``MBHFeedbackJetsThresholdMass``, the jets
    from MBH particles are rendered intermittent, yet dynamically
    important. Default: 10.0
``MBHParticleIO`` (external)
    Set to 1 to print out basic information about MBH particles. Will
    be automatically turned on if ``MBHFeedback`` is set to 2 or 3.
    Default: 0 (FALSE)
``MBHParticleIOFilename`` (external)
    The name of the file used for the parameter above. Default:
    ``mbh_particle_io.dat``

Conduction
~~~~~~~~~~~~~~~~~~~~~~

Isotropic and anisotropic thermal conduction are implemented using the
method of Parrish and Stone: namely, using an explicit, forward
time-centered algorithm.  In the anisotropic conduction, heat can only
conduct along magnetic field lines.  One can turn on the two types of
conduction independently, since there are situations where one might 
want to use both.  The Spitzer fraction can be also set
independently for the isotropic and anisotropic conduction.

``IsotropicConduction`` (external)
    Turns on isotropic thermal conduction using Spitzer conduction.  Default: 0 (FALSE)
``AnisotropicConduction`` (external)
    Turns on anisotropic thermal conduction using Spitzer conduction.
    Can only be used if MHD is turned on (``HydroMethod`` = 4).
    Default: 0 (FALSE)
``IsotropicConductionSpitzerFraction`` (external)
    Prefactor that goes in front of the isotropic Spitzer conduction
    coefficient.  Should be a value between 0 and 1.
    Default: 1.0
``AnisotropicConductionSpitzerFraction`` (external)
    Prefactor that goes in front of the anisotropic Spitzer conduction
    coefficient.  Should be a value between 0 and 1.
    Default: 1.0
``ConductionCourantSafetyNumber`` (external)
    This is a prefactor that controls the stability of the conduction
    algorithm.  In its current explicit formulation, it must be set to
    a value of 0.5 or less.
    Default: 0.5

Shock Finding Parameters
~~~~~~~~~~~~~~~~~~~~~~~~
For details on shock finding in Enzo see :ref:`shock_finding`.

``ShockMethod`` (external)
    This parameter controls the use and type of shock finding. Default: 0
    
    ::

	0 - Off
	1 - Temperature Dimensionally Unsplit Jumps
	2 - Temperature Dimensionally Split Jumps
	1 - Velocity Dimensionally Unsplit Jumps
	2 - Velocity Dimensionally Split Jumps

``ShockTemperatureFloor`` (external)
    When calculating the mach number using temperature jumps, set the
    temperature floor in the calculation to this value.

``StorePreShockFields`` (external)
    Optionally store the Pre-shock Density and Temperature during data output.


.. _testproblem_param:

Test Problem Parameters
-----------------------

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

    ::

              Test  LeftDensity LeftVelocity LeftPressure RightDensity RightVelocity RightPressure
              1.1   1.0         0.0          1.0          0.125        0.0           0.1
              1.2   1.0         -2.0         0.4          1.0          2.0           0.4
              1.3   1.0         0.0          1000.0       1.0          0.0           0.01
              1.4   1.0         0.0          0.01         1.0          0.0           100.0
              1.5   5.99924     19.5975      460.894      5.99242      -6.19633      46.0950


``ShockTubeBoundary`` (external)
    Discontinuity position. Default: 0.5
``ShockTubeDirection`` (external)
    Discontinuity orientation. Type: integer. Default: 0 (shock(s) will
    propagate in x-direction)
``ShockTubeLeftDensity``, ``ShockTubeRightDensity`` (external)
    The initial gas density to the left and to the right of the
    discontinuity. Default: 1.0 and 0.125, respectively
``ShockTubeLeftVelocity``, ``ShockTubeRightVelocity`` (external)
    The same as above but for the velocity component in
    ``ShockTubeDirection``. Default: 0.0, 0.0
``ShockTubeLeftPressure``, ``ShockTubeRightPressure`` (external)
    The same as above but for pressure. Default: 1.0, 0.1

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
``CollapseTestSphereInitialLevel`` (external)
    Failed experiment to try to force refinement to a specified level.
    Not working. Default: 0.

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
    close the universe. Typical value 0.94. Default: 0.0 (no dark
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

Other External Parameters
-------------------------

``huge_number`` (external)
    The largest reasonable number. Rarely used. Default: 1e+20
``tiny_number`` (external)
    A number which is smaller than all physically reasonable numbers.
    Used to prevent divergences and divide-by-zero in C++ functions.
    Modify with caution! Default: 1e-20.

    An independent analog, ``tiny``, defined in ``fortran.def``, does the same
    job for a large family of FORTRAN routines. Modification of ``tiny`` must
    be done with caution and currently requires recompiling the code, since
    ``tiny`` is not a runtime parameter.

``TimeActionParameter[#]``
    Reserved for future use.
``TimeActionRedshift[#]``
    Reserved for future use.
``TimeActionTime[#]``
    Reserved for future use.
``TimeActionType[#]``
    Reserved for future use.

Other Internal Parameters
-------------------------

``TimeLastRestartDump``
    Reserved for future use.
``TimeLastDataDump`` (internal)
    The code time at which the last time-based output occurred.
``TimeLastHistoryDump``
    Reserved for future use.
``TimeLastMovieDump`` (internal)
    The code time at which the last movie dump occurred.
``CycleLastRestartDump``
    Reserved for future use.
``CycleLastDataDump`` (internal)
    The last cycle on which a cycle dump was made
``CycleLastHistoryDump``
    Reserved for future use.
``InitialCPUTime``
    Reserved for future use.
``InitialCycleNumber`` (internal)
    The current cycle
``RestartDumpNumber``
    Reserved for future use.
``DataLabel[#]`` (internal)
    These are printed out into the restart dump parameter file. One
    Label is produced per baryon field with the name of that baryon
    field. The same labels are used to name data sets in HDF files.
``DataUnits[#]``
    Reserved for future use.
``DataDumpNumber`` (internal)
    The identification number of the next output file (the 0000 part of
    the output name). This is used and incremented by both the cycle
    based and time based outputs. Default: 0
``HistoryDumpNumber``
    Reserved for future use.
``MovieDumpNumber`` (internal)
    The identification number of the next movie output file. Default: 0
``VersionNumber`` (internal)
    Sets the version number of the code which is written out to restart
    dumps.

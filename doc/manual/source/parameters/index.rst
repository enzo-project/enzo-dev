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

This list includes parameters for the Enzo 2.3 release.

.. toctree::
   :maxdepth: 2

* `Initialization Parameters`_

* `I/O Parameters`_ 
   
   * `General I/O Parameters`_ 
   
   * `Stopping Parameters`_
   
   * `Streaming Data Format`_
   
   * `Simulation Identifiers and UUIDs`_

* `Hierarchy Control Parameters`_

* `Gravity Parameters`_
   
   * `General Gravity Parameters`_
   
   * `External Gravity Source`_

* `Hydrodynamics Parameters`_
   
   * `General Hydrodynamics Parameters`_
   
   * `Minimum Pressure Support Parameters`_
   
   * `Magnetohydrodynamics (CT) Parameters`_
   
   * `Magnetohydrodynamics (Dedner) Parameters`_

* `Cooling Parameters`_
   
   * `Simple Cooling Options`_
   
   * `Cloudy Cooling`_
   
   * `The Grackle`_


* `Particle Parameters`_

* `Star Formation and Feedback Parameters`_
   
   * `General Star Formation`_
   
   * `Normal Star Formation`_
   
   * `Molecular Hydrogen Regulated Star Formation`_
   
   * `Population III Star Formation`_
   
   * `Radiative Star Cluster Formation`_
   
   * `Massive Black Hole Particle Formation`_
   
   * `Sink Formation and Feedback`_

* `Radiation Parameters`_
   
   * `Background Radiation Parameters`_
   
   * `Radiative Transfer (Ray Tracing) Parameters`_
   
   * `Radiative Transfer (FLD) Parameters`_
   
   * `Radiative Transfer (FLD) Implicit Solver Parameters`_
   
   * `Radiative Transfer (FLD) Split Solver Parameters`_

* `Cosmology Parameters`_

* `Massive Black Hole Physics Parameters`_
   
   * `Accretion Physics`_
   
   * `Feedback Physics`_

* `Shock Finding Parameters`_

* `Cosmic Ray Two-Fluid Model Parameters`_

* `Conduction`_

* `Inline Analysis`_
   
   * `Inline Halo Finding`_
   
   * `Inline Python`_

* `Other Parameters`_
   
   * `Other External Parameters`_
   
   * `Other Internal Parameters`_

* `Problem Type Parameters`_

    * `Shock Tube (1: unigrid and AMR)`_

    * `Wave Pool (2)`_

    * `Shock Pool (3: unigrid 2D, AMR 2D and unigrid 3D)`_

    * `Double Mach Reflection (4)`_

    * `Shock in a Box (5)`_

    * `Implosion (6)`_

    * `Sedov Blast (7)`_

    * `Kelvin-Helmholtz Instability (8)`_

    * `2D/3D Noh Problem (9)`_

    * `Rotating Cylinder (10)`_

    * `Radiating Shock (11)`_

    * `Free Expansion (12)`_

    * `Rotating Sphere (14)`_

    * `Radiating Shock (11)`_
    
    * `Free Expansion (12)`_

    * `Rotating Sphere (14)`_

    * `Zeldovich Pancake (20)`_

    * `Pressureless Collapse (21)`_

    * `Adiabatic Expansion (22)`_

    * `Test Gravity (23)`_

    * `Spherical Infall (24)`_

    * `Test Gravity: Sphere (25)`_

    * `Gravity Equilibrium Test (26)`_

    * `Collapse Test (27)`_

    * `Test Gravity Motion (28)`_

    * `Test Orbit (29)`_

    * `Cosmology Simulation (30)`_

    * `Isolated Galaxy Evolution (31)`_

    * `Shearing Box Simulation (35)`_

    * `Supernova Restart Simulation (40)`_

    * `Photon Test (50)`_

    * `Turbulence Simulation with Stochastic Forcing (59)`_

    * `Turbulence Simulation (60)`_

    * `Protostellar Collapse (61)`_

    * `Cooling Test (62)`_

    * `3D Collapse Test (101)`_

    * `1D Spherical Collapse Test (102)`_

    * `Hydro and MHD Turbulence Simulation (106)`_

    * `Put Sink from Restart (107)`_

    * `Cluster Cooling Flow (108)`_

    * `1D MHD Test (200)`_

    * `2D MHD Test (201)`_

    * `3D MHD Collapse Test (202)`_

    * `MHD Turbulent Collapse Test (203)`_

    * `Galaxy Disk (207)`_

    * `AGN Disk (207)`_

    * `CR Shock Tube (250: unigrid and AMR)`_

    * `Poisson Solver Test (300)`_

    * `Radiation-Hydrodynamics Test 1 - Constant Fields (400)`_

    * `Radiation-Hydrodynamics Test 2 - Streams (401)`_

    * `Radiation-Hydrodynamics Test 3 - Pulse (402)`_

    * `Radiation-Hydrodynamics Test 4 - Grey Marshak Test (403)`_

    * `Radiation-Hydrodynamics Test 5 - Radiating Shock (404/405)`_

    * `Radiation-Hydrodynamics Tests 10 and 11 - I-Front Tests (410/411)`_

    * `Radiation-Hydrodynamics Test 12 - HI ionization of a clump (412)`_

    * `Radiation-Hydrodynamics Test 13 - HI ionization of a steep region (413)`_

    * `Radiation-Hydrodynamics Tests 14/15 - Cosmological HI ionization (414/415)`_

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
    Does not work with isolated gravity.  Option 0 is the slowest of
    the methods.  Option 1 is an aggressive version that is
    memory-intensive.  Option 2 tries to conserve memory at the
    expense of performance.  See also ``Unigrid`` above.  Default: 2.
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

.. _io_parameters:

I/O Parameters
~~~~~~~~~~~~~~

.. _general_io_parameters:

General I/O Parameters 
^^^^^^^^^^^^^^^^^^^^^^

There are three ways to specify the frequency of outputs:
time-based, cycle-based (a cycle is a top-grid timestep), and, for
cosmology simulations, redshift-based. There is also a shortened
output format intended for visualization (movie format). Please
have a look at :ref:`controlling_data_output` for more information.

``dtDataDump`` (external)
    The time interval, in code units, between time-based outputs. A
    value of 0 turns off the time-based outputs. Default: 0
``dtInterpolatedDataDump`` (external)
    The time interval, in code units, between time-based interpolated outputs. A
    value of 0 turns off the time-based outputs. Default: 0
``CycleSkipDataDump`` (external)
    The number of cycles (top grid timesteps) between cycle-based
    outputs. Zero turns off the cycle-based outputs. Default: 0
``SubcycleSkipDataDump`` (external)
    The number of subcycles between subcycle-based
    outputs. Zero turns off the subcycle-based outputs. Default: 0
``dtTracerParticleDump`` (external)
    The time interval, in code units, between time-based tracer particle outputs (defined in ComputeRandomForcingNormalization.C). A
    value of 0 turns off this output. Default: 0
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
``TracerParticleDumpName`` (external)
    The base file name used for tracer particle outputs.
    Default: 
``TracerParticleDumpDir`` (external)
    The dir name used for tracer particle outputs.
    Default: 
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
``FileDirectedOutput``
    If this parameter is set to 1, whenever the finest level has finished
    evolving Enzo will check for new signal files to output.  (See
    :ref:`force_output_now`.)  Default 1.
``TracerParticleOn``
    This parameter is used to set the velocities of the tracer
    particles equal to the gas velocities in the current cells.
    Tracer particles are massless and can be used to output values of
    the gas as they advect with the fluid.  Default: 0
``TracerParticleOutputVelocity``
    This parameter is used to output tracer particle velocity as well
    as position, density, and temperature.  Default: 0
``OutputFirstTimeAtLevel`` (external)
    This forces Enzo to output when a given level is reached, and at
    every level thereafter. Default is 0 (off). User can usefully
    specify anything up to the maximum number of levels in a given
    simulation.
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
    for runs involving > 64 cpus! Default: 0 (FALSE). 
    See ``ParallelParticleIO`` in :ref:`particle_parameters`.    
    See also ``Unigrid`` in :ref:`initialization_parameters`.
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
``SmoothedDarkMatterNeighbors`` (external)
    Number of nearest neighbors to smooth dark matter quantities over.
    Default: 32.
``OutputGriddedStarParticle`` (external)
    Set to 1 or 2 to write out star particle data gridded onto mesh.
    This will be useful e.g. if you have lots of star particles in a
    galactic scale simulation. 1 will output just
    ``star_particle_density``; and 2 will dump
    ``actively_forming_stellar_mass_density``, ``SFR_density``, etc.
    Default: 0.
``PopIIIOutputOnFeedback`` (external)
    Writes an interpolated output when a Pop III is formed or goes
    supernova.  Default: 0
``OutputOnDensity`` (external)
    Should interpolated outputs be generated at varying peak density?
    Default: 0
``StartDensityOutput`` (external)
    The first density (in log g/cc) at which to output.
``CurrentDensityOutput`` (internal)
    The most recent density at which output was generated.
``IncrementDensityOutput`` (external)
    After a density-directed output, how much should the density be increased by?  Default: 999
``ComputePotential`` (external)
    When turned on, the gravitational potential is computed and stored in memory.  Always done when SelfGravity is on.  Default: 0
``WritePotential`` (external)
    When turned on, the gravitational potential is written to file.  Default: 0
``WriteGhostZones`` (external)
    Should ghost zones be written to disk?  Default: 0 
``ReadGhostZones`` (external)
    Are ghost zones present in the files on disk?  Default: 0
``VelAnyl`` (external)
    Set to 1 if you want to output the divergence and vorticity of
    velocity. Works in 2D and 3D.
``BAnyl`` (external)
    Set to 1 if you want to output the divergence and vorticity of
    ``Bfield``. Works in 2D and 3D.
``ExtractFieldsOnly`` (external)
    Used for extractions (enzo -x ...) when only field data are needed
    instead of field + particle data. Default is 1 (TRUE).
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
``ParticleTypeInFile`` (external)    
    Output ParticleType to disk?  Default: 1
``OutputParticleTypeGrouping`` (external)   
    In the grid HDF5 groups, particles are sorted by type, and a reference is created to indicate which particle index range corresponds to each type.  Default: 0
``HierarchyFileInputFormat`` (external) 
    See :ref:`controlling_the_hierarhcy_file_output`.
``HierarchyFileOutputFormat`` (external) 
    See :ref:`controlling_the_hierarhcy_file_output`.
``TimingCycleSkip`` (external)
    Controls how many cycles to skip when timing information is collected, reduced, and written out to performance.out.  Default: 1
``DatabaseLocation`` (external)
    (Not recommended for use at this point)  Where should the SQLite database of outputs be placed?
``CubeDumpEnabled`` (external)
    not recommended for use at this point. Default: 0
``CubeDump[]`` (external)
    not recommended for use at this point
``LocalDir`` (external) 
    See :ref:`controlling_data_output`.
``GlobalDir`` (external) 
    See :ref:`controlling_data_output`.

.. _stopping_parameters:

Stopping Parameters
^^^^^^^^^^^^^^^^^^^

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
``StopFirstTimeAtDensity`` (external)
    Causes the simulation to immediately stop when the maximum gas
    density reaches this value.  In units of proper g/cm^3.  Not used if less
    than or equal to zero. Default: 0.0
``StopFirstTimeAtMetalEnrichedDensity`` (external)
    Causes the simulation to immediately stop when the maximum gas
    density with above some metallicity, specified by
    ``EnrichedMetalFraction``, is reached.  In units of g/cm^3.  Not
    used if less than or equal to zero.  Default: 0.0
``EnrichedMetalFraction`` (external)
    See ``StopFirstTimeAtMetalEnrichedDensity``.  In units of absolute
    metal fraction.  Default: 1e-8
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

.. _streaming_param:

Streaming Data Format
^^^^^^^^^^^^^^^^^^^^^

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

.. _simulation_identifiers_parameters:

Simulation Identifiers and UUIDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    19                 Refine by metal mass
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
``SlopeFlaggingFields`` (external)
    If ``CellFlaggingMethod`` is 1, and you only want to refine on the
    slopes of certain fields then you can enter the
    :ref:`Field Type IDs <Field_List_Reference>` of the fields you want,
    separating the IDs with a space. Up to 7 Field Type IDs can be 
    specified. Default: Refine on slopes of all fields.
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
    For method 13 (metallicity refinement), this is the threshold
    metallicity (in units of solar metallicity) above which cells must
    be refined to a minimum level of
    ``MetallicityRefinementMinLevel``.  For method 19 (metal mass),
    this flags cells for refinement when the metal mass is above the
    necessary baryon mass (method 2) for refinement multiplied by this
    parameter.  Behaves similarly to refinement by baryon mass but
    focuses on metal-enriched regions.  In units of solar metallicity.
    Default: 1.0e-5
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
    is being used, and this parameter is greater than zero, this
    temperature will be used in all cells when calculating the Jeans
    length.  If it is less than or equal to zero, it will be used as a
    temperature floor when calculating the Jeans length. Default: -1.0
``RefineByResistiveLengthSafetyFactor`` (external)
    Resistive length is defined as the curl of the magnetic field over
    the magnitude of the magnetic field. We make sure this length is
    covered by this number of cells. i.w. The resistive length in a MHD simulation should not be smaller than CellWidth * RefineByResistiveLengthSafetyFactor.  Default: 2.0
``MustRefineParticlesCreateParticles`` (external)
    This parameter will flag dark matter particles in cosmological 
    initial conditions as ``MustRefineParticles``.  If ``CellFlaggingMethod`` 
    8 is set, AMR will be restricted to cells surrounding 
    ``MustRefineParticles``.  There are several different modes for creating
    ``MustRefineParticles`` with this parameter described below.  Further 
    information on how to use dark matter ``MustRefineParticles`` in 
    cosmological simulations can be found here (link).  Default: 0

    1. If the user specifies ``MustRefineParticlesLeftEdge`` and 
       ``MustRefineParticlesRightEdge``, dark matter particles within the 
       specified region are flagged.  Otherwise, the code looks for an ascii 
       input file called MustRefineParticlesFlaggingList.in that contains a list
       of particle ids to be flagged.  The ids in this list must be sorted in 
       acending order.
    2. For use with ellipsodial masking in MUSIC inititial conditions.  This
       setting uses traditional static grids for intermediate resolution levels
       MUSIC will generate RefinementMask files and the ``ParticleTypeName``
       parameter should be set to the name of these files.
    3. Same as setting 2, except refinement on intermediate levels is not
       constrained by static grids.  Instead, refinement around dark matter
       particles is allowed down to the level of a particle's generation level.
       Refinement beyond this level is allowed around particles within the MUSIC
       ellipsoidal making region.

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
    This was an experimental parameter to set a minimum for ``MustRefineParticles``. Default: 0.0
``MustRefineParticlesRegionLeftEdge`` (external)
    Bottom-left corner of a region in which dark matter particles are flagged 
    as ``MustRefineParticles`` in nested cosmological simulations.  To be used with 
    ``MustRefineParticlesCreateParticles`` = 1.  Default: 0.0 0.0 0.0
``MustRefineParticlesRegionRightEdge`` (external)
    Top-right corner of a region in which dark matter particles are flagged 
    as ``MustRefineParticles`` in nested cosmological simulations.  To be used with 
    ``MustRefineParticlesCreateParticles`` = 1.  Default: 0.0 0.0 0.0
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
    Default: IND_UNDEFINED
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
``CriticalGridRatio`` (external)
    Critical grid ratio above which subgrids will be split in half along their 
    long axis prior to being split by the second derivative of their 
    signature.  Default: 3.0
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

.. _gravity_parameters:

Gravity Parameters
~~~~~~~~~~~~~~~~~~

.. _general_gravity_parameters:

General Gravity Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

``TopGridGravityBoundary`` (external)
    A single integer which specified the type of gravitational boundary
    conditions for the top grid. Possible values are 0 for periodic and
    1 for isolated (for all dimensions). The isolated boundary
    conditions have not been tested recently, so caveat emptor.
    Default: 0
``SelfGravity`` (external)
    This flag (1 - on, 0 - off) indicates if the baryons and particles
    undergo self-gravity.
``SelfGravityGasOff`` (external)
    This parameter is used in conjuction with SelfGravity so that only particles contribute to potential, not gas. Default = False (i.e. gas does contribute)
``GravitationalConstant`` (external)
    This is the gravitational constant to be used in code units. For cgs units it
    should be 4\*pi\*G. For cosmology, this value must be 1 for the
    standard units to hold. A more detailed decription can be found at :ref:`EnzoInternalUnits`. Default: 4\*pi.
``PotentialIterations`` (external)
    Number of iterations to solve the potential on the subgrids. Values
    less than 4 sometimes will result in slight overdensities on grid
    boundaries. Default: 4.
``MaximumGravityRefinementLevel`` (external)
    This is the lowest (most refined) depth that a gravitational
    acceleration field is computed. More refined levels interpolate
    from this level, providing a mechanism for instituting a minimum
    gravitational smoothing length. Default: ``MaximumRefinementLevel``
``MaximumParticleRefinementLevel`` (external)
    This is the level at which the dark matter particle contribution to
    the gravity is smoothed. This works in an inefficient way (it
    actually smoothes the particle density onto the grid), and so is
    only intended for highly refined regions which are nearly
    completely baryon dominated. It is used to remove the discreteness
    effects of the few remaining dark matter particles. Not used if set
    to a value less than 0. Default: -1
``ParticleSubgridDepositMode`` (external)
    This parameter controls how particles stored in subgrid are deposited
    into the current grid.  Options are:

     0. (CIC_DEPOSIT) - This is a second-order, cloud-in-cell deposition
         method in which the cloud size is equal to the cell size in
         the target grid (particles are in source grid, deposited into
         target grid).  This method preserves the correct center-of-mass
         for a single particle but smears out boundaries and can result
         in small artifacts for smooth particle distributions (e.g.
         nested cosmological simulations with low perturbations).
     1. (CIC_DEPOSIT_SMALL) - This is also a CIC method, but the cloud
         size is taken to be the cell size in the source grid, so for
         subgrids, the cloud is smaller than the grid size.  This
         is an attempt to compromise between the other two methods.
     2. (NGP_DEPOSIT) - This uses a first order, nearest-grid-point
        method to deposit particle mass.  It does not preserve center-
        of mass position and so for single particle results in noisy
        accelerations.  However, it does correctly treat nested
        cosmology simulations with low initial perturbations.

    Default: 1
``BaryonSelfGravityApproximation`` (external)
    This flag indicates if baryon density is derived in a strange,
    expensive but self-consistent way (0 - off), or by a completely
    reasonable and much faster approximation (1 - on). This is an
    experiment gone wrong; leave on. Well, actually, it's important for
    very dense structures as when radiative cooling is turned on, so
    set to 0 if using many levels and radiative cooling is on [ignored
    in current version]. Default: 1

.. _external_gravity_source:

External Gravity Source
^^^^^^^^^^^^^^^^^^^^^^^

These parameters set up an external static background gravity source that is
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
    ``PointSourceGravityCoreRadius``. ``ProblemType`` = 108 (elliptical galaxy and galaxy cluster) also includes the gravity from the stellar component and the SMBH. Default: 1
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
   more aptly named. ``ExternalGravity = 1`` turns on an alternative
   implementation of the NFW profile with properties
   defined via the parameters ``HaloCentralDensity``, ``HaloConcentration`` and ``HaloVirialRadius``. Boxsize is assumed to be 1.0 in this case. ``ExternalGravity = 10`` gives a gravitational field defined by the logarithmic potential in Binney & Tremaine, corresponding to a disk with constant circular velocity.  Default: 0 
``ExternalGravityConstant`` (external)
    If ``ExternalGravity = 10``, this is the circular velocity of the disk in code units. Default: 0.0
``ExternalGravityDensity`` 
   Reserved for future use.
``ExternalGravityPosition`` (external)
    If ``ExternalGravity = 10``, this parameter specifies the center of the gravitational field in code units. Default: 0 0 0
``ExternalGravityOrientation`` (external)
    For ``ExternalGravity = 10``, this is the unit vector of the disk's angular momentum (e.g. a disk whose face-on view is oriented in the x-y plane would have ``ExternalGravityOrientation = 0 0 1``). Default: 0 0 0 
``ExternalGravityRadius`` (external)
   If ``ExternalGravity = 10``, this marks the inner radius of the disk in code units within which the velocity drops to zero. Default: 0.0
``UniformGravity`` (external)
    This flag (1 - on, 0 - off) indicates if there is to be a uniform
    gravitational field. Default: 0
``UniformGravityDirection`` (external)
    This integer is the direction of the uniform gravitational field: 0
    - along the x axis, 1 - y axis, 2 - z axis. Default: 0
``UniformGravityConstant`` (external)
    Magnitude (and sign) of the uniform gravitational acceleration.
    Default: 1
``DiskGravity`` (external)
    This flag (1 - on, 0 - off) indicates if there is to be a
    disk-like gravity field (Berkert 1995; Mori & Burkert 2000).  Default: 0
``DiskGravityPosition`` (external)
    This indicates the position of the center of the disk gravity.
    Default: 0 0 0
``DiskGravityAngularMomentum`` (external)
    Specifies the unit vector of the disk angular momentum.
    Default: 0 0 1
``DiskGravityStellarDiskMass`` (external)
    Total mass of stellar disk (in solar masses)
    Default: 1e11
``DiskGravityDiskScaleHeightR`` (external)
    Disk scale length in radius (in Mpc)
    Default: 4.0e-3
``DiskGravityDiskScaleHeightz`` (external)
    Disk scale height in z (in Mpc)
    Default: 2.5e-4
``DiskGravityStellarBulgeMass`` (external)
    Disk stellar bulge mass (in solar masses)
    Default: 1.0e10
``DiskGravityStellarBulgeR`` (external)
    Disk stellar bulge scalue radius (in Mpc)
    Default: 1.0e-4
``DiskGravityDarkMatterR`` (external)
    Dark matter halo scale radius (in Mpc)
    Default: 2.3e-2
``DiskGravityDarkMatterDensity`` (external)
    Dark matter effective density (in cgs)
    Default: 3.81323e-25

.. _hydrodynamics_parameters:

Hydrodynamics Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

.. _general_hydrodynamics_parameters:

General Hydrodynamics Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

.. _minimum_pressure_support_parameters:

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

.. _mhd_ct_parameters:

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

.. _mhd_dender_parameters:

Magnetohydrodynamics (Dedner) Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following parameters are considered only when ``HydroMethod`` is 3 or 4 (and occasionally only in some test problems).  
Because many of the following parameters are not actively being tested and maintained, users are encouraged to carefully examine the code before using it.

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


.. _cooling_parameters:

Cooling Parameters
~~~~~~~~~~~~~~~~~~

.. _simple_cooling_parameters:

Simple Cooling Options
^^^^^^^^^^^^^^^^^^^^^^

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
``RadiativeCoolingModel`` (external)
    This switches between the tabular look up cooling that is standard (RadiativeCoolingModel=1) and an analytic fit to the Wolfire et al 2003, ApJ, 587, 278 made by Koyama and Inutsuka 2006 (RadiativeCoolingModel = 3, arXiv:astro-ph/0605528).  Default: 1
``GadgetCooling`` (external)
    This flag (1 - on, 0 - off) turns on (when set to 1) a set of
    routines that calculate cooling rates based on the assumption of a
    six-species primordial gas (H, He, no H2 or D) in equilibrium, and
    is valid for temperatures greater than 10,000 K. This requires the
    file ``TREECOOL`` to execute. Default: 0
``GadgetEquilibriumCooling`` (external)
    An implementation of the ionization equilibrium cooling code used
    in the GADGET code which includes both radiative cooling and a
    uniform metagalactic UV background specified by the ``TREECOOL`` file
    (in the ``amr_mpi/exe`` directory). When this parameter is turned on,
    ``MultiSpecies`` and ``RadiationFieldType`` are forced to 0 and
    ``RadiativeCooling`` is forced to 1.
    [Not in public release version]
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
``MultiMetals`` (external)
    This was added so that the user could turn on or off additional
    metal fields - currently there is the standard metallicity field
    (Metal_Density) and two additional metal fields (Z_Field1 and
    Z_Field2). Acceptable values are 1 or 0, Default: 0 (off).
``ThreeBodyRate`` (external)
    Which Three Body rate should be used for H2 formation?: 0 = Abel, Bryan, Norman 2002, 1 = PSS83, 2= CW83, 3 = FH07, 4= G08.  (See `Turk et al 2011 <http://adsabs.harvard.edu/abs/2011ApJ...726...55T>`)
``CIECooling`` (external)
    Should CIE (`Ripamonti & Abel 2004 <http://adsabs.harvard.edu/abs/2004MNRAS.348.1019R>`) cooling be included at high densities?
``H2OpticalDepthApproximation`` (external)
    Should the H2 cooling be attenuated? Taken from `Ripamonti & Abel 2004 <http://adsabs.harvard.edu/abs/2004MNRAS.348.1019R>`. Default: 1?
``H2FormationOnDust`` (external)
    Turns on H2 formation on dust grains and gas-grain heat transfer following `Omukai (2000) <http://adsabs.harvard.edu/abs/2000ApJ...534..809O>`. Default: 0 (OFF)
``NumberOfDustTemperatureBins`` (external)
    Number of dust temperature bins for the dust cooling and H2 formation rates.  Default: 250
``DustTemperatureStart`` (external)
    Minimum dust temperature for dust rates.  Default: 1.0
``DustTemperatureEnd`` (external)
    Maximum dust temperature for dust rates.  Default: 1500
``OutputDustTemperature`` (external)
    Flag to write out the dust temperature field.  Default: 0
``PhotoelectricHeating`` (external)
    If set to be 1, the following parameter will be added uniformly
    to the gas without any shielding (`Tasker & Bryan 2008 <http://adsabs.harvard.edu/abs/2008ApJ...673..810T>`). Default: 0
``PhotoelectricHeatingRate`` (external)
    This is the parameter used as Gamma_pe for uniform photoelectric heating.
    Default: 8.5e-26 erg s^-1 cm^-3

.. _cloudy_cooling:

Cloudy Cooling
^^^^^^^^^^^^^^

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
``CMBTemperatureFloor`` (external)
    An integer (0 or 1) specifying whether a temperature floor is
    created at the temperature of the cosmic microwave background
    (T\ :sub:`CMB`\  = 2.72 (1 + z) K). This is accomplished in the
    code by subtracting the cooling rate at T\ :sub:`CMB`\  such that
    Cooling = Cooling(T) - Cooling(T\ :sub:`CMB`\ ). Default: 1 (on).
``CloudyElectronFractionFactor`` (external)
    A float value to account for additional electrons contributed by
    metals. This is only used with Cloudy datasets with dimension
    greater than or equal to 4. The value of this factor is calculated
    as the sum of (A\ :sub:`i`\  \* i) over all elements i heavier than
    He, where A\ :sub:`i`\  is the solar number abundance relative to
    H. For the solar abundance pattern from the latest version of
    Cloudy, using all metals through Zn, this value is 9.153959e-3.
    Default: 9.153959e-3.

.. _grackle_pars:

The Grackle
^^^^^^^^^^^

The Grackle is an external chemistry and cooling library originally derived from 
Enzo's MultiSpecies chemistry and Cloudy cooling modules.  See :ref:`here <Grackle>` 
for a full description, including why you might use this over Enzo's internal 
chemistry and cooling.  For more information on Grackle parameter, see also the 
`Grackle documentation <https://grackle.readthedocs.org/>`_.  Note, some Grackle 
parameters have been mapped to Enzo parameters for simplicity.

``use_grackle`` (int)
    Flag to use the Grackle machinery (1 - on, 0 - off). Default: 0.

``with_radiative_cooling`` (int)
    Flag to include radiative cooling and actually update the thermal energy during the chemistry solver.  If off, the chemistry species will still be updated.  The most common reason to set this to off is to iterate the chemistry network to an equilibrium state (1 - on, 0 - off).  Default: 1.

``MultiSpecies`` (int) [mapped to Grackle parameter ``primordial_chemistry``]
    Flag to control which primordial chemistry network is used.  Default: 0.

    - 0: no chemistry network.  Radiative cooling for primordial species is solved by interpolating from lookup tables calculated with Cloudy.
    - 1: 6-species atomic H and He.  Active species: H, H\ :sup:`+`, He, He\ :sup:`+`, He\ :sup:`++`, e\ :sup:`-`.
    - 2: 9-species network including atomic species above and species for molecular hydrogen formation.  This network includes formation from the H\ :sup:`-` and H\ :sub:`2`\ :sup:`+` channels, three-body formation (H+H+H and H+H+H\ :sub:`2`), H\ :sub:`2` rotational transitions, chemical heating, and collision-induced emission (optional).  Active species: above + H\ :sup:`-`, H\ :sub:`2`, H\ :sub:`2`\ :sup:`+`.
    - 3: 12-species network include all above plus HD rotation cooling.  Active species: above plus D, D\ :sup:`+`, HD.

``H2FormationOnDust`` (int) [mapped to Grackle parameter ``h2_on_dust``]
    See Enzo equivalent above.  Default: 0.

``MetalCooling`` (int) [mapped to Grackle parameter ``metal_cooling``]
    Flag to enable metal cooling using the Cloudy tables.  If enabled, the cooling table to be used must be specified with the ``grackle_data_file`` parameter (1 - on, 0 - off).  Default: 0.

``CMBTemperatureFloor`` (int) [mapped to Grackle parameter ``cmb_temperature_floor``]
    See Enzo equivalent above.  Default: 1.

``UVbackground`` (int)
    Flag to enable a UV background.  If enabled, the cooling table to be used must be specified with the ``grackle_data_file`` parameter (1 - on, 0 - off).  Default: 0.

``grackle_data_file`` (string)
    Path to the data file containing the metal cooling and UV background tables.  Default: "".

``Gamma`` (float)
    See Enzo equivalent above.  Default:  5/3.

``ThreeBodyRate`` (int) [mapped to Grackle parameter ``three_body_rate``]
    See Enzo equivalent above.  Default: 0.

``CIECooling`` (int) [mapped to Grackle parameter ``cie_cooling``]
    See Enzo equivalent above.  Default: 0.

``H2OpticalDepthApproximation`` (int) [mapped to Grackle parameter ``h2_optical_depth_approximation``]
    See Enzo equivalent above.  Default: 0.

``PhotoelectricHeating`` (int) [mapped to Grackle parameter ``photoelectric_heating``]
    See Enzo equivalent above.  Default: 0.

``PhotoelectricHeatingRate`` (float) [mapped to Grackle parameter ``photoelectric_heating_rate``]
    See Enzo equivalent above.  Default: 8.5e-26.

``Compton_xray_heating`` (int)
   Flag to enable Compton heating from an X-ray background following `Madau & Efstathiou (1999) <http://adsabs.harvard.edu/abs/1999ApJ...517L...9M>`_.  Default: 0.

``LWbackground_intensity`` (float)
    Intensity of a constant Lyman-Werner H\ :sub:`2` photo-dissociating radiation field in units of 10\ :sup:`-21` erg s\ :sup:`-1` cm\ :sup:`-2` Hz\ :sup:`-1` sr\ :sup:`-1`.  Default: 0.

``LWbackground_sawtooth_suppression`` (int)
    Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption (giving a sawtooth pattern), taken from `Haiman & Abel, & Rees (2000) <http://adsabs.harvard.edu/abs/2000ApJ...534...11H>`_.  Default: 0.

.. _particle_parameters:

Particle Parameters
~~~~~~~~~~~~~~~~~~~

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
``ParallelParticleIO`` (external)
    Normally, for the mpi version, the particle data are read into the
    root processor and then distributed to separate processors.
    However, for very large number of particles, the root processor may
    not have enough memory. If this toggle switch is set on (i.e. to
    the value 1), then Ring i/o is turned on and each processor reads
    its own part of the particle data. More I/O is required, but it is
    more balanced in terms of memory. ``ParallelRootGridIO`` and
    ``ParallelParticleIO`` MUST be set for runs involving > 64 cpus!
    See also ``ParallelRootGridIO`` in :ref:`io_parameters`.
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

.. _starparticleparameters:

Star Formation and Feedback Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For details on each of the different star formation methods available in Enzo see :ref:`star_particles`.

.. _general_star_formation_parameters:

General Star Formation
^^^^^^^^^^^^^^^^^^^^^^

``StarParticleCreation`` (external)
    This parameter is bitwise so that multiple types of star formation
    routines can be used in a single simulation. For example if methods
    1 and 3 are desired, the user would specify 10 (2\ :sup:`1`\  +
    2\ :sup:`3`\ ), or if methods 1, 4 and 7 are wanted, this would be
    146 (2\ :sup:`1`\  + 2\ :sup:`4`\  + 2\ :sup:`7`\ ). Default: 0
    
    ::

      0  - Cen & Ostriker (1992)
      1  - Cen & Ostriker (1992) with stocastic star formation
      2  - Global Schmidt Law / Kravstov et al. (2003)
      3  - Population III stars / Abel, Wise & Bryan (2007)
      4  - Sink particles: Pure sink particle or star particle with wind feedback depending on 
           choice for HydroMethod / Wang et al. (2009)
      5  - Radiative star clusters  / Wise & Cen (2009)
      6  - [reserved for future use]
      7  - Cen & Ostriker (1992) with no delay in formation
      8  - Springel & Hernquist (2003)
      9  - Massive Black Hole (MBH) particles insertion by hand / Kim et al. (2010)
      10 - Population III stellar tracers  
      11 - Molecular hydrogen regulated star formation
      13 - Distributed stellar feedback model (So et al. 2014)
      14 - Cen & Ostriker (1992) stochastic star formation with kinetic feedback 
             / Simpson et al. (2015)

``StarParticleFeedback`` (external)
    This parameter works the same way as ``StarParticleCreation`` but only
    is valid for ``StarParticleCreation`` method = 0, 1, 2, 7, 8 and 14 because methods 3, 5 and 9
    use the radiation transport module and ``Star_*.C`` routines to
    calculate the feedback, 4 has explicit feedback and 10 does not use feedback. Default: 0.

``StarFeedbackDistRadius`` (external)
    If this parameter is greater than zero, stellar feedback will be
    deposited into the host cell and neighboring cells within this
    radius.  This results in feedback being distributed to a cube with
    a side of ``StarFeedbackDistRadius+1``. It is in units of cell
    widths of the finest grid which hosts the star particle.  Only
    implemented for ``StarParticleCreation`` method = 0 or 1 with ``StarParticleFeedback`` method =  1. (If ``StarParticleFeedback`` = 0, stellar feedback is only deposited into the cell in which the star particle lives).  Default: 0.

``StarFeedbackDistCellStep`` (external)
    In essence, this parameter controls the shape of the volume where
    the feedback is applied, cropping the original cube.  This volume
    that are within ``StarFeedbackDistCellSteps`` cells from the host
    cell, counted in steps in Cartesian directions, are injected with
    stellar feedback.  Its maximum value is ``StarFeedbackDistRadius``
    * ``TopGridRank``.  Only implemented for ``StarParticleCreation`` method = 0
    or 1  with ``StarParticleFeedback`` method =  1.  See :ref:`distributed_feedback` for an illustration.
    Default: 0.

``StarMakerTypeIaSNe`` (external)
    This parameter turns on thermal and chemical feedback from Type Ia
    supernovae.  The mass loss and luminosity of the supernovae are
    determined from `fits of K. Nagamine
    <http://www.physics.unlv.edu/~kn/SNIa_2/>`_.  The ejecta are
    traced in a separate species field, ``MetalSNIa_Density``.  The
    metallicity of star particles that comes from this ejecta is
    stored in the particle attribute ``typeia_fraction``.  Can be used
    with ``StarParticleCreation`` method = 0, 1, 2, 5, 7, 8, and 13.  Default:
    0.

``StarMakerPlanetaryNebulae`` (external) 
    This parameter turns on thermal and chemical feedback from
    planetary nebulae.  The mass loss and luminosity are taken from
    the same `fits from K. Nagamine
    <http://www.physics.unlv.edu/~kn/SNIa_2/>`_.  The chemical
    feedback injects gas with the same metallicity as the star
    particle, and the thermal feedback equates to a 10 km/s wind.  The
    ejecta are not stored in its own species field.  Can be used
    with ``StarParticleCreation`` method = 0, 1, 2, 5, 7, 8, and 13.  Default: 0.

``StarParticleRadiativeFeedback`` (external)
    By setting this parameter to 1, star particles created with
    methods (0, 1, 2, 5, 7, 8, 13) will become radiation sources with
    the UV luminosity being determined with the parameter
    ``StarEnergyToStellarUV``.  Default: OFF

.. _normal_star_formation_parameters:

Normal Star Formation
^^^^^^^^^^^^^^^^^^^^^

The parameters below are considered in ``StarParticleCreation`` method
0, 1, 2, 7, 8, 13 and 14.

``StarMakerOverDensityThreshold`` (external)
    The overdensity threshold in code units (for cosmological simulations, note that code units are relative to the total mean density, not
    just the dark matter mean density) before star formation will be
    considered. For ``StarParticleCreation`` method = 7 in cosmological
    simulations, however, ``StarMakerOverDensityThreshold`` should be in
    particles/cc, so it is not the ratio with respect to the
    ``DensityUnits`` (unlike most other
    star_makers). This way one correctly represents the Jeans
    collapse and molecular cloud scale physics even in cosmological
    simulations. Default: 100
``StarMakerSHDensityThreshold`` (external)
    The critical density of gas used in Springel & Hernquist star
    formation ( \\rho_{th} in the paper) used to determine the star
    formation timescale in units of g cm\ :sup:`-3`\ . Only valid for ``StarParticleCreation`` method = 8. Default: 7e-26.
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
``StarMakerTimeIndependentFormation`` (external)
    When used, the factor of dt / t_dyn is removed from the calculation of 
    the star particle mass above.  Instead of the local dynamical time, the 
    timescale over which feedback occurs is a constant set by the parameter 
    ``StarMakerMinimumDynamicalTime``.  This is necessary when running with 
    conduction as the timesteps can be very short, which causes the calculated 
    star particle mass to never exceed reasonable values for 
    ``StarMakerMinimumMass``.  This prevents cold, star-forming gas from 
    actually forming stars, and when combined with conduction, results in too 
    much heat being transferred out of hot gas.  When running a cosmological 
    simulation with conduction and star formation, one must use this otherwise 
    bad things will happen.  (1 - ON; 0 - OFF)  Default: 0.
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
    returned as UV radiation with a young star spectrum. This is used
    when calculating the radiation background. Default: 3e-6
``StarEnergyToQuasarUV`` (external)
    The fraction of the rest-mass energy of the stars created which is
    returned as UV radiation with a quasar spectrum. This is used when
    calculating the radiation background. Default: 5e-6
``StarFeedbackKineticFraction`` (external)
    Only valid for ``StarParticleFeedback`` method = 14.  If set to a zero or positive
    value between 0.0 and 1.0, this is the constant fraction of energy injected in kinetic 
    form.  If set to -1, then a variable kinetic fraction is used that depends on local
    gas density, metallicity and resolution.  See Simpson et al. 2015
    for details. Note, some failures may occur in -1 mode.  Default 0.0
``StarMakerExplosionDelayTime`` (external)
    Only valid for ``StarParticleFeedback`` method = 14.  If set to a positive value, energy,
    metals and mass from the particle are injected in a single timestep that is delayed from
    the particle creation time by this amount.  This value is in units of Myrs.  If set
    to a negative value, energy, mass and metals are injected gradually in the same way as is
    done for ``StarParticleFeedback`` method = 1.  Default -1.

.. _molecular_hydrogen_regulated_star_formation_parameters:

Molecular Hydrogen Regulated Star Formation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters below are considered in ``StarParticleCreation`` method 11.

``H2StarMakerEfficiency`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerNumberDensityThreshold`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerMinimumMass`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerMinimumH2FractionForStarFormation`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerStochastic`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerUseSobolevColumn`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerSigmaOverR`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerAssumeColdWarmPressureBalance`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerH2DissociationFlux_MW`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerH2FloorInColdGas`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``H2StarMakerColdGasTemperature`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.
``StarFormationOncePerRootGridTimeStep`` (external)
    See :ref:`molecular_hydrogen_regulated_star_formation`.

.. _popIII_star_formation_parameters:

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
``PopIIIBHLuminosityEfficiency`` (external)
    The radiative efficiency in which the black holes convert accretion
    to luminosity. Default: 0.1.
``PopIIIOverDensityThreshold`` (external)
    The overdensity threshold (relative to the total mean density)
    before Pop III star formation will be considered. Default: 1e6.
``PopIIIH2CriticalFraction`` (external)
    The H_2 fraction threshold before Pop III star formation will be
    considered. Default: 5e-4.
``PopIIIMetalCriticalFraction`` (external)
    The metallicity threshold (relative to gas density, not solar)
    before Pop III star formation will be considered. Note: this should
    be changed to be relative to solar! Default: 1e-4.
``PopIIISupernovaRadius`` (external)
    If the Population III star will go supernova (140<M<260 solar
    masses), this is the radius of the sphere to inject the supernova
    thermal energy at the end of the star's life. Units are in parsecs.
    Default: 1.
``PopIIISupernovaUseColour`` (external)
    Set to 1 to trace the metals expelled from supernovae. If using ``HydroMethod`` 3 or 4, also set ``MixSpeciesAndColors`` to 1 to trace metals. Default: 0.
``PopIIIUseHypernovae`` (external)
    Set to 1 to use the hypernova energies and metal ejecta masses
    from Nomoto et al. (2006).  If set to 0, then the supernova
    energies are always 1e51 erg but use the supernova metal ejecta
    masses from Nomoto et al. (2006).  Default: 1
``PopIIISupernovaExplosions`` (external)
    Set to 1 to consider supernovae from Pop III stars.  Set to 0 to
    neglect all Pop III supernovae, regardless of their masses.
    Default: 1
``PopIIIInitialMassFunction`` (external)
    When turned on, each Pop III stellar mass is randomly drawn from an IMF that is Salpeter above some characteristic mass and exponentially cutoff below this mass.  Default: 0
``PopIIIInitialMassFunctionSeed`` (external)
    Random initial seed for the Pop III stellar mass randomizer.  Default: INT_UNDEFINED
``PopIIILowerMassCutoff`` (external)
    Lower limit of the Pop III IMF.  Default: 1
``PopIIIUpperMassCutoff`` (external)
    Upper limit of the Pop III IMF.  Default: 300
``PopIIIInitialMassFunctionSlope`` (external)
    Slope of the Salpeter (high-mass) portion of the Pop III IMF.  Default: -1.3
``PopIIIInitialMassFunctionCalls`` (internal) 
    Number of times a Pop III mass has been drawn from the IMF.  Used for restarts and reproducibility.  Default: 0
``PopIIISupernovaMustRefine`` (external)
    When turned on, the region around a star about to go supernova is refined to the maximum AMR level.  Experimental.  Default: 0
``PopIIISupernovaMustRefineResolution`` (external)
    Used with PopIIISupernovaMustRefine.  Minimum number of cells across the blastwave.  Default: 32
``PopIIIHeliumIonization`` (external)
    When turned on, Pop III stars will emit helium singly- and doubly-ionizing radiation.  Default: 0
``PopIIIColorDensityThreshold`` (external)
    Above this density, a Pop III "color" particle forms, and it will populate the surrounding region with a color field.  Units: mean density. Default: 1e6
``PopIIIColorMass`` (external)
    A Pop III "color" particle will populate the surrounding region with a mass of PopIIIColorMass.  Units: solar masses.  Default: 1e6

.. _radiative_star_cluster_formation_parameters:

Radiative Star Cluster Formation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters below are considered in ``StarParticleCreation`` method 5.

``StarClusterMinDynamicalTime`` (external)
    When determining the size of a star forming region, one method is
    to look for the sphere with an enclosed average density that
    corresponds to some minimum dynamical time. Observations hint that
    this value should be a few million years. Units are in years.
    Default: 1e7.
``StarClusterIonizingLuminosity`` (external)
    The specific luminosity of the stellar clusters. In units of
    ionizing photons per solar mass. Default: 1e47.
``StarClusterSNEnergy`` (external)
    The specific energy injected into the gas from supernovae in the
    stellar clusters. In units of ergs per solar mass. Default: 6.8e48
    (Woosley & Weaver 1986).
``StarClusterSNRadius`` (external)
    This is the radius of the sphere to inject the supernova thermal
    energy in stellar clusters. Units are in parsecs. Default: 10.
``StarClusterFormEfficiency`` (external)
    Fraction of gas in the sphere to transfer from the grid to the star
    particle. Recall that this sphere has a minimum dynamical time set
    by ``StarClusterMinDynamicalTime``. Default: 0.1.
``StarClusterMinimumMass`` (external)
    The minimum mass of a star cluster particle before the formation is
    considered. Units in solar masses. Default: 1000.
``StarClusterCombineRadius`` (external)
    It is possible to merge star cluster particles together within this
    specified radius. Units in parsecs. This is probably not necessary
    if ray merging is used. Originally this was developed to reduce the
    amount of ray tracing involved from galaxies with hundreds of these
    radiating particles. Default: 10.
``StarClusterUseMetalField`` (external)
    Set to 1 to trace ejecta from supernovae. Default: 0.
``StarClusterHeliumIonization`` (external)
    When turned on, stellar clusters will emit helium singly- and doubly-ionizing radiation.  Default: 0
``StarClusterRegionLeftEdge`` (external)
    Can restrict the region in which star clusters can form.  Origin of this region.  Default: 0 0 0
``StarClusterRegionRightEdge`` (external)
    Can restrict the region in which star clusters can form.  Right corner of this region.  Default: 1 1 1
``StarClusterUnresolvedModel`` (external)
    Regular star clusters live for 20 Myr, but this is only valid when molecular clouds are resolved.  When this parameter is on, the star formation rate is the same as the Cen & Ostriker exponential rate.  Default: 0

.. _massive_black_hole_particle_parameters:

Massive Black Hole Particle Formation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters below are considered in ``StarParticleCreation`` method 9.

``MBHInsertLocationFilename`` (external)
    The mass and location of the MBH particle that has to be inserted.
    For example, the content of the file should be in the following
    form. For details, see ``mbh_maker.src``. Default:
    ``mbh_insert_location.in``
    ::

        #order: MBH mass (in Ms), MBH location[3], MBH creation time
        100000.0      0.48530579      0.51455688      0.51467896      0.0

.. _sink_formation_and_feedback_parameters:

Sink Formation and Feedback
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters below are considered in sink creation routines: sink_maker, star_maker8, star_maker9 (and occasionally only in certain set-ups).  
Because many of the following parameters are not actively being tested and maintained, users are encouraged to carefully examine the code before using it.

``AccretionKernal`` (external)
    While this parameter is used to determine the accretion kernel in star_maker8.C, there is no choice other than 1 at the moment: Ruffert, ApJ (1994) 427 342 (a typo in the parameter name...).  Default: 0
``StellarWindFeedback`` (external)
    This parameter is used to turn on sink particle creation by star_maker8.C and also its feedback.  Currently implemented are: 1 - protostellar jets along the magnetic fields, 2 - protostellar jets along random directions, 3 - isotropic main sequence stellar wind, 4 - not implemented, 5 - not implemented, 6 - methods 2 and 3 combined.  Default: 0
``StellarWindTurnOnMass`` (external)
    This parameter is used to decide whether mass increase reached the ejection threshold for StellarWindFeedback=1, 2, or 6 in star_maker8.C. Default: 0.1
``MSStellarWindTurnOnMass`` (external)
    This parameter is used to decide whether mass increase reached the ejection threshold for StellarWindFeedback = 3 or 6 in star_maker8.C. Default: 10.0
``BigStarFormation`` (external)
    This parameter is used to turn on sink particle creation by star_maker9.C.  
``BigStarFormationDone`` (external)
    In star_maker9.C, this parameter is used when we do not want to form BigStars any more.
``BigStarSeparation`` (external)
    In star_maker[89].C, if the newly-created sink particle is within a certain distance from the closest pre-existing sink, then add to it rather than creating a new one.
``SinkMergeDistance``
    [not used]
``SinkMergeMass``
    [not used]

.. _radiation_parameters:

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
    support). The following values are used.  For field type 15, see
    Table 3 in `Haardt & Madau (2012)
    <http://adsabs.harvard.edu/abs/2012ApJ...746..125H />`_. Default: 0

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
     15 - Haardt & Madau 2012.

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
    automatically turned on when RadiationFieldType is set to 11. When
    set to 1, it calculates shielding for H/He. See
    ``calc_photo_rates.src`` for more details.  When set to 2, it
    shields only H2 with the Sobolev-like approximation from
    Wolcott-Green et al. (2011).  Default: 0
``RadiationFieldRedshift`` (external)
    This parameter specifies the redshift at which the radiation field
    is calculated.  If a UV radiation background is used in a
    non-cosmological simulation, this needs to be defined. Negative
    redshifts are permitted. Default: (undefined)
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
``RadiativeTransferHubbleTimeFraction`` (external)
    Photon packages are deleted when its associated photo-ionization
    timescale, considering the limit when all photons are absorbed in
    one cell, drops below a fraction (this parameter) of a Hubble
    time.  This parameter can be safely set to 0.01 when ray merging
    is used.  Default: 0.1
``RadiativeTransferFluxBackgroundLimit`` (external)
    When the flux of a photon package drops below a fraction (this
    parameter) of the background radiation field, the ray is deleted.
    Only used with ray merging.  Default: 0.01
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
``RadiativeTransferTimestepVelocityLevel`` (external)
    Limit the ray tracing timestep by a sound crossing time (see
    ``RadiativeTransferTimestepVelocityLimit``) across a
    cell on the level specified with this parameter.  Not used if
    equal to INT_UNDEFINED (-99999).  Default: INT_UNDEFINED
``RadiativeTransferHIIRestrictedTimestep`` (external)
    Adaptive ray tracing timesteps will be restricted by a maximum change of 10% in neutral fraction if this parameter is set to 1.  If set to 2, then the incident flux can change by a maximum of 0.5 between cells.  See Wise & Abel (2011) in Sections 3.4.1 and 3.4.4 for more details.  Default: 0
``RadiativeTransferAdaptiveTimestep`` (external)
    Must be 1 when RadiativeTransferHIIRestrictedTimestep is non-zero.  When RadiativeTransferHIIRestrictedTimestep is 0, then the radiative transfer timestep is set to the timestep of the finest AMR level.  Default: 0
``RadiativeTransferLoadBalance`` (external)
    When turned on, the grids are load balanced based on the number of ray segments traced.  The grids are moved to different processors only for the radiative transfer solver.  Default: 0
``RadiativeTransferHydrogenOnly`` (external)
    When turned on, the photo-ionization fields are only created for hydrogen.  Default: 0
``RadiativeTransferRayMaximumLength`` (external)
    The maximum length that a ray is allowed to travel in box units. Thde default value is 1.7320608 (i.e. sqrt(3.0) so a ray covers the entire periodic region with some doubling up inevitably. Setting it to smaller value will reduce the computational cost.
    Default: 1.7320608
``RadiativeTransferUseH2Shielding`` (external)
    Should H2 self-shielding be used.  Default: True
``RadiativeTransferH2ShieldType`` (external)
    If H2 shielding is turned on then which kind should we use. Setting this value to 0 used the self-shielding fit as per
    Draine & Bertoldi (1996). Setting this value to 1 uses the fit as per Wolcott-Green et al. (2011). Default: 0
``RadiativeTransferH2IIDiss`` (external)
    Should we also account for the photo-dissoication of H2II which occurs for radiation between 0.76eV and 13.6 eV.  Default: True
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

.. _radiative_transfer_fld_implicit_solver:

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

.. _radiative_transfer_fld_split_solver:

Radiative Transfer (FLD) Split Solver Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    These parameters should be placed within the file named in
    ``RadHydroParamfile`` in the main parameter file. All are described in
    detail in the User Guide in ``doc/split_fld``.


``RadHydroESpectrum`` (external)
    Type of assumed radiation spectrum for radiation field, Default: 1.

   ::
 
    1  - T=1e5 blackbody spectrum
    0  - power law spectrum, ( nu / nu_{HI})^(-1.5)` 
    -1 - monochromatic spectrum at frequency h nu_{HI}= 13.6 eV
    -2 - monochromatic spectrum at frequency h nu_{HeI}= 24.6 eV
    -3 - monochromatic spectrum at frequency h nu_{HeII}= 54.4 eV

``RadHydroChemistry`` (external)
    Use of primordial chemistry in computing opacities and
    photo-heating/photo-ionization.  Default: 1. 

   ::

    0 no chemistry
    1 hydrogen chemistry
    3 hydrogen and helium chemistry

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
``RadHydroMaxSubcycles`` (external)
    desired number of FLD time steps per hydrodynamics time step (must
    be greater than or equal to 1). This is only recommended if the
    FLD solver is performing chemistry and heating internally, since
    it will only synchronize with the ionization state at each
    hydrodynamic time step.  When using Enzo's chemistry and cooling
    solvers this parameter should be set to 1 to avoid overly
    decoupling radiation and chemistry.  Default: 1.0.
``RadHydroMaxChemSubcycles`` (external)
    desired number of chemistry time steps per FLD time step.  This
    only applies if the FLD solver is performing chemistry and heating
    internally, instead of using Enzo's built-in routines for this
    task. Default: 1.0.
``RadHydroDtNorm`` (external)
    type of p-norm to use in estimating time-accuracy for predicting
    next time step. Default: 2.0.

   ::

    0  - use the max-norm.
    >0 - use the specified p-norm.
    <0 - illegal.

``RadHydroDtGrowth`` (external)
    Maximum growth factor in the FLD time step between successive
    iterations. Default: 1.1 (10% growth).
``RadHydroDtRadFac`` (external)
    Desired time accuracy tolerance for the radiation field. Default:
    1e20 (unused).
``RadHydroDtGasFac`` (external)
    Desired time accuracy tolerance for the gas energy field.  Only
    used if the FLD solver is performing heating internally.  Default:
    1e20 (unused).
``RadHydroDtChemFac`` (external)
    Desired time accuracy tolerance for the hydrogen I number
    density.  Only used if the FLD solver is performing chemistry
    internally.  Default: 1e20 (unused).
``RadiationScaling`` (external)
    Scaling factor for the radiation field, in case standard
    non-dimensionalization fails. Default: 1.0.
``EnergyCorrectionScaling`` (external)
    Scaling factor for the gas energy correction, in case standard
    non-dimensionalization fails. Default: 1.0.
``ChemistryScaling`` (external)
    Scaling factor for the hydrogen I number density, in case standard
    non-dimensionalization fails. Default: 1.0.
``AutomaticScaling`` (external)
    Enables an heuristic approach in the FLD solver to update the
    above scaling factors internally.  Works well for reioniztaion
    calculations, but is not recommended for problems in which the
    optimal unit scaling factor is known a-priori. Default: 1.0.
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
``RadHydroKrylovMethod`` (external)
    Desired outer linear solver algorithm to use.  Default: 1.

    ::

     0 - Preconditioned Conjugate Gradient (PCG)
     1 - Stabilized Bi-Conjugate Gradient (BiCGStab)
     2 - Generalized Minimum Residual (GMRES)

``RadHydroSolTolerance`` (external)
    Desired accuracy for solution to satisfy linear residual (measured
    in the 2-norm). Default: 1e-8.
``RadHydroMaxMGIters`` (external)
    Allowed number of iterations for the inner linear solver (geometric
    multigrid). Default: 50.
``RadHydroMGRelaxType`` (external)
    Relaxation method used by the multigrid solver. Default: 1.

    ::

     0 - Jacobi
     1 - Weighted Jacobi
     2 - Red/Black Gauss-Seidel (symmetric)
     3 - Red/Black Gauss-Seidel (non-symmetric)

``RadHydroMGPreRelax`` (external)
    Number of pre-relaxation sweeps used by the multigrid solver.
    Default: 1.
``RadHydroMGPostRelax`` (external)
    Number of post-relaxation sweeps used by the multigrid solver.
    Default: 1.
``EnergyOpacityC0``, ``EnergyOpacityC1``, ``EnergyOpacityC2`` (external)
    Parameters used in defining the energy-mean opacity used with
    RadHydroModel 10. Default: [1 1 0].

.. _cosmology_parameters:

Cosmology Parameters
~~~~~~~~~~~~~~~~~~~~

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
``CosmologyOmegaRadiationNow`` (external)
    This is the contribution of all relativistic matter to the energy
    density at the current epoch (z=0), in the same units as above.
    Default: 0.0.
``CosmologyHubbleConstantNow`` (external)
    The Hubble constant at z=0, in units of 100 km/s/Mpc. Default:
    0.701
``CosmologyComovingBoxSize`` (external)
    The size of the volume to be simulated in Mpc/h (at z=0). Default:
    64.0
``CosmologyInitialRedshift`` (external)
    The redshift for which the initial conditions are to be generated.
    Default: 20.0
``CosmologyMaxExpansionRate`` (external)
    This float controls the timestep so that cosmological terms are
    accurate followed. The timestep is constrained so that the relative
    change in the expansion factor in a step is less than this value.
    Default: 0.01
``CosmologyTableNumberOfBins`` (external)
    Conversions between time and redshift are computed by interpolating
    from a numerically integrated table of log(scale factor) vs. time.
    This parameter sets the number of bins in the table. Default: 1000.
``CosmologyTableLogaInitial`` (external)
    This sets the lower bound of the table used to convert between time
    and redshift. This is log10 of the lowest value of the scale factor.
    This value will be automatically adjusted if
    ``CosmologyInitialRedshift`` is set to an earlier time.
    Default: -6.0, (i.e., z = 999,999.)
``CosmologyTableLogaFinal`` (external)
    This sets the upper bound of the table used to convert between time
    and redshift. This is log10 of the highest value of the scale factor.
    This value will be automatically adjusted if
    ``CosmologyFinalRedshift`` is set to a later time.
    Default: 0.0, (i.e., z = 0.)
``CosmologyCurrentRedshift`` (information only)
    This is not strictly speaking a parameter since it is never
    interpreted and is only meant to provide information to the user.
    Default: n/a

.. _massive_black_hole_parameters:

Massive Black Hole Physics Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following parameters are for the accretion and feedback from the
massive black hole particle (``PARTICLE_TYPE_MBH``). Details
are described in Kim, Wise, Alvarez, and Abel (2011).

.. _accretion_physics_parameters:

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

.. _Feedback_physics:

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

.. _shock_finding_parameters:

Shock Finding Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

For details on shock finding in Enzo see :ref:`shock_finding`.

``ShockMethod`` (external)
    This parameter controls the use and type of shock finding. Default: 0
    
    ::

      0 - Off
      1 - Temperature Dimensionally Unsplit Jumps
      2 - Temperature Dimensionally Split Jumps
      3 - Velocity Dimensionally Unsplit Jumps
      4 - Velocity Dimensionally Split Jumps

``ShockTemperatureFloor`` (external)
    When calculating the mach number using temperature jumps, set the
    temperature floor in the calculation to this value.

``StorePreShockFields`` (external)
    Optionally store the Pre-shock Density and Temperature during data output.

``FindShocksOnlyOnOutput`` (external)
    0: Finds shocks during Evolve Level and just before writing out data. 1: Only find shocks just before writing out data.  2: Only find shocks during EvolveLevel. Default: 0

.. _cosmic_ray_two_fluid_model_parameters:

Cosmic Ray Two-Fluid Model Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For details on the cosmic ray solver in Enzo see :ref:`cosmic_rays`.

``CRModel`` (external)
    This parameter turns on the model. Default: 0

    0. Off
    1. On

``CRgamma``
    For CR equation of state. Default: 4.0/3.0 (relativistic, adiabatic gas)

``CRDiffusion`` (external)
    Switches on diffusion of the cosmic ray energy density. Default: 0

    0. Off
    1. On with constant coefficient (``CRkappa``)


``CRkappa`` (external)
    Cosmic ray diffusion coefficient in CGS units (cm^2/s), Default: 0.0. For MW-like galaxies: 1E28.

``CRCourantSafetyNumber`` (external)
    Multiplies CR diffusion timestep, for stability should be <= 0.5. Default: 0.5

``CRFeedback`` (external)
    Specify fraction of star formation feedback energy should be diverted into the cosmic
    ray energy density. implemented ONLY for star_maker3 (feedback method 2). Default: 0.0

``CRdensFloor`` (external)
    Floor in gas density, can be imposed, for speed purposes (default 0.0 = off). Any value
    larger than 0.0 is on with that value as the floor in code units.

.. _conduction:

Conduction
~~~~~~~~~~

Isotropic and anisotropic thermal conduction are implemented using the
method of Parrish and Stone: namely, using an explicit, forward
time-centered algorithm.  In the anisotropic conduction, heat can only
conduct along magnetic field lines.  One can turn on the two types of
conduction independently, since there are situations where one might 
want to use both.  The Spitzer fraction can be also set
independently for the isotropic and anisotropic conduction.  Running a 
cosmological simulation with conduction on can be tricky as the timesteps 
can become very short.  It is recommended that you look carefully at all the 
available conduction parameters.  Additionally, if you intend to run with 
star particles, it is highly recommended that you set the parameter, 
``StarMakerTimeIndependentFormation``.  See the description in 
:ref:`starparticleparameters` for more information.

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
``SpeedOfLightTimeStepLimit`` (external)
    When used, this sets a floor for the conduction timestep to be the local light crossing time (dx / c).  This prevents the conduction machinery from prescribing extremely small timesteps.  While this can technically violate the conduction stability criterion, testing has shown that this does not result in notable differences.  (1 - ON; 0 - OFF)  Default: 0 (OFF).
``ConductionDynamicRebuildHierarchy`` (external)
    Using conduction can often result in the code taking extremely short timesteps.  Since the hierarchy is rebuilt each timestep, this can exacerbate memory fragmentation issues and slow the simulation.  In the case where the conduction timestep is the limiter, the hierarchy should not need to be rebuilt every timestep since conduction mostly does not alter the fields which control refinement.  When this option is used, the timestep calculation is carried out as usual, but the hierarchy is only rebuilt on a timescale that is calculated neglecting the conduction timestep.  This results in a decent speedup and reduced memory fragmentation when running with conduction.  (1 - ON; 0 - OFF)  Default: 0 (OFF).
``ConductionDynamicRebuildMinLevel`` (external)
    The minimum level on which the dynamic hierarcy rebuild is performed.  Default: 0.

.. _other_parameters:

Other Parameters
~~~~~~~~~~~~~~~~

.. _other_external_parameters:

Other External Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

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
``StopSteps``
    Reserved for future use
``CoolDataf0to3``
    Reserved for future use
``StageInput``
    Reserved for future use
``LocalPath``
    Reserved for future use
``GlobalPath``
    Reserved for future use

.. _inline_analysis:

Inline Analysis
~~~~~~~~~~~~~~~

.. _inline_halo_finding:

Inline Halo Finding
^^^^^^^^^^^^^^^^^^^

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
``HaloFinderRunAfterOutput`` (external)
    When turned on, the inline halo finder is run after an output is written.  Default: 0
``HaloFinderLastTime`` (internal)
    Last time of a halo find. Default: 0.

.. _inline_python:

Inline Python
^^^^^^^^^^^^^

``PythonTopGridSkip`` (external)
    How many top grid cycles should we skip between calling python at the top of the hierarchy?  Only works with python-yes in compile settings.
``PythonSubcycleSkip`` (external)
    How many subgrid cycles should we skip between calling python at the bottom of the hierarchy?
``PythonReloadScript`` (external)
    Should "user_script.py" be reloaded in between Python calls?
``NumberOfPythonCalls`` (internal)
    Internal parameter tracked by Enzo
``NumberOfPythonTopGridCalls`` (internal)
    Internal parameter tracked by Enzo
``NumberOfPythonSubcycleCalls`` (internal)
    Internal parameter tracked by Enzo

.. _other_internal_parameters:

Other Internal Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

``TimeLastDataDump`` (internal)
    The code time at which the last time-based output occurred.
``TimeLastInterpolatedDataDump`` (internal)
    The code time at which the last interpolated data dump occurred.
``CycleLastDataDump`` (internal)
    The last cycle on which a cycle dump was made
``SubcycleLastDataDump`` (internal)
    The last cycle on which a subcycle dump was made
``TimeLastMovieDump`` (internal)
    The code time at which the last movie dump occurred.
``TimeLastTracerParticleDump`` (internal)
    The code time at which the last tracer particle dump occurred.
``TimeLastRestartDump``
    Reserved for future use.
``TimeLastHistoryDump``
    Reserved for future use.
``CycleLastRestartDump``
    Reserved for future use.
``CycleLastHistoryDump``
    Reserved for future use.
``InitialCPUTime``
    Reserved for future use.
``InitialCycleNumber`` (internal)
    The current cycle
``SubcycleNumber`` (internal)
    The current subcycle
``DataDumpNumber`` (internal)
    The identification number of the next output file (the 0000 part of
    the output name). This is used and incremented by both the cycle
    based and time based outputs. Default: 0
``MovieDumpNumber`` (internal)
    The identification number of the next movie output file. Default: 0
``TracerParticleDumpNumber`` (internal)
    The identification number of the next tracer particle output file. Default: 0    
``RestartDumpNumber``
    Reserved for future use.
``HistoryDumpNumber``
    Reserved for future use.
``DataLabel[#]`` (internal)
    These are printed out into the restart dump parameter file. One
    Label is produced per baryon field with the name of that baryon
    field. The same labels are used to name data sets in HDF files.
``DataUnits[#]`` 
    Reserved for future use.
``VersionNumber`` (internal)
    Sets the version number of the code which is written out to restart
    dumps.

.. _problem_type_parameters:

Problem Type Parameters
~~~~~~~~~~~~~~~~~~~~~~~

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
1            :ref:`shocktube_param`
2            :ref:`wavepool_param`
3            :ref:`shockpool_param`
4            :ref:`doublemach_param`
5            :ref:`shockinabox_param`
6            :ref:`implosion_param`
7            :ref:`sedovblast_param`
8            :ref:`khinstability_param`
9            :ref:`noh_param`
10           :ref:`rotatingcylinder_param`
11           :ref:`radiatingshock_param`
12           :ref:`freeexpansion_param`
14           :ref:`rotatingsphere_param`
20           :ref:`zeldovichpancake_param`
21           :ref:`pressurelesscollapse_param`
22           :ref:`adiabaticexpansion_param`
23           :ref:`testgravity_param`
24           :ref:`sphericalinfall_param`
25           :ref:`testgravitysphere_param`
26           :ref:`gravityequilibriumtest_param`
27           :ref:`collapsetest_param`
28           :ref:`testgravitymotion_param`
29           :ref:`testorbit_param`
30           :ref:`cosmologysimulation_param`
31           :ref:`galaxysimulation_param`
35           :ref:`shearingbox_param`
36           Shearing Box 2D Simulation
37           Stratifeid Shearing Box Simulation
40           :ref:`supernovarestart_param`
50           :ref:`photontest_param`
51           Photon Test Restart
59           :ref:`stochastic_forcing_param`
60           :ref:`turbulence_param` 
61           :ref:`protostellar_param` 
62           :ref:`coolingtest_param`
63           One Zone Free Fall Test
70           Conduction Test with Hydro Off
71           Conduction Test with Hydro On
72           Conduction Bubble Test
73           Conduction Cloud Test
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
209          MHD 1D Waves
210          MHD Decaying Random Magnetic Fields
250          :ref:`cr_shocktube_param`
300          :ref:`poissonsolver_param`
400          :ref:`rhdtest1_param`
401          :ref:`rhdtest2_param`
402          :ref:`rhdtest3_param`
403          :ref:`rhdtest4_param`
404/405      :ref:`rhdtest5_param`
410/411      :ref:`rhdtest10_param`
412          :ref:`rhdtest12_param`
413          :ref:`rhdtest13_param`
414/415      :ref:`rhdtest14_param`
450-452      Free-streaming radiation tests (to be removed)
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
``CollapseTestInitialFractionHII`` (external)
    Initial HII fraction in the domain except for the spheres.
    Default: 1.2e-5
``CollapseTestInitialFractionHeII`` (external)
    Initial HeII fraction in the domain except for the spheres.
    Default: 1e-14
``CollapseTestInitialFractionHeIII`` (external)
    Initial HeIII fraction in the domain except for the spheres.
    Default: 1e-17
``CollapseTestInitialFractionHM`` (external)
    Initial H- fraction in the domain except for the spheres.
    Default: 2e-9
``CollapseTestInitialFractionH2I`` (external)
    Initial H2I fraction in the domain except for the spheres.
    Default: 2e-20
``CollapseTestInitialFractionH2II`` (external)
    Initial H2II fraction in the domain except for the spheres.
    Default: 3e-14
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
``CollapseTestSphereHIIFraction`` (external)
    Initial HII fraction of the sphere.  Default: 1.2e-5
``CollapseTestSphereHeIIFraction`` (external)
    Initial HeII fraction of the sphere.  Default: 1e-14
``CollapseTestSphereHeIIIFraction`` (external)
    Initial HeIII fraction of the sphere.  Default: 1e-17
``CollapseTestSphereHMFraction`` (external)
    Initial H- fraction of the sphere.  Default: 2e-9
``CollapseTestSphereH2IFraction`` (external)
    Initial H2I fraction of the sphere.  Default: 2e-20
``CollapseTestSphereH2IIFraction`` (external)
    Initial H2II fraction of the sphere.  Default: 3e-14
``CollapseTestSphereInitialLevel`` (external)
    Failed experiment to try to force refinement to a specified level.
    Not working. Default: 0
``CollapseTestWind`` (external)
    Boolean flag. Type: integer. This parameter decides if there is wind (inflow boundary). Default: 0 (FALSE)
``CollapseTestWindVelocity`` (external)
    When using inflow boundary, this is the inflow velocity. Default: 0.

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
``GalaxySimulationRPSWind`` (external)
    This flag turns on the ram pressure stripped (RPS) wind in the
    GalaxySimulation problem and sets the mode.  0 = off, 1 = on with
    simple constant wind values, 2 = on with RPS values set from a
    file with the name ICMinflow_data.in.  For the file input case,
    the file should consist of a set of lines with each line
    specifying a 6 columns consisting of time, wind density, wind
    temperature, wind x/y/z velocity.  All units in the file are
    assumed to be CGS and wind values are applied at the time
    indicated to the corner of the box, with linear interpolation
    between key frames.  See Salem et al. (2015) for a worked example.
    Default: 0
``GalaxySimulationRPSWindShockSpeed`` (external)
    This is speed of the RPS driven shock (which differs from the
    wind velocity), to be used to determine where and when to apply
    the appropriate wind boundary condition on the boundary.  Code units.
    Default: 0.0
``GalaxySimulationRPSWindDelay`` (external)
    This is a delay (in code units) for the RPS wind to be applied
    (for example to give time for the galaxy to relax).
    Default: 0.0
``GalaxySimulationRPSWindDensity`` (external)
    For case 1, this is the density of the RPS wind, in code units.
    Default: 1.0
``GalaxySimulationRPSWindtotalEnergy`` (external)
    For case 1, this is the total energy of the RPS wind, in code units.
    Default: 1.0
``GalaxySimulationRPSWindPressure`` (external)
    For case 1, this is the pressutre of the RPS wind (unused).
    Default: 1.0
``GalaxySimulationRPSWindVelocity`` (external)
    For case 1, This is the wind velocity (code units)
    Default: 0 0 0
``GalaxySimulationRPSWindPreWindDensity`` (external)
    This is the density applied to the boundary before the wind arrives.
    Default: 1.0
``GalaxySimulationRPSWindPreWindTotalEnergy`` (external)
    This is the total energy applied to the boundary before the wind arrives.
    Default: 1.0
``GalaxySimulationRPSWindPreWindVelocity`` (external)
    This is the velocity vector applied to the boundary before the
    wind arrives.
    Default:

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
``PhotonTestDensityFilename`` (external)
    Filename of an external density field in HDF5 format.  The file
    should only have one dataset. Default: (undefined)
``PhotonTestHIIFractionFilename`` (external)
    Filename of an external HII fraction field in its own HDF5 format.
    The file should only have one dataset.  Default: (undefined)
``PhotonTestHeIIFractionFilename`` (external)
    Filename of an external HeII fraction field in its own HDF5 format.
    The file should only have one dataset.  Default: (undefined)
``PhotonTestHeIIIFractionFilename`` (external)
    Filename of an external HeIII fraction field in its own HDF5 format.
    The file should only have one dataset.  Default: (undefined)
``PhotonTestTemperatureFilename`` (external)
    Filename of an external temperature field in its own HDF5 format.
    The file should only have one dataset.  Default: (undefined)

.. _stochastic_forcing_param:

Turbulence Simulation with Stochastic Forcing (59)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    Type: integer. 1--Calculate the amount of cold gas around the SMBH and remove it at the rate of 2*Mdot; 2--Calculate Mdot based on the amount of cold gas around the SMBH; 3--Calculate Mdot similar to 2 but change ClusterSMBHJetDim periodically (period = ClusterSMBHJetPrecessionPeriod); 4--Calculate Mdot within Bondi radius (only use this when Bondi radius is resolved); 0--off (do not remove cold gas). Default: 1.
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
``ClusterSMBHDiskRadius`` (external)
    The size of the accretion zone in kpc. Default: 0.5
``ClusterSMBHBCG`` (external)
    The stellar component of the Perseus BCG (in cluster simulations) or the elliptical galaxies (in simulations of isolated elliptical galaxies). Default: 1.0
``ClusterSMBHMass`` (external)
    The mass of the SMBH of the Perseus BCG (in cluster simulations) or the elliptical galaxies (in simulations of isolated elliptical galaxies). Default: 0
``EllipticalGalaxyRe`` (external)
    Re is the radius of the isophote enclosing half of the galaxy's light. In Herquist profile, a=Re/1.8153. Default: 0
``OldStarFeedbackAlpha`` (external)
    Mass ejection rate from evolved stars in the unit of 10^{-19} s^{-1}. It is typically within a factor of 2 of unity. Default: 0
``SNIaFeedbackEnergy`` (external)
    Energy feedback from evolved stars (Type Ia SN). Default: 1.0

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

This problem type sets up many common 2D hydro and MHD problem types.
Many of them can be run also without MHD despite the name. Which problem is done is controled by
MHD2DProblemType which can vary from 0 to 16 so far.

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
    0: Raleigh-Taylor, 1: MHD rotor (Toth 2000, JCompPhys 161, 605.), 2: MHD blast wave (Gardiner and Stone 2005, JCompPhys. 205, 509), 3: MHD Kelvin-Helmholtz (Gardiner & Stone 2005), 4: Another MHD Kelvin Helmholtz, 5: Shock-vortex interaction (Rault, Chiavassa & Donat, 2003, J. Scientific Computing, 19, 1.), 6: Sedov-Taylor Blast Wave (Fryxell et al. 2000, ApJS, 131, 273), 7: Cylindrical Sedov-Taylor Blast Wave (Fryxell et al. 2000), 8: Like MHD2DProblemType = 5 but with a small perturbation upstream of the shock to test odd even coupling of Reimann Solvers, 9: Smoothed Kelvin Helnholtz problem (Robertson, Kravtsov, Gnedin, Abel & Rudd 2010, MNRAS, 401), 10: A modified Raleigh-Taylor problem, 11: Uniform density with sinusoidal shear velocity (Compare to rpSPH tests in Abel 2012), 12: Experimental test, 13: Exploratory blob test, 14: Wengen 2 test to study colliding flows with very soft equations of state, 15: Another experiment with B-fields, 16: A validated non-linear Kelvin Helmholtz test (Lecoanet, McCourt, Quataert, Burns, Vasil, Oishi, Brown, Stone, & O’Leary 2015 preprint)
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

.. _cr_shocktube_param:

CR Shock Tube (250: unigrid and AMR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Very similar to normal shock tube (see problem 1) but includes CR
    component.  See `Salem, Bryan & Hummels (2014)
    <http://adsabs.harvard.edu/abs/2014ApJ...797L..18S>`__ for discussion.

    In addition to the regular shock tube parameters, we add:

``HydroShockTubesLeftCREnDensity``, ``HydroShockTubesRightCREnDensity`` (external)
    The initial CR energy density on the left and right sides.
    Default: 1.0 for each value.

``HydroShockTubesCenterDensity``, ``HydroShockTubesCenterPressure``,
``HydroShockTubesCenterVelocityX``,
``HydroShockTubesCenterVelocityY``,
``HydroShockTubesCenterVelocityZ``,
``HydroShockTubesCenterCREnDensity`` (external)

    In addition to setting a shock tube with two constant regions,
    this version also allows for three constant regions, 
    with a Center region in addition to the Left and Right regions.
    Finally, there are two special cases -- if
    HydroShockTubesCenterCREnDensity is set to 123.4, then the central
    region will be set to a ramp between the left and right regions,
    and if HydroShockTubesCenterCREnDensity is set to 567.8, then a
    gaussian CR energy density is initialized (these problems were set
    up to test the CR diffusion).

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


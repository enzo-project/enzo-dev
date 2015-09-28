.. _io_parameters:

I/O Parameters
--------------

General 
^^^^^^^

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
``OutputOnPop3Feedback`` (external)
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


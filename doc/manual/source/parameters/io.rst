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

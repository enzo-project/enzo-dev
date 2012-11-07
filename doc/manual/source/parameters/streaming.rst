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

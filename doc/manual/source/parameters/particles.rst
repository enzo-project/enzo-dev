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
``AddParticleAttributes`` (external)
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


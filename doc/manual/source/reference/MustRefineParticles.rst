.. _MustRefineParticles:

Must Refine Particles
=====================

Must refine particles within Enzo are dark matter particles that force the code to refine the hydro mesh around their position up to a level specified by ``MustRefineParticlesRefineToLevel``.  They are useful for cosmological zoom simulations and can be used to optimize refinment within the high resoluiton region (see :ref:`MRPNestedGrid`).  They can also be used in more general cases, as demonstrated in the TestOrbitMRP problem within the run directory.

Within the code, must refine particles act identically to dark matter particles, except for the refinement criteria they introduce.  The numerical value for the particle type of must refine particles is 4 (versus 0 for ordinary dark matter particles).  If you wish to use must refine particles in your simulation, particles should be given PARTICLE_TYPE_MUST_REFINE when they are inialized at start-up or when they are created during the simulation. 

The code will flag cells around must refine particles following a cloud-in-cell (CIC) algorithm.  The default cloud size is one cell, so in the general case, each must refine particle will flag a total of 8 hydro cells.  The cloud size can be adjusted by the user with the internal parameter ``ParticleBufferSize``.  Note that if ghost cells are flagged for refinement, this information is not communicated to the grid where the ghost cells are active.  




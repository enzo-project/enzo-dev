Dynamic Refinement Regions in Nested Cosmology Simulations
==============================================

Traditionally, cosomological zoom-in simulations run with Enzo have used nested, rectilinear regions that are static and surround the region of maximum refinement.  These static regions of high and intermediate refinement cover the volume where the highest-resolution mass dark matter particles are placed in the initial conditions at high redshift.  As particles of different masses move during the simulation, they will move to regions with either higher or lower resolution than the region where they were generated during the problem initialization. 

There are several strategies available for combating this problem, some of which are described below.  In addition, it may be sometimes desireable to restrict the region of maximum refinement to a non-rectangular volume or to a subset of the highest resolution particles and some strategies for this are also described below.  Some of these options (but not all) are for use with the initial conditions generator MUSIC (ref).


Non-static Nested Grids
------------

With this method, Enzo no longer uses the static grids that are used during the simulation initialization to control refinement on intermediate nested levels.  Refinement is adjusted to a minimum level appropriate for a dark particle's mass each time the hierarchy is rebuilt.  This option can be activated by setting ``MustRefineParticlesCreateParticles = 3``.  ``CosmologySimulationGrid`` parameters are still necessary because they are used during the problem initialization because particles of different masses can still only be initialized in rectangular regions on start-up.  However, the user should not set ``RefineRegionLeftEdge`` or ``RefineRegionRightEdge``.

This method works by identifying dark matter particles whose particle densities are less than the cosmic dark matter density in the simulation as defined by ``CosmologySimulationOmegaCDMNow`` and ``OmegaMatterNow``.  Users should recall that particle 'masses' stored by Enzo are actually particle densities (see :ref:`EnzoParticleMass`).  When a dark matter particle resides on a grid whose resolution is equal to the resolution of the particle's generation grid, the particle's density will be exactly ``CosmologySimulationOmegaCDMNow/OmegaMatterNow``.  If a particle is on a grid of coarser resolution, the particle's density will be less than this fraction; refinement is triggered around particles when this is the case.  Refinement is done in a CIC (cloud-in-cell) region centered on the particle that is one cell wide and therefore flags a total of 8 cells surrounding the particle.  This flagging continues around the particle until the particle's density equals ``CosmologySimulationOmegaCDMNow/OmegaMatterNow``.

The user should note that this method gaurentees that a particle will never be on a grid coarser than its generation grid, but a consequence of the CIC flagging method and Enzo's grid deconstrution method is that paricles can be on grids with finer resolution.  This method therefore does not prevent contamination of high resolution regions by low resolution paricles.

MUSIC Ellipsoidal Masking
-------------

The nested cosmological initial conditions generator MUSIC has the capability to identify high-resolution particles that lie within an ellipsoidal 'masked' region within the rectangular region containing the highest resolution particles created.  MUSIC does this with files called ``RefinementMask.x`` that contain an integer flag that tells the user whether a particle is in the masked region or not.  AMR can be restricted to the volume containing these particles by setting ``MustRefineParticlesCreateParticles = 2`` or ``3``,``CosmologySimulationParticleTypeName = RefinementMask`` and ``CellFlaggingMethod = 8``.

This refinement restriction is done by flagging particles within the ellipsodial masked region as `must refine particles.'  Must refine particles in Enzo are particles with a special particle type that tells the code to flag cells around the particle for further refinement based on a single-cell CIC cloud centred on the particle.  Therefore each particle will in general flag 8 cells.  If must refine particles are created at start-up with this method, the code also restricts other AMR criteria to only the volume containing these particles, which is not generally a feature of must refine particles that are used for other purposes in Enzo.  
It is important to note that if the CIC-flagging flags cells in the ghost zone of a grid, then these flags are not communicated to the grids where these cells are active.

If ``MustRefineParticlesCreateParticles = 2``, then traditional static nested regions are used on intermediate refinement levels.  If ``MustRefineParticlesCreateParticles = 3``, then non-static nested refinement as described above is conducted.


Non-ellipsoidal Masking
--------------

This method uses traditional static grids on intermediate levels, but within the highest level static region, maximum refinement is restricted to a volume surrounding a subset of the highest resolution dark matter paricles by tagging them as `must refine particles' either from a list of particle IDs or within a specified rectangular box.  It is activated by setting ``MustRefineParticlesCreateParticles = 1``.  

These particles can be selected in one of two ways: (1) by specifiying a rectangular volume with the parameters ``MustRefineParticlesRegionLeftEdge`` and ``MustRefineParticlesRightEdge``; or (2) by providing a list of particle IDs.  If ``MustRefineParticlesRegionLeftEdge`` and ``MustRefineParticlesRegionRightEdge`` are not set, but ``MustRefineParticlesCreateParticles = 1``, then the code looks for an ascii file called ``MustRefineParticlesFlaggingList.in`` within the run directory.

The list of particle IDs can be obtained from a prior dark-matter-only simulation of the refined initial conditions.  The user should be cautitious when doing this, because particle IDs are not contained within the initial conditions; they are assigned during the problem initialization.  Depending on the way in which each simulation is initialized, particle IDs are not gaurenteed to be identical.  These methods do not rely on the MUSIC initial conditions generator.


RefineRegionAutoAdjust
-----------

One simple strategy is to allow the nested rectangular regions to shrink with the comological collapse of structures.  This can be activated with the parameter ``RefineRegionAutoAdjust``.  With this method, however, areas of increased refinement retain an overall rectangular shape.


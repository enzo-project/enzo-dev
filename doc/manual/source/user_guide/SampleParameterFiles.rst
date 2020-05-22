Sample Parameter Files from Publications
========================================

Here is a collection of parameter files that are equivalent or similar to ones
used in publications.  Each example has a description and motivation behind the choice of
important parameters.

Wise et al. (2014)
------------------
(`Publication link
<https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.2560W/abstract>`_)
(:download:`Parameter file <./samples/wise2014_zoomin_update.enzo>`)
(:download:`MUSIC parameter file <./samples/SG64-L2.conf>`)
(`Lagrangian data file <https://www.dropbox.com/s/x53m3a1566vhotq/SG64-L2-Lagrangian.dat.gz?dl=1>`_ 622kB)

This parameter file will set up a simulation **similar** to the paper's
simulation.  Some parameter choices are updated to be consistent with papers
written after the original simulation was run in 2012. It models Pop III star
formation and the transition to galaxy formation.  While the production
simulation was a 256\ :sup:`3` AMR everywhere simulation, this sample file is a
zoom-in calculation (the user needs to generate the ICs with
:ref:`MUSIC <CosmologicalInitialConditions>`) focused on a
single 10\ :sup:`8` solar mass halo at z = 10.

Hierarchy setup
^^^^^^^^^^^^^^^
1 Mpc comoving box with a 64\ :sup:`3` root grid resolution and 2 nested grids

* ``MaximumRefinementLevel = 14``: maximal comoving resolution of ~1 pc
* ``MaximumParticleRefinementLevel = 11``: to remove the discreteness of DM particles, smooth them over ~1 proper pc below which baryons become dominant
* ``CellFlaggingMethod = 2 4 6 8``: Refine by baryon/DM mass, Jeans length, and must-refine cosmological particles
* ``RefineByJeansLengthSafetyFactor = 8``: Resolve the local Jeans length by at least 8 cells
* ``MustRefineParticlesCreateParticles = 3``: Use a :ref:`dynamical nested region <MRPNestedGrid>` that's specified by the particle types given by MUSIC
* ``MustRefineParticlesRefineToLevel = 2``: Refine the innermost region to at least level 2

Hydrodynamics
^^^^^^^^^^^^^
* ``RiemannSolver = 4``: use the HLLC solver for more stability in supernova blastwaves
* ``RiemannSolverFallback = 1``: use the HLL solver when HLLC fails

Cooling and chemistry
^^^^^^^^^^^^^^^^^^^^^
* ``MultiSpecies = 3``: 12-species primordial radiative cooling
* ``MetalCooling = 1``: tabulated metal cooling
* ``grackle_data_file = cloudy_metals_2008_3D-lwb.h5``: metal cooling with a Lyman-Werner (LW) background from `Qin+ (2020) <https://arxiv.org/abs/2003.04442>`_
* ``UVbackground = 1``: Use the LW background the table above
* ``H2_self_shielding = 1``: `Wolcott-Green et al. (2011) <https://ui.adsabs.harvard.edu/abs/2011MNRAS.418..838W/abstract>`_ model

Star formation and feedback
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Metal-free and metal-enriched star formation with radiative and supernova
feedback

* ``StarParticleCreation = 40``: Use model #3 (Pop III stars) and #5 (radiating star clusters)
* ``StarParticleFeedback = 40``: Use model #3 and #5

Pop III particles
"""""""""""""""""
* ``PopIIIOverDensityThreshold = -1e6``: The negative indicates units in particles per cm\ :sup:`3`. This number density is consistent with Pop III formation at the resolution level of 0.1 proper pc.
* ``PopIIIH2CriticalFraction = 1e-3``: This is an appropriate H\ :sub:`2` fraction at the above density threshold
* ``PopIIIMetalCriticalFraction = 1.295e-6``: Use a critical transition metallicity of 10\ :sup:`-4` of solar
* ``PopIIIUseHypernova = 0``: Use a typical 10\ :sup:`51` erg core-collapse supernova energy for stars between 11 and 40 solar masses
* ``PopIIIInitialMassFunction = 1``: Use a Pop III IMF defined in Wise et al. (2014)
* ``PopIIIStarMass = 20``: Use a characteristic mass of 20 solar masses in the above IMF. Motivated by `Hirano et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015MNRAS.448..568H/abstract>`_ and `Susa (2019) <https://ui.adsabs.harvard.edu/abs/2019ApJ...877...99S/abstract>`_

Star cluster particles
""""""""""""""""""""""
* ``StarClusterMinDynamicalTime = 3e+06``: Motivated by `Tan et al. (2006) <https://ui.adsabs.harvard.edu/abs/2006ApJ...641L.121T/abstract>`_ who found that molecular clouds have dynamical times around 700 kyr and form stars over several dynamical times. There could be newer references.
* ``StarClusterIonizingLuminosity = 1.9e+46``: Use the average over 20 Myr from the metal-poor binary model of `Rosdahl et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.479..994R/abstract>`_
* ``StarClusterSNEnergy = 1e49``: Calculated from a Salpeter IMF, i.e. 1 SN per 100 solar masses of stars
* ``StarClusterFormEfficiency = 0.07``: Motivated by `Krumholz et al. (2007) <https://ui.adsabs.harvard.edu/abs/2005ApJ...630..250K/abstract>`_. This converts 7% of the cold gas within the star forming cloud into stars. There could be newer references.
* ``StarClusterMinimumMass = 1000``: (in solar masses) equivalent to the smallest observed star clusters

Radiative transfer
^^^^^^^^^^^^^^^^^^
Uses the adaptive ray tracing machinery that is sourced by both types of star
particles

* ``RadiativeTransferOpticallyThinH2 = 1``: Source a (1/r\ :sup:`2`) LW radiation field from the star particles.  Important for star formation regulation
* ``RadiativeTransferRadiationPressure = 1``: Consider momentum transfer between the ionizing photons and absorbing gas. Plays a role at small distances in dense gas and high flux
* ``RadiativeTransferSourceClustering = 1``: Use ray merging to conserve memory (and maybe speed it up)
* ``RadiativeTransferPhotonMergeRadius = 3.0``: Merge rays at 3 times the separation between sources that are paired in a tree

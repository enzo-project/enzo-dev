.. _starparticleparameters:

Star Formation and Feedback Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For details on each of the different star formation methods available in Enzo see :ref:`star_particles`.

General
^^^^^^^

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
    for details. Default 0.0
``StarMakerExplosionDelayTime`` (external)
    Only valid for ``StarParticleFeedback`` method = 14.  If set to a positive value, energy,
    metals and mass from the particle are injected in a single timestep that is delayed from
    the particle creation time by this amount.  This value is in units of Myrs.  If set
    to a negative value, energy, mass and metals are injected gradually in the same way as is
    done for ``StarParticleFeedback`` method = 1.  Default -1.

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
    Set to 1 to trace the metals expelled from supernovae. Default: 0.
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

Massive Black Hole Physics Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following parameters are for the accretion and feedback from the
massive black hole particle (``PARTICLE_TYPE_MBH``). Details
are described in Kim, Wise, Alvarez, and Abel (2011).

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


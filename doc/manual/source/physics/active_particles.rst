.. _active_particles:


Active Particles
============================================

Active Particles (APs) were introduced into Enzo in version 2.6. They were originally part of the
of the Enzo-3.0 effort and were developed by a number of authors including Nathan Goldbaum,
Matt Turk, John Wise, John Regan, Oliver Hahn and others. Their design philosophy was that the
APs would be highly extensible, have robust feedback capabilities and be very much object orientated.
They were designed to replace the star object class. 

Using Active Particles
______________________

To implement the active particles framework in any Enzo run is very simple. Add the following line to the
parameter file that you are using to (re)start Enzo.
``AppendActiveParticleType = <YourActiveParticleType>``
where ``<YourActiveParticleType>`` is one of the currently available AP types. 

Types of Active Particles
_________________________

At the moment there are a number of AP type implementations. However, with the exception of the "SmartStar"
AP type none of the implementations have been (robustly) tested.

* SmartStar
* Accreting Particle
* Cen-Ostriker
* GalaxyParticle
* Kravtsov
* PopIII
* RadiationParticle
* SpringelHernquist
* Skeleton

The ``Skeleton`` particle is simply an example particle heavily commented so allow a user/developer to develop their
own AP. All the other particles with the exception of the ``AccretingParticle`` and the ``SmartStar`` particle were ported from the
Star Object particle implementation described in :ref:`star_particles`.


Current Limitations of Active Particles
_______________________________________

Most of the AP types have not been tested - though in principle do work correctly. Caution and some testing of the particle types is advisable for now.
Multiple AP types can not currently be run together. 

SmartStar Active Particle Type
______________________________

The ``SmartStar`` particle is built on top of the ``AccretingParticle`` type with additional feedback, accretion and physical
particle type included. The ``SmartStar`` particle was designed to be a single particle type that could adjust to the environment
in which it finds itself. Currently it can represent a PopIII star,
a super-massive star or a black hole. However, there is no inherent limit to the physical object it can represent.

To allow for the creation of a ``SmartStar`` particle in an Enzo simulation the following line must be included in the parameter file:
``AppendActiveParticleType = SmartStar``

Once Enzo reads that line in the parameter file then the AP framework will be engaged and the ``SmartStar`` particle initialised.
The following parameters are currently enabled for the ``SmartStar``

SmartStar Feedback
^^^^^^^^^^^^^^^^^^^^^^

``SmartStarFeedback`` (external)
    This is the master feedback parameter. Set this to 0 and all feedback
    is turned off. It's a master switch. Switch it to 1 and then feedback is on but needs to
    be fine grained by more detailed parameters below. 
    Default: 0

``SmartStarStellarRadiativeFeedback`` (external)
    This parameter controls whether stellar feedback is activated or not. For feedback from PopIII or SMSs then this needs to be on.
    The stellar radiative feedback is divided up into 5 energy bins. The energy bins have energies of 2.0 eV, 12.8 eV, 14.0 eV, 25.0 eV
    and 200 eV. The fraction of energy assigned to each bin is determined using the PopIII tables from Schaerer et. al 2002 Table 4.
    The spectrum for a PopIII star and SMS are different. For a PopIII star a spectrum for a 40 Msolar star is assumed and
    weighted accordingly. For a SMS a 1000 Msolar star is assumed and weighted accordingly.
    Future improvements to the SEDs employed here are under active investigation. 
    Default: 0

``SmartStarBHFeedback`` (external)
    This is a master switch on black hole feedback. Must be turned on if you want black hole feedback. Default: 0

``SmartStarBHRadiativeFeedback`` (external)
    This parameter controls whether black hole radiative feedback gets turned on or not. When turned on the radiative
    feedback from a black hole depends both on the mass of the black hole and the accretion rate onto the black hole. Both of these
    quantities are captured and stored as part of the ``SmartStar``. Details of the SED used can be found in the appendix of
    https://arxiv.org/abs/1811.04953 and is made from assuming a multi-colour disk for the accretion disk and a corona fit by a
    power law.  The radiation emitted by the accretion disk is hard-coded is be emitted by 5 bins with energies of
    2.0 eV, 12.8 eV, 19.1 eV, 217.3 eV and 5190 eV. The fraction of energy assigned to each bin is then determined by the mass of the
    black hole and the associated accretion rate at a given time. The formulism is valid for black hole masses between 1 Msolar and
    1e9 Msolar and for accretion rates between 1e-6 Msolar/yr and 1e3 Msolar/yr.  Default: 0

``SmartStarBHThermalFeedback`` (external)
    This parameter controls whether the black hole thermal feedback gets turned. Thre thermal energy is generated by feedback through the
    accretion process. If this is turned on then the ``SmartStarBHRadiativeFeedback`` should presumably be turned off unless you have a
    good reason to include both. The efficiency, epsilon, depends on both the spin of the black hole and the ISCO oribit. In order to
    calculate this accurately Eqn 32 from Abromowicz & Fragile (2013) is used. See ActiveParticle_SmartStar.h. The feedback is released
    iostropically in a sphere surrounding the SmartStar particle.  Default: 0

``SmartStarBHJetFeedback`` (external)
    The methodology for this algorithm is based on that of Kim et. al (2011) (https://arxiv.org/pdf/1106.4007.pdf). Jets can be activiated
    when a spinning black hole is accreting. The jets are bipolar and are set along the angular momentum vector of the SmartStar. The
    velocity of the jet(s) is set by a separate parameter below. No other parameters need to be set to activate the jet.  Default: 0

``SmartStarEddingtonCap`` (external)
    This parameter allows for accretion onto the SmartStar to be capped at the Eddington limit. I see no good physical reason for doing this in
    general.  Default: 0

``SmartStarSpin`` (external)
    The dimensionless spin of the SmartStar particle. This is a very uncontrained parameter and cannot be readily computed on the fly. This parameter
    should be set if you want to have jet feedback. Setting this is zero and turning on jet feedback wouldn't make sense. The default is set to be 0.7 and
    this is probably reasonable. 
    Default: 0.7

``SmartStarSMSLifetime`` (external)
    This is the lifetime for a supermassive star. After this time has elapsed a SmartStar particle which is behaving like a SMS will collapase
    directly into a black hole with no supernova event. Default: 1e6

``SmartStarJetVelocity`` (external)
    The velocity that the jets are ejected at. Typically jets are observed to travel at a substantial fraction of the speed of light -
    especially those ejected during periods of high accretion. However,
    as mass gets entrained on the jet it slows down. The units of this parameter are as a fraction of the speed of light. Default: 0.1

``SmartStarFeedbackJetsThresholdMass`` (external)
    Jets are only ejected once this amount of mass is available for ejected after an accretion event. Therefore, if there is very limited
    accretion and this parameter is set high then jets will be very infrequent. In units of solar masses. Default: 1.0

``SmartStarSuperEddingtonAdjustment`` (external)
    As accretion rates exceed the canonical Eddington rate the radiative efficiency of the feedback changes. We use the fits from Madau et al. (https://arxiv.org/pdf/1402.6995.pdf)
    to adjust the efficiency when accretion enters the super-critical regime. The fits are based on the slim-disk model of accretion which generate inefficient feedback.
    Default: 1







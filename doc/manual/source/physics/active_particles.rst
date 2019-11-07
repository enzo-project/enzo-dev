.. _active_particles:


Active Particles
================

Active Particles (APs) were introduced into Enzo in version 2.6. They were originally part
of the Enzo-3.0 effort and were developed by a number of authors including Nathan Goldbaum,
Matt Turk, John Wise, John Regan, Oliver Hahn, Greg Meece, Brian Crosby and others. Their design philosophy was that the
APs would be highly extensible, have robust feedback capabilities and be very much object orientated.
They were designed to replace the star object class.

Compared to the star objects the APs can more easily connect with the radiative transfer solver thus making the
framework more extensible in that context. Other feedback mechanisms (e.g. mechanical and thermal) are
also easily implemented. 

Using Active Particles
______________________

To implement the active particles framework in any Enzo run is very simple. Add the following line to the
parameter file that you are using to (re)start Enzo.
``AppendActiveParticleType = <YourActiveParticleType>``
where ``<YourActiveParticleType>`` is one of the currently available AP types.

For relevant parameters, please also see :ref:`active_particles_parameters`.

Types of Active Particles
_________________________

At the moment there are a number of AP type implementations. However, with the exception of the "SmartStar" and the
"Cen-Ostriker" AP type none of the implementations have been (robustly) tested.

* SmartStar
* Accreting Particle
* Cen-Ostriker
* GalaxyParticle
* Kravtsov
* PopIII
* RadiationParticle
* SpringelHernquist
* Skeleton

The ``Skeleton`` particle is simply an example particle heavily commented to allow a user/developer to develop their
own AP. All the other particles with the exception of the ``AccretingParticle`` and the ``SmartStar`` particle were ported from the
Star Object particle implementation described in :ref:`star_particles`.


Current Limitations of Active Particles
_______________________________________

Most of the AP types have not been tested - though in principle do work correctly. Caution and some testing of the particle
types is advisable for now.

**Multiple AP types can not currently be run together.** This isn't a fundamental limitation. In principle multiple APs can work
together without difficulty. Some communication work needs to be undertaken to make this work. 


SmartStar Active Particle Type
______________________________

The ``SmartStar`` particle is built on top of the ``AccretingParticle`` type with additional feedback and accretion protocols attached.
The ``SmartStar`` particle was designed to be a single particle type that could adjust to the environment
in which it finds itself. Currently it can represent a PopIII star,
a super-massive star or a black hole. However, there is no inherent limit to the physical object it can represent. In that sense
it may be suitable to augment the ``SmartStar`` particle with your required feature rather than implmenting a new feature. You
will also be able to build on ``SmartStar`` tests and documentation too rather than starting from scratch. 

To allow for the creation of a ``SmartStar`` particle in an Enzo simulation the following line must be included in the parameter file:
``AppendActiveParticleType = SmartStar``

Once Enzo reads that line in the parameter file then the AP framework will be engaged and the ``SmartStar`` particle initialised.
The following parameters are currently enabled for the ``SmartStar``

The parameters which drive the ``SmartStar`` particle type can be found at  :ref:`active_particles_parameters`.

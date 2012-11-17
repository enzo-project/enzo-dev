Conduction
~~~~~~~~~~

Isotropic and anisotropic thermal conduction are implemented using the
method of Parrish and Stone: namely, using an explicit, forward
time-centered algorithm.  In the anisotropic conduction, heat can only
conduct along magnetic field lines.  One can turn on the two types of
conduction independently, since there are situations where one might 
want to use both.  The Spitzer fraction can be also set
independently for the isotropic and anisotropic conduction.

``IsotropicConduction`` (external)
    Turns on isotropic thermal conduction using Spitzer conduction.  Default: 0 (FALSE)
``AnisotropicConduction`` (external)
    Turns on anisotropic thermal conduction using Spitzer conduction.
    Can only be used if MHD is turned on (``HydroMethod`` = 4).
    Default: 0 (FALSE)
``IsotropicConductionSpitzerFraction`` (external)
    Prefactor that goes in front of the isotropic Spitzer conduction
    coefficient.  Should be a value between 0 and 1.
    Default: 1.0
``AnisotropicConductionSpitzerFraction`` (external)
    Prefactor that goes in front of the anisotropic Spitzer conduction
    coefficient.  Should be a value between 0 and 1.
    Default: 1.0
``ConductionCourantSafetyNumber`` (external)
    This is a prefactor that controls the stability of the conduction
    algorithm.  In its current explicit formulation, it must be set to
    a value of 0.5 or less.
    Default: 0.5


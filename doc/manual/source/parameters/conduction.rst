Conduction
~~~~~~~~~~

Isotropic and anisotropic thermal conduction are implemented using the
method of Parrish and Stone: namely, using an explicit, forward
time-centered algorithm.  In the anisotropic conduction, heat can only
conduct along magnetic field lines.  One can turn on the two types of
conduction independently, since there are situations where one might 
want to use both.  The Spitzer fraction can be also set
independently for the isotropic and anisotropic conduction.  Running a 
cosmological simulation with conduction on can be tricky as the timesteps 
can become very short.  It is recommended that you look carefully at all the 
available conduction parameters.  Additionally, if you intend to run with 
star particles, it is highly recommended that you set the parameter, 
``StarMakerTimeIndependentFormation``.  See the description in 
:ref:`starparticleparameters` for more information.

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
``SpeedOfLightTimeStepLimit`` (external)
    When used, this sets a floor for the conduction timestep to be the local light crossing time (dx / c).  This prevents the conduction machinery from prescribing extremely small timesteps.  While this can technically violate the conduction stability criterion, testing has shown that this does not result in notable differences.  (1 - ON; 0 - OFF)  Default: 0 (OFF).
``ConductionDynamicRebuildHierarchy`` (external)
    Using conduction can often result in the code taking extremely short timesteps.  Since the hierarchy is rebuilt each timestep, this can exacerbate memory fragmentation issues and slow the simulation.  In the case where the conduction timestep is the limiter, the hierarchy should not need to be rebuilt every timestep since conduction mostly does not alter the fields which control refinement.  When this option is used, the timestep calculation is carried out as usual, but the hierarchy is only rebuilt on a timescale that is calculated neglecting the conduction timestep.  This results in a decent speedup and reduced memory fragmentation when running with conduction.  (1 - ON; 0 - OFF)  Default: 0 (OFF).
``ConductionDynamicRebuildMinLevel`` (external)
    The minimum level on which the dynamic hierarcy rebuild is performed.  Default: 0.

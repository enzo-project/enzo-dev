Writing Enzo Parameter Files
============================

Putting together a parameter file for Enzo is possibly the most
critical step when setting up a simulation, and is certainly the step
which is most fraught with peril. There are over 200 parameters that
one can set - see :ref:`parameters` for a complete
listing. For the most part, defaults are set to be sane values for
cosmological simulations, and most physics packages are turned off by
default, so that you have to explicitly turn on modules. All physics
packages are compiled into Enzo (unlike codes such as ZEUS-MP 1.0,
where you have to recompile the code in order to enable new physics).

It is inadvisable for a novice to put together a parameter file from
scratch.  There are many example parameter files in the ``enzo/run``
directory, classified by the kind of test file.  See also the section
:ref:`creating cosmological initial
conditions. <CosmologicalInitialConditions>`.

In order to make the most of this tutorial it is advisable to have
one or more of these parameter files open while reading this page.
For the purposes of this tutorial we assume that the user is
putting together a cosmology simulation and has already generated
the initial conditions files using inits.

All parameters are put into a plain text file (one parameter per
line), the name of which is fed into Enzo at execution time at the
command line. Typically, a parameter is set by writing the
parameter name, an equals sign, and then the parameter value or
values, like this:

.. highlight:: none

::

    NumberOfBufferZones = 3

You must leave at least one space between the parameter, the equals
sign, and the parameter value. It's fine if you use more than one
space - after the first space, whitespace is unimportant. All lines
which start with a # (pound sign) are treated as comments and
ignored. In addition, you can have inline comments by using the
same pound sign, or two forward slashes // after the parameter line.

::

    NumberOfBufferZones = 3 // More may be needed depending on physics used.

Initialization parameters
-------------------------

Complete descriptions of all initialization parameters are given
here. The most fundamental initialization parameter you have to set
is ``ProblemType``, which specifies the type of problem to be run, and
therefore the way that Enzo initiates the data. A cosmology
simulation is problem type 30. As started before, for the purposes
of this introduction I'm assuming that you are generating a
cosmology simulation, so you would put this line in the parameter
file:

::

    ProblemType = 30

``TopGridRank`` specifies the spatial dimensionality of your problem
(1, 2 or 3 dimensions), and must be set. ``TopGridDimensions``
specifies the number of root grid cells along each axis. For a 3D
simulation with 128 grid cells along each axis on the root grid,
put this in the parameter file:

::

    TopGridRank = 3
    TopGridDimensions = 128 128 128

Additionally, you must specify the names of the initial conditions
files with contain the baryon density and velocity information and
the dark matter particle positions and velocities. These are
controlled via the parameters ``CosmologySimulationDensityName``,
``CosmologySimulationVelocity[123]Name`` (where 1, 2 and 3 correspond
to the x, y and z directions, respectively),
``CosmologySimulationParticlePositionName`` and
``CosmologySimulationParticleVelocityName``. Assuming that the baryon
velocity information is all in a single file, and that the baryon
density and velocity file names are ``GridDensity`` and ``GridVelocities``,
and that the particle position and velocity files are named
``ParticlePositions`` and ``ParticleVelocities``, these parameters would be
set as follows:

::

    CosmologySimulationDensityName = GridDensity
    CosmologySimulationVelocity1Name = GridVelocities
    CosmologySimulationVelocity2Name = GridVelocities
    CosmologySimulationVelocity3Name = GridVelocities
    CosmologySimulationParticlePositionName = ParticlePositions
    CosmologySimulationParticleVelocityName = ParticleVelocities

Some more advanced are parameters in the Initialization Parameters
section control domain and boundary value specifications. These
should NOT be altered unless you really, really know what you're
doing!

Cosmology
---------

See also the section on :ref:`creating cosmological initial
conditions. <CosmologicalInitialConditions>`.  

Complete descriptions of all cosmology parameters are given 
:ref:`here <cosmology-parameters>` and 
:ref:`here <cosmologysimulation_param>`. ``ComovingCoordinates`` determines 
whether comoving
coordinates are used or not. In practice, turning this off turns
off all of the cosmology machinery, so you want to leave it set to
1 for a cosmology simulation. ``CosmologyInitialRedshift`` and
``CosmologyFinalRedshift`` control the start and end times of the
simulation, respectively. ``CosmologyHubbleConstantNow`` sets the
Hubble parameter, and is specified at z=0 in units of 100 km/s/Mpc.
``CosmologyComovingBoxSize`` sets the size of the box to be simulated
(in units of Mpc/h) at z=0. ``CosmologySimulationOmegaBaryonNow``,
``CosmologySimulationOmegaCDMNow``, ``CosmologyOmegaMatterNow``, and
``CosmologyOmegaLambdaNow`` set the amounts of baryons, dark matter,
total matter, and vacuum energy (in units of the critical density at
z=0). Setting ``CosmologySimulationUseMetallicityField`` to 1 will 
create an additional tracer field for following metals. This is handy for
simulations with star formation and feedback (described below). For
example, in a cosmology simulation with box size 100 Mpc/h with
approximately the cosmological parameters determined by WMAP, which
starts at z=50 and ends at z=2, and has a metal tracer field, we
put the following into the parameter file:

::

    ComovingCoordinates = 1
    CosmologyInitialRedshift = 50.0
    CosmologyFinalRedshift = 2.0
    CosmologyHubbleConstantNow = 0.7
    CosmologyComovingBoxSize = 100.0
    CosmologyOmegaMatterNow = 0.3
    CosmologyOmegaLambdaNow = 0.7
    CosmologySimulationOmegaBaryonNow = 0.04
    CosmologySimulationOmegaCDMNow = 0.26
    CosmologySimulationUseMetallicityField = 1

Gravity and Particle Parameters
-------------------------------

The parameter list sections on gravity particle positions are here
and here, respectively. The significant gravity-related parameters
are ``SelfGravity``, which turns gravity on (1) or off (0) and
``GravitationalConstant``, which must be 1 in cosmological
simulations. ``BaryonSelfGravityApproximation`` controls whether
gravity for baryons is determined by a quick and reasonable
approximation. It should be left on (1) in most cases. For a
cosmological simulation with self gravity, we would put the
following parameters into the startup file:

::

    SelfGravity = 1
    GravitationalConstant = 1
    BaryonSelfGravityApproximation = 1

We discuss some AMR and parallelization-related particle parameters
in later sections.

Adiabatic hydrodynamics parameters
----------------------------------

The parameter listing section on hydro parameters can be found
here. The most fundamental hydro parameter that you can set is
``HydroMethod``, which lets you decide between the Piecewise Parabolic
Method (aka PPM; option 0), or the finite-difference method used in
the Zeus astrophysics code (option 2). PPM is the more advanced and
optimized method. The Zeus method uses an artificial viscosity-based
scheme and may not be suited for some types of work. When using PPM in
a cosmological simulation, it is important to turn
``DualEnergyFormalism`` on (1), which makes total-energy schemes such
as PPM stable in a regime where there are hypersonic fluid flows,
which is quite common in cosmology. The final parameter that one must
set is ``Gamma``, the ratio of specific heats for an ideal gas. If
``MultiSpecies`` (discussed later in :ref:`RadCooling`) is on, this is
ignored. For a cosmological simulation where we wish to use PPM and
have ``Gamma`` = 5/3, we use the following parameters:

::

    HydroMethod = 0
    DualEnergyFormalism = 1
    Gamma = 1.66667

In addition to these three parameters, there are several others
which control more subtle aspects of the two hydro methods. See the
parameter file listing of hydro parameters for more information on
these.

One final note: If you are interested in performing simulations
where the gas has an isothermal equation of state (gamma = 1), this
can be approximated without crashing the code by setting the
parameter Gamma equal to a number which is reasonably close to one,
such as 1.001.

AMR Hierarchy Control Parameters
--------------------------------

These parameters can be found in the parameter list page here. They
control whether or not the simulation uses adaptive mesh
refinement, and if so, the characteristics of the adaptive meshing
grid creation and refinement criteria. We'll concentrate on a
simulation with only a single initial grid first, and then discuss
multiple levels of initial grids in a subsection.

The most fundamental AMR parameter is ``StaticHierarchy``. When this is
on (1), the code is a unigrid code. When it is off (0), adaptive
mesh is turned on. ``RefineBy`` controls the refinement factor - for
example, a value of 2 means that a child grid is twice as highly
refined as its parent grid. It is important to set ``RefineBy`` to 2
when using cosmology simulations - this is because if you set it to
a larger number (say 4), the ratio of particle mass to gas mass in
a cell grows by a factor of eight during each refinement, causing
extremely unphysical effects.
``MaximumRefinementLevel`` determines how many possible levels of
refinement a given simulation can attain, and
``MaximumGravityRefinementLevel`` defines the maximum level at which
gravitational accelerations are computed. More highly refined
levels have their gravitational accelerations interpolated from
this level, which effectively provides smoothing of the
gravitational force on the spatial resolution of the grids at
``MaximumGravityRefinementLevel``. A simulation with AMR turned on,
where there are 6 levels of refinement (with gravity being smoothed
on level 4) and where each child grid is twice as highly resolved
as its parent grid would have these parameters set as follows:

::

    StaticHierarchy = 0
    RefineBy = 2
    MaximumRefinementLevel = 6
    MaximumGravityRefinementLevel = 4

Once the AMR is turned on, you must specify how and where the
hierarchy
refines. The parameter ``CellFlaggingMethod`` controls the method in
which cells are flagged, and can be set with multiple values. We
find that refining by baryon and dark matter mass (options 2 and 4)
are typically useful in cosmological simulations. The parameter
``MinimumOverDensityForRefinement`` allows you to control the
overdensity at which a given grid is refined, and can is set with
multiple values as well. Another very useful parameter is
``MinimumMassForRefinementLevelExponent``, which modifies the cell
masses/overdensities used for refining grid cells. See the
parameter page for a more detailed explanation. 
Leaving this with a value of 0.0 ensures that gas mass resolution
in dense regions remains more-or-less Lagrangian in nature.
Negative values make the refinement super-Lagrangian (ie, each
level has less gas mass per cell on average than the coarser level
above it) and positive values make the refinement sub-lagrangian.
In an AMR simulation where the AMR triggers on baryon and dark
matter overdensities in a given cell of 4.0 and 8.0, respectively,
where the refinement is slightly super-Lagrangian, these paramaters
would be set as follows:

::

    CellFlaggingMethod = 2 4
    MinimumOverDensityForRefinement = 4.0 8.0
    MinimumMassForRefinementLevelExponent = -0.1

At times it is very useful to constrain your simulation such that
only a small region is adaptively refined (the default is to refine
over an entire simulation volume). For example, if you wish to
study the formation of a particular galaxy in a very large volume,
you may wish to only refine in the small region around where that
galaxy forms in your simulation in order to save on computational
expense and dataset size. Two parameters, ``RefineRegionLeftEdge`` and
``RefineRegionRightEdge`` allow control of this. For example, if we
only want to refine in the inner half of the volume (0.25 - 0.75
along each axis), we would set these parameters as follows:

::

    RefineRegionLeftEdge = 0.25 0.25 0.25
    RefineRegionRightEdge = 0.75 0.75 0.75

This pair of parameters can be combined with the use of nested
initial grids (discussed in the next subsection) to get simulations
with extremely high dark matter mass and spatial resolution in a
small volume at reasonable computational cost.

Multiple nested grids
~~~~~~~~~~~~~~~~~~~~~

At times it is highly advantageous to use multiple nested grids.
This is extremely useful in a situation where you are interested in
a relatively small region of space where you need very good dark
matter mass resolution and spatial resolution while at the same
time still resolving large scale structure in order to preserve
gravitational tidal forces. An excellent example of this is
formation of the first generation of objects in the universe, where
we are interested in a relatively small (10\ :sup:`6`\  solar mass)
halo which is strongly tidally influenced by the large-scale
structure around it. It is important to resolve this halo with a
large number of dark matter particles in order to reduce frictional
heating, but the substructure of the distant large-scale structure
is not necessarily interesting, so it can be resolved by very
massive particles. One could avoid the complication of multiple
grids by using a single very large grid - however, this would be
far more computationally expensive.

Let us assume for the purpose of this example that in addition to
the initial root grid grids (having 128 grid cells along each axis)
there are two subgrids, each of which is half the size of the one
above it in each spatial direction (so subgrid 1 spans from
0.25-0.75 in units of the box size and subgrid 2 goes from
0.375-0.625 in each direction). If each grid is twice as highly
refined spatially as the one above it, the dark matter particles on
that level are 8 times smaller, so the dark matter mass resolution
on grid #2 is 64 times better than on the root grid, while the
total number of initial grid cells only increases by a factor of
three (since each grid is half the size, but twice as highly
refined as the one above it, the total number of grid cells remains
the same). Note: See the page on generating initial conditions for
more information on creating this sort of set of nested grids.

When a simulation with more than one initial grid is run, the total
number of initial grids is specified by setting
``CosmologySimulationNumberOfInitialGrids``. The parameter
``CosmologySimulationGridDimension[#]`` is an array of three integers
setting the grid dimensions of each nested grid, and
``CosmologySimulationGridLeftEdge[#]`` and
``CosmologySimulationGridRightEdge[#]`` specify the left and right
edges of the grid spatially, in units of the box size. In the last
three parameters, "#" is replaced with the grid number. The root
grid is grid 0. None of the previous three parameters need to be
set for the root grid. For the setup described above, the parameter
file would be set as follows:

::

    CosmologySimulationNumberOfInitialGrids = 3
    CosmologySimulationGridDimension[1] = 128 128 128
    CosmologySimulationGridLeftEdge[1] = 0.25 0.25 0.25
    CosmologySimulationGridRightEdge[1] = 0.75 0.75 0.75
    CosmologySimulationGridLevel[1] = 1
    CosmologySimulationGridDimension[2] = 128 128 128
    CosmologySimulationGridLeftEdge[2] = 0.375 0.375 0.375
    CosmologySimulationGridRightEdge[2] = 0.625 0.625 0.625
    CosmologySimulationGridLevel[2] = 2

Multiple initial grids can be used with or without AMR being turned
on. If AMR is used, the parameter ``MinimumOverDensityForRefinement``
must be modified as well. It is advisable to carefully read the
entry for this parameter in the parameter list (in this section).
The minimum overdensity
needs to be divided by r\ :sup:`(d\*l)`\ , where r is the refinement
factor, d is the dimensionality, and l is the zero-based highest
level of the initial grids. So if we wish for the same values for
``MinimumOverDensityForRefinement`` used previous to apply on the most
highly refined grid, we must divide the set values by
2\ :sup:`(3\*2)`\  = 64. In addition, one should only refine on the
highest level, so we must reset ``RefineRegionLeftEdge`` and
``RefineRegionRightEdge``. The parameters would be reset as follows:

::

    RefineRegionLeftEdge = 0.375 0.375 0.375
    RefineRegionRightEdge = 0.625 0.625 0.625
    MinimumOverDensityForRefinement = 0.0625 0.125

A note: When creating multi-level intial conditions, make sure that
the initial conditions files for all levels have the same file name
(ie, ``GridDensity``), but that each file has an extension which is an
integer corresponding to its level. For example, the root grid
``GridDensity`` file would be ``GridDensity.0``, the level 1 file would be
``GridDensity.1``, and so forth. The parameters which describe file
names (discussed above in the section on initialization parameters)
should only have the file name to the left of the period the period
(as in a simulation with a single initial grid), ie,

::

    CosmologySimulationDensityName = GridDensity

Nested Grids and Particles
~~~~~~~~~~~~~~~~~~~~~~~~~~

When initializing a nested grid problem, there can arise an issue of
lost particles as a result of running ring. Please see
:doc:`../reference/NestedGridParticles` for more information.

I/O Parameters
--------------

These parameters, defined in more detail in
:doc:`ControllingDataOutput`, control all aspects of Enzo's data
output. One can output data in a cosmological simulation in both a
time-based and redshift-based manner. To output data regularly in
time, one sets ``dtDataDump`` to a value greater than zero. The size
of this number, which is in units of Enzo's internal time variable,
controls the output frequency.  See the Enzo user's manual section on
output format for more information on physical units. Data can be
output at specific redshifts as controlled by
``CosmologyOutputRedshift[#]``, where # is the number of the output
dump (with a maximum of 10,000 zero-based numbers). The name of the
time-based output files are controlled by the parameter
``DataDumpName`` and the redshift-based output files have filenames
controlled by ``RedshiftDumpName``. For example, if we want to output
data every time the code advances by dt=2.0 (in code units) with file
hierarchiess named ``time_0000``, ``time_0001``, etc., and ALSO output
explicitly at redshifts 10, 5, 3 and 1 with file hierarchy names
``RedshiftOutput0000``, ``RedshiftOutput0001``, etc., we would set
these parameters as follows:

::

    dtDataDump = 2.0
    DataDumpName = time_
    RedshiftDumpName = RedshiftOutput
    CosmologyOutputRedshift[0] = 10.0
    CosmologyOutputRedshift[1] = 5.0
    CosmologyOutputRedshift[2] = 3.0
    CosmologyOutputRedshift[3] = 1.0

Note that Enzo always outputs outputs data at the end of the
simulation, regardless of the settings of ``dtDataDump`` and
``CosmologyOutputRedshift``.

.. _RadCooling:

Radiative Cooling and UV Physics Parameters
-------------------------------------------

Enzo comes with multiple ways to calculate baryon cooling and a
metagalactic UV background, as described in detail here. The
parameter ``RadiativeCooling`` controls whether or not a radiative
cooling module is called for each grid. The cooling is calculated
either by assuming equilibrium cooling and reading in a cooling
curve, or by computing the cooling directly from the species
abundances. The parameter ``MultiSpecies`` controls which cooling
module is called - if ``MultiSpecies`` is off (0) the equilibrium model
is assumed, and if it is on (1 or 2) then nonequilibrium cooling is
calculated using either 6 or 9 ionization states of hydrogen and
helium (corresponding to ``MultiSpecies`` = 1 or 2, respectively). The
UV background is controlled using the parameter ``RadiationFieldType``.
Currently there are roughly a dozen backgrounds to choose from.
``RadiationFieldType`` is turned off by default, and can only be used
when ``Multispecies`` = 1. For example, if we wish to use a
nonequilibrium cooling model with a Haardt and Madau background
with q\ :sub:`alpha`\ = -1.8, we would set these parameters as follows:

::

    RadiativeCooling = 1
    MultiSpecies = 1
    RadiationFieldType = 2

Star Formation and Feedback Physics Parameters
----------------------------------------------

Enzo has multiple routines for star formation and feedback.  Star
particle formation and feedback are controlled separately, by the
parameters ``StarParticleCreation`` and ``StarParticleFeedback``.
Multiple types of star formation and feedback can be used, e.g. models
for Pop III stars for metal-free gas and models for Pop II stars for
metal-enriched gas.  These routines are disabled when these parameters
are set equal to 0.  These parameters are bitwise to allow multiple
types of star formation routines can be used in a single
simulation. For example if methods 1 and 3 are desired, the user would
specify 10 (2\ :sup:`1`\ + 2\ :sup:`3`\ ), or if methods 0, 1 and 4
are wanted, this would be 19 (2\ :sup:`0`\ + 2\ :sup:`1`\ + 2\
:sup:`4`\ ).  See :ref:`StarParticleParameters` for more details.

They are turned on when the i-th bit is flagged.  The value of 2 is
the recommended value. The most commonly used routines (2) are based
upon an algorithm by Cen & Ostriker, and there are a number of free
parameters. Note that it is possible to turn star particle formation
on while leaving feedback off, but not the other way around.

For the star particle creation algorithm, stars are allowed to form
only in cells where a minimum overdensity is reached, as defined by
``StarMakerOverDensityThreshold``. Additionally, gas can only turn into
stars with an efficiency controlled by ``StarMakerMassEfficiency`` and
at a rate limited by ``StarMakerMinimumDynamicalTime``, and the minimum
mass of any given particle is controlled by the parameter
``StarMakerMinimumStarMass``, which serves to limit the number of star
particles. For example, if we wish to use the "standard" star
formation scenario where stars can only form in cells which are at
least 100 times the mean density, with a minimum dynamical time of
10\ :sup:`6`\  years and a minimum mass of 10\ :sup:`7`\  solar
masses, and where only 10% of the baryon gas in a cell can be
converted into stars in any given timestep, we would set these
parameters as follows:

::

    StarParticleCreation = 2
    StarMakerOverDensityThreshold = 100.0
    StarMakerMassEfficiency = 0.1
    StarMakerMinimumDynamicalTime = 1.0e6
    StarMakerMinimumStarMass = 1.0e7

Star particles can provide feedback into the Inter-Galactic Medium via stellar winds,
thermal energy and metal pollution. The parameter
``StarMassEjectionFraction`` controls the fraction of the total initial
mass of the star particle which is eventually returned to the gas
phase. ``StarMetalYield`` controls the mass fraction of metals produced
by each star particle that forms, and ``StarEnergyToThermalFeedback``
controls the fraction of the rest-mass energy of the stars created
which is returned to the gas phase as thermal energy. Note that the
latter two parameters are somewhat constrained by theory and
observation to be somewhere around 0.02 and 1.0e-5, respectively.
The ejection fraction is poorly constrained as of right now. Also,
metal feedback only takes place if the metallicity field is turned
on (``CosmologySimulationUseMetallicityField`` = 1). As an example, if
we wish to use the 'standard' star feedback where 25% of the total
stellar mass is returned to the gas phase, the yield is 0.02 and
10\ :sup:`-5`\  of the rest mass is returned as thermal energy, we
set our parameters as follows:

::

    StarParticleFeedback = 2
    StarMassEjectionFraction = 0.25
    StarMetalYield = 0.02
    StarEnergyToThermalFeedback = 1.0e-5
    CosmologySimulationUseMetallicityField = 1

When using the star formation and feedback algorithms it is
important to consider the regime of validity of our assumptions.
Each "star particle" is supposed to represent an ensemble of stars,
which we can characterize with the free parameters described above.
This purely phenomenological model is only reasonable as long as
the typical mass of the star particles is much greater than the
mass of the heaviest stars so that the assumption of averaging over
a large population is valid. When the typical star particle mass
drops to the point where it is comparable to the mass of a large
star, these assumptions must be reexamined and our algorithms
reformulated.

IO Parallelization Options
--------------------------

One of Enzo's great strengths is that it is possible to do
extremely large simulations on distributed memory machines. For
example, it is possible to intialize a 1024\ :sup:`3`\  root grid
simulation on a linux cluster where any individual node has 1 or 2
GB of memory, which is on the order of 200 times less than the
total dataset size! This is possible because the reading of initial
conditions and writing out of data dumps is fully parallelized - at
startup, when the parameter ``ParallelRootGridIO`` is turned on each
processor only reads the portion of the root grid which is within
its computational domain, and when ``ParallelParticleIO`` is turned on
each processor only reads in the particles within its domain
(though preprocessing is needed - see below). Additionally, the
parameter ``Unigrid`` should be turned on for simulations without AMR,
as it saves roughly a factor of two in memory on startup, allowing
the code to perform even larger simulations for a given computer
size. If we wish to perform an extremely large unigrid simulation
with parallel root grid and particle IO, we would set the following
parameters:

::

    ParallelParticleIO = 1
    ParallelRootGridIO = 1
    Unigrid = 1

AMR simulations can be run with ``ParallelRootGridIO`` and
``ParallelParticleIO`` on, though you must be careful to turn off the
``Unigrid`` parameter. In addition, it is important to note that in the
current version of Enzo you must run the program called "ring" on
the particle position and velocity files before Enzo is started in
order to take advantage of the parallel particle IO. Assuming the
particle position and velocity files are named ``ParticlePositions``
and ``ParticleVelocities``, respectively, this is done by running:

::

    mpirun -np [N] ring ParticlePositions ParticleVelocities

Where mpirun is the executable responsible for running MPI programs
and "-np [N]" tells the machine that there are [N] processors. This
number of processors must be the same as the number which Enzo will
be run with!

Notes
-----

This page is intended to help novice Enzo users put together parameter
files for their first simulation and therefore is not intended to be
an exhaustive list of parameters nor a complete description of each
parameter mentioned. It would be wise to refer to the Enzo user
guide's :ref:`parameters` for a more-or-less complete list of
AMR parameters, some of which may be extremely useful for your
specific application.



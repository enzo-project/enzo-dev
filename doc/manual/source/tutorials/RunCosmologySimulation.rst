.. _RunCosmologySimulation:

How to run a cosmology simulation
=================================

In order to run a cosmology simulation, you'll need to build enzo.exe,
and an initial conditions generator.  We recommend using the `MUSIC
<https://bitbucket.org/ohahn/music>`.

inits.exe and ring.exe (see :doc:`building_enzo`) inits creates the
initial conditions for your simulation, and ring splits up the root
grid which is necessary if you're using parallel IO. Once you have
built the three executables, put them in a common directory where you
will run your test simulation. You will also save the inits and param
files (shown and discussed below) in this directory.

Creating initial conditions
---------------------------

The first step in preparing the simulation is to create the initial
conditions. The file inits uses is a text file which contains a
list of parameters with their associated values. These
values tell the initial conditions generator necessary information
like the simulation box size, the cosmological parameters and the
size of the root grid. The code then takes that information and
creates a set of initial conditions. Here is an example inits
file:

.. highlight:: none

::

    #
    #  Generates initial grid and particle fields for a 
    #    CDM simulation
    #
    #  Cosmology Parameters
    #
    CosmologyOmegaBaryonNow      = 0.044
    CosmologyOmegaMatterNow      = 0.27 
    CosmologyOmegaLambdaNow      = 0.73  
    CosmologyComovingBoxSize     = 10.0    // in Mpc/h
    CosmologyHubbleConstantNow   = 0.71      // in units of 100 km/s/Mpc
    CosmologyInitialRedshift     = 60
    #
    #  Power spectrum Parameters
    #
    
    PowerSpectrumType            = 11
    PowerSpectrumSigma8          = 0.9
    PowerSpectrumPrimordialIndex = 1.0
    PowerSpectrumRandomSeed      = -584783758
    #
    #  Grid info
    #
    Rank                = 3
    GridDims            = 32 32 32
    InitializeGrids     = 1
    GridRefinement      = 1
    #
    #  Particle info
    #
    ParticleDims        = 32 32 32
    InitializeParticles = 1
    ParticleRefinement  = 1
    #
    #  Overall field parameters
    #
    #
    #  Names
    #
    ParticlePositionName = ParticlePositions
    ParticleVelocityName = ParticleVelocities
    GridDensityName      = GridDensity
    GridVelocityName     = GridVelocities

inits is run by typing this command:

::

    ./inits.exe -d Example_Cosmology_Sim.inits

inits will produce some output to the screen to tell you what it is
doing, and will write five files: ``GridDensity``, ``GridVelocities``,
``ParticlePositions``, ``ParticleVelocities`` and ``PowerSpectrum.out``. The
first four files contain information on initial conditions for the
baryon and dark matter componenets of the simulation, and are HDF5
files. The last file is an ascii file which contains information on
the power spectrum used to generate the initial conditions.

It is also possible to run cosmology simulations using initial
nested subgrids.

Parallel IO - the ring tool
---------------------------

This simulation is quite small. The root grid is only 32 cells on a
side and we allow a maximum of three levels of mesh refinement.
Still, we will use the ring tool, since it is important for larger
simulations of sizes typically used for doing science.  Additionally,
if you wish to run with 64 or more processors, you should use
``ParallelRootGridIO``, described in :ref:`ParallelRootGridIO`.

The ring tool is part of the Enzo parallel IO (input-output)
scheme. Examine the last section of the parameter file (see below)
for this example simulation and you will see:

::

    #
    # IO parameters
    #
    ParallelRootGridIO = 1
    ParallelParticleIO = 1

These two parameters turn on parallel IO for both grids and
particles. In a serial IO simulation where multiple processors are
being used, the master processor reads in all of the grid and
particle initial condition information and parcels out portions of
the data to the other processors. Similarly, all simulation output
goes through the master processor as well. This is fine for
relatively small simulations using only a few processors, but slows
down the code considerably when a huge simulation is being run on
hundreds of processors. Turning on the parallel IO options allows
each processor to perform its own IO, which greatly decreases the
amount of time the code spends performing IO.

The process for parallelizing grid and particle information is quite different.
Since it is known exactly where every grid cell in a structured Eulerian grid
is in space, and these cells are stored in a regular and predictable order in
the initial conditions files, turning on ``ParallelRootGridIO`` simply tells
each processor to figure out which portions of the arrays in the GridDensity
and ``GridVelocities`` belong to it, and then read in only that part of the
file. The particle files (``ParticlePositions`` and ``ParticleVelocities``)
store the particle information in no particular order.  In order to efficiently
parallelize the particle IO the ring tool is used.  ring is run on the same
number of processors as the simulation that you intend to run, and is typically
run just before Enzo is called for this reason.  In ring, each processor reads
in an equal fraction of the particle position and velocity information into a
list, flags the particles that belong in its simulation spatial domain, and
then passes its portion of the total list on to another processor. After each
portion of the list has made its way to every processor, each processor then
collects all of the particle and velocity information that belongs to it and
writes them out into files called ``PPos.nnnn`` and ``PVel.nnnn``, where nnnn
is the processor number. Turning on the ``ParallelParticleIO`` flag in the Enzo
parameter file instructs Enzo to look for these files.

For the purpose of this example, you're going to run ring and Enzo on 4
processors (this is a fixed requirement).  The number of processors used in an
MPI job is set differently on each machine, so you'll have to figure out how
that works for you. On some machines, you can request an 'interactive queue' to
run small MPI jobs. On others, you may have to submit a job to the batch queue,
and wait for it to run.

To start an interactive run, it might look something like this:

::

    qsub -I -V -l walltime=00:30:00,size=4

This tells the queuing system that you want four processors total for a
half hour of wall clock time. You may have to wait a bit until
nodes become available, and then you will probably start out back
in your home directory. You then run ring on the particle files by
typing something like this:

::

    mpirun -n 4 ./ring.exe pv ParticlePositions ParticleVelocities

This will then produce some output to your screen, and will
generate 8 files: ``PPos.0000`` through ``PPos.0003`` and ``PVel.0000`` through
``PVel.0003``. Note that the 'mpirun' command may actually be 'aprun'
or something similar. Consult your supercomputer's documentation to
figure out what this command should really be.

Congratulations, you're now ready to run your cosmology
simulation!

Running an Enzo cosmology simulation
------------------------------------

After all of this preparation, running the simulation itself should
be straightforward. First, you need to have an Enzo parameter file.
Here is an example compatible with the inits file above:

::

    #
    # AMR PROBLEM DEFINITION FILE: Cosmology Simulation (AMR version)
    #
    #  define problem
    #
    ProblemType                = 30      // cosmology simulation
    TopGridRank                = 3
    TopGridDimensions          = 32 32 32
    SelfGravity                = 1       // gravity on
    TopGridGravityBoundary     = 0       // Periodic BC for gravity
    LeftFaceBoundaryCondition  = 3 3 3   // same for fluid
    RightFaceBoundaryCondition = 3 3 3
    #
    #  problem parameters
    #
    CosmologySimulationOmegaBaryonNow       = 0.044
    CosmologySimulationOmegaCDMNow      = 0.226 
    CosmologyOmegaMatterNow         = 0.27 
    CosmologyOmegaLambdaNow         = 0.73  
    CosmologySimulationDensityName          = GridDensity
    CosmologySimulationVelocity1Name        = GridVelocities
    CosmologySimulationVelocity2Name        = GridVelocities
    CosmologySimulationVelocity3Name        = GridVelocities
    CosmologySimulationParticlePositionName = ParticlePositions
    CosmologySimulationParticleVelocityName = ParticleVelocities
    CosmologySimulationNumberOfInitialGrids = 1
    #
    #  define cosmology parameters
    #
    ComovingCoordinates        = 1       // Expansion ON
    CosmologyHubbleConstantNow = 0.71    // in km/s/Mpc
    CosmologyComovingBoxSize   = 10.0  // in Mpc/h
    CosmologyMaxExpansionRate  = 0.015   // maximum allowed delta(a)/a
    CosmologyInitialRedshift   = 60.0      // 
    CosmologyFinalRedshift     = 3.0     //
    GravitationalConstant      = 1       // this must be true for cosmology
    #
    #  set I/O and stop/start parameters
    #
    CosmologyOutputRedshift[0] = 25.0 
    CosmologyOutputRedshift[1] = 10.0
    CosmologyOutputRedshift[2] = 5.0  
    CosmologyOutputRedshift[3] = 3.0
    #
    #  set hydro parameters
    #
    Gamma                  = 1.6667
    PPMDiffusionParameter  = 0       // diffusion off
    DualEnergyFormalism    = 1       // use total & internal energy
    InterpolationMethod    = 1     // SecondOrderA
    CourantSafetyNumber    = 0.5
    ParticleCourantSafetyNumber = 0.8
    FluxCorrection         = 1
    ConservativeInterpolation = 0
    HydroMethod            = 0
    #
    #  set cooling parameters
    #
    RadiativeCooling       = 0
    MultiSpecies           = 0
    RadiationFieldType     = 0
    StarParticleCreation   = 0
    StarParticleFeedback   = 0
    #
    #  set grid refinement parameters
    #
    StaticHierarchy           = 0    // AMR turned on!
    MaximumRefinementLevel    = 3
    MaximumGravityRefinementLevel = 3
    RefineBy                  = 2
    CellFlaggingMethod        = 2 4
    MinimumEfficiency         = 0.35
    MinimumOverDensityForRefinement = 4.0 4.0
    MinimumMassForRefinementLevelExponent = -0.1
    MinimumEnergyRatioForRefinement = 0.4 
    
    #
    #  set some global parameters
    #
    GreensFunctionMaxNumber   = 100   // # of greens function at any one time
    
    
    #
    # IO parameters
    #
    
    ParallelRootGridIO = 1
    ParallelParticleIO = 1

Once you've saved this, you start Enzo by typing:

::

    mpirun -n 4 ./enzo.exe -d Example_Cosmology_Sim.param >& output.log

The simulation will now run. The -d flag ensures a great deal of
output, so you may redirect it into a log file called ``output.log``
for later examination. This particular simulation shouldn't take
too long, so you can run this in the same 30 minute interactive job
you started when you ran inits. When the simulation is done, Enzo
will display the message "Successful run, exiting."

Congratulations! If you've made it this far, you have now successfully
run a cosmology simulation using Enzo! 

Cosmology Step by Step
======================

This is a step-by-step guide to getting Enzo compiled and a test problem run on
your computer of choice. We assume no previous experience with the code. The
only assumption is that the computer that you are installing Enzo on meets the
:ref:`CompilationRequirements`, (i.e., MPI and HDF5 are installed), and that
you have GNU Make installed.

This page shows actual inputs and outputs from an Enzo installation performed
on `SDSC's DataStar <http://www.sdsc.edu/us/resources/datastar/>`_ on one of
its final days. Though we show inputs from this machine, they should be
generally applicable to your machine of choice, as long as it runs a `POSIX
<http://en.wikipedia.org/wiki/POSIX>`_ operating system. We have successfully
compiled and run Enzo on Sun, SGI, IBM, Compaq, Hewlett-Packard and Cray
machines, as well as Apple computers running OS X and a wide variety of
machines running the Linux operating system.

Obtaining the Enzo package
--------------------------

Once we're sure that we've met the
`compilation requirements? </wiki/UserGuide/CompilationRequirements>`_,
we can download Enzo and compile it.

Download
~~~~~~~~

::

    ds100 $ wget http://example.edu/TBD.tar

Subversion Checkout
~~~~~~~~~~~~~~~~~~~

::

    ds100 $ svn co http://mngrid.ucsd.edu/svn/Enzo/public/trunk enzo
    ...
    Fetching external item into 'src/yt'
    External at revision 742.
    
    At revision 1788.
    ds100 $

Compiling Enzo
--------------

Before we can run any tests, we need an executable. For detailed
instructions on how to do this, head over to the page on
`building Enzo? </wiki/UserGuide/BuildingEnzo>`_. (If you do follow
the detailed version, be sure you work your way through all the way
to the bottom of the page, and build both ring and inits!)

Here's the abbreviated version, which sets the output style to
packed, builds the excutable, and drops it into a local bin
directory. (The ellipsis (...) indicate additional output.)

::

    ds100 $ cd enzo/
    ds100 $ ./configure 
    ds100 $ cd src/enzo
    ds100 $ make machine-sdsc-datastar
    
     *** Execute 'gmake clean' before rebuilding executables ***
    
       MACHINE: SDSC DataStar
    
    ds100 $ make packed-amr-yes
    
     *** Execute 'gmake clean' before rebuilding executables ***
    
       CONFIG_PACKED_AMR:            yes
    
    ds100 $ make opt-high
    
     *** Execute 'gmake clean' before rebuilding executables ***
    
       CONFIG_OPT:                   high
    
    ds100 $ make clean
    ds100 $ make -j4
    Updating DEPEND
    ...
     Linking
    Success!
    ds100 $ make install
    if [ ! -e ../../bin ]; then mkdir ../../bin; fi
    gmake -s show-flags   >& ../../bin/enzo.show-flags
    gmake -s show-config  >& ../../bin/enzo.show-config
    gmake -s show-version >& ../../bin/enzo.show-version
    gmake -s show-diff    >& ../../bin/enzo.show-diff
    ds100 $

Now we're going to do that in the ring/ and inits/ directories.
First inits:

::

    ds100 $ cd ../inits; make; make install
    Compiling enzo_module.src90
    Updating DEPEND
    ...
    Linking
    Success!
    if [ ! -e ../../bin ]; then mkdir ../../bin; fi
    make show-flags   >& ../../bin/inits.show-flags
    make show-config  >& ../../bin/inits.show-config
    make show-version >& ../../bin/inits.show-version
    ds100 $

Next ring:

::

    ds100 $ cd ../ring; make; make install
    Updating DEPEND
    Compiling Ring_Decomp.C
    ...
    Linking
    Success!
    if [ ! -e ../../bin ]; then mkdir ../../bin; fi
    make show-flags   >& ../../bin/ring.show-flags
    make show-config  >& ../../bin/ring.show-config
    make show-version >& ../../bin/ring.show-version
    ds100 $

We can check that we have the executables we need by looking in the
bin/ directory:

::

    cable:~/tmp/enzo/src/ring rpwagner$ cd ../..     
    cable:~/tmp/enzo rpwagner$ ls bin/
    enzo               enzo.show-version  inits.show-version ring.show-version
    enzo.show-config   inits              ring
    enzo.show-diff     inits.show-config  ring.show-config
    enzo.show-flags    inits.show-flags   ring.show-flags
    cable:~/tmp/enzo rpwagner$ 

Running an Enzo Cosmology Simulation
------------------------------------

After compiling, you should create a directory to run the
simulation in. This is because Enzo cosmology simulations create
quite a few output files, so it's best to store them in their own
directory. For the purposes of this example we assume
that you have created a directory called ``EnzoTestSim`` in your home
directory. You should then download a set of sample parameter
files. The example set used for this tutorial are
`available here <http://lca.ucsd.edu/software/enzo/data/cookbook/>`_.
Download the files called ``Example\_Cosmology\_Sim.inits`` and
``Example\_Cosmology\_Sim.param``, which are the inits and enzo
parameter files, respectively. This tutorial assumes that you have
downloaded these two files and put them on whatever computer you
are using to perform your simulation.

Creating Initial Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step in preparing the simulation is to create the initial
conditions. The file ``Example\_Cosmology\_Sim.inits`` is a text file
that contains a list of parameter file names with their associated
values. These values tell the initial conditions generator necessary
information like the simulation box size, the cosmological
parameters and the size of the root grid. The code then takes that
information and creates a set of initial conditions. ``inits`` is
run by typing this command:

::

    ds100 $ /gpfs/ux455215/Cookbook/enzo/bin/inits -d Example_Cosmology_Sim.inits
    ENZO Inits V64.0 - April 3rd 2006
    
    Reading parameter file
    ...
    successful completion.
    ds100 $ ls
    Example_Cosmology_Sim.inits ParticlePositions
    Example_Cosmology_Sim.param ParticleVelocities
    GridDensity                 PowerSpectrum.out
    GridVelocities

inits will produce some output to the screen to tell you what it is
doing, and will write five files: ``GridDensity``, ``GridVelocities``,
``ParticlePositions``, ``ParticleVelocities`` and ``PowerSpectrum.out``. The
first four files contain information on initial conditions for the
baryon and dark matter componenets of the simulation, and are HDF5
files (formatted binary files). The last file is an ASCII file that
contains information on the power spectrum used to generate the
initial conditions.

Parallel IO Using Ring
~~~~~~~~~~~~~~~~~~~~~~

This example simulation is very small (32\ :sup:`3`\  root grid) so
it is probably not worth using parallel IO. It is definitely
important for larger simulations, though, so we show how to do it
here.
To turn the parallel IO on, these parameters are added into the
Enzo parameter file:

::

    #
    # IO parameters
    #
    ParallelRootGridIO = 1
    ParallelParticleIO = 1

These two parameters turn on parallel IO for both grids and
particles. In a serial IO simulation where multiple
processors are being used, the master processor reads
in all of the grid and particle initial condition information and
parcels out portions of the data to the other processors.
Similarly, all simulation output goes through the master processor
as well.
This is fine for relatively small simulations using only a few
processors, but slows down the code considerably
when a huge simulation is being run on hundreds of processors.
Turning on the parallel IO options allows each processor
to perform its own IO, which greatly decreases the amount of time
the code spends performing IO.

The process for parallelizing grid and particle information is
quite different. Since we know exactly where every
grid cell in a structured Eulerian grid is in space, and these
cells are stored in a regular and predictable order
in the initial conditions files, turning on ParallelRootGridIO
simply tells each processor to
figure out which portions of the arrays in the GridDensity and
GridVelocities belong to it, and
then read in only that part of the file. The particle files
(ParticlePositions and ParticleVelocities)
store the particle information in no particular order, so in order
to efficiently parallelize the particle IO the
ring tool is used. ring is run on the same number of processors as
the simulation that you intend
to run, and can be used right before the simulation itself is run.
In ring, each processor reads in an
equal fraction of the particle position and velocity information
into a list, flags the particles that belong in its
simulation spatial domain,
and then passes its portion of the total list on to another
processor. After each portion of the list has made its
way to every processor, each processor then collects all of the
particle and velocity information that belongs to it
and writes them out into files called PPos.nnnn and PVel.nnnn,
where nnnn is the processor number.
Turning on the ParallelParticleIO flag in the Enzo parameter file
instructs Enzo to look for these files.

For the purpose of this example, I'm going to run ring in parallel,
using four MPI tasks. You run ring on the particle files by
typing:

::

    ds100 $ poe /gpfs/ux455215/Cookbook/enzo/bin/ring pv \
        ParticlePositions ParticleVelocities -nodes 1 -tasks_per_node 4
    Input arg pv should be pv
    PPin = ParticlePositions
    PVin = ParticleVelocities
    Read Position
    Read Velocity
    ...
    Sort completed
    ds100 $ 

This will then produce some output to your screen, and will
generate 8 files:

::

    ds100 $ ls -1 PPos* PVel*
    PPos0000
    PPos0001
    PPos0002
    PPos0003
    PVel0000
    PVel0001
    PVel0002
    PVel0003
    ds100 $ 

Note that if you are using a different machine or platform, you may
use something other than vmirun for MPI-parallel applications.
Consult your system administrator or system documentation for more
information.

Congratulations, you're now ready to run your cosmology
simulation!

Nested Initial Conditions and Particles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When running a
`nested grid cosmology simulation? </wiki/Tutorials/WritingParameterFiles#Multiplenestedgrids>`_,
there can arise an issue of missing particles as a result of
running ring. Please see
`this page? </wiki/Tutorials/NestedGridParticles>`_ for more
information.

Running an Enzo cosmology simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After all of this preparation, running the simulation itself should
be straightforward. You start enzo by typing:

::

    ds100 $ poe /gpfs/ux455215/Cookbook/enzo/bin/enzo -d \
        Example_Cosmology_Sim.param -nodes 1 -tasks_per_node 4 > ExampleSim.log &

The simulation will now run. The -d flag ensures a great deal of
output, so we redirect it into a log file called output.log for
later examination. This particular simulation used to take
approximately two minutes to run on 4 processors on Abe (the NCSA
Linux cluster as of fall 2008). When t



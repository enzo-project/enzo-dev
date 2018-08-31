Deprecated: Using Inits and Other Tools
=======================================

This section has been largely replaced by :ref:`CosmologicalInitialConditions`.

Building inits
--------------

This is the initial conditions generator. See :ref:`using_inits` for more
detailed information. Initial conditions with a single initial grid or multiple
nested grids can be created with this executable.  Output file names are
user-specified, but in a standard cosmology simulation with a single initial
grid there should be a file containing baryon density information, another
containing baryon velocity information, and two more files containing particle
position and velocity information. Simulations with multiple grids will have a
set of these files for each level, appended with numbers to make them unique.

::

    usage: inits [options] param_file
       options are:
          -d(ebug)
          -s(ubgrid) param_file

Building is as simple as the follwing.

::

    ~/enzo/src/ring $ cd ../inits/
    ~/enzo/src/inits $ make
    Compiling enzo_module.src90
    Updating DEPEND
    Compiling acml_st1.src
    ...
    Compiling XChunk_WriteIntField.C
    Linking
    Success!

This will produce ``inits.exe``.


.. _ring:

ring
----

ring must be run on the simulation particle position and velocity
information before a simulation is executed when the Enzo runtime parameter
``ParallelParticleIO`` is set to 1. Running ring generates files called
PPos.nnnn PVel.nnnn where nnnn goes from 0001 to the total number
of processors that are used for the simulation. These files contain
the particle position and velocity information for particles that
belong to each processor individually, and will be read into the
code instead of the monolithic particle position and velocity
files. Note that if ``ParallelParticleIO`` is on and ring is NOT run,
the simulation will crash.

::

    usage:  ring [string] <particle position file> <particle velocity file>

[string] can be one of the following: pv, pvm, pvt, or pvmt. p, v,
m and t correspond to position, velocity, mass, and type,
respectively. The most common [string] choice is 'pv'.
In that case, and if you use the default names for
the particle position and velocity files, your usage will look
like:

::

    ring pv ParticlePositions ParticleVelocities

Building is as simple as the following.

::

    ~/enzo/src/enzo $ cd ../ring/
    ~/enzo/src/ring $ make
    Updating DEPEND
    Compiling Ring_Decomp.C
    Compiling Enzo_Dims_create.C
    Compiling Mpich_V1_Dims_create.c
    Linking
    Success!

This will produce ``ring.exe``.



enzohop
-------

The second (and generally favored) method used for finding density peaks in an
Enzo simulation. More information can be found here. A file called
``HopAnalysis.out`` is output which contains halo position and mass
information.

::

    enzohop [-b #] [-f #] [-t #] [-g] [-d] amr_file
      -b)egin region
      -f)inish region
      -t)hreshold for hop (default 160)
      -g)as particles also used (normally just dm)
      -d)ebug

anyl
----

anyl is the analysis package written in C, previously known as enzo_anyl.
Although the analysis toolkit for enzo that's being constantly updated is YT,
anyl has its own value for some users. It creates radial, disk, vertical
profiles for baryon (each species), dark matter, and star particles. Works with
all AMR formats including HDF4 and packed HDF5.

::

    usage: anyl.exe <amr file> <anyl parameter file>




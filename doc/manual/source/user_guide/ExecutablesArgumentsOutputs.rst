.. _ExecutablesArgumentsOutputs:

Executables, Arguments, and Outputs
===================================

This page is a summary of all of the binaries that are created
after ``make; make install`` is run in the Enzo code bundle. They
should be located in the ``bin`` directory. Links to the various pages
of the manual that describe a particular binary are also included.

enzo
----

This is the main simulation code executable. See :ref:`RunningEnzo`
for more detailed information.

When an Enzo simulation is run, at every datastep several files are
output. There is an ascii file which has no extension (ie, if your
output dumps are named RedshiftOutput then the first parameter file
output will be RedshiftOutput0000). This file contains all of the
parameters that Enzo needs to be able to restart the simulation,
such as cosmology information, redshift, information on box volume,
and which physics modules are turned on. Another file, which has
the same root and the extension '.hierarchy', contains ascii
information on all of the Enzo AMR grids, such as their positions,
sizes, number of particles per grid, etc. There are two files with
extensions '.boundary' and '.boundary.hdf' which contain
information on boundary conditions. And then there will be at least
one file with the extension '.gridNNNN', where NNNN is a number
between 0 and 9999. For simulations with more than 10,000 files,
the numbering will have 5 digits and start from 10000. These files
are where all of the simulation data is actually contained.

::

   usage: ./enzo.exe [options] param_file
      options are:
         -d(ebug)
         -r(estart)
         -x(extract)
            -l(evel_of_extract) level
         -p(roject_to_plane) dimension
         -P(roject_to_plane version 2) dimension
            -m(smooth projection)
         -o(utput as particle data)
         -h(elp)
         -i(nformation output)
         -s(tart  index region) dim0 [dim1] [dim2]
         -e(nd    index region) dim0 [dim1] [dim2]
         -b(egin  coordinate region) dim0 [dim1] [dim2]
         -f(inish coordinate region) dim0 [dim1] [dim2]

inits
-----

This is the initial conditions generator. See :ref:`RunningInits` for more
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

anyl is the analysis package written in C, previously known as enzo\_anyl.
Although the analysis toolkit for enzo that's being constantly updated is YT,
anyl has its own value for some users. It creates radial, disk, vertical
profiles for baryon (each species), dark matter, and star particles. Works with
all AMR formats including HDF4 and packed HDF5.

::

    usage: anyl.exe <amr file> <anyl parameter file>



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

When an Enzo simulation is run, at every datastep several files are output,
inserted into subdirectories.  The most important of these are the files with
no extension and those ending in ``.hierarchy``, of which there will be one of
each for each datadump.  For more information on the format of Enzo output, see
:ref:`EnzoOutputFormats`.

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
      -g (Write Potential field only)
      -M (Write smoothed DM field only)
      -F(riends-of-friends halo finder only)
      -C(ooling time write only)
      -h(elp)
      -i(nformation output)
      -V (show compiler options and flags)
      -s(tart  index region) dim0 [dim1] [dim2]
      -e(nd    index region) dim0 [dim1] [dim2]
      -b(egin  coordinate region) dim0 [dim1] [dim2]
      -f(inish coordinate region) dim0 [dim1] [dim2]

The -g, -M, and -C flags will read in the dataset given on the command
line and write additional data fields to the same data files.  When
running with these flags (or the -F flag), the -r flag must also be
given so that the code knows to read in a dataset.  For example, to
write out the cooling time to the output DD0001, do the following:

::

   enzo.exe -r -C DD0001/DD0001

inits
-----

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



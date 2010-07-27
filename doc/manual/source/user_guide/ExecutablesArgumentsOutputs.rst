Executables, Arguments, and Outputs
===================================

This page is a summary of all of the binaries that are created
after make; make install is run in the Enzo code bundle. They
should be located in the bin/ directory. Links to the various pages
of the manual that describe a particular binary are also included.

enzo
----

This is the main simulation code executable. See the
`page on running Enzo? </wiki/Devel/UserGuide/RunningEnzo>`_ for
more detailed information.

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

    usage: enzo [options] <param_file>
    
       general options:
          -d                            display debug information
          -r                            restart
          -x                            extract
          -l <level>                    level of extract
          -p <dimension>                project to plane
          -m                            smooth projection
          -o                            output as particle data
          -h                            help
          -i                            information output
          -s <dim0> [<dim1> [<dim2>]]   start index region
          -e <dim0> [<dim1> [<dim2>]]   end index region
          -b <dim0> [<dim1> [<dim2>]]   begin coordinates
          -f <dim0> [<dim1> [<dim2>]]   finish coordinate region
    
       performance options:
          -P mode <modeval>             set jbPerf mode
          -P event <eventname>          set jbPerf event
          -P dir <directory>            set jbPerf directory



inits
-----

This is the initial conditions generator. See
`this page? </wiki/Devel/UserGuide/RunningInits>`_ for more
detailed information. Initial conditions with a single initial grid
or multiple nested grids can be created with this executable.
Output file names are user-specified, but in a standard cosmology
simulation with a single initial grid there should be a file
containing baryon density information, another containing baryon
velocity information, and two more files containing particle
position and velocity information. Simulations with multiple grids
will have a set of these files for each level, appended with
numbers to make them unique.

::

    usage: inits [options] param_file
       options are:
          -d(ebug)
          -s(ubgrid) param_file



ring
----

ring must be run on the simulation particle position and velocity
information before a simulation is executed when the parameter
ParallelParticleIO is set to 1. Running ring generates files called
PPos.nnnn PVel.nnnn where nnnn goes from 0001 to the total number
of processors that are used for the simulation. These files contain
the particle position and velocity information for particles that
belong to each processor ind style="text-align: right"> Sedov Blast
(2D unigrid version)
[browser:public/trunk/doc/examples/ShockPool2D.enzo
ShockPool2D.enzo]
2D Shock Propogation test
[browser:public/trunk/doc/examples/ShockPool3D.enzo
ShockPool3D.enzo]
3D Shock Propogation test
[browser:public/trunk/doc/examples/ShockTube.enzo ShockTube.enzo]
ShockTube test (unigrid version)
[browser:public/trunk/doc/examples/SphericalInfall.enzo
SphericalInfall.enzo]
Spherical Infall Test Problem
[browser:public/trunk/doc/examples/StripTest.enzo StripTest.enzo]
Stripping (Collapse) test
[browser:public/trunk/doc/examples/WavePool.enzo WavePool.enzo]
Wave Propogation test
[browser:public/trunk/doc/examples/ZeldovichPancake.enzo
ZeldovichPancake.enzo]
Zeldovich Pancake (unigrid)
Some of these test files include more detailed descriptions as
header comments.

All Enzo test files are also included in the
` lcatest <http://lca.ucsd.edu/projects/lcatest>`_ distribution.
lcatest is a parallel software testing environment, also developed
at LCA but distributed separately.



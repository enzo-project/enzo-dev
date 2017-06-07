.. _RunningEnzo:

Running Enzo
============

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

Once the code is compiled and a parameter file is prepared,
starting the simulation is easy:

::

    mpirun -np 1 enzo [-d] parameter_file

The syntax of the mpirun varies between mpi implementations. The
example given here comes from a machine using a standard MPI
implementation that is initiated by the 'mpirun' command, and implies
the use of a single processors (the argument after the -np flag
indicates the number of processors).




The -d flag triggers a debug option that produces a substantial amount
of output. See :ref:`GettingStarted` for more detailed information on
running simulations. You may also need to use :ref:`ring` if you are
using parallel I/O.

Restarting
----------

During a run, there are a number of forms of output. The largest
will probably be the output of the full dataset as specified by
parameters such as ``dtDataDump`` and the ``CosmologyOutputRedshift``.
Such outputs contain a number of different files (sometimes many files
if there are a large number of grids) and are explained elsewhere.
It is useful to have a fairly large number of such outputs if the
run is a long one, both to provide more information to analyze, but
also in case of an unintended interruption (crash). Fortunately,
any full output can be used to restart the simulation:

::

    mpirun -np 1 enzo [-d] -r output_name

Monitoring information
----------------------

As the simulation runs, at every top grid timestep, it outputs a
line of information to the ascii file OutputLevelInformation (which
is overwritten on restart). The amount of information on this line
can be quite extensive, but here the format is briefly summarized.
The first number is the problem time, while the next 6 relate to
general information about the entire run. Within these six numbers,
the first is the maximum level currently in use, the second is the
number of grids, the third is a number proportional to the memory
used, the fourth is the mean axis ratio of all grids, and the last
two are reserved for future use. Then, there are three spaces,
and another group of numbers, all providing information about the
first (top grid) level. This pattern of three spaces and six
numbers is repeated for every level.  An example of this file is
provided below.

::

      Cycle 151  Time 20.241365  MaxDepth 4  Grids 412  Memory(MB) 53.3117  Ratio 2.22582
         Level 0  Grids 2  Memory(MB) 13.8452  Coverage 1  Ratio 2  Flagged 0  Active 262144
         Level 1  Grids 304  Memory(MB) 31.4977  Coverage 0.166855  Ratio 2.43768  Flagged 0  Active 349920
         Level 2  Grids 76  Memory(MB) 5.81878  Coverage 0.00329208  Ratio 1.66118  Flagged 0  Active 55232
         Level 3  Grids 22  Memory(MB) 1.74578  Coverage 0.000125825  Ratio 1.63561  Flagged 0  Active 16888
         Level 4  Grids 8  Memory(MB) 0.404286  Coverage 2.5034e-06  Ratio 1.21875  Flagged 0  Active 2688

The information for each level is:

#. number of grids on the level
#. memory usage (minus overhead).  Actual memory usage is usually a factor of 10 higher.
#. the volume fraction of the entire region covered by grids on this level,
#. the mean axis ratio of grids on this level
#. the fraction of cells on this level which need refinement (unused)
#. the number of active cells on this level.

Debugging information
---------------------

It is often useful to run with the debug flag turned on,
particularly if the code is crashing for unknown reasons.
However, the amount of output is quite
large so it is useful to redirect this to a log file, such as:

::

    mpirun -np 1 enzo -d -r output_name >& log_file

Some modules (the cooling unit is particularly bad for this),
produce their own debugging logs in the form of fort.?? files.
These can be ignored unless problems occur.

Test Problems
-------------

There are a number of built-in tests, which can be used to debug the
system or characterize how well it solves a particular problem.  (see
:doc:`EnzoTestSuite` for a complete list.) Note that Enzo can run any
problem after compilation, since no compilation flags affect
simulation parameters.  To run a particular test, cd to the
[browser:public/trunk/doc/examples doc/examples] subdirectory of the
Enzo source distribution (after compiling enzo) and use the following
command-line:

::

    mpirun -np 1 enzo [-d] test_name

The syntax of the mpirun various from mpi implementation. The
example given here comes from the Origin2000 and implies a single
processor (the argument after the -np flag indicates the number of
processors).

The parameter test_name corresponds to the parameter file that
specifies the type of test and the test particulars. This file is
ascii, and can be edited.
It consists of a series of lines (and optional comments) each of
which specifies the value of one parameter. The parameters are
discussed in more detail in :ref:`parameters`.

If you just type ``enzo`` without any arguments, or if the number of
arguments is incorrect, the program should respond with a summary
of the command-line usage.

The -d flag turns on a rather verbose debug option.

For example, to run the shock tube test, use:

::

    mpirun -np 1 enzo ShockTube

or

::

    enzo ShockTube

The response should be:

::

    Successfully read in parameter file ShockTube.
    Successful completion...

How do you know if the results are correct?  New for v2.0, we have
added more `regression tests and answer tests
<http://ppcluster.ucsd.edu/lcatest/>`_, using LCAtest.  We hope to
add more answer tests, especially for large production-type
simulations, e.g. a 512\ :sup:`3` cosmology simulation.



.. _RunningEnzo:

Running Enzo
============

Once the code is compiled and a parameter file is prepared,
starting the simulation is easy:

::

    mpirun -np 1 enzo [-d] parameter_file

The syntax of the mpirun varies between mpi implementations. The
example given here comes from a machine using a standard MPI
implementation that is initiated by the 'mpirun' command, and
implies the use of a single processors (the argument after the -np
flag indicates the number of processors).

The -d flag triggers a debug option that produces a substantial
amount of output. See :ref:`Tutorials` for more detailed information on running
simulations. You may also need to use ring (see
:ref:`ExecutablesArgumentsOptions`) if you are using parallel I/O.

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
numbers is repeated for every level. The information for each level
is:


#. number of grids on the level
#. memory usage (roughly)
#. the volume fraction of the entire region covered by grids on
   this level,
#. the mean axis ratio of grids on this level (please don't ask why
   this information is here)
#. the fraction of cells on this level which need refinement (I
   think),
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

There are a number of built-in tests, which can be used to debug the system or
characterize how well it solves a particular problem.  (see
:ref:`EnzoTestSuite` for a complete list.) Note that Enzo can run any problem
after compilation, since no compilation flags affect simulation parameters
(unlike the hydrodynamics code KRONOS for example) To run a particular test, cd
to the [browser:public/trunk/doc/examples doc/examples] subdirectory of the
Enzo source distribution (after compiling enzo) and use the following
command-line:

::

    mpirun -np 1 enzo [-d] test_name

The syntax of the mpirun various from mpi implementation. The
example given here comes from the Origin2000 and implies a single
processor (the argument after the -np flag indicates the number of
processors).

The parameter test\_name corresponds to the parameter file that
specifies the type of test and the test particulars. This file is
ascii, and can be edited.
It consists of a series of lines (and optional comments) each of
which specifies the value of one parameter. The parameters are
discussed in more detail in :ref:`EnzoParameters`.

If you just type enzo without any arguments, or if the number of
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

How do you know if the results are correct? We hope to add the
ability for the code to check against pre-computed results, but for
the moment, if the code doesn't crash, it's probably a reasonable
bet that it is working correctly (the other good check is that the
shock tube run takes 68 or 69 steps). You may also wish to plot the
results (output as HDF files). This section will be expanded in
later editions of this manual.



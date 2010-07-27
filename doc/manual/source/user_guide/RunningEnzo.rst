Running Enzo
============

`TOC(heading=This Page)? </wiki/TOC(heading=This%20Page)>`_

Once the code is compiled and a parameter file is prepared,
starting the simulation is easy:

::

    mpirun -np 1 enzo [-d] parameter_file

The syntax of the mpirun varies between mpi implementations. The
example given here comes from a machine using a standard MPI
implementation that is initiated by the 'mpirun' command, and
implies the use of a single processors (the argument after the -np
flag indicates the number of processors.

Once again, the -d flag triggers a debug option that produces a
substantial amount of output. See the
`Tutorials? </wiki/Tutorials>`_ for more detailed information on
running simulations. You may also need to use
`ring? </wiki/Devel/UserGuide/ExecutablesArgumentsOutputs>`_ if you
are using parallel I/O.

Restarting
----------

During a run, there are a number of forms of output. The largest
will probably be the output of the grid hierarchy as specified by
parameters such as dtDataDump and the CosmologyOutputRedshift. Such
outputs contain a number of different files (sometimes many files
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
two are reserved for future use. Then, there will be three spaces,
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
particularly if the code is crashing for unknown reasons (which
never happens of course). However, the amount of output is quite
large so it is useful to redirect this to a log file, such as:

::

    mpirun -np 1 enzo -d -r output_name > log_file

Some modules (the cooling unit is particularly bad for this),
produce their own debugging logs in the form of fort.?? files.
These can be ignored unless problems occur. We hould probably
stress again that this code is not guaranteed to bug-free and is
offered without any assurances or warranties, implied or
otherwise.

Test Problems
-------------

There are a number of built-in tests, which can be used to debug
the system or characterize how well it solves a particular problem.
(The
a class="missing wiki" href="/wiki/Devel/UserGuide/EnzoTestSu


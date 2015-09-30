.. _PerformanceMeasurement:

Performance Measurement
=======================

EnzoTiming.h and performance_tools
----------------------------------

This framework consists of two pieces -- one that is part of enzo, primarily
contained in src/enzo/EnzoTiming.h, and another which is used to plot and
analyze the performance data, in src/performance_tools/performance_tools.py.


Usage Overview
##############

We have added support for simple, lightweight measurements for the timing and
performance of Enzo.  This allows one to examine which functions are using the
majority of the simulation runtime, and how this varies across multiple
processors. We have built in a number of default timers, such as EvolveLevel for
each level, RebuildHierarchy, SolveHydroEquations, and Group_WriteAllData.
Below, we will outline how to add additional timers and how to generate plots of
the data.

File Format
###########

At each cycle, information is printed out to a file named performance.out.  It
collects the amount of time taken on each of the processors to complete the
listed functions (e.g. Level N EvolveLevel, RebuildHiearchy, etc.) over that
cycle.  Rather than giving all of the values returned by each processor for a
given function, EnzoTiming only outputs the mean amount of time spent per
processor, the maximum & minimum amount of time across processors, and the standard
deviation of this distribution of times.  This is meant to give the user a sense
of how well load-balanced their simulation is across processors, as well as
pinpoint where the majority of the time is being spent.  To explain the output,
we show an example cycle from performance.out:

::

  Cycle_Number 2
  Level_0 6.520748e-05 8.344650e-07 6.389618e-05 6.604195e-05 100 4 3.833916e+05
  Level_1 3.254414e-05 2.804866e-05 1.406670e-05 8.106232e-05 10 1 7.681875e+04
  Level_2 1.159906e-04 2.678922e-05 9.965897e-05 1.623631e-04 14 1 3.017485e+04
  Level_3 2.477765e-04 7.348677e-05 2.028942e-04 3.750324e-04 16 1 1.614358e+04
  Level_4 5.816817e-04 1.630557e-04 4.820824e-04 8.640289e-04 24 1 1.031492e+04
  Level_5 1.266718e-03 3.594168e-04 1.056910e-03 1.889229e-03 26 1 5.131371e+03
  Level_6 2.686501e-03 7.197988e-04 2.262831e-03 3.933191e-03 40 1 3.722315e+03
  RebuildHierarchy 5.715549e-03 1.371242e-04 5.478144e-03 5.801201e-03
  SolveHydroEquations 1.436710e-03 2.407243e-03 4.386902e-05 5.606174e-03
  Total 1.499003e-02 3.440975e-05 1.494408e-02 1.503992e-02 230 10 3.835882e+03

| Each of the Level_N and Total lines have:
| Level_N, mean time, stddev time, min time, max time, number of cell updates, 
| number of grids, mean cell updates/s/processor

| Each non-level line (RebuildHierarchy, SolveHydroEquations, etc.) have:
| Section Name, mean time, stddev time, min time, max time. 

Time is measured in seconds of wall time for each of the processors.

In the example above, we see that more time is being spent in RebuildHierarchy 
than in SolveHydroEquations, and that the load balance is quite poor for the
SolveHydroEquations where the mean is 1.4 ms, with a standard deviation of
2.4 ms. 

At the beginning of each simulation (on Cycle 1), we print out a header to the
performance.out file:

:: 

  # This file contains timing information
  # For instructions on how to decipher this information,
  # see [enzo base directory]/src/performance_tools/README.
  # Times are collected across MPI processes and presented as:
  # Level_N/Total, mean time, std_dev time, min time, max time, cell updates, grids, cell updates/processor/sec
  # Routine, mean time, std_dev time, min time, max time 

Then, at the start of each simulation (whether the beginning or a restart), we
print out the MPI processor count:

::

  # Starting performance log. MPI processes: 4

This is done in case the number of processors changes over time.

Adding New Timers
#################

While there are a number of default timers, it is easy to add new timers to any
section of code in Enzo.

The built-in timers include: EvolveHierarchy (Total), EvolveLevel (for each 
level), SolveHydroEquations, RebuildHierarchy, and Group_WriteAllData.  Adding 
new times should be as simple as doing two things:

1) Add 

.. code-block:: c

  #include "EnzoTiming.h" 

to the top of the file you want to profile,
making sure it is before macros_and_parameters.

2) Add 

.. code-block:: c

  TIMER_START("YourTimerName");
and

.. code-block:: c

  TIMER_STOP("YourTimerName");

around the code you want to time.  And adding an initializer statement to
enzo.C (along with the other timer initializers):

.. code-block:: c

  TIMER_REGISTER("YourTimerName");

The string that you pass in gets collected in a map which is then iterated over
at the end of each evolve hierarchy.  At that time it prints into a file named
performance.out.

Generating Plots
################

performance_tools.py (located in src/performance_tools) is a python module 
for plotting the performance information stored in performance.out.  The easiest
way to generate plots from performance.out is to call performance_tools.py from
the command line:

:: 
    
  python performance_tools.py performance.out
 
| or 

:: 

  python performance_tools.py -s 11 performance.out

to do the same while applying a smoothing kernel to your data 11 cycles in 
width.

By default, performance_tools.py will output 8 plots: 

--p1.png
  Plot the mean time taken per processor on each level and on the 
  simulation as a whole (Total) versus cycle number.  Overplot in 
  lighter tones are the minimum and maximum time taken on a processor 
  for each of these quantities.

--p2.png
  Same as p1.png except scale everything to be as a fraction of the 
  total time taken.

--p3.png
  Plot the mean time taken per processor on each level versus cycle number.  
  Stack each level on the previous layer cumulatively.  

--p4.png
  Plot the mean time taken per processor performing any Non-Level fields versus
  cycle number (e.g. the RebuildHiearchy, SolveHydroEquations, and 
  Group_WriteAllData tasks).  Stack each level on the previous layer 
  cumulatively.  Scale everything to be as a fraction of the total time taken.

--p5.png
  Plot the number of cells updated at each level versus cycle number and 
  stack them cumulatively.

--p6.png
  Plot the efficiency (cell updates/processor/sec) for each level and for
  the simulation as a whole versus cycle number.

--p7.png
  Plot the load balancing (Max Time - Min Time) for all subprocesses and 
  levels of the simulation as a whole versus time.  

--p8.png
  Plot the load balancing (Max Time - Min Time) for all subprocesses and 
  levels of the simulation as a whole versus time.  Normalize them by the 
  mean time taken for each process.

Generating Additional Plots
###########################

If you want to create additional plots of your data beyond the defaults, 
simply add new plot_quantities() and plot_stack() calls to the bottom of 
performance_tools.py.

This can be as simple as adding one of these lines:

.. code-block:: python
  
  # Plot the mean time taken per processor on Level 0 EvolveLevel calls versus
  # Cycle Number.
  p.plot_quantity("Level 0", "Mean Time")

  # Same as above, but stacks the quantity from zero to the mean time.
  p.plot_stack("Level 0", "Mean Time")

  # Plot the mean time take per processor for all defined fields (All levels,
  # All Functions)
  p.plot_quantity([], "Mean Time", repeated_field="All")
  
  # Plot and stack cumulatively on top of each other the number of cell
  # updates for each level versus cycle number.
  p.plot_stack([], "Cell Updates", repeated_field="Level")

  # Plot the mean time taken per processor for all non-level functions versus
  # cycle number (including "Total" time taken by everything).
  p.plot_quantity("Total", "Mean Time", repeated_field="Non-Level")

Full documentation for the plot_quantity and plot_stack functions can
be found in the docstrings for the performance_tools.py module.  You can
view it either by looking at the source code, or by loading it in python:

.. code-block:: python

  import performance_tools as pt
  help(pt.perform)

Additional Performance Tools
############################

An additional performance tool exists in the enzo source which provides
slightly different details about a number of subprocesses in enzo.  By default,
it is turned off, although you can enable it by adding a compiler flag to your 
Makefile.  The downside to this 
performance tool, called MPI Instrumentation, is that it only provides information
if you're running MPI, it only gives you that information at the end of a 
simulation (when it has successfully completed, not when it runs out of time),
and it produces a file for every processor that was used in the simulation 
(which can sometimes crowd your directory).

To enable this feature, you can add a flag to your machine's Makefile to explicitly
set this preprocessor keyword.  Do this by editing your machine's 
Makefile to include this flag to your MACH_DEFINES:

.. code-block:: C++

  MACH_DEFINES = -DMPI_INSTRUMENTATION

In case you want to see what sort of information is provided by MPI Instrumentation,
a sample output file is included below:

.. code-block:: bash

  Elapsed wall time:                   3.582540e+03
  Communication time:                  1.617045e+03
  Global communication time:           9.343419e+02
  Receive communication time:          4.590317e+00
  Waiting communication time:          0.000000e+00
  
  
  Transferring region       ( 1940795 times) 4.588604e+00
  Sending particles         (    1592 times) 6.079674e-04
  Transferring particles    (    9598 times) 5.879667e+01
  Transferring Fluxes       (   32369 times) 9.276295e-02
  ShareGrids                (    5777 times) 8.463278e+01
  Transpose                 ( 1771716 times) 1.597000e+02
  BroadcastValue            (    4915 times) 1.144109e-01
  MinValue                  (   46066 times) 7.745399e+02
  UpdateStarParticleCount   (    5770 times) 1.625819e+01
  
  
  RebuildHierarchy          (    1626 times) 1.555615e+01
  RebuildHierarchy interval (    1626 times) 7.773995e-02
  Load balancing            (       0 times) 0.000000e+00
  Region transfer size      ( 1940795 times) 9.709615e+09
  Particles sent            (    1592 times) 0.000000e+00
  Particle transfer size    (    9598 times) 1.039000e+04
  
    
  Number of load balancing calls 0/0 (LOAD_BALANCE_RATIO=0.000000)
  Number of flagging cells  (    5418 times) 4.116929e+07
  
  
  Average percentage of flagging cells 2.420569e-01(= 1.311464e+03/5418)
  Average percentage of moving cells 0

| Samuel Skillman (samskillman at gmail.com) 
| Cameron Hummels (chummels at gmail.com)


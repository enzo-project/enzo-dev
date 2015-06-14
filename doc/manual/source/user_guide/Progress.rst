.. _Progress:

Measuring Simulation Progress
=============================

Measuring the progress of an Enzo simulation can be tricky as each level of the 
hierarchy has its own timestep and a lot of information is printed.  Fortunately, 
the Enzo source comes with its own progress meter which provides a great deal of 
useful information on the state of a simulation.  The progress meter is called 
**np** and is located in the *bin* directory of the Enzo source.

Running Enzo with the Progress Meter
------------------------------------

To get the most out of the progress meter, simulations should be run with the -d 
flag (for debug output) and have both the standard output and standard error piped 
into a single file.  Note, running Enzo with -d will not slow down the simulation.  
In a bash environment, the standard output and error can be piped into the same 
file in the following way:

::

   [mpirun ...] ./enzo.exe -d AMRCosmology.enzo >& estd.out

Using the Progress Meter
------------------------

To use the progress meter, simply run it from within the simulation directory.

::

   ~ ./np -t "AMRCosmology" -l 5
   +----------------------------------------- AMRCosmology -----------------------------------------+
   | Sat Oct 26 01:20:00 1985                     Status: 30.996% complete.                         |
   +------------------- Time ------------------+ +-------------------- Output ----------------------+
   |  Initial  |  Current  |   Final   | Units | |      |   Time    | Redshift |  Name  | Completed |
   +-----------+-----------+-----------+-------+ +------+-----------+----------+--------+-----------|
   | 8.163e-01 | 7.175e+01 | 2.297e+02 | code  | | Last | 7.082e+01 | 1.558275 | DD0014 |  -------  |
   | 4.911e+07 | 4.316e+09 | 1.382e+10 | years | | Next | 7.582e+01 | 1.438446 | DD0015 | 18.61270% |
   +-------------------------------------------+ +--------------------------------------------------+
   +--------------------------- Hierarchy --------------------------+ +--------- Redshift ----------+
   | L | Grids |  Volume   |    dt     |  Sub  | Completed | Iter |R| | Initial | Current |  Final  |
   +---+-------+-----------+-----------+-------+-----------+------+-+ +---------+---------+---------+
   | 0 |     4 | 1.000e+00 | 2.021e+00 | 1.000 | 1.0000000 |  159 | | | 50.0000 | 1.53496 | 0.00000 |
   | 1 |    49 | 1.500e-01 | 6.820e-01 | 0.688 | 0.6879505 |   54 | | +-----------------------------+
   | 2 |    16 | 4.997e-03 | 2.297e-01 | 0.337 | 0.4641940 |   67 | |
   | 3 |     9 | 1.779e-04 | 9.501e-02 | 1.000 | 0.4641943 |  150 | |
   | 4 |     4 | 3.755e-06 | 9.501e-02 | 1.000 | 0.4641941 |  271 | |
   | 5 |     1 | 2.012e-07 | 2.852e-02 | 0.919 | 0.4603860 |  457 |<|
   +----------------------------------------------------------------+
   +---------------------------------------------
   | TransferSubgridParticles[5]: Moved 0 particles, 0 stars.
   | DetermineSGSize: MaxSubgridSize = 2000, MinSubgridEdge = 4, ncells = 216
   | RebuildHierarchy[5]: Flagged 0/1 grids. 0 flagged cells
   | Level[5]: dt = 0.0285228  0.0285228 (0.0873183/0.0950141)
   | RebuildHierarchy: level = 5
   +---------------------------------------------

The progress meter will continue to update automatically as more information is 
written to the log file.

Progress Meter Output
---------------------

The progress meter has four section: time, output, hierarchy, and redshift 
(if a cosmology simulation.)  The time section gives the initial, current, and 
final time in both code units and years.  The redshift section gives the initial, 
current, and final redshift of the simulation.  The output section gives the time, 
redshift, and names of the previous and next data dump as well as the percentage 
that the simulation is to reaching the next output.  The hierarchy section 
displays, for each level of the hierarchy, the number of grids, the total volume 
of grids, the current timestep, the completion fraction of the level above, the 
completion fraction of the root grid timestep, and the number of iterations taken.  
The far right column shows what level is being computed and the current status.  
See below for an explanation of the symbols.  If the -l flag is given, an 
additional section will appear with the last lines written to the log file.

Additional Options
------------------

Additional options can be seen by running the progress meter with the -h flag.

::

   ~ ./np -h
   np:
           -h: print this help output.
           -d <directory>: simulation directory (default: .).
           -hf <filename>: hierarchy file (default: OutputLevelInformation.out).
           -l <number of output lines>: print enzo standard out lines (default: 0).
           -of <filename>: enzo standard out file (default: estd.out).
           -ol <filename>: enzo output log file (default: OutputLog).
           -pf <filename>: parameter file (default: amr.out).
           -t <title>: title of simulation.
           -w <seconds>: change number of seconds between output (default: 1).
   Status:
           E: Evolve Level
           R: Rebuild Hierarchy
           W: Writing Data
           .: Evolve Level Complete


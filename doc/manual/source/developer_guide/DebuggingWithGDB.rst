Debugging Enzo with GDB
=======================

While it is relatively straightforward to debug enzo in parallel with a
commercial parallel debugger like Totalview or DDT, it is not quite as 
straightforward to debug enzo with a free, open source serial debugger like GDB.
This method works well if you do not have access to a supercomputer or cluster
with a commercial parallel debugger installed, if you would like to run and
debug enzo on a small workstation, or if you prefer to use free and open source
software in your programming life.

There are two general approaches for parallel debugging of Enzo within
GDB, running multiple GDB processes that each run Enzo, or attaching
GDB to an existing Enzo process.


I. Running multiple GDB processes that each run Enzo
----------------------------------------------------

This option works best when running on a single workstation, or on a
cluster to which you have direct access.  The method works best when
running with only a few processors (as will be seen below).

First, build Enzo with debugging symbols enabled and with compiler
optimizations turned off.  This can be accomplished on most systems by
setting ``make opt-debug`` at the command line (see :ref:`MakeOptions`).

Second, launch a number of xterms using ``mpirun`` or ``mpiexec`` that
each internally launch GDB on the Enzo executable::

  18:16:32 [dreynolds@zeno ~]$  mpirun -np 4 xterm -e gdb ./enzo.exe

This will launch 4 xterms, each of which is running a separate gdb
process, that in turn is set to run Enzo.

Within each of these xterms, enter the remaining command-line
arguments needed to run enzo, e.g.::

    (gdb) run -d -r DD0096/DD0096

Once you have hit [enter] in each terminal Enzo will start, with all
process-specific output displayed in it's own xterm.  If you wish to
set breakpoints, these GDB commands should be entered at the various
GDB prompts prior to issuing the ``run`` command.

NOTE: It is possible to insert all of your GDB commands into a GDB
script file, and then have each process run the same script,
eliminating the need to type the commands separately within each
xterm.  To do this, create a file with all of your GDB commands (in
order, one command per line); let's call this file ``gdb.in``.  Then
when you start ``mpirun``, you can specify this script to the GDB
processes::

  18:16:32 [dreynolds@zeno ~]$  mpirun -np 4 xterm -e gdb -x gdb.in ./enzo.exe




II. Attaching GDB to existing Enzo processes
--------------------------------------------


Modify Enzo to allow GDB to attach to a running Enzo process
------------------------------------------------------------

Open ``enzo.C``, located in the main Enzo source directory, and modify the
beginning of the ``MAIN_NAME`` function (the main function where execution
begins) so it looks like the following:

.. code-block:: c

   Eint32 MAIN_NAME(Eint32 argc, char *argv[])
     {
     int i;
     // Initialize Communications
     CommunicationInitialize(&argc, &argv);

   #define DEBUG_MPI
   #ifdef DEBUG_MPI
     if (MyProcessorNumber == ROOT_PROCESSOR) {
       int impi = 0;
       char hostname[256];
       gethostname(hostname, sizeof(hostname));
       printf("PID %d on %s ready for debugger attach\n", getpid(), hostname);
       fflush(stdout);
       while (impi == 0)
         sleep(5);
       }                                                                                                                                  
   #endif

All you should need to do is uncomment the ``#define DEBUG_MPI`` line.  This
code block will make Enzo print the name of the host its being run on and the
process ID number.  You will need both of these pieces of information when you
try to attach to Enzo with GDB.

Once you've modified ``enzo.C``, you will need to rebuild Enzo.  If you haven't
done so already, you should make sure Enzo is built with debugging symbols and
with compiler optimizations turned off.  This can be accomplished on most
systems by setting ``make opt-debug`` at the command line (see
:ref:`MakeOptions`).

Run Enzo
--------

Now you're ready to run a test simulation.  This method works best when using
only a few processors, so don't start a simulation with hundreds of processors
and try to attach to it with GDB unless you know what you're doing.  If you're
running Enzo on a cluster, make sure that you can SSH into the compute nodes.
If not then this debugging method will not work.  Start Enzo normally using
``mpirun``, Enzo should print something like::

  humperdinck:GDB_test goldbaum$ mpirun -np 4 ./enzo.exe -d -r DD0096/DD0096
  MPI_Init: NumberOfProcessors = 4
  PID 34352 on humperdinck.ucolick.org ready for debugger attach

This says that Enzo is running on four cores and has a process ID number of
34352 on the host humperdinck.ucolick.org.  

Attach and Debug With GDB
-------------------------

Next, in a new terminal window, you
should ssh into the appropriate host.  If you're running on your local
workstation there is no need to ssh.  Next, start a GDB session and attach to
the appropriate PID number::

  humperdinck:enzo goldbaum$ gdb
  GNU gdb 6.3.50-20050815 (Apple version gdb-1515) (Sat Jan 15 08:33:48 UTC 2011)
  Copyright 2004 Free Software Foundation, Inc.
  GDB is free software, covered by the GNU General Public License, and you are
  welcome to change it and/or distribute copies of it under certain conditions.
  Type "show copying" to see the conditions.
  There is absolutely no warranty for GDB.  Type "show warranty" for details.
  This GDB was configured as "x86_64-apple-darwin".
  (gdb) attach 34398

GDB should report a long list of warning messages about code in libraries that
enzo links against that was not compiled with debugging symbols.  It's safe to
ignore these errors since we will only be debugging the enzo source.  GDB is now
attached to enzo's process and is probably stuck somewhere in your system's
implimentation of the sleep() function.  To see the execution stack, tell GDB to
print a stack trace::

  0x00007fff8730da6a in __semwait_signal ()
  (gdb) backtrace
  #0  0x00007fff8730da6a in __semwait_signal ()
  #1  0x00007fff8730d8f9 in nanosleep ()
  #2  0x00007fff8735a9ac in sleep ()
  #3  0x0000000100008cee in main (argc=4, argv=0x7fff5fbfef70) at enzo.C:259
  (gdb) 

In this example GDB is stuck three levels down from where we want to be inside
enzo.C.  Move up the stack::

  (gdb) up 3
  #3  0x0000000100008cee in main (argc=4, argv=0x7fff5fbfef70) at enzo.C:259
  259      sleep(5);
  Current language:  auto; currently c++
  (gdb) l
  254    char hostname[256];
  255    gethostname(hostname, sizeof(hostname));
  256    printf("PID %d on %s ready for debugger attach\n", getpid(), hostname);
  257    fflush(stdout);
  258    while (impi == 0)
  259      sleep(5);
  260  }
  261#endif
  262  
  263
  (gdb)

Now GDB is at line 259 of Enzo.C.  To break the infinite loop, you will need to
modify ``impi`` so that it is no longer zero::

  (gdb) set var impi = 7

At this point you can continue execution by typing ``continue`` or ``c``.  If
you want you can also optionally set a breakpoint elsewhere in the enzo source
tree::

  (gdb) break EvolveLevel.C:738

This will pause execution right before Enzo enters ``RebuildHierarchy`` for the
first time.

That should be enough to get you going.  It's also possible to start multiple
GDB processes so you can attach to all of the parallel MPI processes.  See the
GDB docs and the openmpi FAQ page for more information.


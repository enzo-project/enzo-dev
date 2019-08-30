.. _hybrid-parallel:

Running in Hybrid MPI/OpenMP Mode
=================================

Enzo has a functional, if rudimentary, capacity to take advantage of
MPI processes and OpenMP threads simultaneously. 

How To Compile
--------------

First, navigate to your enzo source directory and set your configuration
settings to compile in OpenMP mode with:

::

   make openmp-yes
   
Then, if you're using the Intel compilers add the following to your machine-
specific makefile: 

::

    MACH_OPENMP = -qopenmp

Or if you're using the gcc compilers: 

::

    MACH_OPENMP = -fopenmp

Without this flag, all OpenMP code in the source is ignored as if it were a
comment. 
    
How it Works
------------

OpenMP uses threads within a single process, rather than multiple processes
(like MPI). The threads all share memory, which decreases the communication
overhead and memory footprint of Enzo, which can be very problematic at
larger core counts. Hybrid mode exists to mitigate these issues.


How to Run
----------

On most systems, it is enough to compile according to the instructions above
and set this environment variable:

::
   
   OMP_NUM_THREADS=8

replacing '8' with however many threads you want to use. Leaving this blank
will default to the maximum cores available on your system, so make sure you
set this before running Enzo. If you want to run with only MPI, rather than
re-compiling Enzo, you can just set this variable to 1.

*Quick Note Re: How many Cores You Should Use*
It may be tempting to just use as many cores as your hardware supports.
However, Enzo's OpenMP mode shows diminishing returns in performance past 32 threads,
with the scaling bottoming out at 48 cores. For best results, pick something between
32-48 cores.


Compiler Support
----------------

Hybrid Enzo has only been tested on the GCC and Intel compilers, and relies
on features from OpenMP 4.5, which is officially supported by GCC 6.1, and Intel
17.0, 18.0, and 19.0. However, there have been issues regarding the Intel compilers
failing to compile certain OpenMP pragmas that are legal according to the standard.
Additionally, a bug exists in the 19.0.2 version that causes Enzo to use more stack
memory than it should, which can cause a memory overflow crash. For best results,
always use the latest version of whichever compiler you choose. 

Developing Hybrid Mode Further
------------------------------

At time of writing, only a fraction of Enzo's code has been threaded. Obviously this
means that the rest of the code runs serially (within the process), which prevents Enzo
from achieving ideal performance scaling. While, ideal performance is most likely
impossible, as some sections of the code will have to run serially for a whole host of
reasons, successfully threading more and more of the code can lead to large performance
gains. If you wish to help develop Hybrid Mode further, there are certain pitfalls that
you should be aware of before you get started.

1. Avoid placing parallel sections inside for loops. Every OpenMP parallel section starts
   and ends by spawning or destroying a bunch of threads, which creates overhead. Placing
   these sections inside of a for loop causes those threads to be spawned and destroyed
   everytime the loop iterates, which can waste a lot of processor time if you have an
   especially long list. Best practice is typically to place the start of your parallel
   section outside any for loops, and place an ``omp for`` right before your loops
   start.

Bad practice, don't do this

::
   
   for(int i=0;i<10000;i++){
   #pragma omp parallel for
       // threads created here everytime outer loop runs
       for (int j=0;j<10000;j++){
	  doSomething(i,j);
       }
       // threads destroyed here everytime outer loop runs
   #pragma omp parallel for
       // threads created here everytime outer loop runs      
       for (int j=0;j<10000;j++){
           doSomethingElse(i,j);
       }
       // threads destroyed here everytime outer loop runs
    }
				

This is fine
::

      #pragma omp parallel
      // threads created here once
      for(int i=0;i<10000;i++){
      #pragma omp for
	  for (int j=0;j<10000;j++)
	      doSomething(i,j);
	  

      #pragma omp for
	  for (int j=0;j<10000;j++)
	      doSomethingElse(i,j);
              
      } // threads destroyed here once
	      


This is also fine
::

      #pragma omp parallel for
      // threads created here once
      for(int i=0;i<10000;i++){       
	  for (int j=0;j<10000;j++)       
	      doSomething(i,j);
	  for (int j=0;j<10000;j++)       
	      doSomethingElse(i,j);
      }// threads destroyed here once		       
       


2. OpenMP and MPI do not mix well. DO NOT call MPI routines from inside OpenMP
   parallel sections unless you know exactly what you are doing. Sends and Recieves
   will be sent and recieved out of order, which can cause a lot of issues for
   MPI. If you are getting a ``MPI_ERR_TRUNCATE``, it is very likely that an MPI_Send
   and/or MPI_Recv is being made from inside an OpenMP parallel section. Unless
   you have a specific plan for getting OpenMP and MPI to play nice, stick to
   keeping communication routines and thread parallelism separate. 

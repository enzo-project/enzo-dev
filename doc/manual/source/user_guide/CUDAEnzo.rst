.. _CUDAEnzo:

Running Enzo with CUDA
======================

Enzo contains a CUDA MHD solver. Supported parameters include

    - Riemann Solver: HLL-PLM
    - Gravity
    - Color fields: chemistry, etc
    - Driving field
    - Comoving coordinates

How to compile with CUDA
------------------------

In order to make Enzo compatible with CUDA, a few changes need to be
made in the settings. 

The first thing is to locate the Makefiles located within the ``src/`` 
directory in the Enzo repository. Make sure then that ``CONFIG_ECUDA=yes``

Then locate your machine specific Makefile, e.g. Make.mach.mymach, and 
then set the following variables:

::

    MACH_OPT_CUDA = -arch=sm_20 -03
    MACH_LIBS_CUDA = -L/your/cuda/install/dir/lib64 -lcudart

If using the Kepler K20 GPU, then 

::

    MACH_OPT_CUDA = -arch=sm_35 -03

Last thing to note is that the CUDA solver is single precision only. 
Check to make sure that all the precision flags are set correctly.
For example,

::

    MACH_FFLAGS_INTEGER_32 = 
    MACH_FFLAGS_INTEGER_64 = -i8
    MACH_FFLAGS_REAL_32 = 
    MACH_FFLAGS_REAL_64 = -r8


How to run with CUDA
--------------------

The only thing to do is to set ``UseCUDA=1`` in whichever parameter
file. That's all!

Be sure that each node has at least 1 NVIDIA GPU. Also note that 
although each GPU can be running multiple MPI processes, the 
performance will not increase beyond  1 MPI process per GPU.



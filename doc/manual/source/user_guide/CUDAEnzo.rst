.. _CUDAEnzo:

Running Enzo with CUDA
======================

Enzo contains CUDA version of PPM and MHD solver. Supported parameters include

    - PPM: TwoShock, HLL, HLLC, PPMFlatteningParameter, PPMSteepeningParameter, DualEnergyFormalism, RiemannSolverFallback
    - MHD: HLL-PLM
    - Gravity
    - Color fields: chemistry, etc
    - Driving field
    - Comoving coordinates
    - EOSType: 0

How to compile with CUDA
------------------------

In order to make Enzo compatible with CUDA, a few changes need to be
made in the settings. 

Currently Enzo CUDA is only compatible with 32-bit precision.  In order to
correctly set this, make sure that in src/enzo/ you run the following 
configuration commands:

::

    make cuda-yes
    make integers-32
    make precision-32
    make particles-32
    make particle-id-32
    make inits-32
    make io-32

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
each GPU can be running multiple MPI processes, the 
performance will typically increase with mulitple MPI processes per GPU. 
So it's recommended to set the number of MPI processes per node to be the number 
of CPU cores to fully use both the CPU and GPU resources.
Furthermore, on Kepler K20 GPU, it's recommended to turn on CUDA MPS (Multi-Process Service),
which enables concurrent running of multiple MPI processes on the GPU.


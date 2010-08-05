.. _CompilationRequirements:

Enzo Compilation Requirements
=============================

Enzo can be compiled on any POSIX-compatible operating system, such as Linux,
BSD (including Mac OS X), and AIX.  In addition to a C/C++ and Fortran-90
compiler, the following libraries are necessary:

   * `HDF5 <http://hdf.ncsa.uiuc.edu/HDF5/>`_, the hierarchical data format.
     Note that HDF5 also may require the szip and zlib libraries, which can be
     found at the HDF5 website.  Note that compiling with HDF5 1.8 or greater
     requires that the compiler directive ``H5_USE_16_API`` be specified;
     typically this is done with ``-DH5_USE_16_API`` and it's set in most of
     the provided makefiles.
   * `MPI <http://www-unix.mcs.anl.gov/mpi/>`_, for multi-processor parallel
     jobs.  Note that Enzo will compile without MPI, but it's fine to compile
     with MPI and only run oon a single processor.

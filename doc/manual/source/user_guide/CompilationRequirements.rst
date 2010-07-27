Enzo Compilation Requirements
=============================

Enzo can be compiled on any POSIX-compatible operating system, such as Linux,
BSD (including Mac OS X), and AIX.  In addition to a C/C++ and Fortran-90
compiler, the following libraries are necessary:

   * `HDF5 <http://hdf.ncsa.uiuc.edu/HDF5/>`_, the hierarchical data format.
     Note that HDF5 also may require the szip and zlib libraries, which can be
     found at the HDF5 website.  (See :ref:`hdf5_versions` for a note about
     compiling with HDF5 1.8 or greater.)
   * `MPI <http://www-unix.mcs.anl.gov/mpi/>`_, for multi-processor parallel
     jobs.

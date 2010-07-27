Enzo compilation requirements
=============================

Enzo can be compiled on any POSIX-compatible operating system, such
as Linux, BSD (including Mac OS X), and AIX. In addition to a C/C++
and Fortran-90 compiler, the following libraries are necessary:


-  ` HDF5 <http://hdf.ncsa.uiuc.edu/HDF5/>`_, the hierarchical data
   format. Note that HDF5 also may require the szip and zlib
   libraries, which can be found at the HDF5 website. (Please read
   `this note? </wiki/Devel/UserGuide/BuildingEnzo#HDF5Versions>`_
   about how to handle compiling with HDF5 versions 1.8+.)
-  ` MPI <http://www-unix.mcs.anl.gov/mpi/>`_, the Message-Passing
   Interface



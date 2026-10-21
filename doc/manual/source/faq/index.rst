Frequently Asked Questions
==========================

Building Enzo
-------------


**Q: I’m getting a compilation error related to HDF5. What is HDF5 and how to I get it?**

A: HDF5 is a data format with accompanying library for writing very large
data sets. Enzo uses HDF5 for data output. If you do not have a version of HDF5
available on your machine, you can download binaries or source code for HDF5
from https://www.hdfgroup.org/downloads/hdf5/. CMake will automatically locate HDF5 if installed in standard system locations or module paths. If HDF5 is installed in a custom directory (e.g. ``/home/enzo-user/local/hdf5/``), you can specify its path to CMake:
:: 

  $ cmake -B build -S . -DHDF5_ROOT=/home/enzo-user/local/hdf5
  $ cmake --build build

When running Enzo, make sure that the HDF5 library is in ``LD_LIBRARY_PATH``:
::

  $ export LD_LIBRARY_PATH=/home/enzo-user/local/hdf5/lib/:$LD_LIBRARY_PATH


Running Simulations
-------------------

Common Crashes
--------------


Misc.
-----


**Q: What is the difference between enzo-dev (week-of-code) and the stable
branch? Should I only use the stable branch?**

A:

The "week-of-code" branch of enzo-dev is the primary development branch, which
is updated on a fairly regular basis (the name "week-of-code" is historical).
Changes are migrated into the stable branch on a roughly annual basis. In
general, if you want code that is somewhat more reliable but may be
significantly behind the cutting-edge Enzo version, you should use the 'stable'
branch. If you are comfortable with more recent (and thus possibly less
reliable) code, you should use the "week-of-code" branch.



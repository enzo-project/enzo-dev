Frequently Asked Questions
==========================

This document lists a few frequently asked questions about Enzo.

**WORK IN PROGRESS: Please add questions new Enzo users may ask**

Building Enzo
-------------


**Q: I’m getting a compilation error that looks something like:**
::

 Grid_ComputeCoolingTime.C(201): error: expression must have class type
    if (grackle_data.use_grackle == TRUE)

**what does this mean?**


A: Check what branch you are on using with Enzo (``$ hg branch``). This issue is
likely related to a version conflict between Enzo and Grackle 2.0 and Grackle
3.0. Currently the stable branch only works with Grackle 2.0, and the
development branch (``$ hg checkout week-of-code``) works with Grackle 3.0. We
recommend switching to using the development branch and using the most recent
version of  Grackle 3.0 instead of using Grackle 2.0 if you need Grackle for
your simulations.


**Q: I’m getting a compilation error related to HDF5. What is HDF5 and how to I get it?**

A: HDF5 is a data format with accompanying library for writing very large
data sets. Enzo uses HFD5 for data output. If you do not have a version of HDF5
available on your machine, you can download binaries or source code for HDF5
from https://www.hdfgroup.org/downloads/hdf5/. Once you have a version of HDF5
installed on your machine, you need to notify Enzo where it is located for the
build process in the Makefile (eg. ``Make.mach.linux-gnu`` or
``Make.mach.my-machine``). For example, if HDF5 was installed in
``/home/enzo-user/local/hdf5/``, you would edit the line
::

  LOCAL_HDF5_INSTALL = /home/enzo-user/local/hdf5

then run
:: 

  $ make machine-linux-gnu
  $ make clean
  $ make

to rebuild enzo.exe with your HDF5 installation. When running enzo.exe, make
sure that the HDF5 library is in ``LD_LIBRARY_PATH``. In this example, if you
are running bash, run the command
::

  $ export LD_LIBRARY_PATH=/home/enzo-user/local/hdf5/lib/:$LD_LIBRARY_PATH 

to put the HDF5 library in the library path before running Enzo.


Running Simulations
-------------------

Common Crashes
--------------


Misc.
-----


**Q: What is the difference between enzo-dev (week-of-code) and the stable
branch? Should I only use the stable branch?**

A:



.. _obtaining_and_building_enzo:

Obtaining and Building Enzo
===========================


.. _CompilationRequirements:

Enzo Compilation Requirements
-----------------------------

Enzo can be compiled on any POSIX-compatible operating system, such as Linux,
BSD (including Mac OS X), and AIX. In addition to a C/C++ and Fortran-90
compiler, the following libraries and tools are necessary:

   * `CMake <https://cmake.org/>`_ (version 3.20 or newer), the cross-platform build system.
   * `HDF5 <http://www.hdfgroup.org/HDF5/>`_, the hierarchical data format.
     Note that HDF5 also may require the szip and zlib libraries, which can be
     found at the HDF5 website.
   * `MPI <http://www.mcs.anl.gov/research/projects/mpi/>`_, for multi-processor parallel
     jobs. Note that Enzo will compile without MPI if desired, but it is recommended to compile
     with MPI.
   * `git <https://git-scm.org/>`__, a free, distributed source control management tool.
   * `yt <http://yt-project.org>`_, the yt visualization and analysis suite.
     While it is not required to run enzo, ``yt`` enables the easiest analysis
     of its outputs, as well as the ability to run the enzo testing tools.
   * `libyt <https://libyt.readthedocs.io/en/latest/>`__, a C library for in situ analysis using Python and ``yt``.
     This is optional, see the details in :ref:`in_situ_python_analysis`.
 

Downloading Enzo
----------------

We encourage anyone who uses Enzo to sign up for the `Enzo Users'
List <http://groups.google.com/group/enzo-users>`_, where one can ask questions
to the community of enzo users and developers.  

Please visit the `Enzo Project home page <http://enzo-project.org>`_ to learn
more about the code and different installation methods. To directly access the source
code, you can visit the `Enzo Github page <https://github.com/enzo-project>`_.

If you already have Fortran, C, C++ compilers, CMake, MPI, and HDF5 installed,
then installation of Enzo should be straightforward. Simply run the following at the command line 
to get the latest stable version of the Enzo source using git:

.. highlight:: none

::

    ~ $ git clone https://github.com/enzo-project/enzo-dev.git

Enzo development continues regularly, and if you wish to use the
latest changes, you can update the code as follows:

::

    ~/enzo $ git fetch
    ~/enzo $ git merge
    Already up to date.

This covers the basics, but for more information about interacting with the
git version control system please peruse the :ref:`developers_guide`,
and any git tutorial such as `this one <https://git-scm.com/docs/gittutorial>`_.


Building Enzo
-------------

Enzo uses CMake to configure and build all targets across platforms. A comprehensive list of the CMake build options can be found in :ref:`CMakeOptions`.

.. note::
   If you are compiling a historical version of Enzo that uses the legacy ``Make.mach.*`` Makefile system, please consult the :ref:`LegacyMakeOptions` reference.

Quick Start: Building with Default Options
++++++++++++++++++++++++++++++++++++++++++

To build Enzo with default options (64-bit precision, MPI enabled, HDF5 enabled):

::

    ~ $ cd enzo-dev/
    ~/enzo-dev $ make

This wrapper ``Makefile`` automatically invokes CMake with parallel build jobs matching your available CPU cores.

Alternatively, you can call CMake directly from the top-level directory:

::

    ~/enzo-dev $ cmake -B build -S .
    ~/enzo-dev $ cmake --build build

All output executables will be generated in the ``bin/`` directory:

::

    ~/enzo-dev $ ls bin/
    enzo  inits  enzohop  ring  anyl  P-GroupFinder

Configuring Custom Build Options
++++++++++++++++++++++++++++++++

You can customize compilation options (such as precision, external libraries, or physics features) by passing ``-DOPTION=VALUE`` flags during CMake configuration.

For example, to build with 32-bit floating point precision and enable Grackle chemistry:

::

    ~/enzo-dev $ cmake -B build -S . -DENZO_PRECISION=32 -DENZO_USE_GRACKLE=ON
    ~/enzo-dev $ cmake --build build

Common Configuration Options
++++++++++++++++++++++++++++

Below are some frequently used CMake options:

- ``-DENZO_PRECISION=32|64``: Floating-point precision (default: 64)
- ``-DENZO_INTEGERS=32|64``: Grid integer size (default: 64)
- ``-DENZO_PARTICLES=32|64|128``: Particle position precision (default: 64)
- ``-DENZO_USE_MPI=ON|OFF``: Enable/disable MPI parallel support (default: ON)
- ``-DENZO_USE_GRACKLE=ON|OFF``: Enable/disable Grackle cooling library (default: OFF)
- ``-DENZO_USE_HYPRE=ON|OFF``: Enable/disable HYPRE linear solvers (default: OFF)

For a complete reference of all available options, see :ref:`CMakeOptions`.

Building Individual Tools
+++++++++++++++++++++++++

To compile a specific tool without building the entire suite, pass the target name to ``make`` or ``cmake --build``:

::

    ~/enzo-dev $ make enzo        # Builds bin/enzo only
    ~/enzo-dev $ make inits       # Builds bin/inits only
    ~/enzo-dev $ make enzohop     # Builds bin/enzohop only
    ~/enzo-dev $ make anyl        # Builds bin/anyl only

Cleaning the Build
++++++++++++++++++

To remove built binaries and object files:

::

    ~/enzo-dev $ make clean

Or remove the ``build/`` directory:

::

    ~/enzo-dev $ rm -rf build/

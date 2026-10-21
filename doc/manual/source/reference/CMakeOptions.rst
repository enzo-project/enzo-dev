.. _CMakeOptions:

The ``Enzo`` CMake Build System Reference
==========================================

``Enzo`` uses `CMake <https://cmake.org/>`_ (version 3.20 or newer) as its build system. CMake automatically locates system libraries (MPI, HDF5, Grackle, HYPRE, libuuid) and generates native Makefiles or build scripts for your platform.

Quickstart
----------

To configure and build Enzo with default settings (64-bit precision, MPI enabled, HDF5 enabled):

.. highlight:: none

::

       cd enzo/
       cmake -B build -S .
       cmake --build build

Or simply use the top-level wrapper ``Makefile``:

::

       cd enzo/
       make

The compiled executables will be placed in ``bin/`` (e.g. ``bin/enzo``, ``bin/inits``, ``bin/enzohop``, ``bin/ring``, ``bin/anyl``, ``bin/P-GroupFinder``).

Configuring Build Options
-------------------------

Build options are passed to CMake using ``-DOPTION=VALUE`` during the configuration step:

::

       cmake -B build -S . -DENZO_PRECISION=64 -DENZO_USE_GRACKLE=ON

Or interactively using ``ccmake build`` or ``cmake-gui build``.

Precision & Bit-Width Settings
------------------------------

======================  ===========  ==============================================================
Option                  Default      Description
======================  ===========  ==============================================================
**ENZO_PRECISION**      ``64``       Floating-point precision for grid calculations (32 or 64).
**ENZO_INTEGERS**       ``64``       Integer size for grid indices (32 or 64).
**ENZO_PARTICLES**      ``64``       Particle position precision (32, 64, or 128).
**ENZO_PARTICLE_IDS**   ``64``       Particle ID integer size (32 or 64).
**ENZO_INITS**          ``64``       Initial conditions precision (32 or 64).
**ENZO_IO**             ``32``       HDF5 I/O precision (32 or 64).
======================  ===========  ==============================================================

External Feature & Library Options
----------------------------------

======================  ===========  ==============================================================
Option                  Default      Description
======================  ===========  ==============================================================
**ENZO_USE_MPI**        ``ON``       Enable MPI parallel support.
**ENZO_USE_HDF4**       ``OFF``      Use HDF4 instead of HDF5.
**ENZO_USE_GRACKLE**    ``OFF``      Enable Grackle chemistry and cooling library integration.
**ENZO_USE_HYPRE**      ``OFF``      Enable HYPRE linear solver library (for implicit FLD RT).
**ENZO_USE_CUDA**       ``OFF``      Enable CUDA GPU acceleration.
**ENZO_USE_UUID**       ``ON``       Enable libuuid unique run identifier generation.
======================  ===========  ==============================================================

Physics & Algorithmic Options
-----------------------------

======================  ===========  ==============================================================
Option                  Default      Description
======================  ===========  ==============================================================
**ENZO_TRANSFER**       ``ON``       Enable radiative transfer (adaptive ray-tracing / ray-casting).
**ENZO_FAST_SIB**       ``ON``       Enable fast sibling grid search algorithm.
**ENZO_SAB**            ``ON``       Enable SAB (Subgrid Aggregation for Boundaries).
======================  ===========  ==============================================================

Free Parameters
---------------

===========================  ===========  ==============================================================
Option                       Default      Description
===========================  ===========  ==============================================================
**ENZO_MAX_SUBGRIDS**        ``100000``   Maximum number of subgrids.
**ENZO_MAX_BARYONS**         ``30``       Maximum number of baryon fields.
**ENZO_MAX_TASKS_PER_NODE**  ``8``        Maximum parallel tasks per node.
**ENZO_MEMORY_POOL_SIZE**    ``100000``   Initial memory pool size (in photon packages).
===========================  ===========  ==============================================================

Building Specific Targets
-------------------------

You can build individual targets using ``cmake --build build --target <target_name>`` or ``make <target_name>``:

- ``make enzo``: Main Enzo simulation executable
- ``make inits``: Initial conditions generator
- ``make enzohop``: Hop halo finder
- ``make ring``: Ring analysis tool
- ``make anyl``: Analysis package
- ``make P-GroupFinder``: Parallel group finder

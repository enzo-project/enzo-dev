.. _LegacyMakeOptions:

The Legacy ``Enzo`` Makefile System (Historical Archive)
=========================================================

.. note::
   This page documents the historical ``Makefile`` build system (``Make.mach.*``, ``./configure``, etc.) used prior to the transition to CMake. It is preserved here for historical reference when working with legacy Enzo releases. For current build instructions, please see :ref:`obtaining_and_building_enzo` and :ref:`CMakeOptions`.

The makefile system in ``Enzo`` was organized into separate files summarized below. Note that the files discussed on this page were found in the ``src/enzo`` subdirectory in legacy releases.

==================  ============
**Makefile**        The main makefile for compiling the ``Enzo`` executable ``enzo.exe``
**Make.mach.\***    These files contained all machine-dependent settings
**Make.config.\***  These files contained all compile-time configuration settings
==================  ============

For example, to compile ``Enzo`` on NICS's Kraken platform (starting from the top-level ``Enzo`` directory in legacy versions):

.. highlight:: none

::

       ./configure
       cd src/enzo
       gmake machine-nics-kraken
       gmake

Machine settings
----------------

These machine-specific configuration files were named ``Make.mach.machinename``.

General variables:

================ ============
**MACH_FILE**    Name of the make include file for the machine, e.g. ``Make.mach.nics-kraken``
**MACH_TEXT**    Description of the platform, e.g. ``"NICS Kraken"``
**MACH_VALID**   Should be set to 1, though not currently accessed
================ ============

Paths to compilers:

===================== ============
**MACH_CPP**          The C preprocessor
**MACH_CC_MPI**       The MPI C compiler
**MACH_CC_NOMPI**     The C compiler
**MACH_CXX_MPI**      The MPI C++ compiler
**MACH_CXX_NOMPI**    The C++ compiler
**MACH_F90_MPI**      The MPI F90 compiler
**MACH_F90_NOMPI**    The F90 compiler
**MACH_FC_MPI**       The MPI F77 compiler
**MACH_FC_NOMPI**     The F77 compiler
**MACH_CUDACOMPILER** The CUDA compiler
**MACH_LD_MPI**       The MPI linker (typically the MPI C++ compiler)
**MACH_LD_NOMPI**     The linker (typically the C++ compiler)
===================== ============

Compiler flags:

================== ============
**MACH_CPPFLAGS**  Machine-dependent flags for the C preprocessor, e.g.  ``-P -traditional``
**MACH_CFLAGS**    Machine-dependent flags for the C compiler
**MACH_CXXFLAGS**  Machine-dependent flags for the C++ compiler
**MACH_F90FLAGS**  Machine-dependent flags for the F90 compiler
**MACH_FFLAGS**    Machine-dependent flags for the F77 compiler
**MACH_LDFLAGS**   Machine-dependent flags for the linker
================== ============

Machine-specific flags:

============================== ============
**MACH_DEFINES**               Machine-specific defines, e.g. ``-DLINUX``, ``-DIBM``, ``-DIA64``, etc.
============================== ============

Paths to include header files:

========================= ============
**MACH_INCLUDES**         All required machine-dependent includes--should at least include HDF5.
**MACH_INCLUDES_HYPRE**   Includes for optional Hypre linear solver package
**MACH_INCLUDES_MPI**     Includes for MPI if needed
**MACH_INCLUDES_CUDA**    Includes for CUDA if needed
**MACH_INCLUDES_PYTHON**  Includes for Python if needed
========================= ============

Paths to library files:

====================== ============
**MACH_LIBS**          All required machine-dependent libraries--should at least include HDF5.
**MACH_LIBS_HYPRE**    Libraries for optional Hypre linear solver package
**MACH_LIBS_MPI**      Libraries for MPI if needed
**MACH_LIBS_PAPI**     Libraries for optional PAPI performance package (optionally called by ``lcaperf``)
**MACH_LIBS_CUDA**     Libraries for CUDA if needed
**MACH_LIBS_PYTHON**   Libraries for Python if needed
====================== ============

Optimization flags:

========================= ============
**MACH_OPT_AGGRESSIVE**   Compiler/link flags for "aggressive" optimization
**MACH_OPT_DEBUG**        Compiler/link flags for debugging
**MACH_OPT_HIGH**         Compiler/link flags for standard optimizations
**MACH_OPT_WARN**         Compiler/link flags to generate verbose warning messages
========================= ============

Makefile commands
-----------------

===============  ==============================================
**gmake**        Compile and generate the executable ``enzo.exe``
**gmake help**   Display help information
**gmake clean**  Remove object files, executable, etc.
===============  ==============================================

Configuration options
---------------------

Precision settings
~~~~~~~~~~~~~~~~~~

============================   =====================================
**integers-[32\|64]**          Set integer size to 32- or 64-bits.
**precision-[32\|64]**         Set floating-point precision to 32- or 64-bits.
**particles-[32\|64\|128]**    Set particle position precision to 32-, 64-, or 128-bits. 
**inits-[32\|64]**             Set inits precision to 32- or 64-bits.
**io-[32\|64]**                Set IO precision to 32- or 64-bits.
**particle-id-[32\|64]**       Set integer size for particle IDs
============================   =====================================

Algorithmic & Library settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

==============================  =====================================
**use-mpi-[yes\|no]**            Set whether to use MPI.
**fastsib-[no\|yes]**	         Include fast sibling search
**photon-[no\|yes]**	         Include radiative transfer (adaptive ray tracing)
**hypre-[no\|yes]**              Include HYPRE libraries (implicit RT solvers)
**cuda-[no\|yes]**               Set whether to use CUDA (GPU-computing)
**use-hdf4-[no\|yes]**           Set whether to use HDF4
==============================  =====================================

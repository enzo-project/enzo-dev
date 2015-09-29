.. _MakeOptions:

The ``Enzo`` Makefile System
============================

The makefile system in ``Enzo`` is a bit complicated, because it's
designed to work on many different platforms, allow many different
compile-time configuration settings, and be usable by automated
systems such as the ``lcatest`` parallel program testing
environment. 

To decouple machine-specific settings from configuration-specific
settings, it's organized into separate files summarized below.  Note
that the files discussed on this page are found in the ``src/enzo``
subdirectory.

==================  ============
**Makefile**        The main makefile for compiling the ``Enzo`` executable ``enzo.exe``
**Make.mach.\***    These files contain all machine-dependent settings
**Make.config.\***  These files contain all compile-time configuration settings
==================  ============

If there is already a ``Make.mach.*`` file present for the particular
machine you want to compile on, and you just want to compile ``Enzo``
with the default configuration, then compiling is relatively
straightforward. For example, to compile ``Enzo`` on NICS's Kraken
platform (starting from the top-level ``Enzo`` directory):

.. highlight:: none

::

       ./configure
       cd src/enzo
       gmake machine-nics-kraken
       gmake

If all goes well, this should create the ``enzo.exe`` executable in
the ``src/enzo`` subdirectory.  Also, note that ``gmake`` is required,
though ``make`` may work on your system as well.

Machine settings
----------------

If there is not already a ``Make.mach.*`` file present for your
platform, you will need to create one.  The easiest way to port
``Enzo`` to a new platform is to copy an existing ``Make.mach.*`` file
to a new one and edit it accordingly.  Generally, all variables
prefixed by ``MACH_`` in ``Make.mach.*`` files should be assigned a
value (even if that value is an empty string), and all variables that
begin with ``LOCAL_`` (or anything else) are optional and only
accessed within the ``Make.mach.*`` file itself.

The list of ``MACH_`` variables that can be set are listed below.

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
**MACH_INCLUDES**         All required machine-dependent includes--should at least include    HDF5.
**MACH_INCLUDES_HYPRE**   Includes for optional Hypre linear solver package
**MACH_INCLUDES_MPI**     Includes for MPI if needed
**MACH_INCLUDES_CUDA**    Includes for CUDA if needed
**MACH_INCLUDES_PYTHON**  Includes for Python if needed
========================= ============

Paths to library files:

====================== ============
**MACH_LIBS**          All required machine-dependent libraries--should at least include    HDF5.
**MACH_LIBS_HYPRE**    Libraries for optional Hypre linear solver package
**MACH_LIBS_MPI**      Libraries for MPI if needed
**MACH_LIBS_PAPI**     Libraries for optional PAPI performance package (optionally called    by ``lcaperf``)
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

Although it breaks from the ``MACH_*`` naming convention, there is
also a **MACHINE_NOTES** variable for machine-specific information
that is displayed whenever ``Enzo`` is compiled.



Makefile commands
-----------------

The default action of typing ``gmake`` without a target is to attempt
to compile ``Enzo``.  Other high-level makefile targets are ``help``,
and ``clean``:

===============  ==============================================
**gmake**        Compile and generate the executable ``enzo.exe``
**gmake help**   Display this help information
**gmake clean**  Remove object files, executable, etc.
===============  ==============================================

(For brevity we'll omit the ``gmake`` portion for the remainder of the
discussion.)

Configuration-related targets are ``help-config``, ``show-config``,
``show-flags``, and ``default``:

=================  ======================================================
**help-config**    Display detailed help on configuration make targets
**show-config**    Display the current configuration settings
**show-flags**     Display the current compilers and compilation flags
**default**        Reset the configuration to the default values
=================  ======================================================

Note that ``gmake default`` may also clear your machine setting, in
which case you will need to rerun gmake machine-*platform*.

Configuration options
---------------------


Other configuration targets, set using e.g. ``gmake integers-32``,
are listed below:

Free parameters
~~~~~~~~~~~~~~~

========================= ============
**max-subgrids-N**        Set the maximum number of subgrids to *N*.
**max-baryons-N**         Set the maximum number of baryon fields to *N*.
**max-tasks-per-node-N**  Set the number of tasks per node to *N*.
**memory-pool-N**         Set initial memory pool size (in number of photons).
========================= ============

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

Global settings
~~~~~~~~~~~~~~~

============================   =====================================
**object-mode-[32\|64]**       Set address/pointer size to 32-bit or 64-bit object files.  This is an    obsolete setting and is no longer used.
**testing-[yes\|no]**          Include hooks for the lcatest regression tests
============================   =====================================

Algorithmic settings
~~~~~~~~~~~~~~~~~~~~

========================   =====================================
**bitwise-[no\|yes]**      Turn on blocking-gravity for bitwise identical runs
**emissivity-[no\|yes]**   Include emissivity field
**fastsib-[no\|yes]**	   Include fast sibling search
**fluxfix-[no\|yes]**	   Include sibling subgrid boundary fix
**newgridio-[no\|yes]**	   Use the new Grid IO routines
**photon-[no\|yes]**	   Include radiative transfer (adaptive ray tracing)
========================   =====================================

External libraries
~~~~~~~~~~~~~~~~~~

===========================   =====================================
**use-mpi-[yes\|no]**         Set whether to use MPI.
**isolated-bcs-[yes\|no]**    Set whether to compile in isolated boundary conditions code
**tpvel-[yes\|no]**           Set whether to compile in tracer particle velocity information
**lcaperf-[yes\|no]**         Set whether to call the optional lcaperf performance tool
**papi-[yes\|no]**            Set whether to link in the PAPI library if required by lcaperf
**hypre-[no\|yes]**           Include HYPRE libraries (implicit RT solvers)
**cuda-[no\|yes]**            Set whether to use CUDA (GPU-computing)
**python-[no\|yes]**          Set whether to use inline python
**use-hdf4-[no\|yes]**        Set whether to use HDF4
===========================   =====================================

Performance settings
~~~~~~~~~~~~~~~~~~~~

================================= ============================
**opt-VALUE**                     Set optimization/debug/warning levels, where VALUE = [warn\|debug\|high\|aggressive\|cudadebug]
**taskmap-[yes\|no]**             Set whether to use unigrid taskmap performance modification
**packed-amr-[yes\|no]**          Set whether to use 'packed AMR' disk performance modification.
**packed-mem-[yes\|no]**          Set whether to use 'packed memory' option: requires packed AMR.
**unigrid-transpose-[yes\|no]**   Set whether to perform unigrid communication transpose performance   optimization
**ooc-boundary-[yes\|no]**        Set whether to use out-of-core handling of the boundary
**log2alloc-[yes\|no]**           Set whether to compile with grid/particle arrays allocated in sizes of powers of 2
================================= ============================


The ``Make.config.*`` Files
---------------------------

The ``Make.config.settings`` and ``Make.config.override`` files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default configuration settings and current configuration
settings are stored in the two files ``Make.config.settings`` and
``Make.config.override``.

The **Make.config.settings** file consists of assignments to the
``CONFIG_*`` make variables that define the default configuration
settings in ``Enzo``'s makefile. This file should not be modified
lightly.  If you type ``gmake default``, then these will become the
currently active settings.

The **Make.config.override** file, together with the
``Make.config.settings`` file, define the current configuration
settings. This file should also not be edited (since misspelled
configuration variable names may not be detected, leading to behavior
that is unexpected and difficult to locate), though it will be modified
indirectly through ``gmake`` when setting new configuration
values. For example, if you were to type ``gmake integers-32``, then
the ``Make.config.override`` file would contain ``CONFIG_INTEGERS =
32``.  The values in the ``Make.config.override`` file essentially
override the settings in ``Make.config.settings``.

In summary:

    default settings = **Make.config.settings**


    current settings =
    **Make.config.settings + Make.config.override**


Typing ``gmake default`` will clear the ``Make.config.override``
file entirely, making the default settings in ``Make.config.settings``
the current settings.

The ``Make.config.objects`` file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file is used simply to define the list of all object files,
excluding the file containing ``main()``. Only one variable needs to
be set.

======================  ==============
**OBJS_CONFIG_LIB**     List of all object files excluding the file containing ``main()``
======================  ==============

Dependencies are generated automatically using the makedepend
command and stored in the ``DEPEND`` file, so dependencies don't need
to be explicitly included.  If it complains about missing files,
such as ``DEPEND`` or ``Make.config.override``, then try (re)-running
the ``./configure`` script in the top-level ``Enzo`` subdirectory.

The ``Make.config.targets`` file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file contains rules for all configuration-related make
targets. It exists mainly to reduce the size of the top-level
Makefile. When adding new configuration settings, this file will
need to be modified.

The ``Make.config.assemble`` file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file contains all the makefile magic to convert configuration
settings (defined by ``$(CONFIG_*)`` make variables) into appropriate
compiler flags (such as ``$(DEFINES)``, ``$(INCLUDES)``, etc.). When
adding a new configuration setting, this file will need to be
modified.

James Bordner (jobordner at ucsd.edu)



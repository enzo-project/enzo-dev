The Enzo Makefile System
========================

The makefile in Enzo is a bit complicated, because it's designed to
work on many different platforms and allow many different
compile-time configuration settings. To decouple machine-specific
settings from configuration-specific settings, it's organized into
separate files as follows:

Makefile.config
    The main Makefile. Currently this is is sym-linked to Makefile by
    the configure script.
Make.mach.\*
    All machine-dependent settings are defined in these files.
Make.config.\*
    All compile-time configuration settings related to Enzo are in
    these files.

If you have a Make.mach.\* file for the particular machine you want
to compile on, and you just want to compile Enzo with the default
configuration, then compiling is straightforward. For example, to
compile Enzo on SDSC's `DataStar? </wiki/DataStar>`_ platform:

::

       ./configure
       cd src/enzo
       gmake machine-sdsc-datastar
       gmake

Note that gmake is required--standard make will generally not
work.

Makefile commands
-----------------

The default action of typing **gmake** without a target is to
attempt to compile Enzo.
Other high-level makefile targets are **install**, **help**, and
**clean**:

gmake
    Compile and generate the executable 'enzo.exe'
gmake install
    Copy the executable to bin/enzo
gmake help
    Display this help information
gmake clean
    Remove object files, executable, etc.

Configuration-related targets are **help-config**, **show-config**,
**show-flags**, and **default**:
gmake help-config
    Display detailed help on configuration make targets
gmake show-config
    Display the configuration settings
gmake show-flags
    Display specific compilation flags
gmake default
    Reset the configuration to the default values

Machine settings
----------------

All machine-dependent settings are defined in Make.mach.\* files,
and all settings in Make.mach.\* files should be machine-dependent.
The easiest way to port Enzo to a new platform is to copy an
existing Make.mach.\* file to a new one (preferably using the
Make.mach.ORG-MACHINE convention), and editing it accordingly.
Generally, all variables prefixed by MACH\_ in Make.mach.\* files
should be set (even if they are set to an empty string), and all
variables that begin with LOCAL\_ are optional and for convenience
only.

The list of MACH\_ variables that can be set are listed below:

MACH\_FILE
    Name of the makefile, e.g. *Make.mach.sdsc-datastar*
MACH\_TEXT
    Description of the platform, e.g. *SDSC Datastar*
MACH\_VALID
    Should be set to 1, though not currently accessed

MACH\_CPP
    The C preprocessor
MACH\_CC\_MPI
    The MPI C compiler
MACH\_CC\_NOMPI
    The C compiler
MACH\_CXX\_MPI
    The MPI C++ compiler
MACH\_CXX\_NOMPI
    The C++ compiler
MACH\_F90\_MPI
    The MPI F90 compiler
MACH\_F90\_NOMPI
    The F90 compiler
MACH\_FC\_MPI
    The MPI F77 compiler
MACH\_FC\_NOMPI
    The F77 compiler
MACH\_LD\_MPI
    The MPI linker (typically the MPI C++ compiler)
MACH\_LD\_NOMPI
    The linker (typically the C++ compiler)

MACH\_CPPFLAGS
    Machine-dependent flags for the C preprocessor, e.g.
    *-P -traditional*
MACH\_CFLAGS
    Machine-dependent flags for the C compiler
MACH\_CXXFLAGS
    Machine-dependent flags for the C++ compiler
MACH\_F90FLAGS
    Machine-dependent flags for the F90 compiler
MACH\_FFLAGS
    Machine-dependent flags for the F77 compiler
MACH\_LDFLAGS
    Machine-dependent flags for the linker

MACH\_DEFINES
    Machine-specific defines, e.g. *-DLINUX*, *-DIBM*, *-DIA64*, etc.
MACH\_FFLAGS\_INTEGER\_32
    Fortran flags for specifying 32-bit integers
MACH\_FFLAGS\_INTEGER\_64
    Fortran flags for specifying 64-bit integers
MACH\_FFLAGS\_REAL\_32
    Fortran flags for specifying 32-bit reals
MACH\_FFLAGS\_REAL\_64
    Fortran flags for specifying 64-bit reals

MACH\_INCLUDES
    All required machine-dependent includes--should at least include
    HDF5.
MACH\_INCLUDES\_HYPRE
    Includes for optional Hypre linear solver package
MACH\_INCLUDES\_JBPERF
    Includes for optional jbPerf (lcaperf) performance package
MACH\_INCLUDES\_MPI
    Includes for MPI if needed
MACH\_INCLUDES\_PAPI
    Includes for optional PAPI performance package (optionally called
    by jbPerf)

MACH\_LIBS
    All required machine-dependent libraries--should at least include
    HDF5.
MACH\_LIBS\_HYPRE
    Libraries for optional Hypre linear solver package
MACH\_LIBS\_JBPERF
    Libraries for optional jbPerf (lcaperf) performance package
MACH\_LIBS\_MPI
    Libraries for MPI if needed
MACH\_LIBS\_PAPI
    Libraries for optional PAPI performance package (optionally called
    by jbPerf)

MACH\_OPT\_AGGRESSIVE
    Compiler/link flags for "aggressive" optimization
MACH\_OPT\_DEBUG
    Compiler/link flags for debugging
MACH\_OPT\_HIGH
    Compiler/link flags for standard optimizations
MACH\_OPT\_WARN
    Compiler/link flags to generate verbose warning messages

Configuration options
---------------------

Use gmake help-config for online help about configuration settings.
Use gmake show-config for a summary of current settings in effect.
Use gmake default to set default settings (this may also clear your
machine setting, so you may need to rerun gmake machine-*platform*
to use settings in the corresponding Make.mach.*platform* machine
file.)

The configuration targets, set using e.g. gmake integers-32, are
listed below:

Free parameters
~~~~~~~~~~~~~~~

max-subgrids-*N*
    Set the maximum number of subgrids to *N*.
max-baryons-*N*
    Set the maximum number of baryon fields to *N*.
max-tasks-per-node-*N*
    Set the number of tasks per node to *N*.

Precision settings
~~~~~~~~~~~~~~~~~~

integers-[32\|64]
    Set integer size to 32- or 64-bits.
precision-[32\|64]
    Set floating-point precision to 32- or 64-bits.
particles-[32\|64\|128]
    Set particle position precision to 32-, 64-, or 128-bits. This
    should be 64.
inits-[32\|64]
    Set inits precision to 32- or 64-bits.
io-[32\|64]
    Set IO precision to 32- or 64-bits.

Global settings
~~~~~~~~~~~~~~~

object-mode-[32\|64]
    Set address/pointer size to 32-bit or 64-bit object files [NOT
    IMPLEMENTED]
testing-[yes\|no]
    Include hooks for the lcatest regression tests

External libraries
~~~~~~~~~~~~~~~~~~

use-mpi-[yes\|no]
    Set whether to use MPI. [REQUIRED FOR ENZO]
isolated-bcs-[yes\|no]
    Set whether to compile in isolated boundary conditions code
tpvel-[yes\|no]
    Set whether to compile in tracer particle velocity information
jbperf-[yes\|no]
    Set whether to call the optional jbPerf (lcaperf) performance tool
papi-[yes\|no]
    Set whether to link in the PAPI library if required by jbPerf

Performance settings
~~~~~~~~~~~~~~~~~~~~

opt-[warn\|debug\|high\|aggressive
    Set optimization/debug/warning levels
taskmap-[yes\|no]
    Set whether to use unigrid taskmap performance modification
packed-amr-[yes\|no]
    Set whether to use 'packed AMR' disk performance modification.
packed-mem-[yes\|no]
    Set whether to use 'packed memory' option: requires packed AMR.
unigrid-transpose-[yes\|no]
    Set whether to perform unigrid communication transpose performance
    optimization
ooc-boundary-[yes\|no]
    Set whether to use out-of-core handling of the boundary
load-balance-[yes\|no]
    Set whether to use load balancing of grids.

The Make\* Files
----------------

The Make.config.settings and Make.config.override files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default configuration settings and current configuration
settings are stored in the two files Make.config.settings and
Make.config.override.

The **Make.config.settings** file consists of assignments to the
CONFIG\_\* make variables that define
the default configuration settings in Enzo's makefile. Generally
this file should never be modified.
If you type "gmake default", then these will become the currently
active settings.

The **Make.config.override** file, together with the
Make.config.settings file, define the current configuration
settings. This file should also not be edited, though it may be
modified indirectly when setting new configuration settings. For
example, if you were to type "gmake integers-32", then the
Make.config.override file would contain "CONFIG\_INTEGERS = 32".
The values in the Make.config.override file essentially override
the settings in Make.config.settings.

In summary:

    default settings = **Make.config.settings**


    current settings =
    **Make.config.settings + Make.config.override**


Typing "gmake default" will clear the Make.config.override file
entirely, making the default settings in Make.config.settings the
current settings.

The Make.config.objects file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file is used simply to define the list of all object files,
excluding the file containing "main()". Only one variable needs to
be set.

OBJS\_CONFIG\_LIB
    List of all object files excluding the file containing "main()"

Dependencies are generated automatically using the makedepend
command and stored in the DEPEND file, so dependencies don't need
to be explicitly included.

The Make.config.targets file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file contains rules for all configuration-related make
targets. It exists mainly to reduce the size of the top-level
Makefile. When adding new configuration settings, this file will
need to be modified.

The Make.config.assemble file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file contains all the makefile magic to convert configuration
settings (defined by $(CONFIG\_\*) make variables) into appropriate
compiler flags (such as $(DEFINES), $(INCLUDES), etc.). When adding
a new configuration setting, this file will need to be modified.

James Bordner (jobordner at ucsd.edu)



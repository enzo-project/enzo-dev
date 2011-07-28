.. _obtaining_and_building_enzo:

Obtaining and Building Enzo
===========================


.. _CompilationRequirements:

Enzo Compilation Requirements
-----------------------------

Enzo can be compiled on any POSIX-compatible operating system, such as Linux,
BSD (including Mac OS X), and AIX.  In addition to a C/C++ and Fortran-90
compiler, the following libraries are necessary:

   * `HDF5 <http://www.hdfgroup.org/HDF5/>`_, the hierarchical data format.
     Note that HDF5 also may require the szip and zlib libraries, which can be
     found at the HDF5 website.  Note that compiling with HDF5 1.8 or greater
     requires that the compiler directive ``H5_USE_16_API`` be specified;
     typically this is done with ``-DH5_USE_16_API`` and it's set in most of
     the provided makefiles.
   * `MPI <http://www.mcs.anl.gov/research/projects/mpi/>`_, for multi-processor parallel
     jobs.  Note that Enzo will compile without MPI, but it's fine to compile
     with MPI and only run oon a single processor.

Mercurial Check Out Instructions
--------------------------------

Enzo is provided in both a stable and an unstable form.  **It is highly
recommended that for any production run the stable version is used.**
Additionally, we encourage anyone who uses Enzo to sign up for the `Enzo Users'
List <http://groups.google.com/group/enzo-users>`_.  A source
browser is also available.

Please visit the Google Code project website to access the Enzo source tree and
read the latest source checkout instructions.

http://enzo.googlecode.com/

Updating a source tree with Mercurial is beyond the scope of this document; for
more information, please peruse :ref:`developers_guide` and the Mercurial
documentation.  The `mercurial <http://mercurial.selenic.com/>`_ commands of
most use are ``pull``, ``update`` and ``incoming``.

Building Enzo
-------------

This is a quick, line by line example of checking out and building
Enzo using current build system. A comprehensive list of the make
system arguments can be found in :ref:`MakeOptions`.

This assumes that we're working from a checkout from the Enzo project page,
located at http://enzo.googlecode.com/ .  Checkout instructions can be found
there, and for more detailed information about the structure of the Enzo source
control repository, see :ref:`enzo_modification`.

Initializing the Build System
+++++++++++++++++++++++++++++

This just clears any existing configurations left over from a previous machine,
and creates a couple of files for building.

.. highlight:: none

::

    ~ $ cd enzo/
    ~/enzo $ ./configure 

This should output a brief message saying that the build system has been
initialized.  To confirm that it ran, there should be a file called
Make.config.machine in the src/enzo subdirectory.

Go to the Source Directory
++++++++++++++++++++++++++

The source code for the various Enzo components are laid out in the
src/ directory.

::

    ~/enzo/src $ cd src/
    ~/enzo/src $ ls
    Makefile      P-GroupFinder anyl          enzo          enzohop       inits
    lcaperf       mpgrafic      ring
    ~/enzo/src $ 

Right now, we're just building the main executable (the one that
does the simulations), so we need the ``enzo/`` directory.

::

    ~/enzo/src $ cd enzo/

Find the Right Machine File
+++++++++++++++++++++++++++

We've chosen to go with configurations files based on specific
machines. This means we can provide configurations files for most
of the major NSF resources, and examples for many of the one-off
(clusters, laptops, etc.).

These machine-specific configuration files are named ``Make.mach.machinename``.

::

    ~/enzo/src/enzo $ ls Make.mach.*
    Make.mach.darwin          Make.mach.nasa-discover   Make.mach.ncsa-cobalt
    Make.mach.ornl-jaguar-pgi Make.mach.tacc-ranger     Make.mach.unknown
    Make.mach.kolob           Make.mach.nasa-pleiades   Make.mach.nics-kraken
    Make.mach.scinet          Make.mach.triton
    Make.mach.linux-gnu       Make.mach.ncsa-abe        Make.mach.orange
    Make.mach.sunnyvale       Make.mach.triton-intel
    ~/enzo/src/enzo $ 

In this example, we choose ``Make.mach.darwin``, which is appropriate for Mac
OS X machines.

Porting
+++++++

If there's no machine file for the machine you're on, you will have
to do a small amount of porting. However, we have attempted to
provide a wide base of Makefiles, so you should be able to find one
that is close, if not identical, to the machine you are attempting
to run Enzo on. The basic steps are as follows:


#. Find a Make.mach file from a similar platform.
#. Copy it to Make.mach.site-machinename (site = sdsc or owner,
   machinename = hostname).
#. Edit the machine-specific settings (compilers, libraries, etc.).
#. Build and test.

If you expect that you will have multiple checkouts of the Enzo source code,
you should feel free to create the directory $HOME/.enzo/ and place your custom
makefiles there, and Enzo's build system will use any machine name-matching
Makefile in that directory to provide or override Make settings.

Make sure you save your configuration file! If you're on a big system (multiple
Enzo users), please post your file to `the Enzo mailing list
<http://groups.google.com/group/enzo-users>`_, and it will be
considered for inclusion with the base Enzo distribution.

HDF5 Versions
+++++++++++++

If your system uses a version of HDF5 greater than or equal to 1.8, you
probably need to add a flag to your compile settings, unless your HDF5 library
was compiled using --with-default-api-version=v16. The simplest thing to do is
to find the line in your Make.mach file that sets up MACH_DEFINES, which may
look like this

::

    MACH_DEFINES   = -DLINUX # Defines for the architecture; e.g. -DSUN, -DLINUX, etc.

and change it to

::

    MACH_DEFINES   = -DLINUX -DH5_USE_16_API # Defines for the architecture; e.g. -DSUN, -DLINUX, etc.

This will ensure that the HDF5 header files expose the correct API
for Enzo.

Build the Makefile
++++++++++++++++++

Now that you have your configuration file, tell the build system to
use it:

::

    ~/enzo/src/enzo $ make machine-darwin
    
     *** Execute 'gmake clean' before rebuilding executables ***
    
       MACHINE: Darwin (OSX Leopard)
    
    ~/enzo/src/enzo $ 

You may also to know the settings (precision, etc.) that's being
use. You can find this out using ``make show-config``. For a detailed
explanation of what these mean, see :ref:`MakeOptions`.

::

    ~/enzo/src/enzo $ make show-config
    
    MACHINE: Darwin (OSX Leopard)
    MACHINE-NAME: darwin
    
    PARAMETER_MAX_SUBGRIDS:       100000
    PARAMETER_MAX_BARYONS:        20
    PARAMETER_MAX_TASKS_PER_NODE: 8
    PARAMETER_MEMORY_POOL_SIZE:   100000
    
    CONFIG_PRECISION:             64
    CONFIG_PARTICLES:             64
    CONFIG_INTEGERS:              64
    CONFIG_PARTICLE_IDS:          64
    CONFIG_INITS:                 64
    CONFIG_IO:                    32
    CONFIG_USE_MPI:               yes
    CONFIG_OBJECT_MODE:           64
    CONFIG_TASKMAP:               no
    CONFIG_PACKED_AMR:            yes
    CONFIG_PACKED_MEM:            no
    CONFIG_LCAPERF:               no
    CONFIG_PAPI:                  no
    CONFIG_PYTHON:                no
    CONFIG_ECUDA:                 no
    CONFIG_OOC_BOUNDARY:          no
    CONFIG_OPT:                   debug
    CONFIG_TESTING:               no
    CONFIG_TPVEL:                 no
    CONFIG_PHOTON:                yes
    CONFIG_HYPRE:                 no
    CONFIG_EMISSIVITY:            no
    CONFIG_USE_HDF4:              no
    CONFIG_NEW_GRID_IO:           yes
    CONFIG_BITWISE_IDENTICALITY:  yes
    CONFIG_FAST_SIB:              yes
    CONFIG_FLUX_FIX:              yes
    
    ~/enzo/src/enzo $ 

Build Enzo
++++++++++

The default build target is the main executable, Enzo.

::

    ~/enzo/src/enzo $ make
    Updating DEPEND
    pdating DEPEND
    Compiling enzo.C
    Compiling acml_st1.src
    ... (skipping) ...
    Compiling Zeus_zTransport.C
    Linking
    Success!
    ~/enzo/src/enzo $ 

After compiling, you will have ``enzo.exe`` in the current directory.

Building other Tools
++++++++++++++++++++

Building other tools is typically very straightforward; they rely on the same
Makefiles, and so should require no porting or modifications to configuration.

Inits
~~~~~

::

    ~/enzo/src/ring $ cd ../inits/
    ~/enzo/src/inits $ make
    Compiling enzo_module.src90
    Updating DEPEND
    Compiling acml_st1.src
    ...
    Compiling XChunk_WriteIntField.C
    Linking
    Success!

This will produce ``inits.exe``.

Ring
~~~~

::

    ~/enzo/src/enzo $ cd ../ring/
    ~/enzo/src/ring $ make
    Updating DEPEND
    Compiling Ring_Decomp.C
    Compiling Enzo_Dims_create.C
    Compiling Mpich_V1_Dims_create.c
    Linking
    Success!

This will produce ``ring.exe``.

.. _build_yt:

YT
~~

To install yt, you can use the installation script provided with the yt source
distribution.  See `the yt homepage <http://yt.enzotools.org/>`_ for more
information.

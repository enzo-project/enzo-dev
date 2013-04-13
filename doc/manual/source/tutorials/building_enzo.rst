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
     with MPI and only run on a single processor.
   * `yt <http://yt-project.org>`_, the yt visualization and analysis suite.  
     While it is not required to run enzo, ``yt`` enables the easiest analysis
     of its outputs, as well as the ability to run the enzo testing tools.  It
     also provides an easy way to download enzo as part of its installation script.
     See the `Enzo Project home page <http://enzo-project.org/>`_ for more 
     information.

Downloading Enzo
----------------

Enzo is provided in both a stable and an unstable (developer's) form.  **It is highly
recommended that for any production run the stable version is used.**
Additionally, we encourage anyone who uses Enzo to sign up for the `Enzo Users'
List <http://groups.google.com/group/enzo-users>`_, where one can ask questions
to the community of enzo users and developers.  

Please visit the `Enzo Project home page <http://enzo-project.org>`_ to learn
more about the code and different installation methods.  To directly access the source
code, you can visit the `Enzo Bitbucket page <https://bitbucket.org/enzo>`_.

If you already have Fortran, C, C++ compilers, 
`Mercurial <http://mercurial.selenic.com>`_, 
`MPI <http://www.mcs.anl.gov/research/projects/mpi/>`_, and 
`HDF5 <http://www.hdfgroup.org/HDF5/>`_ installed, then installation of
Enzo should be straightforward.  Simply run the following at the command line 
to get the latest stable version of the Enzo source using Mercurial. This 
command makes a copy of the existing enzo source code repository on your local 
computer in the current directory:

.. highlight:: none

::

    ~ $ hg clone https://bitbucket.org/enzo/enzo-stable

Or if you're feeling fiesty, you can get the developer's version which is 
more up to date, but may also have some experimental code mixed in (which 
might affect stability):

.. highlight:: none

::

    ~ $ hg clone https://bitbucket.org/enzo/enzo-dev

Later on, if you want to update your code and get any additional modifications 
which may have occurred since you originally cloned the source repository, 
you will have to ``pull`` them from the server and then ``update`` your 
local copy (in this example, no new changes have occurred):

.. highlight:: none

::

    ~ $ cd enzo-stable
    ~/enzo-stable $ hg pull
    pulling from https://bitbucket.org/enzo/enzo-stable
    searching for changes
    no changes found

    ~/enzo-stable $ hg update
    0 files updated, 0 files merged, 0 files removed, 0 files unresolved

    ~/enzo-stable $ 

This covers the basics, but for more information about interacting with the
mercurial version control system please peruse the :ref:`developers_guide`,
the `Mercurial Documentation <http://mercurial.selenic.com/>`_, and/or 
this entertaining `tutorial on Mercurial <http://hginit.com>`_.

Building Enzo
-------------

This is a quick, line by line example for building
Enzo using the current build system. A comprehensive list of the make
system arguments can be found in :ref:`MakeOptions`.

This assumes that we're working from a checkout (or download) of the source
after following instructions on the `Enzo Project home page <http://enzo-project.org>`_, or the instructions in the last section.  For more detailed information 
about the structure of the Enzo source control repository, see 
:ref:`enzo_modification`.

Initializing the Build System
+++++++++++++++++++++++++++++

This just clears any existing configurations left over from a previous machine,
and creates a couple of files for building.

::

    ~ $ cd enzo-stable/
    ~/enzo-stable $ ./configure 
    Configure complete.

    ~/enzo-stable $ 

This message just confirms that the build system has been
initialized.  To further confirm that it ran, there should be a file called
Make.config.machine in the src/enzo subdirectory.

Go to the Source Directory
++++++++++++++++++++++++++

The source code for the various Enzo components are laid out in the
src/ directory.

::

    ~/enzo-stable $ cd src
    ~/enzo-stable/src $ ls
    Makefile      P-GroupFinder      TREECOOL      anyl      enzo      enzohop
    inits         lcaperf            mpgrafic      performance_tools   ring

    ~/enzo-stable/src $ 

Right now, we're just building the main executable (the one that
does the simulations), so we need the ``enzo/`` directory.

::

    ~/enzo-stable/src $ cd enzo/

Find the Right Machine File
+++++++++++++++++++++++++++

We've chosen to go with configurations files based on specific
machines. This means we can provide configurations files for most
of the major NSF resources, and examples for many of the one-off
(clusters, laptops, etc.).

These machine-specific configuration files are named ``Make.mach.machinename``.

::

    ~/enzo-stable/src/enzo $ ls Make.mach.*
    Make.mach.arizona               Make.mach.darwin                
    Make.mach.hotfoot-condor        Make.mach.kolob                 
    Make.mach.linux-gnu             Make.mach.nasa-discover         
    Make.mach.nasa-pleiades         Make.mach.ncsa-bluedrop         
    Make.mach.ncsa-bluewaters-gnu   Make.mach.ncsa-cobalt           
    Make.mach.nics-kraken           Make.mach.nics-kraken-gnu       
    Make.mach.nics-kraken-gnu-yt    Make.mach.nics-nautilus
    Make.mach.orange                Make.mach.ornl-jaguar-pgi
    Make.mach.scinet                Make.mach.sunnyvale
    Make.mach.tacc-ranger           Make.mach.trestles
    Make.mach.triton                Make.mach.triton-gnu
    Make.mach.triton-intel          Make.mach.unknown

    ~/enzo-stable/src/enzo $ 

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
use it (remember to ``make clean`` if you change any previous settings):

::

    ~/enzo-stable/src/enzo $ make machine-darwin
    
     *** Execute 'gmake clean' before rebuilding executables ***
    
       MACHINE: Darwin (OSX Leopard)
    
    ~/enzo-stable/src/enzo $ 

You may also want to know the settings (precision, etc.) that are being
use. You can find this out using ``make show-config``. For a detailed
explanation of what these mean, see :ref:`MakeOptions`.

::

    ~/enzo-stable/src/enzo $ make show-config
    
    MACHINE: Darwin (OSX Leopard)
    MACHINE-NAME: darwin

    PARAMETER_MAX_SUBGRIDS  [max-subgrids-###]                : 100000
    PARAMETER_MAX_BARYONS  [max-baryons-###]                  : 30
    PARAMETER_MAX_TASKS_PER_NODE  [max-tasks-per-node-###]    : 8
    PARAMETER_MEMORY_POOL_SIZE  [memory-pool-###]             : 100000
 
    CONFIG_PRECISION  [precision-{32,64}]                     : 64
    CONFIG_PARTICLES  [particles-{32,64,128}]                 : 64
    CONFIG_INTEGERS  [integers-{32,64}]                       : 64
    CONFIG_PARTICLE_IDS  [particle-id-{32,64}]                : 64
    CONFIG_INITS  [inits-{32,64}]                             : 64
    CONFIG_IO  [io-{32,64}]                                   : 32
    CONFIG_USE_MPI  [use-mpi-{yes,no}]                        : yes
    CONFIG_OBJECT_MODE  [object-mode-{32,64}]                 : 64
    CONFIG_TASKMAP  [taskmap-{yes,no}]                        : no
    CONFIG_PACKED_AMR  [packed-amr-{yes,no}]                  : yes
    CONFIG_PACKED_MEM  [packed-mem-{yes,no}]                  : no
    CONFIG_LCAPERF  [lcaperf-{yes,no}]                        : no
    CONFIG_PAPI  [papi-{yes,no}]                              : no
    CONFIG_PYTHON  [python-{yes,no}]                          : no
    CONFIG_NEW_PROBLEM_TYPES  [new-problem-types-{yes,no}]    : no
    CONFIG_ECUDA  [cuda-{yes,no}]                             : no
    CONFIG_OOC_BOUNDARY  [ooc-boundary-{yes,no}]              : no
    CONFIG_ACCELERATION_BOUNDARY  [acceleration-boundary-{yes,no}]    : yes
    CONFIG_OPT  [opt-{warn,debug,cudadebug,high,aggressive}]  : debug
    CONFIG_TESTING  [testing-{yes,no}]                        : no
    CONFIG_TPVEL  [tpvel-{yes,no}]]                           : no
    CONFIG_PHOTON  [photon-{yes,no}]                          : yes
    CONFIG_HYPRE  [hypre-{yes,no}]                            : no
    CONFIG_EMISSIVITY  [emissivity-{yes,no}]                  : no
    CONFIG_USE_HDF4  [use-hdf4-{yes,no}]                      : no
    CONFIG_NEW_GRID_IO  [newgridio-{yes,no}]                  : yes
    CONFIG_BITWISE_IDENTICALITY  [bitwise-{yes,no}]           : no
    CONFIG_FAST_SIB  [fastsib-{yes,no}]                       : yes
    CONFIG_FLUX_FIX  [fluxfix-{yes,no}]                       : yes
    CONFIG_GRAVITY_4S  [gravity-4s-{yes,no}]                  : no
    CONFIG_ENZO_PERFORMANCE  [enzo-performance-{yes,no}]      : yes
    
    ~/enzo-stable/src/enzo $ 

Build Enzo
++++++++++

The default build target is the main executable, Enzo.

::

    ~/enzo-stable/src/enzo $ make
    Updating DEPEND
    pdating DEPEND
    Compiling enzo.C
    Compiling acml_st1.src
    ... (skipping) ...
    Compiling Zeus_zTransport.C
    Linking
    Success!

    ~/enzo-stable/src/enzo $ 

After compiling, you will have ``enzo.exe`` in the current directory.
If you have a failure during the compiler process, you may get enough of
an error message to track down what was responsible.  If there is a failure
during linking, examine the ``compile.out`` file to learn more about 
what caused the problem.  A common problem is that you forgot to include the 
current location of the HDF5 libraries in your machine-specific makefile.

Building other Tools
++++++++++++++++++++

Building other tools is typically very straightforward; they rely on the same
Makefiles, and so should require no porting or modifications to configuration.

Inits
~~~~~

::

    ~/enzo-stable/src/ring $ cd ../inits/
    ~/enzo-stable/src/inits $ make
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

    ~/enzo-stable/src/enzo $ cd ../ring/
    ~/enzo-stable/src/ring $ make
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

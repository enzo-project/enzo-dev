Building enzo
=============

`TOC? </wiki/TOC>`_

This is a quick, line by line example of checking out and building
Enzo using current build system. A comprehensive list of the make
system arguments can be found on the
`MakeOptions? </wiki/Devel/UserGuide/MakeOptions>`_ page.

This assumes that we're working from the public version, either the
repository or the release. To get a copy, either see the
`CheckOutInstructions? </wiki/Devel/UserGuide/CheckOutInstructions>`_,
or `DownloadInstructions? </wiki/DownloadInstructions>`_ page.

Checkout the Code
-----------------

For the impatient, you can get a copy by checking out the latest
public version.

::

    ~ $ svn co http://mngrid.ucsd.edu/svn/Enzo/public/trunk enzo

Initializing the Build System
-----------------------------

This just clears any existing configurations left over from a
previous machine,
and creates a couple of files for building.

::

    ~ $ cd enzo/
    ~/enzo $ ./configure 

Don't be worried if it doesn't output anything--it's not supposed
to. To confirm that it ran, there should be a file called
Make.config.machine in the src/enzo subdirectory.

Go to the Source Directory
--------------------------

The source code for the various Enzo components are laid out in the
src/ directory.

::

    ~/enzo/src $ cd src/
    ~/enzo/src $ ls
    Makefile enzo     enzohop  inits    ring     yt
    ~/enzo/src $ 

Right now, we're just building the main executable (the one that
does the simulations), so we need the enzo/ directory.

::

    ~/enzo/src $ cd enzo/

Find the Right Machine File
---------------------------

We've chosen to go with configurations files based on specific
machines. This means we can provide configurations files for most
of the major NSF resources, and examples for many of the one-off
(clusters, laptops, etc.).

These machine-specific configuration files are named:
Make.mach.machinename. For this example, I'm working on my laptop,
which is still suitable for small test problems.

::

    ~/enzo/src/enzo $ ls Make.mach.*
    Make.mach.bordner-krummhorn   Make.mach.ornl-jaguar-pgi
    Make.mach.bwoshea-fnord       Make.mach.padoan-cluster
    Make.mach.bwoshea-thunderhead Make.mach.psc-bigben
    Make.mach.gso-mac             Make.mach.rpwagner-cable
    Make.mach.ncsa-abe            Make.mach.sdsc-datastar
    Make.mach.ncsa-cobalt         Make.mach.sdsc-teragrid
    Make.mach.nics-kraken         Make.mach.unknown
    Make.mach.ornl-jaguar-gnu
    ~/enzo/src/enzo $ 

As you can see, we already have a makefile:
Make.mach.rpwagner-cable. Additional example Makefiles reside in
doc/example\_makefiles.

Porting
-------

If there's no machine file for the machine you're on, you will have
to do a small amount of porting. However, we have attempted to
provide a wipe base of Makefiles, so you should be able to find one
that is close, if not identical, to the machine you are attempting
to run Enzo on. The basic steps are as follows:


#. Find a Make.mach file from a similar platform.
#. Copy it to Make.mach.site-machinename (site = sdsc or owner,
   machinename = hostname).
#. Edit the machine-specific settings (compilers, libraries, etc.).
#. Build and test.

If you expect that you will have multiple checkouts of the enzo
source code, you should feel free to create the directory
$HOME/.enzo/ and place your custom makefiles there, as Enzo's build
system will use any machine name-matching Makefile in that
directory to provide or override Make settings.

Make sure you save your configuration file! If you're on a big
system (multiple Enzo users), please post your file to
` the Enzo mailing list <http://mailman.ucsd.edu/mailman/listinfo/enzo-users-l>`_,
and it will be considered for inclusion with the base Enzo
distribution.

HDF5 Versions
-------------

If your system uses a version of HDF5 greater than or equal to 1.8,
you probably need to add a flag to your compile settings, unless
your HDF5 library was compiled using
--with-default-api-version=v16. The simplest thing to do is to find
the line in your Make.mach file that sets up MACH\_DEFINES, which
may look like this

::

    MACH_DEFINES   = -DLINUX # Defines for the architecture; e.g. -DSUN, -DLINUX, etc.

and change it to

::

    MACH_DEFINES   = -DLINUX -DH5_USE_16_API # Defines for the architecture; e.g. -DSUN, -DLINUX, etc.

This will ensure that the HDF5 header files expose the correct API
for Enzo.

Build the Makefile
------------------

Now that you have your configuration file, tell the build system to
use it:

::

    ~/enzo/src/enzo $ make machine-rpwagner-cable
    
     *** Execute 'gmake clean' before rebuilding executables ***
    
       MACHINE: Rick's Laptop (Make.mach.rpwagner-cable)
    
    ~/enzo/src/enzo $ 

You may also to know the settings (precision, etc.) that's being
use. You can find this out using make show-config. For a detailed
explanation of what these mean, head over to the
`MakeOptions? </wiki/Devel/UserGuide/MakeOptions>`_ page.

::

    ~/enzo/src/enzo $ make show-config
    
       MACHINE: Rick's Laptop (Make.mach.rpwagner-cable)
    
       PARAMETER_MAX_SUBGRIDS:       100000
       PARAMETER_MAX_BARYONS:        20
       PARAMETER_MAX_TASKS_PER_NODE: 8
    
       CONFIG_PRECISION:             64
       CONFIG_PARTICLES:             64
       CONFIG_INTEGERS:              64
       CONFIG_INITS:                 64
       CONFIG_IO:                    32
       CONFIG_USE_MPI:               yes
       CONFIG_OBJECT_MODE:           64
       CONFIG_TASKMAP:               no
       CONFIG_PACKED_AMR:            yes
       CONFIG_PACKED_MEM:            no
       CONFIG_JBPERF:                no
       CONFIG_PAPI:                  no
       CONFIG_UNIGRID_TRANSPOSE:     yes
       CONFIG_OOC_BOUNDARY:          no
       CONFIG_OPT:                   debug
       CONFIG_TESTING:               no
       CONFIG_ISOBCS:                no
       CONFIG_TPVEL:                 no
    
    ~/enzo/src/enzo $ 

Build Enzo
----------

The default build target is the main executable, enzo.

::

    ~/enzo/src/enzo $ make
    awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_config(FILE *fp) {"}; {print "   fprintf (fp,\""$0"\\n\");"}; END {print "}"}' < temp.show-config > auto_show_config.C
    awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_flags(FILE *fp) {"}; {print "   fprintf (fp,\""$0"\\n\");"}; END {print "}"}' < temp.show-flags > auto_show_flags.C
    awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_version(FILE *fp) {"}; {print "   fprintf (fp,\""$0"\\n\");"}; END {print "}"}' < temp.show-version > auto_show_version.C
    Updating DEPEND
    pdating DEPEND
    Compiling enzo.C
    Compiling acml_st1.src
    ...
    Compiling Zeus_zTransport.C
    Linking
    Success!
    ~/enzo/src/enzo $ 

After compiling, you can have the build system copy the executable
to a bin/ directory at the top level.

::

    ~/enzo/src/enzo $ make install 
    if [ ! -e ../../bin ]; then mkdir ../../bin; fi
    make -s show-flags   >& ../../bin/enzo.show-flags
    make -s show-config  >& ../../bin/enzo.show-config
    make -s show-version >& ../../bin/enzo.show-version
    make -s show-diff    >& ../../bin/enzo.show-diff
    ~/enzo/src/enzo $

Now that you've got things build, maybe you'll want to check out
some
`Tutorials on running simulations? </wiki/Tutorials#ControllingEnzoSimulations>`_.

Building other Tools
--------------------

Here's the quick steps to building ring, inits and
[http;*yt.enzotools.org/ yt].*

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
    ~/enzo/src/inits $ make install 
    if [ ! -e ../../bin ]; then mkdir ../../bin; fi
    make show-flags   >& ../../bin/inits.show-flags
    make show-config  >& ../../bin/inits.show-config
    make show-version >& ../../bin/inits.show-version
    ~/enzo/src/inits $

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
    ~/enzo/src/ring $ make install 
    if [ ! -e ../../bin ]; then mkdir ../../bin; fi
    make show-flags   >& ../../bin/ring.show-flags
    make show-config  >& ../../bin/ring.show-config
    make show-version >& ../../bin/ring.show-version

YT
~~

YT comes with an installer script which will run from within the
Enzo source distribution, obtaining all needed dependencies. If you
are comfortable with installing software, you should feel free to
follow the standard installation instructions (which are available
in the ` documentation <http://yt.enzotools.org/doc/>`_ or
` wiki <http://yt.enzotools.org/wiki/InstallationInstructions>`_)
but otherwise you should be able to set the variable DEST\_DIR
inside the installer script and execute it to have it handle all of
those steps for you. (If you're going to run on OSX, the
instructions are
` slightly different <http://yt.enzotools.org/wiki/OSXInstallation>`_
and the installation script is not recommended.)

::

    ~/enzo/src/yt $ nano doc/install_script.sh  # Or your favorite editor!
    ~/enzo/src/yt $ bash doc/install_script.sh

If you run into problems with linking or compilation, common
solutions include requesting YT to install zlib and HDF5. Please
also feel free to post requests for help with installation or usage
on the yt-users
` mailing list <http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org>`_!



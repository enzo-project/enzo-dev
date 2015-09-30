Programming Guide
=================

There are several coding practices that we should adhere to when
programing for Enzo. Some are style, some are more important for
the health of the code (and other Enzo users' projects/sanity).

Remember that other programmers will read your code
---------------------------------------------------

    "Everyone knows that debugging is twice as hard as writing a
    program in the first place. So if you're as clever as you can be
    when you write it, how will you ever debug it?"
    --Brian Kernighan "The Elements of Programming Style", 2nd edition,
    chapter 2


Comments
--------

At the very least, put the following in things at the top of each
of your functions:


-  Your name
-  The date you wrote it. If you modified it in a significant way,
   note the date and modification.
-  Effects on global or class member variables.
-  Variable names that are not obvious. As a rule of thumb, the
   name is not obvious.
-  Primary references where applicable.

Two more rules:


-  Write your comments now. You will not have time to come back and
   clean it up later.
-  If you change something, change the comments. Now. Wrong
   comments are worse than no comments.

float is double
---------------

One must constantly be wary of the possibility of built in C types
to be re-defined to higher precision types. This is outlined
in :ref:`FloatIsDouble`.

Allocating arrays with ``new``
------------------------------

Enzo has a significant issue with memory fragmentation.  This is due
to grid and particle storage arrays being constantly created,
destroyed, and re-created with slightly different sizes.  The most
successful solution to this problem has been to limit the sizes of
arrays to powers of 2.  When compiling Enzo with log2alloc-yes,
``new`` is overloaded such that arrays are always created with sizes
that are the nearest power of 2.  Thus, it is important to keep in
mind that arrays created with ``new`` will usually be slightly
larger than you think.

Fortran types
-------------

Unlike Enzo's C and C++ routines, Fortran files (.F and .F90) do not
re-define the built-in 'integer' and 'real' types, so all variables
and constants must be defined with the appropriate precision.  There
are pre-defined type specifiers that will match Enzo's C and C++
precision re-definitions, which should be used for all variables that
pass through the C/Fortran interface.  This is discussed in detail in 
:ref:`FloatIsDouble`.

Header Files
------------

Header files must be included in the correct order. This is due, among other
things, to the redefinition of float which is done in
``macros_and_parameters.h``. This must be done before Enzo headers, but after
external libraries. The order should be as follows:

.. code-block:: c

    #include "ErrorExceptions.h"
    #include "svn_version.def"
    #include "EnzoTiming.h"
    #include "performance.h"
    #include "macros_and_parameters.h"
    #include "typedefs.h"
    #include "global_data.h"
    #include "units.h"
    #include "flowdefs.h"
    #include "Fluxes.h"
    #include "GridList.h"
    #include "ExternalBoundary.h"
    #include "Grid.h"
    #include "Hierarchy.h"
    #include "LevelHierarchy.h"
    #include "TopGridData.h"
    #include "communication.h"
    #include "CommunicationUtilities.h"

Accessing BaryonField
---------------------

Access data in the BaryonField array as is described in the page on
:ref:`BaryonFieldAccess`.

Accessing the Hierarchy
-----------------------

The hierarchy should be traversed as described in
:ref:`LinkedLists`.

enum
----

The enum construct in C has no standardized size, which can cause
problems when using 64 bit integers. Direct integer assignment
should be used instead. This also has the added advantage of making
explicit the values of parameters that are also used in parameter
files. The typical idiom should be:

.. code-block:: c

    #ifdef SMALL_INTS
    typedef int hydro_method;
    #endif
    #ifdef LARGE_INTS
    typedef long_int hydro_method;
    #endif
    const hydro_method
      PPM_DirectEuler      = 0,
      PPM_LagrangeRemap    = 1,
      Zeus_Hydro           = 2,
      HD_RK                = 3,
      MHD_RK               = 4,
      HydroMethodUndefined = 5;



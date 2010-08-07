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


File Naming Convention
----------------------

With very few exceptions, Enzo has a one function per file layout, with the
file name being the function name. Object methods have the object name
prepended to the beginning, such as the member of the grid class
``SolveHydroEquations`` lives in the file ``Grid_SolveHydroEquations.C``.

This does create a large number of files. Familiarity with ``grep`` or ``ack``
and pipes like ``ls -1 |grep`` are essential.

Internal capitalization is used for C files, all lowercase with underscores for
fortran files and header files. All Fortran files end with .F.

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

Header Files
------------

Header files must be included in the correct order. This is due, among other
things, to the redefinition of float which is done in
``macros_and_parameters.h``. This must be done before Enzo headers, but after
external libraries. The order should be as follows:

::

    #include "ErrorExceptions.h"
    #include "svn_version.def"
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

::

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



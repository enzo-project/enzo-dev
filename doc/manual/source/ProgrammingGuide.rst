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


`Funciton/File? </wiki/Funciton/File>`_ Names
---------------------------------------------

With very few exceptions, Enzo has a one function per file layout,
with the file name being the function name. Object methods have the
object name prepended to the beginning, such as the member of the
grid class SolveHydroEquations lives in the file
Grid\_SolveHydroEquations.

This does create a large number of files. Familiarity with grep or
ack and pipes like ls -1 \|grep are essential.

Internal capitalization is used for C files, all lowercase with
underscores for fortran files and header files. Due to the
preprocessing method for fortran used in Enzo, all fortran files
must be named .src or .src90

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
`in this page? </wiki/Tutorials/FloatIsDouble>`_

Header Files
------------

Header files must be included in the correct order. This is due,
among other things, to the redefinition of float which is done in
macros\_and\_parameters.h. This must be done before Enzo headers,
but after external libraries. The order should be as follows:

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
`Accessing the Baryon Field? </wiki/Tutorials/BaryonFieldAccess>`_.

Accessing the Hierarchy
-----------------------

The hierarchy should be traversed as described in
`this page on the linked lists? </wiki/Tutorials/LinkedLists>`_

enum
----

The enum construct in C has no standardized size, which can cause
problems when using 64 bit integers. Direct integer assignment
should be used instean
Here we plot the number of "cell/particle operations" per processor
per second, which is proportional to FLOPS per processor. This
quantity is roughly sum\_over\_level [ (number of timesteps) x
(number\_of\_cells + number\_of\_particles) ]. It is computed in
EvolveLevel as

::

    if (level == 0) CPOP = 0;
    while (TimeThisLevel < TimeLevelAbove) {
      for (grid = FirstGrid; grid != NULL; grid = grid->NextGrid)
         CPOP += grid->NumberOfCells + grid->NumberOfParticles;
      TimeThisLevel += dt;
    }

`Image(allmachines.png, 25%)? </wiki/Image(allmachines.png,%2025%)>`_

Because enzo parallelizes over grids, it is also useful to plot
this quantity versus the maximum number of grids on a given level.

`Image(allmachines\_grids.png, 25%)? </wiki/Image(allmachines_grids.png,%2025%)>`_



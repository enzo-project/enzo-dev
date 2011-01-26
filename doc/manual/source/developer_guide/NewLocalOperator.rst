Adding a new Local Operator.
============================

If you're adding new physics to Enzo, chances are you'll need to
add some kind of new operator.

This page is only to describe new physics that is


-  Operator split from everything else
-  Completely local, so depends on the field value and their
   derivatives in a cell.
-  Doesn't depend on the grid's position in the hierarchy.

Global operators, such as solution to Poisson's equation, are much
more significant undertakings, and should be discussed with the
Enzo development team.


#. `Read all the supporting documents found here.? </wiki/EnzoPrimaryReferences>`_
   This is not a simple piece of software.

It's really in your best interest to understand the basic
algorithms before trying to write code to extend it. It's much more
complex than Gadget or Zeus, and much much easier to break.


2. Open ``FastSib\_EvolveHierarchy.C``


3. Read it, and understand the structure. The flowcharts can help,
   they can be found HEY FIND THE LINK.


4. `Add a parameter to drive your code ? </wiki/Tutorials/AdddingNewParameters>`_


5. Write your new routine. This can either be a grid member
   function (old style) or a non-member function that accesses the
   enzo data using the `Field Array? </wiki/GridFieldArrays>`_ objects
   (prefered method.)


6. Locate this block of code:
   ::

             if (Grids[grid1]->GridData->SolveHydroEquations(LevelCycleCount[level],
                NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1], level) == FAIL) {
               fprintf(stderr, "Error in grid->SolveHydroEquations.\n");
               return FAIL;
             }
       
             JBPERF_STOP("evolve-level-13"); // SolveHydroEquations()
       
       //      fprintf(stderr, "%"ISYM": Called Hydro\n", MyProcessorNumber);
       
             /* Solve the cooling and species rate equations. */


This is in the primary grid loop on this level.


6. Insert your new grid operation right before the last comment. It
   should look something like this:
   ::

             if (Grids[grid1]->GridData->SolveHydroEquations(LevelCycleCount[level],
                NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1], level) == FAIL) {
               fprintf(stderr, "Error in grid->SolveHydroEquations.\n");
               return FAIL;
             }
       
             JBPERF_STOP("evolve-level-13"); // SolveHydroEquations()
       
       //      fprintf(stderr, "%"ISYM": Called Hydro\n", MyProcessorNumber);
       
             /* Solve the cooling and species rate equations. */
       
             if( YourFlag ){
               if( Grids[grid1]->GridData->YourRoutine(YourArguments) == FAIL ){
                 fprintf(stderr,"Error in grid->YourRoutine\n");
                 return FAIL;
               }


If your code isn't a grid member, you can omit the
``Grids[grid1]->GridData->`` part.



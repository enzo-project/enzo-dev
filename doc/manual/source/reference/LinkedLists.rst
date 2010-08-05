.. _LinkedLists:

Getting Around the Hierarchy: Linked Lists in Enzo
==================================================

There are two primary linked lists in Enzo; HierarchyEnetry and
LevelHierarchyEntry. They're both used to traverse the hierarchy,
but in very different ways. HierarchyEntry is used to traverse
**down** the hierarchy, from a parent to its children.
LevelHierarchyEntry is used to traverse **across** the hierarchy,
on a single level.

One of the primary things to note about the two lists is that
NextGridThisLevel (which exists in both) serve different purposes.

In LevelHierarchyEntry, NextGridThisLevel links all the grids on a
given level together.

In HierarchyEntry, NextGridThisLevel only counts things on a given
level that share a parent.

Below I will present a description of the structures and they're
creation and usage in Enzo.

HierarchyEntry
--------------

The HierarchyEntry linked list is used for traversing *down* the
hierarchy, from parents to children.

This is the contents of the definition of the structure, which you
can find in [source:/public/trunk/src/enzo/Hierarchy.h
HierarchyEntry.h]

::

    struct HierarchyEntry
    {
      HierarchyEntry *NextGridThisLevel; /* pointer to the next grid on level */
      HierarchyEntry *NextGridNextLevel; /* pointer to first child of this grid */
      HierarchyEntry *ParentGrid;        /* pointer to this grid's parent */
      grid           *GridData;          /* pointer to this grid's data */
    };

NextGridThisLevel connects all children of a parent.
NextGridNextLevel points to the first child of the given grid.
ParentGrid connects to the parent, and GridData points to the
actual grid structure.

Usage of HierarchyEntry lists
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The HierarchyEntry list is used (among other things) whenever
communication between
child and parent grids needs to be done. The typical pattern for
looping over all the children of a parent grid is as following:

::

    1.)    HierarchyEntry * NextGrid = ParentGrid->NextGridNextLevel;
    2.)    while (NextGrid != NULL ){
    3.)      if (NextGrid->GridData->SomeFunctionOnChildren(args) == FAIL )
    4.)        fprintf(stderr, "Error in your function\n");
    5.)        return FAIL;
    6.)      }
    7.)      NextGrid = NextGrid->NextGridThisLevel;
    8.)    }

Line 1 sets the pointer NextGrid to the "first" child of the parent
grid.

Line 2 starts the while loop.

Lines 3-6 is the standard function call pattern in Enzo.

Line 7 advances the pointer to the next child on the *child*
level.

This loop stops once all the children of ParentGrid have been
accessed, because the last child grid
of a given parent has NULL as NextGridThisLevel.

Generation of HierarchyEntry lists
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The HierarchyEntry linked list is generated in several different
points in the code. The details are slightly different for each
place it's used, depending on the details of what that linked list
is used for and the assumed structure of the hierarchy at that
point. The list most used in the code is the one generated in
[source:public/trunk/src/enzo/FindSubgrids.C FindSubgrids.C],
called in [source:public/trunk/src/enzo/RebuildHierarchy.C
RebuildHierarchy.C]. This code is called on a single 'Parent Grid'
at a time. Paraphrased and annotated:

::

     1   HierarchyEntry *, *ThisGrid;
     2   PreviousGrid = &ParentGrid;
     3    for (i = 0; i < NumberOfSubgrids; i++) {
     4
     5      ThisGrid = new HierarchyEntry;
     6
     7      if (PreviousGrid == &ParentGrid)
     8        ParentGrid.NextGridNextLevel = ThisGrid;
     9      else
    10        PreviousGrid->NextGridThisLevel = ThisGrid;
    11      ThisGrid->NextGridNextLevel = NULL;
    12      ThisGrid->NextGridThisLevel = NULL;
    13      ThisGrid->ParentGrid        = &ParentGrid;
    14
    15      ThisGrid->GridData = new grid;
    16      ThisGrid->GridData = Setup Functions Skipped for clarity;
    17
    18      PreviousGrid = ThisGrid;
    19   }

Line 1 starts the HierarchyEntry list with ParentGrid. (Called
simply Grid in the source, changed here for clarity.)

Line 5 creates the next HierarchyEntry to be added to the list.

Line 7-8 attaches the new subgrid, and the ensuing subgrid chain,
to the parent grid (note that this is only done for the first new
subgrid)

line 10 attaches all subsequent new subgrids to the
NextGridThisLevel chain.

Lines 11 and 12 ensure that both lists terminate with this new
grid. NextGridThisLevel will be replaced if there is in fact a next
grid. Since this routine is called only on a single Parent at a
time, one can now see that for HierarchyEntry, the
NextGridThisLevel list only links children that belong to the same
Parent Grid.

Lines 13-17 finish setting up this grid.

If you're writing a new problem generator, and have been brought
here by the AMR problem generation page, we advise that you examine
one of the other code patterns that are used in Enzo. They look
fairly similar to the above code, though have some details
different. Some suggestions are:

For adding a single subgrid, visit
[source:public/trunk/src/enzo/SphericalInfallInitialize.C
SphericalInfallInitialize.C]

For adding a single stack of nested subgrids, see
[source:public/trunk/src/enzo/ProtostellarCollapseInitialize.C
ProtostellarCollapseInitialize.C ]

For a completely general, though more complex setup, see
[source:public/trunk/src/enzo/CosmologySimulationInitialize.C
CosmologySimulationInitialize.C]

Another notable routine that generates HierarchyEntry lists is
[source:public/trunk/src/enzo/CommunicationPartitionGrid.C], which
breaks the TopGrid pointer across multiple processors.

LevelHierarchyEntry and LevelArray
----------------------------------

The LevelHierarchyEntry Linked List is used for traversing all the
grids on a given level. It's a simpler structure than
HierarchyEntry. The source can be found in
[source:/public/trunk/src/enzo/LevelHierarchy.h LevelHierarchy.h]

::

    struct LevelHierarchyEntry
    {
      LevelHierarchyEntry *NextGridThisLevel;  /* next entry on this level */
      grid                *GridData;           /* pointer to this entry's grid */
      HierarchyEntry      *GridHierarchyEntry; /* pointer into hierarchy */
    };

NextGridThisLevel connects all grids on a given level. GridData
points to the actual grid object, and GridHierarchyEntry points to
the (unique) HierarchyEntry node discussed above.

The LevelHierarchyEntry lists, one for each populated level, are
all bundled together in the LevelArray object. Both data structures
will be discussed presently.

Usage of LevelHierarchyEntry and LevelArray
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main usage of the LevelHierarchyEntry list is quite similar to
the main loop for HierarchyEntry lists.

::

      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
        if (Temp->GridData->MyCode(MyArgs) == FAIL) {
          fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
          return FAIL;
        }
        Temp = Temp->NextGridThisLevel;
      }

This calls MyCode for each grid on level.

Generation of LevelHierarchyEntry and LevelArray
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is done in two places in the code: in
[source:/public/trunk/src/enzo/main.C main.C] and
[source:/public/trunk/src/enzo/RebuildHierarchy.C
RebuildHierarchy.C]. It's done by the code
[source:/public/trunk/src/enzo/LevelHierarchy\_AddLevel.C
LevelHierarchy\_AddLevel.C], which I'll spend a moment with.

The setup:

::

    Prep in main.C:
      for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
        LevelArray[level] = NULL;

The call in main:

::

    AddLevel(LevelArray, &TopGrid, 0);

The fill:

::

     1   void AddLevel(LevelHierarchyEntry *LevelArray[], HierarchyEntry *Grid,
     2                 int level)
     3   {
     4      LevelHierarchyEntry *ThisLevel;
     5
     6     /* create a new LevelHierarchyEntry for the HierarchyEntry Grid                                          
     7        and insert it into the head of the linked list (LevelArray[level]). */
     8
     9     ThisLevel = new LevelHierarchyEntry;
    10     ThisLevel->GridData = Grid->GridData;
    11     ThisLevel->NextGridThisLevel = LevelArray[level];
    12     ThisLevel->GridHierarchyEntry = Grid;
    13     LevelArray[level] = ThisLevel;
    14
    15     /* recursively call this for the next grid on this level. */
    16
    17     if (Grid->NextGridThisLevel != NULL)
    18       AddLevel(LevelArray, Grid->NextGridThisLevel, level);
    19  
    20     /* ... and then descend the tree. */
    21
    22     if (Grid->NextGridNextLevel != NULL)
    23       AddLevel(LevelArray, Grid->NextGridNextLevel, level+1);
    }

This is a recursive function that takes LevelArray that's to be
filled, the HierarchyEntry list that fills it, and a counter for
the level. It's recursive in both HierarchyEntry's lists, both
NextGridNextLevel and NextGridThisLevel. The most notable lines are
11, 13, and 17. In lines 11 and 13, one can see that the current
HierarchyEntry is attached to the HEAD of the list, but line 17
shows that the HierarchyEntry list is traversed from its head to
its tail: so the LevelArray list is backwards from the
HierarchyEntry. This is only really needed information on the top
grid.

Traversing the Entire Hierarchy
-------------------------------

Sometimes the user needs to traverse the entire hierarchy. This is
done with a recursive function call on the HierarchyEntry. This
should be done in a manner akin to the AddLevel code above.



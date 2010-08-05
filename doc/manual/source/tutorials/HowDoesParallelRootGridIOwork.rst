.. _ParallelRootGridIO:

ParallelRootGridIO
==================

ParallelRootGridIO is a set of Enzo behaviors that allow the user
to run problems who's root grids are larger than the available
memory on a single node.

This page is intended for developers that need to write new problem
generators that will be run at extremely large scale. Large problem
size will need to utilize the ParallelRootGridIO machinery in Enzo.
As this brings a significant amount of added complexity, it isn't
recommended for smaller problems. It is also recommended that you
write the problem generator without this machinery first, and test
on smaller problems, before adding the additional complexity. If
you don't intend to write your own problem generator, this page is
basically irrelevant.

Background: why it is how it is
-------------------------------

Parallel Root Grid IO (PRGIO) is an essential component of doing
any simulations at large scale. In its initial inception, Enzo
worked on shared memory machines. This meant that the total
computer memory available dictated the problem size. Enzo would
allocate the root grid on the root processor, then distribute it to
the other processors. When it came time to write the data, the root
grid was collected back to the root processor, and written in a
single file.

This worked fine until the distributed machine showed up at the
party. This also came with a growth of the desired root grid size
for the enzo simulation. Now, the total aggregate memory of the
computer and the maximum allocatable array were vastly different.
So the old model broke down, because you simply can't fit the
15x512\ :sup:`3`\  arrays you need on a 512 Mb of ram, but you can
on 64 nodes. So out of necessity, PRGIO was born.

Short version
-------------

Essentially, PRGIO has three components (though not called in this
order)


-  `Input/Restart? </wiki/Input/Restart>`_
-  Output
-  Initialization

`Input/Restart? </wiki/Input/Restart>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is pretty straight forward: each grid is read by the processor
that wants it, out of the hdf5 file containing it.

Output
~~~~~~

This is pretty straight forward: it simply doesn't do the
collection back to the root processor of the root grid tiles. Each
processor writes an hdf5 file for each grid it owns. In the packed
AMR version, each processor writes one file, and in it go all the
grids it owns.

Initialization
~~~~~~~~~~~~~~

This is the part that needs attention, because the details are not
obvious from the code itself.

Initialization **BEFORE** PRGIO happened in three steps:


-  Set up grid
-  Allocate Data on the TopGrid object, on the Root Processor
-  Partition TopGrid across processors.

**WITH** PRGIO, the order is different:


-  Set up grid
-  Partition TopGrid
-  Allocate Data on the *working* grids.

Setup and Allocation
~~~~~~~~~~~~~~~~~~~~

This is pretty straight forward in principle, but the
implementation is a little confusing.

First, grids need to be set up. There aren't very many things you
need to do. See
`Part 1 of Making a New Test Problem.? </wiki/NewTestProblem/Part1_SerialUnigrid>`_
(which currently doesn't exist.) Basically, it needs to count the
NumberOfBaryonFields and record which one goes where in the
FieldType array.

After the Partition (next section), you need to allocate the data.

The confusing bits are in the implementation. We'll describe this
by way of example, using Cosmology simulations as our descriptor.
CosmologySimulationInitialize.C contains two routines:
CosmologySimulationInitialize (from here out CSI) and
CosmologySimulation*Re*Initialize (CSRI). These are both called in
InitializeNew. The job of the first routine is to set up the
hierarchy of grids and subgrids you'll need for your cosmology
simulation, and call CosmologySimulationInitialize*Grid* (CSIG).
Both CSI and CSIG are called whether or not PRGIO is on. CSRI is
called from InitializeNew after the Top Grid is partitioned. It is
*only* called when PRGIO is on.

Stated a different way:

::

    InitializeNew: Reads the parameter file, then calls
        CosmologySimulationInitialize: sets up the grid hierarchy.  On each of those grids gets called
            CosmologySimulationInitializeGrid: which sets NumberOfBaryonFields, and may allocate data.
        PartitionGrid: breaks the root grid into parts, and sends those parts to the other processors.
        CosmologySimulationReInitialize: If PRGIO is on, this is called.  It loops over grids and calls
            CosmologySimulationInitializeGrid: again, which allocates and defines the data.      

CSI passes a flag, TotalRefinement to CSIG for each grid you
initialize. This is equal (refinement
factor)\ :sup:`(refinement level of this grid.)`\  So for the Top
grid, this is equal to 1, and something that isn't 1 on all other
grids.

Inside of CSIG: if PRGIO is on **and** TotalRefinement == 1, then
statements relating to reading data from disk, allocating memory,
and accessing memory are **skipped.** (this is done by setting
ReadData = FALSE) In all other cases, it's left on. (So if PRGIO is
off, or **this grid** is not on the root level.) Thus at the first
pass at initialization, the TopGrid doesn't get it's BaryonFields
allocated.

Author's Note: It's not entirely obvious to me right now how this
works on Nested Sub Grids: It looks to me like they get allocated
on the root processor regardless.

CSRI is called AFTER the root grid has been partitioned and sent
off to the other processors. It does very little except call CSIG
again. This time when CSIG is called, TotalRefinement = -1. This
allows the data to be allocated.

Author's Second Note: Near as I can tell, TotalRefinement serves no
purpose other than to force the root grid to not be allocated
before the TopGrid is Partitioned.

Partition TopGrid and /\* bad kludge \*/
----------------------------------------

The other confusing part the partition, specifically a line in
ExternalBoundary::Prepare.

::

    if (ParallelRootGridIO == TRUE)
        TopGrid->NumberOfBaryonFields = 0; /* bad kludge! */

More on that in a moment.

CommunicationPartitionGrid is the routine that takes the TopGrid
(or, any grid) and breaks it across the processors. It first sorts
out the layout of the processors with MPI\_Dims\_create. It then
evenly splits the initial grid over those processors by first
creating a new grid on each tile, linking them to the Hierarchy
linked list. It then (and here's the part that caught me off guard)
allocates each grid on the Root processor and copies data from the
Initial Grid to the new tile. Finally, it take these freshly minted
root grid tiles and sends them to their new processor home.

Here's where the **bad kludge!** comes in. You'll note that in the
above description, there's an allocate on each of the newly minted
tiles *on the root processor*, which will allocate more than the
root grid data. This is the problem we were trying to avoid. So
ExternalBoundary::Prepare is nice enough to set
NumberOfBaryonFields to zero, so when the allocate comes around
it's allocating Zero fields.

Why is it in ExternalBoundary::Prepare? A look at the lines
immediately preceding the 'kludge' help:

::

      BoundaryRank = TopGrid->GridRank;
      NumberOfBaryonFields = TopGrid->NumberOfBaryonFields;
      if (ParallelRootGridIO == TRUE)
        TopGrid->NumberOfBaryonFields = 0; /* bad kludge! */

In order to do its job properly, the ExternalBoundary objects need
to know how many BaryonFields there are in the simulation. So
ExternalBoundary::Prepare records the data, and because that's the
last place NumberOfBaryonFields is needed, sets it to zero.

So now, when CommunicationPartitionGrid gets to the point where it
allocates the data, NumberOfBaryonFields is now zero, so it
allocates no data. These empty root grid tiles are then distributed
to the other processors.

Finally, CosmologyReInitialize is called, which calls
CosmologyInitializeGrid. This code then resets NumberOfBaryonFields
to its proper value, and since TotalRefinement = -1 allocates all
the data.

Then the simulation continues on, only aware of PRGIO when it comes
time to not collect the data again.



.. _ParallelRootGridIO:

Parallel Root Grid IO
=====================

Parallel Root Grid IO (PRGIO) is a set of Enzo behaviors that allow
the user to run problems that has a root grid larger than the
available memory on a single node.

This page is intended for developers that need to write new problem
generators that will be run at extremely large scale. Large problem
size will need to utilize the PRGIO machinery in Enzo.  As this brings
a significant amount of added complexity, it isn't recommended for
smaller problems. It is also recommended that you write the problem
generator without this machinery first, and test on smaller problems,
before adding the additional complexity. If you don't intend to write
your own problem generator, this page is basically irrelevant.

Background: why it is how it is
-------------------------------

PRGIO is an essential component of doing any simulations at large
scale. In its initial inception, Enzo worked on shared memory
machines. This meant that the total computer memory available dictated
the problem size. Enzo would allocate the root grid on the root
processor, then distribute spatially decomposed parts of the root grid
to the other processors. When it came time to write the data, the root
grid was collected back to the root processor, and written in a single
file.

This worked fine until distributed computers were deployed in response
to the limitations of a shared memory computer.  This coincided with a
growth of the desired root grid size for the Enzo simulation. Now, the
total aggregate memory of a single shared memory computer and the
memory required were vastly different.  The old model broke down
because you simply can't fit the 15x512\ :sup:`3`\ arrays you need in
512 Mb of RAM, but you can on 64 nodes if the memory is taken as an
aggregate total.  So out of necessity, PRGIO was born.

Short version
-------------

Essentially, PRGIO has three components (though not called in this
order)


-  :doc:`Input/Restart </user_guide/RunningEnzo>`
-  :ref:`PRGIO_Output`
-  :ref:`PRGIO_Initialization`

Input and Restarting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During initialization, the root grid is partitioned into tiles, and
each processor reads the part, i.e. a HDF5 hyperslab, of the initial
data files.  For restarts, each grid is read by one processor that
owns the data (``ProcessorNumber == MyProcessorNumber``) from the HDF5
file containing it.

.. _PRGIO_Output:

Output
~~~~~~

Unlike early versions of Enzo that collected all the grid data on one
processor before writing to disk, with PRGIO each processor writes an
HDF5 file for each grid it owns. In the packed AMR output mode, each
processor writes one HDF5 file, and in it go all the grids it owns.

.. _PRGIO_Initialization:

Initialization
~~~~~~~~~~~~~~

This is the part that needs attention, because the details are not
obvious from the code itself.

Initialization **BEFORE** PRGIO happens in three steps:


-  Set up grid
-  Allocate Data on the ``TopGrid`` object, on the Root Processor
-  Partition ``TopGrid`` across processors.

**WITH** PRGIO, the order is different:


-  Set up grid
-  Partition ``TopGrid``
-  Allocate Data on the *working* grids.

Setup and Allocation
~~~~~~~~~~~~~~~~~~~~

This is pretty straightforward in principle, but the
implementation is a little confusing.

First grids need to be set up. There aren't very many things you need
to do. See :ref:`InitializeGrid` for a more
comprehensive overview.  Simplified, a count of the
``NumberOfBaryonFields`` is made and a record of which field is which
goes in the ``FieldType`` array.

After the Partition (next section), you need to allocate the data.

The confusing bits are in the implementation. We'll describe this by
way of example, using Cosmology simulations as our descriptor.
``CosmologySimulationInitialize.C`` contains two routines:
``CosmologySimulationInitialize()`` (CSI) and
``CosmologySimulationReInitialize()`` (CSRI). These are both called in
``InitializeNew()``. The job of the first routine is to set up the
hierarchy of grids and subgrids you'll need for your cosmology
simulation, and call ``CosmologySimulationInitializeGrid`` (CSIG).
Both CSI and CSIG are called whether or not PRGIO is on. CSRI is
called from ``InitializeNew()`` after the Top Grid is partitioned. It
is *only* called when PRGIO is on.

Stated a different way:

#. InitializeNew: reads the parameter file, then calls
#. CosmologySimulationInitialize: sets up the grid hierarchy.  On each of those grids gets called
#. CosmologySimulationInitializeGrid: which sets NumberOfBaryonFields, and may allocate data.
#. PartitionGrid: breaks the root grid into parts, and sends those parts to the other processors.
#. CosmologySimulationReInitialize: If PRGIO is on, this is called. It loops over grids and calls CosmologySimulationInitializeGrid again, which allocates and defines the data.

CSI passes a flag, ``TotalRefinement`` to CSIG for each grid you
initialize. This is equal to (refinement factor)\ :sup:`(refinement
level of this grid)`. So for the Top grid, this is equal to 1, and
something that is greater than 1 on all other grids.

Inside of CSIG: if PRGIO is on **and** ``TotalRefinement`` == 1, then
statements relating to reading data from disk, allocating memory,
and accessing memory are **skipped.** (this is done by setting
``ReadData = FALSE``) In all other cases, it's left on. (So if PRGIO is
off, or **this grid** is not on the root level.) Thus at the first
pass at initialization, the ``TopGrid`` doesn't get it's ``BaryonFields``
allocated.

The same procedure is done on the nested initial grids if
``PartitionNestedGrids`` == 1.  If not, the root processor will read
the entire nested grid, partition it into smaller subgrids, and
finally send the data to different processors if ``LoadBalancing >
0``.  Regardless of the value of ``PartitionNestedGrids``, the
partitions of the static nested grids will never be re-combined for
I/O, unlike the behavior of the root grid when PRGIO is off.

CSRI is called AFTER the root grid has been partitioned and sent
off to the other processors. It does very little except call CSIG
again. This time when CSIG is called, ``TotalRefinement = -1``. This
allows the data to be allocated.

Partition TopGrid and /\* bad kludge \*/
----------------------------------------

The other confusing part the partition, specifically a line in
``ExternalBoundary::Prepare()``.

::

    if (ParallelRootGridIO == TRUE)
        TopGrid->NumberOfBaryonFields = 0; /* bad kludge! */

More on that in a moment.

``CommunicationPartitionGrid()`` is the routine that takes the ``TopGrid``
(or, any grid) and breaks it across the processors. It first sorts
out the layout of the processors with ``MPI_Dims_create()``. It then
evenly splits the initial grid over those processors by first
creating a new grid on each tile, linking them to the Hierarchy
linked list. It then (and here's the tricky part)
allocates each grid on the Root processor and copies data from the
Initial Grid to the new tile. Finally, it take these freshly created
root grid tiles and sends them to their new processor home.

Here's where the **bad kludge!** comes in. You'll note that in the
above description, there's an allocate on each of the newly created
tiles *on the root processor*, which will allocate more than the root
grid data. This is the problem we were trying to avoid. So
``ExternalBoundary::Prepare()`` sets ``NumberOfBaryonFields`` to zero,
so when the allocate comes around it's allocating Zero fields.

Why is it in ``ExternalBoundary::Prepare()``? A look at the lines
immediately preceding the 'kludge' help:

::

      BoundaryRank = TopGrid->GridRank;
      NumberOfBaryonFields = TopGrid->NumberOfBaryonFields;
      if (ParallelRootGridIO == TRUE)
        TopGrid->NumberOfBaryonFields = 0; /* bad kludge! */

In order to do its job properly, the ``ExternalBoundary`` objects need
to know how many ``BaryonFields`` there are in the simulation. So
``ExternalBoundary::Prepare()`` records the data, and because that's
the last place ``NumberOfBaryonFields`` is needed, sets it to zero.

When ``CommunicationPartitionGrid()`` gets to the point where it
allocates the data, ``NumberOfBaryonFields`` is now zero, so it
allocates no data. These empty root grid tiles are then distributed to
the other processors.

Finally, ``CosmologyReInitialize()`` is called, which calls
``CosmologyInitializeGrid()``. This code then resets
``NumberOfBaryonFields`` to its proper value, and since
``TotalRefinement = -1`` allocates all the data.

Then the simulation continues on, only aware of PRGIO when it comes
time to not collect the data again.



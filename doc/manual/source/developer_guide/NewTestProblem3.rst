Using Parallel Root Grid IO
===========================

First, read :doc:`../reference/HowDoesParallelRootGridIOwork`.  Come
back when you're finished.

Parallel root grid IO (PRGIO) is necessary when initializing problems
that don't fit in memory on one machine.  A PRGIO problem generator
needs to function in two passes.  First it needs to set up the basic
problem (see :ref:`UnigridInitialize`) *without* allocating any data.
This will create a temporary root grid that covers the entire domain.
Then ``CommunicationPartitionGrid`` splits this grid into several
pieces.  Usually there is one partition per MPI process unless the
parameter ``NumberOfRootGridTilesPerDimensionPerProcessor`` is greater
than 1.  The temporary root grid is then deleted, leaving only the
empty level-0 grids.  Finally each processor **re**-initializes the newly
created subgrids, this time allocating the data only when the grid
belongs to it, i.e. ``MyProcessorNumber == ProcessorNumber``.  Both
passes are done in ``InitializeNew.C``.

For an example, see either the

* ``CosmologySimulationInitialize`` and ``CosmologySimulationReInitialize``
* ``TurbulenceSimulationInitialize`` and ``TurbulenceSimulationReInitialize`` 

routines in ``InitializeNew.C``.


.. _running_large_simulations:

Running Large Simulations
=========================

Here we describe how to efficiently run a large simulation on a high
number of processors, such as particular parameters to set and
suggested number of MPI tasks for a given problem size.  For a problem
to be scalable, most of the code must be parallel to achieve high
performance numbers on large MPI process counts (see `Amdahl's
Law`__).  In general, the user wants to pick the number of processors
so that computation is still dominant over communication time.  If the
processor count is too high, communication time will become too large
and might even *slow* down the simulation!

For picking the number of processors for an Enzo run, a good starting
point is putting a 64\ :sup:`3` box on each processor for both AMR and
unigrid setups.  For example, a 256\ :sup:`3` simulation would run
well on (256/64)\ :sup:`3` = 64 processors.  For nested grid
simulations, the outer boxes usually require little computation
compared to the "zoom-in" region, so the processor count should be
based on the inner-most nested grid size.  The user can experiment
with increasing the processor count from this suggestion, but strong
scaling (i.e. linear speedup with processor count) is not to be
expected.  Little performance gains (as of v2.0) can be expected
beyond assigning a 32\ :sup:`3` cube per processor.

.. note:: 

   The level-0 grid is only partitioned during the problem
   initialization.  It will *never* be re-partitioned if the user
   restarts with a different number of processors.  However, some
   performance gains can be expected even if a processor does not
   contain a level-0 grid because of the work on finer levels.

Important Parameters
--------------------

* ``LoadBalancing``: Default is 1, which moves work from overloaded to
  underutilized processes, regardless of the grid position.  **New for
  v2.1**: In some cases but not always, speedups can be found in load
  balancing on a `space filling curve`_ (``LoadBalancing = 4``).  Here
  the grids on each processor will be continuous on the space filling
  curve.  This results in a grouped set of grids, requiring less
  communication from other processors (and even other compute nodes).

* ``SubgridSizeAutoAdjust`` and ``OptimalSubgridsPerProcessor``: **New for
  v2.1** Default is ON and 16, respectively.  The maximum subgrid size
  and edge length will be dynamically adjusted on each AMR level
  according to the number of cells on the level and number of
  processors.  The basic idea behind increasing the subgrid sizes
  (i.e. coalescing grids) is to reduce communication between grids.

* ``MinimumSubgridEdge`` and ``MaximumSubgridSize``: *Unused if
  SubgridAutoAdjust is ON*.  Increase both of these parameters to
  increase the average subgrid size, which might reduce communication
  and speedup the simulation.

* ``UnigridTranspose``: Default is 0, which is employs blocking MPI
  communication to transpose the root grid before and after the FFT.
  In level-0 grids >= 1024\ :sup:`3`, this becomes the most
  expense part of the calculation.  In these types of large runs,
  Option 2 is recommended, which uses non-blocking MPI calls; however
  it has some additional memory overhead, which is the reason it is
  not used by default.

Compile-time options
--------------------

* ``max-subgrids``: If the number of subgrids in a single AMR level
  exceeds this value, then the simulation will crash.  Increase as
  necessary.  Default: 100,000

* ``ooc-boundary-yes``: Stores the boundary conditions out of core,
  i.e. on disk.  Otherwise, each processor contains a complete copy of
  the external boundary conditions.  This becomes useful in runs with
  large level-0 grids.  For instance in a 1024\ :sup:`3` simulation
  with 16 baryon fields, each processor will contain a set of
  boundary conditions on 6 faces of 1024\ :sup:`2` with 16 baryon
  fields.  In single precision, this requires 402MB!  Default: OFF

* ``fastsib-yes``: Uses a chaining mesh to help locate sibling grids
  when constructing the boundary conditions.  Default: ON

.. |ge| unicode:: 0x2265

.. _space filling curve: http://en.wikipedia.org/wiki/Hilbert_curve

.. __: http://en.wikipedia.org/wiki/Amdahl's_law

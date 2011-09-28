Estimated Simulation Resource Requirements
==========================================

Estimating problem sizes for most Enzo calculations is at best an
inexact science, given the nature of Adaptive Mesh Refinement (AMR)
simulations. The fundamental issue with an AMR calculation in
cosmology or in many astrophysical situations where gravitational
collapse is important has to do with memory. The amount of memory
used at the beginning of the simulation (when you have a single
grid or a handful of grids) is far, far less than the memory
consumption at the end of the simulation, when there can be
hundreds of grids per processor. The amount of memory required can
easily grow by an order of magnitude over the course of a
cosmological simulation, so it is very important to make sure to
take this into account to ensure that enough memory is available in
later stages of your simulation. It is also important to realize
that in general one should try to keep the largest amount of data
per processing core that you can so that individual cores are never
data-starved. Data-starved processing units cause poor scaling, as
your CPUs will then be sitting idle while waiting for data from
other computing nodes. Computational fluid dynamics simulations are
notoriously communication-heavy, making this a challenging corner
of parameter space to operate in.

This page contains some rules of thumb that will help you along
your way, based on data collected up to the release of
Enzo v1.5 (so up to Fall 2008), when
supercomputers typically have 1GB-2GB of memory per processing unit
(a dual-processor node with two cores per processor would have 4-8
GB of memory, for example).

Cosmology or non-cosmology unigrid (non-AMR) simulations
--------------------------------------------------------

These are actually quite straightforward to predict, given that in
a unigrid simulation the grid is partitioned up in an approximately
equal fashion and then left alone. Experimentation shows that, for
machines with 1-2 GB of memory per core, one gets near-ideal
scaling with 128\ :sup:`3`\  cells per core (so a 512\ :sup:`3`\ 
cell calculations should be run on 64 processors, and a
1024\ :sup:`3`\  cell run should be done on 512 processors). This
is comfortably within memory limits for non-cosmology runs, and
there is no danger of running up against a node's memory ceiling
(which causes tremendous slowdown, if not outright program
failure). Unigrid cosmology runs have a further complication due to
the dark matter particles - these move around in space, and thus
move from processor to processor. Areas where halos and other
cosmological structures form will correspond to regions with
greater than average memory consumption. Keeping 128\ :sup:`3`\ cells
and particles per core seems to scale extremely efficiently
up to thousands of processors, though if one is using a machine
like an
`IBM Blue Gene <http://domino.research.ibm.com/comm/research_projects.nsf/pages/bluegene.index.html>`_,
which typically has far less memory per core than other computers,
one might have to go to 64\ :sup:`3`\  cells/particles per core so
that nodes corresponding to dense regions of the universe don't run
out of memory.

Cosmology adaptive mesh simulations
-----------------------------------

Scaling and problem size is much more difficult to predict for an
AMR cosmology run than for its unigrid equivalent. As discussed
above, the amount of memory consumed can grow strongly over time.
For example, a 512\ :sup:`3`\  root grid simulation with seven
levels of adaptive mesh refinement started out with 512 root grid
tiles, and ended up with over 400,000 grids! This calculation was
run on 512 processors, though memory consumption grew to the point
that it had to be run on a system where half of the cores per node
were kept this particle mass field over these processors. For each
grid, only processors with particles contribute to this sum to
reduce the amount of computation and communication. In short, this
routine performs a non-blocking ``MPI_SUM`` over a select number of
processors.

``CommunicationCollectParticles(SUBGRIDS_LOCAL)`` -- This routine
replaces ``grid::MoveSubgridParticlesFast()``. It keeps the particles on
the same processor, but this doesn't matter here because the
children grids are always created on the same processor as its
parent and then moved to another processor during load balancing.
``CommunicationCollectParticles(SIBLINGS_ONLY)`` -- After load
balancing is complete on level L\ :sub:`sub`\, we can safely move the
particles to their host processor without the worry of running out
of memory.


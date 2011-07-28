SetAccelerationBoundary (SAB)
=============================

One of the minor bugs in Enzo that was uncovered by the addition of
MHD-CT is the boundary on the gravitational acceleration field.

Enzo currently solves gravity in two phases: first by Fast Fourier
Transform on the root grid, then by multigrid relaxation on the
subgrids. Unfortunately, each subgrid is solved as an individual
problem, and is not very concious of its neighbours.

The problem with this is the ghost zones. Enzo MHD-CT is not a
divergence *free* method, but a divergence *preserving* method.
There isn't a mechanism that reduces the divergence of the magnetic
field. Unfortunately, inconsistencies in *any* fluid quantity can
lead to divergence in the magnetic field. The magnetic field is
stored on the faces of each computational zone, and are updated by
an electric field that is stored on the edges. Since this data sits
in the face of the zone, whenever two grids abut, they share a
face, so it is vital that both grids describe everything in the
stencil of the face centered fields identically, otherwise they
will get different results for the magnetic field on that face, and
divergence will be generated. It was noticed that in the case of
the ``AccelerationField`` that due to the isolated nature of the
gravity solver, the ghost zones of a subgrid didn't necessarily
equal the active zones of grids that were next to it. Thus the
Magnetic fields in the shared face would ultimately be computed
slightly differently, and divergence would show up.

The proper fix for this is replacing the gravity solver with one
that is aware of the entire subgrid hierarchy at once, but this is
quite costly in both programmer time and in compute time. Work has
begun on this project at the LCA, but has not yet been finished.

As an intermediate step, Enzo was hacked a little bit. Initially,
the main loop in ``EvolveLevel.C`` looked like this:

.. code-block:: c

     for( grid=0, grid< NumberOfGrids, grid++){
        Grid[grid]->SolvePotential
        Grid[grid]->SolveHydroEquations
     }

Among, of course, many other physics and support routines. This was
broken into two loops, and a call to ``SetBoundaryConditions()`` as
inserted between the two.

.. code-block:: c

     for( grid=0, grid< NumberOfGrids, grid++){
        Grid[grid]->SolvePotential
     }
     SetBoundaryConditions
     for( grid=0, grid< NumberOfGrids, grid++){
        Grid[grid]->SolveHydroEquations
     }

However, since ``SetBoundaryConditions()`` doesn't natively know about
the ``AccelerationField``, another kludge was done. A new set of
pointers ``ActualBaryonField`` was added to ``Grid.h``, and the true
pointers are saved here, while the ``BaryonField`` array is temporarily
pointed to ``AccelerationField``. This saved a substantial rewrite of
the boundary setting routines, at the expense of some
less-than-ideal code.

This is not a bug that makes much difference overall in cosmology
simulations, and it does not solve the problem of artificial
fragmentation that has been noticed by some groups. Cosmology tests
have been done that compare solutions both with and without this
fix, and only negligible changes appear. So for most runs, it
simply adds the expense of an extra boundary condition set.
However, with MHD-CT runs it is absolutely necessary, for explosive
divergence will show up.  Additionally, and other simulations that 
are extremely sensitive to overall conservation or consistency will require
this flag.  In any condition where the user is potentially concerned about 
we suggest running a test both with and without ``SAB``, and comparing the answers.
``SAB`` brings the compuational expense of an additional boundary condition call, and 
the memory expense of three global fields, since without it the ``AccelerationField`` exists
only on a single grid at a time, while with it all three fields must be created on the entire hierarchy
at once.  This is not a major expense on either count for most simulations.

This is controled by the preprocessor directive ``SAB``. If this is
defined, the necessary steps are taken to call the acceleration
boundary.  In the file machine make file, ``Make.mach.machine-name``, this should be
added to the variable ``MACH_DEFINES``



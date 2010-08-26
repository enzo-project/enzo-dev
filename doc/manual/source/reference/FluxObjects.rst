.. _TheFluxObject:

The Flux Object
===============

This page is intended to document the creation, use, and
destruction of the Fluxes object in Enzo. This will not be a
complete description of the Flux Correction algorithm, see the
primary references for that.

Purpose
-------

In order to keep the change in zone size across grid boundaries
consistent with the underlying conservation law, Flux Correction is
used. Basically, it makes sure that the change in Total Energy inside
a subgrid (or mass, momentum, or any other conserved quantitiy) is
equal to the flux across the boundary **as seen by both levels**. This
means that the coarse grid, which gets its solution in that space
replaced by the fine grid data, also needs to have the zones right
outside that space updated so they also see that same flux.

To facilitate this operation, the Fluxes object is used.

For each subgrid, there are two Fluxes objects, that store the flux
computed in the solver (typically ``Grid_[xyz]EulerSweep``).  One
stores the fluxes that the fine grid computes, and one stores the
fluxes that the coarse grid computes. These are stored in two objects:
a grid member fluxes ``BoundaryFluxes`` for the fine data, and fluxes
``***SubgridFluxesEstimate`` for the coarse data.

Fluxes.h
--------

The actual object can be found in
``src/enzo/Fluxes.h``.

::

    struct fluxes
    {
      long_int LeftFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
      long_int LeftFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
      long_int RightFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
      long_int RightFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
      float *LeftFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
      float *RightFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
    };

This contains two sets of arrays for the actual flux values, and 4
arrays to describe the position of the flux in the computational
domain. There is a flux on each face of the subgrid, and each flux
has a vector describing its start and end. For instance,
``LeftFluxStartGlobalIndex[0][dim]`` describes the starting index for
the X face left flux. ``LeftFluxes[densNum][0]`` describes the flux of
density across the left x face.

SubgridFluxesEstimate
---------------------

``SubgridFluxesEstimate`` is a 2 dimensional array of pointers to
Fluxes objects that a given grid patch will fill. Its indexing is like
``*SubgridFluxesEstimate[Grid][Subgrid]`` , where ``Grid`` goes over
all the grids on a level, and ``Subgrid`` goes over that grid's
subgrids PLUS ONE for the grid itself, as each grid needs to keep
track of its own boundary flux for when it communicates with the
parent. (This last element is used in conjunction with the
``BoundaryFluxes`` object, as we'll see later)

Allocation
~~~~~~~~~~

Allocation of the pointer array for the grids on this level happens
at the beginning of ``EvolveLevel``:

::

    fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];

At the beginning of the time loop, each grid has its subgrid fluxes
array allocated, and a fluxes object is allocated for each subgrid
(plus one for the grid itself)

::

     while (dtThisLevelSoFar < dtLevelAbove) {
      ... timestep computation  ...
       for (grid = 0; grid < NumberOfGrids; grid++) {
    
          // The array for the subgrids of this grid
          SubgridFluxesEstimate[grid] = new fluxes *[NumberOfSubgrids[grid]];
    
          if (MyProcessorNumber ==
              Grids[grid]->GridData->ReturnProcessorNumber()) {
    
            for( Subgrids of grid ){
              SubgridFluxesEstimate[grid][counter] = new fluxes;
              ... Setup meta data ...
            }
    
            /* and one for the grid itself */
            SubgridFluxesEstimate[grid][counter] = new fluxes;
            ... and some meta data ...
    
          }
        } // end loop over grids (create Subgrid list)

Note that in older versions of Enzo are missing the processor
check, so fluxes objects are allocated for each grid and subgrid on
each processor, causing a bit of waste. This has been fixed since Enzo
1.5.

The ``LeftFluxes`` and ``RightFluxes`` are allocated in
``Grid_SolveHydroEquations.C``

Assignment
~~~~~~~~~~

After the ``LeftFluxes`` and ``RightFluxes`` are allocated in
``Grid_SolveHydroEquations.C``, they are filled with fluxes from the
solver.  In v2.0, the C++ and FORTRAN interface with the hydrodynamics
solver was improved to avoid the previous method that juggled pointers
to a temporary array for the fluxes returned from the FORTRAN hydro
solver.  Now ``Grid_[xyz]EulerSweep.C`` allocates memory for each of
the flux variables and passes them into each of the FORTRAN hydro
routines.  This removes any size limitations that the old wrappers had
when the temporary array was too large.

.. This should be included when the EnzoMHD is implemented in the
.. public version
..
.. The MHD solver Grid_SolveMHDEquations and mhd_li.src allocates the
.. flux data in a single array for each fluid quantity, and addresses the
.. 7 dimensional (non-rectangular) array in much the same way that
.. BaryonField is accessed. This doesn't resort to any questionable
.. pointer arithmetic like the original ppm_de.src, though it hasn't been
.. as stringently vetted.  For more details, one should refer to the
.. source code.

Flux Correction
~~~~~~~~~~~~~~~

After being filled with coarse grid fluxes, ``SubgridFluxesEstimate``
is then passed into ``UpdateFromFinerGrids``, where it is used to
correct the coarse grid cells and boundary fluxes. For each
grid/subgrid, ``SubgridFluxesEstimate`` is passed into
``Grid_CorrectForRefinedFluxes`` as ``InitialFluxes``. The difference of
``InitialFluxes`` and ``RefinedFluxes`` is used to update the appropriate
zones. (Essentially, the coarse grid flux is removed from the
update of those zones ex post facto, and replaced by the average of
the (more accurate) fine grid fluxes.

See the section below for the details of ``SubgridFluxesRefined`` and
``RefinedFluxes``.

AddToBoundaryFluxes
~~~~~~~~~~~~~~~~~~~

The last thing to be done with ``SubgridFluxesEstimate`` is to update
the ``BoundaryFluxes`` object for each grid on the current level. Since
multiple fine grid timesteps are taken for each parent timestep,
the **total** flux must be stored on the grids boundary. This is
done in ``Grid_AddToBoundaryFluxes``, at the end of the EvolveLevel
timestep loop.

Deallocation
~~~~~~~~~~~~

In the same grid loop that ``BoundaryFluxes`` is updated, the
``SubgridFluxesEstimate`` object is destroyed with ``DeleteFluxes``, and
the pointers themselves are freed.

::

       for (grid = 0; grid < NumberOfGrids; grid++) {
          if (MyProcessorNumber == Grids[grid]->GridData->ReturnProcessorNumber()) {
    
           Grids[grid]->GridData->AddToBoundaryFluxes
               (SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1])
    
    
           for (subgrid = 0; subgrid < NumberOfSubgrids[grid]; subgrid++) {
    
            DeleteFluxes(SubgridFluxesEstimate[grid][subgrid]);
    
            delete SubgridFluxesEstimate[grid][subgrid];
           } 
          delete [] SubgridFluxesEstimate[grid];
         }

grid.BoundaryFluxes
-------------------

Each instance of each grid has a fluxes ``BoundaryFluxes`` object that
stores the flux across the surface of that grid. It's used to
correct it's Parent Grid.

Allocation
~~~~~~~~~~

``BoundaryFluxes`` is allocated immediately *before* the timestep loop
in ``EvolveLevel`` by the routine ``ClearBoundaryFluxes``.

Usage
~~~~~

For each grid, ``BoundaryFluxes`` is filled at the end of the
``EvolveLevel`` timestep loop by the last element of the array
``SubgridFluxesEstimate[grid]`` for that grid. This is additive,
since each grid will have multiple timesteps that it must correct
its parent for. This is done by ``AddToBoundaryFluxes``, as
described above.

``BoundaryFluxes`` is used in ``UpdateFromFinerGrids`` to populate another
fluxes object, ``SubgridFluxesRefined``. This is done in
``GetProjectedBoundaryFluxes``. The values in ``SubgridFluxesRefined`` are
area weighted averages of the values in ``BoundaryFluxes``, coarsened
by the refinement factor of the simulation. (So for factor of 2
refinement, ``SubgridFluxesRefined`` has half the number of zones in
each direction than ``BoundaryFluxes``, and matches the cell width of
the parent grid.)

``BoundaryFluxes`` is also updated from subgrids in
``CorrectForRefinedFluxes``. This happens when a subgrid boundary lines
up exactly with a parent grid boundary. However, in many versions
of Enzo, this is deactivated by the following code:

::

            CorrectLeftBoundaryFlux = FALSE;
            CorrectRightBoundaryFlux = FALSE;
    #ifdef UNUSED
            if (Start[dim] == GridStartIndex[dim]-1)
              CorrectLeftBoundaryFlux = TRUE;
            if (Start[dim] + Offset == GridEndIndex[dim]+1)
              CorrectRightBoundaryFlux = TRUE;
    #endif /* UNUSED */

It is unclear why this is, but removal of the UNUSED lines restores
conservation in the code, and is essential for proper functioning
of the MHD version of the code (which will be released in the
future.) I have seen no problems from removing this code.

Many implementations of block structured AMR require a layer of
zones between parent and subgrid boundaries. Enzo is not one of
these codes.

Deallocation
~~~~~~~~~~~~

``BoundaryFluxes`` is only deleted once the grid itself is deleted.
This happens mostly in ``RebuildHierarchy``.

SubgridFluxesRefined
--------------------

The final instance of a fluxes object is fluxes
``SubgridFluxesRefined``. This object takes the fine grid fluxes,
resampled to the coarse grid resolution, and is used to perform the
flux correction itself. This section is short, as its existance has
been largely documented in the previous sections.

Allocation
~~~~~~~~~~

``SubgridFluxesRefined`` is declared in ``UpdateFromFinerGrids``. The
actual allocation occurs in ``Grid_GetProjectedBoundaryFluxes``, where
it's passed in as ``ProjectedFluxes``.

Usage
~~~~~

``SubgridFluxesRefined`` is also filled in
``Grid_GetProjectedBoundaryFluxes``, as the area weighted average of
the subgrid boundary flux.

It is then passed into ``Grid_CorrectForRefinedFluxes``, Here, it is
used to update the coarse grid zones that need updating.

Deallocation
~~~~~~~~~~~~

``SubgridFluxesRefined`` is deleted after it is used in
``Grid_CorrectForRefinedFluxes``.


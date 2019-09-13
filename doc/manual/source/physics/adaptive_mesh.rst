.. _adaptive_mesh:

Adaptive Mesh Refinement
========================

Enzo's adaptive mesh implements a version of the `Berger and Colella (1989)
<https://ui.adsabs.harvard.edu/abs/1989JCoPh..82...64B/abstract>`_
block-structured AMR algorithm on a Cartesian mesh.  In this
algorithm, cells that are flagged for refinement are combined into
rectangular solid "child grid" patches with cells whose spatial resolution are
an integer multiple more finely resolved than their coarser "parent
grids" (with the ratio of resolutions typically, but not necessarily,
being 2).  This method is distinct from cell-based AMR codes in that
cells are aggregated into grids, and distinct from oct-tree block
structured AMR in that the grids that are created can be of arbitrary
size and aspect ratio (i.e., each grid dimension can differ, rather
than being a cube) and can be located at arbitrary locations within
the parent grid rather than in octants of the parent grid.  Enzo does
not implement the full Berger and Colella method - for the sake of
efficiency, it is restricted in the following ways:

* Higher-resolution child grids must be contained completely within
  their parent grids, rather than spanning multiple parent grids.
  (Parent grids of level L are, however, allowed to have
  multiple child grids of level L+1.)
* The edges of child grids must align with the cell edges of parent
  grids (which also implicitly requires that child grids align
  with the edges of their parent grids).
* All grids are aligned with the principle axes (x,y,z), and may not
  be arbitrarily rotated with respect to those axes.

Enzo implements time subcycling within its AMR, with every grid level L
determining its own timestep and all grids at that level taking the
same timestep.  This timestep is restricted so that the timestep at
level L may not exceed the timestep at level L-1.
The `Enzo method paper <https://ui.adsabs.harvard.edu/abs/2014ApJS..211...19B/abstract>`_ 
describes the algorithm in more detail.

The parameters controlling Enzo's grid hierarchy can be found in :ref:`hierarchy_control_parameters`.


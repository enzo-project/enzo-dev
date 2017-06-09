The Enzo Hierarchy File - Explanation and Usage
===============================================

The Enzo Hierarchy file is a representation of the internal memory
state of the entire hierarchy of grids. As such, its format --
while somewhat obtuse at first -- reflects that context. Each grid
entry has a set number of fields that describe its position in
space, as well as the fields that are affiliated with that grid:

Note: We are in the process of transitioning to an `HDF5-formatted
Hierarchy File`_.

.. highlight:: none

::

    Grid = 1
    Task              = 4
    GridRank          = 3
    GridDimension     = 38 22 22 
    GridStartIndex    = 3 3 3 
    GridEndIndex      = 34 18 18 
    GridLeftEdge      = 0 0 0 
    GridRightEdge     = 1 0.5 0.5 
    Time              = 646.75066015177
    SubgridsAreStatic = 0
    NumberOfBaryonFields = 8
    FieldType = 0 1 4 5 6 19 20 21 
    BaryonFileName = ./RD0005/RedshiftOutput0005.cpu0000
    CourantSafetyNumber    = 0.300000
    PPMFlatteningParameter = 0
    PPMDiffusionParameter  = 0
    PPMSteepeningParameter = 0
    NumberOfParticles   = 20
    ParticleFileName = ./RD0005/RedshiftOutput0005.cpu0000
    GravityBoundaryType = 0
    Pointer: Grid[1]->NextGridThisLevel = 2

The final field, starting with "Pointer", is slightly more
complicated and will be discussed below.

``Grid = 1``

  This is the ID of the grid. Enzo grids are indexed internally
  starting at 1.

``Task = 3``

  This grid was written by processor 3 and will be read in by it if
  restarting more than 4 processors.

``GridRank = 3``

  This is the dimensionality of the grid.

``GridDimension = 38 22 22``

  Dimensions, *including* ghost zones.

``GridStartIndex = 3 3 3``

  The first index of data values *owned* by this grid.

``GridEndIndex = 34 18 18``

  The final index *owned* by this grid. The active zones have
  dimensionality of GridEndIndex - GridStartIndex + 1.

``GridLeftEdge = 0 0 0``

  In code units, between ``DomainLeftEdge`` and ``DomainRightEdge``,
  the origin of this grid.

``GridRightEdge = 1 0.5 0.5``

  In code units, between ``DomainLeftEdge`` and ``DomainRightEdge``,
  the right-edge of this grid. ``dx = (GridRightEdge -
  GridLeftEdge)/(GridEndIndex - GridStartIndex + 1)``.


``Time = 646.75066015177``

  The current time (in code units) to which the baryon values in this
  grid have been evolved.


``SubgridsAreStatic = 0``

  Whether refinement can occur in the subgrids.

``NumberOfBaryonFields = 8``

  The number of data fields associated with this grid.

``FieldType = 0 1 4 5 6 19 20 21``

  The integer identifiers of each field, in order, inside this grid.

``BaryonFileName = ./RD0005/RedshiftOutput0005.cpu0000``

  The HDF5 file in which the baryons fields are stored.

``CourantSafetyNumber = 0.300000``

  Courant safety number for this grid (governs timestepping.)

``PPMFlatteningParameter = 0``

  Flattening parameter for this grid (governs PPM hydro.)

``PPMDiffusionParameter = 0``

  Diffusion parameter for this grid (governs PPM hydro.)

``PPMSteepeningParameter = 0``

  Steepening parameter for this grid (governs PPM hydro.)

``NumberOfParticles = 20``

  How many particles are located in this grid at this timestep.

``ParticleFileName = ./RD0005/RedshiftOutput0005.cpu0000``

  The HDF5 file in which the baryon fields and particle data are
  stored.  This field will not exist if there aren't any particles in
  the grid.

``GravityBoundaryType = 0``

   Boundary type inside gravity solver.




HDF5-formatted Hierarchy File
-----------------------------

We are transitioning to an HDF5-formatted hierarchy file. This is an
improvement because reading a large (many thousand grid) ASCII
hierarchy file take a long time, and be a possible cause of precision
errors in deep hierarchies.

The structure of the file:

Although HDF5 tools like 'h5ls' and 'h5dump' can be used to explore
the structure of the file, it's probably easiest to use python and
h5py. This is how to open an example hierarchy file (from
run/Cosmology/Hydro/AMRCosmologySimulation) in python.

::

     >>> import h5py
     >>> f = h5py.File('RD0007/RedshiftOutput0007.hierarchy.hdf5','r')

The root group ('/') contains a number of attributes.

::

     >>> f.attrs.keys()
     ['Redshift', 'NumberOfProcessors', 'TotalNumberOfGrids']
     >>> f.attrs['Redshift']
     0.0
     >>> f.attrs['NumberOfProcessors']
     1
     >>> f.attrs['TotalNumberOfGrids']
     44

So we see that this is a z=0 output from a simulation run on a single
core and it contains a total of 44 grids.

Now let's look at the groups contained in this file.

::

     >>> f.keys()
     ['Level0', 'Level1', 'Level2', 'LevelLookupTable']

The simulation has two levels of refinement, so there are a total of
three HDF5 groups that contain information about the grids at each
level. Additionally, there is one more dataset ('LevelLookupTable')
that is useful for finding which level a given grid belongs to. Let's
have a closer look.

::

     >>> level_lookup = f['LevelLookupTable']
     >>> level_lookup.shape
     (44,)
     >>> level_lookup[:]
     array([0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

This shows you that the first grid is on level 0, the second on level
1, and all the remaining grids on level 2. Let's have a look at the
'Level2' group.

::

     >>> g = f['Level2']
     >>> g.keys()
    ['Grid00000003', 'Grid00000004', 'Grid00000005', ..., 'Grid00000043', 'Grid00000044']

Each level group also has one attribute, 'NumberOfGrids'.

::

     >>> g.attrs['NumberOfGrids']
     42

The hierarchy information about each of the grids is stored as both
attributes and datasets.

::

     >>> grid = g['Grid00000003']
     >>> grid.attrs.keys()
     ['Task', 'GridRank', 'Time', 'OldTime', 'SubgridsAreStatic', 'NumberOfBaryonFields', 'FieldType',
      'BaryonFileName', 'CourantSafetyNumber', 'PPMFlatteningParameter', 'PPMDiffusionParameter',
      'PPMSteepeningParameter', 'ParticleFileName', 'GravityBoundaryType', 'NumberOfDaughterGrids',
      'NextGridThisLevelID', 'NextGridNextLevelID']
     >>> grid.keys()
     ['GridDimension', 'GridEndIndex', 'GridGlobalPosition',
      'GridLeftEdge', 'GridRightEdge', 'GridStartIndex', 'NumberOfParticles']

Besides the parameters that have been described above, there are few
new elements:

``GridGlobalPosition`` is LeftGridEdge[] expressed in integer indices
of this level, i.e. running from 0 to RootGridDimension[] *
RefinementFactors[]**level - 1. This may be useful for re-calculating
positions in long double precision (which is not universally supported
by HDF5) at runtime.


``NumberOfDaughterGrids`` gives you the number of daughter grids.


``DaughterGrids`` is a group that contains HDF5-internal soft links to
the daugher datasets. Example:

::

     >>> daughters = grid['DaughterGrids']
     >>> daughters.keys()
     ['DaughterGrid0000', 'DaughterGrid0001', 'DaughterGrid0002', ..., 'DaughterGrid0041']
     >>> daughters.get('DaughterGrid0000', getlink=True)
     <SoftLink to "/Level2/Grid00000003">

In this case there are 42 daughter grids.


``ParentGrids`` is a group that contains HDF5-internal soft links to
parent grids on all levels above the present grid's level. Example for
a level 2 grid:

::

     >>> grid = f['Level2']['Grid00000044']
     >>> parents = grid['ParentGrids']
     >>> parents.keys()
     ['ParentGrid_Level0', 'ParentGrid_Level1']
     >>> parents.get('ParentGrid_Level0', getlink=True)
     <SoftLink to "/Level0/Grid00000001">

Lastly, there's one additional (experimental) feature that is
available only if you've compiled with verson 1.8+ of HDF5. In that
case you can set '#define HAVE_HDF5_18' in
Grid_WriteHierarchyInformationHDF5.C [perhaps this should become a
Makefile configuration option?], and then there will be an external
HDF5 link to the HDF5 file containing the actual data for that grid. Example:

::

     >>> grid.get('GridData', getlink=True)
     >>> <ExternalLink to "Grid00000002" in file "./RD0007/RedshiftOutput0007.cpu0000"


.. _controlling_the_hierarhcy_file_output: 

Controlling the Hierarchy File Output Format
--------------------------------------------

There are two new parameters governing the format of the hierarchy
format:

``[OutputControl.]HierarchyFileInputFormat = 0, 1``

  This specifies the format of the hierarchy file to be read in: 0 =
  ASCII, 1 = HDF5. Default set to 0 for now, but will change to 1 in the
  future.

``[OutputControl.]HierarchyFileOutputFormat = 0, 1, 2``  [OutputControl.HierarchyFileOutputFormat in new-config]

  This specifies the format of the hierarchy file to be written out: 0
  = ASCII, 1 = HDF5, 2 = both. Default set to 2 for now, but will change
  to 1 in the future.

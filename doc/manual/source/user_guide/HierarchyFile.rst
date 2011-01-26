The Enzo Hierarchy File - Explanation and Usage
===============================================

The Enzo Hierarchy file is a representation of the internal memory
state of the entire hierarchy of grids. As such, its format --
while somewhat obtuse at first -- reflects that context. Each grid
entry has a set number of fields that describe its position in
space, as well as the fields that are affiliated with that grid:

::

    Grid = 1
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
    BaryonFileName = /scratch/username/RD0005/RedshiftOutput0005.cpu0000
    CourantSafetyNumber    = 0.300000
    PPMFlatteningParameter = 0
    PPMDiffusionParameter  = 0
    PPMSteepeningParameter = 0
    NumberOfParticles   = 20
    ParticleFileName = /scratch/username/RD0005/RedshiftOutput0005.cpu0000
    GravityBoundaryType = 0
    Pointer: Grid[1]->NextGridThisLevel = 2

The final field, starting with "Pointer", is slightly more
complicated and will be discussed below.

    Grid = 1

        This is the ID of the grid. Enzo grids are indexed internally
        starting at 1.


    GridRank = 3

        This is the dimensionality of the grid.


    GridDimension = 38 22 22

        Dimensions, *including* ghost zones.


    GridStartIndex = 3 3 3

        The first index of data values *owned* by this grid.


    GridEndIndex = 34 18 18

        The final index *owned* by this grid. The active zones have
        dimensionality of GridEndIndex - GridStartIndex + 1.


    GridLeftEdge = 0 0 0

        In code units, between `DomainLeftEdge? </wiki/DomainLeftEdge>`_
        and `DomainRightEdge? </wiki/DomainRightEdge>`_, the origin of this
        grid.


    GridRightEdge = 1 0.5 0.5

        In code units, between `DomainLeftEdge? </wiki/DomainLeftEdge>`_
        and `DomainRightEdge? </wiki/DomainRightEdge>`_, the right-edge of
        this grid. dx = (GridRightEdge - GridLeftEdge)/(GridEndIndex -
        GridStartIndex + 1).


    Time = 646.75066015177

        The current time to which the baryon values in this grid have been
        evolved.


    SubgridsAreStatic = 0

        Whether refinement can occur in the subgrids.


    NumberOfBaryonFields = 8

        The number of data fields associated with this grid.


    FieldType = 0 1 4 5 6 19 20 21

        The integer identifiers of each field, in order, inside this
        grid.


    BaryonFileName =
    /scratch/username/RD0005/RedshiftOutput0005.cpu0000

        The HDF5 file in which the baryons fields are stored.


    CourantSafetyNumber = 0.300000

        Courant safety number for this grid (governs timestepping.)


    PPMFlatteningParameter = 0

        Flattening parameter for this grid (governs PPM hydro.)


    PPMDiffusionParameter = 0

        Diffusion parameter for this grid (governs PPM hydro.)


    PPMSteepeningParameter = 0

        Steepening parameter for this grid (governs PPM hydro.)


    NumberOfParticles = 20

        How many particles are located in this grid at this timestep.


    ParticleFileName =
    /scratch/username/RD0005/RedshiftOutput0005.cpu0000

        The HDF5 file in which the baryon fields and particle data are stored.


    GravityBoundaryType = 0

        Boundary type inside gravity solver.





.. _analyzing_with_yt:

Analyzing Enzo Data
=================

Using YT
--------

The Enzo community predominantly makes use of the Python package
YT to inspect and analyze data.  YT
is completely free and open source, with an active and expanding
development community, and it presents to the user both high-level
and low-level APIs. The
`documentation <http://yt-project.org/doc/>`_ contains a
tutorial as well as an API reference, but here we will step through
some simple steps toward creating script to make simple plots of a
cosmological simulation.

For more information on yt, see the `yt website <http://yt-project.org>`_,
where you will find mailing lists, documentation, API documentation, a cookbook
and even a gallery of images.

HDF5 Tools
----------

Enzo reads in initial conditions files and outputs simulation data using the
`HDF5 <http://www.hdfgroup.org/>`_ structured data format (created and
maintained by the NCSA HDF group). Though this format takes a bit more effort
to code than pure C/C++ binary output, we find that the advantages are worth
it. Unlike raw binary, HDF5 is completely machine-portable and the HDF5
library takes care of error checking. There are many useful standalone
utilities included in the HDF5 package that allow a user to examine the
contents and structure of a dataset. In addition, there are several
visualization and data analysis packages that are HDF5-compatible. See the
page on Data Vizualization for more information about this. The NCSA HDF group
has an excellent tutorial on working with HDF5.

Note that as of the Enzo 2.0 code release, Enzo still supports reading the HDF4
data format, but not writing to it. We strongly suggest that new users
completely avoid this and use the HDF5 version instead. Enzo's parallel IO
only works with HDF5, and we are encouraging users migrate as soon as is
feasible.

Writing your own tools, I - the Enzo Grid Hierarchy
---------------------------------------------------

Enzo outputs each individual adaptive mesh block as its own grid
file. Each of these files is completely self-contained, and has
information about all of the grid cells that are within that volume
of space. Information on the size and spatial location of a given
grid file can be obtained from the hierarchy file, which has the
file extension ".hierarchy". This ascii file has a listing for each
grid that looks something like this:

.. highlight:: none

::

    Grid = 26
    GridRank          = 3
    GridDimension     = 34 22 28 
    GridStartIndex    = 3 3 3 
    GridEndIndex      = 30 18 24 
    GridLeftEdge      = 0.5 0.28125 0.078125 
    GridRightEdge     = 0.71875 0.40625 0.25 
    Time              = 101.45392321467
    SubgridsAreStatic = 0
    NumberOfBaryonFields = 5
    FieldType = 0 1 4 5 6 
    BaryonFileName = RedshiftOutput0011.grid0026
    CourantSafetyNumber    = 0.600000
    PPMFlatteningParameter = 0
    PPMDiffusionParameter  = 0
    PPMSteepeningParameter = 0
    NumberOfParticles   = 804
    ParticleFileName = RedshiftOutput0011.grid0026
    GravityBoundaryType = 2
    Pointer: Grid[26]->NextGridThisLevel = 27

``GridRank`` gives the dimensionality of the grid (this one is 3D),
``GridDimension`` gives the grid size in grid cells, including ghost
zones. ``GridStartIndex`` and ``GridEndIndex`` give the starting and ending
indices of the non-ghost zone cells, respectively. The total size
of the baryon datasets in each grid along dimension i is (1+
``GridEndIndex[i]`` - ``GridStartIndex[i]``). ``GridLeftEdge`` and
``GridRightEdge`` give the physical edges of the grids (without ghost
zones) in each dimension. ``NumberOfParticles`` gives the number of
dark matter particles (and/or star particles, for simulations
containing star particles) in a given grid. Note that when there
are multiple grids covering a given region of space at various
levels of resolution, particles are stored in the most highly
refined grid. ``BaryonFileName`` is the name of the actual grid file,
and should be the same as ``ParticleFileName``. ``Time`` is the simulation
time, and should be the same as ``InitialTime`` in the parameter file
for the same data dump. The other parameters for each entry are
more advanced and probably not relevant for simple data analysis.

Possibly the greatest source of potential confusion in Enzo's
datasets is the overlap of grid cells. In a simulation, when a
given grid is further refined, the coarse cells which have not been
refined are still kept. The solution to the hydro and gravity
equations are still calculated on that level, but are updated with
information from more highly refined levels. What this is means is
that a volume of space which has been refined beyond the root grid
is covered by multiple grid patches at different levels of
resolution. Typically, when doing analysis you only want the most
highly refined information for a given region of space (or the most
highly refined up to a certain level) so that you don't
double-count (or worse) the gas in a given cell. Look at this
example analysis code.

.. _EnzoPhysicalUnits:

Writing your own tools, II - Enzo Physical Units
------------------------------------------------

Yet another significant source of confusion is the units that Enzo
uses. When doing a cosmology simulation, the code uses a set of
units that make most quantities on the order of unity (in
principle). The Enzo manual section on
the code output format :ref:`EnzoOutputFormats`
explains how to convert code units to cgs units. However, there are
some subtleties:

**Density fields**
    All density fields are in the units described in the AMR guide
    **except** electron density. Electron density is only output when
    ``MultiSpecies`` is turned on, and in order to convert the electron
    density to cgs it must be multiplied by the code density conversion
    factor and then (m\:sub:`e`\/m\:sub:`p`\), where
    m\:sub:`e`\ and m\:sub:`p`\ are the electron
    and proton rest masses (making electron density units different
    from the other fields by a factor of m\:sub:`e`\/m\:sub:`p`\).
    The reason this is
    done is so that in the code the electron density can be computed
    directly from the abundances of the ionized species.
**Energy fields**
    There are two possible energy fields that appear in the code - Gas
    energy and total energy. Both are in units of **specific energy**,
    ie, energy per unit mass. When Zeus hydro is being used
    (``HydroMethod`` = 2, there should be only one energy field - "total
    energy". This is a misnomer - the Zeus hydro method only follows
    the specific internal (ie, thermal) energy of the gas explicitly.
    When the total energy is needed, it is calculated from the
    velocities. When PPM is used (``HydroMethod`` = 0) the number of energy
    fields depends on whether or not ``DualEnergyFormalism`` is turned on
    or off. If it is ON (1), there is a "gas energy" field and a "total
    energy" field, where "gas energy" is the specific internal energy
    and "total energy" is "gas energy" plus the specific kinetic energy
    of the gas in that cell. If ``DualEnergyFormalism`` is OFF (0), there
    should only be "total energy", which is kinetic+internal specific
    energies. Confused yet?
**Particle mass field**
    Particle "masses" are actually stored as densities. This is to
    facilitate calculation of the gravitational potential. The net
    result of this is that, in order to calculate the stored particle
    "mass" to a physical mass, you must first multiply this field by the volume of
    a cell in which the particle resides.
    Remember that particle data is only stored in the most refined grid that
    covers that portion of the simulational volume.
    
    
When the simulation is done, Enzo will display the message
"Successful run, exiting."
Enzo is a complicated code, with a similarly complicated output
format. See the Enzo User Guide page on
the Enzo output format :ref:`EnzoOutputFormats` for
more information on the data outputs.

Congratulations! If you've made it this far, you have now
successfully run a simulation using Enzo!



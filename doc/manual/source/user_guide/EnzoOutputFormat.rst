.. _EnzoOutputFormats:

Enzo Output Formats
===================

Although there are a number of ways of specifying when (and how
often) Enzo outputs information, there is only one type of output
'dump' (well, not quite -- there are now movie dumps, see below),
which can also be used to restart the simulation. The output format
uses the following files, each of which begins with the output
name, here we use the example base_name, and are then followed by
the output number, ranging from 0000 to 9999 (if more than 10000
grids are generated then the number goes to 10000, etc.). When
restarting, or other times when an output filename needs to be
specified, use the name without any extension (e.g. ``enzo -r
base_name0000``).

Summary of Files
----------------

**base_name0000**
    This ascii file contains a complete listing of all the parameter
    settings, both those specified in the initial parameter file, as
    well as all those for which default values were assumed. The
    parameters (see :ref:`parameters`) are in the same format
    as that used in the input file: parameter_name = value. This file
    is modifiable if you would like to restart from a certain point
    with different parameter values.
**base_name0000.hierarchy**
    This ascii file specifies the hierarchy structure as well as the
    names of the grid files, their sizes, and what they contain. It
    should not be modified.
**base_name0000.cpu00001**
    The field information for each cpu (padded with zeros) is contained
    in separate files with a root 'Node' for each grid, padded with
    zeros to be eight digits. The format is the Hierarchy Data Format
    (HDF) version 5, a self-describing machine-independent data format
    developed and supported by the National Center for Supercomputing
    Applications (NCSA). More information can be found on their
    `home page <http://www.hdfgroup.org/>`_. Most scientific
    visualization packages support this format. Each field is stored as
    it's own one-, two- or three-dimensional Scientific Data Set (SDS),
    and is named for identification. Particles (if any) are included
    with a set of one-dimensional datasets under the top 'grid' node.
**base_name0000.boundary**
    An ascii file which specifies boundary information. It is not
    generally useful to modify.
**base_name0000.boundary.hdf**
    Contains field-specific boundary information, in HDF format.
**base_name0000.radiation**
    This ascii file is only generated if using the self-consistent
    radiation field.

Output Units
------------

The units of the physical quantities in the grid SDS's are depend
on the problem being run. For most test problems there is no
physical length or time specified, so they can be be simply scaled.
For cosmology there are a set of units designed to make most
quantities of order unity (so single precision variables can be
used). These units are defined below (rho0 =
3\*OmegaMatterNow\*(100\*HubbleConstantNow
km/s/Mpc)\ :sup:`2`\ /(8\*Pi\*G)).


-  length: ComovingBoxSize/HubbleConstantNow \* Mpc / (1+z)
-  density: rho0 \* (1+z)\ :sup:`3`\ 
-  time: 1/sqrt(4\*Pi\*G\*rho0\*(1+InitialRedshift)\ :sup:`3`\ )
-  temperature: K
-  velocity: (length/time)\*(1+z)/(1+InitialRedshift) (this is z
   independent)

The conversion factor is also given in the ascii output file
(base_name0000): search for ``DataCGSConversionFactor``. Each field
has its own conversation factor, which converts that field to cgs
units. Users can also set completely arbitrary internal units, as
long as they are self-consistent: to see how to do this, go to
:ref:`EnzoInternalUnits`.

.. _StreamingDataFormat:

Streaming Data Format
---------------------

**Purpose:** To provide data on every N-th timestep of each AMR
level.

Method
~~~~~~

We keep track of the elapsed timesteps on every AMR level.  Every N-th
timestep on a particular level L, all grids on levels >= L are written
for the baryon fields (specified by the user in ``MovieDataField``)
and particles. The integers in ``MovieDataField`` correspond to the
field element in ``BaryonField``, i.e. 0 = Density, 7 = HII
density. Temperature has a special value of 1000.

See :ref:`streaming_param` for a full description of the streaming
data format parameters.

File format
~~~~~~~~~~~

All files are written in HDF5 with one file per processor per
top-level timestep. The filename is named AmiraDataXXXX_PYYY.hdf5
where XXXX is the file counter, which should equal the cycle
number, and YYY is the processor number. Each file has a header
indicating


-  whether the data are cell-centered (1) or vertex-centered (0)
   [int]
-  number of baryon fields written [int]
-  number of particle fields written [int]
-  field names with the baryon fields first, followed by the
   particle fields [array of variable-length strings]

The group names (grid-%d) are unique only in the file. Unique grids
are identified by their timestep number attribute and position.
Each
grid has the following attributes:


-  AMR level [int]
-  Timestep [int]
-  Code time [double]
-  Redshift [double]
-  Ghost zones flag for each grid face [6 x int]
-  Number of ghost zones in each dimension [3 x int]
-  Cell width [3 x double]
-  Grid origin in code units [3 x double]
-  Grid origin in units of cell widths [3 x long long]

In addition to the HDF5 files, a binary index file is created for
fast I/O in post-processing. The filenames of the these files are the
same as the main data files but with the extension .idx. The header
consists of


-  pi (to indicate endianness) [float]
-  cell width on the top level [float]
-  number of fields [char]
-  cell-centered (1) or vertex-centered (0) [char]
-  field names [number of fields x (64 char)]

For every grid written, an index entry is created with


-  grid ID [int]
-  code time [double]
-  timestep [int]
-  redshift [double]
-  level [char]
-  grid origin in units of cell widths [long long]
-  grid dimensions [short]
-  number of particles [int]

Lastly, we output an ASCII file with the code times and redshifts of every top
level timestep for convenience when choosing files to read afterwards.


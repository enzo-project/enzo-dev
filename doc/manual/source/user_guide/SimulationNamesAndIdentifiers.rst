Simulation Names and Identifiers
================================

To help track and identify simulations and datasets, a few new
lines have been added to the parameter file:

::

    MetaDataIdentifier                = <short string persisted across datasets>
    MetaDataSimulationUUID            = <uuid persisted across datasets>
    MetaDataDatasetUUID               = <unique dataset uuid>
    MetaDataRestartDatasetUUID        = <input dataset uuid>
    MetaDataInitialConditionsUUID     = <initial conditions uuid>

The parameters stored during a run are members of the TopGridData
struct.

::

      char *MetaDataIdentifier;     // A name (string) that will be persisted between datasets
      char *SimulationUUID;         // Unique identifier for the simulation
      char *RestartDatasetUUID;     // Identifier of the dataset restarting from
      char *InitialConditionsUUID;  // Identifier of the initial conditions used

MetaDataIdentifier
------------------

This is a character string without spaces (specifically, something
that can be picked by "%s"), that can be defined in a parameter
file, and will be written out in every following output. It's
intended to be a human-friendly way of tracking datasets. For
example

Example:

::

    MetaDataIdentifier   = Cosmology512_Mpc_run4

MetaDataSimulationUUID
----------------------

The ``MetaDataSimulationUUID`` is a globally unique identifier for a
collection of datasets.
`Universally Unique Identifiers <http://en.wikipedia.org/wiki/Universally_Unique_Identifier>`_
(UUIDs) are opaque identifiers using random 128-bit numbers, with
an extremely low chance of collision. Therefore, they are very
useful when trying to label data coming from multiple remote
resources (say, computers distributed around the world).

Example:

::

    MetaDataSimulationUUID             = e5f72b77-5258-45ba-a376-ffe11907fae1

Like the ``MetaDataIdentifier``, the ``MetaDataSimulationUUID`` is read in
at the beginning of a run, and then re-written with each output.
However, if one is not found initially, a new one will be
generated, using code from the
`ooid library <http://sourceforge.net/projects/ooid/>`_ included
in Enzo.

If you want to define one up front, there are several ways. The
tool uuidgen is on most POSIX systems, and UUID libraries exist for
the languages we commonly use.

Some examples:

::

    namor:~ rpwagner$ uuidgen 
    2521B497-CBE7-4849-8642-2766793241EF
    namor:~ rpwagner$ python
    Python 2.6.4 (r264:75706, Jan 26 2010, 00:29:28) 
    [GCC 4.2.1 (Apple Inc. build 5646)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import uuid
    >>> str(uuid.uuid4()).lower()
    'd731bdce-f993-4a91-aa2a-a1b1af2d85e0'
    >>> 

MetaDataDatasetUUID
-------------------

A ``MetaDataDatasetUUID`` is created at each output.

Example:

::

    MetaDataDatasetUUID                = b9d78cc7-2ecf-4d66-a23c-a1dcd40e7955

MetaDataRestartDatasetUUID
--------------------------

While reading the parameter file, if a ``MetaDataDatasetUUID`` line is
found, it is stored, and re-written as ``MetaDataRestartDatasetUUID``.
The intention of this is help track datasets across restarts and
parameter tweaks.

Example:

::

    MetaDataRestartDatasetUUID   = b9d78cc7-2ecf-4d66-a23c-a1dcd40e7955

MetaDataInitialConditionsUUID
-----------------------------

This is simih3>

nprocs
    number of processors
boxsize
    box size in comoving Mpc (not Mpc/h)
resolution
    top-level grid resolution
n\_levels
    level of the finest nested grid
inner\_width
    width of the finest nested grid
buffer\_cells
    number of cells separating nested grids
seed
    random seed (must be 9 digits)
name
    name of the data directory (saved in mpgrafic/data/name/)
center
    how much to shift the data in order to center on a particular
    region.
LargeScaleCorrection
    whether to use a noise file from a lower-resolution run
LargeScaleFile
    noise file from that lower-resolution run
OneDimPerFile
    whether we're using one file per velocity component
omega\_m
    Omega matter
omega\_v
    Omega lambda
omega\_b
    Omega baryon
h0
    Hubble constant in units of [km/s/Mpc]
sigma8
    sigma\_8
n\_plawslope
    slope of power spectrum

After you set your parameters, run this script with

::

    python make_ic.py 

and it will re-compile mpgrafic and (for nested grids)
degraf. Then it will run mpgrafic for the full-resolution box. If
the user wants nested grids, it will copy the data files to
mpgrafic/degraf and create the set of nested grid files.

The user cannot specify the initial redshift because mpgrafic
determines it from the parameter sigstart that is the maximum
initial
density fluctuation. From this, mpgrafic calculates the initial
redshift. This file is overwritten by the python script, so if you
want to change this parameter, change it in the python script
(routine
write\_grafic1inc).

The noise file is always kept in mpgrafic/mpgrafic-0.2/src and is
named $seed\_$resolution.dat, where $resolution is the top-level
grid
resolution. It can be re-used with
`LargeScaleFile? </wiki/LargeScaleFile>`_ if the user wants
to re-simulate the volume at a higher resolution.

The data files are moved to mpgrafic/data/$name. If nested grids
were
created, degraf writes a set of parameters in enzo.params for
copy-pasting into an Enzo parameter file. Now you can move the
files
to the simulation directory and start your Enzo cosmology
simulation!



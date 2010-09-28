.. _SimulationNamesAndIdentifiers:

Simulation Names and Identifiers
================================

To help track and identify simulations and datasets, a few new
lines have been added to the parameter file:

``MetaDataIdentifier``
   short string persisted across datasets
``MetaDataSimulationUUID``
   uuid persisted across datasets
``MetaDataDatasetUUID``
   unique dataset uuid
``MetaDataRestartDatasetUUID``
   input dataset uuid
``MetaDataInitialConditionsUUID``
   initial conditions uuid


The parameters stored during a run are members of the
TopGridData struct.

MetaDataIdentifier
------------------

This is a character string without spaces (specifically, something
that can be picked by "%s"), that can be defined in a parameter
file, and will be written out in every following output. It's
intended to be a human-friendly way of tracking datasets. For
example

Example:

::

    MetaDataIdentifier = Cosmology512_Mpc_run4


MetaDataSimulationUUID
----------------------

The MetaDataSimulationUUID is a globally unique identifier for a collection of
datasets.   `Universally Unique Identifiers
<http://en.wikipedia.org/wiki/Universally_Unique_Identifier>`_ (UUIDs) are
opaque identifiers using random 128-bit numbers, with an extremely low chance
of collision. Therefore, they are very useful when trying to label data coming
from multiple remote resources (say, computers distributed around the world).

Example:

::

    MetaDataSimulationUUID = e5f72b77-5258-45ba-a376-ffe11907fae1


Like the ``MetaDataIdentifier``, the ``MetaDataSimulationUUID`` is read in at
the beginning of a run, and then re-written with each output.  However, if one
is not found initially, a new one will be generated, using code from the `ooid
library <http://sourceforge.net/projects/ooid/>`_ included in Enzo.

UUIDs can be generated with a variety of tools, including the python standard
library.

MetaDataDatasetUUID
-------------------

A MetaDataDatasetUUID is created at each output.

Example:

::

    MetaDataDatasetUUID = b9d78cc7-2ecf-4d66-a23c-a1dcd40e7955


MetaDataRestartDatasetUUID

--------------

While reading the parameter file, if a MetaDataDatasetUUID line is
found, it is stored, and re-written as MetaDataRestartDatasetUUID.
The intention of this is help track datasets across restarts and
parameter tweaks.

Example:

::

    MetaDataRestartDatasetUUID = b9d78cc7-2ecf-4d66-a23c-a1dcd40e7955

MetaDataInitialConditionsUUID
-----------------------------

This is similar to ``MetaDataRestartDatasetUUID``, except it's intended for tracking which initial conditions were used for a simulation.

Example:

::

    MetaDataInitialConditionsUUID   = 99f71bdf-e56d-4daf-88f6-1ecd988cbc9f

Still to be done
----------------

 * Add UUID generation to ``inits`` store it in the HDF5 output.
 * Preserve the UUID when using ``ring``.
 * Have Enzo check for the UUID in both cases.

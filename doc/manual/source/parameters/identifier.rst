Simulation Identifiers and UUIDs
--------------------------------

These parameters help to track, identify and group datasets. For reference,
`Universally Unique Identifiers
<http://en.wikipedia.org/wiki/Universally_Unique_Identifier>`_ (UUIDs) are
opaque identifiers using random 128-bit numbers, with an extremely low chance
of collision. (See :ref:`SimulationNamesAndIdentifiers` for a longer
description of these parameters.)

``MetaDataIdentifier`` (external)
    This is a character string without spaces (specifically, something
    that can be picked by "%s"), that can be defined in a parameter
    file, and will be written out in every following output, if it is
    found.
``MetaDataSimulationUUID`` (internal)
    A UUID that will be written out in all of the following outputs.
    Like ``MetaDataIdentifier``, an existing UUID will be kept, but if one
    is not found, and new one will be generated.
``MetaDataDatasetUUID`` (internal)
    A UUID created for each specific output.
``MetaDataRestartDatasetUUID`` (internal)
    If a ``MetaDataDatasetUUID`` UUID is found when the parameter file is
    read in, it will written to the following datasets. This is used to
    track simulations across restarts and parameter adjustments.
``MetaDataInitialConditionsUUID`` (internal)
    This is similar to ``MetaDataRestartDatasetUUID``, except it's used to
    track which initial conditions were used.


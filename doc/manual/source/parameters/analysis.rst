Inline Analysis
~~~~~~~~~~~~~~~

Inline Halo Finding
^^^^^^^^^^^^^^^^^^^

Enzo can find dark matter (sub)halos on the fly with a
friends-of-friends (FOF) halo finder and a subfind method,
originally written by Volker Springel. All output files will be
written in the directory FOF/.

``InlineHaloFinder`` (external)
    Set to 1 to turn on the inline halo finder. Default: 0.
``HaloFinderSubfind`` (external)
    Set to 1 to find subhalos inside each dark matter halo found in the
    friends-of-friends method. Default: 0.
``HaloFinderOutputParticleList`` (external)
    Set to 1 to output a list of particle positions and IDs for each
    (sub)halo. Written in HDF5. Default: 0.
``HaloFinderMinimumSize`` (external)
    Minimum number of particles to be considered a halo. Default: 50.
``HaloFinderLinkingLength`` (external)
    Linking length of particles when finding FOF groups. In units of
    cell width of the finest static grid, e.g. unigrid -> root cell
    width. Default: 0.1.
``HaloFinderCycleSkip`` (external)
    Find halos every N\ :sup:`th`\  top-level timestep, where N is this
    parameter. Not used if set to 0. Default: 3.
``HaloFinderTimestep`` (external)
    Find halos every dt = (this parameter). Only evaluated at each
    top-level timestep. Not used if negative. Default: -99999.0
``HaloFinderRunAfterOutput`` (external)
    When turned on, the inline halo finder is run after an output is written.  Default: 0
``HaloFinderLastTime`` (internal)
    Last time of a halo find. Default: 0.

Inline Python
^^^^^^^^^^^^^

``PythonTopGridSkip`` (external)
    How many top grid cycles should we skip between calling python at the top of the hierarchy?  Only works with python-yes in compile settings.
``PythonSubcycleSkip`` (external)
    How many subgrid cycles should we skip between calling python at the bottom of the hierarchy?
``PythonReloadScript`` (external)
    Should "user_script.py" be reloaded in between Python calls?
``NumberOfPythonCalls`` (internal)
    Internal parameter tracked by Enzo
``NumberOfPythonTopGridCalls`` (internal)
    Internal parameter tracked by Enzo
``NumberOfPythonSubcycleCalls`` (internal)
    Internal parameter tracked by Enzo

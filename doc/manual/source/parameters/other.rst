Other Parameters
~~~~~~~~~~~~~~~~

Other External Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

``huge_number`` (external)
    The largest reasonable number. Rarely used. Default: 1e+20
``tiny_number`` (external)
    A number which is smaller than all physically reasonable numbers.
    Used to prevent divergences and divide-by-zero in C++ functions.
    Modify with caution! Default: 1e-20.

    An independent analog, ``tiny``, defined in ``fortran.def``, does the same
    job for a large family of FORTRAN routines. Modification of ``tiny`` must
    be done with caution and currently requires recompiling the code, since
    ``tiny`` is not a runtime parameter.

``TimeActionParameter[#]``
    Reserved for future use.
``TimeActionRedshift[#]``
    Reserved for future use.
``TimeActionTime[#]``
    Reserved for future use.
``TimeActionType[#]``
    Reserved for future use.
``StopSteps``
    Reserved for future use
``CoolDataf0to3``
    Reserved for future use
``StageInput``
    Reserved for future use
``LocalPath``
    Reserved for future use
``GlobalPath``
    Reserved for future use

Other Internal Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

``TimeLastDataDump`` (internal)
    The code time at which the last time-based output occurred.
``TimeLastInterpolatedDataDump`` (internal)
    The code time at which the last interpolated data dump occurred.
``CycleLastDataDump`` (internal)
    The last cycle on which a cycle dump was made
``SubcycleLastDataDump`` (internal)
    The last cycle on which a subcycle dump was made
``TimeLastMovieDump`` (internal)
    The code time at which the last movie dump occurred.
``TimeLastTracerParticleDump`` (internal)
    The code time at which the last tracer particle dump occurred.
``TimeLastRestartDump``
    Reserved for future use.
``TimeLastHistoryDump``
    Reserved for future use.
``CycleLastRestartDump``
    Reserved for future use.
``CycleLastHistoryDump``
    Reserved for future use.
``InitialCPUTime``
    Reserved for future use.
``InitialCycleNumber`` (internal)
    The current cycle
``SubcycleNumber`` (internal)
    The current subcycle
``DataDumpNumber`` (internal)
    The identification number of the next output file (the 0000 part of
    the output name). This is used and incremented by both the cycle
    based and time based outputs. Default: 0
``MovieDumpNumber`` (internal)
    The identification number of the next movie output file. Default: 0
``TracerParticleDumpNumber`` (internal)
    The identification number of the next tracer particle output file. Default: 0    
``RestartDumpNumber``
    Reserved for future use.
``HistoryDumpNumber``
    Reserved for future use.
``DataLabel[#]`` (internal)
    These are printed out into the restart dump parameter file. One
    Label is produced per baryon field with the name of that baryon
    field. The same labels are used to name data sets in HDF files.
``DataUnits[#]`` 
    Reserved for future use.
``VersionNumber`` (internal)
    Sets the version number of the code which is written out to restart
    dumps.

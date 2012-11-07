Stopping Parameters
-------------------

``StopTime`` (external)
    This parameter specifies the time (in code units) when the
    calculation will halt. For cosmology simulations, this variable is
    automatically set by ``CosmologyFinalRedshift``. *No default.*
``StopCycle`` (external)
    The cycle (top grid timestep) at which the calculation stops. A
    value of zero indicates that this criterion is not be used.
    *Default: 100,000*
``StopFirstTimeAtLevel`` (external)
    Causes the simulation to immediately stop when a specified level is
    reached. Default value 0 (off), possible values are levels 1
    through maximum number of levels in a given simulation.
``NumberOfOutputsBeforeExit`` (external)
    After this many datadumps have been written, the code will exit.  If 
    set to 0 (default), this option will not be used.  Default: 0.
``StopCPUTime`` (external)
    Causes the simulation to stop if the wall time exceeds ``StopCPUTime``.
    The simulation will output if the wall time after the next
    top-level timestep will exceed ``StopCPUTime``, assuming that the wall
    time elapsed during a top-level timestep the same as the previous
    timestep. In units of seconds. Default: 2.592e6 (30 days)
``ResubmitOn`` (external)
    If set to 1, the simulation will stop if the wall time will exceed
    ``StopCPUTime`` within the next top-level timestep and run a shell
    script defined in ``ResubmitCommand`` that should resubmit the job
    for the user. Default: 0.
``ResubmitCommand`` (external)
    Filename of a shell script that creates a queuing (e.g. PBS)
    script from two arguments, the number of processors and parameter
    file.  This script is run by the root processor when stopping with
    ``ResubmitOn``. An example script can be found in
    input/resubmit.sh. Default: (null)

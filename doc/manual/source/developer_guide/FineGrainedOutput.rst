Fine Grained Output
===================

When making significant changes to Enzo that have non-local impact, such as adding a new accretion
mechanism for sink particles or face centered magnetic fields, there are many
places to introduce errors.  In order to examine the effect of changes at
specific points in the code, one can use ``ExtraOutputs``.  This run time
parameter makes a call to ``WriteAllData`` at various points in ``EvolveLevel``.
For instance, putting::

    ExtraOutput = 1
    StopCycle = 2
    MaximumRefinementLevel = 0

will cause::

    ED01_0000
    ED01_0001

to be written, along with your regular outputs.  With one level of refinement,
six outputs will be written.  The relation between output number and position is
below.  Up to 10 output points can be specified.  


Unraveling what output gets written when can be a challenge.  One technique is
to run with ``-d``, and use the following command::

    egrep "^Level||ExtraOutput" output.log

on the output log, which will show what output gets called on which level, and a string indicating
at which point in ``EvolveLevel`` it was called.

It should be noted that ``ExtraOutputs`` is not written into parameter files on
data dumps, though it can be added to restart parameter files.  This is to prevent absurd amounts of data being written to disk.
By design, this technique outputs many data dumps for each root grid timestep,
following the W cycle.
This has the added disadvantage of making the code slower, as disk access is
rarely the fastest part of any machine.   

In the code, overhead is minimized by wrapping the full function signature in a macro.  New calls can be added with::

    EXTRA_OUTPUT_MACRO(42, "After My Special Purpose")

where, of course, 42 is replaced by an integer not used by another output, and the string represents the location in the code.  It is often instructive to include this output mechanism in ``EvolveHierarchy`` as well, though this has not been done in the current checkin.



Here's a table of output number vs. position in ``EvolveLevel``.  Please refer to the :ref:`FlowChart` 
to understand each entry. The non-continuity represents some outputs that will be introduced when 
MHDCT is merged, but not relevant for pure hydro.

===== ===================================
Index Position in ``EvolveLevel``
----- -----------------------------------
1     Before time loop
2     After SolveHydroEquations grid loop
25     After SetBoundaryConditions
3     Before UpdateFromFinerGrids
4     After UpdateFromFinerGrids
6     After the time loop
===== ===================================

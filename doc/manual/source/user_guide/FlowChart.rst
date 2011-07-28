Enzo Flow Chart, Source Browser
===============================

`Here's a cartoon of Enzo. <http://lca.ucsd.edu/software/enzo/v1.5/flowchart/>`_
This was written as a first look as the details of how enzo works.
Black arrows indicate further flow charts. Grey boxes (usually)
indicate direct links to the source code.

No guarantees are made regarding the correctness of this flowchart
-- it's meant to help get a basic understanding of the flow of Enzo
before extensive code modifications.
`Also see the Enzo Source Browser. <http://lca.ucsd.edu/software/enzo/v1.0.1/source_browser/>`_
This is a second attempt at the same thing in a more dynamic way.
It allows one to (in principle) see all the routines called from a
function, in order, and jump to the source showing the call. It
also allows you to see a reverse call stack of every routine that
calls a particular function. It runs relatively well on the newer
versions of Firefox. Older browsers have a harder time with the
Javascript. This is also somewhat buggy, for a host of reasons.

Flow Chart Errors
-----------------

As mentioned above, this was written with an old version of the
code and has not been made perfect. Issues may remain.

Return CPU Time
~~~~~~~~~~~~~~~

This used to be an Enzo routine. Now, it's an MPI call, MPI_Wtime



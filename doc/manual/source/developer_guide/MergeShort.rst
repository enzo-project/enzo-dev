Branching and Merging: The Short Version
========================================

This is short version of branching and merging for experienced
developers with no long term memory. (And can't keep trunk/devel
and devel/trunk straight.)

Short, single line copy and commit:

::

    svn copy -m "branch for testing NewBaryon, NewProblemGen" http://lca.ucsd.edu/svn/Enzo/devel/trunk http://lca.ucsd.edu/svn/Enzo/devel/branches/TestingForDocument

Then check out:

::

    svn co http://lca.ucsd.edu/svn/Enzo/devel/branches/TestingForDocument

Then merge back (into trunk sandbox)

::

    svn merge -r1157:HEAD   http://lca.ucsd.edu/svn/Enzo/devel/branches/TestingForDocument .

Make sure it builds and whatnot.

Then check in.



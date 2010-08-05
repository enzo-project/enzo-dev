Diffing Copies
==============

Set Up: Add and remove a file
-----------------------------

::

    cable:~/Projects/tmp/Enzo rpwagner$ emacs -nw TESTFILE
    cable:~/Projects/tmp/Enzo rpwagner$ svn add TESTFILE
    A         TESTFILE
    cable:~/Projects/tmp/Enzo rpwagner$ svn commit TESTFILE 
    Adding         TESTFILE
    Transmitting file data .
    Committed revision 279.
    cable:~/Projects/tmp/Enzo rpwagner$ svn rm TESTFILE
    D         TESTFILE
    cable:~/Projects/tmp/Enzo rpwagner$ svn commit TESTFILE 
    Deleting       TESTFILE
    
    Committed revision 280.
    cable:~/Projects/tmp/Enzo rpwagner$

Revision vs. Working Directory
------------------------------

Diff the local directory against the previous revision, using the
revision number, and see the change.

::

    cable:~/Projects/tmp/Enzo rpwagner$ svn diff -r 279
    Index: TESTFILE
    ===================================================================
    --- TESTFILE    (revision 279)
    +++ TESTFILE    (working copy)
    @@ -1 +0,0 @@
    -Repeating add/remove file.
    \ No newline at end of file

For fun, look two revisions back, and see that it's identical to
the working copy.

::

    cable:~/Projects/tmp/Enzo rpwagner$ svn diff -r 278
    cable:~/Projects/tmp/Enzo rpwagner$ 

Diffing Between Copies
----------------------

**NOTE:** It turns out tags and branches in Subversion are just
copies. So when cvs2svn imported the CVS tags, all it did was
create a copy for each tag.

To use svn diff to compare different copies, use the --old and
--new flags.

Copy vs. Working Directory
--------------------------

This compares ones of the "tags" (read "copy") made at import,
against the working directory.

::

    cable:~/Projects/tmp/Enzo rpwagner$ cd doc
    cable:~/Projects/tmp/Enzo/doc rpwagner$ svn diff --old http://lca.ucsd.edu/svn/Enzo/tags/v1_3_4/Enzo/doc --new . | more
    Index: README.v130-rad-hydro
    ===================================================================
    --- README.v130-rad-hydro       (.../http://lca.ucsd.edu/svn/Enzo/tags/v1_3_4/Enzo/doc)   (revision 280)
    +++ README.v130-rad-hydro       (working copy)
    @@ -328,9 +328,13 @@
        http://lca.ucsd.edu/projects/LLNL/development/imp_rad_system-v4.pdf
     
        The updated self-gravity solver solves the basic Poisson equation
    -   for the gravitational potential, using the standard Inexact Newton
    -   solver framework.
    +   for the gravitational potential.

Copy vs. Copy
-------------

This compares two different tags (copies) against each other.

::

    cable:~/Projects/tmp/Enzo/doc rpwagner$ svn diff --old http://lca.ucsd.edu/svn/Enzo/tags/v1_3_0/Enzo/doc \
      --new http://lca.ucsd.edu/svn/Enzo/tags/v1_3_5/Enzo/doc | more
    Index: README.v130-rad-hydro
    ===================================================================
    --- README.v130-rad-hydro       (.../v1_3_0/Enzo/doc)   (revision 280)
    +++ README.v130-rad-hydro       (.../v1_3_5/Enzo/doc)   (revision 280)
    @@ -328,9 +328,13 @@
        http://lca.ucsd.edu/projects/LLNL/development/imp_rad_system-v4.pdf
     
        The updated self-gravity solver solves the basic Poisson equation
    -   for the gravitational potential, using the standard Inexact Newton
    -   solver framework.
    +   for the gravitational potential.

More Diff
=========

Local vs Repository
-------------------

If you have multiple "branches" of a code on a machine, you can
compare either of the working copies with either of the repository
copies. (this only makes a difference if you've changed something
in one of the copies.)

Brief Diffing
-------------

If you want to just see the *files* that differ between two
Subversion Revisions, you're SOL. Subversion doesn't have that
functionality. The dumbest workaround I came up with is piping the
output through grep:

::

    ds011:~/Microscope/trunk>svn diff -r20:99 | grep "^Index"



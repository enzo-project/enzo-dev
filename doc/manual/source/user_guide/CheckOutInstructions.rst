Subversion Check Out Instructions
=================================

The latest public version of Enzo is available for anonymous
checkout using ` Subversion <http://subversion.tigris.org/>`_. This
is where bug fixes and new features will appear between releases.

You also browse the source tree, either through the
` default Subversion HTTP interface <http://mngrid.ucsd.edu/svn/Enzo/public>`_,
or the nicer [browser:public Trac browser].

*A hint*: Before you try to build Enzo, you might want to make sure
you meet the
`compilation requirements? </wiki/Devel/UserGuide/CompilationRequirements>`_.

Subversion Clients
------------------

To check out a local copy of the Enzo source, you need a
` Subversion <http://subversion.tigris.org/>`_ client. These are
available as part of all recent Linux distributions; for other
operating systems (OS X, AIX, etc.), binaries are available
(there's a list of third-party clients on the
` Subversion front page <http://subversion.tigris.org/>`_), or the
client can be built from source. GUI clients are available, but
these instructions assume you're using a command line client.

To see if you already have the client installed, use which.

::

    $ which svn
    /usr/local/bin/svn

If you don't have the client, and are a shared shared resource (a
cluster, or supercomputer), please ask the system administrator to
install Subversion for everyone. This will make things easier for
the next user who comes along.

Getting a Copy
--------------

FYI: Checking out Enzo will also get you a copy
` YT <http://yt.enzotools.org/>`_, the
` Python <http://www.python.org>`_ based analysis toolkit. Check
the ` YT website <http://yt.enzotools.org/>`_ for instructions on
compiling and using ` YT <http://yt.enzotools.org/>`_.

Once you have the client, you can use it *checkout* a local copy.

::

    $ svn checkout http://mngrid.ucsd.edu/svn/Enzo/public/trunk enzo
    A    enzo/configure
    A    enzo/doc
    A    enzo/doc/flowchart
    ...
    A    enzo/src/enzo/Grid_FastSiblingLocatorFindSiblings.C
    A    enzo/src/enzo/AnalyzeClusters.h
    A    enzo/src/Makefile
    A    enzo/bin
    
    Fetching external item into 'enzo/src/yt'
    A    enzo/src/yt/LICENSE.txt
    A    enzo/src/yt/epydoc.cfg
    A    enzo/src/yt/tests
    ...
    A    enzo/src/yt/examples/test_parallel_projection.py
    A    enzo/src/yt/setup.cfg
     U   enzo/src/yt
    Checked out external at revision 731.
    
    Checked out revision 1761.

And now you have a copy of the latest public version. Time to work
on `building Enzo? </wiki/Devel/UserGuide/BuildingEnzo>`_.

Updating
--------

To update your local copy, you can use svn to pull down only the
latest changes.

*Note*: If you've modified any files in your copy, this will merge
changes from the trunk in to your working copy, which may generate
conflicts. If you're doing development on Enzo itself, you may want
to check the [log:public/trunk revision log] before doing an
update.

::

    $ cd enzo/
    $ svn update
    A    README
    
    Fetching external item into 'src/yt'
    Updated external to revision 731.
    
    Updated to revision 1762.

Now, you can do a make clean; make and get back to work.



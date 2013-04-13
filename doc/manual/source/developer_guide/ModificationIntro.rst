.. _enzo_modification:

Introduction to Enzo Modification
=================================

.. note:: This is not a comprehensive document, but it does cover some of the
          grounds of modifying Enzo.  Please don't hesitate to email the users'
          mailing list with any further questions about Enzo, Mercurial, or how
          to write and execute new test problems.

If this is the first time you've opened the hood to Enzo, welcome.  If you're
an old hand and have already added new physics to it, welcome back.

Enzo is an extremely powerful piece of software, but by no means a complete
representation of the observable universe. It's quite likely that there will be
some piece of physics that you'll want to model, and these span a broad range
of software complexities. In all cases, whether it's a mildly invasive change
such as a new background heating model or extremely invasive like adding
relativistic non-neutral multi-fluid plasma physics, we strongly recommend
taking advantage of some basic tools. These are outlined in the sections that
follow.  These tools prevent the developer from breaking existing features
(which is far easier than one would expect), keeping track of your changes, and
sharing those changes with others. We strongly recommend you start with getting
LCATest running before you start programming, so mistakes can be caught early.

Additionally in the Tutorials section you'll see a pair of flow chart tools
that are intended as educational tools, and several descriptions on how to
actually add various components to the code.  It is intended that these will be
at least read in order, as doing complex things with the code require the
ability to do the simpler things.

We are very happy to accept patches, features, and bugfixes from any member of
the community!  Enzo is developed using mercurial, primarily because it enables
very easy and straightforward submission of changesets.  We're eager to hear
from you, and if you are developing Enzo, please subscribe to the users'
mailing list:

http://groups.google.com/group/enzo-users

This document describes how to use Mercurial to make changes to Enzo, how to
send those changes upstream, and how to navigate the Enzo source tree.

.. highlight:: none

Mercurial Introduction
----------------------

If you're new to Mercurial, these three resources are pretty great for learning
the ins and outs:

   * http://hginit.com
   * http://hgbook.red-bean.com/read/
   * http://mercurial.selenic.com/

The major difference between Mercurial (and other distributed version control
systems) and centralized version control systems (like CVS, RCS, SVN) is that
of the directed acyclic graph (DAG).  Rather than having a single timeline of
modifications, Mercurial (or "hg") can have multiple, independent streams of
development.

There are a few concepts in Mercurial to take note of:

Changesets
   Every point in the history of the code is referred to as a changeset.  These
   are specific states of the code, which can be recovered at any time in *any*
   checkout of the repository.  These are analagous to revisions in Subversion.
Children
   If a changeset has changesets that were created from its state, those are
   called children.  A changeset can have many children; this is how the graph
   of development branches.
Heads
   Every changeset that has no children is called a head.
Branches
   Every time the DAG branches, these are branches.  Enzo also uses "named
   branches," where the branches have specific identifiers that refer to the
   feature under development or some other characteristic of a line of
   development.

When you check out the Enzo repository, you receive a full and complete copy of
the entire history of that repository; you can update between revisions at
will without ever touching the network again.  This allows not only for
network-disconnected development, but it also means that if you are creating
some new feature on top of Enzo you can (and should!) conduct local version
control on your development.  Until you choose explicitly to share changes,
they will remain private to your checkout of the repository.

Enzo Source Trees
-----------------

Enzo has two primary repositories, the "stable" repository which is curated and
carefully modified, and the "development" repository which is where active
development occurs.  Please note that while we test and verify the results of
the "stable" repository, the "unstable" repository is not guaranteed to be
tested, verified, or even to provide correct answers.

.. note:: The "stable" Enzo source tree is *not* for general development.  If
   you want to contribute to Enzo, make your changes to a fork of the
   development repository!

To conceptually -- and technically! -- separate these two repositories, they
also live in different places.  We keep both the stable repository 
and the development repository at BitBucket.  Enzo is (as of 2.2) developed in
a relatively simple fashion:

  #. On BitBucket, developers "fork" the primary development repository.
  #. When a piece of work is ready to be shared, a "pull request" is issued.
     This notifies the current set of Enzo curators that a new feature has been
     suggested for inclusion.
  #. After these features have been accepted, they are pulled into the
     development branch.  New features will be aggregated into patch
     releases on the "stable" branch.
  #. When a new patch release is issued, the current development branch is
     pushed to the "enzo-stable" repository on Bitbucket.

The idea here is that there is a double firewall: the development repository is
very high-cadence and with high-turnover, but the stable repository is much
slower, more carefully curated, and inclusions in it are well-tested.

 * Stable code lives at: https://bitbucket.org/enzo/enzo-stable
 * Development code lives at: http://bitbucket.org/enzo/enzo-dev

How To Share Changes
--------------------

Sharing your changes to Enzo is easy with Mercurial and the BitBucket
repository.

Go here:

http://bitbucket.org/enzo/enzo-dev/fork

Now, clone your new repository.  Make your changes there.  Now go back and
issue a pull request.  For instance, you might do something like this:

 #. Clone Enzo, make a few changes, commit them, and decide you want to share.
 #. Fork the main enzo repository at that link.
 #. Now, edit ``.hg/hgrc`` to add a new path, and push to that path.
 #. Go to the BitBucket URL for your new repository and click "Pull Request".
    Fill it out, including a summary of your changes, and then submit.  It will
    get evaluted -- and it might not get accepted right away, but the response
    will definitely include comments and suggestions.

That's it!  If you run into any problems, drop us a line on the `Enzo Users'
Mailing List <http://groups.google.com/group/enzo-users>`_.

How To Use Branching
--------------------

.. warning:: In most cases, you do *not* need to make a new named branch!  Do
   so with care, as it lives forever.

If you are planning on making a large change to the code base that may not be
ready for many, many commits, or if you are planning on breaking some
functionality and rewriting it, you can create a new named branch.  You can
mark the current repository as a new named branch by executing: ::

   $ hg branch new_feature_name

To merge changes in from another branch, you would execute: ::

   $ hg merge some_other_branch

Note also that you can use revision specifiers instead of "some_other_branch".
When you are ready to merge back into the main branch, execute this process: ::

   $ hg merge name_of_main_branch
   $ hg commit --close-branch
   $ hg up -C name_of_main_branch
   $ hg merge name_of_feature_branch
   $ hg commit

When you execute the merge you may have to resolve conflicts.  Once you resolve
conflicts in a file, you can mark it as resolved by doing: ::

   $ hg resolve -m path/to/conflicting/file.py

Please be careful when resolving conflicts in files.

Once your branch has been merged in, mark it as closed on the wiki page.

The Patch Directory
--------------------

If you are experimenting with a code change or just debugging, then
the patch directory, found in the top level of your Enzo directory,
may be of use. Files put in here are compiled in preference to those
in ``/src/enzo``, so you can implement changes without overwriting the
original code. To use this feature, run ``make`` from inside
``/patch``. You may need to add ``-I../src/enzo`` to the
``MACH_INCLUDES`` line of your machine makefile
(e.g. ``Make.mach.triton``) to ensure the .h files are found when compiling.

As an example, suppose you wish to check the first few values of the acceleration field as Enzo runs through ``EvolveLevel.C``. Copy ``EvolveLevel.C`` from ``/src/enzo`` into ``/patch`` and put the appropriate print statements throughout that copy of the routine. Then recompile Enzo from inside the patch directory. When you no longer want those changes, simply delete EvolveLevel.C from ``/patch`` and the next compile of the code will revert to using the original ``/src/enzo/EvolveLevel.C``. If you make adjustments you wish to keep, just copy the patch version of the code into ``/src/enzo`` to replace the original.


How To Include Tests
--------------------

If you have added any new functionality, you should add it as a test in the
directory tree ``run/`` under the (possibly new!) appropriate directory.  Your
test file should consist of:

 * A parameter file, ending in the extension ``.enzo``
 * A file of ``notes.txt``, describing the problem file, the expected results,
   and how to verify correctness
 * A test file, using the yt extension ``enzo_test``, which verifies
   correctness.  (For more information on this, see some of the example test
   files.)
 * (optional) Scripts to plot the output of the new parameter file.

Please drop a line to the mailing list if you run into any problems!

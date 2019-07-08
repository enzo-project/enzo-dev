.. _enzo_modification:

Introduction to Enzo Modification
=================================

.. note:: This is not a comprehensive document, but it does cover some of the 
          grounds of modifying Enzo. Please don't hesitate to email the users' 
          mailing list with any further questions about Enzo, Git, or how
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

.. _contributing-code:

How to Develop Enzo
===================

Enzo is a community project!

We are very happy to accept patches, features, and bugfixes from any member of
the community!  Enzo is developed using Git, primarily because of how well
it enables open-source, community contribution. We're eager to hear from you.

.. note:: If you are already familiar with Git and `GitHub <https://github.com>`_,
   the best way to contribute is to fork the `main Enzo repository
   <https://github.com/enzo-project/enzo-dev.git>`__, make your changes, push them
   to your fork, and issue a pull request. The rest of this document is just an
   explanation of how to do that. Enzo was previously hosted on BitBucket and 
   maintained with mercurial but has since transitioned to git and GitHub. 
   Please do not fork or clone the repository on BitBucket. The code on BitBucket
   is no longer maintained.

Keep in touch, and happy hacking!

.. _open-issues:

Open Issues
-----------

If you're interested in participating in Enzo development, take a look at the
`issue tracker on GitHub <https://github.com/enzo-project/enzo-dev/issues>`_.
If you are encountering a bug that is not already tracked there, please `open a
new issue <https://github.com/enzo-project/enzo-dev/issues/new>`__.

Contributing to Enzo with Git and Github
----------------------------------------

We provide a brief introduction to submitting changes here.  We encourage
contributions from any user. If you are new to Git and/or GitHub, there are
excellent guides available at `guides.github.com <https://guides.github.com/>`_,
specifically the `Git Handbook
<https://guides.github.com/introduction/git-handbook/>`__, and the `GitHub
Hello World <https://guides.github.com/activities/hello-world/>`__. We are also
happy to provide guidance on the mailing list or in our slack channel.

Licensing
+++++++++

Enzo is under a BSD-like license.

All contributed code must be BSD-compatible.  If you'd rather not license in
this manner, but still want to contribute, please consider creating an external
package, which we'll happily link to in the Enzo documentation.

How To Get The Source Code For Editing
++++++++++++++++++++++++++++++++++++++

Enzo is hosted on GitHub. In order to modify the source code for Enzo,
we ask that you make a "fork" of the main Enzo repository on GitHub.  A
fork is simply an exact copy of the main repository (along with its history)
that you will now own and can make modifications as you please.  You can create
a personal fork by visiting the Enzo GitHub webpage at
https://github.com/enzo-project/enzo-dev/.  After logging in, you should see an
option near the top right labeled "fork".  You now have a forked copy of
the Enzo repository for your own personal modification.

This forked copy exists on GitHub under your username, so in order to access
it locally, follow the instructions at the top of that webpage for that
forked repository:

.. code-block:: bash

   $ git clone http://github.com/<USER>/<REPOSITORY_NAME>

This downloads that new forked repository to your local machine, so that you can
access it, read it, make modifications, etc.  It will put the repository in a
local directory of the same name as the repository in the current working
directory.

.. code-block:: bash

   $ cd enzo_dev

Verify that you are on the master branch of Enzo by running:

.. code-block:: bash

   $ git branch

If you're not on the master branch, you can get to it with:

.. code-block:: bash

   $ git checkout master

You can see any past state of the code by using the git log command.
For example, the following command would show you the last 5 revisions
(modifications to the code) that were submitted to that repository.

.. code-block:: bash

   $ git log -n 5

Using the revision specifier (the number or hash identifier next to each
changeset), you can update the local repository to any past state of the
code (a previous changeset or version) by executing the command:

.. code-block:: bash

   $ git checkout revision_specifier

.. _sharing-changes:

Making and Sharing Changes
--------------------------

The simplest way to submit changes to Enzo is to do the following:

#. Fork the main repository.
#. Clone your fork.
#. Make some changes and commit them.
#. Push the changesets to your fork.
#. Issue a pull request.

Here's a more detailed flowchart of how to submit changes.

#. Fork Enzo on GitHub.  (This step only has to be done once.)  You can do
   this by clicking on the **fork** button in the top-right corner of `the main
   repository <https://github.com/enzo-project/enzo-dev>`__.
#. Create a new branch in which to make your changes by doing ``git
   checkout -b <new branch name>``. This will make it easy to move back and
   forth between the main branch of the code and your changes.
#. Edit the source file you are interested in and test your changes.
#. Use ``git add <files>`` to stage files to be committed.
#. Commit your changes with ``git commit``. This will open a text editor so you
   can write a commit message. To add your message inline, do
   ``git commit -m "<commit message>"``. You can list specific file to be
   committed.
#. Remember that this is a large development effort and to keep the code
   accessible to everyone, good documentation is a must.  Add in source code
   comments for what you are doing.  Add documentation to the appropriate
   section of the online docs so that people other than yourself know how
   to use your new code.
#. If your changes include new functionality or cover an untested area of the
   code, add a test. Commit these changes as well.
#. Push your changes to your new fork using the command::

      $ git push origin <branch name>

   .. note::
     Note that the above approach uses HTTPS as the transfer protocol
     between your machine and GitHub.  If you prefer to use SSH - or
     perhaps you're behind a proxy that doesn't play well with SSL via
     HTTPS - you may want to set up an `SSH key
     <https://help.github.com/articles/connecting-to-github-with-ssh/>`__
     on GitHub.  Then, you use
     the syntax ``ssh://git@github.com/<USER>/enzo_dev``, or equivalent, in
     place of ``https://github.com/<USER>/enzo_dev`` in git commands.
     For consistency, all commands we list in this document will use the HTTPS
     protocol.

#. Issue a pull request by going to the main repository and clicking on the
   green button that says **Compare & pull request**. This will open up a page
   that will allow you to enter a description of the changes to be merged. Once
   submitted, a series of automated tests will run and their status will be
   reported on the pull request page.

During the course of your pull request you may be asked to make changes.  These
changes may be related to style issues, correctness issues, or requesting
tests.  The process for responding to pull request code review is relatively
straightforward.

#. Make requested changes, or leave a comment on the pull request page on
   GitHub indicating why you don't think they should be made.
#. Commit those changes to your local repository.
#. Push the changes to your fork::

      $ git push origin <branch name>

#. Your pull request will be automatically updated.

Once your pull request has been accepted, you can safely delete your
branch::

      $ git branch --delete <branch name>

Updating Your Branch
++++++++++++++++++++

If your branch or pull request has been open for some time, it may be useful
to keep it up to date with the latest changes from the main repository. This
can be done by `rebasing your changes <https://git-scm.com/docs/git-rebase>`__.
Before doing this, you will need to be able to pull the latest changes from
the main repository.

#. Add the main repository as a remote::

      $ git remote add enzo_dev https://github.com/enzo-project/enzo-dev

   You can verify that it has been added by doing ``git remote -v``. This
   only needs to be done once.

#. Go back to the master branch and pull the changes::

      $ git checkout master
      $ git pull enzo_dev master

#. Return to your branch and rebase your changes onto the head of the master
   branch::

      $ git checkout <branch name>
      $ git rebase master

This should go smoothly unless changes have been made to the same lines in
the source, in which case you will need to fix conflicts. After rebasing,
you will get an error when trying to push your branch to your fork. This is
because you have changed the order of commits and git does not like that.
In this case, you will need to add "-f" to your push command to force
the changes to be accepted.::

      $ git push -f origin <branch name>

If you envision a decent number of conflicts, instead of rebasing, it may
be better to do::

      $ git merge master

commit the merge, then push the changes to your fork.


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


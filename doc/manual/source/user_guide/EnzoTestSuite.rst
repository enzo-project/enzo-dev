.. _EnzoTestSuite:

Enzo Test Suite
===============

The Enzo test suite is a set of tools whose purpose is to perform
regression tests on the Enzo codebase, in order to help developers
discover bugs that they have introducted, to verify that the code is
producing correct results on new computer systems and/or compilers,
and, more generally, to demonstrate that Enzo is behaving as expected
under a wide variety of conditions.

What's in the test suite?
-------------------------

The suite is composed of a large number of individual test problems
that are designed to span the range of physics and dimensionalities
that are accessible using the Enzo test code, both separately and in
various permutations.  Tests can be selected based on a variety of
criteria, including (but not limited to) the physics included, the
estimated runtime of the test, and the dimensionality.  For
convenience, three pre-created, overlapping sets of tests are
provided:

1.  The "quick suite" (``--quicksuite=True``).  This is composed of
small calculations that test critical physics packages both
alone and in combination.  The intent of this package is to be run
relatively frequently (multiple times a day) to ensure that bugs have
not been introduced during the code development process.  All runs 
in the quick suite use no more than a single processor.  The total 
run time should be about 25 minutes.  The gold standard results for 
the quick suite alone can be downloaded 
`here <http://enzo-project.org/tests/gold_standard_quick.tar.gz>`_.

2.  The "push suite" (``--pushsuite=True``).  This is a slightly 
large set of tests, encompassing all of the quick suite and 
some additional larger simulations that test a wider variety of physics 
modules.  The intent of this package is to provide a thorough validation 
of the code prior to changes being pushed to the main repository.  The 
total run time is roughly 90 minutes and all simulations use only a single 
processor.  The gold standard results for the push suite can be downloaded 
`here <http://enzo-project.org/tests/gold_standard_push.tar.gz>`_.

3.  The "full suite" (``--fullsuite=True``).  This encompasses essentially 
all of test simulations contained within the run directory.  This suite 
provides the most rigorous possible validation of the code in many different 
situations, and is intended to be run prior to major changes being pushed 
to the stable branch of the code.  A small number of simulations in the full 
suite are designed to be run on 2 processors and will take multiple hours to 
complete.  The total run time is roughly 36 hours.  The gold standard results
for the full suite can be downloaded 
`here <http://enzo-project.org/tests/gold_standard_full.tar.gz>`_.

How to run the test suite
-------------------------

The Enzo test suite is run within the ``run/`` subdirectory of the
Enzo source distribution, using the ``test_runner.py`` file.  To
run the test suite, follow these instructions:

1.  Before running the test suite, you should download the "gold
standard" results for the 
`quick <http://enzo-project.org/tests/gold_standard_quick.tar.gz>`_, 
`push <http://enzo-project.org/tests/gold_standard_push.tar.gz>`_, or 
`full <http://enzo-project.org/tests/gold_standard_full.tar.gz>`_ 
suites and untar that file into a convenient directory.

2.  Compile Enzo.  The gold standard calculations use the default 
compiler settings that can be restored with ``make default``.  
If you use significantly different compilation options
(higher-level optimization in particular) you may see somewhat
different outputs that will result in failed tests.

3.  Go into the ``run/`` subdirectory in the Enzo repository and
type the following command:

::

    ./test_runner.py --quicksuite=True  --compare-dir=/path/to/gold_standard \
            --output-dir=/enzo/test/directory

In this comand, ``--quicksuite=True`` instructs the test runner to
use the quick suite (other possible keyboards here are
``--pushsuite=True`` and ``--fullsuite=True``).
``--output-dir=/enzo/test/directory`` instructs the test runner to
write output to the user-specified directory, and
``--compare-dir=/path/to/gold_standard`` instructs the test runner
to use the set of data files in the listed directory as a gold
standard for comparison.  It is also possible to choose sets of tests
that are sorted by dimensionality, physics modules, runtime, number of
processors required, and other criteria.  A single named test can be run 
by giving ``--name=<name of test>``.  Type ``./test_runner.py
--help`` for a more complete listing.


How to add a new test to the library
------------------------------------

It is hoped that any newly-created or revised physics module will be
accompanied by one or more test problems, which will ensure the
continued correctness of the code.  This sub-section explains the
structure of the test problem system as well as how to add a new test
problem to the library.

Test problems are contained within the ``run/`` directory in the
Enzo repository.  This subdirectory contains a tree of directories
where test problems are arranged by the primary physics used in that
problem (e.g., Cooling, Hydro, MHD).  These directories may be further
broken down into sub-directories (Hydro is broken into Hydro-1D,
Hydro-2D, and Hydro-3D), and finally into individual directories
containing single problems.  A given directory contains, at minimum,
the Enzo parameter file (having extension ``.enzo``, described in
detail elsewhere in the manual) and the Enzo test suite parameter file
(with extension ``.enzotest``).  The latter contains a set of
parameters that specify the properties of the test.  Consider the test
suite parameter file for InteractingBlastWaves, which can be found in the
``run/Hydro/Hydro-1D/InteractingBlastWavest`` directory:

::

    name = 'InteractingBlastWaves'
    answer_testing_script = None
    nprocs = 1
    runtime = 'short'
    critical = True
    cadence = 'nightly'
    hydro = True
    gravity = False
    dimensionality = 1
    max_time_minutes = 1

This allows the user to specify the dimensionality, physics used, the
runtime (both in terms of 'short', 'medium', and 'long' calculations,
and also in terms of an actual wall clock time), and whether the test
problem is critical (i.e., tests a fundamental piece of the code) or
not.  A general rule for choosing the runtime value is 'short' for runs 
taking less than 5 minutes, 'medium' for run taking between 5 and 30 minutes, 
and 'long' for runs taking more than 30 minutes.  A full listing of options 
can be found in the ``run/README`` file.

Once you have created a new problem type in Enzo and thoroughly
documented the parameters in the Enzo parameter list, you should
follow these steps to add it as a test problem:

1.  Create a new subdirectory in the appropriate place in the
``run/`` directory.  If your test problem uses multiple pieces of
physics, put it under the most relevant one.

2.  Add an Enzo parameter file, ending in the extension ``.enzo``,
for your test problem to that subdirectory.

3.  Add an Enzo test suite parameter file, ending in the extension
``.enzotest``.  In that file, add any relevant parameters (as
described in the ``run/README`` file).

4.  Create a "gold standard" set of data for your test problem, by
running with the default compile options. Contact Britton Smith 
(brittonsmith@gmail.com) and arrange 
to send him this data.  Please try to minimize the quantity of data
generated by your calculation by only writing out data at the end of
the calculation, not during the interim (unless evolution of a
quantity or quantities is important).

If you want to examine the output of your test problem for something
specific, you can optionally add a script that is indicated by the
``answer_testing_script`` parameter.  Look in the directory
``run/Hydro/Hydro-3D/RotatingCylinder`` for an example of how this
is done.

Congratulations, you've created a new test problem!


What to do if you fix a bug in Enzo
-----------------------------------

It's inevitable that bugs will be found in Enzo, and that some of
those bugs will affect the actual simulation results (and thus the
test problems used in the problem suite).  If you fix a bug that
results in a change to some or all of the test problems, the gold
standard solutions will need to be updated.  Here is the procedure for
doing so:

1.  Run the "push suite" of test problems (``--pushsuite=True``)
for your newly-revised version of Enzo, and determine which test
problems now fail.

2.  Visually inspect the failed solutions, to ensure that your new
version is actually producing the correct results!

3.  Email the enzo-developers mailing list at
enzo-dev@googlegroups.com to explain your bug fix, and to show the
results of the now-failing test problems.

4.  Once the denizens of the mailing list concur that you have
correctly solved the bug, create a new set of gold standard test
problem datasets, following the instructions in the next section.

5.  After these datasets are created, send the new gold standard
datasets to Britton Smith (brittonsmith@gmail.com), who will update
the gold standards.

6.  Push your Enzo changes to the repository.


How to create a new set of reference calculations
-------------------------------------------------

It may be necessary for you to generate a set of reference
calculations for some reason.  If so, here is how you do this.

1.  First, build Enzo using the default set of compile options.  
Type ``make default`` to restore the defaults.  You will 
now have an enzo binary in the ``src/enzo`` directory.

2.  Go into the ``run/`` directory and call test_runner.py without the ``--compare-dir`` directory.  If you
are have multiple Enzo repositories, you can specify the one you want:

::

    ./test_runner.py --repo=/path/to/desired/enzo/repo \
         --output-dir=/path/to/new/reference/directory

Note that you should only use the top-level directory in the
repository, not src/enzo, and if you simply want to use the current
repository (that is, the one your run directory is located in) you can
leave out the ``--repo`` option.  Once this step is completed, you should
have a full set of test problems.

3.  If you then want to compare against this set of test problems, use
the following command:

::

    ./test_runner.py --repo=/path/to/desired/enzo/repo  \
         --compare-dir=/path/to/new/reference/directory \
         --output-dir=/path/to/output/directory





.. _EnzoTestSuite:

Enzo Test Suite
===============

The Enzo test suite is a set of tools whose purpose is to perform
regression tests on the Enzo codebase, in order to help developers
discover bugs that they have introducted, to verify that the code is
producing correct results on new computer systems and/or compilers, and,
more generally, to demonstrate that Enzo is behaving as expected under
a wide variety of conditions.

What's in the test suite?
--------------------

The suite is composed of a large number of individual test problems that are
designed to span the range of physics and dimensionalities that are
accessible using the Enzo test code, both separately and in various
permutations.  Tests can be selected based on a variety of criteria,
including (but not limited to) the physics included, the estimated
runtime of the test, and the dimensionality.  For convenience, 
three pre-created, overlapping sets of tests are provided:

1.  The "quick suite" (keyword: quick).  This set of tests should run
in a few minutes on a single core of a laptop.  It is composed of
one-dimensional calculations that test critical physics packages both
alone and in combination.  The intent of this package is to be run
relatively frequently (multiple times a day) to ensure that bugs have
not been introduced during the code development process.  A
full list of the test problems included in the "quick suite" can be
found at PUT LINK HERE.

2.  The "push suite" (keyword: push).  This set of tests should run in
roughly an hour or two on a single core of a laptop or desktop
machine.  It is composed of one-, two- and three-dimensional
calculations that test a wider variety of physics modules, both alone
and in combination, than the "quick suite" described above.  The
intent of this package is to provide a thorough validation of the code
prior to changes being pushed to the Google Code mercurial
repository.  A
full list of the test problems included in the "push suite" can be
found at PUT LINK HERE.

3.  The "full suite" (keyword: full).  This set of tests should run in
no more than 24 hours on 8-16 cores of a Linux cluster, and includes a
variety of 3D, multiphysics tests that are not in the "quick" and
"push" suites.  This suite provides the most rigorous possible
validation of the code in many different situations, and is intended
to be run prior to major changes being pushed to the Google Code
repository and prior to public releases of the code.   A
full list of the test problems included in the "full suite" can be
found at PUT LINK HERE.


How to run the test suite
------------------

The Enzo test suite is run within the :file:`run/` subdirectory of the
Enzo source distribution, using the :file:`test_runner.py` file.  To
run the test suite, follow these instructions:

1.  Before running the test suite, you should download the "gold standard"
datasets from INSERT URL HERE, and untar that file into a convenient
directory.

2.  Compile Enzo.  The gold standard calculations use opt-debug and
64-bit precision everywhere (:file:`make opt-debug`, :file:`make
precision-64`, :file:`make particles-64`, and :file:`make
integers-64`).

3.  Go into the :file:`run/` subdirectory in the Enzo repository and type the following
command:

::

    ./test_runner.py --suite=quick -o <enzo/test/directory> --compare-dir=<path/to/gold_standard>

In this comand, :file:`--suite=quick` instructs the test runner to use
the quick suite (other possible keyboards here are 'push' and 'full').
:file:`-o <enzo/test/directory>` instructs the test runner to write
output to the user-specified directory, and
:file:`--compare-dir=<path/to/gold_standard>` instructs the test
runner to use the set of data files in the listed directory as a gold
standard for comparison. It is also possible to choose sets of tests
that are sorted by dimensionality, physics modules, runtime, number of
processors required, and other criteria.  Type :file:`./test_runner.py
--help` for a more complete listing.


How to add a new test to the library
----------------------------

It is hoped that any new physics module will be accompanied by one or
more new test problems, which will ensure the continued correctness of
the code.  This sub-section explains the structure of the test problem
system as well as how to add a new test problem to the library.

Test problems are contained within the :file:`run/` subdirectory in
the Enzo repository.  This subdirectory contains a tree of further
directories, where test problems are arranged by the primary physics
used by a problem (e.g., Cooling, Hydro, MHD).  These directories may
be further broken down into sub-categories (Hydro is broken into
Hydro-1D, Hydro-2D, and Hydro-3D), and finally into individual
directories containing single problems.  A given directory contains,
at minimum, the Enzo parameter file (having extension  :file:`.enzo`,
described in detail elsewhere in the manual)
and the Enzo test suite parameter file (with extension
:file:`.enzotest`).  The latter contains a set of parameters that
specify the properties of the test.  Consider the test suite parameter
file for 
InteractingBlastWaves, found at
:file:`run/Hydro/Hydro-1D/InteractingBlastWaves/InteractingBlastWaves.enzotest`:

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
not.  A full listing of options can be found in the :file:`run/README`
file.

Once you have created a new problem type in Enzo and thoroughly
documented the parameters in the Enzo parameter list, you should
follow these steps to add it as a test problem:

1.  Create a new subdirectory in the appropriate place in the
:file:`run/` directory.  If your test problem uses multiple pieces of
physics, put it under the most relevant one.

2.  Add an Enzo parameter file, ending in the extension :file:`enzo`,
for your test problem to that subdirectory.

3.  Add an Enzo test suite parameter file, ending in the extension
:file:`enzotest`.  In that file, add any relevant parameters (as
described in the :file:`run/README file).

4.  Create a "gold standard" set of data for your test problem, by
running with opt-debug and 64-bit precision for floats and
integers. Contact Britton Smith (brittonsmith@gmail.com) and arrange
to send him this data.  Please try to minimize the quantity of data
generated by your calculation by only writing out data at the ending
of the calculation.

If you want to examine the output of your test problem for something
specific, you can optionally add a script that
is indicated by the answer_testing_script parameter.  Look in the
directory :file:`run/Hydro/Hydro-3D/RotatingCylinder` for an example
of how this is done.

Congratulations, you have created a new test problem!

What to do if you fix a bug in Enzo
---------------------------

1. talk about it on the mailing list, get people to sign off

2. create a new set of gold standard datasets 

3. send the new gold standard datasets to britton (who will update the
gold standard)

How to create a new set of reference calculations
---------------------------------------
includes how to add a new gold standard! (64-bit everything, opt-debug)
- britton is our test czar: email him! 

::

    ./test_runner.py --repo=../../enzo-2.0 -o /Users/bwoshea/Desktop/testrun

and

::

    ./test_runner.py --repo=myrepodir  --compare-dir=$HOME/testrun/goldstandarddir -o $HOME/testrun/myoutputdir


and that's it.



Old text
-------

Sample Enzo parameter files are located in the :file:`run/`
subdirectory of the Enzo source distribution. There are a couple
examples of how to run these and some sample output on
:doc:`../tutorials/RunTestProblem`, and each test problem has
:ref:`problem-specific parameters<testproblem_param>`.  Below is the
list of Enzo parameter files in the :file:`run/` subdirectory,
together with an illuminating description of each test:

.. include:: TestProblems.rst

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
estimated runtime of the test, and the dimensionality.  The 
testing suite runs enzo on each selected test problem, then compares
the outputs against those generated from a "good" version of the 
enzo code run on the same problems.  In a few cases, the outputs
are compared directly against analytical solutions.  Lastly, the 
results of these tests are returned to the user for interpretation.

For convenience, three pre-created, overlapping sets of tests are
provided.  For each set of tests, the test suite can automatically
pull the "gold standard" results from a remote server; or one
can generate their own standard locally against which she can compare
different builds of the code.

1.  The "quick suite" (``--suite=quick``).  This is composed of
small calculations that test critical physics packages both
alone and in combination.  The intent of this package is to be run
automatically and relatively frequently (multiple times a day) on 
a remote server to ensure that bugs have not been introduced during the code 
development process.  All runs in the quick suite use no more than 
a single processor.  The total run time should be about 15 minutes.  

2.  The "push suite" (``--suite=push``).  This is a slightly 
large set of tests, encompassing all of the quick suite and 
some additional larger simulations that test a wider variety of physics 
modules.  The intent of this package is to provide a thorough validation 
of the code prior to changes being pushed to the main repository.  The 
total run time is roughly 30 minutes and all simulations use only a single 
processor.  

3.  The "full suite" (``--suite=full``).  This encompasses essentially 
all of test simulations contained within the run directory.  This suite 
provides the most rigorous possible validation of the code in many different 
situations, and is intended to be run prior to major changes being pushed 
to the stable branch of the code.  A small number of simulations in the full 
suite are designed to be run on 2 processors and will take multiple hours to 
complete.  The total run time is roughly 60 hours.  

How to run the test suite
-------------------------

1.  **Compile Enzo.**  The gold standard calculations use the default 
compiler settings that can be restored with ``make default``.  
If you use significantly different compilation options
(high-level optimization in particular) you may see somewhat
different outputs that will result in failed tests.  To compile 
enzo with the standard settings, complete these commands:

::

    cd <enzo_root>/src/enzo
    make default
    make clean
    make

2.  **Get/update yt.**  The enzo tests are generated and compared using the 
yt analysis suite.  If you do not yet have yt, visit 
http://yt-project.org/#getyt for installation instructions.  
If you already have yt and yt is in your path, make sure you're using
the most up-to-date version by running the following command:

::

    yt update

3.  **Generate the standard test files.**  The testing suite operates by 
running a series of enzo test files throughout the ``run/`` subdirectory.
Some unique test files are already generated for specific test problems, 
but the standard generic tests to be run on each test problem need to be 
created by you with the following command: 

::

    cd <enzo_root>/run
    python make_new_tests.py

4.  **Run the test suite.** While remaining in the ``run/`` 
subdirectory, you can initiate the quicksuite test simulations and 
their comparison against the gold standard by running the following 
commands:

::

    python test_runner.py --suite=quick -o <output_dir>

In this comand, ``--suite=quick`` instructs the test runner to
use the quick suite. ``--output-dir=<output_dir>`` instructs the 
test runner to output its results to a user-specified directory 
(preferably outside of the enzo file hierarchy).  For a full 
description of the many flags associated with test_runner.py, see 
below.

5.  **Review the results.**  While the test_runner is executing, you should 
see the results coming up at the terminal in real time, but you can review 
these results in a file output at the end of the run.  The test_runner 
creates a subdirectory in the output directory you provided it, as shown
in the example below.  

::

    ls <output_dir>

    fe7d4e298cb2    


    ls <output_dir>/fe7d4e298cb2    

    Cooling                 GravitySolver           MHD                     test_results.txt 
    Cosmology               Hydro                   RadiationTransport      version.txt

The name of this directory will be the unique hash of the version of 
enzo you chose to run with the testing suite.  Within this directory 
are all of the test problems that you ran along with their simulation
outputs, organized based on test type (e.g. ``Cooling``, ``AMR``, 
``Hydro``, etc.)  Additionally, you should see a file called
``test_results.txt``, which contains a summary of the test runs.  

The testing suite does not expect bitwise agreement between the gold standard
and your results, due to compiler, architecture and operating system
differences between versions of enzo.  There must be a significant 
difference between your result and the gold standard for you to fail 
any tests, thus you should be passing all of the tests.  If you are not, 
then examine more closely what modifications you made to the enzo source
which caused the test failure.  If this is a fresh version of enzo that 
you grabbed and compiled, then you should write the enzo-users email 
list with details of your test run (computer os, architecture, version of enzo, 
version of yt, what test failed, what error message you received), so that
we can address this issue.

For more details about the results of an individual test, examine the
``estd.out`` file in the test problem within this directory hierarchy,
as it contains the stderr and stdout for each test simulation.

How to generate your own set of reference results and compare against them
--------------------------------------------------------------------------

Explanation of all the ``test_runner.py`` flags
-----------------------------------------------

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
``run/Hydro/Hydro-1D/InteractingBlastWaves`` directory:

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

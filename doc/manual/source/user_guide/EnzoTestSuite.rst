.. _EnzoTestSuite:

Enzo Test Suite
===============

The Enzo test suite is a set of tools whose purpose is to perform
regression tests on the Enzo codebase, in order to help developers
discover bugs that they have introduced, to verify that the code is
producing correct results on new computer systems and/or compilers,
and, more generally, to demonstrate that Enzo is behaving as expected
under a wide variety of conditions.

What's in the test suite?
-------------------------

The suite is composed of a large number of individual test problems
that are designed to span the range of physics and dimensionalities
that are accessible using the Enzo code, both separately and in
various permutations.  Tests can be selected based on a variety of
criteria, including (but not limited to) the physics included, the
estimated runtime of the test, and the dimensionality.  The 
testing suite runs enzo on each selected test problem, produces 
a series of outputs, and then uses yt to process these outputs
in a variety of different ways (making projections, looking at
fields, etc.).  The results of these yt analyses are then compared
against similarly generated results from an earlier "good" version 
of the enzo code run on the same problems.  In test problems where
we have them, analytical solutions are compared against the test
results (e.g. shocktubes).  Lastly, a summary of these test results 
are returned to the user for interpretation.

One can run individual tests or groups of tests using the various run time
flags_.  For convenience, three pre-created, overlapping sets of tests are
provided.  For each set of tests, one must generate their own standard locally
against which she can compare different builds of the code.

1.  The "quick suite" (``--suite=quick``).  This is composed of
small calculations that test critical physics packages both
alone and in combination.  The intent of this package is to be run
automatically and relatively frequently (multiple times a day) on 
a remote server to ensure that bugs have not been introduced during the code 
development process.  All runs in the quick suite use no more than 
a single processor.  The total run time should be about 15 minutes 
on the default lowest level of optimization..  

2.  The "push suite" (``--suite=push``).  This is a slightly 
large set of tests, encompassing all of the quick suite and 
some additional larger simulations that test a wider variety of physics 
modules.  The intent of this package is to provide a thorough validation 
of the code prior to changes being pushed to the main repository.  The 
total run time is roughly 60 minutes for default optimization, and 
all simulations use only a single processor.  

3.  The "full suite" (``--suite=full``).  This encompasses essentially 
all of test simulations contained within the run directory.  This suite 
provides the most rigorous possible validation of the code in many different 
situations, and is intended to be run prior to major changes being pushed 
to the stable branch of the code.  A small number of simulations in the full 
suite are designed to be run on 2 processors and will take multiple hours to 
complete.  The total run time is roughly 60 hours for the default lowest
level of optimization.

.. _running:

How to run the test suite
-------------------------


1.  **Compile Enzo.** If you have already built enzo, you can skip this step and
the test will use your existing enzo executable.  To compile enzo with the
standard settings, complete these commands:

::

    $ cd <enzo_root>
    $ ./configure
    $ cd ./src/enzo
    $ make load-config-allphysics
    $ make clean
    $ make

Note that you need not copy the resulting enzo executable to your path,
since the enzo.exe will be symbolically linked from the src/enzo directory
into each test problem directory before tests are run.

This build configuration requires that the Hypre and Grackle libraries are
installed and visible in your compiler's search paths. If you do not have these
libraries available, then you can set:

::

    $ make grackle-no
    $ make hypre-no

.. note::

  If Enzo is compiled without support for the grackle and hypre libraries, tests
  of Enzo modules that depend on these libraries will likely fail.


2.  **Install the necessary Python libraries**  The test suite works
    with both Python 2.x and Python 3.x, but requires python-hglib
    (https://pypi.python.org/pypi/python-hglib) to access Mercurial.
    This should be installable via pip.
    
   
3.  **Get the correct yt version** The enzo tests are generated and compared
using the yt analysis suite.  You must be using yt 3.3.0 or newer in order for
the test suite to work.  If you do not yet have yt, visit
http://yt-project.org/#getyt for installation instructions.  If you already have
yt and yt is in your path, make sure you are using the latest verion of yt by
running the following commands:

::

    $ cd /path/to/yt_mercurial_repository
    $ hg update yt
    $ python setup.py develop

4. **Generate answers to test with.** Run the test suite with these flags within
the ``run/`` subdirectory in the enzo source hierarchy:

::

    $ cd <enzo_root>/run
    $ ./test_runner.py --suite=quick -o <output_dir> --answer-store
                        --answer-name=<test_name> --local 
 
Note that we're creating test answers in this example with the quick suite, but
we could just as well create a reference from any number of test problems using
other test problem flags_.

Here, we are storing the results from our tests locally in a file called
<test_name> which will now reside inside of the ``<output_dir>``.  If you want
to, you can leave off ``--answer-name`` and get a sensible default.

.. _directory layout:

::

    $ ls <output_dir>
    fe7d4e298cb2    <test_name>        

    $ ls <output_dir>/<test_name>
    <test_name>.db

When we inspect this directory, we now see that in addition to the subdirectory
containing the simulation results, we also have a <test_name> subdirectory which
contains python-readable shelve files, in this case a dbm file.  These are the
files which actually contain the reference standard.  You may have a different
set of files or extensions depending on which OS you are using, but don't worry
Python can read this no problem.  Congratulations, you just produced your own
reference standard.  Feel free to test against this reference standard or tar
and gzip it up and send it to another machine for testing.


5.  **Run the test suite using your local answers.** The testing suite operates
by running a series of enzo test files throughout the ``run`` subdirectory.
Note that if you want to test a specific enzo changeset, you must update to it
and recompile enzo. You can initiate the quicksuite test simulations and their
comparison against your locally generated answers by running the following
commands:

::

    $ cd <enzo_root>/run
    $ ./test_runner.py --suite=quick -o <output_dir> --answer-name=<test_name>
                       --local --clobber

In this command, ``--output-dir=<output_dir>`` instructs the test runner to
output its results to a user-specified directory (preferably outside of the enzo
file hierarchy).  Make sure this directory is created before you call
test_runner.py, or it will fail.  The default behavior is to use the quick
suite, but you can specify any set of tests using the ``--suite`` or ``--name``
flags_. We are comparing the simulation results against a local (``--local``)
reference standard which is named ``<test_name>`` also located in the
``<output_dir>`` directory.  Note, we included the ``--clobber`` flag to rerun
any simulations that may have been present in the ``<output_dir>`` under the
existing enzo version's files, since the default behavior is to not rerun
simulations if their output files are already present.  Because we didn't set
the ``--answer-store`` flag, the default behavior is to compare against the
``<test_name>``.


5.  **Review the results.** While the test_runner is executing, you should see
the results coming up at the terminal in real time, but you can review these
results in a file output at the end of the run.  The test_runner creates a
subdirectory in the output directory you provided it, as shown in the example
below.

::

    $ ls <output_dir>
    fe7d4e298cb2    

    $ ls <output_dir>/fe7d4e298cb2    
    Cooling        GravitySolver    MHD                    test_results.txt 
    Cosmology      Hydro            RadiationTransport     version.txt

The name of this directory will be the unique hash of the version of
enzo you chose to run with the testing suite.  In this case it is
``fe7d4298cb2``, but yours will likely be different, but equally
unintelligible.  You can specify an optional additional suffix to be
appended to this directory name using ``--run-suffix=<suffix>``. This
may be useful to distinguish multiple runs of a given version of enzo,
for example with different levels of optimization. Within this
directory are all of the test problems that you ran along with their
simulation outputs, organized based on test type (e.g.  ``Cooling``,
``AMR``, ``Hydro``, etc.)  Additionally, you should see a file called
``test_results.txt``, which contains a summary of the test runs and
which ones failed and why.  

My tests are failing and I don't know why
-----------------------------------------

A variety of things cause tests to fail: differences in compiler,
optimization level, operating system, MPI submission method, 
and of course, your modifications to the code.  Go through your 
``test_results.txt`` file for more information about which tests 
failed and why.  You could try playing with the relative tolerance 
for error using the ``--tolerance`` flag as described in the flags_ 
section.  For more information regarding the failures of a specific 
test, examine the ``estd.out`` file in that test problem's subdirectory
within the ``<output_dir>`` directory structure, as it contains the 
``STDERR`` and ``STDOUT`` for that test simulation.

If you are receiving ``EnzoTestOutputFileNonExistent`` errors, it
means that your simulation is not completing.  This may be due to
the fact that you are trying to run enzo with MPI which your 
system doesn't allow you to initiate from the command line.
(e.g. it expects you to submit mpirun jobs to the queue).  
You can solve this problem by recompiling your enzo executable with
MPI turned off (i.e. ``make use-mpi-no``), and then just pass the 
local_nompi machine flag (i.e. ``-m local_nompi``) to your 
test_runner.py call to run the executable directly without MPI support.  
Currently, only a few tests use multiple cores, so this is not a 
problem in the quick or push suites.

If you see a lot of ``YTNoOldAnswer`` errors, it may mean that your simulation
is running to a different output than what was reached for your locally
generated answers does, and the test suite is trying to compare your last output
file against a non-existent file in the answers.  Look carefully at the
results of your simulation for this test problem using the provided python file
to determine what is happening.  Or it may simply mean that you specified the
wrong answer name.

.. _flags:

Descriptions of all the testing suite flags
-------------------------------------------

You can type ``./test_runner.py --help`` to get a quick summary of all 
of the command line options for the testing suite.  Here is a more 
thorough explanation of each.

**General flags**

``-h, --help``
    list all of the flags and their argument types (e.g. int, str, etc.)

``-o str, --output-dir=str`` default: None
    Where to output the simulation and results file hierarchy.  Recommended
    to specify outside of the enzo source hierarchy.

``-m str, --machine=str`` default: local
    Specify the machine on which you're running your tests.  This loads 
    up a machine-specific method for running your tests.  For instance,
    it might load qsub or mpirun in order to start the enzo executable
    for the individual test simulations.  You can only use machine
    names of machines which have a corresponding machine file in the 
    ``run/run_templates`` subdirectory (e.g. nics-kraken). *N.B.*
    the default, ``local``, will attempt to run the test simulations using
    mpirun, so if you are required to queue on a machine to execute 
    mpirun, ``test_runner.py`` will silently fail before finishing your
    simulation.  You can avoid this behavior by compiling enzo without
    MPI and then setting the machine flag to ``local_nompi``.

``--repo=str`` default: current directory
    Path to repository being tested.

``--interleave`` default: False
    Interleaves preparation, running, and testing of each 
    individual test problem as opposed to default batch
    behavior.

``--clobber`` default: False
    Rerun enzo on test problems which already have 
    results in the destination directory

``--tolerance=int`` default: see ``--strict``
    Sets the tolerance of the relative error in the 
    comparison tests in powers of 10.  

    Ex: Setting ``--tolerance=3`` means that test results
    are compared against the standard and fail if
    they are off by more than 1e-3 in relative error.
    
``--bitwise`` default: see ``--strict``
    Declares whether or not bitwise comparison tests
    are included to assure that the values in output
    fields exactly match those in the reference standard.

``--strict=[high, medium, low]`` default: low
    This flag automatically sets the ``--tolerance``
    and ``--bitwise`` flags to some arbitrary level of
    strictness for the tests.  If one sets ``--bitwise``
    or ``--tolerance`` explicitly, they trump the value
    set by ``--strict``.  When testing enzo general 
    functionality after an installation, ``--strict=low``
    is recommended, whereas ``--strict=high`` is suggested
    when testing modified code against a local reference 
    standard.

    ``high``: tolerance = 13, bitwise = True
    ``medium``: tolerance = 6, bitwise = False
    ``low``: tolerance = 3, bitwise = False

``--sim-only`` default: False
    Only run simulations, do not store the tests or compare them against a 
    standard.

``--test-only`` default: False
    Only perform tests on existing simulation outputs, do not rerun the simulations.

``--time-multiplier=int`` default: 1
    Multiply simulation time limit by this factor.  Useful if you're on a slow
    machine or you cannot finish the specified tests in their allocated time.

``--run-suffix=str`` default: None
    An optional suffix to append to the test run directory. Useful 
    to distinguish multiple runs of a given changeset.

``-v, --verbose`` default: False
    Verbose output in the testing sequence.  Very good for tracking down
    specific test failures.

``--pdb`` default: False
    When a test fails a pdb session is triggered.  Allows interactive inspection
    of failed test data.

``--changeset=str`` default: latest
    Changeset to use in simulation repo.  If supplied,
    make clean && make is also run

    
**Flags for storing, comparing against different standards**

``--answer-store`` default: False
    Should we store the results as a reference or just compare
    against an existing reference?

``--answer-name=str`` default: latest gold standard
    The name of the file where we will store our reference results,
    or if ``--answer-store`` is false, the name of the reference against 
    which we will compare our results. 

``--local`` default: False
    Store/Compare the reference standard locally (i.e. not on the cloud)


**Flags not used**

``--with-answer-testing`` default: False
    DO NOT USE.  This flag is used in the internal yt answer testing
    and has no purpose in the enzo testing infrastructure.

``--answer-big-data`` default: False
    DO NOT USE.  This flag is used in the internal yt answer testing
    and has no purpose in the enzo testing infrastructure.

**Flags for specifying test problems**

These are the various means of specifying which test problems you want
to include in a particular run of the testing suite.

``--suite=[quick, push, full]`` default: None
    A precompiled collection of several different test problems.
    quick: 37 tests in ~15 minutes, push: 48 tests in ~30 minutes, 
    full: 96 tests in ~60 hours.

``--answer_testing_script=str`` default: None

``--AMR=bool`` default: False         
    Test problems which include AMR

``--author=str`` default: None
    Test problems authored by a specific person

``--chemistry=bool`` default: False
    Test problems which include chemistry

``--cooling=bool`` default: False
    Test problems which include cooling

``--cosmology=bool`` default: False   
    Test problems which include cosmology

``--dimensionality=[1, 2, 3]``
    Test problems in a particular dimension

``--gravity=bool`` default: False        
    Test problems which include gravity

``--hydro=bool`` default: False          
    Test problems which include hydro

``--max_time_minutes=float``
    Test problems which finish under a certain time limit

``--mhd=bool`` default: False            
    Test problems which include MHD

``--name=str`` default: None
    A test problem specified by name

``--nprocs=int`` default: 1
    Test problems which use a certain number of processors

``--problematic=bool`` default: False 
    Test problems which are deemed problematic

``--radiation=[None, fld, ray]`` default: None    
    Test problems which include radiation

``--runtime=[short, medium, long]`` default: None
    Test problems which are deemed to have a certain predicted runtime


.. _bisect:

How to track down which changeset caused your test failure
----------------------------------------------------------

In order to identify changesets that caused problems, we have 
provided the ``--bisect`` flag.  This runs hg bisect on revisions 
between those which are marked as --good and --bad.

hg bisect automatically manipulates the repository as it runs its 
course, updating it to various past versions of the code and 
rebuilding.  In order to keep the tests that get run consistent through 
the course of the bisection, we recommend having two separate enzo
installations, so that the specified repository (using ``--repo``) where 
this rebuilding occurs remains distinct from the repository where the 
testing is run.  

To minimize the number of tests run, bisection is only run on tests 
for which ``problematic=True``.  This must be set by hand by the user 
before running bisect.  It is best that this is a single test problem, 
though if multiple tests match that flag, failures are combined with "or"


An example of using this method is as follows:

::

    $ echo "problematic = True" >> Cosmology/Hydro/AdiabaticExpansion/AdiabaticExpansion.enzotest
    $ ./test_runner.py  --output-dir=/scratch/dcollins/TESTS --repo=/SOMEWHERE_ELSE 
                        --answer-compare-name=$mylar/ac7a5dacd12b --bisect --good=ac7a5dacd12b 
                        --bad=30cb5ff3c074 -j 8

To run preliminary tests before bisection, we have also supplied the 
``--changeset`` flag.  If supplied, ``--repo`` is updated to 
``--changeset`` and compiled.  Compile errors cause ``test_runner.py`` 
to return that error, otherwise the tests/bisector is run. 

.. _new_test:

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
    hydro = True
    gravity = False
    AMR = True
    dimensionality = 1
    max_time_minutes = 1
    fullsuite = True
    pushsuite = True
    quicksuite = True

This allows the user to specify the dimensionality, physics used, the
runtime (both in terms of 'short', 'medium', and 'long' calculations,
and also in terms of an actual wall clock time).  A general rule for 
choosing the runtime value is 'short' for runs taking less than 5 minutes, 
'medium' for run taking between 5 and 30 minutes, and 'long' for runs taking 
more than 30 minutes.  If the test problem runs successfully in any amount 
of time, it should be in the full suite, selected by setting 
``fullsuite=True``.  If the test runs in a time that falls under 'medium' 
or 'short', it can be added to the push suite (``pushsuite=True``).  If 
the test is 'short' and critical to testing the functionality of the code, 
add it to the quick suite (``quicksuite=True``).

Once you have created a new problem type in Enzo and thoroughly
documented the parameters in the Enzo parameter list, you should
follow these steps to add it as a test problem:

1.  Create a fork of Enzo.

2.  Create a new subdirectory in the appropriate place in the
``run/`` directory.  If your test problem uses multiple pieces of
physics, put it under the most relevant one.

3.  Add an Enzo parameter file, ending in the extension ``.enzo``,
for your test problem to that subdirectory.

4.  Add an Enzo test suite parameter file, ending in the extension
``.enzotest``.  In that file, add any relevant parameters as described 
above.

5.  By default, the final output of any test problem will be tested by 
comparing the min, max, and mean of a set of fields.  If you want to 
have additional tests performed, create a script in the problem type 
directory and set the ``answer_testing_script`` parameter in the 
``.enzotest`` file to point to your test script.  For an example of 
writing custom tests, see 
``run/Hydro/Hydro-3D/RotatingCylinder/test_rotating_cylinder.py``.

6.  Submit a Pull Request with your changes and indicate that you have 
created a new test to be added to the testing suites.

Congratulations, you've created a new test problem!


What to do if you fix a bug in Enzo
-----------------------------------

It's inevitable that bugs will be found in Enzo, and that some of
those bugs will affect the actual simulation results (and thus the
test problems used in the problem suite).  Here is the procedure for
doing so:

1.  Run the "push suite" of test problems (``--pushsuite=True``)
for your newly-revised version of Enzo, and determine which test
problems now fail.

2.  Visually inspect the failed solutions, to ensure that your new
version is actually producing the correct results!

3.  Email the enzo-developers mailing list at
enzo-dev@googlegroups.com to explain your bug fix, and to show the
results of the now-failing test problems.

4.  Create a pull request for your fix.

.. _http://yt-project.org/#getyt: http://yt-project.org/#getyt

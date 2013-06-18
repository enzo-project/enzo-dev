How to run an Enzo test problem
===============================

Enzo comes with a set of pre-written parameter files which are used
to test Enzo. This is useful when migrating to a new machine with
different compilers, or when new versions of compilers and
libraries are introduced. Also, all the test problems should run to
completion, which is generally not a guarantee!

At the top of each Enzo parameter file is a line like ``ProblemType =
23``, which tells Enzo the type of problem. You can see how this
affects Enzo by inspecting ``InitializeNew.C``. In this
example, this gets called:

.. code-block:: c

      if (ProblemType == 23)
        ret = TestGravityInitialize(fptr, Outfptr, TopGrid, MetaData);

which then calls the routine in ``TestGravityInitialize.C``,
and so on. By inspecting the initializing routine for each kind of
problem, you can see what and how things are being included in the
simulation.

The test problem parameter files are inside the ``run`` subdirectory.
Please see :ref:`EnzoTestSuite` for a full list of test
problems. The files that end in .enzo are the Enzo parameter files,
and .inits are inits parameter files. inits files are only used for
cosmology simulations, and you can see an example of how to run
that in :ref:`RunCosmologySimulation`. Let's try a
couple of the non-cosmology test problems.

ShockPool3D test
----------------

The ShockPool3D is a purely hydrodynamical simulation testing a
shock with non-periodic boundary conditions. Once you've
built enzo (:ref:`obtaining_and_building_enzo`), make a directory
to run the test problem in. Copy enzo.exe and ShockPool3D.enzo into
that directory.
This example test will be run using an interactive session.
On `Kraken <http://www.nics.tennessee.edu/computing-resources/kraken>`_,
to run in an interactive queue, type:

.. highlight:: none

::

    qsub -I -V -q debug -lwalltime=2:00:00,size=12

12 cores (one node) is requested for two hours. Of course, this
procedure may differ on your machine. Once you're in the
interactive session, inside your test run directory, enter:

::

    aprun -n 12 ./enzo.exe -d ShockPool3D.enzo > 01.out

The test problem is run on 12 processors, the debug flag (-d) is
on, and the standard output is piped to a file (01.out). This took
about an hour and twenty minutes to run on Kraken. When it's
finished, you should see ``Successful run, exiting.`` printed to
stderr. Note that if you use other supercomputers, ``aprun`` may be
replaced by 'mpirun', or possibly another command. Consult your
computer's documentation for the exact command needed.

If you want to keep track of the progress of the run, in another
terminal type:

::

    tail -f 01.out
    tail -f 01.out | grep dt

The first command above gives too verbose output to keep track of
the progress. The second one will show what's more interesting,
like the current cycle number and how deep in the AMR hierarchy the
run is going (look for Level[n] where n is the zero-based AMR level
number). This command is especially useful for batch queue jobs
where the standard out always goes to a file.

GravityTest test
----------------

The GravityTest.enzo problem only tests setting up the gravity
field of 5000 particles. A successful run looks like this and
should take less than a second, even on one processor:

::

    test2> aprun -n 1 ./enzo.exe GravityTest.enzo > 01.out
    ****** GetUnits:  1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 *******
    CWD test2
    Global Dir set to test2
    Successfully read in parameter file GravityTest.enzo.
    INITIALIZATION TIME =   6.04104996e-03
    Successful run, exiting.

Other Tests & Notes
-------------------

All the outputs of the tests have been linked to on this page,
below. Some of the tests were run using only one processor, and
others that take more time were run using 16. All tests were run
with the debug flag turned on (which makes the output log, 01.out
more detailed). Enzo was compiled in debug mode without any
optimization turned on (gmake opt-debug). The tests that produce
large data files have only the final data output saved. If you wish
to do analysis on these datasets, you will have to change the
values of GlobalDir, BoundaryConditionName, BaryonFileName and
ParticleFileName in the restart, boundary and hierarchy files to
match where you've saved the data.

PressurelessCollapse
~~~~~~~~~~~~~~~~~~~~

The PressurelessCollapse test required isolated boundary
conditions, so you need to compile Enzo with that turned on (gmake
isolated-bcs-yes). You will also need to turn off the top grid
bookkeeping (gmake unigrid-transpose-no).

Input Files
~~~~~~~~~~~

A few of the test require some input files to be in the run
directory. They are kept in input:

::

    > ls input/
    ATOMIC.DAT  cool_rates.in  lookup_metal0.3.data

You can either copy the files into your run directory as a matter
of habit, or copy them only if they're needed.

Outputs
-------


-  ` AMRCollapseTest.tar.gz <http://lca.ucsd.edu/software/enzo/data/AMRCollapseTest.tar.gz>`_
   - 24 MB
-  ` AMRShockPool2D.tar.gz <http://lca.ucsd.edu/software/enzo/data/AMRShockPool2D.tar.gz>`_
   - 35 KB
-  ` AMRShockTube.tar.gz <http://lca.ucsd.edu/software/enzo/data/AMRShockTube.tar.gz>`_
   - 23 KB
-  ` AMRZeldovichPancake.tar.gz <http://lca.ucsd.edu/software/enzo/data/AMRZeldovichPancake.tar.gz>`_
   - 72 KB
-  ` AdiabaticExpansion.tar.gz <http://lca.ucsd.edu/software/enzo/data/AdiabaticExpansion.tar.gz>`_
   - 31 KB
-  ` CollapseTest.tar.gz <http://lca.ucsd.edu/software/enzo/data/CollapseTest.tar.gz>`_
   - 5.4 MB
-  ` CollideTest.tar.gz <http://lca.ucsd.edu/software/enzo/data/CollideTest.tar.gz>`_
   - 7.6 MB
-  ` DoubleMachReflection.tar.gz <http://lca.ucsd.edu/software/enzo/data/DoubleMachReflection.tar.gz>`_
   - 2.1 MB
-  ` ExtremeAdvectionTest.tar.gz <http://lca.ucsd.edu/software/enzo/data/ExtremeAdvectionTest.tar.gz>`_
   - 430 KB
-  ` GravityStripTest.tar.gz <http://lca.ucsd.edu/software/enzo/data/GravityStripTest.tar.gz>`_
   - 12 MB
-  ` GravityTest.tar.gz <http://lca.ucsd.edu/software/enzo/data/GravityTest.tar.gz>`_
   - 99 KB
-  ` GravityTestSphere.tar.gz <http://lca.ucsd.edu/software/enzo/data/GravityTestSphere.tar.gz>`_
   - 4.6 MB
-  ` Implosion.tar.gz <http://lca.ucsd.edu/software/enzo/data/Implosion.tar.gz>`_
   - 5.6 MB
-  ` ImplosionAMR.tar.gz <http://lca.ucsd.edu/software/enzo/data/ImplosionAMR.tar.gz>`_
   - 3.5 MB



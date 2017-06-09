Particles in Nested Grid Cosmology Simulations
==============================================

When running a nested grid cosmology simulation, not all the
particles created by inits necessarily lie inside of the intended
grid. This has to do with they way particle positions are
calculated from the velocity field. This problem is not a flaw in
the way inits makes initial conditions, but it can lead to
unreliable results if it is not addressed.

**Note:** This effect does not always occur. But it should be
checked for when doing nested initial conditions.

The Problem
-----------

Following the
:doc:`cosmology tutorial </user_guide/CosmologicalInitialConditions>` for
:doc:`nested grids </user_guide/WritingParameterFiles>`,
first inits is run, and then ring is run on the output of inits to
prepare data for the Parallel Root Grid IO mode of Enzo. The contents of the
initial conditions are easily inspected:

.. highlight:: none

::

    $ h5ls Particle*
    ParticleMasses.0         Dataset {1, 2064384}
    ParticleMasses.1         Dataset {1, 262144}
    ParticlePositions.0      Dataset {3, 2064384}
    ParticlePositions.1      Dataset {3, 262144}
    ParticleVelocities.0     Dataset {3, 2064384}
    ParticleVelocities.1     Dataset {3, 262144}

In this example, there are two initial grids. The root grid has
2,064,384 particles, and the nested grid has 262,144. After ring is
run, a number of files with prefixes ``PPos``, ``PVel`` and ``PMass`` are
created. Using eight tasks, here are the contents of the PPos files
for the top grid:

::

    $ h5ls PP*0
    PPos0000.0               Dataset {3, 258304}
    PPos0001.0               Dataset {3, 258304}
    PPos0002.0               Dataset {3, 257792}
    PPos0003.0               Dataset {3, 257792}
    PPos0004.0               Dataset {3, 258304}
    PPos0005.0               Dataset {3, 258304}
    PPos0006.0               Dataset {3, 257792}
    PPos0007.0               Dataset {3, 257792}

And the nested grid:

::

    $ h5ls PP*1
    PPos0000.1               Dataset {3, 32743}
    PPos0001.1               Dataset {3, 32665}
    PPos0002.1               Dataset {3, 32767}
    PPos0003.1               Dataset {3, 32844}
    PPos0004.1               Dataset {3, 32715}
    PPos0005.1               Dataset {3, 32151}
    PPos0006.1               Dataset {3, 32749}
    PPos0007.1               Dataset {3, 32692}

The sum of the particles in the top grid files is 2,064,384
particles, but in the nested grid files it is only 261,326, a
deficit of 818 particles. The missing particles have been thrown
out by ring because they lie outside the nested grid boundaries.

If the sum of the particles in the files after ring has been run is
equal to the original total, the problem is not extant in the
dataset.

The Solution
------------

The solution to this problem is to introduce an extra step between
inits and ring, where particles are moved to the correct grid.
However, when a particle is moved to a grid with a different
refinement, the mass of the particle must be modified. During this
step, when a particle changes grid, this move must be tracked and
its mass updated to reflect the different grid refinement. Please
see :ref:`EnzoPhysicalUnits` 
for more on why the particle mass must be changed when moving
between grids.

One wrinkle to this solution is the ``ParticleMasses`` file *must* be
created by inits, for all grids, along with the ``ParticlePositions``
and ``ParticleVelocities files``. ``CosmologySimulationParticleMassName``
must therefore also be specified as an input in the Enzo parameter
file.

`Linked here <http://barn.enzotools.org/inits_sort/>`_
is a simple `Python <http://python.org/>`_ script
that will fix the initial condition files. After running the
script, run ring on the new initial condition files. The script
requires a Python installation that has both
`Numpy <http://numpy.scipy.org/>`_ and
`h5py <http://code.google.com/p/h5py/>`_. A simple way to gain an
installation of Python with these modules is to install
`yt <http://yt.enzotools.org/>`_, which is one of the
:doc:`data analysis tools </user_guide/AnalyzingWithYT>`
available for Enzo.

Procedure
---------

Save a copy of the script to the same directory as your nested
initial condition files. Edit the top of the file, where noted, to
match your setup. Please note the order items should be entered.
Once the settings are correct, invoke ``python inits_sort.py``. The
updated initial condition files will be placed inside the directory
``new_ICs``. Then run ring on the new initial condition files, and use
the results with Enzo.



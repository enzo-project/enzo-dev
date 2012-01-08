.. _sample-parameter-files:

Sample inits and Enzo parameter files
=====================================

This page contains a large number of example inits and Enzo parameter
files that should cover any possible kind of Enzo cosmology simulation
that you are interested in doing. All should run with minimal
tinkering. They can be downloaded separately below, or as a single
tarball.

Note: unless otherwise specified, inits is run by calling

.. highlight:: none

::

      inits -d <name of inits parameter file>

and Enzo is run by calling

::

      [mpirun ...] enzo -d <name of enzo parameter file>

In both cases, the -d flag displays debugging information, and can be
omitted. Leaving out the -d flag can significantly speed up Enzo
calculations. Also note that Enzo is an MPI-parallel program, whereas
inits is not.

**Unigrid dark matter-only cosmology simulation**.  This is the
simplest possible Enzo cosmology simulation - a dark matter-only
calculation (so no baryons at all) and no adaptive mesh. See the
:download:`inits parameter file <./samples/dmonly.inits>` and
:download:`Enzo parameter file <./samples/dmonly_unigrid.enzo>`.

**AMR dark matter-only cosmology simulation**.  This is a dark
matter-only cosmology calculation (using the same initial conditions
as the previous dm-only run) but with adaptive mesh refinement turned
on.  See the :download:`inits parameter file <./samples/dmonly.inits>`
and :download:`Enzo parameter file <./samples/dmonly_amr.enzo>`.

**Unigrid hydro+dark matter cosmology simulation (adiabatic)**.  This
is a dark matter plus hydro cosmology calculation **without** adaptive
mesh refinement and no additional physics.  See the :download:`inits
parameter file <./samples/gas_plus_dm.inits>` and :download:`Enzo
parameter file <./samples/gas_plus_dm_uni_adia.enzo>`.

**AMR hydro+dark matter cosmology simulation (adiabatic)**.  This is a
dark matter plus hydro cosmology calculation (using the same initial
conditions as the previous dm+hydro run)**with** adaptive mesh
refinement (refining everywhere in the simulation volume) and no
additional physics.  See the :download:`inits parameter file
<./samples/gas_plus_dm.inits>` and :download:`Enzo parameter file
<./samples/gas_plus_dm_amr_adia.enzo>`.

**AMR hydro+dark matter cosmology simulation (lots of physics)**.
This is a dark matter plus hydro cosmology calculation (using the same
initial conditions as the previous two dm+hydro runs) **with**
adaptive mesh refinement (refining everywhere in the simulation
volume) and including radiative cooling, six species primordial
chemistry, a uniform metagalactic radiation background, and
prescriptions for star formation and feedback.  See the
:download:`inits parameter file <./samples/gas_plus_dm.inits>` and
:download:`Enzo parameter file <./samples/gas_plus_dm_amr_multiphys.enzo>`.

**AMR hydro+dark matter nested-grid cosmology simulation (lots of
physics)**.  This is a dark matter plus hydro cosmology calculation
with two static nested grids providing excellent spatial and dark
matter mass resolution for a single Local Group-sized halo and its
progenitors. This simulation only refines in a small subvolume of the
calculation, and includes radiative cooling, six species primordial
chemistry, a uniform metagalactic radiation background, and
prescriptions for star formation and feedback. All parameter files can
be downloaded in :download:`one single tarball
<./samples/example_ics.tar.gz>`. Note that inits works differently
for multi-grid setups. Instead of calling inits one time, it is called
N times, where N is the number of grids. For this example, where there
are three grids total (one root grid and two nested subgrids), the
procedure would be::

  inits -d -s SubGridFile.inits  TopGridFile.inits
  inits -d -s SubSubGridFile.inits SubGridFile.inits 
  inits -d SubSubGridFile.inits 

(but note that there is now an easier way to do multiple-grid
initialization with inits -- see :ref:`using_inits`).
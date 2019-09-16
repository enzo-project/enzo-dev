.. _test_problems:

Enzo Test Problems
===================

The following is a list of the test problems that are available in the 
``run`` subdirectory found in the Enzo repository. These problems are also the
same set of problems used in the test suite. They are listed in broad
categories following the directory structure.

This list however is not a complete list of all the ProblemTypes that 
are available in Enzo.

The test problem specific  parameters can be found in the :ref:`parameters`.

.. toctree::
   :maxdepth: 2


* `Cooling`_ 

   * `CoolingTest_Cloudy`_

   * `CoolingTest_Grackle`_ 

   * `CoolingTest_JHW`_  

   * `OneZoneFreefallTest`_


* `Cosmology`_

   * `AdiabaticExpansion`_

   * `AMRZeldovichPancake`_  
   
   * `AMRZeldovichPancake_Streaming`_

   * `MHDAdiabaticExpansion_CT`_       
   
   * `MHDAdiabaticExpansion_Dedner`_  
   
   * `MHDZeldovichPancake`_           

   * `MHDZeldovichPancake_2_CT`_      

   * `MHDZeldovichPancake_2_Dedner`_

   * `SphericalInfall`_

   * `ZeldovichPancake`_


* `Cosmology Simulation`_

   * `amr_cosmology`_

   * `amr_nested_cosmology`_

   * `dm_only`_

   * `ReionizationHydro`_

   * `ReionizationRadHydro`_


* `DrivenTurbulence3D`_

* `FLD`_

* `FuzzyDarkMatter`_

* `Gravity Solver`_

   * `BinaryCollapse`_  

   * `BinaryCollapseMHDCT`_  

   * `GravityStripTest`_  

   * `GravityTest`_  

   * `GravityTestSphere`_  

   * `MaximumGravityRefinementTest`_  

   * `TestOrbit`_  

   * `TestOrbitMRP`_

* `Hydro/Hydro-1D`_

   * `FreeExpansion`_

   * `InteractingBlastWaves`_

   * `PressurelessCollapse`_

   * `ShockInABox`_

   * `SodShockTube`_

   * `Toro-1-ShockTube`_

   * `Toro-1-ShockTubeAMR`_

   * `Toro-2-ShockTube`_

   * `Toro-2-ShockTubeAMR`_

   * `Toro-3-ShockTube`_

   * `Toro-3-ShockTubeAMR`_

   * `Toro-4-ShockTube`_

   * `Toro-4-ShockTubeAMR`_

   * `Toro-5-ShockTube`_

   * `Toro-5-ShockTubeAMR`_

   * `Toro-6-ShockTube`_

   * `Toro-7-ShockTube`_

   * `WavePool`_

* `Hydro/Hydro-2D`_

   * `AMRShockPool2D`_

   * `Athena-RayleighTaylor`_

   * `DoubleMachReflection`_

   * `FreeExpansionAMR`_

   * `HDMHD2DCheckOddEvenCouplingOfRiemannSolver`_

   * `Implosion`_

   * `ImplosionAMR`_

   * `KelvinHelmholtz`_

   * `KelvinHelmholtzAMR`_

   * `NohProblem2D`_

   * `NohProblem2DAMR`_

   * `RadiatingShockWave`_

   * `RampedKelvinHelmholtz2D`_

   * `SedovBlast`_

   * `SedovBlastAMR`_

   * `ShockPool2D`_

   * `ValidatedNonlinearKelvinHelmholtz`_

* `Hydro/Hydro-3D`_

   * `AgoraGalaxy`_

   * `Athena-RayleighTaylor3D`_

   * `CollapseTestNonCosmological`_

   * `CollideTest`_

   * `ExtremeAdvectionTest`_

   * `GalaxySimulation`_

   * `NFWCoolCoreCluster`_

   * `NohProblem3D`_

   * `NohProblem3DAMR`_

   * `ProtostellarCollapse_Std`_

   * `RotatingCylinder`_

   * `RotatingSphere`_

   * `ShockPool3D`_

   * `StripTest`_

* `MHD/1D`_

   * `BrioWu-MHD-1D`_

   * `BrioWu-MHD-1D-MHDCT`_

   * `CR-ShockTube`_

   * `MHD_Metal_Advection_CT`_

   * `MHD_Metal_Advection_Dedner`_

* `MHD/2D`_

   * `LoopAdvection_CT`_

   * `LoopAdvection_Dedner`_

   * `MHD2DRotorTest`_

   * `MHDCTOrszagTang`_

   * `MHDCTOrszagTangAMR`_

   * `MHDDednerOrszagTang`_

   * `RayleighTaylor_CT_Suppressed`_

   * `SedovBlast-MHD-2D-Fryxell`_

   * `SedovBlast-MHD-2D-Gardiner`_

   * `Wengen2-CollidingFlow`_

* `MHD/3D`_

   * `ShearingBox`_

   * `StochasticForcing`_

* `RadiationTransport`_

   * `PhotonShadowing`_

   * `PhotonTest`_
   
   
   * `PhotonTestAMR`_  

   * `PhotonTestMultiFrequency`_

* `RadiationTransportFLD`_

   * `CosmoIonization_q05z10`_

   * `CosmoIonization_q05z10_sp`_
   
   * `CosmoIonization_q05z4`_

   * `CosmoIonization_q05z4_sp`_

   * `CosmoIonization_q5z10`_

   * `CosmoIonization_q5z10_sp`_

   * `CosmoIonization_q5z4`_

   * `Grey_Enzochem`_

   * `Grey_Split`_

   * `RadiatingShockLab`_

   * `RadiatingShockLab1D`_

   * `RadiatingShockLab1D_sp`_

   * `RadiationStream1D`_

   * `RadiationStream1D_sp`_

   * `RadiationStreamX0`_

   * `RadiationStreamX1`_

   * `RadiationStreamX1_sp`_

   * `RadiationStreamY0`_

   * `RadiationStreamY0_sp`_

   * `RadiationStreamY1`_

   * `RadiationStreamY1_sp`_
 
   * `RadiationStreamZ0`_
   
   * `RadiationStreamZ0_sp`_

   * `RadiationStreamZ1`_

   * `RHIonization1`_

   * `RHIonization1_sp`_

   * `RHIonization2`_

   * `RHIonization2_sp`_

   * `TurnerStoneEquil1`_

   * `TurnerStoneEquil2`_

* `StarParticle`_


.. _CoolingProblems:

Cooling
~~~~~~~

.. _CoolingTest_Cloudy:

CoolingTest_Cloudy
^^^^^^^^^^^^^^^^^^
This test problem will set up a single grid that varies smoothly in density, 
metallicity, and temperature, then iterate the rate equations in the chemisty 
module for 50,000 years with hydro deactivated.  The code will make an 
output at the end of the run that includes the cooling time.  This problem 
type did not exist in Enzo 1.5, so there is no comparison.

The cooling tests will run in a few minutes on a single processor.

The three parameter files are:
CoolingTest_Cloudy.enzo - uses Cloudy cooling along with the MultiSpecies = 1 chemistry.  
The input data provided is a three dimensional table that varies in density, metallicity, 
and temperature.

Cooling data files:
primordial_cie.dat - CIE cooling rates for atomic H and He taken from Black (1981).
solar_2008_3D_metals.h5 - input data for Cloudy cooling.

The script plot.py will plot cooling rates from the cooling test 
along with the H/He cooling rate from Black (1981) and the Z = Zsun 
rate from Sarazin & White (1987)


.. _CoolingTest_Grackle:

CoolingTest_Grackle
^^^^^^^^^^^^^^^^^^^
This test problem will set up a single grid that varies smoothly in density, 
metallicity, and temperature, then iterate the rate equations in the chemisty 
module for 50,000 years with hydro deactivated.  The code will make an 
output at the end of the run that includes the cooling time.

The cooling tests will run in a few minutes on a single processor.

The three parameter files are:
CoolingTest_Grackle.enzo - uses Grackle cooling along with the 
non-equilibrium atomic H/He chemistry.

Cooling data files:
primordial_cie.dat - CIE cooling rates for atomic H and He taken from 
Black (1981).
CloudyData_UVB=HM2012.h5 - input data for Grackle cooling.

The script plot.py will plot cooling rates from the cooling test 
along with the H/He cooling rate from Black (1981) and the Z = Zsun 
rate from Sarazin & White (1987)


.. _CoolingTest_JHW:

CoolingTest_JHW
^^^^^^^^^^^^^^^
This test problem will set up a single grid that varies smoothly in density, 
metallicity, and temperature, then iterate the rate equations in the chemisty 
module for 50,000 years with hydro deactivated.  The code will make an 
output at the end of the run that includes the cooling time.  This problem 
type did not exist in Enzo 1.5, so there is no comparison.

The cooling tests will run in a few minutes on a single processor.

The three parameter files are:
CoolingTest_JHW.enzo - uses John Wise's metal cooling along with 
the MultiSpecies = 1 chemitry.

Cooling data files:
primordial_cie.dat - CIE cooling rates for atomic H and He taken from 
Black (1981).
cool_rates.in - analytic cooling rates for Z = 0.5 and 1 Zsun from 
Sarazin & White (1987).
metal_cool.dat - input data for John Wise's metal cooling.

The script plot.py will plot cooling rates from the cooling test
along with the H/He cooling rate from Black (1981) and the Z = Zsun 
rate from Sarazin & White (1987)
.. _OneZoneFreefallTest:

OneZoneFreefallTest
^^^^^^^^^^^^^^^^^^^
This test problem will set up a 2D grid varying in energy and metallicity.  
All points have the same density, which evolves according to the analytical 
solution for free-fall collapse.  The timestep is calculated as a fraction of 
the free-fall time.  Since the timestep continually decreases, outputs are 
done on cycles.  This test problem can be used to test chemistry and 
cooling routines.

The script, plot.py, will create plots of n vs. T (and Tdust), n
vs. f_H2, and n vs. t_cool/t_dyn.  If using H2 formation on 
dust grains, set dust=True on line 10.  Run this script like this:

python plot.py OneZoneFreefallTest.enzo

.. _CosmologyProblems:

Cosmology
~~~~~~~~~

.. _AdiabaticExpansion:

AdiabaticExpansion
^^^^^^^^^^^^^^^^^^

A test for time-integration accuracy of the expansion terms (Bryan thesis 1996, Sect. 3.3.3).

.. _AMRZeldovichPancake:  

AMRZeldovichPancake  
^^^^^^^^^^^^^^^^^^^

This test simulates a collapsing sinusoidal cosmological pertubation
in one-dimension, commonly known as a Zel'dovich pancake.  This
problem tests both the hydrodynamics and gravity solvers and the
implementation of cosmological expansion.  The system will form a
caustic in the center of the domain with a density and temperature
peak.  There should be a small dip in temperature at the center of the
broad peak.  In flat cosmology, there exists an analytical solution in
the linear phase of collapse and is given in Zel'dovich (1970).

This test runs to completion and creates 2 outputs -- the initial
(z=20) and final (z=0) states.  There are two levels of refinement by
factors of 4.  The finest resolution element is the same as the
non-AMR version of this test.  There are some small differences between
Enzo v1.5 and v2.0 at the parent-child grid boundaries.

.. _AMRZeldovichPancake_Streaming:

AMRZeldovichPancake_Streaming
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(John Wise, July 2010)

This test simulates a collapsing sinusoidal cosmological pertubation
in one-dimension, commonly known as a Zel'dovich pancake.  This
problem tests both the hydrodynamics and gravity solvers and the
implementation of cosmological expansion.  The system will form a
caustic in the center of the domain with a density and temperature
peak.  There should be a small dip in temperature at the center of the
broad peak.  In flat cosmology, there exists an analytical solution in
the linear phase of collapse and is given in Zel'dovich (1970).

This test runs to completion and creates 2 outputs -- the initial
(z=20) and final (z=0) states.  There are two levels of refinement by
factors of 4.  The finest resolution element is the same as the
non-AMR version of this test.  There are some small differences between
Enzo v1.5 and v2.0 at the parent-child grid boundaries.

(John Wise, February 2019)

Modified to add DM particles and a bulk gas velocity

.. _MHDAdiabaticExpansion_CT:

MHDAdiabaticExpansion_CT
^^^^^^^^^^^^^^^^^^^^^^^^

Adiabatic expansion test for MHD, using Athena CT.

This test differs from the PPM Adiabatic Expansion by the increased initial
temperature.  The extremely low floor on the PPM version makes it a poor test,
since the thermal expansion is dominated by the temperature floor, rather than
physical integration. 

This test is not entirely uniform for two reasons.  
First is self gravity (actually off in this version) which causes issues at the
corners of the domain as well as subgrid boundary.
The second is the time-interpolation in the boundary for the subgrid, which
causes slight acceleration due to slightly different expansion of the fluid in
the boundary of the subgrid relative to the "ideal" solution.

.. _MHDAdiabaticExpansion_Dedner:

MHDAdiabaticExpansion_Dedner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Adiabatic expansion test for MHD, using Dedner.

This test differs from the PPM Adiabatic Expansion by the increased initial
temperature.  The extremely low floor on the PPM version makes it a poor test,
since the thermal expansion is dominated by the temperature floor, rather than
physical integration. 

This test is not entirely uniform for two reasons.  
First is self gravity (actually off in this version) which causes issues at the
corners of the domain as well as subgrid boundary.
The second is the time-interpolation in the boundary for the subgrid, which
causes slight acceleration due to slightly different expansion of the fluid in
the boundary of the subgrid relative to the "ideal" solution.

.. _MHDZeldovichPancake:

MHDZeldovichPancake
^^^^^^^^^^^^^^^^^^^


.. _MHDZeldovichPancake_2_CT:   

MHDZeldovichPancake_2_CT
^^^^^^^^^^^^^^^^^^^^^^^^
This is another iteration of Zel'Dovich pancake.  This is tuned to almost 
reproduce the result from Collins et al 2010, as well as the Dedner run in 
run/Cosmology/MHDZeldovichPancake_2_Dedner

Slight differences with the method paper exist due to the uniform initial
distribution in the method paper.  

Neither this nor the Dedner version use Dual Energy Formalism, in order to match
the temperature field as well as possible.
.. _MHDZeldovichPancake_2_Dedner:

MHDZeldovichPancake_2_Dedner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is another iteration of Zel'Dovich pancake.  This is tuned to almost 
reproduce the result from Collins et al 2010, as well as the Dedner run in 
run/Cosmology/MHDZeldovichPancake_2_Dedner

Slight differences with the method paper exist due to the uniform initial
distribution in the method paper.  

Neither this nor the Dedner version use Dual Energy Formalism, in order to match
the temperature field as well as possible.
.. _SphericalInfall:

SphericalInfall
^^^^^^^^^^^^^^^


.. _ZeldovichPancake:

ZeldovichPancake
^^^^^^^^^^^^^^^^
(John Wise, July 2010)

This test simulates a collapsing sinusoidal cosmological pertubation
in one-dimension, commonly known as a Zel'dovich pancake.  This
problem tests both the hydrodynamics and gravity solvers and the
implementation of cosmological expansion.  The system will form a
caustic in the center of the domain with a density and temperature
peak.  There should be a small dip in temperature at the center of the
broad peak.  In flat cosmology, there exists an analytical solution in
the linear phase of collapse and is given in Zel'dovich (1970).

This test runs to completion and creates 2 outputs -- the initial
(z=20) and final (z=0) states.  There is no refinement.  Enzo v1.5 and
v2.0 produce exactly the same results.
.. _Cosmology Simulation:

Cosmology Simulation
~~~~~~~~~~~~~~~~~~~~

.. _amr_cosmology:

amr_cosmology
^^^^^^^^^^^^^

This is a cosmology simulation with Cen & Ostriker star
formation/feedback, cooling with Grackle, and initial conditions made
with MUSIC. It consists of a 32 Mpc/h box with 32^3 root grid cells
and dark matter particles and 5 levels of AMR.

This simulation will run to z = 0 in about 10 minutes on a single
core. To run this, you will need to copy the file,
CloudyData_UVB=HM2012.h5, from the input directory of your Grackle
source, to the run directory of the simulation.

To generate the ICs with MUSIC:
./MUSIC amr_cosmology.conf

Initial conditions can also be downloaded from the "Enzo test data"
collection on the yt Hub (hub.yt), or do:
pip install girder-client
girder-cli --api-url https://girder.hub.yt/api/v1 download 5afb0040ec1bd30001fcd002

To run the simulation:
./enzo.exe -d amr_cosmology.enzo
.. _amr_nested_cosmology:

amr_nested_cosmology
^^^^^^^^^^^^^^^^^^^^
This is a variation of the amr_cosmology simulation with 1 level of
initial refinement and must-refine-particles. The simulation uses Cen
& Ostriker star formation/feedback, cooling with Grackle, and initial
conditions made with MUSIC. The nested refinement zooms in on the most
massive halo in the box at z = 1.8. Must-refine-particles are used to
allow AMR only for particles that end up within 3 virial radii of the
most massive halo at z = 1.8. The simulation will run in about 1
minute on a single core.

The zoom initial conditions with must-refine-particle flagging was
created using the method outlined in
https://bitbucket.org/jwise77/enzo-mrp-music

Initial conditions can be downloaded from the "Enzo test data"
collection on the yt Hub (hub.yt), or do:
pip install girder-client
girder-cli --api-url https://girder.hub.yt/api/v1 download 5afef79bec1bd30001fcd07e

To run the simulation:
./enzo.exe -d amr_nested_cosmology.enzo

.. _dm_only:

dm_only
^^^^^^^
This is a dark-matter-only version of the amr_cosmology simulation. It
consists of a 32 Mpc/h box with 32^3 dark matter particles and 5
levels of AMR. This simulation will run to z = 0 in less than a minute
on a single core.

To generate the ICs with MUSIC:
./MUSIC dm_only.conf

Initial conditions can also be downloaded from the "Enzo test data"
collection on the yt Hub (hub.yt), or do:
pip install girder-client
girder-cli --api-url https://girder.hub.yt/api/v1 download 5afb0145ec1bd30001fcd024

To run the simulation:
./enzo.exe -d dm_only.enzo

.. _ReionizationHydro:

ReionizationHydro
^^^^^^^^^^^^^^^^^
This is a cosmology simulation that simulates reionization using the 
convention, non-radiative star formation and feedback and a Haardt &
Madau background.  It will run on 2 processors in about 20 minutes.

Usage:
./inits.exe -d ReionizationHydro.inits
mpirun -np 2 ./ring.exe pv ParticlePositions ParticleVelocities
mpirun -np 2 ./enzo.exe -d ReionizationHydro.enzo

.. _ReionizationRadHydro:

ReionizationRadHydro
^^^^^^^^^^^^^^^^^^^^
This is a cosmology simulation that simulates reionization using the 
ray tracing radiation transfer method with radiating star particles
and a Haardt & Madau background.  It will run on 2 processors in about
40 minutes.

Usage:
./inits.exe -d ReionizationRadHydro.inits
mpirun -np 2 ./ring.exe pv ParticlePositions ParticleVelocities
mpirun -np 2 ./enzo.exe -d ReionizationRadHydro.enzo

.. _DrivenTurbulence3D:

DrivenTurbulence3D
~~~~~~~~~~~~~~~~~~

Unless hydromethod == 4  this will do hydrodynamic turbulence. 
 Tom Abel 2009

This can do fixed force pattern driving as well as decaying turbulence set ups.
Set UseDrivingField = 1 to use the driving with HydroMethod 3 or 4 (hydro/MHD)
Only decaying is implemeted for HydroMethod < 3 (Zeus & standard PPM)

.. _FLD:

FLD
~~~

.. _FuzzyDarkMatter:

FuzzyDarkMatter
~~~~~~~~~~~~~~~

.. _Gravity Solver:

Gravity Solver
~~~~~~~~~~~~~~

.. _BinaryCollapse:

BinaryCollapse
^^^^^^^^^^^^^^
(Stephen Skory, June 2010)


This test runs to completion. There are 8 data dumps, DD0000 to DD0007.
This test does not work in 1.5 because problem 202 is not defined there.

I ran my tests with opt-debug using the intel/openmpi stack on Triton on 3 June 2010.

2.0 test:

Running with the 2.0 enzo produces a binary collapse of
matter that begins to form two dense clumps best viewed along the z-axis.
The yt plotting script included will plot this nicely, zoomed in on the inner
third of the volume.

I have made some changes to this file from how it was previously.
I removed the magnetic field, and I changed the stop time so this test takes
roughly 1 hour on 8 cores. As it was written before, it would have taken many
hours to reach completion, and I am fundamentally against test problems take more than 1 hour.
.. _BinaryCollapseMHDCT:

BinaryCollapseMHDCT
^^^^^^^^^^^^^^^^^^^

(Stephen Skory, June 2010, David Collins, January 2015)

This is an update to the test in run/GravitySolver/BinaryCollapse, to include
MHDCT.  A shorter test to ensure MHDCT runs.  This runs in less than 5 minutes.

.. _GravityStripTest:

GravityStripTest
^^^^^^^^^^^^^^^^

.. _GravityTest:

GravityTest
^^^^^^^^^^^
(Greg Bryan, July 2010)

This test places a single, massive particle at the center of the box
and then places 5000 nearly massless particles throughout the rest
of the box, randomly spaced in log radius from the center, and
randomly placed in angle.  A single small step is taken, and the
velocity is then divided by the timestep to measure the acceleration.
A single subgrid is placed in the center from 0.4375 to 0.5625
(in units where the box size is 1.0). 

This tests the acceleration of the particles from a single point
mass and so can be directed compared to the r^-2 expected
result.  An output file is generated, called TestGravityCheckResults.out,
which contains four columns with one entry for each of the
5000 particles.  The columns are the radius (in units of the
cell length of the most refined grid), the tangential component
of the measured force, as computed by the code, 
the radial component of the computed force, and finally
the "true" (unsoftened) force.  

The tangential component of the force should be zero; the
radial component should follow the r^-2 law, but is softened
for radii less than about one cell length (or slightly larger).
It also falls below r^-2 are large distances because this
problem uses periodic boundary conditions (The code has
been modified since this problem was originally written
to use isolated boundary conditions, but this problem
has not been changed).

The test for this problem is to compute the rms force
error between 1 and 8 cell distances from the center.
For CIC, this is measured to be about 5%.  The test
checks to make sure this has not changed by more than 1%.
This is not an ideal check, as many errors could
conceivably escape undetected (e.g. those having to do
with force errors at small and large radii); however, the problem
with a bitwise comparison is that the positions of
the 5000 particles are random (with no setable seed).

.. _GravityTestSphere:

GravityTestSphere
^^^^^^^^^^^^^^^^^
(Greg Bryan, July 2010)

This test places a sphere with radius 0.01 and a very
large overdensity in an AMR box with size 1.0 and
computes the resulting velocity change due to gravity.
It is similar to TestGravity but uses a sphere and
baryons instead of a single particle.  No analytic
solution is computed.  More work on this test should
be done to compare to an analytic solution.


.. _MaximumGravityRefinementTest:

MaximumGravityRefinementTest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This test must be run on 2 cores, otherwise the problem isn't triggered.


This test is to ensure that MaximumGravityRefinement is functional.  As of this
writing (2014-07-12) when MaximumGravityRefinement < MaximumRefinementLevel, the
wrong SiblingList is used in PrepareDensityField.  

The current test creates 2 grids on Level 1, and one grid on Level 2, in which
case the code tries to iterate over the SiblingList for level=1 when on level=2.
This causes a seg fault.

.. _TestOrbit:

TestOrbit
^^^^^^^^^
(Greg Bryan, July 2010)

This test places two particles in an otherwise empty box, one at the
center (0.5 0.5 0.5) with mass 1.0, and a second at position 0.2 0.5
0.5 with mass 1.0e-6.  The initial velocity is set such that the
second particle is in a circular orbit around the first.  The system
is evolved for 2.0 time units without AMR but with isolated gravity
boundary conditions.  This makes it complete a bit over 1 orbit.
The test is to check the final particle positions.

For the test to run with enzo 1.5 compile with:
  make isolated-bcs-yes
  make unigrid-transpose-no

.. _TestOrbitMRP:

TestOrbitMRP
^^^^^^^^^^^^
Orbit Test Problem with MRPs (MustRefineParticles)

NOTE: This version of the parameter file writes
out the gravitational potential, and is intended
to be used with run/TestOrbit/make_TE_plot.py 
to make a graph of total particle energy as a 
function of orbit!


.. _Hydro/Hydro-1D:

Hydro/Hydro-1D
~~~~~~~~~~~~~~

.. _FreeExpansion:

FreeExpansion
^^^^^^^^^^^^^
(John Wise, July 2010)

This test simulates a blastwave in the free expansion phase.  In the
initial setup, the interior region has a uniform density and a
linearly increasing radial velocity.  The blastwave should advect
outwards, and create a high entropy shell and have little oscillations
in the shock.  

This test runs to completion and creates 41 outputs.  This test
problem is new for version 2.0.  It uses the new hydro_rk solver.

The initial setup is taken from Truelove & McKee, 1999, ApJS, 120,
299.

.. _InteractingBlastWaves:

InteractingBlastWaves
^^^^^^^^^^^^^^^^^^^^^
Two interacting blast waves

This is the first test problem in Woodward & Colella (1984), 
JCP, 54, 115.  With the outer tenths of the domain 
overpressurized, two blast waves move in toward the center.  The
boundaries are reflecting.  One can see the solution by ATHENA in 
Stone et al. (2008), ApJS, 178, 137.

.. _PressurelessCollapse:

PressurelessCollapse
^^^^^^^^^^^^^^^^^^^^

.. _ShockInABox:

ShockInABox
^^^^^^^^^^^


.. _SodShockTube:

SodShockTube
^^^^^^^^^^^^
This is a fairly mild test. The solution consists of left
rarefaction wave, a contact discontinuity, and a right shock.

.. _Toro-1-ShockTube:

Toro-1-ShockTube
^^^^^^^^^^^^^^^^
This is Problem #1 from Chapter 10.8 in Toro's "Riemann Solvers and
Numerical Methods for Fluid Dynamics" (2nd edition).

The solution to this test consists of a left sonic rarefaction
wave, a right travelling contact discontinuity, and a right
shock. It is useful for assessing the entropy satisfaction property
of numerical methods.

.. _Toro-1-ShockTubeAMR:

Toro-1-ShockTubeAMR
^^^^^^^^^^^^^^^^^^^
AMR Version of Toro Problem #1

.. _Toro-2-ShockTube:

Toro-2-ShockTube
^^^^^^^^^^^^^^^^

This is Problem #2 from Chapter 10.8 in Toro's "Riemann Solvers and
Numerical Methods for Fluid Dynamics" (2nd edition).

The solution to this test consists of two symmetric strong
rarefaction waves and a trivial contact discontinuity. The region
between the two non-linear waves is close to vacuum, thus
testing the numerical performance for low density flows.


.. _Toro-2-ShockTubeAMR:

Toro-2-ShockTubeAMR
^^^^^^^^^^^^^^^^^^^
AMR Version of Toro Problem #2

.. _Toro-3-ShockTube:

Toro-3-ShockTube
^^^^^^^^^^^^^^^^
This is Problem #3 from Chapter 10.8 in Toro's "Riemann Solvers and
Numerical Methods for Fluid Dynamics" (2nd edition).

This problem is the left half of Woodward & Colella's blast wave
problem. The solution consists of a left rarefaction wave, a
contact discontinuity, and a strong right shock wave (shock Mach
number 198).

.. _Toro-3-ShockTubeAMR:

Toro-3-ShockTubeAMR
^^^^^^^^^^^^^^^^^^^

AMR Version of Toro Problem #3

.. _Toro-4-ShockTube:

Toro-4-ShockTube
^^^^^^^^^^^^^^^^
This is Problem #4 from Chapter 10.8 in Toro's "Riemann Solvers and
Numerical Methods for Fluid Dynamics" (2nd edition).

This very severe test is made up of the right and left shocks
emerging from the solution to the left and right half of Woodward &
Colella's blast wave test problem. The collision of these two
strong shocks results in three right travelling discontinuities: a
slow left shock, a contact discontinuity, and a right shock.

.. _Toro-4-ShockTubeAMR:

Toro-4-ShockTubeAMR
^^^^^^^^^^^^^^^^^^^

AMR Version of Toro Problem #4

.. _Toro-5-ShockTube:

Toro-5-ShockTube
^^^^^^^^^^^^^^^^
This is Problem #5 from Chapter 10.8 in Toro's "Riemann Solvers and
Numerical Methods for Fluid Dynamics" (2nd edition).

This test is designed to assess the code's ability to resolve
slowly-moving contact discontinuities. Its solution consists of a
left rarefaction wave, a right travelling shock, and a stationary
contact discontinuity.

.. _Toro-5-ShockTubeAMR:

Toro-5-ShockTubeAMR
^^^^^^^^^^^^^^^^^^^

AMR Version of Toro Problem #5

.. _Toro-6-ShockTube:

Toro-6-ShockTube
^^^^^^^^^^^^^^^^
This is Problem #6 from Chapter 10.8 in Toro's "Riemann Solvers and
Numerical Methods for Fluid Dynamics" (2nd edition).

This test consists of an isolated stationary contact wave. It
demonstrates the advantage of the HLLC Riemann solver over the HLL
solver in capturing stationary and slowly moving contact
waves (see also Toro-7-ShockTube).

.. _Toro-7-ShockTube:

Toro-7-ShockTube
^^^^^^^^^^^^^^^^
This is Problem #7 from Chapter 10.8 in Toro's "Riemann Solvers and
Numerical Methods for Fluid Dynamics" (2nd edition).

This test consists of an isolated contact wave moving slowly to the
right. It demonstrates the advantage of the HLLC Riemann solver
over the HLL solver in capturing stationary and slowly moving
contact waves (see also Toro-6-ShockTube).

.. _WavePool:

WavePool
^^^^^^^^


.. _Hydro/Hydro-2D:

Hydro/Hydro-2D
~~~~~~~~~~~~~~

.. _AMRShockPool2D:

AMRShockPool2D
^^^^^^^^^^^^^^

2D Shock Propagation Test (AMR Version)

.. _Athena-RayleighTaylor:

Athena-RayleighTaylor
^^^^^^^^^^^^^^^^^^^^^
 
classic Raleigh Taylor setup with sharp contact 
this file should work with all hydro methods

compare to 
http://www.astro.princeton.edu/~jstone/tests/rt/rt.html

.. _DoubleMachReflection:

DoubleMachReflection
^^^^^^^^^^^^^^^^^^^^
Double Mach Reflection test (see WC84, Section IVc)

Mach = 10 shock in air at 60 degrees angle with respect to a reflecting wall.
The gas density ahead of the shock is 1.4, and the pressure is 1.0; density
behind the shock is 8.0.

Compare with WC84 Figure 4 at t = 2.0. Note, the diffusion is OFF here;
also note the figure shows only a part of the domain: 0 < x < 3.

Most of the required parameters were hardwired, see DoubleMachInitialize.C.


.. _FreeExpansionAMR:

FreeExpansionAMR
^^^^^^^^^^^^^^^^
(John Wise, July 2010)

This test simulates a blastwave in the free expansion phase.  In the
initial setup, the interior region has a uniform density and a
linearly increasing radial velocity.  The blastwave should advect
outwards, and create a high entropy shell and have little oscillations
in the shock.  It is exactly the same as the FreeExpansion 1D test,
but in two dimensions and with 2 levels of AMR.  There is usually some
non-spherical in the reverse shock, but the main outer shock should
always be nearly spherical.  The artifacts in the reverse shock
originate from the initial discretization and decreases if the
resolution is increased.

This test runs to completion and creates 21 outputs.  This test
problem is new for version 2.0.  It uses the new hydro_rk solver.

The initial setup is taken from Truelove & McKee, 1999, ApJS, 120,
299.

.. _HDMHD2DCheckOddEvenCouplingOfRiemannSolver:

HDMHD2DCheckOddEvenCouplingOfRiemannSolver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Standing shock with slight density perturbation. 
Look at the y-velocity to check for the effect. 
It is subtle because we use a small perturbation.
define problem

HydroMethod = 3/RiemanSolver = 4 (HLLC) shows this strongly, so does PPM (HM=0)
HM=3 wit RiemanSolver 1 (HLL) or 3 (LLF) pass this better and so does Zeus

.. _Implosion:

Implosion
^^^^^^^^^

A 2D converging shock test problem.

Liska & Wendroff, 2003, SIAM J. Sci. Comp., 25, N3, 995-1017
http://www-troja.fjfi.cvut.cz/~liska/CompareEuler/compare8

Jim Stone's Athena test page
http://www.astro.princeton.edu/~jstone/tests/implode/Implode.html

.. _ImplosionAMR:

ImplosionAMR
^^^^^^^^^^^^
AMR Version of Implosion test

.. _KelvinHelmholtz:

KelvinHelmholtz
^^^^^^^^^^^^^^^
The KH Test problem creates two fluids moving antiparallel to each other
in a periodic 2D grid (inner fluid and outer fluid).  The inside fluid
has a higher density than the outside fluid.  There is a slight ramp region
in density and x-velocity connecting the two regions so there are no 
discontinuities in the flow.  The y-velocity is perturbed with small sinusoidal
perturbation.  As the flows shear past each other, the KH instability
is excited, which develops over time.  This test watches the evolution of
those instabilities.  --Cameron Hummels, 2013



.. _KelvinHelmholtzAMR:

KelvinHelmholtzAMR
^^^^^^^^^^^^^^^^^^
This version incorporates 1 level of AMR using the shear criterion on the
interface between the two fluids.  --Cameron Hummels, 2013


.. _NohProblem2D:

NohProblem2D
^^^^^^^^^^^^
NohProblem2D (See Noh (1987) J. Comp. Phys. 72, 78)

The Noh Problem test sets up a a uniform gas of density of 1.0 that
has a uniform inward radial velocity of 1.0

Noh (1987) J. Comp. Phys. 72, 78 introduced an infinite shock
reflection problem that has an exact analytical solution. Gas with
initially uniform density of 1 and zero pressure converges onto the
origin with a uniform radial velocity of 1 in a unit domain x, y \in
[0, 1].

In cylindrical geometry (2D) the solution is an infinite strength
circularly symmetric shock reflecting from the origin. In the
postshock region the density is 16, the velocity is zero, and the
pressure is 16/3. The shock speed is 1/3, and in the preshock region
the density varies as a function of time and radius as (1 + t/sqrt(x^2
+ y^2)) while velocity and pressure keep their initial values.

We set the initial pressure to be 10^-6 instead of zero for numerical
reasons. We use reflecting boundaries at x=0 and at y=0 and set up the
outer boundaries at x=1 and y=1 based on the exact analytical
solution.  We follow the propagation of the shock until it reaches a
radius of 2/3 at t=2. At this point we compare our results with
similar tests performed with other Eulerian numerical schemes, see
Liska & Wendroff (2003), Section 4.5 in this PDF document or in SIAM
J. Sci. Comput. 25, 995, 2003. See also Rider (2000),
J. Comp. Phys. 162, 395 for a discussion of "wall heating" phenomenon
near the origin that seriously affects the results obtained with
Lagrangian schemes. Numerically, this is a difficult problem.

.. _NohProblem2DAMR:

NohProblem2DAMR
^^^^^^^^^^^^^^^

AMR Version of NohProblem2D

.. _RadiatingShockWave:

RadiatingShockWave
^^^^^^^^^^^^^^^^^^

A 2D explosion test problem which includes radiative cooling.

.. _RampedKelvinHelmholtz2D:

RampedKelvinHelmholtz2D
^^^^^^^^^^^^^^^^^^^^^^^

Kelvin Helmholtz with a ramp

.. _SedovBlast:

SedovBlast
^^^^^^^^^^
A 2d explosion test problem

.. _SedovBlastAMR:

SedovBlastAMR
^^^^^^^^^^^^^
AMR version of SedovBlast 
.. _ShockPool2D:

ShockPool2D
^^^^^^^^^^^

2D Shock Propagation Test

.. _ValidatedNonlinearKelvinHelmholtz:

ValidatedNonlinearKelvinHelmholtz
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Implements initial condition generator for 
http://arxiv.org/abs/1509.03630
A Validated Nonlinear Kelvin-Helmholtz Benchmark for Numerical Hydrodynamics
Daniel Lecoanet, Michael McCourt, Eliot Quataert, Keaton J. Burns, 
Geoffrey M. Vasil, Jeffrey S. Oishi, Benjamin P. Brown, James M. Stone, Ryan M. O'Leary (2015)

Note, this has not been tested much. 
It uses the machinery for the MHD2D tests in hydro_rk. 
It should also work with PPM (i.e. HydroMethod = 0) but is not setup to do Zeus (HydroMethod=1 tests). 
However, HydroMethod=3 offers many options for Riemann Solvers and Slope limiting which may be interesting to test. 

To run in parallel keep the ParallelRootgridIO = 1 on.

This is a code test for a limit high resolution shock captruing codes are not often optimized for.
It is useful to see whether your choice of hydro method has sufficient diffusivity to give behave sensibly in
a convergence study. 


.. _Hydro/Hydro-3D:

Hydro/Hydro-3D
~~~~~~~~~~~~~~

.. _AgoraGalaxy:

AgoraGalaxy
^^^^^^^^^^^
The initial conditions for this problem depend on star and dark matter particle
data files that must be generated using the MakeGalaxy code.

Example star particle data files corresponding to the AGORA isolated galaxy
initial conditions can be found at the following URL:

https://www.dropbox.com/sh/1xzt1rysy9v3a9l/AAAHZyjrfTz88aG12H0Q_Rqla

The included test problem parameter file corresponds to the LOW initial
conditions (80 pc resolution).  Follow the readme.txt and use the IPython
notebook to adjust the parameter file if you want to simulate a galaxy at
different resolution.

Enzo must be compiled with new-problem-types-yes and grackle-yes for this test
problem to function correctly. The bash script prepare_sim.sh downloads and
extracts the necessary grackle data file and initial conditions. Run it in the
directory you would like to run the test simulation.

If run to completion, this simulation needs about 200 gigabytes of free space.
If you are constrained by hard disk space, consider decreasing the output
interval (1 Myr of physical time per output by default).  The simulation should
finish in approximately 12 hours when run on 16 cores.  Of course, your mileage
may vary for the amount of wallclock time necessary.
.. _Athena-RayleighTaylor3D:

Athena-RayleighTaylor3D
^^^^^^^^^^^^^^^^^^^^^^^
classic Raleigh Taylor setup with sharp contact 
this file should work with all hydro methods

compare to 
http://www.astro.princeton.edu/~jstone/tests/rt/rt.html

.. _CollapseTestNonCosmological:

CollapseTestNonCosmological
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test problem initializes a constant density sphere at the center of the box 
with radius of 0.15 of the box width.  The sphere is in pressure equilibrium with 
its surroundings.  The sphere will collapse dynamically and reach its peak density 
at t ~ 5.2 in code units.  Radiative cooling is turned off, so the sphere will bounce 
and reach an equilibrium state.

This test runs to completion (t = 7) and produces 71 outputs in roughly an hour 
on a single processor.  In plot.py, we plot the evolution of the peak density versus 
time create density projections for each output.

Note for comparing with Enzo 1.5 - By default, Enzo 2.0 performs 4 iterations of 
the multigrid gravity solver (configurable with the PotentialIterations parameter).  
However, the default number of potential iteration in Enzo 1.5 is 0 and there is 
no parameter to change this.  For a direct comparison, Enzo 1.5 must be compiled 
with MAX_POTENTIAL_ITERATIONS (in macros_and_parameters.h) set to 4.
.. _CollideTest:

CollideTest
^^^^^^^^^^^


.. _ExtremeAdvectionTest:

ExtremeAdvectionTest
^^^^^^^^^^^^^^^^^^^^

.. _GalaxySimulation:

GalaxySimulation
^^^^^^^^^^^^^^^^
This problem sets up a galaxy disk using code written by Stephanie
Tonnesen and Munier Salem.  The gravity is modelled with a static
potential that includes a dark matter halo, bulge, and stellar disk.
The gas disk is modelled with self gravity and the mass, scale height
and length can be set with parameters.  A full description is in Salem
et al 2014 (ApJ, 815, 77; http://arxiv.org/pdf/1507.07935.pdf).

This version also includes a ram pressure wind that uses a lookup
table (included here as ICMinflow_data.in) to set the density,
temperature and velocity of the inflowing gas as a function of time
(see documentation and the above paper).

This is not really intended as a test (and so one should not trust the
enzotest settings), but more as an example run file for this kind of
problem.


.. _NFWCoolCoreCluster:

NFWCoolCoreCluster
^^^^^^^^^^^^^^^^^^
The NFW Cool-Core Cluster is a simulation of the cooling flow in an idealized
cool-core cluster that resembles Perseus cluster. It can be a test for cooling, and
maybe gravity too. 

The default set up is with a root grid of 64^3, a maximum refinement level of 12,
and MinimumMassForRefinementLevelExponent of -0.2 (for better resolution, use -1.2)
which can be changed based on the resolution one needs.

The default set up has a static gravity and no self-gravity of the gas since the latter
is much smaller than the gravity of the dark matter and does not change much during the 
early stage of the cooling flow evolution.

As the cooling catastrophe happens, the temperature drops to the bottom of the cooling 
function ~ 10^4 K in the center within ~1kpc with the default resolution. The size of 
the region becomes smaller with higher resolution.

The projected gas density shows a disk of size ~ 1kpc (inside the cooling catastrophe
region) at late times along z axis which is the direction of 
the initial angular momentum of the gas.  

.. _NohProblem3D:

NohProblem3D
^^^^^^^^^^^^
The Noh Problem test sets up a a uniform gas of density of 1.0 that
has a uniform inward radial velocity of 1.0

Noh (1987) J. Comp. Phys. 72, 78 introduced an infinite shock
reflection problem that has an exact analytical solution. Gas with
initially uniform density of 1 and zero pressure converges onto the
origin with a uniform radial velocity of 1 in a unit domain x, y \in
[0, 1].

In spherical geometry (3D) the postshock density is 64 and in the
preshock region the density varies as (1 + t/sqrt(x^2 + y^2))^2. All
other dependencies remain the same as in the 2D case.

We set the initial pressure to be 10^-6 instead of zero for numerical
reasons. We use reflecting boundaries at x=0 and at y=0 and set up the
outer boundaries at x=1 and y=1 based on the exact analytical
solution.  We follow the propagation of the shock until it reaches a
radius of 2/3 at t=2. At this point we compare our results with
similar tests performed with other Eulerian numerical schemes, see
Liska & Wendroff (2003), Section 4.5 in this PDF document or in SIAM
J. Sci. Comput. 25, 995, 2003. See also Rider (2000),
J. Comp. Phys. 162, 395 for a discussion of "wall heating" phenomenon
near the origin that seriously affects the results obtained with
Lagrangian schemes. Numerically, this is a difficult problem.


.. _NohProblem3DAMR:

NohProblem3DAMR
^^^^^^^^^^^^^^^
AMR version of NohProblem3D

.. _ProtostellarCollapse_Std:

ProtostellarCollapse_Std
^^^^^^^^^^^^^^^^^^^^^^^^

.. _RotatingCylinder:

RotatingCylinder
^^^^^^^^^^^^^^^^
The rotating cylinder is a test of the conservation of the
angular momentum in Enzo.

The test sets up a rotating gas cyclinder in the center of the box at 0.5 0.5
0.5 with an overdensity of 200. The default set up is with a root grid
of 32^3 and a single level of refinement. 

The cyclinder collapses towards its center and oscillates before
settling. 

The results signify correctness if the net change in the total angular
momentum of the system is low (< 5%). The percentage change per output
should be less than 1% and should decrease over time as the collapse
reaches equilibrium. 

In plots.py, we image slices in the x-direction and plot the angular
momentum evolution of the system. 

.. _RotatingSphere:

RotatingSphere
^^^^^^^^^^^^^^

.. _ShockPool3D:

ShockPool3D
^^^^^^^^^^^

A 3D Shock Propagation Test

.. _StripTest:

StripTest
^^^^^^^^^

.. _MHD/1D:

MHD/1D
~~~~~~~~~~

.. _BrioWu-MHD-1D:

BrioWu-MHD-1D
^^^^^^^^^^^^^
From
 Brio, M., & Wu, C. C. 1988, J. Comput. Phys., 75, 400
 Wang, P., & Abel, T. 2009, Astrophysical Journal, 696:96-109

Run: Ji-hoon Kim, July 2010

This test sets up an one-dimensional Riemann problem for MHD, and has become a useful 
standard test for any MHD solver.  Detailed description of the initial set up can be 
found in the papers above.  This test problem is new for enzo2.0.  It uses the new 
Stanford hydro_rk solver.

This test runs to completion while generating 12 outputs, and scripts.py will 
produce the plots for Density, x-velocity, By, Internal Energy for the last snapshot 
(t=0.1).  This last snapshot should be compared to figure 18 from Figure 2 of Brio & 
Wu (1988) or Figure 15 of Wang & Abel (2009)

Success in test_briowu.py is determined by nearly exact match (5e-3) in Density and By. 

.. _BrioWu-MHD-1D-MHDCT:

BrioWu-MHD-1D-MHDCT
^^^^^^^^^^^^^^^^^^^

This also serves as an Example of how to do 1D HD/MHD tests with the myriad 
of shock tube problems defined in the literature

.. _CR-ShockTube:

CR-ShockTube
^^^^^^^^^^^^

This is a fairly mild test. The solution consists of left
rarefaction wave, a contact discontinuity, and a right shock.

See Pfrommer et al 2006 for information on the analytic sol'tn

.. _MHD_Metal_Advection_CT:

MHD_Metal_Advection_CT
^^^^^^^^^^^^^^^^^^^^^^
Square wave advection with a single species field.
The metals are offset by 0.25 from the density. 
python plot.py will make a plot.  

.. _MHD_Metal_Advection_Dedner:

MHD_Metal_Advection_Dedner
^^^^^^^^^^^^^^^^^^^^^^^^^^
Square wave advection with a single species field.
The metals are offset by 0.25 from the density.

.. _MHD/2D:

MHD/2D
~~~~~~~~~~

.. _LoopAdvection_CT:

LoopAdvection_CT
^^^^^^^^^^^^^^^^

Advection of a magnetic field loop.
Originally described by Gardiner & Stone 2005 (Journal of Computational Physics,
205,509).  The multidimensional MHD analogue of a square wave advection test,
the field loop severely deforms for many CT schemes.  

.. _LoopAdvection_Dedner:

LoopAdvection_Dedner
^^^^^^^^^^^^^^^^^^^^
Advection of a magnetic field loop.
Originally described by Gardiner & Stone 2005 (Journal of Computational Physics,
205,509).  The multidimensional MHD analogue of a square wave advection test,
the field loop severely deforms for many CT schemes. 

.. _MHD2DRotorTest:

MHD2DRotorTest
^^^^^^^^^^^^^^
From
 G. Toth, J. Comput. Phys. 161 (2000) 605

Initially discussed in 
 D. Balsara & D. Spicer, J. Comput. Phys. 149, 270292 (1999) 

Run: dcollins, July 2010

A cylander of gas is set rotating in a uniform medium, with uniform magnetic
field perpandicular to the rotation axis.  A torsional Alfven wave is launched
from the surface, and propogates to the simulation boundary at the end of the
simulation.

This test is useful as it can and will cause negative densities for some
solvers; it also demonstrates physics unique to the MHD system, providing a
useful comparison.

Visual comparison to the plots in Toth (2000) shows that the general morphology
and extent of the Alfven wave is similar.  Round contours in the mach number
field indicate solid body rotation, without the artifacts seen in some other
solvers.  

Success in test_rotor.py is determined by nearly exact match (1e-12) in L1 norm between
Density, Bx, P, and Mach number. 

This test generates 11 outputs, and snapshots for the 4 above fields for each
snapshot.  The 11th snapshot should be compared to figure 18 from Toth (2000)

.. _MHDCTOrszagTang:

MHDCTOrszagTang
^^^^^^^^^^^^^^^

.. _MHDCTOrszagTangAMR:

MHDCTOrszagTangAMR
^^^^^^^^^^^^^^^^^^

.. _MHDDednerOrszagTang:

MHDDednerOrszagTang
^^^^^^^^^^^^^^^^^^^

.. _RayleighTaylor_CT_Suppressed:

RayleighTaylor_CT_Suppressed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MHD suppresses the Rayleigh Taylor instability (in 2D)
This is set with a magnetic field that is just slightly below the critical value
for complete suppression, but one does see a severe reduction in growth rate and
secondary instability.

Fun things to try include reducing the field strength and changing the
direction!  

.. _SedovBlast-MHD-2D-Fryxell:

SedovBlast-MHD-2D-Fryxell
^^^^^^^^^^^^^^^^^^^^^^^^^
From
 Fryxell et al, 2000, ApJS, 131, 273 (Section 7.4)

Run: Ji-hoon Kim, July 2010


This test sets up a two-dimensional blast wave problem for MHD.  While the initial 
condition  essentially describes a circular overpressurized region in a low-pressure 
magnetized medium, more detailed description of the initial set up can be found in 
the papers above.  This test problem is new for enzo2.0.  It uses the new Stanford 
hydro_rk solver.  

Unfortunately for some users, some of the key parameters are hard-coded in 
Grid_MHD2DTestInitializeGrid.C.  Depending on B-field values, this test can be a 
pure Sedov blast wave problem, or a Gardiner blast wave problem.  The current 
parameter will give Sedov blast wave solution. 

(1) Setting LowerBx=LowerBy=0 will give the traditional Sedov test 
    with no B-field, essentially very similar to SedovBlast-AMR.enzo; 
    but with different hydro solver (Stanford HD/MHD)
(2) Setting LowerBx=LowerBy=5 will give the Sedov blast with the 
    presence of B-field, very similar to SedovBlast-MHD-2D-Gardiner.enzo;
    but the exact setup is somewhat different (see the code)

This test runs to completion while generating 12 outputs, and scripts.py will 
produce the plots and slices for Density and Pressure of the last snapshot (t=0.05).  
This last snapshot should be compared to Figure 29 and 30 of Fryxell (2000).

Success in test_fryxell.py is determined by nearly exact match (3e-5) in Density and 
Pressure. 

.. _SedovBlast-MHD-2D-Gardiner:

SedovBlast-MHD-2D-Gardiner
^^^^^^^^^^^^^^^^^^^^^^^^^^
From
 Gardiner, T. A., & Stone, J. M. 2005, J. Comput. Phys., 205, 509
 Wang, P., & Abel, T. 2009, Astrophysical Journal, 696:96-109

Run: Ji-hoon Kim, July 2010

This test sets up a two-dimensional blast wave problem for MHD, and has become a useful 
standard test for any MHD solver.  While the initial condition  essentially describes 
a circular overpressurized region in a low-pressure magnetized medium, more detailed 
description of the initial set up can be found in the papers above.  This test problem 
is new for enzo2.0.  It uses the new Stanford hydro_rk solver.  

Unfortunately for some users, most of the key parameters are hard-coded in 
Grid_MHD2DTestInitializeGrid.C.  With zero B-field, this test should be a pure Sedov 
blast wave problem.    

This test runs to completion while generating 12 outputs, and scripts.py will 
produce the plots for Density, x-velocity, By, Internal Energy for the last snapshot 
(t=0.02).  This last snapshot should be compared to Figure 13 of Gardiner & Stone (2005)
or Figure 16 of Wang & Abel (2009)

Success in test_gardiner.py is determined by nearly exact match (3e-5) in Density, 
Pressure, Bx, and By. 

.. _Wengen2-CollidingFlow:

Wengen2-CollidingFlow
^^^^^^^^^^^^^^^^^^^^^
Wengen 2 colliding flow
Reference: http://www-theorie.physik.unizh.ch/~agertz/Wengen_2/Code_tests.html
Tom Abel September 2010

Also works with magnetic fields. 

.. _MHD/3D:

MHD/3D
~~~~~~~~~~

.. _ShearingBox:

ShearingBox
^^^^^^^^^^^

.. _StochasticForcing:

StochasticForcing
^^^^^^^^^^^^^^^^^
MHD/HD turbulence problem with stochastic forcing with subgrid-scale (SGS) turbulence model
Philipp Grete 2014

Typical "turbulence-in-a-box" problem with non-static driving field.
For details on stochastic forcing, see Schmidt et al. 2009 A&A 494, 127-145 
http://dx.doi.org/10.1051/0004-6361:200809967
For details on the SGS model, see Grete et al. (2017) Phys. Rev. E. 95 033206
https://dx.doi.org/10.1103/PhysRevE.95.033206

Works/properly tested only on 3D uniform grids with MUSCL type solvers and MHDCT at this point.
For hydro use HydroMethod 3
For MHD use HydroMethod 4
For MHDCT use HydroMethod 6

.. _RadiationTransport:

RadiationTransport
~~~~~~~~~~~~~~~~~~

.. _PhotonShadowing:

PhotonShadowing
^^^^^^^^^^^^^^^
This problem tests shadowing capabilities of the ray tracing module
with an optically-thick clump absorbing radiation from a source at the
opposite edge.  It has the same parameters (at lower resolution but
with AMR) as Test 3 in Iliev et al. (2006) MNRAS, 371, 1057.  The
ambient medium is optically-thin and the radiation should immediately
hit the clump and start ionizing and heating it.  The ionization front
should reach slightly past halfway in the clump at the end of the
simulation.  In this setup, there is a analytic solution and can be
found in Shapiro et al. (2004), MNRAS, 348, 753.


This test produces 15 outputs at intervals of 1 Myr.  The analysis
script (1) calculates the average temperature and ionized fraction of the
clump in each output, (2) gives line cuts through the clump center at
1,3,5,10,15 Myr, and (3) produces slices through the clump center at
the final time.

.. _PhotonTest:

PhotonTest
^^^^^^^^^^
** Test 1 from Iliev et al. (2006), MNRAS, 371, 1057 **

 - Source at the origin
 - Luminosity = 5e48 ph/s
 - Fixed temperature, no hydro
 - Density = 1e-3 cm^-3

.. _PhotonTestAMR:

PhotonTestAMR
^^^^^^^^^^^^^
This problem is the classical HII region expansion in an isothermal
sphere with a constant density core with a point source at the origin.
The source has a T=10^5 K blackbody spectrum with 1 energy group.  It
has the same parameters (at a lower base resolution, but with AMR) as
Test 6 in Iliev et al. (2009) MNRAS, 400, 1283.  The source should
heat and ionize a spherical region around it and drive a shock
outwards.  There is no analytical solution to this problem but has
been studied in great detail.

This test produces 25 outputs at intervals of 1 Myr.  The analysis
script (1) finds the ionization front radius for each output, (2)
creates radial profiles at a few times, (3) produces slices at the
origin in the final output, and (4) computes the deviation in the
photo-ionization rates from 1/r^2.

.. _PhotonTestMultiFrequency:

PhotonTestMultiFrequency
^^^^^^^^^^^^^^^^^^^^^^^^
This test is derived from the PhotonTest problem which is itself
based on the classical HII region expansion in a uniform
static medium similar to Test 1 in Iliev et al. (2006). Rather than 
a mono-chromatic source this source contains 7 frequency bins from 0.5 eV 
up to 100 eV. It can be used to test the multi-frequency photon solver. In a separate
subdirectory is the optically thin version (OT) which runs with H2 and HM dissociation
run in optically thin mode.

This test runs with Grackle. 

.. _RadiationTransportFLD:

RadiationTransportFLD
~~~~~~~~~~~~~~~~~~~~~

.. _CosmoIonization_q05z10:

CosmoIonization_q05z10
^^^^^^^^^^^^^^^^^^^^^^

.. _CosmoIonization_q05z10_sp:

CosmoIonization_q05z10_sp
^^^^^^^^^^^^^^^^^^^^^^^^^

.. _CosmoIonization_q05z4:

CosmoIonization_q05z4
^^^^^^^^^^^^^^^^^^^^^

.. _CosmoIonization_q05z4_sp:

CosmoIonization_q05z4_sp
^^^^^^^^^^^^^^^^^^^^^^^^

.. _CosmoIonization_q5z10:

CosmoIonization_q5z10
^^^^^^^^^^^^^^^^^^^^^

.. _CosmoIonization_q5z10_sp:

CosmoIonization_q5z10_sp
^^^^^^^^^^^^^^^^^^^^^^^^

.. _CosmoIonization_q5z4:

CosmoIonization_q5z4
^^^^^^^^^^^^^^^^^^^^

.. _Grey_Enzochem:

Grey_Enzochem
^^^^^^^^^^^^^
.. _Grey_Split:

Grey_Split
^^^^^^^^^^

.. _RadiatingShockLab:

RadiatingShockLab
^^^^^^^^^^^^^^^^^

.. _RadiatingShockLab1D:

RadiatingShockLab1D
^^^^^^^^^^^^^^^^^^^
.. _RadiatingShockLab1D_sp:

RadiatingShockLab1D_sp
^^^^^^^^^^^^^^^^^^^^^^
.. _RadiationStream1D:

RadiationStream1D
^^^^^^^^^^^^^^^^^

.. _RadiationStream1D_sp:

RadiationStream1D_sp
^^^^^^^^^^^^^^^^^^^^
.. _RadiationStreamX0:

RadiationStreamX0
^^^^^^^^^^^^^^^^^
.. _RadiationStreamX1:

RadiationStreamX1
^^^^^^^^^^^^^^^^^
.. _RadiationStreamX1_sp:

RadiationStreamX1_sp
^^^^^^^^^^^^^^^^^^^^
.. _RadiationStreamY0:

RadiationStreamY0
^^^^^^^^^^^^^^^^^

.. _RadiationStreamY0_sp:

RadiationStreamY0_sp
^^^^^^^^^^^^^^^^^^^^
.. _RadiationStreamY1:

RadiationStreamY1
^^^^^^^^^^^^^^^^^
.. _RadiationStreamY1_sp:

RadiationStreamY1_sp
^^^^^^^^^^^^^^^^^^^^
.. _RadiationStreamZ0:

RadiationStreamZ0
^^^^^^^^^^^^^^^^^
.. _RadiationStreamZ0_sp:

RadiationStreamZ0_sp
^^^^^^^^^^^^^^^^^^^^
.. _RadiationStreamZ1:

RadiationStreamZ1
^^^^^^^^^^^^^^^^^
.. _RHIonization1:

RHIonization1
^^^^^^^^^^^^^
.. _RHIonization1_sp:

RHIonization1_sp
^^^^^^^^^^^^^^^^
.. _RHIonization2:

RHIonization2
^^^^^^^^^^^^^
.. _RHIonization2_sp:

RHIonization2_sp
^^^^^^^^^^^^^^^^

.. _TurnerStoneEquil1:

TurnerStoneEquil1
^^^^^^^^^^^^^^^^^

.. _TurnerStoneEquil2:

TurnerStoneEquil2
^^^^^^^^^^^^^^^^^

.. _StarParticle:

StarParticle
~~~~~~~~~~~~

This test places a single star particle in the center of a box with uniform
gas density and thermal energy.  The gas is initially at rest.  The particle 
will then produce feedback according to the method set for 
StarParticleFeedback.  

By default, the star particle produces feedback with method 14, kinetic 
feedback.  An inital timestep is set to account for the large bulk velocities 
created by the feedback event in the first timestep.  

The particle is also set to be at rest by default, but it can be given a motion
relative to the gas with the parameter TestStarParticleStarVelocity = vx vy vz, 
where vx, vy and vz have units of km/s.



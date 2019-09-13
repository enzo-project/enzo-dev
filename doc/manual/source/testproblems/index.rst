.. _parameters:

Enzo Test Problems
===================

The following is a list of the test problems that are available in the 
run directory found in the Enzo repository. These problems are also the
same set of problems used in the test suite. They are listed in broad
categories following the directory structure.

This list however is not a complete list of all the ProblemTypes that 
are available in Enzo.

The test problem specific  parameters can be found in the Parameter List.

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

   * `SedovBlasti`_

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

* `MHD/MHD-1D`_

   * `BrioWu-MHD-1D`_

   * `BrioWu-MHD-1D-MHDCT`_

   * `CR-ShockTube`_

   * `MHD_Metal_Advection_CT`_

   * `MHD_Metal_Advection_Dedner`_

* `MHD/MHD-2D`

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

* `MHD/MHD-3D`

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


.. _Cooling:

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
CoolingTest_Cloudy.enzo - uses Cloudy cooling along with the 
			       	      	    MultiSpecies = 1 chemistry.  The input data 
					    provided is a three dimensional table that 
					    varies in density, metallicity, and temperature.

Cooling data files:
primordial_cie.dat - CIE cooling rates for atomic H and He taken from 
		     	 	 Black (1981).
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

.. _Cosmology:

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

.. _MHDAdiabaticExpansion_Dedner:

MHDAdiabaticExpansion_Dedner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. _MHDZeldovichPancake:

MHDZeldovichPancake
^^^^^^^^^^^^^^^^^^^


.. _MHDZeldovichPancake_2_CT:   

MHDZeldovichPancake_2_CT
^^^^^^^^^^^^^^^^^^^^^^^^

.. _MHDZeldovichPancake_2_Dedner:

MHDZeldovichPancake_2_Dedner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _SphericalInfall:

SphericalInfall
^^^^^^^^^^^^^^^


.. _ZeldovichPancake:

ZeldovichPancake
^^^^^^^^^^^^^^^^

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

   * `SedovBlasti`_

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

* `MHD/MHD-1D`_

   * `BrioWu-MHD-1D`_

   * `BrioWu-MHD-1D-MHDCT`_

   * `CR-ShockTube`_

   * `MHD_Metal_Advection_CT`_

   * `MHD_Metal_Advection_Dedner`_

* `MHD/MHD-2D`

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

* `MHD/MHD-3D`

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






File naming conventions and routine locations
=============================================

The large number of source files can be intimidating even to the
experienced Enzo developer, and this page describes some of the naming
conventions.  Familiarity with ``grep`` or ``ack`` and pipes like ``ls
-1 | grep`` are essential.  Also routines with a similar functionality
are grouped together with a common name.

Here are some file naming rules that are used.

1. Internal capitalization is used for C files, all lowercase with
   underscores for fortran files and header files. All Fortran files
   end with .F.

2. With very few exceptions, Enzo has a one function per file layout, with the
   file name being the function name. 

3. Object methods have the object name prepended to the beginning,
   such as the member of the grid class ``SolveHydroEquations`` lives
   in the file ``Grid_SolveHydroEquations.C``.

Below we list some examples of filenames, grouped by functionality.
*This is not a complete list of files in Enzo.*

Grid methods
------------

Initializers
^^^^^^^^^^^^
::

	Grid_NohInitializeGrid.C
	Grid_OneZoneFreefallTestInitializeGrid.C
	Grid_NestedCosmologySimulationInitializeGrid.C
	Grid_RHIonizationTestInitializeGrid.C
	Grid_RadHydroStreamTestInitializeGrid.C
	Grid_RadiatingShockInitializeGrid.C


Particles
^^^^^^^^^
::

	Grid_AddParticlesFromList.C
	Grid_MoveAllParticles.C
	Grid_MoveSubgridParticlesFast.C
	Grid_TracerParticleCreateParticles.C
	Grid_TracerParticleOutputData.C
	Grid_TracerParticleSetVelocity.C
	Grid_TransferSubgridParticles.C
	Grid_TransferSubgridStars.C
	Grid_UpdateParticlePosition.C
	Grid_UpdateParticleVelocity.C

Solvers
^^^^^^^
::

	Grid_ComputeCoolingTime.C
	Grid_ComputeDustTemperatureField.C
	Grid_ComputeGammaField.C
	Grid_ComputePressure.C
	Grid_MultiSpeciesHandler.C
	Grid_SolveHydroEquations.C
	Grid_SolvePPM_DE.C
	Grid_SolvePrimordialChemistryCVODE.C
	Grid_SolveRadiativeCooling.C
	Grid_SolveRateAndCoolEquations.C
	Grid_SolveRateEquations.C
	Grid_ZeusSolver.C
	Grid_xEulerSweep.C
	Grid_yEulerSweep.C
	Grid_zEulerSweep.C

Gravity and acceleration
^^^^^^^^^^^^^^^^^^^^^^^^
::

	Grid_AddBaryonsToGravitatingMassField.C
	Grid_AddExternalAcceleration.C
	Grid_AddExternalPotentialField.C
	Grid_ComputeAccelerationField.C
	Grid_ComputeAccelerationFieldExternal.C
	Grid_ComputeAccelerations.C
	Grid_ComputeAccelerationsFromExternalPotential.C
	Grid_DepositBaryons.C
	Grid_DepositMustRefineParticles.C
	Grid_DepositParticlePositions.C
	Grid_PrepareFFT.C
	Grid_PrepareGreensFunction.C
	Grid_PreparePotentialField.C
	Grid_SolveForPotential.C

Hierarchy control
^^^^^^^^^^^^^^^^^
::

	Grid_AddFieldMassToMassFlaggingField.C
	Grid_AddOverlappingParticleMassField.C
	Grid_AllocateGrids.C
	Grid_CopyZonesFromGrid.C
	Grid_FlagCellsToAvoidRefinement.C
	Grid_FlagCellsToAvoidRefinementRegion.C
	Grid_FlagCellsToBeRefinedByCoolingTime.C
	Grid_FlagCellsToBeRefinedByJeansLength.C
	Grid_FlagCellsToBeRefinedByMass.C
	Grid_SetFlaggingField.C
	Grid_SetFlaggingFieldStaticRegions.C

Utilities
^^^^^^^^^
::

	Grid_AccessBaryonFields.C
	Grid_ComputeTemperatureField.C
	Grid_IdentifyColourFields.C
	Grid_IdentifyGloverSpeciesFields.C
	Grid_IdentifyNewSubgrids.C
	Grid_IdentifyNewSubgridsSmall.C
	Grid_IdentifyPhysicalQuantities.C
	Grid_IdentifyRadiationPressureFields.C
	Grid_IdentifyRadiativeTransferFields.C
	Grid_IdentifyShockSpeciesFields.C
	Grid_IdentifySpeciesFields.C

Conduction
^^^^^^^^^^
::

	Grid_ConductHeat.C
	Grid_ConductionBubbleInitialize.C
	Grid_ConductionCloudInitialize.C
	Grid_ConductionTestInitialize.C


Radiation
^^^^^^^^^
::

	Grid_AddH2Dissociation.C
	Grid_AddRadiationImpulse.C
	Grid_AddRadiationPressureAcceleration.C
	Grid_AllocateInterpolatedRadiation.C
	Grid_ComputePhotonTimestep.C
	Grid_ComputePhotonTimestepHII.C
	Grid_ComputePhotonTimestepTau.C
	Grid_FinalizeRadiationFields.C
	Grid_PhotonPeriodicBoundary.C
	Grid_PhotonSortLinkedLists.C
	Grid_SetSubgridMarkerFromParent.C
	Grid_SetSubgridMarkerFromSibling.C
	Grid_SetSubgridMarkerFromSubgrid.C
	Grid_Shine.C

I/O
^^^
::

	New_Grid_ReadGrid.C
	New_Grid_WriteGrid.C
	Grid_WriteNewMovieData.C
	Grid_WriteNewMovieDataSeparateParticles.C


Communcation
^^^^^^^^^^^^
::

	Grid_CommunicationMoveGrid.C
	Grid_CommunicationReceiveRegion.C
	Grid_CommunicationSendParticles.C
	Grid_CommunicationSendPhotonPackages.C
	Grid_CommunicationSendRegion.C
	Grid_CommunicationSendStars.C
	Grid_CommunicationTransferParticlesOpt.C
	Grid_CommunicationTransferStarsOpt.C


Feedback
^^^^^^^^
::

	Grid_ChangeParticleTypeBeforeSN.C
	Grid_AddFeedbackSphere.C
	Grid_FindNewStarParticles.C

Analysis
^^^^^^^^
::

	Grid_CalculateAngularMomentum.C
	Grid_ConvertToNumpy.C

Turbulence
^^^^^^^^^^
::

	Grid_AddRandomForcing.C
	Grid_AppendForcingToBaryonFields.C
	Grid_ComputeRandomForcingFields.C
	Grid_DetachForcingFromBaryonFields.C
	Grid_PrepareRandomForcingNormalization.C
	Grid_ReadRandomForcingFields.C
	Grid_RemoveForcingFromBaryonFields.C


Communication methods
---------------------
::

Core methods
------------
::

External boundary methods
-------------------------
::

Halo finder methods
-------------------
::

Gravity methods
---------------
::

Hierarchy methods
-----------------
::

Radiation methods
-----------------
::

Flux limited diffusion
^^^^^^^^^^^^^^^^^^^^^^
::

Adaptive ray tracing
^^^^^^^^^^^^^^^^^^^^
::

I/O
---
::


Star formation methods
----------------------
::

Utilities
---------
::

Headers
-------
::

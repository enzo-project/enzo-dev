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

.. contents::

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

	CommunicationBroadcastValue.C
	CommunicationBufferedSend.C
	CommunicationCollectParticles.C
	CommunicationCombineGrids.C
	CommunicationInitialize.C
	CommunicationLoadBalanceGrids.C
	CommunicationLoadBalancePhotonGrids.C
	CommunicationLoadBalanceRootGrids.C
	CommunicationMergeStarParticle.C
	CommunicationNonblockingRoutines.C
	CommunicationParallelFFT.C
	CommunicationPartitionGrid.C
	CommunicationReceiveFluxes.C
	CommunicationReceiveHandler.C
	CommunicationReceiverPhotons.C
	CommunicationSendFluxes.C
	CommunicationShareGrids.C
	CommunicationShareParticles.C
	CommunicationShareStars.C
	CommunicationSyncNumberOfParticles.C
	CommunicationSyncNumberOfPhotons.C
	CommunicationTransferParticlesOpt.C
	CommunicationTransferPhotons.C
	CommunicationTransferStarsOpt.C
	CommunicationTransferSubgridParticles.C
	CommunicationTranspose.C
	CommunicationUpdateStarParticleCount.C
	CommunicationUtilities.C


Core methods
------------
::

	EvolveLevel.C
	EvolveHierarchy.C
	enzo.C

External boundary methods
-------------------------
::

	ExternalBoundary_AddField.C
	ExternalBoundary_AppendForcingToBaryonFields.C
	ExternalBoundary_DeleteObsoleteFields.C
	ExternalBoundary_DetachForcingFromBaryonFields.C
	ExternalBoundary_IdentifyPhysicalQuantities.C
	ExternalBoundary_InitializeExternalBoundaryFaceIO.C
	ExternalBoundary_Prepare.C
	ExternalBoundary_ReadExternalBoundary.C
	ExternalBoundary_SetDoubleMachBoundary.C
	ExternalBoundary_SetExternalBoundary.C
	ExternalBoundary_SetExternalBoundaryIO.C
	ExternalBoundary_SetExternalBoundaryParticles.C
	ExternalBoundary_SetShockPoolBoundary.C
	ExternalBoundary_SetWavePoolBoundary.C
	ExternalBoundary_SetWengenCollidingFlowBoundary.C
	ExternalBoundary_WriteExternalBoundary.C
	

Halo finder methods
-------------------
::

	FOF.C
	FOF_Finalize.C
	FOF_Initialize.C
	FOF_allocate.C
	FOF_cmpfunc.C
	FOF_density.C
	FOF_forcetree.C
	FOF_iindexx.C
	FOF_indexx.C
	FOF_ngbtree.C
	FOF_nrutil.C
	FOF_potential.C
	FOF_properties.C
	FOF_selectb.C
	FOF_sort2_flt_int.C
	FOF_sort2_int.C
	FOF_sort_int.C
	FOF_subfind.C
	FOF_subgroups.C
	FOF_unbind.C

Hydrodynamics methods
---------------------
::

	pgas2d.F
	pgas2d_dual.F
	twoshock.F
	inteuler.F
	intlgrg.F
	intpos.F
	intprim.F
	intrmp.F
	intvar.F
	calc_eigen.F
	calcdiss.F
	euler.F
	flux_hll.F
	flux_hllc.F
	flux_twoshock.F

Chemistry and energy solvers
----------------------------
::

	solve_cool.F
	solve_rate.F
	solve_rate_cool.F
	calc_photo_rates.F
	calc_rad.F
	calc_rates.F
	calc_tdust_1d.F
	calc_tdust_3d.F
	cool1d.F
	cool1d_cloudy.F
	cool1d_koyama.F
	cool1d_multi.F
	cool1d_sep.F
	cool_multi_lum.F
	cool_multi_time.F
	cool_time.F


Gravity methods
---------------
::

	mg_calc_defect.F
	mg_prolong.F
	mg_prolong2.F
	mg_relax.F
	mg_restrict.F
	FastFourierTransform.C
	FastFourierTransformPrepareComplex.C
	FastFourierTransformSGIMATH.C
	PrepareDensityField.C
	PrepareGravitatingMassField.C
	PrepareIsolatedGreensFunction.C

Hierarchy methods
-----------------
::

	RebuildHierarchy.C
	CopyZonesFromOldGrids.C
	CreateSUBlingList.C
	CreateSiblingList.C
	DepositParticleMassFlaggingField.C
	FastSiblingLocatorFinalize.C
	FastSiblingLocatorInitialize.C
	FastSiblingLocatorInitializeStaticChainingMesh.C
	FindSubgrids.C
	HilbertCurve3D.C
	LoadBalanceHilbertCurve.C
	LoadBalanceHilbertCurveRootGrids.C
	LoadBalanceSimulatedAnnealing.C
	TransposeRegionOverlap.C
	UpdateFromFinerGrids.C
	

Radiation methods
-----------------

Flux limited diffusion
^^^^^^^^^^^^^^^^^^^^^^
::

	RadiativeTransferCallFLD.C
	RHIonizationClumpInitialize.C
	RHIonizationSteepInitialize.C
	RHIonizationTestInitialize.C
	RadHydroConstTestInitialize.C
	RadHydroGreyMarshakWaveInitialize.C
	RadHydroPulseTestInitialize.C
	RadHydroRadShockInitialize.C
	RadHydroStreamTestInitialize.C
	gFLDProblem_ComputeRHS.C
	gFLDProblem_ComputeRadiationIntegrals.C
	gFLDProblem_ComputeTemperature.C
	gFLDProblem_ComputeTimeStep.C
	gFLDProblem_CrossSections.C
	gFLDProblem_Dump.C
	gFLDProblem_EnforceBoundary.C
	gFLDProblem_Evolve.C
	gFLDProblem_FInterface.C
	gFLDProblem_InitialGuess.C
	gFLDProblem_Initialize.C
	gFLDProblem_LocRHS.C
	gFLDProblem_RadiationSpectrum.C
	gFLDProblem_SetupBoundary.C
	gFLDProblem_UpdateBoundary.C
	gFLDProblem_WriteParameters.C
	gFLDProblem_constructor.C
	gFLDProblem_destructor.C
	gFLDProblem_lsetup.C
	gFLDProblem_lsolve.C
	gFLDProblem_nlresid.C
	gFLDSplit_ComputeRadiationIntegrals.C
	gFLDSplit_ComputeTemperature.C
	gFLDSplit_ComputeTimeStep.C
	gFLDSplit_CrossSections.C
	gFLDSplit_Dump.C
	gFLDSplit_EnforceBoundary.C
	gFLDSplit_Evolve.C
	gFLDSplit_FInterface.C
	gFLDSplit_InitialGuess.C
	gFLDSplit_Initialize.C
	gFLDSplit_RadiationSpectrum.C
	gFLDSplit_SetupBoundary.C
	gFLDSplit_WriteParameters.C
	gFLDSplit_constructor.C
	gFLDSplit_destructor.C

Adaptive ray tracing
^^^^^^^^^^^^^^^^^^^^
::

	EvolvePhotons.C
	RadiativeTransferComputeTimestep.C
	RadiativeTransferHealpixRoutines.C
	RadiativeTransferInitialize.C
	RadiativeTransferLoadBalanceRevert.C
	RadiativeTransferMoveLocalPhotons.C
	RadiativeTransferPrepare.C
	RadiativeTransferReadParameters.C
	RadiativeTransferWriteParameters.C
	SetSubgridMarker.C
	FindSuperSource.C
	FindSuperSourceByPosition.C
	

I/O
---
::

	OutputAsParticleData.C
	OutputCoolingTimeOnly.C
	OutputFromEvolveLevel.C
	OutputLevelInformation.C
	OutputPotentialFieldOnly.C
	OutputSmoothedDarkMatterOnly.C
	ReadAllData.C
	ReadAttr.C
	ReadDataHierarchy.C
	ReadEvolveRefineFile.C
	ReadFile.C
	ReadGridFile.C
	ReadIntFile.C
	ReadMetalCoolingRates.C
	ReadMetalCoolingRatios.C
	ReadParameterFile.C
	ReadPhotonSources.C
	ReadRadiationData.C
	ReadRadiativeTransferSpectrumTable.C
	ReadStarParticleData.C
	ReadUnits.C
	WriteAllData.C
	WriteAllDataCubes.C
	WriteConfigure.C
	WriteDataCubes.C
	WriteDataHierarchy.C
	WriteHDF5HierarchyFile.C
	WriteMemoryMap.C
	WriteParameterFile.C
	WritePhotonSources.C
	WriteRadiationData.C
	WriteStarParticleData.C
	WriteStreamData.C
	WriteStringAttr.C
	WriteTaskMap.C
	WriteTracerParticleData.C
	WriteUnits.C


Star formation methods
----------------------
::

	StarParticleAccretion.C
	StarParticleAddFeedback.C
	StarParticleCountOnly.C
	StarParticleDeath.C
	StarParticleFinalize.C
	StarParticleFindAll.C
	StarParticleInitialize.C
	StarParticleMergeMBH.C
	StarParticleMergeNew.C
	StarParticlePopIII_IMFInitialize.C
	StarParticleRadTransfer.C
	StarParticleSetRefinementLevel.C
	StarParticleSubtractAccretedMass.C
	StarRoutines.C
	Star_Accrete.C
	Star_AccreteAngularMomentum.C
	Star_ActivateNewStar.C
	Star_ApplyFeedbackTrue.C
	Star_AssignAccretedAngularMomentum.C
	Star_AssignFinalMassFromIMF.C
	Star_CalculateFeedbackParameters.C
	Star_CalculateMassAccretion.C
	Star_ComputePhotonRates.C
	Star_DeleteCopyInGridGlobal.C
	Star_DeleteParticle.C
	Star_DisableParticle.C
	Star_FindFeedbackSphere.C
	Star_HitEndpoint.C
	Star_IsARadiationSource.C
	Star_MirrorToParticle.C
	Star_MultiplyAccretionRate.C
	Star_RemoveMassFromStarAfterFeedback.C
	Star_SetFeedbackFlag.C
	Star_SphereContained.C
	Star_SubtractAccretedMassFromCell.C
	cluster_maker.F
	star_feedback_pn_snia.F
	star_maker1.F
	star_maker2.F
	star_maker3.F
	star_maker4.F
	star_maker5.F
	star_maker7.F
	star_maker8.C
	star_maker9.C
	star_maker_h2reg.F
	sink_maker.C
	pop3_color_maker.F
	pop3_maker.F
	pop3_properties.F

Utilities
---------
::

	cic_deposit.F
	cic_flag.F
	cic_interp.F
	cicinterp.F
	smooth.F
	smooth_deposit.F
	rotate2d.F
	rotate3d.F
	int_lin3d.F
	int_spline.F
	interp1d.F
	interp2d.F
	interp3d.F
	interpolate.F
	utilities.F
	MemoryAllocationRoutines.C
	MemoryPoolRoutines.C
	SortCompareFunctions.C


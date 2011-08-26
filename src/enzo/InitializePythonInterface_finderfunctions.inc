#ifdef MEMORY_POOL
#endif /* MEMORY_POOL */ 
#ifdef DEFINE_STORAGE
#endif /* DEFINE_STORAGE */ 
#ifdef NEW_PROBLEM_TYPES
#endif /* NEW_PROBLEM_TYPES */ 

    if (strncmp(parameter_name, "LoadBalancing", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadBalancing);
        return rv;
    }

    if (strncmp(parameter_name, "LoadBalancingCycleSkip", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadBalancingCycleSkip);
        return rv;
    }

    if (strncmp(parameter_name, "ResetLoadBalancing", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ResetLoadBalancing);
        return rv;
    }

    if (strncmp(parameter_name, "CoresPerNode", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) CoresPerNode);
        return rv;
    }

    if (strncmp(parameter_name, "PreviousMaxTask", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) PreviousMaxTask);
        return rv;
    }

    if (strncmp(parameter_name, "LoadBalancingMinLevel", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadBalancingMinLevel);
        return rv;
    }

    if (strncmp(parameter_name, "LoadBalancingMaxLevel", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadBalancingMaxLevel);
        return rv;
    }

    if (strncmp(parameter_name, "FileDirectedOutput", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) FileDirectedOutput);
        return rv;
    }

    if (strncmp(parameter_name, "WriteBinaryHierarchy", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) WriteBinaryHierarchy);
        return rv;
    }

    if (strncmp(parameter_name, "debug", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) debug);
        return rv;
    }

    if (strncmp(parameter_name, "debug1", 7) == 0) {
        PyObject *rv = PyLong_FromLong((long) debug1);
        return rv;
    }

    if (strncmp(parameter_name, "debug2", 7) == 0) {
        PyObject *rv = PyLong_FromLong((long) debug2);
        return rv;
    }

    if (strncmp(parameter_name, "extract", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) extract);
        return rv;
    }

    if (strncmp(parameter_name, "CheckpointRestart", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) CheckpointRestart);
        return rv;
    }

    if (strncmp(parameter_name, "ProblemType", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) ProblemType);
        return rv;
    }
#ifdef NEW_PROBLEM_TYPES
#endif /* NEW_PROBLEM_TYPES */ 

    if (strncmp(parameter_name, "huge_number", 12) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) huge_number);
        return rv;
    }

    if (strncmp(parameter_name, "tiny_number", 12) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) tiny_number);
        return rv;
    }

    if (strncmp(parameter_name, "Gamma", 6) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) Gamma);
        return rv;
    }

    if (strncmp(parameter_name, "PressureFree", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) PressureFree);
        return rv;
    }

    if (strncmp(parameter_name, "RefineBy", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) RefineBy);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumRefinementLevel", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) MaximumRefinementLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumGravityRefinementLevel", 30) == 0) {
        PyObject *rv = PyLong_FromLong((long) MaximumGravityRefinementLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumParticleRefinementLevel", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) MaximumParticleRefinementLevel);
        return rv;
    }

    if (strncmp(parameter_name, "FastSiblingLocatorEntireDomain", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) FastSiblingLocatorEntireDomain);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineRegionMinRefinementLevel", 35) == 0) {
        PyObject *rv = PyLong_FromLong((long) MustRefineRegionMinRefinementLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MetallicityRefinementMinLevel", 30) == 0) {
        PyObject *rv = PyLong_FromLong((long) MetallicityRefinementMinLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MetallicityRefinementMinMetallicity", 36) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MetallicityRefinementMinMetallicity);
        return rv;
    }

    if (strncmp(parameter_name, "MetallicityRefinementMinDensity", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MetallicityRefinementMinDensity);
        return rv;
    }

    if (strncmp(parameter_name, "TimestepSafetyVelocity", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) TimestepSafetyVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "FluxCorrection", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) FluxCorrection);
        return rv;
    }

    if (strncmp(parameter_name, "ConservativeInterpolation", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) ConservativeInterpolation);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumEfficiency", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumEfficiency);
        return rv;
    }

    if (strncmp(parameter_name, "SubgridSizeAutoAdjust", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) SubgridSizeAutoAdjust);
        return rv;
    }

    if (strncmp(parameter_name, "OptimalSubgridsPerProcessor", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) OptimalSubgridsPerProcessor);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumSubgridEdge", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) MinimumSubgridEdge);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumSubgridSize", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) MaximumSubgridSize);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfBufferZones", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfBufferZones);
        return rv;
    }

    if (strncmp(parameter_name, "DimUnits", 9) == 0) {
        PyObject *rv = PyString_FromString((char *) DimUnits);
        return rv;
    }

    if (strncmp(parameter_name, "DimLabels", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) DimLabels);
        return rv;
    }

    if (strncmp(parameter_name, "DataLabel", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) DataLabel);
        return rv;
    }

    if (strncmp(parameter_name, "DataUnits", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) DataUnits);
        return rv;
    }

    if (strncmp(parameter_name, "RefineRegionAutoAdjust", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) RefineRegionAutoAdjust);
        return rv;
    }

    if (strncmp(parameter_name, "UniformGravity", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) UniformGravity);
        return rv;
    }

    if (strncmp(parameter_name, "UniformGravityDirection", 24) == 0) {
        PyObject *rv = PyLong_FromLong((long) UniformGravityDirection);
        return rv;
    }

    if (strncmp(parameter_name, "UniformGravityConstant", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) UniformGravityConstant);
        return rv;
    }

    if (strncmp(parameter_name, "PointSourceGravity", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) PointSourceGravity);
        return rv;
    }

    if (strncmp(parameter_name, "PointSourceGravityConstant", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PointSourceGravityConstant);
        return rv;
    }

    if (strncmp(parameter_name, "PointSourceGravityCoreRadius", 29) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PointSourceGravityCoreRadius);
        return rv;
    }

    if (strncmp(parameter_name, "SelfGravity", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) SelfGravity);
        return rv;
    }

    if (strncmp(parameter_name, "SelfGravityGasOff", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) SelfGravityGasOff);
        return rv;
    }

    if (strncmp(parameter_name, "AccretionKernal", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) AccretionKernal);
        return rv;
    }

    if (strncmp(parameter_name, "CopyGravPotential", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) CopyGravPotential);
        return rv;
    }

    if (strncmp(parameter_name, "PotentialIterations", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) PotentialIterations);
        return rv;
    }

    if (strncmp(parameter_name, "BaryonSelfGravityApproximation", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) BaryonSelfGravityApproximation);
        return rv;
    }

    if (strncmp(parameter_name, "GravitationalConstant", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GravitationalConstant);
        return rv;
    }

    if (strncmp(parameter_name, "S2ParticleSize", 15) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) S2ParticleSize);
        return rv;
    }

    if (strncmp(parameter_name, "GravityResolution", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GravityResolution);
        return rv;
    }

    if (strncmp(parameter_name, "ComputePotential", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) ComputePotential);
        return rv;
    }

    if (strncmp(parameter_name, "WritePotential", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) WritePotential);
        return rv;
    }

    if (strncmp(parameter_name, "GreensFunctionMaxNumber", 24) == 0) {
        PyObject *rv = PyLong_FromLong((long) GreensFunctionMaxNumber);
        return rv;
    }

    if (strncmp(parameter_name, "GreensFunctionMaxSize", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) GreensFunctionMaxSize);
        return rv;
    }

    if (strncmp(parameter_name, "DualEnergyFormalism", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) DualEnergyFormalism);
        return rv;
    }

    if (strncmp(parameter_name, "DualEnergyFormalismEta1", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DualEnergyFormalismEta1);
        return rv;
    }

    if (strncmp(parameter_name, "DualEnergyFormalismEta2", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DualEnergyFormalismEta2);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleCourantSafetyNumber", 28) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ParticleCourantSafetyNumber);
        return rv;
    }

    if (strncmp(parameter_name, "RootGridCourantSafetyNumber", 28) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RootGridCourantSafetyNumber);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeCooling", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeCooling);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeCoolingModel", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeCoolingModel);
        return rv;
    }

    if (strncmp(parameter_name, "GadgetEquilibriumCooling", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) GadgetEquilibriumCooling);
        return rv;
    }

    if (strncmp(parameter_name, "RandomForcing", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) RandomForcing);
        return rv;
    }

    if (strncmp(parameter_name, "RandomForcingEdot", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RandomForcingEdot);
        return rv;
    }

    if (strncmp(parameter_name, "RandomForcingMachNumber", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RandomForcingMachNumber);
        return rv;
    }

    if (strncmp(parameter_name, "MultiSpecies", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) MultiSpecies);
        return rv;
    }

    if (strncmp(parameter_name, "PrimordialChemistrySolver", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) PrimordialChemistrySolver);
        return rv;
    }

    if (strncmp(parameter_name, "ThreeBodyRate", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) ThreeBodyRate);
        return rv;
    }

    if (strncmp(parameter_name, "GloverChemistryModel", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) GloverChemistryModel);
        return rv;
    }

    if (strncmp(parameter_name, "GloverRadiationBackground", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) GloverRadiationBackground);
        return rv;
    }

    if (strncmp(parameter_name, "GloverOpticalDepth", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) GloverOpticalDepth);
        return rv;
    }

    if (strncmp(parameter_name, "MultiMetals", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) MultiMetals);
        return rv;
    }

    if (strncmp(parameter_name, "ShockMethod", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShockMethod);
        return rv;
    }

    if (strncmp(parameter_name, "ShockTemperatureFloor", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ShockTemperatureFloor);
        return rv;
    }

    if (strncmp(parameter_name, "StorePreShockFields", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) StorePreShockFields);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationFieldType", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationFieldType);
        return rv;
    }

    if (strncmp(parameter_name, "AdjustUVBackground", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) AdjustUVBackground);
        return rv;
    }

    if (strncmp(parameter_name, "SetUVBAmplitude", 16) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SetUVBAmplitude);
        return rv;
    }

    if (strncmp(parameter_name, "SetHeIIHeatingScale", 20) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SetHeIIHeatingScale);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationFieldLevelRecompute", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationFieldLevelRecompute);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationXRaySecondaryIon", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationXRaySecondaryIon);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationXRayComptonHeating", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationXRayComptonHeating);
        return rv;
    }

    if (strncmp(parameter_name, "TabulatedLWBackground", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) TabulatedLWBackground);
        return rv;
    }

    if (strncmp(parameter_name, "PhotoelectricHeating", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) PhotoelectricHeating);
        return rv;
    }

    if (strncmp(parameter_name, "OutputCoolingTime", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputCoolingTime);
        return rv;
    }

    if (strncmp(parameter_name, "OutputTemperature", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputTemperature);
        return rv;
    }

    if (strncmp(parameter_name, "OutputSmoothedDarkMatter", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputSmoothedDarkMatter);
        return rv;
    }

    if (strncmp(parameter_name, "SmoothedDarkMatterNeighbors", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmoothedDarkMatterNeighbors);
        return rv;
    }

    if (strncmp(parameter_name, "OutputGriddedStarParticle", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputGriddedStarParticle);
        return rv;
    }

    if (strncmp(parameter_name, "ZEUSLinearArtificialViscosity", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ZEUSLinearArtificialViscosity);
        return rv;
    }

    if (strncmp(parameter_name, "ZEUSQuadraticArtificialViscosity", 33) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ZEUSQuadraticArtificialViscosity);
        return rv;
    }

    if (strncmp(parameter_name, "UseMinimumPressureSupport", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseMinimumPressureSupport);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumPressureSupportParameter", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumPressureSupportParameter);
        return rv;
    }

    if (strncmp(parameter_name, "RefineRegionTimeType", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) RefineRegionTimeType);
        return rv;
    }

    if (strncmp(parameter_name, "MyProcessorNumber", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) MyProcessorNumber);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfProcessors", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfProcessors);
        return rv;
    }

    if (strncmp(parameter_name, "CommunicationTime", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CommunicationTime);
        return rv;
    }

    if (strncmp(parameter_name, "ParallelRootGridIO", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParallelRootGridIO);
        return rv;
    }

    if (strncmp(parameter_name, "ParallelParticleIO", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParallelParticleIO);
        return rv;
    }

    if (strncmp(parameter_name, "Unigrid", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) Unigrid);
        return rv;
    }

    if (strncmp(parameter_name, "CubeDumpEnabled", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) CubeDumpEnabled);
        return rv;
    }

    if (strncmp(parameter_name, "PartitionNestedGrids", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) PartitionNestedGrids);
        return rv;
    }

    if (strncmp(parameter_name, "StaticPartitionNestedGrids", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) StaticPartitionNestedGrids);
        return rv;
    }

    if (strncmp(parameter_name, "ExtractFieldsOnly", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExtractFieldsOnly);
        return rv;
    }

    if (strncmp(parameter_name, "First_Pass", 11) == 0) {
        PyObject *rv = PyLong_FromLong((long) First_Pass);
        return rv;
    }

    if (strncmp(parameter_name, "UnigridTranspose", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) UnigridTranspose);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfRootGridTilesPerDimensionPerProcessor", 46) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfRootGridTilesPerDimensionPerProcessor);
        return rv;
    }

    if (strncmp(parameter_name, "CosmologySimulationNumberOfInitialGrids", 40) == 0) {
        PyObject *rv = PyLong_FromLong((long) CosmologySimulationNumberOfInitialGrids);
        return rv;
    }

    if (strncmp(parameter_name, "OutputOnDensity", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputOnDensity);
        return rv;
    }

    if (strncmp(parameter_name, "StartDensityOutputs", 20) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StartDensityOutputs);
        return rv;
    }

    if (strncmp(parameter_name, "CurrentDensityOutput", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CurrentDensityOutput);
        return rv;
    }

    if (strncmp(parameter_name, "CurrentMaximumDensity", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CurrentMaximumDensity);
        return rv;
    }

    if (strncmp(parameter_name, "IncrementDensityOutput", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) IncrementDensityOutput);
        return rv;
    }

    if (strncmp(parameter_name, "PythonTopGridSkip", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) PythonTopGridSkip);
        return rv;
    }

    if (strncmp(parameter_name, "PythonSubcycleSkip", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) PythonSubcycleSkip);
        return rv;
    }

    if (strncmp(parameter_name, "InlineHaloFinder", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) InlineHaloFinder);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderSubfind", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderSubfind);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderOutputParticleList", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderOutputParticleList);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderMinimumSize", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderMinimumSize);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderCycleSkip", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderCycleSkip);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderRunAfterOutput", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderRunAfterOutput);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderLinkingLength", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloFinderLinkingLength);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderTimestep", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloFinderTimestep);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderLastTime", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloFinderLastTime);
        return rv;
    }

    if (strncmp(parameter_name, "DepositPositionsParticleSmoothRadius", 37) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DepositPositionsParticleSmoothRadius);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumPressureJumpForRefinement", 33) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumPressureJumpForRefinement);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumEnergyRatioForRefinement", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumEnergyRatioForRefinement);
        return rv;
    }

    if (strncmp(parameter_name, "RefineByJeansLengthSafetyFactor", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RefineByJeansLengthSafetyFactor);
        return rv;
    }

    if (strncmp(parameter_name, "JeansRefinementColdTemperature", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) JeansRefinementColdTemperature);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineParticlesRefineToLevel", 33) == 0) {
        PyObject *rv = PyLong_FromLong((long) MustRefineParticlesRefineToLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineParticlesRefineToLevelAutoAdjust", 43) == 0) {
        PyObject *rv = PyLong_FromLong((long) MustRefineParticlesRefineToLevelAutoAdjust);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineParticlesMinimumMass", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MustRefineParticlesMinimumMass);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumShearForRefinement", 26) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumShearForRefinement);
        return rv;
    }

    if (strncmp(parameter_name, "RefineByResistiveLengthSafetyFactor", 36) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RefineByResistiveLengthSafetyFactor);
        return rv;
    }

    if (strncmp(parameter_name, "ShockwaveRefinementMinMach", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ShockwaveRefinementMinMach);
        return rv;
    }

    if (strncmp(parameter_name, "ShockwaveRefinementMinVelocity", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ShockwaveRefinementMinVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "ShockwaveRefinementMaxLevel", 28) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ShockwaveRefinementMaxLevel);
        return rv;
    }

    if (strncmp(parameter_name, "NohProblemFullBox", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) NohProblemFullBox);
        return rv;
    }

    if (strncmp(parameter_name, "ComovingCoordinates", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) ComovingCoordinates);
        return rv;
    }

    if (strncmp(parameter_name, "StarParticleCreation", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) StarParticleCreation);
        return rv;
    }

    if (strncmp(parameter_name, "StarParticleFeedback", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) StarParticleFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfParticleAttributes", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfParticleAttributes);
        return rv;
    }

    if (strncmp(parameter_name, "AddParticleAttributes", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) AddParticleAttributes);
        return rv;
    }

    if (strncmp(parameter_name, "BigStarFormation", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) BigStarFormation);
        return rv;
    }

    if (strncmp(parameter_name, "BigStarFormationDone", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) BigStarFormationDone);
        return rv;
    }

    if (strncmp(parameter_name, "BigStarSeparation", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) BigStarSeparation);
        return rv;
    }

    if (strncmp(parameter_name, "SimpleQ", 8) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SimpleQ);
        return rv;
    }

    if (strncmp(parameter_name, "SimpleRampTime", 15) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SimpleRampTime);
        return rv;
    }

    if (strncmp(parameter_name, "CubeDumps", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) CubeDumps);
        return rv;
    }

    if (strncmp(parameter_name, "TracerParticleOn", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) TracerParticleOn);
        return rv;
    }

    if (strncmp(parameter_name, "TracerParticleCreationSpacing", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) TracerParticleCreationSpacing);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleTypeInFile", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParticleTypeInFile);
        return rv;
    }

    if (strncmp(parameter_name, "OutputParticleTypeGrouping", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputParticleTypeGrouping);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalBoundaryIO", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalBoundaryIO);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalBoundaryTypeIO", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalBoundaryTypeIO);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalBoundaryValueIO", 24) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalBoundaryValueIO);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalBoundaryField", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalBoundaryField);
        return rv;
    }

    if (strncmp(parameter_name, "SimpleConstantBoundary", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) SimpleConstantBoundary);
        return rv;
    }

    if (strncmp(parameter_name, "LoadGridDataAtStart", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadGridDataAtStart);
        return rv;
    }

    if (strncmp(parameter_name, "PrevParameterFileName", 22) == 0) {
        PyObject *rv = PyString_FromString((char *) PrevParameterFileName);
        return rv;
    }
#ifdef MPI_INSTRUMENTATION

    if (strncmp(parameter_name, "tracename", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) tracename);
        return rv;
    }

    if (strncmp(parameter_name, "flagging_count", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) flagging_count);
        return rv;
    }

    if (strncmp(parameter_name, "in_count", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) in_count);
        return rv;
    }

    if (strncmp(parameter_name, "out_count", 10) == 0) {
        PyObject *rv = PyLong_FromLong((long) out_count);
        return rv;
    }

    if (strncmp(parameter_name, "moving_count", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) moving_count);
        return rv;
    }

    if (strncmp(parameter_name, "flagging_pct", 13) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) flagging_pct);
        return rv;
    }

    if (strncmp(parameter_name, "moving_pct", 11) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) moving_pct);
        return rv;
    }
#endif /* MPI_INSTRUMENTATION */ 

    if (strncmp(parameter_name, "name", 5) == 0) {
        PyObject *rv = PyString_FromString((char *) name);
        return rv;
    }

    if (strncmp(parameter_name, "traceMPI", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) traceMPI);
        return rv;
    }
#ifdef MEM_TRACE

    if (strncmp(parameter_name, "traceMEM", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) traceMEM);
        return rv;
    }

    if (strncmp(parameter_name, "memtracename", 13) == 0) {
        PyObject *rv = PyString_FromString((char *) memtracename);
        return rv;
    }
#endif /* MEM_TRACE */ 

    if (strncmp(parameter_name, "MovieSkipTimestep", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) MovieSkipTimestep);
        return rv;
    }

    if (strncmp(parameter_name, "Movie3DVolumes", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) Movie3DVolumes);
        return rv;
    }

    if (strncmp(parameter_name, "MovieVertexCentered", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) MovieVertexCentered);
        return rv;
    }

    if (strncmp(parameter_name, "NewMovieDumpNumber", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) NewMovieDumpNumber);
        return rv;
    }

    if (strncmp(parameter_name, "NewMovieParticleOn", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) NewMovieParticleOn);
        return rv;
    }

    if (strncmp(parameter_name, "UseHydro", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseHydro);
        return rv;
    }

    if (strncmp(parameter_name, "Coordinate", 11) == 0) {
        PyObject *rv = PyLong_FromLong((long) Coordinate);
        return rv;
    }

    if (strncmp(parameter_name, "NSpecies", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) NSpecies);
        return rv;
    }

    if (strncmp(parameter_name, "NColor", 7) == 0) {
        PyObject *rv = PyLong_FromLong((long) NColor);
        return rv;
    }

    if (strncmp(parameter_name, "Theta_Limiter", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) Theta_Limiter);
        return rv;
    }

    if (strncmp(parameter_name, "RKOrder", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) RKOrder);
        return rv;
    }

    if (strncmp(parameter_name, "UsePhysicalUnit", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) UsePhysicalUnit);
        return rv;
    }

    if (strncmp(parameter_name, "iden", 5) == 0) {
        PyObject *rv = PyLong_FromLong((long) iden);
        return rv;
    }

    if (strncmp(parameter_name, "ietot", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) ietot);
        return rv;
    }

    if (strncmp(parameter_name, "ivx", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) ivx);
        return rv;
    }

    if (strncmp(parameter_name, "ivy", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) ivy);
        return rv;
    }

    if (strncmp(parameter_name, "ivz", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) ivz);
        return rv;
    }

    if (strncmp(parameter_name, "iBx", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iBx);
        return rv;
    }

    if (strncmp(parameter_name, "iBy", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iBy);
        return rv;
    }

    if (strncmp(parameter_name, "iBz", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iBz);
        return rv;
    }

    if (strncmp(parameter_name, "iPhi", 5) == 0) {
        PyObject *rv = PyLong_FromLong((long) iPhi);
        return rv;
    }

    if (strncmp(parameter_name, "ieint", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) ieint);
        return rv;
    }

    if (strncmp(parameter_name, "iD", 3) == 0) {
        PyObject *rv = PyLong_FromLong((long) iD);
        return rv;
    }

    if (strncmp(parameter_name, "iEtot", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) iEtot);
        return rv;
    }

    if (strncmp(parameter_name, "iS1", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iS1);
        return rv;
    }

    if (strncmp(parameter_name, "iS2", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iS2);
        return rv;
    }

    if (strncmp(parameter_name, "iS3", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iS3);
        return rv;
    }

    if (strncmp(parameter_name, "iEint", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) iEint);
        return rv;
    }

    if (strncmp(parameter_name, "SmallRho", 9) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmallRho);
        return rv;
    }

    if (strncmp(parameter_name, "SmallP", 7) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmallP);
        return rv;
    }

    if (strncmp(parameter_name, "SmallEint", 10) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmallEint);
        return rv;
    }

    if (strncmp(parameter_name, "SmallT", 7) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmallT);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumAlvenSpeed", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MaximumAlvenSpeed);
        return rv;
    }

    if (strncmp(parameter_name, "NEQ_HYDRO", 10) == 0) {
        PyObject *rv = PyLong_FromLong((long) NEQ_HYDRO);
        return rv;
    }

    if (strncmp(parameter_name, "NEQ_MHD", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) NEQ_MHD);
        return rv;
    }

    if (strncmp(parameter_name, "ReconstructionMethod", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) ReconstructionMethod);
        return rv;
    }

    if (strncmp(parameter_name, "PositiveReconstruction", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) PositiveReconstruction);
        return rv;
    }

    if (strncmp(parameter_name, "RiemannSolverFallback", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) RiemannSolverFallback);
        return rv;
    }

    if (strncmp(parameter_name, "RiemannSolver", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) RiemannSolver);
        return rv;
    }

    if (strncmp(parameter_name, "ConservativeReconstruction", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) ConservativeReconstruction);
        return rv;
    }

    if (strncmp(parameter_name, "EOSType", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) EOSType);
        return rv;
    }

    if (strncmp(parameter_name, "EOSSoundSpeed", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) EOSSoundSpeed);
        return rv;
    }

    if (strncmp(parameter_name, "EOSCriticalDensity", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) EOSCriticalDensity);
        return rv;
    }

    if (strncmp(parameter_name, "EOSGamma", 9) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) EOSGamma);
        return rv;
    }

    if (strncmp(parameter_name, "C_h", 4) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) C_h);
        return rv;
    }

    if (strncmp(parameter_name, "C_p", 4) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) C_p);
        return rv;
    }

    if (strncmp(parameter_name, "DivBDampingLength", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DivBDampingLength);
        return rv;
    }

    if (strncmp(parameter_name, "UseConstantAcceleration", 24) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseConstantAcceleration);
        return rv;
    }

    if (strncmp(parameter_name, "Mu", 3) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) Mu);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalGravity", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalGravity);
        return rv;
    }

    if (strncmp(parameter_name, "StringKick", 11) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StringKick);
        return rv;
    }

    if (strncmp(parameter_name, "StringKickDimension", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) StringKickDimension);
        return rv;
    }

    if (strncmp(parameter_name, "UseFloor", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseFloor);
        return rv;
    }

    if (strncmp(parameter_name, "UseViscosity", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseViscosity);
        return rv;
    }

    if (strncmp(parameter_name, "ViscosityCoefficient", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ViscosityCoefficient);
        return rv;
    }

    if (strncmp(parameter_name, "UseAmbipolarDiffusion", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseAmbipolarDiffusion);
        return rv;
    }

    if (strncmp(parameter_name, "UseResistivity", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseResistivity);
        return rv;
    }

    if (strncmp(parameter_name, "UseH2OnDust", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseH2OnDust);
        return rv;
    }

    if (strncmp(parameter_name, "CoolingCutOffDensity1", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CoolingCutOffDensity1);
        return rv;
    }

    if (strncmp(parameter_name, "CoolingCutOffDensity2", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CoolingCutOffDensity2);
        return rv;
    }

    if (strncmp(parameter_name, "CoolingPowerCutOffDensity1", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CoolingPowerCutOffDensity1);
        return rv;
    }

    if (strncmp(parameter_name, "CoolingPowerCutOffDensity2", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CoolingPowerCutOffDensity2);
        return rv;
    }

    if (strncmp(parameter_name, "CoolingCutOffTemperature", 25) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CoolingCutOffTemperature);
        return rv;
    }

    if (strncmp(parameter_name, "CoolingModel", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) CoolingModel);
        return rv;
    }

    if (strncmp(parameter_name, "HaloConcentration", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloConcentration);
        return rv;
    }

    if (strncmp(parameter_name, "HaloRedshift", 13) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloRedshift);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalGravityDensity", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ExternalGravityDensity);
        return rv;
    }

    if (strncmp(parameter_name, "UseDivergenceCleaning", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseDivergenceCleaning);
        return rv;
    }

    if (strncmp(parameter_name, "DivergenceCleaningBoundaryBuffer", 33) == 0) {
        PyObject *rv = PyLong_FromLong((long) DivergenceCleaningBoundaryBuffer);
        return rv;
    }

    if (strncmp(parameter_name, "DivergenceCleaningThreshold", 28) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DivergenceCleaningThreshold);
        return rv;
    }

    if (strncmp(parameter_name, "PoissonApproximationThreshold", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PoissonApproximationThreshold);
        return rv;
    }

    if (strncmp(parameter_name, "PoissonBoundaryType", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) PoissonBoundaryType);
        return rv;
    }

    if (strncmp(parameter_name, "ShiningParticleID", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShiningParticleID);
        return rv;
    }

    if (strncmp(parameter_name, "SinkMergeDistance", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SinkMergeDistance);
        return rv;
    }

    if (strncmp(parameter_name, "SinkMergeMass", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SinkMergeMass);
        return rv;
    }

    if (strncmp(parameter_name, "TotalSinkMass", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) TotalSinkMass);
        return rv;
    }

    if (strncmp(parameter_name, "StellarWindFeedback", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) StellarWindFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "StellarWindTurnOnMass", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StellarWindTurnOnMass);
        return rv;
    }

    if (strncmp(parameter_name, "MSStellarWindTurnOnMass", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MSStellarWindTurnOnMass);
        return rv;
    }

    if (strncmp(parameter_name, "NBodyDirectSummation", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) NBodyDirectSummation);
        return rv;
    }

    if (strncmp(parameter_name, "UseDrivingField", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseDrivingField);
        return rv;
    }

    if (strncmp(parameter_name, "DrivingEfficiency", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DrivingEfficiency);
        return rv;
    }

    if (strncmp(parameter_name, "UseCUDA", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseCUDA);
        return rv;
    }

    if (strncmp(parameter_name, "ran1_init", 10) == 0) {
        PyObject *rv = PyLong_FromLong((long) ran1_init);
        return rv;
    }

    if (strncmp(parameter_name, "rand_init", 10) == 0) {
        PyObject *rv = PyLong_FromLong((long) rand_init);
        return rv;
    }
#ifdef STAGE_INPUT

    if (strncmp(parameter_name, "StageInput", 11) == 0) {
        PyObject *rv = PyLong_FromLong((long) StageInput);
        return rv;
    }

    if (strncmp(parameter_name, "LocalPath", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) LocalPath);
        return rv;
    }

    if (strncmp(parameter_name, "GlobalPath", 11) == 0) {
        PyObject *rv = PyString_FromString((char *) GlobalPath);
        return rv;
    }
#endif /* STAGE_INPUT */ 
#ifdef USE_PYTHON

    if (strncmp(parameter_name, "NumberOfPythonCalls", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfPythonCalls);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfPythonTopGridCalls", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfPythonTopGridCalls);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfPythonSubcycleCalls", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfPythonSubcycleCalls);
        return rv;
    }
#endif /* USE_PYTHON */ 

    if (strncmp(parameter_name, "MetalCooling", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) MetalCooling);
        return rv;
    }

    if (strncmp(parameter_name, "CIECooling", 11) == 0) {
        PyObject *rv = PyLong_FromLong((long) CIECooling);
        return rv;
    }

    if (strncmp(parameter_name, "H2OpticalDepthApproximation", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) H2OpticalDepthApproximation);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransfer", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransfer);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferHydrogenOnly", 30) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferHydrogenOnly);
        return rv;
    }
#ifdef TRANSFER

    if (strncmp(parameter_name, "PhotonTime", 11) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PhotonTime);
        return rv;
    }

    if (strncmp(parameter_name, "dtPhoton", 9) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) dtPhoton);
        return rv;
    }
#endif /* TRANSFER */ 
#ifdef DEFINE_STORAGE
#endif /* DEFINE_STORAGE */ 

    if (strncmp(parameter_name, "RadiativeTransferSourceRadius", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferSourceRadius);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPropagationSpeedFraction", 42) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferPropagationSpeedFraction);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPropagationDistance", 37) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferPropagationDistance);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferSplitPhotonRadius", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferSplitPhotonRadius);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferRaysPerCell", 29) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferRaysPerCell);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferInitialHEALPixLevel", 37) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferInitialHEALPixLevel);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPhotonEscapeRadius", 36) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferPhotonEscapeRadius);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferInterpolateField", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferInterpolateField);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferTimestepVelocityLimit", 39) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferTimestepVelocityLimit);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferSourceClustering", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferSourceClustering);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPhotonMergeRadius", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferPhotonMergeRadius);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationPressure", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationPressure);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationPressureScale", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiationPressureScale);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferOpticallyThinH2", 33) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferOpticallyThinH2);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPeriodicBoundary", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferPeriodicBoundary);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferFLDCallOnLevel", 32) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferFLDCallOnLevel);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferHIIRestrictedTimestep", 39) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferHIIRestrictedTimestep);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferAdaptiveTimestep", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferAdaptiveTimestep);
        return rv;
    }

    if (strncmp(parameter_name, "GlobalMaximumkphIfront", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GlobalMaximumkphIfront);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferSourceBeamAngle", 33) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferSourceBeamAngle);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferTraceSpectrum", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferTraceSpectrum);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferLoadBalance", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferLoadBalance);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferFluxBackgroundLimit", 37) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferFluxBackgroundLimit);
        return rv;
    }
#ifdef MEMORY_POOL
#endif /* MEMORY_POOL */ 

    if (strncmp(parameter_name, "RadiativeTransferCoupledRateSolver", 35) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferCoupledRateSolver);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferFLD", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferFLD);
        return rv;
    }

    if (strncmp(parameter_name, "ImplicitProblem", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) ImplicitProblem);
        return rv;
    }

    if (strncmp(parameter_name, "StarMakerEmissivityField", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) StarMakerEmissivityField);
        return rv;
    }

    if (strncmp(parameter_name, "uv_param", 9) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) uv_param);
        return rv;
    }

    if (strncmp(parameter_name, "AngularVelocity", 16) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) AngularVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "VelocityGradient", 17) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) VelocityGradient);
        return rv;
    }

    if (strncmp(parameter_name, "ShearingBoundaryDirection", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShearingBoundaryDirection);
        return rv;
    }

    if (strncmp(parameter_name, "ShearingVelocityDirection", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShearingVelocityDirection);
        return rv;
    }

    if (strncmp(parameter_name, "ShearingOtherDirection", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShearingOtherDirection);
        return rv;
    }

    if (strncmp(parameter_name, "useMHD", 7) == 0) {
        PyObject *rv = PyLong_FromLong((long) useMHD);
        return rv;
    }

    if (strncmp(parameter_name, "ShearingBoxProblemType", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShearingBoxProblemType);
        return rv;
    }

    if (strncmp(parameter_name, "IsothermalSoundSpeed", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) IsothermalSoundSpeed);
        return rv;
    }

    if (strncmp(parameter_name, "RefineByJeansLengthUnits", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) RefineByJeansLengthUnits);
        return rv;
    }

    if (strncmp(parameter_name, "MoveParticlesBetweenSiblings", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) MoveParticlesBetweenSiblings);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleSplitterIterations", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParticleSplitterIterations);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleSplitterChildrenParticleSeparation", 43) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ParticleSplitterChildrenParticleSeparation);
        return rv;
    }

    if (strncmp(parameter_name, "ResetMagneticField", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ResetMagneticField);
        return rv;
    }

    if (strncmp(parameter_name, "MBHParticleIO", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) MBHParticleIO);
        return rv;
    }

    if (strncmp(parameter_name, "OutputWhenJetsHaveNotEjected", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputWhenJetsHaveNotEjected);
        return rv;
    }

    if (strncmp(parameter_name, "VelAnyl", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) VelAnyl);
        return rv;
    }

    if (strncmp(parameter_name, "BAnyl", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) BAnyl);
        return rv;
    }

    if (strncmp(parameter_name, "current_error", 14) == 0) {
        PyObject *rv = PyString_FromString((char *) current_error);
        return rv;
    }

    if (strncmp(parameter_name, "IsotropicConduction", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) IsotropicConduction);
        return rv;
    }

    if (strncmp(parameter_name, "AnisotropicConduction", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) AnisotropicConduction);
        return rv;
    }

    if (strncmp(parameter_name, "IsotropicConductionSpitzerFraction", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) IsotropicConductionSpitzerFraction);
        return rv;
    }

    if (strncmp(parameter_name, "AnisotropicConductionSpitzerFraction", 37) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) AnisotropicConductionSpitzerFraction);
        return rv;
    }

    if (strncmp(parameter_name, "ConductionCourantSafetyNumber", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ConductionCourantSafetyNumber);
        return rv;
    }

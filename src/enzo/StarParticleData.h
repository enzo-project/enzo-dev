/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE STAR PARTICLES
/
/  written by: Greg Bryan
/  date:       February, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#ifndef STAR_PARTICLE_DATA_DEFINED__
#define STAR_PARTICLE_DATA_DEFINED__

#include "global_data.h"
#ifdef DEFINE_STORAGE
# define SPEXTERN
#else /* DEFINE_STORAGE */
# define SPEXTERN extern
#endif /* DEFINE_STORAGE */

#define STAR_PARTICLE_NUMBER_START 1000000000

/* #include "macros_and_parameters.h" */

struct ParticleEntry {
  FLOAT Position[3];
  float Mass;
  float Velocity[3];
  float Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  PINT Number;
  int Type;
  float InitialMass = 0;
};


/* Number of Star particles. */

SPEXTERN int NumberOfStarParticles;
SPEXTERN int NumberOfActiveParticles;
SPEXTERN int NumberOfDeletedParticles;
SPEXTERN PINT NumberOfOtherParticles; //all the particles other than type=2
SPEXTERN PINT NextActiveParticleID;
SPEXTERN int G_TotalNumberOfStars;

/* Star particle parameters. */

SPEXTERN int StarFeedbackType;
SPEXTERN int StarMakerTypeIaSNe;
SPEXTERN int StarMakerTypeIISNeMetalField;
SPEXTERN int StarMakerPlanetaryNebulae;
SPEXTERN int StarMakerTimeIndependentFormation;
SPEXTERN float StarMakerOverDensityThreshold;
SPEXTERN int StarMakerUseOverDensityThreshold;
SPEXTERN float StarMakerMaximumFractionCell;
SPEXTERN float StarMakerSHDensityThreshold;
SPEXTERN float StarMakerMassEfficiency;
SPEXTERN float StarMakerMinimumMass;
SPEXTERN float StarMakerMinimumDynamicalTime;
SPEXTERN float StarMassEjectionFraction;
SPEXTERN float StarMetalYield;
SPEXTERN float StarEnergyToThermalFeedback;
SPEXTERN float StarFeedbackAdditionalThermalEnergy;
SPEXTERN float MomentumMultiplier;
SPEXTERN int MomentumCancellationToThermal;
SPEXTERN int WriteFeedbackLogFiles;
SPEXTERN float StarEnergyFeedbackRate;
SPEXTERN float StarEnergyToStellarUV;
SPEXTERN float StarEnergyToQuasarUV;
SPEXTERN int StarFeedbackDistRadius;
SPEXTERN int StarFeedbackDistCellStep;
SPEXTERN int StarFeedbackDistTotalCells;
SPEXTERN float StarFeedbackKineticFraction;
SPEXTERN int   StarFeedbackUseTabularYields; 
SPEXTERN char* StarFeedbackTabularFilename;
SPEXTERN float StarFeedbackTabularSNIIEnergy;
SPEXTERN float StarFeedbackTabularSNIaEnergy;
SPEXTERN int   StarFeedbackTrackMetalSources;
SPEXTERN float StarMakerExplosionDelayTime;
SPEXTERN int   StarMakerUseJeansMass;
SPEXTERN int   StarMakerVelDivCrit;
SPEXTERN int   StarMakerSelfBoundCrit;
SPEXTERN int   StarMakerThermalCrit;
SPEXTERN int   StarMakerH2Crit;
SPEXTERN int   StarMakerStochasticStarFormation;
SPEXTERN float StarMakerTemperatureThreshold;
SPEXTERN int   StarMakerStoreInitialMass;

SPEXTERN float PopIIIStarMass;
SPEXTERN int   PopIIIInitialMassFunction;
SPEXTERN int   PopIIIInitialMassFunctionSeed;
SPEXTERN int   PopIIIInitialMassFunctionCalls;
SPEXTERN float PopIIILowerMassCutoff;
SPEXTERN float PopIIIUpperMassCutoff;
SPEXTERN float PopIIIInitialMassFunctionSlope;
SPEXTERN int   PopIIIBlackHoles;
SPEXTERN float PopIIIBHLuminosityEfficiency;
SPEXTERN float PopIIIOverDensityThreshold;
SPEXTERN float PopIIIH2CriticalFraction;
SPEXTERN float PopIIIMetalCriticalFraction;
SPEXTERN int   PopIIIHeliumIonization;
SPEXTERN float PopIIISupernovaRadius;
SPEXTERN int   PopIIISupernovaUseColour;
SPEXTERN int   PopIIISupernovaMustRefine;
SPEXTERN int   PopIIISupernovaMustRefineResolution;
SPEXTERN float PopIIIColorDensityThreshold;
SPEXTERN float PopIIIColorMass;
SPEXTERN int   PopIIIUseHypernova;
SPEXTERN int   PopIIISupernovaExplosions;
SPEXTERN int   PopIIIOutputOnFeedback;

SPEXTERN int    StarClusterUseMetalField;
SPEXTERN int    StarClusterHeliumIonization;
SPEXTERN float  StarClusterMinDynamicalTime;
SPEXTERN double StarClusterIonizingLuminosity;
SPEXTERN double StarClusterSNEnergy;
SPEXTERN float  StarClusterSNRadius;
SPEXTERN float  StarClusterFormEfficiency;
SPEXTERN float  StarClusterMinimumMass;
SPEXTERN float  StarClusterCombineRadius;
SPEXTERN int    StarClusterUnresolvedModel;
SPEXTERN float  StarClusterRegionLeftEdge[3];
SPEXTERN float  StarClusterRegionRightEdge[3];

SPEXTERN float  MBHMinDynamicalTime;
SPEXTERN float  MBHMinimumMass;
SPEXTERN int    MBHAccretion;
SPEXTERN float  MBHAccretionRadius;
SPEXTERN float  MBHAccretingMassRatio;
SPEXTERN float  MBHAccretionFixedTemperature;
SPEXTERN float  MBHAccretionFixedRate;
SPEXTERN int    MBHTurnOffStarFormation;
SPEXTERN float  MBHCombineRadius;

SPEXTERN float UnfulfilledStarFormationMass;

SPEXTERN int    MBHFeedback;
SPEXTERN float  MBHFeedbackRadiativeEfficiency;
SPEXTERN float  MBHFeedbackEnergyCoupling;
SPEXTERN float  MBHFeedbackMassEjectionFraction;
SPEXTERN float  MBHFeedbackMetalYield;
SPEXTERN float  MBHFeedbackThermalRadius;
SPEXTERN float  MBHFeedbackJetsThresholdMass;

SPEXTERN int    H2StarMakerH2FractionMethod;
SPEXTERN float  H2StarMakerEfficiency;
SPEXTERN float  H2StarMakerNumberDensityThreshold;
SPEXTERN float  H2StarMakerMinimumMass;
SPEXTERN float  H2StarMakerMinimumH2FractionForStarFormation;
SPEXTERN int    H2StarMakerStochastic;
SPEXTERN int    H2StarMakerUseSobolevColumn;
SPEXTERN float  H2StarMakerSigmaOverR;
SPEXTERN int    H2StarMakerAssumeColdWarmPressureBalance;
SPEXTERN int    H2StarMakerUseLocalDensityMax;
SPEXTERN float  H2StarMakerH2DissociationFlux_MW;
SPEXTERN float  H2StarMakerH2FloorInColdGas;
SPEXTERN float  H2StarMakerColdGasTemperature;
SPEXTERN int    H2StarMakerWriteStarLogFiles;

SPEXTERN int AccretingParticleRadiation;
SPEXTERN double AccretingParticleLuminosity;

SPEXTERN float minStarLifetime;
SPEXTERN FLOAT LastSupernovaTime;
SPEXTERN float *IMFData;

/* for star particle minimum mass ramp */
SPEXTERN int StarMakerMinimumMassRamp;
SPEXTERN float StarMakerMinimumMassRampStartTime;
SPEXTERN float StarMakerMinimumMassRampStartMass;
SPEXTERN float StarMakerMinimumMassRampEndTime;
SPEXTERN float StarMakerMinimumMassRampEndMass;

SPEXTERN int StarFeedbackThermalEfficiencyRamp;
SPEXTERN float StarFeedbackThermalEfficiencyRampStartTime;
SPEXTERN float StarFeedbackThermalEfficiencyRampStartValue;
SPEXTERN float StarFeedbackThermalEfficiencyRampEndTime;
SPEXTERN float StarFeedbackThermalEfficiencyRampEndValue;

#endif

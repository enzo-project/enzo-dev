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

#ifdef DEFINE_STORAGE
# define SPEXTERN
#else /* DEFINE_STORAGE */
# define SPEXTERN extern
#endif /* DEFINE_STORAGE */

#define STAR_PARTICLE_NUMBER_START 1000000000

struct ParticleEntry {
  FLOAT Position[3];
  float Mass;
  float Velocity[3];
  float Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  int Number;
  int Type;
};


/* Number of Star particles. */

SPEXTERN int NumberOfStarParticles;
SPEXTERN int NumberOfDeletedParticles;
SPEXTERN int NumberOfOtherParticles; //all the particles other than type=2
SPEXTERN int G_TotalNumberOfStars;

/* Star particle parameters. */

SPEXTERN int StarFeedbackType;
SPEXTERN float StarMakerOverDensityThreshold;
SPEXTERN float StarMakerMassEfficiency;
SPEXTERN float StarMakerMinimumMass;
SPEXTERN float StarMakerMinimumDynamicalTime;
SPEXTERN float StarMassEjectionFraction;
SPEXTERN float StarMetalYield;
SPEXTERN float StarEnergyToThermalFeedback;
SPEXTERN float StarEnergyFeedbackRate;
SPEXTERN float StarEnergyToStellarUV;
SPEXTERN float StarEnergyToQuasarUV;

SPEXTERN float PopIIIStarMass;
SPEXTERN int   PopIIIBlackHoles;
SPEXTERN float PopIIIBHLuminosityEfficiency;
SPEXTERN float PopIIIOverDensityThreshold;
SPEXTERN float PopIIIH2CriticalFraction;
SPEXTERN float PopIIIMetalCriticalFraction;
SPEXTERN float PopIIISupernovaRadius;
SPEXTERN int   PopIIISupernovaUseColour;
SPEXTERN float PopIIIColorDensityThreshold;
SPEXTERN float PopIIIColorMass;

SPEXTERN int    StarClusterUseMetalField;
SPEXTERN float  StarClusterMinDynamicalTime;
SPEXTERN double StarClusterIonizingLuminosity;
SPEXTERN double StarClusterSNEnergy;
SPEXTERN float  StarClusterSNRadius;
SPEXTERN float  StarClusterFormEfficiency;
SPEXTERN float  StarClusterMinimumMass;
SPEXTERN float  StarClusterCombineRadius;
SPEXTERN float  StarClusterRegionLeftEdge[3];
SPEXTERN float  StarClusterRegionRightEdge[3];

SPEXTERN float  MBHMinDynamicalTime;
SPEXTERN float  MBHMinimumMass;
SPEXTERN int    MBHAccretion;
SPEXTERN float  MBHAccretingMassRatio;
SPEXTERN int    MBHFeedbackThermal;
SPEXTERN float  MBHFeedbackRadius;
SPEXTERN float  MBHFeedbackRadiativeEfficiency;
SPEXTERN float  MBHFeedbackThermalCoupling;
SPEXTERN float  MBHFeedbackMassEjectionFraction;
SPEXTERN float  MBHFeedbackMetalYield;
SPEXTERN float  MBHCombineRadius;

SPEXTERN float minStarLifetime;
SPEXTERN FLOAT LastSupernovaTime;

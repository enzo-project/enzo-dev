/***********************************************************************
/
/  COMPUTES THE CLUMPING FACTOR (TO COMPARE WITH BEN MATHIESON)
/
/  written by: Greg Bryan
/  date:       March, 2000
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../enzo/macros_and_parameters.h"
#include "../enzo/typedefs.h"
#include "../enzo/global_data.h"
#include "../enzo/Fluxes.h"
#include "../enzo/GridList.h"
#include "../enzo/ExternalBoundary.h"
#include "../enzo/Grid.h"
#include "../enzo/Hierarchy.h"
#include "../enzo/LevelHierarchy.h"
#include "../enzo/TopGridData.h"
#include "../enzo/CosmologyParameters.h"

/* function prototypes */

void QuickSortAndDragFloat(float List[], int left, int right, int NumberToDrag,
			   float *DragList[]);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


#define NUMBER_OF_BINS 18


int AnalyzeClusterComputeClumpingFactor(LevelHierarchyEntry *LevelArray[],
					TopGridData *MetaData,
					int NumberOfGridPoints,
					int NumberOfParticles,
					FLOAT SphereCenter[], 
					float SphereRadius, char *Name)
{

  /* Declarations */

  int i, j, level;
  char ClumpingName[MAX_LINE_LENGTH];
  const float pi = 3.14159;

  if (ComovingCoordinates != 1) {  
    InitialRedshift = 0; 
    //FinalRedshift = 0;
    HubbleConstantNow = 0.7; 
    OmegaMatterNow = 0.3;
    OmegaLambdaNow = 0.7;
    //    float ComovingBoxSize = 1;
    //    float MaxExpansionRate = 1;
  }  

  /* First, allocate space. */

  float *ParticleRadius, *ParticleDensity, *ParticleVolume;
  int TotalNumber = NumberOfGridPoints + NumberOfParticles;

  ParticleRadius = new float[TotalNumber];
  ParticleDensity = new float[TotalNumber];
  ParticleVolume = new float[TotalNumber];
  if (ParticleVolume == NULL) {
    fprintf(stderr, "malloc problem (out of memory?)\n");
    return FAIL;
  }

  /* Next, loop through hierarchy and collect particles & grid points. */

  LevelHierarchyEntry *Temp;
  int ParticleCount = 0;
  if (debug) printf("->computing special clumping factor\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      if (Temp->GridData->CollectParticleInfo(SphereCenter, SphereRadius, 
		                   &ParticleCount, ParticleRadius, 
				   ParticleDensity, ParticleVolume) == FAIL) {
	fprintf(stderr, "Error in grid->AnalyzeClusterCollectParticleInfo\n");
	return FAIL;
      }
      Temp = Temp->NextGridThisLevel;
    }
  }

  if (debug)
    printf("TotalNumberOfParticles = %d ParticleCount = %d\n", 
	   TotalNumber, ParticleCount);

  /* Sort particles by radius. */

  float *DragList[2];
  DragList[0] = ParticleDensity;
  DragList[1] = ParticleVolume;
  QuickSortAndDragFloat(ParticleRadius, 0, TotalNumber-1, 2, DragList);

  /* Compute ratio of mean density to crital density. */

  FLOAT a = 1, dadt, CurrentRedshift = 0;
  float OmegaCurvatureNow = 0.0, Esquared, MeanToCriticalDensityRatio = 1;

  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(LevelArray[0]->GridData->ReturnTime(),
					&a, &dadt) == FAIL) {
      fprintf(stderr, "Error in ComputeExpansionFactor.\n");
      return FAIL;
    }

    CurrentRedshift = (1.0+InitialRedshift)/a - 1.0;

    OmegaCurvatureNow = 1.0 - OmegaMatterNow - OmegaLambdaNow;
    Esquared = OmegaMatterNow    * pow(1.0+CurrentRedshift, 3) +
               OmegaCurvatureNow * pow(1.0+CurrentRedshift, 2) +
               OmegaLambdaNow;
    MeanToCriticalDensityRatio = OmegaMatterNow * pow(1.0+CurrentRedshift, 3)
                                 / Esquared;
  }

  /* Compute run of density as a function of radius. */

  float BinDelta[NUMBER_OF_BINS], BinRadius[NUMBER_OF_BINS], 
        BinDensity[NUMBER_OF_BINS], BinDensity2[NUMBER_OF_BINS],
        BinVolume[NUMBER_OF_BINS];
  int   BinNumber[NUMBER_OF_BINS];
  for (j = 0; j < NUMBER_OF_BINS; j++) {
    BinDelta[j] = pow(10.0, 4.3-float(j)*0.2);
    BinNumber[j] = 0;
    BinDensity[j] = BinDensity2[j] = BinVolume[j] = 0.0;
  }

  float Mass = 0, Volume, delta;
  j = 0;
  for (i = 0; i < TotalNumber; i++) {

    /* Compute delta: the mean density interior to the point divided by
       the critical density. Mass and volume are in code units, so their
       ratio is the density relative to the mean density. */

    Mass += fabs(ParticleDensity[i]*ParticleVolume[i]);
    Volume = 4.0/3.0*pi*pow(ParticleRadius[i], 3);
    delta = Mass/Volume*MeanToCriticalDensityRatio;
    while (delta < BinDelta[j])
      BinRadius[j++] = ParticleRadius[max(i-1, 0)];

    /* if the particle is a gas paraticle (i.e. mass is positive), then
       collect information in this delta bin. */

    if (ParticleDensity[i] > 0) {
      BinNumber[j]++;
      BinDensity[j] += ParticleDensity[i]*ParticleVolume[i];
      BinDensity2[j] += ParticleDensity[i]*ParticleDensity[i]*
	                ParticleVolume[i];
      BinVolume[j] += ParticleVolume[i];

    }

  }

  /* Compute and output clumping factors. */

  strcpy(ClumpingName, Name);
  strcat(ClumpingName, ".ClumpingFactor");
  FILE *fptr = fopen(ClumpingName, "w");
  float VolumeShell, Clump, MeanClump = 0, TotalMass = 0;
  fprintf(fptr, "# bin #   radius    delta    density  density^2  volume1  volume2  f_clump  # points\n");
  for (j = 0; j < NUMBER_OF_BINS; j++) 
    if (BinNumber[j] > 0) {
      VolumeShell = 4.0/3.0*pi*pow(BinRadius[j], 3);
      if (j > 0)
	VolumeShell -= 4.0/3.0*pi*pow(BinRadius[j-1], 3);
      Mass = BinDensity[j];
      BinDensity[j] /= VolumeShell;
      BinDensity2[j] /= VolumeShell;
      Clump = BinDensity2[j]/(BinDensity[j]*BinDensity[j]);
      if (BinDelta[j] > 500) {
	MeanClump += Clump * Mass;
	TotalMass += Mass;
      }
      fprintf(fptr, "%d %g %g %g %g %g %g %g %d\n", j, BinRadius[j], 
	      BinDelta[j], BinDensity[j]*MeanToCriticalDensityRatio,
	      BinDensity2[j]*pow(MeanToCriticalDensityRatio, 2),
	      BinVolume[j], VolumeShell, Clump, BinNumber[j]);
    }
  MeanClump = MeanClump/TotalMass;
  fprintf(fptr, "# mean clumping within r500 = %g\n", MeanClump);
  fclose(fptr);

  return SUCCESS;
}

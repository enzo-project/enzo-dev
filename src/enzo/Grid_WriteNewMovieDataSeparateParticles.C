/*------------------------------------------------------------------------
  WRITE MOVIE DATAFILE (particles only in a separate hdf5 file)
  
  Created : 05 Apr 2010, Ji-hoon Kim
  Modified: 

  Purpose : To output particle data into a separate single hdf5 file, 
            while being grouped by WriteTime.

  History : 

------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"
#include "AMRH5writer.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);

int grid::WriteNewMovieDataSeparateParticles(FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[], 
					     FLOAT StopTime, AMRHDF5Writer &AmiraGrid,
					     int lastMovieStep, int WriteMe, 
					     FLOAT WriteTime, 
					     int alreadyopened[], int NumberOfStarParticlesOnProc[])
{

  /* Declarations */
  int i, j, k, dim;
  int dataWritten;
  FLOAT CurrentRedshift, Left[MAX_DIMENSION], Right[MAX_DIMENSION], a = 1, dadt;

  double doubleTime, doubleRedshift;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NewMovieParticleOn != NON_DM_PARTICLES_MERGED_ALL)
    return SUCCESS;

  /* Exit if last timestep, but not called from WriteDataHierarchy */
  if (!lastMovieStep && Time >= StopTime)
    return SUCCESS;

  // Flag to write data if it's the last timestep or the n-th timestep
  dataWritten = (lastMovieStep || WriteMe);

  // If not n-th timestep, exit
  if (!dataWritten) return SUCCESS;

  /* If outside the region, skip the grid */
  
  for (dim = 0; dim < GridRank; dim++) {
    Left[dim] = max(RegionLeftEdge[dim], GridLeftEdge[dim]);
    Right[dim] = min(RegionRightEdge[dim], GridRightEdge[dim]);
    if (Left[dim] >= Right[dim])
      return SUCCESS;
  }

  /* Get expansion factor */

  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      return FAIL;
    }
    CurrentRedshift = (1 + InitialRedshift)/a - 1;
  } else
    CurrentRedshift = -1.0;

  /* Create attributes for this grid's entry */

  doubleTime = (double) WriteTime;
  doubleRedshift = (double) CurrentRedshift;

  /* Output particles */

  if (NewMovieParticleOn == NON_DM_PARTICLES_MERGED_ALL) {

    /* Search for non dark matter particles and record their array
       element.  We don't store the data first because we need the
       particle number first. */

    int *NonDMParticleIndices = new int[NumberOfParticles];
    int NumberOfNonDMParticles = 0;
    int ii, iattr;

    FLOAT *TempPosition[3];
    float *TempVelocity[3], *TempMass;
    float *TempAttr[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    int *TempType;
    PINT *TempNumber;
    
    for (i = 0; i < NumberOfParticles; i++)
      NonDMParticleIndices[i] = -1;
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleType[i] != PARTICLE_TYPE_DARK_MATTER)
	NonDMParticleIndices[NumberOfNonDMParticles++] = i;

    /* Allocate memory */

    if (NumberOfNonDMParticles > 0) {
      for (dim = 0; dim < GridRank; dim++) {
	TempPosition[dim] = new FLOAT[NumberOfNonDMParticles];
	TempVelocity[dim] = new float[NumberOfNonDMParticles];
      }
      TempMass = new float[NumberOfNonDMParticles];
      for (i = 0; i < NumberOfParticleAttributes; i++)
	TempAttr[i] = new float[NumberOfNonDMParticles];
      TempType = new int[NumberOfNonDMParticles];
      TempNumber = new PINT[NumberOfNonDMParticles];
    } // ENDIF non-DM particles > 0

    /* Move non-DM particles into temp arrays */

    for (i = 0; i < NumberOfNonDMParticles; i++) {
      j = NonDMParticleIndices[i];
      for (dim = 0; dim < GridRank; dim++) {
	TempPosition[dim][i] = ParticlePosition[dim][j];
	TempVelocity[dim][i] = ParticleVelocity[dim][j];
      }
      TempMass[i] = ParticleMass[j];
      for (iattr = 0; iattr < NumberOfParticleAttributes; iattr++)
	TempAttr[iattr][i] = ParticleAttribute[iattr][j];
      TempType[i] = ParticleType[j];
      TempNumber[i] = ParticleNumber[j];
    } // ENDFOR non-DM particles

    //    fprintf(stdout, "alreadyopened[] = %d\n", alreadyopened[ProcessorNumber]);  

    if (AmiraGrid.writeSeparateParticles(NumberOfNonDMParticles, NumberOfParticleAttributes,
					 GridRank, 
					 (void **) TempPosition, 
					 (void **) TempVelocity,
					 (void *) TempType, (void *) TempNumber, 
					 (void *) TempMass,
					 (void **) TempAttr, 
					 doubleTime, doubleRedshift, 
					 alreadyopened[ProcessorNumber],
					 NumberOfStarParticlesOnProc[ProcessorNumber]) 
	!= 0) {
      fprintf(stderr, "Error in AMRHDF5Writer->writeSeparateParticles\n");
      return FAIL;
    }

    /* Free memory */

    if (NumberOfNonDMParticles > 0) {
      for (dim = 0; dim < GridRank; dim++) {
	delete [] TempPosition[dim];
	delete [] TempVelocity[dim];
      }
      for (i = 0; i < NumberOfParticleAttributes; i++)
	delete [] TempAttr[i];
      delete [] TempMass;
      delete [] TempType;
      delete [] TempNumber;
    } // ENDIF non-DM particles > 0

    delete [] NonDMParticleIndices;

  } // ENDIF Non-DM particles

  AmiraGrid.IncreaseGridCount();

  return SUCCESS;

}

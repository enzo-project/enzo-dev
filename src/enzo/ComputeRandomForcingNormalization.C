/***********************************************************************
/
/  COMPUTE RANDOM FORCING NORMALIZATION
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
 
/* ======================================================================= */
/* Function prototypes. */
 
int   GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
			HierarchyEntry **Grids[]);
 
/* This routine calculates the normalization for Random Forcing on
   level 0 (since this involves communication). */
 
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
				      int level, TopGridData *MetaData,
				      float * norm, float * pTopGridTimeStep)
{
  /* Return if this does not concern us */
  if (!RandomForcing) return SUCCESS;
 
  LCAPERF_START("ComputeRandomForcingNormalization");

  /* If level is above 0 then complain: forcing will only work on level 0
     grid(s). */

  if ((MetaData->CycleNumber <= 0) || (level != 0)) {
      LCAPERF_STOP("ComputeRandomForcingNormalization");
      return SUCCESS;
  }
 
  /* Create an array (Grids) of all the grids on level 0. */
 
  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
 
  /* Loop over level 0 grids and compute sums over each of the grids. */
 
  int grid, GlobNum = 9;
  float * GlobVal = new float[GlobNum];
  for (int num = 0; num < GlobNum; num++)
    GlobVal[num] = 0.0;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->PrepareRandomForcingNormalization(GlobVal,
								 GlobNum)
	== FAIL) {
      ENZO_FAIL("Error in grid->PrepareRandomForcingNormalization.");
    }
 
  /* Communicate grid-specific sums and compute global sums;
     also communicate min/max Density values. */
 
  CommunicationAllSumValues(GlobVal, GlobNum-2);
  float minDens = CommunicationMinValue(*(GlobVal+GlobNum-2));
  float maxDens = CommunicationMaxValue(*(GlobVal+GlobNum-1));
 
  /* Compute normalization. */
 
  int numberOfGridZones = 1;
  for (int dim = 0; dim <  MetaData->TopGridRank; dim++)
    numberOfGridZones *= MetaData->TopGridDims[dim];
 
  float dt = LevelArray[level]->GridData->ReturnTimeStep();
  *pTopGridTimeStep = dt;
  if (RandomForcingEdot == 0.0)
    *norm = 0.0;
  else
    *norm = ( sqrt(GlobVal[0]*GlobVal[0] + GlobVal[1]*dt*RandomForcingEdot*2.0*
		 numberOfGridZones) - GlobVal[0] )/GlobVal[1];
 
  if (debug) printf("RandomForcingNormalization %.10"GSYM"\n", *norm);
  if (debug) printf("RandomForcingGlobals: E_k %"GSYM" M_m %"GSYM" M_v %"GSYM" \n",
		    0.5*GlobVal[4]/numberOfGridZones,
		    sqrt(GlobVal[2]/numberOfGridZones),
		    sqrt(GlobVal[3]/numberOfGridZones));
 
 
  /* Output results for level 0 grids. */
 
  if (MetaData->CycleSkipGlobalDataDump != 0)
    if (MyProcessorNumber == ROOT_PROCESSOR &&
	MetaData->CycleNumber % MetaData->CycleSkipGlobalDataDump == 0.0) {
      FILE * Fptr;
      if ((Fptr = fopen("randomForcing.out", "a")) == NULL)
	ERROR_MESSAGE;
      fprintf( Fptr, "%"ISYM" %9.6"FSYM" ", MetaData->CycleNumber, MetaData->Time);
      fprintf( Fptr, " %9.3"GSYM, GlobVal[1]); //gv1
      fprintf( Fptr, " %9.3"GSYM, GlobVal[2]); //gv1
      fprintf( Fptr, " %9.3"GSYM, GlobVal[3]); //gv1
      fprintf( Fptr, " %9.3"GSYM, GlobVal[4]); //gv1
      fprintf( Fptr, " %9.3"GSYM, GlobVal[5]); //gv1
      fprintf( Fptr, " %9.3"GSYM, GlobVal[6]); //gv1
      fprintf( Fptr, " %9.3"GSYM, GlobVal[7]); //gv1
      fprintf( Fptr, " %9.3"GSYM, GlobVal[8]); //gv1
      fprintf( Fptr, " %9.3"GSYM, GlobVal[9]); //gv1
      fprintf( Fptr, " %9.6"FSYM, 0.50*GlobVal[4]/numberOfGridZones);   // kinetic energy
      fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[2]/numberOfGridZones));  // mass weighted rms Mach
      fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[3]/numberOfGridZones));  // volume weighed rms Mach
      fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[5]/numberOfGridZones));  // rms Velocity
      fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[6]/numberOfGridZones));  // Density variance
      fprintf( Fptr, " %9.6"FSYM" %10.5"FSYM"\n", minDens, maxDens );                  // min/max Density
      fclose(Fptr);
    }
 
  /* clean up. */
 
  delete [] GlobVal;
 
  LCAPERF_STOP("ComputeRandomForcingNormalization");
  return SUCCESS;
}

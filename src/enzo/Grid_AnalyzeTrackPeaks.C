/***********************************************************************
/
/  GRID CLASS (ANALYZE DENSITY FIELD FOR PEAKS AND RECORD)
/
/  written by: Greg Bryan
/  date:       December, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "CosmologyParameters.h"
#include "Grid.h"
 
#define MAX_PEAKS 200
 
int grid::AnalyzeTrackPeaks(int level, int ReportLevel)
{
 
  static int PeakActive[MAX_DEPTH_OF_HIERARCHY][MAX_PEAKS];
  static float PeakDensity[MAX_DEPTH_OF_HIERARCHY][MAX_PEAKS];
  static FLOAT PeakLastUpdateTime[MAX_DEPTH_OF_HIERARCHY][MAX_PEAKS];
  static FLOAT PeakPosition[MAX_DEPTH_OF_HIERARCHY][MAX_PEAKS][MAX_DIMENSION];
 
  /* Return if this grid is not on this processor
     (this routine is not multi-processor safe). */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* ReportLevel == 0 means no reporting. */
 
  if (ReportLevel == 0 || NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* For ReportLevel > 1, ignore level 0. */
 
//  if (level == 0 && ReportLevel > 1)
//    return SUCCESS;
 
  /* declarations */
 
  int dim, n, i, j, k, i1, j1, k1, l, index, index1,
      FoundIt, OutputAlready[MAX_PEAKS];
  FLOAT pos[MAX_DIMENSION], radius;
  FILE *fptr;
  char PeakOutputName[MAX_LINE_LENGTH];
 
  for (n = 0; n < MAX_PEAKS; n++)
    OutputAlready[n] = FALSE;
 
  /* Allocate field and compute temperature. */
 
  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  float *temperature = new float[size];
  if (this->ComputeTemperatureField(temperature) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");
  }
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* Open output file. */
 
  sprintf(PeakOutputName, "%s.L%1.1"ISYM, "PeakData", level);
 
  if ((fptr = fopen(PeakOutputName, "a")) == FAIL) {
    ENZO_VFAIL("Error opening %s.\n", PeakOutputName)
  }
 
  /* Compute the MinimumPeakDensity. */
 
  float MinimumPeakDensity = 0.9*MinimumMassForRefinement[0]/
    (CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0]);
  if (level > 0)
    MinimumPeakDensity /= POW(float(RefineBy), GridRank);
 
  /* Remove old peaks (more than 30 timesteps ago). */
 
#ifdef UNUSED
  for (n = 0; n < MAX_PEAKS; n++)
    if (PeakActive[level][n] == TRUE)
      if (PeakLastUpdateTime[level][n] < Time - 30.0*dtFixed)
	PeakActive[level][n] = FALSE;
#endif /* UNUSED */
 
  /* Loop over grid looking for new peaks. */
 
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
 
	index = (k*GridDimension[1] + j)*GridDimension[0] + i;
 
	/* If the density is too low, ignore it. */
 
	if (BaryonField[DensNum][index] < MinimumPeakDensity)
	   continue;
 
	/* For level 0, make sure it is in the allowed refined region. */
 
	if (level == 0)
	  if (CellLeftEdge[0][i] < RefineRegionLeftEdge[0] ||
	      CellLeftEdge[1][j] < RefineRegionLeftEdge[1] ||
	      CellLeftEdge[2][k] < RefineRegionLeftEdge[2] ||
	      CellLeftEdge[0][i+1] > RefineRegionRightEdge[0] ||
	      CellLeftEdge[1][j+1] > RefineRegionRightEdge[1] ||
	      CellLeftEdge[2][k+1] > RefineRegionRightEdge[2]   )
	    continue;
 
	/* Check to see if it is a local max. */
 
	for (k1 = -1; k1 <= 1; k1++)
	  for (j1 = -1; j1 <= 1; j1++)
	    for (i1 = -1; i1 <= 1; i1++) {
 
	      index1 = index +
		(k1*GridDimension[1] + j1)*GridDimension[0] + i1;
 
	      if (BaryonField[DensNum][index1] >= BaryonField[DensNum][index]
		  && index != index1)
		goto NotLocalMax;
	    }
 
	pos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	pos[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
 
	/* Found a local maxima, check if we already have it in the list */
 
	FoundIt = FALSE;
	for (n = 0; n < MAX_PEAKS; n++)
	  if (PeakActive[level][n] == TRUE) {
	    radius = 0;
	    for (dim = 0; dim < GridRank; dim++)
	      radius += (pos[dim] - PeakPosition[level][n][dim])*
		        (pos[dim] - PeakPosition[level][n][dim]);
	    radius = sqrt(radius);
 
	    /* Found it */
 
	    if (radius < 1.9*CellWidth[0][i]) {
	      FoundIt = TRUE;
	      break;
	    }
 
	  } // end: loop over peaks
 
	/* If we didn't find it, add it (if possible) to the list. */
 
	if (FoundIt == FALSE)
	  for (n = 0; n < MAX_PEAKS; n++)
	    if (PeakActive[level][n] == FALSE) {
	      PeakActive[level][n] = TRUE;
	      FoundIt = TRUE;
	      break;
	    }
	
	/* If we found a spot, Update time and output to file. */
 
	if (FoundIt == TRUE && OutputAlready[n] == FALSE) {
	  OutputAlready[n] = TRUE;
	  PeakLastUpdateTime[level][n] = Time;
	  for (dim = 0; dim < GridRank; dim++)
	    PeakPosition[level][n][dim] = pos[dim];
	  PeakDensity[level][n] = BaryonField[DensNum][index];
 
	  fprintf(fptr, "peak=%"ISYM" ", n);
 
	  for (l = 0; l < level; l++)
	    for (int n1 = 0; n1 < MAX_PEAKS; n1++)
	      if (PeakActive[l][n1] == TRUE) {
		radius = 0;
		for (dim = 0; dim < GridRank; dim++)
		  radius += (pos[dim] - PeakPosition[l][n1][dim])*
		            (pos[dim] - PeakPosition[l][n1][dim]);
		radius = sqrt(radius);
 
		if (radius < CellWidth[0][0]*POW(float(RefineBy), level-l))
		  fprintf(fptr, "%"ISYM"l%"ISYM" ", l, n1);
	      }
 
	  fprintf(fptr, "! %"GOUTSYM" ", Time);
 
	  fprintf(fptr, "%"GSYM" %"GSYM" %"GSYM" ", PeakDensity[level][n],
		  temperature[index], temperature[index]/
		  POW(BaryonField[DensNum][index], Gamma-1));
 
/*
	  if (SelfGravity && GravityResolution == 1 &&
	      GravitatingMassFieldParticles != NULL && level > 0)
	    fprintf(fptr, "%"GSYM" ", GravitatingMassFieldParticles[index]);
	  else
	    fprintf(fptr, "%"GSYM" ", tiny_number);
*/
	
	  if (SelfGravity && GravitatingMassFieldParticles != NULL) {

	    int gdims[3], gindex;
	    for (dim = 0; dim < 3; dim++)
	      gdims[dim] =
		nint((pos[dim] - GravitatingMassFieldParticlesLeftEdge[dim])/
		     GravitatingMassFieldParticlesCellSize);
	    gindex = (gdims[2]*GravitatingMassFieldParticlesDimension[1] +
		      gdims[1])*GravitatingMassFieldParticlesDimension[0] +
		      gdims[0];
	    fprintf(fptr, "%"GSYM" ", GravitatingMassFieldParticles[gindex]);
	  }
	  else
	    fprintf(fptr, "%"GSYM" ", tiny_number);
 
	  fprintf(fptr, "%"GOUTSYM" %"GOUTSYM" %"GOUTSYM" ", PeakPosition[level][n][0],
		  PeakPosition[level][n][1], PeakPosition[level][n][2]);
 
	  for (l = 1; l < NumberOfBaryonFields; l++)
	    fprintf(fptr, "%"GSYM" ", BaryonField[l][index]);
 
	  fprintf(fptr, "\n");
 
	}
 
      NotLocalMax:
	FoundIt = FALSE;
	
      } // end: loop over grid
 
  /* Close file. */
 
  fclose(fptr);
  delete temperature;
 
  return SUCCESS;
}

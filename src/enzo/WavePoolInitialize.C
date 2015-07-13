/***********************************************************************
/
/  INITIALIZE A WAVE POOL SIMULATION
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:
/
/  PURPOSE:
/    The wave pool sets up a system which allows a 1d sinusoidal wave to
/    enter from the left boundary.  The initial active region
/    is completely uniform, and wave enters via inflow boundary conditions.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
 
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#define DEFINE_STORAGE
#include "WavePoolGlobalData.h"
#undef DEFINE_STORAGE
 
int WavePoolInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
       SubgridDims[MAX_DIMENSION];
  float WavePoolTotalEnergy;
  FLOAT WavePoolSubgridLeft, WavePoolSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  float ZeroBField[3] = {0.0, 0.0, 0.0};
 
  /* set default parameters */
 
  WavePoolAmplitude     = 0.01;   // linear wave
  WavePoolWavelength    = 0.1;    // one-tenth of the box
  WavePoolNumberOfWaves = 1;      // just one wave
  WavePoolAngle         = 0.0;    // direction of wave propogation wrt x-axis
 
  WavePoolDensity       = 1.0;  // uniform pool
  WavePoolPressure      = 1.0;
  WavePoolVelocity[0]   = 0.0;  // x, y and z velocities
  WavePoolVelocity[1]   = 0.0;
  WavePoolVelocity[2]   = 0.0;
 
  WavePoolSubgridLeft   = 0.0;   // start of subgrid
  WavePoolSubgridRight  = 0.0;  // end of subgrid
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "WavePoolAmplitude = %"FSYM, &WavePoolAmplitude);
    ret += sscanf(line, "WavePoolWavelength = %"FSYM, &WavePoolWavelength);
    ret += sscanf(line, "WavePoolNumberOfWaves = %"FSYM, &WavePoolNumberOfWaves);
    ret += sscanf(line, "WavePoolAngle = %"FSYM, &WavePoolAngle);
 
    ret += sscanf(line, "WavePoolDensity = %"FSYM, &WavePoolDensity);
    ret += sscanf(line, "WavePoolPressure = %"FSYM, &WavePoolPressure);
    ret += sscanf(line, "WavePoolVelocity1 = %"FSYM, &WavePoolVelocity[0]);
    ret += sscanf(line, "WavePoolVelocity2 = %"FSYM, &WavePoolVelocity[1]);
    ret += sscanf(line, "WavePoolVelocity3 = %"FSYM, &WavePoolVelocity[2]);
 
    ret += sscanf(line, "WavePoolSubgridLeft = %"PSYM, &WavePoolSubgridLeft);
    ret += sscanf(line, "WavePoolSubgridRight = %"PSYM, &WavePoolSubgridRight);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "WavePool"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  } // end input from parameter file
 
  /* set the inflow boundary on the left, otherwise leave things alone. */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    MetaData.LeftFaceBoundaryCondition[dim] = inflow;
 
  /* compute total energy */
 
  WavePoolTotalEnergy = WavePoolPressure/((Gamma - 1.0)*WavePoolDensity);
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    WavePoolTotalEnergy += 0.5*WavePoolVelocity[dim]*WavePoolVelocity[dim];
 
  /* set up grid */
 
  if (TopGrid.GridData->InitializeUniformGrid(WavePoolDensity,
					      WavePoolTotalEnergy,
					      WavePoolTotalEnergy,
					      WavePoolVelocity,
                          ZeroBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }
 
  /* If requested, create a subgrid */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((WavePoolSubgridRight - WavePoolSubgridLeft)/
	   ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	    FLOAT(MetaData.TopGridDims[dim])))
	*RefineBy;
 
  if (NumberOfSubgridZones[0] > 0) {
 
    /* create a new HierarchyEntry, attach to the top grid and fill it out */
 
    HierarchyEntry *Subgrid    = new HierarchyEntry;
    TopGrid.NextGridNextLevel  = Subgrid;
    Subgrid->NextGridNextLevel = NULL;
    Subgrid->NextGridThisLevel = NULL;
    Subgrid->ParentGrid        = &TopGrid;
 
    /* Compute the dimensions and left/right edges for the subgrid. */
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*NumberOfGhostZones;
      LeftEdge[dim]    = WavePoolSubgridLeft;
      RightEdge[dim]   = WavePoolSubgridRight;
    }
 
    /* create a new subgrid and initialize it */
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->InitializeUniformGrid(WavePoolDensity,
						 WavePoolTotalEnergy,
						 WavePoolTotalEnergy,
						 WavePoolVelocity,
                         ZeroBField) == FAIL) {
            ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
    }			
  }
 
  /* set up field names and units */
 
  DataLabel[0] = DensName;
  DataLabel[1] = TEName;
  DataLabel[2] = Vel1Name;
  DataLabel[3] = Vel2Name;
  DataLabel[4] = Vel3Name;
 
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "WavePoolAmplitude     = %"FSYM"\n", WavePoolAmplitude);
    fprintf(Outfptr, "WavePoolWavelength    = %"FSYM"\n", WavePoolWavelength);
    fprintf(Outfptr, "WavePoolNumberOfWaves = %"FSYM"\n", WavePoolNumberOfWaves);
    fprintf(Outfptr, "WavePoolAngle         = %"FSYM"\n\n", WavePoolAngle);
 
    fprintf(Outfptr, "WavePoolDensity       = %"FSYM"\n", WavePoolDensity);
    fprintf(Outfptr, "WavePoolPressure      = %"FSYM"\n", WavePoolPressure);
    fprintf(Outfptr, "WavePoolVelocity1     = %"FSYM"\n", WavePoolVelocity[0]);
    fprintf(Outfptr, "WavePoolVelocity2     = %"FSYM"\n", WavePoolVelocity[1]);
    fprintf(Outfptr, "WavePoolVelocity3     = %"FSYM"\n\n", WavePoolVelocity[2]);
 
    fprintf(Outfptr, "WavePoolSubgridLeft   = %"GOUTSYM"\n", WavePoolSubgridLeft);
    fprintf(Outfptr, "WavePoolSubgridRight  = %"GOUTSYM"\n\n", WavePoolSubgridRight);
  }
 
  return SUCCESS;
 
}

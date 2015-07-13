/***********************************************************************
/
/  Rotating Disk Test Problem
/
/  written by:  Elizabeth Tasker
/  date:        May 2012
/  modified1:  
/
/  PURPOSE: sets up a simple galaxy disk in an NFW profile
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
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
#include "TopGridData.h"
#include "phys_constants.h"

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);


 
int RotatingDiskInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			   TopGridData &MetaData)
{

  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";

  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, ret, level;
     
  float RotatingDiskScaleRadius           = 2265.0; // [pc]
  float RotatingDiskScaleHeight           = 100.0;  // [pc]
  float RotatingDiskOuterRadius           = 12e3;   // [pc]
  float RotatingDiskCentralDensity        = 0.67;   // [Msun/pc3]
  float RotatingDiskTemperature           = 1e4;    // [K]
  float RotatingDiskTotalDMMass           = 1e11;   // [Msun]
  float RotatingDiskDMConcentration       = 10.0; 
  int RotatingDiskRefineAtStart           = 1;
 
  /* read input from file */
	
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {	
    ret = 0;
    
    ret += sscanf(line, "RotatingDiskScaleRadius = %"FSYM,
		  &RotatingDiskScaleRadius);
    ret += sscanf(line, "RotatingDiskScaleHeight = %"FSYM,
		  &RotatingDiskScaleHeight);
    ret += sscanf(line, "RotatingDiskOuterRadius = %"FSYM,
		  &RotatingDiskOuterRadius);
    ret += sscanf(line, "RotatingDiskCentralDensity = %"FSYM,
		  &RotatingDiskCentralDensity);
    ret += sscanf(line, "RotatingDiskDMConcentration = %"FSYM,
		  &RotatingDiskDMConcentration);
    ret += sscanf(line, "RotatingDiskTotalDMMass = %"FSYM,
		  &RotatingDiskTotalDMMass);
    ret += sscanf(line, "RotatingDiskTemperature = %"FSYM,
		  &RotatingDiskTemperature);
    ret += sscanf(line, "RotatingDiskRefineAtStart = %"ISYM,
		  &RotatingDiskRefineAtStart);
    
    /* if the line is suspicious, issue a warning */	
    if (ret == 0 && strstr(line, "=") && strstr(line, "RotatingDisk") && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
		
  } // end input from parameter file

  
  // Initialize the top (root) grid with a uniform gas
  if (TopGrid.GridData->RotatingDiskInitializeGrid(RotatingDiskScaleRadius,
						   RotatingDiskScaleHeight, 
						   RotatingDiskTemperature,
						   RotatingDiskDMConcentration, 
						   RotatingDiskTotalDMMass,
						   RotatingDiskCentralDensity,
						   RotatingDiskOuterRadius) == FAIL) {	 
    ENZO_FAIL("Error in root grid initialise.");
  }

  // Set up initial AMR levels

  if (RotatingDiskRefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */
 
    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */
 
    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	ENZO_FAIL("Error in RebuildHierarchy.");
      }

      if (LevelArray[level+1] == NULL)
	break;
      
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      
       while (Temp != NULL) {
	if (Temp->GridData->RotatingDiskInitializeGrid(RotatingDiskScaleRadius,
						       RotatingDiskScaleHeight, 
						       RotatingDiskTemperature,
						       RotatingDiskDMConcentration, 
						       RotatingDiskTotalDMMass,
						       RotatingDiskCentralDensity,
						       RotatingDiskOuterRadius) == FAIL) {
	  ENZO_FAIL("Error in RotatingDiskInitializeGrid [subgrid]");
	}
	Temp = Temp->NextGridThisLevel;
      }
    } // end loop over levels

    /* Loop back from the bottom, restoring the consistency amoung levels. */
 
    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
		 *Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL) {
	  ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.");
	}
	Temp = Temp->NextGridThisLevel;
      }
    }
 
  } // end: if (RotatingDiskRefineAtStart)

  /* set up field names and units */

  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if(DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;

  for(j=0; j < i; j++)
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */
  
  if (MyProcessorNumber == ROOT_PROCESSOR) {
	
    fprintf(Outfptr, "RotatingDiskRefineAtStart      = %"ISYM"\n",
	    RotatingDiskRefineAtStart);
    fprintf(Outfptr, "RotatingDiskScaleRadius        = %"FSYM"\n",
	    RotatingDiskScaleRadius);
    fprintf(Outfptr, "RotatingDiskScaleHeight        = %"FSYM"\n",
	    RotatingDiskScaleHeight);
    fprintf(Outfptr, "RotatingDiskCentralDensity     = %"FSYM"\n",
	    RotatingDiskCentralDensity);
    fprintf(Outfptr, "RotatingDiskDMConcentration    = %"FSYM"\n",
	    RotatingDiskDMConcentration);
    fprintf(Outfptr, "RotatingDiskTotalDMMass        = %"FSYM"\n",
	    RotatingDiskTotalDMMass);
    fprintf(Outfptr, "RotatingDiskTemperature        = %"FSYM"\n",
	    RotatingDiskTemperature);
    fprintf(Outfptr, "RotatingDiskOuterRadius        = %"FSYM"\n",
	    RotatingDiskOuterRadius);
    
  }

 
#ifdef USE_MPI
	
  /* In Grid_RotatingDiskInitializeGrid.C, we set the values of these two global variables,
     PointSourceGravityConstant and PointSourceGravityCoreRadius. When running in parallel, these
     values need to be shares with all processors. */
	
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Bcast(&PointSourceGravityConstant,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(&PointSourceGravityCoreRadius,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(&PointSourceGravityPosition,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);
 	
#endif
	
  return SUCCESS;

}

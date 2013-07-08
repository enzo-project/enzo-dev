/***********************************************************************
/
/  INITIALIZE A GRAVITY TEST
/
/  written by: Greg Bryan
/  date:       September, 1995
/  modified1:
/
/  PURPOSE:
/    Set up a spherical mass distribution to test the gravity solver.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
 
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
#include "LevelHierarchy.h"
#include "TopGridData.h"
#define DEFINE_STORAGE
#include "TestGravitySphereGlobalData.h"
#undef DEFINE_STORAGE
 
 
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
 
int TestGravitySphereInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *GPotName  = "Grav_Potential";
 
  /* declarations */
 
  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level;
  int   NumberOfSubgridZones[MAX_DIMENSION], SubgridDims[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  /* Error check. */
 
  if (!SelfGravity)
    fprintf(stderr, "TestGravitySphere: gravity is off!?!");
 
  /* set default parameters */
 
        TestGravitySphereInteriorDensity   = 1.0;  // density inside sphere
        TestGravitySphereExteriorDensity   = tiny_number;  // outside sphere
        TestGravitySphereRadius            = 0.1;
        TestGravitySphereType              = 0;    // uniform density
  FLOAT TestGravitySphereSubgridLeft       = 0.0;  // start of subgrid
  FLOAT TestGravitySphereSubgridRight      = 0.0;  // end of subgrid
  int   TestGravitySphereUseBaryons        = TRUE;
  int   TestGravitySphereRefineAtStart     = FALSE;
  FLOAT TestGravitySphereCenter[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    TestGravitySphereCenter[dim] = 0.5*(DomainLeftEdge[dim] +
					DomainRightEdge[dim]);
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "TestGravitySphereInteriorDensity = %"FSYM,
		  &TestGravitySphereInteriorDensity);
    ret += sscanf(line, "TestGravitySphereExteriorDensity = %"FSYM,
		  &TestGravitySphereExteriorDensity);
    ret += sscanf(line, "TestGravitySphereRadius = %"FSYM,
		  &TestGravitySphereRadius);
    ret += sscanf(line, "TestGravitySphereType = %"ISYM,
		  &TestGravitySphereType);
    ret += sscanf(line, "TestGravitySphereSubgridLeft = %"PSYM,
		  &TestGravitySphereSubgridLeft);
    ret += sscanf(line, "TestGravitySphereSubgridRight = %"PSYM,
		  &TestGravitySphereSubgridRight);
    ret += sscanf(line, "TestGravitySphereUseBaryons = %"ISYM,
		  &TestGravitySphereUseBaryons);
    ret += sscanf(line, "TestGravitySphereRefineAtStart = %"ISYM,
		  &TestGravitySphereRefineAtStart);
    ret += sscanf(line, "TestGravitySphereCenter = %"PSYM" %"PSYM" %"PSYM,
		  TestGravitySphereCenter, TestGravitySphereCenter+1,
		  TestGravitySphereCenter+2);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "TestGravitySphere")
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  } // end input from parameter file
 
  /* set up grid */
 
  if (TopGrid.GridData->TestGravitySphereInitializeGrid(
                                          TestGravitySphereInteriorDensity,
                                          TestGravitySphereExteriorDensity,
					  TestGravitySphereRadius,
					  TestGravitySphereType,
					  TestGravitySphereUseBaryons,
					  TestGravitySphereCenter
						  ) == FAIL){
    ENZO_FAIL("Error in TestGravitySphereInitializeGrid.");
  }
 
  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */
 
  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }
 
  /* If requested, refine the grid to the desired level. */
 
  if (TestGravitySphereRefineAtStart) {
 
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
	if (Temp->GridData->TestGravitySphereInitializeGrid(
                                          TestGravitySphereInteriorDensity,
                                          TestGravitySphereExteriorDensity,
					  TestGravitySphereRadius,
					  TestGravitySphereType,
					  TestGravitySphereUseBaryons,
					  TestGravitySphereCenter
							        ) == FAIL) {
	  ENZO_FAIL("Error in TestGravitySphereInitializeGrid.");
	}
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels
 
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
 
  } // end: if (TestGravitySphereRefineAtStart)
 
  /* If requested, create a subgrid */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((TestGravitySphereSubgridRight - TestGravitySphereSubgridLeft)/
	   ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	    FLOAT(MetaData.TopGridDims[dim])))
	*RefineBy;
 
  if (NumberOfSubgridZones[0] > 0) {
 
    /* Error check */
 
    if (TestGravitySphereRefineAtStart) {
      ENZO_FAIL("Cannot RefineAtStart AND create subgrid.");
    }
 
    /* create a new HierarchyEntry, attach to the top grid and fill it out */
 
    HierarchyEntry *Subgrid    = new HierarchyEntry;
    TopGrid.NextGridNextLevel  = Subgrid;
    Subgrid->NextGridNextLevel = NULL;
    Subgrid->NextGridThisLevel = NULL;
    Subgrid->ParentGrid        = &TopGrid;
 
    /* compute the dimensions and left/right edges for the subgrid */
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*NumberOfGhostZones;
      LeftEdge[dim]    = TestGravitySphereSubgridLeft;
      RightEdge[dim]   = TestGravitySphereSubgridRight;
    }
 
    /* create a new subgrid and initialize it */
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->TestGravitySphereInitializeGrid(
                                          TestGravitySphereInteriorDensity,
                                          TestGravitySphereExteriorDensity,
					  TestGravitySphereRadius,
					  TestGravitySphereType,
					  TestGravitySphereUseBaryons,
					  TestGravitySphereCenter)
	== FAIL) {
      ENZO_FAIL("Error in TestGravitySphereInitializeGrid.");
    }			
  }
 
  /* set up field names and units */
 
  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
 
  if (WritePotential){
    DataLabel[count++] = GPotName;
  }

  for (int j = 0; j<count; j++){
    DataUnits[j] = NULL;
  }

  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "TestGravitySphereInteriorDensity   = %"FSYM"\n",
	    TestGravitySphereInteriorDensity);
    fprintf(Outfptr, "TestGravitySphereExteriorDensity   = %"FSYM"\n",
	    TestGravitySphereExteriorDensity);
    fprintf(Outfptr, "TestGravitySphereRadius            = %"FSYM"\n",
	    TestGravitySphereRadius);
    fprintf(Outfptr, "TestGravitySphereType              = %"ISYM"\n",
	    TestGravitySphereType);
    fprintf(Outfptr, "TestGravitySphereSubgridLeft       = %"GOUTSYM"\n",
	    TestGravitySphereSubgridLeft);
    fprintf(Outfptr, "TestGravitySphereSubgridRight      = %"GOUTSYM"\n",
	    TestGravitySphereSubgridRight);
    fprintf(Outfptr, "TestGravitySphereUseBaryons        = %"ISYM"\n",
	    TestGravitySphereUseBaryons);
    fprintf(Outfptr, "TestGravitySphereRefineAtStart     = %"ISYM"\n",
	    TestGravitySphereRefineAtStart);
    fprintf(Outfptr, "TestGravitySphereCenter            = ");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank, TestGravitySphereCenter);
  }
 
  return SUCCESS;
 
}

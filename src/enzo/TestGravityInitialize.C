/***********************************************************************
/
/  INITIALIZE A GRAVITY TEST
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/    We set up a system in which there is one grid point with mass in order
/     to see the resulting acceleration field.  If finer grids are specified,
/     the mass is one grid on the subgrid as well.
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
#include "TopGridData.h"
 
int TestGravityInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char  line[MAX_LINE_LENGTH];
  int   dim, ret;
  int   NumberOfSubgridZones[MAX_DIMENSION], SubgridDims[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  /* Error check. */
 
  if (!SelfGravity)
    fprintf(stderr, "TestGravity: gravity is off!?!");
 
  /* set default parameters */
 
  float TestGravityDensity           = 1.0;  // density of central peak
  FLOAT TestGravitySubgridLeft       = 0.0;  // start of subgrid
  FLOAT TestGravitySubgridRight      = 0.0;  // end of subgrid
  int   TestGravityNumberOfParticles = 0;    // number of test particles
  int   TestGravityUseBaryons        = FALSE;
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "TestGravityDensity = %"FSYM, &TestGravityDensity);
    ret += sscanf(line, "TestGravitySubgridLeft = %"PSYM,
		  &TestGravitySubgridLeft);
    ret += sscanf(line, "TestGravitySubgridRight = %"PSYM,
		  &TestGravitySubgridRight);
    ret += sscanf(line, "TestGravityNumberOfParticles = %"ISYM,
		  &TestGravityNumberOfParticles);
    ret += sscanf(line, "TestGravityUseBaryons = %"ISYM,
		  &TestGravityUseBaryons);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "TestGravity")
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  } // end input from parameter file
 
  /* set up grid */
 
  if (TopGrid.GridData->TestGravityInitializeGrid(TestGravityDensity,
				                  TestGravityNumberOfParticles,
						  TestGravityUseBaryons
						  ) == FAIL){
    ENZO_FAIL("Error in TestGravityInitializeGrid.\n");
  }
 
  /* If requested, create a subgrid */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((TestGravitySubgridRight - TestGravitySubgridLeft)/
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
 
    /* compute the dimensions and left/right edges for the subgrid */
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*NumberOfGhostZones;
      LeftEdge[dim]    = TestGravitySubgridLeft;
      RightEdge[dim]   = TestGravitySubgridRight;
    }
 
    /* create a new subgrid and initialize it */
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->TestGravityInitializeGrid(TestGravityDensity*
	 POW(float(RefineBy), MetaData.TopGridRank), 0, TestGravityUseBaryons)
	== FAIL) {
      ENZO_FAIL("Error in TestGravityInitializeGrid.\n");
    }			
 
    /* Generate a static refine region. */
 
    StaticRefineRegionLevel[0] = 0;
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      StaticRefineRegionLeftEdge[0][dim] = TestGravitySubgridLeft;
      StaticRefineRegionRightEdge[0][dim] = TestGravitySubgridRight;
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
 
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  DataUnits[5] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "TestGravityDensity           = %"FSYM"\n",
	    TestGravityDensity);
    fprintf(Outfptr, "TestGravitySubgridLeft       = %"GOUTSYM"\n",
	    TestGravitySubgridLeft);
    fprintf(Outfptr, "TestGravitySubgridRight      = %"GOUTSYM"\n",
	    TestGravitySubgridRight);
    fprintf(Outfptr, "TestGravityNumberOfParticles = %"ISYM"\n",
	    TestGravityNumberOfParticles);
    fprintf(Outfptr, "TestGravityUseBaryons        = %"ISYM"\n\n",
	    TestGravityUseBaryons);
  }
 
  return SUCCESS;
 
}

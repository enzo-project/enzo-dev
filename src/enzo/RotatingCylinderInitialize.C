/***********************************************************************
/
/  Rotating Cylinder Test Problem
/
/  written by: Brian O'Shea
/  date:       February 2008
/  modified1:  
/
/  PURPOSE:
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
 
int RotatingCylinderInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData)
{
  if(debug){
    printf("Entering RotatingCylinderInitialize\n");
    fflush(stdout);
  }

  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *MetalName = "Metal_Density";

  /* parameter declarations */
 
  FLOAT RotatingCylinderSubgridLeft[MAX_DIMENSION], RotatingCylinderSubgridRight[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT RotatingCylinderCenterPosition[MAX_DIMENSION];

  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];
 
  /* make sure it is 3D */
 
  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do RotatingCylinder in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }
 
  for(i=0; i<MAX_DIMENSION; i++)
    RotatingCylinderCenterPosition[i] = 0.5;  // right in the middle of the box

  float RotatingCylinderVelocity[3]   = {0.0, 0.0, 0.0};   // gas initally at rest
  float RotatingCylinderBField[3]   = {0.0, 0.0, 0.0};   // gas initally at rest
  FLOAT RotatingCylinderRadius = 0.3;
  float RotatingCylinderLambda = 0.05;
  float RotatingCylinderOverdensity = 20.0;
  float RotatingCylinderDensity = 1.0;
  float RotatingCylinderTotalEnergy = 1.0;
  float Pi                      = 3.14159;

  /* set no subgrids by default. */
 
  RotatingCylinderSubgridLeft[0] = RotatingCylinderSubgridLeft[1] = 
    RotatingCylinderSubgridLeft[2] = 0.0;    // start of subgrid(s)

  RotatingCylinderSubgridRight[0] = RotatingCylinderSubgridRight[1] = 
    RotatingCylinderSubgridRight[2] = 0.0;    // end of subgrid(s)

  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters specifically for radiating shock problem*/

    ret += sscanf(line, "RotatingCylinderOverdensity  = %"FSYM, &RotatingCylinderOverdensity);
    ret += sscanf(line, "RotatingCylinderSubgridLeft = %"PSYM" %"PSYM" %"PSYM,
		  RotatingCylinderSubgridLeft,RotatingCylinderSubgridLeft+1,RotatingCylinderSubgridLeft+2);
    ret += sscanf(line, "RotatingCylinderSubgridRight = %"PSYM" %"PSYM" %"PSYM,
		  RotatingCylinderSubgridRight,RotatingCylinderSubgridRight+1,RotatingCylinderSubgridRight+2);
    ret += sscanf(line, "RotatingCylinderLambda = %"FSYM,
		        &RotatingCylinderLambda);

    ret += sscanf(line, "RotatingCylinderTotalEnergy = %"FSYM,
		        &RotatingCylinderTotalEnergy);

    ret += sscanf(line, "RotatingCylinderRadius = %"PSYM,
		        &RotatingCylinderRadius);
    ret += sscanf(line, "RotatingCylinderCenterPosition = %"PSYM" %"PSYM" %"PSYM,
		  RotatingCylinderCenterPosition, RotatingCylinderCenterPosition+1,
		  RotatingCylinderCenterPosition+2);

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && (strstr(line, "RotatingCylinder") || strstr(line, "TestProblem")) &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,
	 "*** warning: the following parameter line was not interpreted:\n%s\n",
	      line);
 
  } // end input from parameter file
 
 
  if (TopGrid.GridData->InitializeUniformGrid(RotatingCylinderDensity,
					      RotatingCylinderTotalEnergy,
					      RotatingCylinderTotalEnergy,
					      RotatingCylinderVelocity,
					      RotatingCylinderBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }
 
  /* Create as many subgrids as there are refinement levels
     needed to resolve the initial explosion region upon the start-up. */
 
  HierarchyEntry ** Subgrid;
  if (MaximumRefinementLevel > 0)
    Subgrid   = new HierarchyEntry*[MaximumRefinementLevel];
 
  /* Create new HierarchyEntries. */
 
  int lev;
  for (lev = 0; lev < MaximumRefinementLevel; lev++)
    Subgrid[lev] = new HierarchyEntry;
 
  for (lev = 0; lev < MaximumRefinementLevel; lev++) {
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      NumberOfSubgridZones[dim] =
	nint((RotatingCylinderSubgridRight[dim] - RotatingCylinderSubgridLeft[dim])/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *int(POW(RefineBy, lev + 1));
 
    if (debug)
      printf("RotatingCylinder:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
	     NumberOfSubgridZones[0]);
 
    if (NumberOfSubgridZones[0] > 0) {
 
      /* fill them out */
 
      if (lev == 0)
	TopGrid.NextGridNextLevel  = Subgrid[0];
      Subgrid[lev]->NextGridThisLevel = NULL;
      if (lev == MaximumRefinementLevel-1)
	Subgrid[lev]->NextGridNextLevel = NULL;
      else
	Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];
      if (lev == 0)
	Subgrid[lev]->ParentGrid        = &TopGrid;
      else
	Subgrid[lev]->ParentGrid        = Subgrid[lev-1];
 
      /* compute the dimensions and left/right edges for the subgrid */
 
      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*NumberOfGhostZones;
	LeftEdge[dim]    = RotatingCylinderSubgridLeft[dim];
	RightEdge[dim]   = RotatingCylinderSubgridRight[dim];
      }
 
      /* create a new subgrid and initialize it */
 
      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				     LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(RotatingCylinderDensity,
						   RotatingCylinderTotalEnergy,
						   RotatingCylinderTotalEnergy,
							RotatingCylinderVelocity,
							RotatingCylinderBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }
 
      /* set up the initial explosion area on the finest resolution subgrid */
 
      if (lev == MaximumRefinementLevel - 1)
	if (Subgrid[lev]->GridData->RotatingCylinderInitializeGrid(RotatingCylinderRadius,
								   RotatingCylinderCenterPosition,
								   RotatingCylinderLambda,
								   RotatingCylinderOverdensity) 
	    == FAIL) {
	  	  ENZO_FAIL("Error in RotatingCylinderInitialize[Sub]Grid.");
	}

    }
    else{
      printf("RotatingCylinder: single grid start-up.\n");
    }
  }

 
  /* set up subgrids from level 1 to max refinement level -1 */
 
  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
				       *(Subgrid[lev-1]->GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
 
  /* set up the root grid */
 
  if (MaximumRefinementLevel > 0) {
    if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
  }
  else
    if (TopGrid.GridData->RotatingCylinderInitializeGrid(RotatingCylinderRadius,
							 RotatingCylinderCenterPosition,
							 RotatingCylinderLambda,
							 RotatingCylinderOverdensity) == FAIL) {
            ENZO_FAIL("Error in RotatingCylinderInitializeGrid.");
    }

 
  /* set up field names and units -- NOTE: these absolutely MUST be in 
     the same order that they are in Grid_InitializeUniformGrids.C, or 
     else you'll find out that data gets written into incorrectly-named
     fields.  Just FYI. */

  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if(DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;

  if(MetaData.TopGridRank > 1)
    DataLabel[i++] = Vel2Name;

  if(MetaData.TopGridRank > 2)
    DataLabel[i++] = Vel3Name;

  if (TestProblemData.UseMetallicityField)
    DataLabel[i++] = MetalName;

  for(j=0; j < i; j++)
    DataUnits[j] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "RotatingCylinderOverdensity         = %"FSYM"\n"  , RotatingCylinderOverdensity);
    fprintf(Outfptr, "RotatingCylinderLambda         = %"FSYM"\n"  , RotatingCylinderLambda);
    fprintf(Outfptr, "RotatingCylinderTotalEnergy         = %"FSYM"\n"  , RotatingCylinderTotalEnergy);
    fprintf(Outfptr, "RotatingCylinderRadius         = %"PSYM"\n"  , RotatingCylinderRadius);
    fprintf(Outfptr, "RotatingCylinderCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
		  RotatingCylinderCenterPosition, RotatingCylinderCenterPosition+1,
		  RotatingCylinderCenterPosition+2);
    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

  } //   if (MyProcessorNumber == ROOT_PROCESSOR) 


  if(debug){

    printf("Exiting RotatingCylinderInitialize\n");
    fflush(stdout);
  }
 
  return SUCCESS;
 
}

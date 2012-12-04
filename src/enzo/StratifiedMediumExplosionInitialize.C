////////////////////////////////////////////////////////////////////////////////
//
//  Conduction Test Problem
//
//  written by: David A. Ventimiglia, Brian O'Shea
//  date:       June 2009
//  modified:  February 10, 2011 by BWO
//
//  PURPOSE: This initializes a bubble with user-controlled entropy in an 
//     ambient medium that is initially in hydrostatic equilibrium
//     (with density, temperature profile controlled by the user).  The
//     bubble then does its thing based on its entropy relative to the 
//     ambient medium.
//
//  RETURNS: SUCCESS or FAIL
//
////////////////////////////////////////////////////////////////////////////////
 
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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

// Problem Initializer
int StratifiedMediumExplosionInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData){

  if(debug){
    printf("Entering StratifiedMediumExplosionInitialize\n");
    fflush(stdout);
  }

  char line[MAX_LINE_LENGTH];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  int i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION], SubgridDims[MAX_DIMENSION];

  float StratifiedMediumExplosionDensity = 1.0;
  float StratifiedMediumExplosionTotalEnergy = 1.0;
  float StratifiedMediumExplosionGasEnergy = 1.0;
  float StratifiedMediumExplosionEnergy = 1.0;
  float StratifiedMediumExplosionVelocity[3] = {0.0,0.0,0.0};
  float StratifiedMediumExplosionInitialUniformBField[3] = {0.0,0.0,0.0};  // in Gauss

  FLOAT StratifiedMediumExplosionRadiusOfBubble = 0.1;  // units of box size
  int   StratifiedMediumExplosionPulseType = 1;  // pulse type
  FLOAT StratifiedMediumExplosionCenter[MAX_DIMENSION] = {0.5,0.5,0.5};

  FLOAT StratifiedMediumExplosionSubgridLeft[3]={0.0,0.0,0.0}; 
  FLOAT StratifiedMediumExplosionSubgridRight[3]={0.0,0.0,0.0};

  // Read parameters

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;
    ret += sscanf(line, "StratifiedMediumExplosionRadiusOfBubble = %"PSYM, &StratifiedMediumExplosionRadiusOfBubble);
    ret += sscanf(line, "StratifiedMediumExplosionEnergy = %"FSYM, &StratifiedMediumExplosionEnergy);
    ret += sscanf(line, "StratifiedMediumExplosionPulseType = %"ISYM, &StratifiedMediumExplosionPulseType);
    ret += sscanf(line, "StratifiedMediumExplosionCenter = %"PSYM" %"PSYM" %"PSYM, &StratifiedMediumExplosionCenter[0],
		  &StratifiedMediumExplosionCenter[1],&StratifiedMediumExplosionCenter[2]);
    ret += sscanf(line, "StratifiedMediumExplosionSubgridLeft = %"PSYM" %"PSYM" %"PSYM,
		  StratifiedMediumExplosionSubgridLeft,StratifiedMediumExplosionSubgridLeft+1,
		  StratifiedMediumExplosionSubgridLeft+2);
    ret += sscanf(line, "StratifiedMediumExplosionSubgridRight = %"PSYM" %"PSYM" %"PSYM,
		  StratifiedMediumExplosionSubgridRight,StratifiedMediumExplosionSubgridRight+1,
		  StratifiedMediumExplosionSubgridRight+2);

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    if (ret == 0 && 
	strstr(line, "=") && strstr(line, "StratifiedMediumExplosion") &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "*** warning: the following parameter line was not interpreted:\n%s\n", line);
    }
  }

  StratifiedMediumExplosionGasEnergy = StratifiedMediumExplosionTotalEnergy;

  float DensityUnits=1.0, LengthUnits=1.0, TemperatureUnits=1.0, TimeUnits=1.0,
    VelocityUnits=1.0;
  double MassUnits=1.0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, 0.0) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  // Create a uniform grid
  if (TopGrid.GridData->InitializeUniformGrid(StratifiedMediumExplosionDensity,
					      StratifiedMediumExplosionTotalEnergy,
					      StratifiedMediumExplosionGasEnergy,
					      StratifiedMediumExplosionVelocity,
					      StratifiedMediumExplosionInitialUniformBField) == FAIL) {
    ENZO_FAIL("Error in InitializeUniformGrid.");
  }

  // first pass puts in the stratified medium
  if (TopGrid.GridData->StratifiedMediumExplosionInitialize(StratifiedMediumExplosionRadiusOfBubble, 
							    StratifiedMediumExplosionPulseType,
							    StratifiedMediumExplosionEnergy,
							    StratifiedMediumExplosionCenter) == FAIL) {
    ENZO_FAIL("Error in StratifiedMediumExplosionInitialize.");
  }
  
  /*------------------------------------------------------*/


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
	nint((StratifiedMediumExplosionSubgridRight[dim] - StratifiedMediumExplosionSubgridLeft[dim])/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *int(POW(RefineBy, lev + 1));
 
    if (debug)
      printf("StratifiedMediumExplosion:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
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
	LeftEdge[dim]    = StratifiedMediumExplosionSubgridLeft[dim];
	RightEdge[dim]   = StratifiedMediumExplosionSubgridRight[dim];
      }
 
      /* create a new subgrid and initialize it */
 
      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				     LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(StratifiedMediumExplosionDensity,
					      StratifiedMediumExplosionTotalEnergy,
					      StratifiedMediumExplosionGasEnergy,
					      StratifiedMediumExplosionVelocity,
					      StratifiedMediumExplosionInitialUniformBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }
 
      /* set up the initial explosion area on the finest resolution subgrid */
 
      if (lev == MaximumRefinementLevel - 1)
	if (Subgrid[lev]->GridData->StratifiedMediumExplosionInitialize(StratifiedMediumExplosionRadiusOfBubble, 
							    StratifiedMediumExplosionPulseType,
							    StratifiedMediumExplosionEnergy,
							    StratifiedMediumExplosionCenter) 
	    == FAIL) {
	  	  ENZO_FAIL("Error in StratifiedMediumExplosionInitialize[Sub]Grid.");
	}

    }
    else{
      printf("StratifiedMediumExplosion: single grid start-up.\n");
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
    if (TopGrid.GridData->StratifiedMediumExplosionInitialize(StratifiedMediumExplosionRadiusOfBubble, 
							    StratifiedMediumExplosionPulseType,
							    StratifiedMediumExplosionEnergy,
							    StratifiedMediumExplosionCenter) == FAIL) {
            ENZO_FAIL("Error in StratifiedMediumExplosionInitializeGrid.");
    }


  /*------------------------------------------------------*/

  // Then perturb it

  // set up field names and units
  i = 0;
  DataLabel[i++] = "Density";
  DataLabel[i++] = "TotalEnergy";
  if (DualEnergyFormalism) {DataLabel[i++] = "GasEnergy";}
  if (MetaData.TopGridRank > 0) {DataLabel[i++] = "x-velocity";}
  if (MetaData.TopGridRank > 1 || HydroMethod > 2) {DataLabel[i++] = "y-velocity";}
  if (MetaData.TopGridRank > 2 || HydroMethod > 2) {DataLabel[i++] = "z-velocity";}

  if (TestProblemData.UseMetallicityField)
    DataLabel[i++] = "Metal_Density";

  for (j=0; j < i; j++) 
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "StratifiedMediumExplosionRadiusOfBubble = %"PSYM"\n", StratifiedMediumExplosionRadiusOfBubble);
    fprintf(Outfptr, "StratifiedMediumExplosionPulseType = %"ISYM"\n", StratifiedMediumExplosionPulseType);
    fprintf(Outfptr, "StratifiedMediumExplosionEnergy = %"FSYM"\n", StratifiedMediumExplosionEnergy);
    fprintf(Outfptr, "StratifiedMediumExplosionCenter = %"PSYM" %"PSYM" %"PSYM"\n", StratifiedMediumExplosionCenter,
		  StratifiedMediumExplosionCenter+1,StratifiedMediumExplosionCenter+2);
    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);
  }

  if(debug){
    printf("Exiting StratifiedMediumExplosionInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;

}

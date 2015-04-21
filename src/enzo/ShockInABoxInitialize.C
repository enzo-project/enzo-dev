/***********************************************************************
/
/  INITIALIZE A SHOCK IN A BOX
/
/  written by: Greg Bryan
/  date:       November, 1994
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
 
int ShockInABoxInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData, ExternalBoundary &Exterior)
{
  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  char *MachName   = "Mach";
  char *PSTempName = "PreShock_Temperature";
  char *PSDenName  = "PreShock_Density";
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
       SubgridDims[MAX_DIMENSION], ShockInABoxDirection;
  float ShockInABoxDensity[2], ShockInABoxPressure[2], ShockInABoxVelocity[2];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  /* set default parameters */
 
  FLOAT ShockInABoxBoundary = 0.5;     //
 
  float d1 = 1, v1 = 0, m = 2, p1 = 1, d2, p2, v2, c1, shockspeed = 0;
 
  d2 = d1*((Gamma+1)*m*m)/((Gamma-1)*m*m + 2);
  p2 = p1*(2.0*Gamma*m*m - (Gamma-1))/(Gamma+1);
  c1 = sqrt(Gamma*p1/d1);
  v2 = m*c1*(1-d1/d2);
 
  shockspeed = 0.9*c1 * m;
 
  ShockInABoxDirection   = 0;
  ShockInABoxDensity[0]  = d1;
  ShockInABoxVelocity[0] = shockspeed-v1;
  ShockInABoxPressure[0] = p1;
 
  ShockInABoxDensity[1]  = d2;
  ShockInABoxPressure[1] = p2;
  ShockInABoxVelocity[1] = shockspeed-v2;
 
  FLOAT ShockInABoxSubgridLeft  = 0.0;
  FLOAT ShockInABoxSubgridRight = 0.0;
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "ShockInABoxBoundary = %"PSYM, &ShockInABoxBoundary);
 
    ret += sscanf(line, "ShockInABoxLeftDensity = %"ESYM,
		  &ShockInABoxDensity[0]);
    ret += sscanf(line, "ShockInABoxLeftPressure = %"ESYM,
		  &ShockInABoxPressure[0]);
    ret += sscanf(line, "ShockInABoxLeftVelocity = %"ESYM,
		  &ShockInABoxVelocity[0]);
 
    ret += sscanf(line, "ShockInABoxRightDensity = %"ESYM,
		  &ShockInABoxDensity[1]);
    ret += sscanf(line, "ShockInABoxRightPressure = %"ESYM,
		  &ShockInABoxPressure[1]);
    ret += sscanf(line, "ShockInABoxRightVelocity = %"ESYM,
		  &ShockInABoxVelocity[1]);
 
    ret += sscanf(line, "ShockInABoxSubgridLeft = %"PSYM,
		  &ShockInABoxSubgridLeft);
    ret += sscanf(line, "ShockInABoxSubgridRight = %"PSYM,
		  &ShockInABoxSubgridRight);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "ShockInABox"))
      fprintf(stderr, "ShockInABox warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* set up grid */
 
  if( ShockInABoxDirection != 0 )
    ENZO_FAIL("Only ShockInABoxDirection=0 supported at the moment!");

  if (TopGrid.GridData->
      HydroShockTubesInitializeGrid(ShockInABoxBoundary,
				    ShockInABoxDensity[0], ShockInABoxDensity[1],
				    ShockInABoxVelocity[0], ShockInABoxVelocity[1],
				    0.0, 0.0,
				    0.0, 0.0, 
				    ShockInABoxPressure[0], ShockInABoxPressure[1]) == FAIL) {
    ENZO_FAIL("Error in HydroShockTubesInitializeGrid (called from ShockInABoxInitialize).\n");
  }
 
  /* If requested, create a subgrid */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((ShockInABoxSubgridRight - ShockInABoxSubgridLeft)/
	   ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	    float(MetaData.TopGridDims[dim])))
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
      LeftEdge[dim]    = ShockInABoxSubgridLeft;
      RightEdge[dim]   = ShockInABoxSubgridRight;
    }
 
    /* create a new subgrid and initialize it */
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->
	HydroShockTubesInitializeGrid(ShockInABoxBoundary,
				      ShockInABoxDensity[0], ShockInABoxDensity[1],
				      ShockInABoxVelocity[0], ShockInABoxVelocity[1],
				      0.0, 0.0,
				      0.0, 0.0, 
				      ShockInABoxPressure[0], ShockInABoxPressure[1]) == FAIL) {
      ENZO_FAIL("Error in HydroShockTubesInitializeGrid (called from ShockInABoxInitialize).\n");
    }
  }
  
 
  /* set up field names and units */
  int labelCounter = 0;

  // Density
  int DensNum = labelCounter;
  DataUnits[labelCounter] = NULL;
  DataLabel[labelCounter++] = DensName;
  // Velocity 1
  int Vel1Num = labelCounter;
  DataUnits[labelCounter] = NULL; 
  DataLabel[labelCounter++] = Vel1Name;
  // Velocity 2
  int Vel2Num = -1;
  int Vel3Num = -1;
  if (TopGrid.GridData->GetGridRank() > 1 || HydroMethod > 2) {
    Vel2Num = labelCounter;
    DataUnits[labelCounter] = NULL; 
    DataLabel[labelCounter++] = Vel2Name;
    // Velocity 3
    if (TopGrid.GridData->GetGridRank() > 2 || HydroMethod > 2) {
      Vel3Num = labelCounter;
      DataUnits[labelCounter] = NULL; 
      DataLabel[labelCounter++] = Vel3Name;
    }
  }
  // Total Energy
  int TENum = labelCounter;
  DataUnits[labelCounter] = NULL; 
  DataLabel[labelCounter++] = TEName;
 
  if (ShockMethod) {
    DataUnits[labelCounter] = NULL; 
    DataLabel[labelCounter++] = MachName;
    if (StorePreShockFields) {
      DataUnits[labelCounter] = NULL; 
      DataLabel[labelCounter++] = PSTempName;
      DataUnits[labelCounter] = NULL; 
      DataLabel[labelCounter++] = PSDenName;
    }
  } 

  /* Initialize the exterior. */
 
  Exterior.Prepare(TopGrid.GridData);
 
  float InflowValue[labelCounter], Dummy[labelCounter];
  for (int j = 0; j < labelCounter; j++){
    InflowValue[j] = 0.0;
    Dummy[j] = 0.0;
  }
  InflowValue[DensNum] = ShockInABoxDensity[0];
  InflowValue[Vel1Num] = ShockInABoxVelocity[0];
  InflowValue[TENum] = ShockInABoxPressure[0]/(Gamma-1.0)/ShockInABoxDensity[0]
                   + 0.5*POW(ShockInABoxVelocity[0], 2);
 
  if (Exterior.InitializeExternalBoundaryFace(0, inflow, outflow, InflowValue,
					      Dummy) == FAIL) {
      ENZO_FAIL("Error in InitializeExternalBoundaryFace.\n");
    }
 
  if (MetaData.TopGridRank > 1)
    Exterior.InitializeExternalBoundaryFace(1, reflecting, reflecting,
					    Dummy, Dummy);
  if (MetaData.TopGridRank > 2)
    Exterior.InitializeExternalBoundaryFace(2, reflecting, reflecting,
					    Dummy, Dummy);
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "ShockInABoxDirection     = %"ISYM"\n", ShockInABoxDirection);
    fprintf(Outfptr, "ShockInABoxBoundary      = %"GOUTSYM"\n\n",
	    ShockInABoxBoundary);
 
    fprintf(Outfptr, "ShockInABoxLeftDensity   = %"ESYM"\n", ShockInABoxDensity[0]);
    fprintf(Outfptr, "ShockInABoxLeftPressure  = %"ESYM"\n",
	    ShockInABoxPressure[0]);
    fprintf(Outfptr, "ShockInABoxLeftVelocity  = %"ESYM"\n\n",
	    ShockInABoxVelocity[0]);
 
    fprintf(Outfptr, "ShockInABoxRightDensity  = %"ESYM"\n", ShockInABoxDensity[1]);
    fprintf(Outfptr, "ShockInABoxRightPressure = %"ESYM"\n",
	    ShockInABoxPressure[1]);
    fprintf(Outfptr, "ShockInABoxRightVelocity = %"ESYM"\n\n",
	    ShockInABoxVelocity[1]);
  }
 
  return SUCCESS;
}

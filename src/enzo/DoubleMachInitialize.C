/***********************************************************************
/
/  INITIALIZE A DOUBLE-MACH REFLECTION TEST
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:
/
/  PURPOSE:
/    Initializes the double Mach reflection of a strong shock
/           test problem from Woodward and Collela.
/
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
 
int DoubleMachInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData, ExternalBoundary &Exterior)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
       SubgridDims[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  /* set default parameters */
 
  float d0 = 8.0, e0 = 291.25, u0 = 8.25*sqrt(3.0)/2.0, v0 = -8.25*0.5, w0 = 0;
  FLOAT DoubleMachSubgridLeft   = 0.0;    // start of subgrid
  FLOAT DoubleMachSubgridRight  = 0.0;    // end of subgrid
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "DoubleMachSubgridLeft = %"PSYM,
		  &DoubleMachSubgridLeft);
    ret += sscanf(line, "DoubleMachSubgridRight = %"PSYM,
		  &DoubleMachSubgridRight);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "DoubleMach"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  } // end input from parameter file
 
  /* set up grid */
 
  if (TopGrid.GridData->DoubleMachInitializeGrid(d0, e0, u0, v0, w0) == FAIL) {
    ENZO_FAIL("Error in DoubleMachInitializeGrid.\n");
  }
 
  /* If requested, create a subgrid */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((DoubleMachSubgridRight - DoubleMachSubgridLeft)/
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
      LeftEdge[dim]    = DoubleMachSubgridLeft;
      RightEdge[dim]   = DoubleMachSubgridRight;
    }
 
    /* create a new subgrid and initialize it */
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->DoubleMachInitializeGrid(d0, e0, u0, v0, w0)
	== FAIL) {
      ENZO_FAIL("Error in DoubleMachInitializeGrid (subgrid).\n");
    }			
  }
 
  /* Initialize the exterior. */
 
  Exterior.Prepare(TopGrid.GridData);
 
  float InflowValue[5], Dummy[5];
  InflowValue[0] = d0;
  InflowValue[1] = e0/d0 + 0.5*(u0*u0 + v0*v0 + w0*w0);
  InflowValue[2] = u0;
  InflowValue[3] = v0;
  InflowValue[4] = w0;
 
  if (Exterior.InitializeExternalBoundaryFace(0, inflow, outflow, InflowValue,
					      Dummy) == FAIL) {
      ENZO_FAIL("Error in InitializeExternalBoundaryFace.\n");
    }
 
  Exterior.InitializeExternalBoundaryFace(1, inflow, inflow, InflowValue,
					  InflowValue);
  if (MetaData.TopGridRank > 2)
    Exterior.InitializeExternalBoundaryFace(2, reflecting, reflecting,
					    Dummy, Dummy);
 
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

    fprintf(Outfptr, "DoubleMachSubgridLeft  = %"GOUTSYM"\n"  ,DoubleMachSubgridLeft);
    fprintf(Outfptr, "DoubleMachSubgridRight = %"GOUTSYM"\n\n",DoubleMachSubgridRight);
  }
 
  return SUCCESS;
 
}

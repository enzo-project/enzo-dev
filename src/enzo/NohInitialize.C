/***********************************************************************
/
/  INITIALIZE A NOH TEST PROBLEM
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:  Alexei Kritsuk, May 2005
/
/  PURPOSE:
/    Initializes the Noh Problem in 2D or 3D. 
/
/    Liska & Wendroff, 2003, SIAM J. Sci. Comput. 25, 995,
/    Section 4.5, Fig. 4.4.
/
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.

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

int NohInitialize(FILE *fptr, 
		  FILE *Outfptr, 
		  HierarchyEntry &TopGrid,
		  TopGridData &MetaData)
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

  /* The setup & parameters:
   *
   * There is an exact analytical solution for an ideal gas with
   * Gamma = 5/3.
   *
   *  The initial density = 1, pressure = 0, radial velocity = 1 (towards 
   * the center. Result -- circularly symmetric shock reflection from the
   * origin. The density behind the shock = 16, velocity = 0, pressure = 16/3.
   * The shock speed = 1/3. Ahead of the shock, at sqrt(x^2 + y^2) > t/3, the
   * density is (1 + t/sqrt(x^2 + y^2)), while the velocity and pressure are 
   * the same as initially.
   * 
   * In 3D spherical case the preshock conditions for density at R>t/3
   * (1 + t/sqrt(x^2 + y^2))^2. The density behind the shock is 64.
   * 
   * Domain [0,1]x[0,1]. Reflecting boundaries at x=y=z=0, exact solution at 
   * x=y=z=1.
   */

  /* Input parameters */

  /* Default parameters */

  NohProblemFullBox      = 0;
  float NohDensity       = 1.0000000000;
  float NohPressure      = 1.0000000000e-6;
  float NohVelocity      = - 1.0000000000;
  float NohSubgridLeft   = 0.0;    // start of subgrid
  float NohSubgridRight  = 0.0;    // end of subgrid

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "NohProblemFullBox = %"ISYM,
		  &NohProblemFullBox);
    ret += sscanf(line, "NohSubgridLeft = %"FSYM, 
		  &NohSubgridLeft);
    ret += sscanf(line, "NohSubgridRight = %"FSYM,
		  &NohSubgridRight);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "Noh"))
      fprintf(stderr, 
	  "warning: the following parameter line was not interpreted:\n%s\n", 
	      line);

  } // end input from parameter file

  /* set up grid */

  if (TopGrid.GridData->NohInitializeGrid(NohDensity,
					  NohPressure,
					  NohVelocity) == FAIL) {
    ENZO_FAIL("Error in NohInitializeGrid.\n");
  }

  /* If requested, create a subgrid */

  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] = 
      nint((NohSubgridRight - NohSubgridLeft)/
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
      LeftEdge[dim]    = NohSubgridLeft;
      RightEdge[dim]   = NohSubgridRight;
    }

    /* create a new subgrid and initialize it */

    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->NohInitializeGrid(NohDensity,
					  NohPressure,
					  NohVelocity)
	== FAIL) {
      ENZO_FAIL("Error in NohInitializeGrid (subgrid).\n");
    }			   
  }

  /* set up boundary types. */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = reflecting;
    MetaData.RightFaceBoundaryCondition[dim] = BoundaryUndefined;
  }

  if (NohProblemFullBox == 1)
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      MetaData.LeftFaceBoundaryCondition[dim]  = BoundaryUndefined;

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

    fprintf(Outfptr, "NohSubgridLeft  = %"GOUTSYM"\n"  , NohSubgridLeft);
    fprintf(Outfptr, "NohSubgridRight = %"GOUTSYM"\n\n", NohSubgridRight);
  }

  return SUCCESS;

}

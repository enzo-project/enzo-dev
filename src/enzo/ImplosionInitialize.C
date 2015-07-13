/***********************************************************************
/
/  INITIALIZE AN IMPLOSION SIMULATION
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:  Alexei Kritsuk, December 2004.
/
/  PURPOSE:
/    The implosion test sets up a converging shock problem in a square domain
/    (x,y) \in (0, 0.3)x(0, 0.3) with gas initially at rest. Initial
/    pressure and density is 1 everywhere except for a triangular region
/    (0.15,0)(0.15,0) where d=0.125 and p=0.14. Reflecting boundary conditions
/    at all boundaries. Adiabatic index gamma=1.4.
/
/    If AMR is used, a hierarchy of subgrids (one per level) will be generated
/    at start-up to properly resolve the initial discontinuity.
/
/         0.3
/           __________________________________
/           |                |               |
/           |                |    Domain     |
/           |           0.15 |               |
/           |               /|\              |
/           |              / | \             |
/           |             /  |  \            |
/           |            /  0|___\___________|
/           |            \   0   /0.15       |
/           |             \     /            |
/           |              \ Diamond         |
/           |               \ /              |
/           |                                |
/           |                                |
/           |                                |
/      -0.3 ---------------------------------- 0.3
/          -0.3
/
/
/   REFERENCE: Hui Li and Z. Li, JCP 153, 596, 1999.
/              Chang et al. JCP 160, 89, 1999.
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
#define DEFINE_STORAGE
#include "ImplosionGlobalData.h"
#undef DEFINE_STORAGE
 
int ImplosionInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* parameter declarations */
 
  FLOAT ImplosionSubgridLeft, ImplosionSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
                          SubgridDims[MAX_DIMENSION];
 
  /* set default parameters */
 
  float ImplosionVelocity[3]     = {0.0, 0.0, 0.0};   // gas initally at rest
  float ImplosionBField[3]      = {0.0, 0.0, 0.0};   // no magnetic field
  float ImplosionPressure        = 1.0;
  float ImplosionDiamondPressure = 0.14;
  ImplosionDensity               = 1.0;
  ImplosionDiamondDensity        = 0.125;
  ImplosionSubgridLeft           = 0.0;    // start of subgrid(s)
  ImplosionSubgridRight          = 0.0;    // end of subgrid(s)
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "ImplosionDensity  = %"FSYM, &ImplosionDensity);
    ret += sscanf(line, "ImplosionPressure = %"FSYM, &ImplosionPressure);
    ret += sscanf(line, "ImplosionDiamondDensity  = %"FSYM,
		        &ImplosionDiamondDensity);
    ret += sscanf(line, "ImplosionDiamondPressure = %"FSYM,
		        &ImplosionDiamondPressure);
    ret += sscanf(line, "ImplosionSubgridLeft = %"FSYM,
		        &ImplosionSubgridLeft);
    ret += sscanf(line, "ImplosionSubgridRight = %"FSYM,
		        &ImplosionSubgridRight);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "Implosion") &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,
	 "warning: the following parameter line was not interpreted:\n%s\n",
	      line);
 
  } // end input from parameter file
 
 
  /* Compute total energies */
 
  ImplosionTotalEnergy = ImplosionPressure/((Gamma - 1.0)*ImplosionDensity);
  ImplosionDiamondTotalEnergy = ImplosionDiamondPressure/((Gamma - 1.0)*
						   ImplosionDiamondDensity);
 
  /* set the reflecting boundaries */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = reflecting;
    MetaData.RightFaceBoundaryCondition[dim] = reflecting;
  }
 
  /* set up uniform grid without a "diamond" */
 
  if (TopGrid.GridData->InitializeUniformGrid(ImplosionDensity,
					      ImplosionTotalEnergy,
					      ImplosionTotalEnergy,
					      ImplosionVelocity, ImplosionBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }
 
  /* set up the diamond */
 
  if (TopGrid.GridData->ImplosionInitializeGrid(ImplosionDiamondDensity,
					      ImplosionDiamondTotalEnergy)
      == FAIL) {
        ENZO_FAIL("Error in ImplosionInitializeGrid.");
  }
 
 
  /* Create as many subgrids as refinement levels to resolve
     the initial discontinuity upon the start-up.            */
 
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
	nint((ImplosionSubgridRight - ImplosionSubgridLeft)/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *POW(RefineBy, lev + 1);
 
    if (debug)
      printf("Implosion:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
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
	LeftEdge[dim]    = ImplosionSubgridLeft;
	RightEdge[dim]   = ImplosionSubgridRight;
      }
 
      /* create a new subgrid and initialize it */
 
      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				     LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(ImplosionDensity,
						   ImplosionTotalEnergy,
						   ImplosionTotalEnergy,
					        ImplosionVelocity, ImplosionBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }
 
      /* set up the diamond */
 
      if (Subgrid[lev]->GridData->ImplosionInitializeGrid(
				  ImplosionDiamondDensity,
				  ImplosionDiamondTotalEnergy)
	  == FAIL) {
		ENZO_FAIL("Error in ImplosionInitialize[Sub]Grid.");
      }
    }
    else
      printf("Implosion: single grid start-up.\n");
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
    fprintf(Outfptr, "ImplosionDensity         = %"FSYM"\n"  , ImplosionDensity);
    fprintf(Outfptr, "ImplosionPressure        = %"FSYM"\n"  , ImplosionPressure);
    fprintf(Outfptr, "ImplosionDiamondDensity  = %"FSYM"\n",
	    ImplosionDiamondDensity);
    fprintf(Outfptr, "ImplosionDiamondPressure = %"FSYM"\n",
	    ImplosionDiamondPressure);
  }
 
  return SUCCESS;
 
}

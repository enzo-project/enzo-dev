/***********************************************************************
/
/  INITIALIZE PROTOSTELLAR CORE COLLAPSE
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:  Alexei Kritsuk, June 2005. 
/  modified2:  Robert Harkness, August 12th 2006 (for 64b ints)
/
/  PURPOSE:
/
/  REFERENCE: Bate 1998, ApJL 508, L95-L98
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
#define DEFINE_STORAGE
#undef DEFINE_STORAGE

int ProtostellarCollapseInitialize(FILE *fptr, FILE *Outfptr, 
				   HierarchyEntry &TopGrid,
				   TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *GPotName = "Grav_Potential";

  /* parameter declarations */

  FLOAT ProtostellarCollapseSubgridLeft, ProtostellarCollapseSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  
  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
                          SubgridDims[MAX_DIMENSION],
                                  xyz[MAX_DIMENSION];
  float TopCell[MAX_DIMENSION];

  /* make sure it is 3D */

  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do ProtostellarCollapse in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }    

  /* Setup and parameters:

   1. ambient density (should be very small) - free parameter
   2. core density = 10^6
   3. core radius - free parameter
   4. core angular velocity - free parameter
   5. Gamma as a function of density. 
      Gamma = 1.001 within the initial core radius
      Gamma = 5/3 outside the initial core radius
   6. pressure = 1 everywhere at t=0.
   7. ambient gas velocity = 0.

  */

  float ProtostellarCollapseVelocity[3]     = {0.0, 0.0, 0.0}; // ambient gas initally at rest
  float ProtostellarCollapseBField[3]     = {0.0, 0.0, 0.0}; // ambient gas initally at rest
  float ProtostellarCollapseCoreDensity     = 2000.0; // 10^6/500
  float ProtostellarCollapseCoreEnergy      = 1e3;  // thermal energy, assumes Gamma = 1.001, p=1, d=1
  float ProtostellarCollapseOuterDensity    = 1.0;
  float ProtostellarCollapseOuterEnergy     = 1e3;
  float ProtostellarCollapseCoreRadius      = 0.005;
  float ProtostellarCollapseAngularVelocity = 0.0;
  float dx = (DomainRightEdge[0] - DomainLeftEdge[0])/
                                                   MetaData.TopGridDims[0];

  /* set no subgrids by default. */

  ProtostellarCollapseSubgridLeft           = 0.0;    // start of subgrid(s)
  ProtostellarCollapseSubgridRight          = 0.0;    // end of subgrid(s)

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "ProtostellarCollapseCoreRadius      = %"FSYM, 
		  &ProtostellarCollapseCoreRadius);
    ret += sscanf(line, "ProtostellarCollapseOuterDensity    = %"FSYM, 
		  &ProtostellarCollapseOuterDensity);
    ret += sscanf(line, "ProtostellarCollapseAngularVelocity = %"FSYM, 
		  &ProtostellarCollapseAngularVelocity);
    ret += sscanf(line, "ProtostellarCollapseSubgridLeft     = %"FSYM, 
		  &ProtostellarCollapseSubgridLeft);
    ret += sscanf(line, "ProtostellarCollapseSubgridRight    = %"FSYM, 
		  &ProtostellarCollapseSubgridRight);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "ProtostellarCollapse") && 
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, 
	 "warning: the following parameter line was not interpreted:\n%s\n", 
	      line);

  } // end input from parameter file


  /* set the periodic boundaries */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }

  /* set up uniform grid without the core */

  if (TopGrid.GridData->InitializeUniformGrid(ProtostellarCollapseOuterDensity, 
					      ProtostellarCollapseOuterEnergy,
					      ProtostellarCollapseOuterEnergy,
					      ProtostellarCollapseVelocity,
                          ProtostellarCollapseBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }

  /* Create as many subgrids as there are refinement levels 
     needed to resolve the initial core upon the start-up. 

     Synchronize Left and Right edges of the refined region to match 
     cell boundaries on the top level and the number of procs. 
     XYZ symmetry is assumed nearly everywhere.
  */

  int startWithSubgrids = FALSE;
  if (ProtostellarCollapseSubgridRight > ProtostellarCollapseSubgridLeft) {
    startWithSubgrids = TRUE;
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      TopCell[dim] = (DomainRightEdge[dim] - DomainLeftEdge[dim])/
	(float(MetaData.TopGridDims[dim]));

    int nzones = nint( (ProtostellarCollapseSubgridRight - ProtostellarCollapseSubgridLeft)/
		       TopCell[0]/2.0 );
    ProtostellarCollapseSubgridLeft  = (DomainRightEdge[0] + DomainLeftEdge[0])/2.0 
      - TopCell[0]*nzones;
    ProtostellarCollapseSubgridRight = (DomainRightEdge[0] + DomainLeftEdge[0])/2.0 
      + TopCell[0]*nzones;
  }

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
	nint((ProtostellarCollapseSubgridRight - ProtostellarCollapseSubgridLeft)/
	     TopCell[dim])*POW(RefineBy, lev + 1);

    if (debug)
      printf("PCI: Level[%"ISYM"]: NumberOfSubgridZones = [%"ISYM",%"ISYM",%"ISYM"]\n", lev+1, 
	     NumberOfSubgridZones[0], NumberOfSubgridZones[1], NumberOfSubgridZones[2]);

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
        LeftEdge[dim]    = ProtostellarCollapseSubgridLeft;
        RightEdge[dim]   = ProtostellarCollapseSubgridRight;
      }

      if (debug)
	printf("PCI: Level[%"ISYM"]: SubgridDims = [%"ISYM",%"ISYM",%"ISYM"] [%"FSYM"-%"FSYM", %"FSYM"-%"FSYM", %"FSYM"-%"FSYM"]\n", lev+1, 
	       SubgridDims[0], SubgridDims[1], SubgridDims[2],
	       LeftEdge[0],RightEdge[0],LeftEdge[1],RightEdge[1],LeftEdge[2],RightEdge[2]);

      /* create a new subgrid and initialize it */
	
      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
					    LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(ProtostellarCollapseOuterDensity,
							ProtostellarCollapseOuterEnergy,
							ProtostellarCollapseOuterEnergy,
							ProtostellarCollapseVelocity,
                            ProtostellarCollapseBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }

      /* set up the dense core on the finest resolution subgrid */
      
      if (lev == MaximumRefinementLevel - 1)
	if (Subgrid[lev]->GridData->ProtostellarCollapseInitializeGrid(
				    ProtostellarCollapseCoreDensity,
				    ProtostellarCollapseCoreEnergy,
				    ProtostellarCollapseCoreRadius,
				    ProtostellarCollapseAngularVelocity) 
	    == FAIL) {
	  	  ENZO_FAIL("Error in ProtostellarCollapseInitialize[Sub]Grid.");
	}
    }
    else
      printf("ProtostellarCollapse: single grid start-up.\n");
  }


  /* set up subgrids from level 1 to max refinement level -1 */

  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
      if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
			     *(Subgrid[lev-1]->GridData))
	  == FAIL) {
		ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
      }
  
  /* set up the root grid */

  if (MaximumRefinementLevel > 0 && startWithSubgrids)
    if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }

  else
    if (TopGrid.GridData->ProtostellarCollapseInitializeGrid(
                          ProtostellarCollapseCoreDensity,
			  ProtostellarCollapseCoreEnergy,
			  ProtostellarCollapseCoreRadius,
			  ProtostellarCollapseAngularVelocity) == FAIL) {
            ENZO_FAIL("Error in ProtostellarCollapseInitializeGrid.");
    }

  /* set up field names and units */

  DataLabel[0] = DensName;
  DataLabel[1] = TEName;
  DataLabel[2] = Vel1Name;
  DataLabel[3] = Vel2Name;
  DataLabel[4] = Vel3Name;
  //  if (WritePotential)   //DR: not using this of Alexei's for now
  //    DataLabel[5] = GPotName;

  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  //  if (WritePotential)   //DR: not using this of Alexei's for now
  //    DataUnits[5] = NULL;


  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "ProtostellarCollapseCoreDensity     = %"FSYM"\n"  , 
	              ProtostellarCollapseCoreDensity);
    fprintf(Outfptr, "ProtostellarCollapseCoreEnergy      = %"FSYM"\n"  , 
	              ProtostellarCollapseCoreEnergy);
    fprintf(Outfptr, "ProtostellarCollapseOuterDensity    = %"FSYM"\n"  , 
	              ProtostellarCollapseOuterDensity);
    fprintf(Outfptr, "ProtostellarCollapseCoreRadius      = %"FSYM"\n"  , 
	              ProtostellarCollapseCoreRadius);
    fprintf(Outfptr, "ProtostellarCollapseAngularVelocity = %"FSYM"\n"  , 
	              ProtostellarCollapseAngularVelocity);
  }

  return SUCCESS;

}

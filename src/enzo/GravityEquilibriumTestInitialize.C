/***********************************************************************
/
/  INITIALIZE PLANE PARALLEL ATMOSPHERE TO TEST GRAVITY-PRESSURE EQUILIBRIUM
/
/  written by: Greg Bryan
/  date:       June, 1996
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
 
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
 
 
int GravityEquilibriumTestInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char  line[MAX_LINE_LENGTH];
  int   ret;
 
  /* Error check. */
 
  if (!UniformGravity)
    fprintf(stderr, "warning: UniformGravity is off, should be on");
 
  /* set default parameters */
 
  float GravityEquilibriumTestScaleHeight = 0.1;
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "GravityEquilibriumTestScaleHeight = %"FSYM,
		  &GravityEquilibriumTestScaleHeight);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") &&
	strstr(line, "GravityEquilibriumTest") && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  } // end input from parameter file
 
  /* Set up grid. */
 
  if (TopGrid.GridData->GravityEquilibriumTestInitializeGrid(
				GravityEquilibriumTestScaleHeight
						  ) == FAIL){
    ENZO_FAIL("Error in GravityEquilibriumTestInitializeGrid.\n");
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
 
  if (MyProcessorNumber == ROOT_PROCESSOR)

    fprintf(Outfptr, "GravityEquilibriumTestScaleHeight = %"GSYM"\n",
	    GravityEquilibriumTestScaleHeight);
 
  return SUCCESS;
 
}

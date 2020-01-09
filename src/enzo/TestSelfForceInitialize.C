/***********************************************************************
/
/  INITIALIZE A SELFFORCE TEST
/
/  written by: Jean-Claude Passy
/  date:       June 2013
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
 
int TestSelfForceInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
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
 
  float TestSelfForceDensity           = 1.0;  // density of central peak
 
  FLOAT TestSelfForcePartciclePositionX = 0.5;
  FLOAT TestSelfForcePartciclePositionY = 0.5;
  FLOAT TestSelfForcePartciclePositionZ = 0.5;

  FLOAT TestSelfForcePartcicleVelocityX = 0.0;
  FLOAT TestSelfForcePartcicleVelocityY = 0.0;
  FLOAT TestSelfForcePartcicleVelocityZ = 0.0;

  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "TestSelfForceDensity = %"FSYM, &TestSelfForceDensity);
    ret += sscanf(line, "TestSelfForcePartciclePosition = %"PSYM" %"PSYM" %"PSYM,
                  &TestSelfForcePartciclePositionX,
                  &TestSelfForcePartciclePositionY,
                  &TestSelfForcePartciclePositionZ);
    ret += sscanf(line, "TestSelfForcePartcicleVelocity  = %"PSYM" %"PSYM" %"PSYM,
                  &TestSelfForcePartcicleVelocityX,
                  &TestSelfForcePartcicleVelocityY,
                  &TestSelfForcePartcicleVelocityZ);

    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "TestSelfForce")
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  } // end input from parameter file
 
  /* set up grid */
 
  if (TopGrid.GridData->TestSelfForceInitializeGrid(TestSelfForceDensity,1,
						    TestSelfForcePartciclePositionX,
						    TestSelfForcePartciclePositionY,
						    TestSelfForcePartciclePositionZ,
						    TestSelfForcePartcicleVelocityX,
						    TestSelfForcePartcicleVelocityY,
						    TestSelfForcePartcicleVelocityZ
						    ) == FAIL)
    ENZO_FAIL("Error in TestSelfForceInitializeGrid.\n");
  
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

    fprintf(Outfptr, "TestSelfForceDensity           = %"FSYM"\n",
	    TestSelfForceDensity);
  }
 
  return SUCCESS;
 
}

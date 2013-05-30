/***********************************************************************
/
/  INITIALIZE A STAR PARTICLE TEST
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/    Initialize a star particle test in a uniform medium.
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

int TestStarParticleInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			       TopGridData &MetaData,float *Initialdt)
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

  /* Error check. */


  /* set default parameters */

  float TestStarParticleDensity     = 1.0;
  float TestStarParticleEnergy      = 1.0;
  float TestStarParticleVelocity[3] = {0.0, 0.0, 0.0};
  FLOAT TestStarParticleStarVelocity[3] = {0.0, 0.0, 0.0};
  FLOAT TestStarParticleStarPosition[3] = {0.5, 0.5, 0.5};
  float TestStarParticleBField[3]   = {0.0, 0.0, 0.0};
  float TestStarParticleStarMass    = 100.0;

  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "TestStarParticleDensity = %"FSYM,
		  &TestStarParticleDensity);
    ret += sscanf(line, "TestStarParticleEnergy = %"FSYM,
		  &TestStarParticleEnergy);
    ret += sscanf(line, "TestStarParticleStarMass = %"FSYM,
		  &TestStarParticleStarMass);
    ret += sscanf(line,"TestStarParticleStarVelocity = %"PSYM" %"PSYM" %"PSYM, 
		  &TestStarParticleStarVelocity[0],
		  &TestStarParticleStarVelocity[1],
		  &TestStarParticleStarVelocity[2]);
    ret += sscanf(line,"TestStarParticleStarPosition = %"PSYM" %"PSYM" %"PSYM, 
		  &TestStarParticleStarPosition[0],
		  &TestStarParticleStarPosition[1],
		  &TestStarParticleStarPosition[2]);
    
    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "TestStarParticle") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up uniform grid as of before explosion */
  
  if (TopGrid.GridData->InitializeUniformGrid(TestStarParticleDensity, 
					      TestStarParticleEnergy,
					      TestStarParticleEnergy,
					      TestStarParticleVelocity,
					      TestStarParticleBField) == FAIL)
    ENZO_FAIL("Error in InitializeUniformGrid.");
 
  if (TopGrid.GridData->
      TestStarParticleInitializeGrid(TestStarParticleStarMass,
				     Initialdt, 
				     TestStarParticleStarVelocity,
				     TestStarParticleStarPosition) == FAIL)
    ENZO_FAIL("Error in TestStarParticleInitializeGrid.\n");

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

    fprintf(Outfptr, "TestStarParticleDensity = %"FSYM"\n",
	    TestStarParticleDensity);
    fprintf(Outfptr, "TestStarParticleEnergy = %"FSYM"\n",
	    TestStarParticleEnergy);
  }

  return SUCCESS;

}


/***********************************************************************
/
/  INITIALIZE A GRAVITY TEST
/
/  written by: Greg Bryan
/  date:       August, 2005
/  modified1:
/
/  PURPOSE:
/    Initialize a particle orbit test.  Create two particles, one in the
/     center of the grid and the other a much less massive test particle
/     in order to test how well the gravity and integration works.
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

int TestOrbitInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
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

  int   TestOrbitNumberOfParticles = 1;    // number of test particles
  FLOAT TestOrbitRadius            = 0.2;  
  float TestOrbitCentralMass       = 1.0;
  float TestOrbitTestMass          = 1.0e-6;
  int   TestOrbitUseBaryons        = FALSE; // not implemented 

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "TestOrbitNumberOfParticles = %"ISYM,
		  &TestOrbitNumberOfParticles);
    ret += sscanf(line, "TestOrbitRadius = %"PSYM,
		  &TestOrbitRadius);
    ret += sscanf(line, "TestOrbitCentralMass = %"FSYM,
		  &TestOrbitCentralMass);
    ret += sscanf(line, "TestOrbitTestMass = %"FSYM,
		  &TestOrbitTestMass);
    ret += sscanf(line, "TestOrbitUseBaryons = %"ISYM,
		  &TestOrbitUseBaryons);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "TestOrbit") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up grid */

  if (TopGrid.GridData->TestOrbitInitializeGrid(TestOrbitNumberOfParticles,
						TestOrbitRadius,
						TestOrbitCentralMass,
						TestOrbitTestMass,
						TestOrbitUseBaryons
						  ) == FAIL){
    ENZO_FAIL("Error in TestOrbitInitializeGrid.\n");
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
  
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "TestOrbitNumberOfParticles = %"ISYM"\n",
	    TestOrbitNumberOfParticles);
    fprintf(Outfptr, "TestOrbitRadius            = %"GOUTSYM"\n",
	    TestOrbitRadius);
    fprintf(Outfptr, "TestOrbitCentralMass       = %"GOUTSYM"\n",
	    TestOrbitCentralMass);
    fprintf(Outfptr, "TestOrbitTestMass          = %"GOUTSYM"\n",
	    TestOrbitTestMass);
    fprintf(Outfptr, "TestOrbitUseBaryons        = %"ISYM"\n\n",
	    TestOrbitUseBaryons);
  }

  return SUCCESS;

}


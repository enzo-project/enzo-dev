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
  char *GPotName = "Grav_Potential";

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
  int   TestOrbitUseBaryons        = FALSE; // not on by default

  // Quantities for the baryon field if we turn it on (necessary for the potential)
  float TestOrbitBackgroundVelocity[3]   = {0.0, 0.0, 0.0};   // gas initally at rest
  float TestOrbitBackgroundDensity;
  float TestOrbitBackgroundEnergy = 1.0;
  float ZeroBField[3] = {0.0, 0.0, 0.0};

  /* The test orbit background density needs to be small, otherwise bad things will 
     happen to the orbit.  The reason for this is that the test orbit requires 
     non-periodic (isolated) gravity boundary conditions, and it's hard to maintain
     a circular orbit inside of a cube of gas.  BUT, we need to have baryon fields
     enabled in order to write out the gravitational potential, so we compromise by
     setting this to the smallest practical number.  As long as 
     (bayron density) * (simulation volume) is much smaller than the mass of 
     the big particle at the center of the volume, everything ought to be fine.
     --BWO, May 2013 */

  TestOrbitBackgroundDensity = tiny_number;   
  

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


  if(TestOrbitUseBaryons == TRUE && UseHydro == TRUE){
    UseHydro = FALSE;
    fprintf(stderr,"TestOrbit: UseBaryons = TRUE, turning off hydro!\n");
  }

  /* set up uniform top grid if we're using baryon fields (this
     must be done if we want to write out the gravitational potential). */
  if(TestOrbitUseBaryons == TRUE)
    if (TopGrid.GridData->InitializeUniformGrid(TestOrbitBackgroundDensity,
						TestOrbitBackgroundEnergy,
						TestOrbitBackgroundEnergy,
						TestOrbitBackgroundVelocity,
						ZeroBField
						) == FAIL) {
      ENZO_FAIL("Error in InitializeUniformGrid.");
    }

  /* put particles in the grid - doesn't need grid baryon quantities initialized */
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

  if (WritePotential)
    DataLabel[count++] = GPotName;  

  for (int j = 0; j < count; j++)
    DataUnits[j] = NULL;

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


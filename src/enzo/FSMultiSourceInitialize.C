/***********************************************************************
/
/  INITIALIZE FREE-STREAMING RADIATION TEST -- INITIALLY ZERO FIELD 
/  WITH MULTIPLE RANDOMLY-LOCATED SOURCES IN EACH SUBDOMAIN
/
/  written by: Daniel Reynolds
/  date:       June 2009
/
/  PURPOSE:
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


/* default constants */
#define DEFAULT_MU 0.6       // mean molecular mass
#define MIN_TEMP 1.0         // minimum temperature [K]


// function prototypes
int InitializeRateData(FLOAT Time);




int FSMultiSourceInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			    TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stdout,"Entering FSMultiSourceInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "TotalEnergy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName   = "FS_Radiation";

  // local declarations
  char line[MAX_LINE_LENGTH];
  int  dim, ret;

  // Setup and parameters:
  float Density              = 10.0;
  float X0Velocity           = 0.0;
  float X1Velocity           = 0.0;
  float X2Velocity           = 0.0;
  float TEnergy              = 1.0;
  float RadiationEnergy      = 10.0;

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
	ret += sscanf(line, "FSProbVelocity = %"FSYM" %"FSYM" %"FSYM,
		      &X0Velocity, &X1Velocity, &X2Velocity);
	ret += sscanf(line, "FSProbDensity = %"FSYM, &Density);
	ret += sscanf(line, "FSProbTEnergy = %"FSYM, &TEnergy);
	ret += sscanf(line, "FSProbRadiationEnergy = %"FSYM, 
		      &RadiationEnergy);
      } // end input from parameter file
      fclose(RHfptr);
    }
  }


  /* error checking */
  if (Mu != DEFAULT_MU) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: mu =%f assumed in initialization; setting Mu = %f for consistency.\n", DEFAULT_MU);
    Mu = DEFAULT_MU;
  }

  // set up the grid(s) on this level
  if (debug)
    printf("FSMultiSourceInitialize: calling grid initializer\n");
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->FSMultiSourceInitializeGrid(Density, X0Velocity, 
			X1Velocity, X2Velocity, TEnergy, 
			RadiationEnergy, local) == FAIL) {
      fprintf(stderr, "Error in Grid_FSMultiSourceInitializeGrid.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }

  // set up field names and units
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  DataLabel[BaryonField++] = Vel0Name;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  DataLabel[BaryonField++] = RadName;

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;
 
#endif

}

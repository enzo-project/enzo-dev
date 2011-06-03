/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- EVOLVE A CONSTANT FIELD
/
/  written by: Daniel Reynolds
/  date:       November, 2006
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


// function prototypes
int InitializeRateData(FLOAT Time);




int RadHydroPulseTestInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stdout,"Entering RadHydroPulseTestInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "TotalEnergy";
  char *IEName    = "Internal_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName   = "Grey_Radiation_Energy";

  // local declarations
  int dim;

  // make sure it is 3D
  if (MetaData.TopGridRank != 3) {
    printf("Cannot do Rad-Hydro Tests in %"ISYM" dimension(s)\n", 
	   MetaData.TopGridRank);
    return FAIL;
  }    


  // Setup and parameters:
  //  1. ambient density (should be very small) - free parameter
  //  2. ambient radiation energy
  //  3. initial time step size
  //  4. dimension for streaming test problem {0,1,2}
  //  5. direction for streaming radiation {0:l->r, 1:r->l}
  float RadHydroDensity   = 1.0;
  float RadHydroRadEnergy = 1.0e-10;
  int RadPulseDim         = 0;

  // overwrite parameters from RadHydroParamFile file, if it exists
  char line[MAX_LINE_LENGTH];
  int  ret;
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
	ret += sscanf(line, "RadHydroDensity = %"FSYM, &RadHydroDensity);
	ret += sscanf(line, "RadHydroRadEnergy = %"FSYM, &RadHydroRadEnergy);
	ret += sscanf(line, "RadPulseDim = %"ISYM, &RadPulseDim);
      } // end input from parameter file
      fclose(RHfptr);
    }
  }


  /* error checking */
  if (Mu != 0.6) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: mu = 0.6 assumed in initialization; setting Mu = 0.6 for consistency.\n");
    Mu = 0.6;
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }

  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->RadHydroPulseTestInitializeGrid(RadHydroDensity, 
			RadHydroRadEnergy, RadPulseDim, local) == FAIL) {
      fprintf(stderr, "Error in RadHydroPulseTestInitializeGrid.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }

  // set up field names and units
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  if (DualEnergyFormalism) 
    DataLabel[BaryonField++] = IEName;
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


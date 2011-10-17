/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- Grey Su-Olson Marshak Wave
/
/  written by: Daniel Reynolds and John Hayes
/  date:       June, 2007
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




int RadHydroGreyMarshakWaveInitialize(FILE *fptr, FILE *Outfptr, 
				      HierarchyEntry &TopGrid,
				      TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stdout,"Entering RadHydroGreyMarshakWaveInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "TotalEnergy";
  char *IEName    = "Internal_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName   = "Grey_Radiation_Energy";

  // Setup and parameters:
  //  1. ambient density (should be very small) - free parameter
  //  2. ambient radiation energy
  //  3. initial time step size
  //  4. propagation coordinate for marshak problem {0,1,2}
  float RadHydroDensity   = 1.0;
  float RadHydroGasEnergy = 1.0;
  float RadHydroRadEnergy = 1.0;
  int   GreyMarshDir      = 0;

  // overwrite parameters from RadHydroParamFile file, if it exists
  char line[MAX_LINE_LENGTH];
  int  dim, ret;
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
	ret += sscanf(line, "RadHydroDensity = %"FSYM, &RadHydroDensity);
	ret += sscanf(line, "RadHydroGasEnergy = %"FSYM, &RadHydroGasEnergy);
	ret += sscanf(line, "RadHydroRadEnergy = %"FSYM, &RadHydroRadEnergy);
	ret += sscanf(line, "GreyMarshDir = %"ISYM, &GreyMarshDir);
      } // end input from parameter file
      fclose(RHfptr);
    }
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
    if (Temp->GridData->RadHydroGreyMarshakWaveInitializeGrid(RadHydroDensity, 
	  	        RadHydroGasEnergy, RadHydroRadEnergy, GreyMarshDir, 
			local) == FAIL) {
      fprintf(stderr, "Error in RadHydroGreyMarshakWaveInitializeGrid.\n");
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


/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- IONIZATION TEST IN A 
/  R^{-2} DENSITY PROFILE
/
/  written by: Daniel Reynolds
/  date:       December 2007
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




int RHIonizationSteepInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stdout,"Entering RHIonizationSteepInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "TotalEnergy";
  char *IEName    = "Internal_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName   = "Grey_Radiation_Energy";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *DeName    = "Electron_Density";

  // local declarations
  char line[MAX_LINE_LENGTH];
  int  dim, ret;

  // Setup and parameters:
  float RadHydroX0Velocity           = 0.0;
  float RadHydroX1Velocity           = 0.0;
  float RadHydroX2Velocity           = 0.0;
  float RadHydroNumDensity           = 3.2;           // [cm^{-3}]
  float RadHydroDensityRadius        = 2.8234155e+20; // 91.5 pc [cm]
  float DensityCenter0               = 0.0;
  float DensityCenter1               = 0.0;
  float DensityCenter2               = 0.0;
  float RadHydroTemperature          = 100.0;         // [K]
  float RadHydroRadiationEnergy      = 1.0e-20;
  float RadHydroInitialFractionHII   = 0.0;
  int   RadHydroChemistry            = 1;
  int   RadHydroModel                = 1;

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
	ret += sscanf(line, "RadHydroVelocity = %"FSYM" %"FSYM" %"FSYM,
		      &RadHydroX0Velocity, &RadHydroX1Velocity, 
		      &RadHydroX2Velocity);
	ret += sscanf(line, "RadHydroChemistry = %"ISYM, 
		      &RadHydroChemistry);
	ret += sscanf(line, "RadHydroModel = %"ISYM, 
		      &RadHydroModel);
	ret += sscanf(line, "RadHydroNumDensity = %"FSYM, 
		      &RadHydroNumDensity);
	ret += sscanf(line, "RadHydroDensityRadius = %"FSYM, 
		      &RadHydroDensityRadius);
	ret += sscanf(line, "RadHydroTemperature = %"FSYM, 
		      &RadHydroTemperature);
	ret += sscanf(line, "RadHydroRadiationEnergy = %"FSYM, 
		      &RadHydroRadiationEnergy);
	ret += sscanf(line, "RadHydroInitialFractionHII = %"FSYM, 
		      &RadHydroInitialFractionHII);
	ret += sscanf(line, "EtaCenter = %"FSYM" %"FSYM" %"FSYM, 
		      &DensityCenter0, &DensityCenter1, &DensityCenter2);
      } // end input from parameter file
      fclose(RHfptr);
    }
  }

  // ensure that we're performing only Hydrogen chemistry
  if (RadHydroChemistry != 1) 
    ENZO_FAIL("RHIonizationSteepInitialize error: RadHydroChemistry must equal 1!");

  /* error checking */
  if (Mu != DEFAULT_MU) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: mu =%f assumed in initialization; setting Mu = %f for consistency.\n", DEFAULT_MU);
    Mu = DEFAULT_MU;
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }

  // convert input temperature to internal energy
  RadHydroTemperature = max(RadHydroTemperature,MIN_TEMP); // enforce minimum
  float mp = 1.67262171e-24;    // proton mass [g]
  float kb = 1.3806504e-16;     // boltzmann constant [erg/K]
  float HI = 1.0 - RadHydroInitialFractionHII;
  float HII = RadHydroInitialFractionHII;
  float ne = HII;
  float num_dens = HI + HII + ne;
  float mu = 1.0/num_dens;
  // correct mu if using a special model
  if ((RadHydroModel == 4) || (RadHydroModel == 5)) 
    mu = DEFAULT_MU;
  // compute the internal energy
  float RadHydroIEnergy = kb*RadHydroTemperature/mu/mp/(Gamma-1.0);	

  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->RHIonizationSteepInitializeGrid(
                        RadHydroChemistry, RadHydroNumDensity, 
			RadHydroDensityRadius, DensityCenter0, 
			DensityCenter1, DensityCenter2, RadHydroX0Velocity, 
			RadHydroX1Velocity, RadHydroX2Velocity, 
			RadHydroIEnergy, RadHydroRadiationEnergy, 
			RadHydroInitialFractionHII, local) == FAIL) {
      fprintf(stderr, "Error in RHIonizationSteepInitializeGrid.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }

  // set up field names and units
  // note: we must set up He species fields as well since Enzo 
  //       requires them for H chemistry (initialized to zero)
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  if (DualEnergyFormalism) 
    DataLabel[BaryonField++] = IEName;
  DataLabel[BaryonField++] = Vel0Name;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  DataLabel[BaryonField++] = RadName;
  DataLabel[BaryonField++] = DeName;
  DataLabel[BaryonField++] = HIName;
  DataLabel[BaryonField++] = HIIName;

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;
 
#endif

}

/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- CLUMP IONIZATION TEST
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


int RHIonizationClumpInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("Entering RHIonizationClumpInitialize routine\n");

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
  //  1. ambient density (should be very small) - free parameter
  //  2. ambient gas velocity - free parameter
  //  3. ambient gas temperature
  //  4. ambient radiation energy
  //  5. Hydrogen mass fraction 
  //  6. initial fraction HII
  //  7. Number of chemical species
  //  8. mesh spacing
  float RadHydroX0Velocity           = 0.0;
  float RadHydroX1Velocity           = 0.0;
  float RadHydroX2Velocity           = 0.0;
  float RadHydroNumDensityIn         = 0.04;
  float RadHydroNumDensityOut        = 0.0002;
  float RadHydroTemperatureIn        = 40.0;
  float RadHydroTemperatureOut       = 8000.0;
  float RadHydroRadiationEnergy      = 1.0e-20;
  float RadHydroHydrogenMassFraction = 1.0;
  float RadHydroInitialFractionHII   = 0.0;
  int   RadHydroChemistry            = 1;
  int   RadHydroModel                = 1;
  float ClumpCenterX                 = 1.54285e22;  // cm (5 kpc)
  float ClumpCenterY                 = 1.018281e22; // cm (3.3 kpc)
  float ClumpCenterZ                 = 1.018281e22; // cm (3.3 kpc)
  float ClumpRadius                  = 2.46856e21;  // cm (0.8 kpc)

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
	ret += sscanf(line, "RadHydroNumDensityIn = %"FSYM, 
		      &RadHydroNumDensityIn);
	ret += sscanf(line, "RadHydroNumDensityOut = %"FSYM, 
		      &RadHydroNumDensityOut);
	ret += sscanf(line, "RadHydroTemperatureIn = %"FSYM, 
		      &RadHydroTemperatureIn);
	ret += sscanf(line, "RadHydroTemperatureOut = %"FSYM, 
		      &RadHydroTemperatureOut);
	ret += sscanf(line, "RadHydroRadiationEnergy = %"FSYM, 
		      &RadHydroRadiationEnergy);
	ret += sscanf(line, "RadHydroInitialFractionHII = %"FSYM, 
		      &RadHydroInitialFractionHII);
	
	ret += sscanf(line, "ClumpCenter = %"FSYM" %"FSYM" %"FSYM,
		      &ClumpCenterX, &ClumpCenterY, &ClumpCenterZ);
	ret += sscanf(line, "ClumpRadius = %"FSYM, &ClumpRadius);
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

  // ensure that we're performing only Hydrogen chemistry
  if (RadHydroChemistry != 1) 
    ENZO_FAIL("RHIonizationClumpInitialize error: RadHydroChemistry must equal 1!");

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }

  // if temperature specified and not internal energy, perform conversion here
  RadHydroTemperatureIn = max(RadHydroTemperatureIn,MIN_TEMP); // enforce minimum
  RadHydroTemperatureOut = max(RadHydroTemperatureOut,MIN_TEMP); // enforce minimum
  float mp = 1.67262171e-24;    // proton mass [g]
  float kb = 1.3806504e-16;     // boltzmann constant [erg/K]
  float nH, HI, HII, ne, num_dens, mu;
  HI = 1.0 - RadHydroInitialFractionHII;
  HII = RadHydroInitialFractionHII;
  ne = HII;
  num_dens = HI + HII + ne;
  mu = 1.0/num_dens;
  // correct mu if using a special model
  if ((RadHydroModel == 4) || (RadHydroModel == 5)) 
    mu = DEFAULT_MU;
  // compute the internal energy
  float RadHydroIEnergyIn  = kb*RadHydroTemperatureIn/mu/mp/(Gamma-1.0);
  float RadHydroIEnergyOut = kb*RadHydroTemperatureOut/mu/mp/(Gamma-1.0);
  
  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->RHIonizationClumpInitializeGrid(
			RadHydroChemistry, RadHydroNumDensityIn, 
			RadHydroNumDensityOut, RadHydroX0Velocity, 
			RadHydroX1Velocity, RadHydroX2Velocity, 
			RadHydroIEnergyIn, RadHydroIEnergyOut, 
			RadHydroRadiationEnergy, RadHydroHydrogenMassFraction, 
			RadHydroInitialFractionHII, ClumpCenterX, ClumpCenterY, 
			ClumpCenterZ, ClumpRadius, local) == FAIL) {
      fprintf(stderr, "Error in RHIonizationClumpInitializeGrid.\n");
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

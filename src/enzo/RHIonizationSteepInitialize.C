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
#include "phys_constants.h"

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

  const char *kphHIName    = "HI_kph";
  const char *kphHeIName   = "HeI_kph";
  const char *kphHeIIName  = "HeII_kph";
  const char *gammaName    = "PhotoGamma";
  const char *kdissH2IName = "H2I_kdiss";
  const char *DensName  = "Density";
  const char *TEName    = "TotalEnergy";
  const char *IEName    = "Internal_Energy";
  const char *Vel0Name  = "x-velocity";
  const char *Vel1Name  = "y-velocity";
  const char *Vel2Name  = "z-velocity";
  const char *RadName   = "Grey_Radiation_Energy";
  const char *HIName    = "HI_Density";
  const char *HIIName   = "HII_Density";
  const char *HeIName   = "HeI_Density";
  const char *HeIIName  = "HeII_Density";
  const char *HeIIIName = "HeIII_Density";
  const char *DeName    = "Electron_Density";

  // local declarations
  char line[MAX_LINE_LENGTH];

  // Setup and parameters:
  float RadHydroX0Velocity           = 0.0;
  float RadHydroX1Velocity           = 0.0;
  float RadHydroX2Velocity           = 0.0;
  float RadHydroNumDensity           = 3.2;           // [cm^{-3}]
  float RadHydroDensityRadius        = 1.14375e-1;    // [code units]
  float DensityCenter0               = 0.0;
  float DensityCenter1               = 0.0;
  float DensityCenter2               = 0.0;
  float RadHydroTemperature          = 100.0;         // [K]
  float RadHydroRadiationEnergy      = 1.0e-20;
  float RadHydroHydrogenMassFraction = 1.0;
  float RadHydroInitialFractionHII   = 0.0;
  float RadHydroInitialFractionHeII  = 0.0;
  float RadHydroInitialFractionHeIII = 0.0;
  int   RadHydroChemistry            = 1;

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	// read relevant problem parameters
	sscanf(line, "RadHydroVelocity = %" FSYM" %" FSYM" %" FSYM,
		      &RadHydroX0Velocity, &RadHydroX1Velocity, 
		      &RadHydroX2Velocity);
	sscanf(line, "RadHydroChemistry = %" ISYM, 
		      &RadHydroChemistry);
	sscanf(line, "RadHydroNumDensity = %" FSYM, 
		      &RadHydroNumDensity);
	sscanf(line, "RadHydroDensityRadius = %" FSYM, 
		      &RadHydroDensityRadius);
	sscanf(line, "RadHydroTemperature = %" FSYM, 
		      &RadHydroTemperature);
	sscanf(line, "RadHydroRadiationEnergy = %" FSYM, 
		      &RadHydroRadiationEnergy);
	sscanf(line, "RadHydroInitialFractionHII = %" FSYM, 
		      &RadHydroInitialFractionHII);
	sscanf(line, "RadHydroHFraction = %" FSYM, 
		      &RadHydroHydrogenMassFraction);
	sscanf(line, "RadHydroInitialFractionHeII = %" FSYM, 
		      &RadHydroInitialFractionHeII);
	sscanf(line, "RadHydroInitialFractionHeIII = %" FSYM, 
		      &RadHydroInitialFractionHeIII);
	sscanf(line, "EtaCenter = %" FSYM" %" FSYM" %" FSYM, 
		      &DensityCenter0, &DensityCenter1, &DensityCenter2);
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

  // convert input temperature to internal energy
  RadHydroTemperature = max(RadHydroTemperature,MIN_TEMP); // enforce minimum
  float nH, HI, HII, nHe, HeI, HeII, HeIII, ne, num_dens, mu;
  if (RadHydroChemistry == 1) {
    HI = 1.0 - RadHydroInitialFractionHII;
    HII = RadHydroInitialFractionHII;
    ne = HII;
    num_dens = HI + HII + ne;
    mu = 1.0/num_dens;
  }
  else if (RadHydroChemistry == 3) {
    nH = RadHydroHydrogenMassFraction;
    nHe = (1.0 - RadHydroHydrogenMassFraction);
    HI = nH*(1.0 - RadHydroInitialFractionHII);
    HII = nH*RadHydroInitialFractionHII;
    HeII = nHe*RadHydroInitialFractionHeII;
    HeIII = nHe*RadHydroInitialFractionHeIII;
    HeI = nHe - HeII - HeIII;
    ne = HII + HeII/4.0 + HeIII/2.0;
    num_dens = 0.25*(HeI + HeII + HeIII) + HI + HII + ne;
    mu = 1.0/num_dens;
  }
  // compute the internal energy
  float RadHydroIEnergy = kboltz*RadHydroTemperature/mu/mh/(Gamma-1.0);	

  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->RHIonizationSteepInitializeGrid(
                        RadHydroChemistry, RadHydroNumDensity, 
			RadHydroDensityRadius, DensityCenter0, 
			DensityCenter1, DensityCenter2, RadHydroX0Velocity, 
			RadHydroX1Velocity, RadHydroX2Velocity, 
			RadHydroIEnergy, RadHydroRadiationEnergy, 
			RadHydroHydrogenMassFraction, 
			RadHydroInitialFractionHII, 
			RadHydroInitialFractionHeII, 
			RadHydroInitialFractionHeIII, local) == FAIL) {
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
  if ((RadHydroChemistry == 3) || (MultiSpecies > 0)) {
    DataLabel[BaryonField++] = HeIName;
    DataLabel[BaryonField++] = HeIIName;
    DataLabel[BaryonField++] = HeIIIName;
  }

  // if using external chemistry/cooling, set rate labels and update params
  if (RadiativeCooling) {
    DataLabel[BaryonField++] = kphHIName;
    DataLabel[BaryonField++] = gammaName;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      DataLabel[BaryonField++] = kphHeIName;
      DataLabel[BaryonField++] = kphHeIIName;
    }
    if (MultiSpecies > 1)
      DataLabel[BaryonField++] = kdissH2IName;
  }

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;
 
#endif

}

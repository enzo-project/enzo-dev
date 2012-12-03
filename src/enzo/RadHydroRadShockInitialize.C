/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- Grey Radiating Shock Test
/
/  written by: Daniel Reynolds and John Hayes
/  date:       July, 2007
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

#define DEF_MU 0.5  // fully ionized hydrogen gas

// function prototypes
int InitializeRateData(FLOAT Time);




int RadHydroRadShockInitialize(FILE *fptr, FILE *Outfptr, 
		               HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stdout,"Entering RadHydroRadShockInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "TotalEnergy";
  char *IEName    = "Internal_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName   = "Grey_Radiation_Energy";

  // local declarations
  int dim;

  // Setup and parameters:
  //  1. ambient density (g/cc)
  //  2. ambient gas temperature (K or eV, depending CGSType)
  //  3. ambient rad temperature (K or eV, depending CGSType)
  //  4. imposed fluid velocity  (cm/sec)
  //  5. coordinate along which to propagate shock {0,1,2}
  //  6. Problem Type (1 = astrophysical setup parameters;
  //                   2=  "lab" setup parameters, after Lowrie)

  float DensityConstant  = 1.0;
  float GasTempConstant  = 1.0;
  float RadTempConstant  = 1.0;
  float VelocityConstant = 1.0;
  int   ShockDir    = 0;
  int   CGSType     = 1;

  // overwrite parameters from RadHydroParamFile file, if it exists
  char line[MAX_LINE_LENGTH];
  int  ret;
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
	ret += sscanf(line, "DensityConstant = %"FSYM, &DensityConstant);
	ret += sscanf(line, "GasTempConstant = %"FSYM, &GasTempConstant);
	ret += sscanf(line, "RadTempConstant = %"FSYM, &RadTempConstant);
	ret += sscanf(line, "VelocityConstant = %"FSYM, &VelocityConstant);
	ret += sscanf(line, "ShockDir = %"ISYM, &ShockDir);
	ret += sscanf(line, "CGSType = %"ISYM, &CGSType);
      } // end input from parameter file
      fclose(RHfptr);
    }
  }

  /* error checking */
  if (Mu != DEF_MU) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: mu =%f assumed in initialization; setting Mu = %f for consistency.\n", DEF_MU);
    Mu = DEF_MU;
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }

  // set some physical constants
  float mp     = 1.67262171e-24;    // proton mass [g]
  float kb     = 1.3806504e-16;     // boltzmann constant [erg/K]
  float rc_lab = 137.20172;         // radiation constant [erg/cc/eV^4] 
  float rc_ast = 7.56e-15;          // radiation constant [erg/cc/K^4]
  float Cv     = 2.218056e+12;      // specific heat [erg/g/eV]  (From Lowrie setup)

  // Compute total gas energy (erg/gm) and radiation energy 
  // density (erg/cc) from input temperatures
  float gas_pressure;
  if ( CGSType == 1 ) 
    gas_pressure  = DensityConstant * kb * GasTempConstant / DEF_MU / mp;
  if ( CGSType == 2 ) 
    gas_pressure  = (Gamma - 1.0) * Cv * DensityConstant * GasTempConstant;

  float TEConstant = gas_pressure/(Gamma-1.0)/DensityConstant 
                   + 0.5 * VelocityConstant * VelocityConstant;

  float rc;
  if ( CGSType == 1 ) rc = rc_ast;
  if ( CGSType == 2 ) rc = rc_lab;

  float REConstant = rc * POW(RadTempConstant,4.0);

  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->RadHydroRadShockInitializeGrid(DensityConstant, 
                        TEConstant, REConstant, VelocityConstant, 
			ShockDir, local) == FAIL) {
      fprintf(stderr, "Error in RadHydroRadShockInitializeGrid.\n");
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


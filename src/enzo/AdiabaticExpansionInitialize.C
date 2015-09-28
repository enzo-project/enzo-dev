/***********************************************************************
/
/  INITIALIZE AN ADIABATIC EXPANSION TEST
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
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
#include "CosmologyParameters.h"
 
/* Set the mean molecular mass as in Grid_ComputeTemperatureField.C */
 
#define DEFAULT_MU 0.6
 
/* function prototypes */
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int CosmologyComputeTimeFromRedshift(FLOAT Redshift, FLOAT *TimeCodeUnits);

int AdiabaticExpansionInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *DebugName = "Debug";
  char *Phi_pName = "Phip";
  char *CRName = "CREnergyDensity";
  char *GPotName  = "Grav_Potential";

  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int ret;
 
  /* Error check. */
 
  if (!ComovingCoordinates) {
        ENZO_FAIL("ComovingCoordinates must be TRUE!");
  }
 
  /* set default parameters */
 
  float AdiabaticExpansionOmegaBaryonNow     = 1.0;  // standard
  float AdiabaticExpansionOmegaCDMNow        = 0.0;  // no dark matter
  float AdiabaticExpansionInitialTemperature = 200;  // degrees K
  float AdiabaticExpansionInitialUniformBField[MAX_DIMENSION];  // in Gauss
  float AdiabaticExpansionInitialVelocity    = 100;  // km/s
	float AdiabaticExpansionInitialCR          = 1.8e-15; // ergs/cm^3
 
  float InitialVels[MAX_DIMENSION];  // Initialize arrays
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    AdiabaticExpansionInitialUniformBField[dim] = 0.0;
    InitialVels[dim] = 0.0;
  }

  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "AdiabaticExpansionOmegaBaryonNow = %"FSYM,
		  &AdiabaticExpansionOmegaBaryonNow);
    ret += sscanf(line, "AdiabaticExpansionOmegaCDMNow = %"FSYM,
		  &AdiabaticExpansionOmegaCDMNow);
    ret += sscanf(line, "AdiabaticExpansionInitialTemperature = %"FSYM,
		  &AdiabaticExpansionInitialTemperature);
    ret += sscanf(line, "AdiabaticExpansionInitialUniformBField = %"FSYM" %"FSYM" %"FSYM,
		  AdiabaticExpansionInitialUniformBField,
		  AdiabaticExpansionInitialUniformBField+1,
		  AdiabaticExpansionInitialUniformBField+2);
    ret += sscanf(line, "AdiabaticExpansionInitialVelocity = %"FSYM,
		  &AdiabaticExpansionInitialVelocity);
    ret += sscanf(line, "AdiabaticExpansionInitialCR = %"FSYM,
      &AdiabaticExpansionInitialCR);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "AdiabaticExpansion"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }

  /* error checking */
  if (Mu != DEFAULT_MU) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: mu =%f assumed in initialization; setting Mu = %f for consistency.\n", DEFAULT_MU);
    Mu = DEFAULT_MU;
  }

 
  /* Get the units so we can convert temperature later. */
 
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1, PressureUnits=1.,MagneticUnits=1., a=1,dadt=0,CRUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, InitialTimeInCodeUnits) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  PressureUnits = DensityUnits * (LengthUnits/TimeUnits)*(LengthUnits/TimeUnits);
  MagneticUnits = sqrt(PressureUnits*4.0*M_PI);
	CRUnits = DensityUnits * (LengthUnits/TimeUnits)*(LengthUnits/TimeUnits);

  for (int dim = 0; dim < MAX_DIMENSION; dim++) 
    AdiabaticExpansionInitialUniformBField[dim] /= MagneticUnits;

  AdiabaticExpansionInitialCR /= CRUnits;

  /* Put inputs in a form that will be understood by InitializeUniformGrid. */
 
  float InitialTotalEnergy, InitialGasEnergy;
  InitialGasEnergy = AdiabaticExpansionInitialTemperature/TemperatureUnits /
    (Gamma - 1.0);
  if (MultiSpecies == FALSE) InitialGasEnergy /=  DEFAULT_MU;
  InitialTotalEnergy = InitialGasEnergy;
  InitialVels[0] = AdiabaticExpansionInitialVelocity/VelocityUnits*1.0e5;
  InitialTotalEnergy += 0.5*POW(InitialVels[0],2);
  for (int dim = 1; dim < MAX_DIMENSION; dim++)
    InitialVels[dim] = 0.0;
 
  /* set up grid */
 
  if (TopGrid.GridData->InitializeUniformGrid(
					      AdiabaticExpansionOmegaBaryonNow,
					      InitialTotalEnergy,
					      InitialGasEnergy, InitialVels,
					      AdiabaticExpansionInitialUniformBField,
					      AdiabaticExpansionInitialCR
					      ) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }
 
  /* set up field names and units */
 
  int i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;
  if (HydroMethod == MHD_RK) {
    DataLabel[i++] = BxName;
    DataLabel[i++] = ByName;
    DataLabel[i++] = BzName;
    DataLabel[i++] = PhiName;
  }
  if(UseDivergenceCleaning){
    DataLabel[i++] = Phi_pName;
    DataLabel[i++] = DebugName;
  }
	if(CRModel) DataLabel[i++] = CRName;
  if (WritePotential)
    DataLabel[i++] = GPotName;  
 

  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "AdiabaticExpansionOmegaBaryonNow     = %"FSYM"\n",
	    AdiabaticExpansionOmegaBaryonNow);
    fprintf(Outfptr, "AdiabaticExpansionOmegaCDMNow        = %"FSYM"\n",
	    AdiabaticExpansionOmegaCDMNow);
    fprintf(Outfptr, "AdiabaticExpansionInitialTemperature = %"FSYM"\n",
	    AdiabaticExpansionInitialTemperature);
    fprintf(Outfptr, "AdiabaticExpansionInitialUniformBField = ");
    WriteListOfFloats(Outfptr, 3, AdiabaticExpansionInitialUniformBField);

    fprintf(Outfptr, "AdiabaticExpansionInitialVelocity    = %"FSYM"\n\n",
	    AdiabaticExpansionInitialVelocity);
    fprintf(Outfptr, "AdiabaticExpansionInitialCR          = %"FSYM"\n\n",
      AdiabaticExpansionInitialCR);
  }
 
  return SUCCESS;
}

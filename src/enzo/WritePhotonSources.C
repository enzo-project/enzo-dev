#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);

int WritePhotonSources(FILE *fptr, FLOAT CurrentTime)
{

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, CurrentTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");

  }

  RadiationSourceEntry *RS = GlobalRadiationSources->NextSource;
  int i, nSources = 0;

  while (RS != NULL) {
    nSources++;
    RS = RS->NextSource;
  }

  RS = GlobalRadiationSources->NextSource;

  fprintf(fptr, "PhotonTestNumberOfSources       = %"ISYM"\n", nSources);
  for (i = 0; i < nSources; i++) {
    fprintf(fptr, "PhotonTestSourceType[%"ISYM"]         = %"ISYM"\n", i, RS->Type);
    fprintf(fptr, "PhotonTestSourcePosition[%"ISYM"]     = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
	    i, RS->Position[0], RS->Position[1], RS->Position[2]);
    if (RS->Orientation != NULL) {
      fprintf(fptr, "PhotonTestSourceOrientation[%"ISYM"]     = %"GSYM" %"GSYM" %"GSYM"\n",
	      i, RS->Orientation[0], RS->Orientation[1], RS->Orientation[2]);
    }
    fprintf(fptr, "PhotonTestSourceLuminosity[%"ISYM"]   = %"GSYM"\n",
	    i, RS->Luminosity/TimeUnits*pow(LengthUnits,3));
    fprintf(fptr, "PhotonTestSourceLifeTime[%"ISYM"]     = %"GSYM"\n", i, RS->LifeTime);
    fprintf(fptr, "PhotonTestSourceCreationTime[%"ISYM"] = %"GSYM"\n", i, 
	    RS->CreationTime);
    fprintf(fptr, "PhotonTestSourceRampTime[%"ISYM"]     = %"GSYM"\n", i, RS->RampTime);
    fprintf(fptr, "PhotonTestSourceEnergyBins[%"ISYM"]   = %"ISYM"\n", i, RS->EnergyBins);
    fprintf(fptr, "PhotonTestSourceSED[%"ISYM"]          = ", i);
    WriteListOfFloats(fptr, RS->EnergyBins, RS->SED);
    fprintf(fptr, "PhotonTestSourceEnergy[%"ISYM"]       = ", i);
    WriteListOfFloats(fptr, RS->EnergyBins, RS->Energy);
    fprintf(fptr, "\n");
    RS = RS->NextSource;
  }

  return SUCCESS;

}

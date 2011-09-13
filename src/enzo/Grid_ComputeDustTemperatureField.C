/***********************************************************************
/
/  GRID CLASS (COMPUTE THE TEMPERATURE FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the pressure at the requested time.  The pressure here is
//   just the ideal-gas equation-of-state.
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
#include "CosmologyParameters.h"
 
/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int f, int farray[], int n);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

extern "C" void FORTRAN_NAME(calc_tdust_3d)(
	float *d, float *de, float *HI, float *HII, 
	float *HeI, float *HeII, float *HeIII,
	float *HM, float *H2I, float *H2II, 
	int *in, int *jn, int *kn, 
	int *nratec, int *iexpand,
	int *ispecies, int *idim,
	int *is, int *js, int *ks, 
	int *ie, int *je, int *ke, 
	float *aye, float *temstart, float *temend,
	float *gasgra,
	float *utem, float *uxyz, float *uaye,
	float *urho, float *utim,
	float *gas_temp, float *dust_temp);

int grid::ComputeDustTemperatureField(float *temperature, float *dust_temperature)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Return if not using MultiSpecies chemistry. */
  if (!MultiSpecies) {
    ENZO_FAIL("Dust temperature calculation requires MultiSpecies > 0.\n");
  }

  int DensNum;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  
  /* Compute the size of the fields. */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0)
    ENZO_FAIL("Cannot find density.");

  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
 
  /* Find the units. */
 
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
 
    aUnits = 1.0/(1.0 + InitialRedshift);
  }
  float afloat = float(a);

  /* Find Multi-species fields. */

  if (MultiSpecies) {  
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }
  }


  /* Call the appropriate FORTRAN routine to do the work. */

    FORTRAN_NAME(calc_tdust_3d)(
       BaryonField[DensNum], BaryonField[DeNum], BaryonField[HINum], BaryonField[HIINum],
       BaryonField[HeINum], BaryonField[HeIINum], BaryonField[HeIIINum],
       BaryonField[HMNum], BaryonField[H2INum], BaryonField[H2IINum],
       GridDimension, GridDimension+1, GridDimension+2,
       &CoolData.NumberOfTemperatureBins, &ComovingCoordinates,
       &MultiSpecies, &GridRank, 
       GridStartIndex, GridStartIndex+1, GridStartIndex+2,
       GridEndIndex, GridEndIndex+1, GridEndIndex+2,
       &afloat, &CoolData.TemperatureStart, &CoolData.TemperatureEnd, 
       CoolData.gas_grain, 
       &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
       temperature, dust_temperature);

  return SUCCESS;
}

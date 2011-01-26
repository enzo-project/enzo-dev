/***********************************************************************
/
/  GRID CLASS (COMPUTE THE X-RAY EMISSIVITY FIELD)
/
/  written by: Greg Bryan
/  date:       February, 2000
/  modified1:
/
/  PURPOSE:  Compute X-ray emissivity in given band.  Input field must
/            be the temperature.  Units are in 10^-23 ergs/cm^3/s.
/            Does not allocate xray_emissivity field.
/
/  RETURNS:
/
************************************************************************/
 
// Use a lookup table to compute the X-ray emissivity in the requested
//   band.  Current version assumes complete ionization.
 
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
 
int FindField(int f, int farray[], int n);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
/* The following variables are used only in computing the emissivity,
   and are read in from the specified file which is created with a
   Raymond-Smith code.
   SpectrumEnergykeV is the energy of each bin (in keV),
   SpectrumEmissivity is the emissivity in that bin, at that temperature,
   temp1, temp2 are log10(K) of the first and last temperature bins. */
 
static int FirstTimeCalled = TRUE;
static float *SpectrumEnergykeV, *SpectrumEmissivity, *TotalEmissivity;
static int NumberOfSpectralBins;
static int NumberOfTemperatureBins;
static float temp1, temp2;
 
 
 
int grid::ComputeXrayEmissivity(float *temperature,
				float *xray_emissivity, float keV1, float keV2,
				char *XrayTableFileName)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* Compute the size of the fields. */
 
  int DensNum, i, j, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    ENZO_FAIL("Cannot find density.\n");
  }
 
  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1, CurrentRedshift = 1.0;
  FLOAT a, dadt;
 
  /* Find the temperature units. */
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    CurrentRedshift = (1.0+InitialRedshift)/a - 1.0;
  }
 
  /* If not already done so, read in spectral table. */
 
  if (FirstTimeCalled) {
 
    /* open file */
 
    FILE *fptr = fopen(XrayTableFileName, "r");
    char *dummy = new char[MAX_LINE_LENGTH];
    if (fptr == NULL) {
      ENZO_VFAIL("Error in spectral file %s\n", XrayTableFileName)
    }
 
    /* read in table description: ignore first two lines, then read
       # of temperature bins and min/max temperature (log10(K)), then
       skip one line and read in number of spectral bins.  */
 
    fgets(dummy, MAX_LINE_LENGTH, fptr);
    fgets(dummy, MAX_LINE_LENGTH, fptr);
    if (fscanf(fptr, "%"ISYM" %"FSYM" %"FSYM, &NumberOfTemperatureBins, &temp1, &temp2)
	!= 3) {
      ENZO_FAIL("Error reading temperature info\n");
    }
    if (debug)
      printf("NumberOfTemperatureBins = %"ISYM" (%"GSYM"-%"GSYM")\n", NumberOfTemperatureBins,
	     temp1, temp2);
    fgets(dummy, MAX_LINE_LENGTH, fptr);
    fgets(dummy, MAX_LINE_LENGTH, fptr);
    fscanf(fptr, "%"ISYM, &NumberOfSpectralBins);
    if (debug)
      printf("NumberOfSpectralBins = %"ISYM"\n", NumberOfSpectralBins);
 
    /* Skip five lines and read in spectral energy bins and then emissivity
       table (one per temperature bin). */
 
    for (j = 0; j < 5; j++)
      fgets(dummy, MAX_LINE_LENGTH, fptr);
    SpectrumEnergykeV = new float[NumberOfSpectralBins];
    SpectrumEmissivity =
              new float[NumberOfTemperatureBins*NumberOfSpectralBins];
 
    for (j = 0; j < NumberOfSpectralBins; j++)
      if (fscanf(fptr, "%"FSYM, SpectrumEnergykeV+j) != 1) {
	ENZO_VFAIL("Error reading energy: %s\n", XrayTableFileName)
      }
 
    int n = 0;
    for (i = 0; i < NumberOfTemperatureBins; i++)
      for (j = 0; j < NumberOfSpectralBins; j++, n++) {
	if (fscanf(fptr, "%"FSYM, SpectrumEmissivity+n) != 1) {
	  ENZO_VFAIL("Error reading file %s\n", XrayTableFileName)
	}
      }
 
    /* Create a single look-up table (to speed integration over bins). */
 
    float keV1_shifted = keV1*(1.0+CurrentRedshift),
          keV2_shifted = keV2*(1.0+CurrentRedshift);
    TotalEmissivity = new float[NumberOfTemperatureBins];
    for (i = 0; i < NumberOfTemperatureBins; i++) {
      TotalEmissivity[i] = 0;
      for (j = 0; j < NumberOfSpectralBins; j++)
	if (SpectrumEnergykeV[j] > keV1_shifted &&
	    SpectrumEnergykeV[j] < keV2_shifted)
	  TotalEmissivity[i] += SpectrumEmissivity[i*NumberOfSpectralBins+j];
    }
 
    FirstTimeCalled = FALSE;
    delete [] dummy;
  }
 
  /* Loop over grid and compute emissivity
     (fh is hydrogen fraction by mass). */
 
  float temp, ne, nH, frac, fh = 0.76;
  float ConvertToNumberDensity = DensityUnits/1.67e-24;
  float deltemp = (temp2-temp1)/float(NumberOfTemperatureBins-1);
  for (i = 0; i < size; i++) {
 
    /* Get temperature and electron number density (assume fully ionized). */
 
    temp = log10(temperature[i]);
    nH   = fh*BaryonField[DensNum][i] * ConvertToNumberDensity;
    ne   = nH + 0.5*(1.0-fh)*BaryonField[DensNum][i] * ConvertToNumberDensity;
 
#define NO_ENTROPY_FLOOR
#ifdef ENTROPY_FLOOR
 
    /* If the density is greater than 200 times the mean (note: the problem
       is that we don't have access to the mean omega baryon, so just guess
       it -- fix this), then set the density such that the "entropy" is no
       greater than 100 keV cm^2. */
 
    if (BaryonField[DensNum][i] > 0.04/OmegaMatterNow*200.0) {
      //      ne = min(ne, POW(temperature[i]/(100*1.16e7), 1.5));
      ne = POW((POW(ne, -0.666666) + 200*1.16e7/temperature[i]), -1.5);
      nH = ne*2.0*fh/(1.0+fh);
    }
 
#endif /* ENTROPY_FLOOR */
 
    /* Look-up temperature and compute temperature bins to interpolate from. */
 
    temp = min(max(temp, temp1), temp2*0.999);
    j = int((temp-temp1)/deltemp);
    frac = (temp - (temp1+j*deltemp))/deltemp;
 
    if (frac < -0.01 || frac > 1.01) {
      ENZO_VFAIL("prob: %"GSYM" %"ISYM" %"ISYM" %"GSYM"\n", frac, j, i, temp)

    }
    frac = min(max(frac, 0), 1);
 
    /* Add up the emissivity in selected (redshifted) band and multiply
       by n_e * n_H (assuming complete ionization!),
       in units of 10^-23 erg cm^-3 s^-1. */
 
    xray_emissivity[i] = TotalEmissivity[j  ]*     frac +
                         TotalEmissivity[j+1]*(1.0-frac);
    xray_emissivity[i] *= ne*nH;
    xray_emissivity[i] = max(xray_emissivity[i], 1e-20);
 
  }
 
  return SUCCESS;
}

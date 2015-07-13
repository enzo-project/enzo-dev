/***********************************************************************
/
/  GRID CLASS (COMPUTE IONIZED ELEMENTAL DENSITY FIELD)
/
/  written by: Greg Bryan
/  date:       February, 2000
/  modified1:
/
/  PURPOSE:  This routine computes the elemental number density due to
/            collisional ionized of one of two species: OVII, OVIII.
/            Some code originally from Taotao Fang.
/            It also can use a generalized density, temperature look-up table
/
/            Temperature field must be passed in; elemental density
/            field must be allocated before calling this routine.
/             Type refers to element type:
/               0 - OVII             1 - OVIII
/               2 - SiXIII           3 - SIXIV
/               4 - FeXXV            5 - FeXXVI
/               6 - Oxygen using lookup table
/
/            Metallicity is 0.3 unless there is a metallcity field
/
/  RETURNS:
/
************************************************************************/
 
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
#include "Grid.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
static int TableRead = FALSE;
static float *LookupTable, *TableDensity, *TableTemperature;
static int TableSize[2];
 
/* Look-up table for collisional ionization; the first column is log T, in K,
   the second is fraction abundance, in log. */
 
float OVIILookupTable[] = {
5.2,     4.845,
5.3,     2.554,
5.4,     0.937,
5.5,     0.200,
5.6,     0.044,
5.7,     0.013,
5.8,     0.005,
5.9,     0.003,
6.0,     0.004,
6.1,     0.019,
6.2,     0.077,
6.3,     0.232,
6.4,     0.547,
6.5,     1.038,
6.6,     1.610,
6.7,     2.175,
6.8,     2.699,
6.9,     3.178,
7.0,     3.616,
7.1,     4.017,
7.2,     4.387,
7.3,     4.730
};
 
float OVIIILookupTable[] = {
5.8,     4.656,
5.9,     3.310,
6.0,     2.234,
6.1,     1.403,
6.2,     0.807,
6.3,     0.452,
6.4,     0.355,
6.5,     0.498,
6.6,     0.766,
6.7,     1.057,
6.8,     1.331,
6.9,     1.578,
7.0,     1.799,
7.1,     1.997,
7.2,     2.175,
7.3,     2.336,
7.4,     2.487,
7.5,     2.625,
7.6,     2.755,
7.7,     2.687
};
 
float SiXIIILookupTable[] = {
6.00,  4.675,
6.10,  2.376,
6.20,  0.926,
6.30,  0.295,
6.40,  0.106,
6.50,  0.049,
6.60,  0.029,
6.70,  0.029,
6.80,  0.055,
6.90,  0.130,
7.00,  0.283,
7.10,  0.539,
7.20,  0.887,
7.30,  1.283,
7.40,  1.686,
7.50,  2.074,
7.60,  2.440,
7.70,  2.783,
7.80,  3.103,
7.90,  3.405,
8.00,  3.690
};
 
 
float SiXIVLookupTable[] = {
6.40,  4.447,
6.50,  3.268,
6.60,  2.342,
6.70,  1.614,
6.80,  1.055,
6.90,  0.658,
7.00,  0.424,
7.10,  0.354,
7.20,  0.420,
7.30,  0.566,
7.40,  0.744,
7.50,  0.926,
7.60,  1.102,
7.70,  1.267,
7.80,  1.421,
7.90,  1.566,
8.00,  1.702
};
 
 
float FeXXVLookupTable[] = {
6.9,  4.755,
7.0,  2.865,
7.1,  1.601,
7.2,  0.842,
7.3,  0.453,
7.4,  0.269,
7.5,  0.180,
7.6,  0.139,
7.7,  0.132,
7.8,  0.160,
7.9,  0.227,
8.0,  0.336,
8.1,  0.487,
8.2,  0.671,
8.3,  0.879,
8.4,  1.101,
8.5,  1.328,
8.6,  1.557,
8.7,  1.787,
8.8,  2.015,
8.9,  2.242,
9.0,  2.469
};
 
 
float FeXXVILookupTable[] = {
7.10,  4.891,
7.20,  3.523,
7.30,  2.620,
7.40,  1.964,
7.50,  1.458,
7.60,  1.061,
7.70,  0.757,
7.80,  0.540,
7.90,  0.407,
8.00,  0.350,
8.10,  0.355,
8.20,  0.405,
8.30,  0.482,
8.40,  0.574,
8.50,  0.675,
8.60,  0.780,
8.70,  0.891,
8.80,  0.999,
8.90,  1.105,
9.00,  1.208
};
 
 
int ElementalTableSizes[6] = {
                     sizeof(OVIILookupTable)/(2*sizeof(float)),
                     sizeof(OVIIILookupTable)/(2*sizeof(float)),
                     sizeof(SiXIIILookupTable)/(2*sizeof(float)),
                     sizeof(SiXIVLookupTable)/(2*sizeof(float)),
                     sizeof(FeXXVLookupTable)/(2*sizeof(float)),
                     sizeof(FeXXVILookupTable)/(2*sizeof(float)) };
 
/* Solar abundances ratio of O/H by number (first line),
   Si/H by number (second line), Fe/H by number (third line). */
 
float ElementalBaseAbundance[7] = {8.53e-4, 8.53e-4,
			           3.58e-5, 3.58e-5,
			           3.23e-5, 3.23e-5,
                                   8.53e-4};
			
 
int grid::ComputeElementalDensity(float *temperature,
				  float *elemental_density, int Type)
{
 
  const float DefaultMetallicity = 0.3;
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* Compute the size of the fields. */
 
  int DensNum, i, j, n, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    ENZO_FAIL("Cannot find density.\n");
  }
 
  /* Find metallicity field and set flag. */
 
  int MetallicityField = FALSE, MetalNum;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;
 
  /* Find the temperature units if we are using comoving coordinates. */
 
  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
 
  /* Set lookup table pointer and size. */
 
  switch (Type) {
  case 0: LookupTable = OVIILookupTable; break;
  case 1: LookupTable = OVIIILookupTable; break;
  case 2: LookupTable = SiXIIILookupTable; break;
  case 3: LookupTable = SiXIVLookupTable; break;
  case 4: LookupTable = FeXXVLookupTable; break;
  case 5: LookupTable = FeXXVILookupTable; break;
  case 6: break; /* general Oxygen */
  default:
    ENZO_VFAIL("Unrecognized element type: %"ISYM"\n", Type)
  }
 
  if (Type < 6) {
 
    /* Set up for a table which is already there. */
 
    TableSize[0] = ElementalTableSizes[Type];
    TableSize[1] = 1;
    TableRead = FALSE;
 
  } else {
 
    /* Read in general table if not done so already. */
 
    if (TableRead == FALSE) {
 
      char filename[] = "fraction_data.out";
      FILE *fptr;
      float DummyFloat;
 
      printf("reading %s\n", filename);
      if ((fptr = fopen(filename, "r")) == NULL) {
	ENZO_VFAIL("Erroring opening %s.\n", filename)
      }
      if (fscanf(fptr, "NumberOfDensityPoints = %"ISYM"\n", TableSize) != 1) {
	ENZO_FAIL("Erroring reading number of density points\n");
      }
      if (fscanf(fptr, "NumberOfTemperaturePoints = %"ISYM"\n", TableSize+1) != 1) {
	ENZO_FAIL("Erroring reading number of temperature points\n");
      }
      printf("NumberOfBins (temp,dens) = %"ISYM",%"ISYM"\n", TableSize[0], TableSize[1]);
 
      TableDensity     = new float[TableSize[0]];
      TableTemperature = new float[TableSize[1]];
      LookupTable      = new float[TableSize[0]*TableSize[1]];
      for (n = 0, j = 0; j < TableSize[1]; j++)
	for (i = 0; i < TableSize[0]; i++, n++)
	  if (fscanf(fptr, "%"FSYM" %"FSYM" %"FSYM" %"FSYM,
		     TableDensity+i, TableTemperature+j,
		     LookupTable+n, &DummyFloat) != 4) {
	    ENZO_VFAIL("Error reading table %"ISYM" %"ISYM"\n", i, j)
	  }
      fclose(fptr);
      TableRead = TRUE;
 
    }
  }
 
  /* Loop over grid and compute emissivity
     (fh is hydrogen fraction by mass; 1.67e-24 is proton mass). */
 
  float fh = 0.76, ConvertToNumberDensity = DensityUnits/1.67e-24;
  float nElement, nH, LogFraction, Metallicity, Fraction,
        deld, delt, logd, logt, dd, dt;
  int ilogd, ilogt;
 
  deld = (TableDensity[1]    -TableDensity[0]    );
  delt = (TableTemperature[1]-TableTemperature[0]);
 
  for (i = 0; i < size; i++)
 
    if (BaryonField[DensNum][i] > 0) {
 
    /* Get metallicity for this cell. */
 
    if (MetallicityField)
      Metallicity = BaryonField[MetalNum][i]/BaryonField[DensNum][i];
    else
      Metallicity = DefaultMetallicity;
 
#define DONT_USE_DENSITY_DEPENDENT_METALLICITY
#ifdef USE_DENSITY_DEPENDENT_METALLICITY
 
    /* This is a fit to Cen & Ostriker metallicity function from simul. */
 
    float baryonfraction = 0.1333;
    Metallicity = POW(10.0,
              0.34*log10(BaryonField[DensNum][i]/baryonfraction) - 1.62);
 
#endif /* USE_DENSITY_DEPENDENT_METALLICITY */
 
    /* Get elemental number density (of all ionization states). */
 
    nH   = fh*BaryonField[DensNum][i] * ConvertToNumberDensity;
    nElement = nH * ElementalBaseAbundance[Type] * Metallicity;
 
    /* Convert log(temperature.) into an integer in the table */
 
    logt = log10(max(temperature[i], 1));
    logt = min(max(logt, TableTemperature[0]),
	       TableTemperature[TableSize[1]-1]);
 
    ilogt = int((logt-TableTemperature[0])/delt);
    ilogt = min(ilogt, TableSize[1]-1);
 
    dt = (logt-TableTemperature[ilogt])/delt;
 
    if (TableSize[1] > 1) {
 
      /* Convert log(densit) into an integer in the table */
 
      logd = log10(max(nH, tiny_number));
      logd = min(max(logd, TableDensity[0]), TableDensity[TableSize[0]-1]);
 
      ilogd = int((logd-TableDensity[0])/deld);
      ilogd = min(ilogd, TableSize[0]-1);
 
      dd = (logd-TableDensity[ilogd])/deld;
 
      /* Lookup in 2-dimensional table (linear fraction). */
 
      Fraction =
       LookupTable[ ilogt   *TableSize[0]+ilogd  ]*(1.0-dd)*(1.0-dt) +
       LookupTable[ ilogt   *TableSize[0]+ilogd+1]*     dd *(1.0-dt) +
       LookupTable[(ilogt+1)*TableSize[0]+ilogd  ]*(1.0-dd)*     dt  +
       LookupTable[(ilogt+1)*TableSize[0]+ilogd+1]*     dd *     dt;
 
    } else {
 
      /* Lookup in 1-dimensional (temperature only) table (log fraction). */
 
      LogFraction =
       LookupTable[ilogt  ]*(1.0-dt) +
       LookupTable[ilogt+1]*     dt;
 
      Fraction = POW(10, -LogFraction);
 
    }
 
    /* Set density based on fraction. */
 
    elemental_density[i] = nElement * Fraction;
 
  } // end: loop over grid
 
#ifdef UNUSED
 
    if (temperature[i] > 0) {
 
      /* Get temperature and oxygen/Si/Fe number density. */
 
 
      /* Based on density and temperature,
	 look up relative fraction from table. */
 
      for (j = 0; j < TableSize; j++)
	if (LookupTable[j*2] > temp)
	  break;
 
      LogFraction = 20.0; /* something very small. */
      if (j > 0 && j < TableSize) {

	delTemp = LookupTable[j*2] - LookupTable[j*2-2];
	LogFraction = LookupTable[j*2-1]+(temp - LookupTable[j*2-2])/delTemp*
	  (LookupTable[j*2+1]-LookupTable[j*2-1]);
      }
 
      /* Set density based on fraction. */
 
      elemental_density[i] = nElement * POW(10.0, -LogFraction);
 
    } else
      elemental_density[i] = tiny_number;
 
#endif /* UNUSED */
 
  return SUCCESS;
}

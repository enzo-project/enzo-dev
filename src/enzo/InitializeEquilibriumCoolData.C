/***********************************************************************
/
/  INITIALIZE THE EQUILIBRIUM COOL RATES
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
int InitializeEquilibriumCoolData(FLOAT Time)
{
 
  /* Declarations. */
 
  FLOAT a = 1, dadt;
  float temp, dummy;
 
  /* Set relevant CoolData parameters. */
 
  CoolData.HydrogenFractionByMass   = 0.76;
 
  /* Open input file for data. */
 
  if( RadiativeCoolingModel == 1){
      
  FILE *fptr = fopen("cool_rates.in", "r");
  if (fptr == NULL) {
    ENZO_FAIL("Error opening cool_rates.in\n");
  }
 
  /* Read rate data, skipping over comments (count first). */
 
  int index = 0;
  char line[MAX_LINE_LENGTH];
  CoolData.NumberOfTemperatureBins = 0;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    if (line[0] != '#') CoolData.NumberOfTemperatureBins++;
  CoolData.EquilibriumRate = new float[CoolData.NumberOfTemperatureBins];
  rewind(fptr);
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    if (line[0] != '#')
      if (sscanf(line, "%"FSYM" %"FSYM" %"FSYM, &temp, &dummy,
		 &CoolData.EquilibriumRate[index]) == 3) {
	CoolData.EquilibriumRate[index] =
	  POW(10, CoolData.EquilibriumRate[index]);
	if (index == 0)
	  CoolData.TemperatureStart = POW(10, temp);
	index++;
      }
  CoolData.NumberOfTemperatureBins = index;
  CoolData.TemperatureEnd = POW(10, temp);
  fclose(fptr);
 
  if (debug) {
    printf("InitializeEqCoolData: NumberOfTemperatureBins = %"ISYM"\n",
		    CoolData.NumberOfTemperatureBins);
    printf("InitializeEqCoolData: TemperatureStart = %"GSYM"\n",
		    CoolData.TemperatureStart);
    printf("InitializeEqCoolData: TemperatureEnd = %"GSYM"\n",
		    CoolData.TemperatureEnd);
  }
 
  /* If using cosmology, compute the expansion factor and get units. */
 
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
 
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt)
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");

    }
 
    aUnits = 1.0/(1.0 + InitialRedshift);
 
  }
 
  /* Get conversion units */
 
  /* t/x/dbase1 is the number (z dependant) that converts from the
       dimensionless code units to physical units.  Also, in the
       code aye = 1 at z=zinit, so to convert the usual a (=1 at z=0)
       to a' (written in the code as aye), we use a = a'*[a]  */
 
  /* Set the dimension of the cooling coefficients (including constants)
     (this equation has a rho because e is the specific energy, not
      energy/unit volume).
       delta(e)  = L     * n1        * n2        * dt     / dens   / a^3
       [e]       = L     * [dens]/mh * [dens]/mh * [time] / [dens] / [a]^3
       delta(e') = L'    * n1'       * n2'       * dt'    / dens'  / a'^3 [']
     so L = [L] * L' where [L] = [e] * mh**2 * [a]^3 / ([dens] * [time]) [']
       but [e] = ([a]*[x])**2 / [time]**2 and ([a] = 1 / (1 + zri) )
      [L] = ([a]**5 * [x]**2 * mh**2) / ([dens] * [time]**3)  */
 
  double tbase1 = TimeUnits;
  double xbase1 = LengthUnits/(a*aUnits);
  double dbase1 = DensityUnits * POW(a*aUnits, 3);
  double mh = 1.67e-24;
  double CoolUnit = (POW(aUnits,5) * POW(xbase1,2) * POW(mh,2)) /
                    (POW(tbase1,3) * dbase1);
  for (index = 0; index < CoolData.NumberOfTemperatureBins; index++)
    CoolData.EquilibriumRate[index] /= CoolUnit;
 
  /*  Photoelectric heating by UV-irradiated dust (Wolfire 1995)
      Default is 8.5e-26 for epsilon=0.05, G_0=1.7 (rate in erg s^-1 cm^-3) 
      The extra factor of dbase1/mh is taken out when used in cool1d.F */

  CoolData.gammah = PhotoelectricHeatingRate / CoolUnit;


  /* Output cooling rate in code units. */

  if( MyProcessorNumber == ROOT_PROCESSOR ){
    fptr = fopen("cool_rates.out", "w");
    for (index = 0; index < CoolData.NumberOfTemperatureBins; index++)
      fprintf(fptr, "%"GSYM" %"GSYM"\n", log10(CoolData.TemperatureStart) +
       (log10(CoolData.TemperatureEnd)-
        log10(CoolData.TemperatureStart)) * float(index)/
       float(CoolData.NumberOfTemperatureBins-1),
       CoolData.EquilibriumRate[index]);
    fclose(fptr);
  } // end rootgrid if

  }else if( RadiativeCoolingModel == 3){
      CoolData.NumberOfTemperatureBins = 2;
      CoolData.EquilibriumRate = new float[CoolData.NumberOfTemperatureBins];
      CoolData.TemperatureEnd = 1e4  ;
      CoolData.TemperatureStart = 18;
  }
 
  return SUCCESS;
}

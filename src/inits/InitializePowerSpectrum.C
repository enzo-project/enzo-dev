/***********************************************************************
/
/  INITIALIZE POWER SPECTRUM AND CREATE LOOK-UP TABLE
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:  Robert Harkness
/  date:       November, 2003
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
 
#include "macros_and_parameters.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "PowerSpectrumParameters.h"
 
// function prototypes
 
extern "C" void FORTRAN_NAME(qromo)(FLOAT (*func)(FLOAT*), FLOAT *a, FLOAT *b,
				FLOAT *ss, void (*choose)());
extern "C" void FORTRAN_NAME(midpnt)();
extern "C" void FORTRAN_NAME(midinf)();
extern "C" void FORTRAN_NAME(set_common)(FLOAT *lam0_in, FLOAT *omega0_in,
					 FLOAT *zri_in, FLOAT *hub_in);
extern "C" FLOAT FORTRAN_NAME(calc_growth)(FLOAT *z);
FLOAT ps_tophat(FLOAT *kptnr);
FLOAT EvaluatePowerSpectrum(FLOAT k, int Species);
 
  
 
int InitializePowerSpectrum()
{
  // Set amplitude based on sigma8
 
  FLOAT Zero = 0.0, Huge = 1.0e2;
  FLOAT Pi = 3.14159265358979324;
  FLOAT s1, s2, kbreak;

  if (debug) printf("Initializing power spectrum\n");
  
  if (debug) printf("Redshift = %e\n",Redshift);

  // Set-up for integration
 
  Normalization = 1.0;
  GrowthFactor  = 1.0;
  TophatRadius = 8.0/HubbleConstantNow;
 
  // Integrate smoothed density field
 
  kbreak = 2.0*Pi/TophatRadius;

  FORTRAN_NAME(qromo)(&ps_tophat, &Zero, &kbreak, &s1, &FORTRAN_NAME(midpnt));
  FORTRAN_NAME(qromo)(&ps_tophat, &kbreak, &Huge, &s2, &FORTRAN_NAME(midinf));
 
  // Set Normalization by sigma8
 
  Normalization = sigma8*sigma8/(s1 + s2);
  if (debug) printf("Found a normalization at 8/h Mpc as %g in input power spectrum.\n", 
		    sqrt(s1+s2));
  if (debug) printf("Set it to input sigma8 at 8/h Mpc of %g.\n", sigma8);
 
  return SUCCESS;
}
 
 
 
// mqk: Pass in the name of the power spectrum output file. This
// allows separate outputs of initial and z=0 power spectra, which
// differ by more than just the growth factor in the CMBFAST case.

int MakePowerSpectrumLookUpTable(char *PowerSpectrumFilename)
{
 
  int i, species;
  FLOAT k, k1, k2, delk, Zero = 0;
 
  if (debug) printf("Generating look-up table\n");
 
  /* Compute the linear growth rate from z=init to 0. */


  // mqk: Allow arbitrary redshift (i.e. not necessarily equal to
  // InitialRedshift).

  // FORTRAN_NAME(set_common)(&OmegaLambdaNow, &OmegaMatterNow, &InitialRedshift,
  // 			   &HubbleConstantNow);
  // GrowthFactor = FORTRAN_NAME(calc_growth)(&InitialRedshift) /
  //                FORTRAN_NAME(calc_growth)(&Zero);

  FORTRAN_NAME(set_common)(&OmegaLambdaNow, &OmegaMatterNow, &Redshift,
			   &HubbleConstantNow);
  GrowthFactor = FORTRAN_NAME(calc_growth)(&Redshift) /
                 FORTRAN_NAME(calc_growth)(&Zero);
  if (debug) printf("GrowthFactor = %"GSYM"\n", GrowthFactor);
 
  /* Prepare PS lookup table. */
 
  k1 = log(kmin);
  k2 = log(kmax);
  delk = (k2 - k1)/(NumberOfkPoints-1);
  for (species = 0; species < MAX_SPECIES; species++)
    PSLookUpTable[species] = new FLOAT[NumberOfkPoints];
 
  /* Generate table for the mean mass power spectrum. */
 
  for (i = 0; i < NumberOfkPoints; i++) {
    k = exp(k1 + i*delk);
    PSLookUpTable[0][i] = EvaluatePowerSpectrum(k, 0);
  }
 
  /* If using psfunc, then just duplicate the table for the other species. */
 
  if (PowerSpectrumType < 20)
    for (species = 1; species < 3; species++)
      PSLookUpTable[species] = PSLookUpTable[0];
 
  /* If using CMBFast, then get the power spectrum for the other species. */
 
  if (PowerSpectrumType == 20) {
    for (species = 1; species < 3; species++)
      for (i = 0; i < NumberOfkPoints; i++) {
	k = exp(k1 + i*delk);
	PSLookUpTable[species][i] = EvaluatePowerSpectrum(k, species);
      }
  }
 
  /* Output power spectrum */
 
  FILE *fptr;
  //  if ((fptr = fopen("PowerSpectrum.out", "w")) == NULL) {
  if ((fptr = fopen(PowerSpectrumFilename, "w")) == NULL) {
    fprintf(stderr, "Error opening %s\n",PowerSpectrumFilename);
    return FAIL;
  }
  for (i = 0; i < NumberOfkPoints; i += 50) {
    fprintf(fptr, "%"GSYM" ", exp(k1 + i*delk));
    for (species = 0; species < 3; species++)
      fprintf(fptr, "%"GSYM" ", PSLookUpTable[species][i]);
    fprintf(fptr, "\n");
  }
  fclose(fptr);
 
  return SUCCESS;
}
 
 
 
FLOAT ps_tophat(FLOAT *kpntr)
{
 
  FLOAT psval, x, tophat, k = *kpntr;
 
  psval  = EvaluatePowerSpectrum(k, 0);
  x      = TophatRadius * k;
  tophat = 3.0/(x*x*x)*(sin(x) - x*cos(x));
  return 4.0*3.14159 * psval * (k*k*tophat*tophat);
}
 

/***********************************************************************
/
/  GRID CLASS (COMPUTE ELECTRON FRACTION ESTIMATE)
/
/  written by: John H. Wise
/  date:       April, 2007
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#define EFRAC_LOWERLIMIT 1e-2
#define MIN_TEMP 50
#define MAX_IONTEMP 30000

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

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::ElectronFractionEstimate(float dt)
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Find necessary fields. */

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum, gammaNum, kphHMNum, kdissH2IINum;

  IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
			     Vel3Num, TENum);

  IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
			HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, 
        VelocityUnits, TimeUnits, aUnits = 1;
  FLOAT a = 1.0, dadt;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {

    if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }

    aUnits = 1.0/(1.0 + InitialRedshift);

  }

  float afloat = float(a);

  /* For cells with photo-ionization rates and low e-fractions (shell
     of HII regions), estimate e-fraction. */

  float mh = 1.673e-24;
  float alpha_recombination = 2.59e-13;  // Assume T = 1e4 K
  alpha_recombination *= (TimeUnits * (DensityUnits / mh));

  int i, j, k, index;
  float efrac, t_i, x_eq, x_estimate;
  float total, total_h, t_frac, new_hii;

  // Compton cooling
  float comp1 = 1e-20, comp2 = 1e-20, zr;
  if (ComovingCoordinates) {
    zr = 1.0 / (afloat * aUnits) - 1.0;
    comp1 = CoolData.comp * (1.0 + zr);
    comp2 = 2.73 * (1.0 + zr);
  }

  double CoolUnit, xbase1, dbase1;
  xbase1 = LengthUnits / (afloat * aUnits);
  dbase1 = DensityUnits * pow((afloat*aUnits), 3);
  CoolUnit = (pow(aUnits,5) * pow(xbase1,2) * pow(mh,2)) /
    (pow(TimeUnits,3) * dbase1);
  float a3 = afloat*afloat*afloat;
  float dom = DensityUnits * a3 / mh;
  double rtunits = 1.602e-12 / TimeUnits / CoolUnit;

  float proper_d, proper_de, proper_hi, proper_hii, proper_hei, proper_heii,
    proper_heiii, pressure, temperature, max_edotplus;
  float logtem, logtem0, logtem9, dlogtem, t1, t2, tdef;
  float ceHI, ceHeI, ceHeII, ceHeIII, ciHI, ciHeI, ciHeII, ciHeIS, reHII,
    reHeII1, reHeII2, reHeIII;
  float edot, edotplus, brem;
  int cindex;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	efrac = BaryonField[DeNum][index] / BaryonField[DensNum][index];

	if (BaryonField[kphHINum][index] > 0 && efrac < EFRAC_LOWERLIMIT) {

	  /* Make an estimate of temperature only using the 6-species
	     model (H2 cooling) will be unimportant in the shells). */

#ifdef UNUSED
	  proper_d     = BaryonField[DensNum][index] 	/ a3;
	  proper_de    = BaryonField[DeNum][index] 	/ a3;
	  proper_hi    = BaryonField[HINum][index] 	/ a3;
	  proper_hii   = BaryonField[HIINum][index] 	/ a3;
	  proper_hei   = BaryonField[HeINum][index] 	/ a3;
	  proper_heii  = BaryonField[HeIINum][index] 	/ a3;
	  proper_heiii = BaryonField[HeIIINum][index] 	/ a3;

	  if (DualEnergyFormalism)
	    pressure = (Gamma - 1.0) * proper_d * BaryonField[GENum][index];
	  else {
	    pressure = BaryonField[TENum][index]
	      - 0.5 * (BaryonField[Vel1Num][index]*BaryonField[Vel1Num][index]
		       +BaryonField[Vel2Num][index]*BaryonField[Vel2Num][index]
		       +BaryonField[Vel3Num][index]*BaryonField[Vel3Num][index]);
	    pressure = max((Gamma - 1.0) * proper_d * pressure, tiny_number);
	  }

	  temperature = 0.25*(proper_hei + proper_heii + proper_heiii) 
	    + proper_hi + proper_hii + proper_de;
	  max_edotplus = (temperature * MAX_IONTEMP) 
	    / ((Gamma - 1.0) * proper_d * TemperatureUnits);
	  temperature = max(pressure*TemperatureUnits/temperature, MIN_TEMP);

	  logtem = log(temperature);
	  logtem0 = log(CoolData.TemperatureStart);
	  logtem9 = log(CoolData.TemperatureEnd);
	  dlogtem = (logtem9 - logtem0) / 
	    float(CoolData.NumberOfTemperatureBins);
	  cindex = (logtem - logtem0) / dlogtem + 1;

	  t1 = logtem0 + (cindex - 1) * dlogtem;
	  t2 = logtem0 + (cindex    ) * dlogtem;
	  tdef = 1.0 / (t2-t1);

	  // Lookup cooling rates and linearly interpolate
	  ceHI = CoolData.ceHI[cindex] + (logtem - t1) * 
	    (CoolData.ceHI[cindex+1] - CoolData.ceHI[cindex]) * tdef;
	  ceHeI = CoolData.ceHeI[cindex] + (logtem - t1) * 
	    (CoolData.ceHeI[cindex+1] - CoolData.ceHeI[cindex]) * tdef;
	  ceHeII = CoolData.ceHeII[cindex] + (logtem - t1) * 
	    (CoolData.ceHeII[cindex+1] - CoolData.ceHeII[cindex]) * tdef;
	  ciHI = CoolData.ciHI[cindex] + (logtem - t1) * 
	    (CoolData.ciHI[cindex+1] - CoolData.ciHI[cindex]) * tdef;
	  ciHeI = CoolData.ciHeI[cindex] + (logtem - t1) * 
	    (CoolData.ciHeI[cindex+1] - CoolData.ciHeI[cindex]) * tdef;
	  ciHeIS = CoolData.ciHeIS[cindex] + (logtem - t1) * 
	    (CoolData.ciHeIS[cindex+1] - CoolData.ciHeIS[cindex]) * tdef;
	  ciHeII = CoolData.ciHeII[cindex] + (logtem - t1) * 
	    (CoolData.ciHeII[cindex+1] - CoolData.ciHeII[cindex]) * tdef;
	  reHII = CoolData.reHII[cindex] + (logtem - t1) * 
	    (CoolData.reHII[cindex+1] - CoolData.reHII[cindex]) * tdef;
	  reHeII1 = CoolData.reHeII1[cindex] + (logtem - t1) * 
	    (CoolData.reHeII1[cindex+1] - CoolData.reHeII1[cindex]) * tdef;
	  reHeII2 = CoolData.reHeII2[cindex] + (logtem - t1) * 
	    (CoolData.reHeII2[cindex+1] - CoolData.reHeII2[cindex]) * tdef;
	  reHeIII = CoolData.reHeIII[cindex] + (logtem - t1) * 
	    (CoolData.reHeIII[cindex+1] - CoolData.reHeIII[cindex]) * tdef;
	  brem = CoolData.brem[cindex] + (logtem - t1) * 
	    (CoolData.brem[cindex+1] - CoolData.brem[cindex]) * tdef;

	  edot = 
	    - ceHI * proper_hi * proper_de
	    - ceHeI * proper_heii * proper_de * proper_de * dom * 0.25
	    - ceHeII * proper_heii * proper_de * 0.25
	    - ciHI * proper_hi * proper_de
	    - ciHeI * proper_hei * proper_de * 0.25
	    - ciHeII * proper_heii * proper_de * 0.25
	    - ciHeIS * proper_heii * proper_de * proper_de * dom * 0.25
	    - reHII * proper_hii * proper_de
	    - reHeII1 * proper_heii * proper_de * 0.25
	    - reHeII2 * proper_heii * proper_de * 0.25
	    - reHeIII * proper_heiii * proper_de * 0.25
	    - comp1 * (temperature - comp2) * proper_de / dom
	    - CoolData.comp_xray * (temperature - CoolData.temp_xray)
	      * proper_de / dom
	    - brem * (proper_hii + 0.25*proper_heii + proper_heiii) * proper_de;

	  //	  printf("edot[0] = %"GSYM"\n", edot);

	  edotplus = CoolData.ipiht * BaryonField[gammaNum][index] * rtunits
	    * proper_hi / dom;
	  edotplus = min(edotplus, max_edotplus);

	  edot += edotplus;

	  //	  printf("edotplus = %"GSYM" (%"GSYM")\n", edotplus, max_edotplus);

	  if (DualEnergyFormalism)

	    BaryonField[GENum][index] += edot / proper_d * dt;
	  BaryonField[TENum][index] += edot / proper_d * dt;
#endif /* UNUSED */

//	  printf("(%"ISYM" %"ISYM" %"ISYM") tem=%.2g, gg=%.2g, edot=%.2g, "
//		 "ge=%.2g\n",
//		 i, j, k, temperature, BaryonField[gammaNum][index],
//		 edot, BaryonField[GENum][index]);

	  /* Estimate ionization fraction */

	  total = BaryonField[kphHINum][index] + 
	    (alpha_recombination * BaryonField[DeNum][index]);
	  t_i = 1.0 / total;
	  x_eq = BaryonField[kphHINum][index] / total;
	  t_frac = dt / t_i;
	  x_estimate = x_eq + (efrac - x_eq) * (1 - exp(-t_frac)) / t_frac;

	  /* Correct electron and hydrogen species for this estimate */

	  total_h = BaryonField[HINum][index] + BaryonField[HIINum][index];
	  new_hii = x_eq * total_h;

//	  printf("(%"ISYM" %"ISYM" %"ISYM") t_i=%.2g, x_eq=%.2g, t_frac=%.2g, "
//		 "x_est=%.2g, x0=%.2g, kph=%.2g\n",
//		 i, j, k, t_i, x_eq, t_frac, x_estimate, efrac, 
//		 BaryonField[kphHINum][index]);

	  BaryonField[DeNum][index] += new_hii - BaryonField[HIINum][index];
	  BaryonField[HINum][index] = (1.0 - x_eq) * total_h;
	  BaryonField[HIINum][index] = new_hii;

	} // ENDIF shell cell
	
      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  return SUCCESS;

}

/***********************************************************************
/
/  INITIALIZE HII REGION AROUND SOURCE
/
/  written by: John H. Wise
/  date:       January, 2006
/  modified1:  
/
/ PURPOSE: To avoid geometrical artifacts in the HII region, we
/          initialize the region around the source with an ionized
/          medium.
/
***********************************************************************/

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

/* function prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

#define MAX_RADIUS 5

int grid::InitializeSource(RadiationSourceEntry *RS)
{

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }

  /* Get units. */

  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits; 

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  int pos[3], rmin[3], radius = 1000000, dim, i, j, k;;

  /* Determine the radius.  The sphere shouldn't go outside of the host grid */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {

    pos[dim] = (int) ((RS->Position[dim] - GridLeftEdge[dim]) / 
		      CellWidth[dim][0]) + GridStartIndex[dim];
    rmin[dim] = min(pos[dim]-GridStartIndex[dim], GridEndIndex[dim]-pos[dim]);
    radius = min(radius, rmin[dim]);

  } /* ENDFOR dim */

  radius = min(radius, MAX_RADIUS);
  printf("InitializeSource: radius = %"ISYM"\n", radius);

  /* Select the correct field if we're using the coupled transfer/rate solver */

  float *density, *HI, *HII, *f_e, *te, *ge;

  //  if (RadiativeTransferCoupledRateSolver) {
//  if (0) {
//    density = NewBaryonField[DensNum];
//    HI      = NewBaryonField[HINum];
//    HII     = NewBaryonField[HIINum];
//    f_e     = NewBaryonField[DeNum];
//    te      = NewBaryonField[TENum];
//    ge      = NewBaryonField[GENum];
//  } else {
    density = BaryonField[DensNum];
    HI      = BaryonField[HINum];
    HII     = BaryonField[HIINum];
    f_e     = BaryonField[DeNum];
    te      = BaryonField[TENum];
    ge      = BaryonField[GENum];
    //  }

  float HIRecombinationRate = 2.6e-13;
  float HICrossSection = 6.8e-18;

  int delk, delj, deli, index;
  FLOAT delr, delCell;
  double Luminosity= (double) RS->Luminosity * (double) pow(LengthUnits,3) / 
    (double) TimeUnits;
  double ConvertToProperDensity = double(DensityUnits)/double(1.67e-24);
  float fH, fHI, Old_fHII;

  for (k = pos[2]-radius; k <= pos[2]+radius; k++) {
    delk = k - pos[2];
    for (j = pos[1]-radius; j <= pos[1]+radius; j++) {
      delj = j - pos[1];
      index = (GridDimension[1]*k + j) * GridDimension[0];
      for (i = pos[0]-radius; i <= pos[0]+radius; i++) {

	deli = i - pos[0];
	delCell = sqrt(delk*delk + delj*delj + deli*deli);
	delCell = max(delCell, 0.1);
	delr = delCell * CellWidth[0][0] * LengthUnits;  // in cm
	
	if (delCell > radius) continue;

	/* f_HI = (n/Flux) * (k_rec / sigma) :: Equation (5.10) in
	   "This Physics of the ISM" by Dyson & Williams */

	fH = CoolData.HydrogenFractionByMass;
	fHI = (ConvertToProperDensity * 
	       (double) density[index+i] * fH / 
	       (Luminosity / (4*M_PI*(double)pow(delr,2)))) *
	  (HIRecombinationRate / HICrossSection);

//	printf("rho = %"GSYM", SA = %"GSYM", ratio = %"GSYM"\n", 
//	       ConvertToProperDensity * (double)BaryonField[DensNum][index+i]*fH,
//	       (4*M_PI*(double)pow(delr,2)),
//	       (HIRecombinationRate / HICrossSection));
//	printf("fHI = %"GSYM", delr = %"GSYM", delC = %"GSYM", L = %"GSYM"\n", fHI, delr, delCell,
//	       Luminosity);
//	printf("delta(i,j,k) = %"ISYM" %"ISYM" %"ISYM"\n", deli, delj, delk);	

	Old_fHII = HII[index+i];

	HI[index+i ] = fHI      * density[index+i];
	HII[index+i] = (fH-fHI) * density[index+i];
	f_e[index+i]+= HII[index+i] - Old_fHII;

//	te[index+i] = 1e4 / TemperatureUnits;
//
//	if (DualEnergyFormalism)

//	  ge[index+i] = te[index+i];

      } /* ENDFOR i */
    } /* ENDFOR j */
  } /* ENDFOR k */

  return SUCCESS;

}

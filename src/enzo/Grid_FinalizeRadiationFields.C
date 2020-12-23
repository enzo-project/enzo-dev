#define DEBUG 0
#define MYPROC MyProcessorNumber == ProcessorNumber

/***********************************************************************
/
/  GRID CLASS (FINALIZE THE PHOTO-RATES BY DIVIDING BY NUMBER OF PARTICLES)
/
/  written by: John Wise
/  date:       June, 2009
/  modified1:
/
/  PURPOSE:
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
#include "fortran.def"
#include "Grid.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int FindField(int f, int farray[], int n);

float ComputeHeatingRateFromDustModel(const float &n_H, const float &n_e,
                                      const float &Z,
                                      const float &G, const float &dx);



int grid::FinalizeRadiationFields(void)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

#ifdef TRANSFER

  int i, j, k, index, dim;
  float CellVolume = 1;

  for (dim = 0; dim < GridRank; dim++)
    CellVolume *= CellWidth[dim][0];

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }

  /* Find radiative transfer fields. */

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum,
				  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  /* Get units. */

  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits,
    DensityUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, PhotonTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  float DensityConversion = DensityUnits / mh;
  float factor = DensityConversion * CellVolume;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	BaryonField[kphHINum][index] /= factor * BaryonField[HINum][index];
	BaryonField[gammaNum][index] /= factor * BaryonField[HINum][index]; //divide by N_HI = n_HI*(dx)^3
      } // ENDFOR i
    } // ENDFOR j

  if (RadiativeTransferHydrogenOnly == FALSE)
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  BaryonField[kphHeINum][index] /=
	    0.25 * factor * BaryonField[HeINum][index];
	  BaryonField[kphHeIINum][index] /=
	    0.25 * factor * BaryonField[HeIINum][index];
	} // ENDFOR i
      } // ENDFOR j

   if (IndividualStarFUVHeating && !RadiativeTransferOpticallyThinFUV){
     const int FUVRateNum = FindField(FUVRate, this->FieldType, this->NumberOfBaryonFields);
     const int PeNum      = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);
     const int DensNum    = FindField(Density, this->FieldType, this->NumberOfBaryonFields);
     const int MetalNum   = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);
     const double MassUnits = DensityUnits*LengthUnits*LengthUnits*LengthUnits;

     float *temperature;
     int size = 1;
     for (dim=0;dim<GridRank;dim++) size *= this->GridDimension[dim];
     temperature = new float[size];

     if(this->ComputeTemperatureField(temperature) == FAIL){
         ENZO_FAIL("Error in compute temperature in Grid_GrackleWrapper");
     }

     const float EnergyUnits = (MassUnits*VelocityUnits*VelocityUnits);
     const float EnergyUnits_inv = 1.0 / EnergyUnits;

     const float FluxConv = EnergyUnits / TimeUnits / LengthUnits / LengthUnits;
     const float FluxConv_inv = 1.0 / FluxConv;
     const float PeConversion = 1.0 / (EnergyUnits/TimeUnits/POW(LengthUnits,3)); 


     for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++){
       for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++){
         index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
         for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++){

           BaryonField[FUVRateNum][index] *= (EnergyUnits_inv * POW(LengthUnits,3));

           float n_H, n_e, Z;

           n_H = BaryonField[HINum][index]+BaryonField[HIINum][index];
           if(MultiSpecies>1){
	     n_H += BaryonField[HMNum][index]+
	            0.5*(BaryonField[H2INum][index]+BaryonField[H2IINum][index]);
           }
           n_H *= DensityConversion;
           n_e  = BaryonField[DeNum][index]*DensityUnits/me;
           Z    = BaryonField[MetalNum][index]/BaryonField[DensNum][index];

//         This is now done in initialization of RT fields
//           BaryonField[FUVRateNum][index] += (G_background * FluxConv_inv);

	   if (temperature[index] >= IndividualStarFUVTemperatureCutoff){
	     BaryonField[PeNum][index] = 0.0;
           } else {
             BaryonField[PeNum][index] += ComputeHeatingRateFromDustModel(n_H, n_e, Z,
                                                                         BaryonField[FUVRateNum][index]*FluxConv,
                                                                         -1.0) * PeConversion; // set dx < 0 to turn off self-shielding approx
           }

         }
       }
     } // end k

     delete [] temperature;

   } // end FUV PE heating


   if (MultiSpecies > 1 && !RadiativeTransferFLD && !RadiativeTransferOpticallyThinH2)
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  if(!RadiativeTransferUseH2Shielding) {
	    BaryonField[kdissH2INum][index] /=
	     1.0 * factor * BaryonField[H2INum][index];
	  }
	  BaryonField[kphHMNum][index] /=
	    1.0 * factor * BaryonField[HMNum][index];
	  BaryonField[kdissH2IINum][index] /=
	    1.0 * factor * BaryonField[H2IINum][index];
	} // ENDFOR i
      } // ENDFOR j

   if(this->IndexOfMaximumkph >= 0)
     this->MaximumkphIfront /= (factor * BaryonField[HINum][IndexOfMaximumkph]);

#endif /* TRANSFER */

  return SUCCESS;
}

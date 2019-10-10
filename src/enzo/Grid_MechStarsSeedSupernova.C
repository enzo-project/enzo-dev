/*
    If the region can form stars but has no metals, 
    We seed the region with a Pop III supernova.
    This takes mass from the host cell, converting it
    into ejecta mass, energy, and ejecta metals and 
    deposits it using the MechStars infrastructure.

    07/2019: Azton Wells
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#include "phys_constants.h"


    int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, float Time);
extern "C"  void FORTRAN_NAME(cic_deposit)(float* xPosition, float* yPosition,
        float* zPosition, int* gridRank, int* nParticles, 
        float* DepositQuantity, float* FieldToDepositTo,
        float* leftEdge, int* xGridDim, int* yGridDim, 
        int* zGridDim, float* gridDx, float* cloudsize);
int search_lower_bound(float *arr, float value, int low, int high, 
		       int total);
unsigned_long_int mt_random(void);
int StarParticlePopIII_IMFInitialize(void);

int grid::MechStars_SeedSupernova(float* totalMetal, float* temperature, int* seedIndex){
    debug = true;
    /* Initialize the IMF lookup table if requested and not defined */
    if (debug) fprintf(stdout, "setting IMF\n");
    if (PopIIIInitialMassFunction)
    StarParticlePopIII_IMFInitialize();
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    /* Find fields: density, total energy, velocity1-3. */
    int CRNum;
    if (debug) printf("IMF Set!\n");
    /* Find Multi-species fields. */
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
        DINum, DIINum, HDINum; 

    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                        Vel3Num, TENum) == FAIL) {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        return FAIL;
    }
    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
        TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                &TimeUnits, &VelocityUnits, &MassUnits, this->Time) == FAIL) {
        fprintf(stderr, "Error in GetUnits.\n");
        return FAIL;
    }
    int size = 1;
    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];

    /* set feedback to random cell in grid*/
    int ip = seedIndex[0];//rand() % (GridDimension[0]-2*NumberOfGhostZones)+NumberOfGhostZones ;
    int jp = seedIndex[1];//rand() % (GridDimension[1]-2*NumberOfGhostZones)+NumberOfGhostZones;
    int kp = seedIndex[2];//rand() % (GridDimension[2]-2*NumberOfGhostZones)+NumberOfGhostZones;
    float position[3] = {((float)ip+0.5)*CellWidth[0][0]+CellLeftEdge[0][0], 
                        ((float)jp+0.5)*CellWidth[0][0]+CellLeftEdge[1][0],
                        ((float)kp+0.5)*CellWidth[0][0]+CellLeftEdge[2][0]};
    int index = ip + jp*GridDimension[0]+kp*GridDimension[0]*GridDimension[1];
    /* Get mass of star in exact method used by Star_AssignFinalMassFromIMF.C */
    if (debug) fprintf (stdout, "Setting final mass\n");
    unsigned_long_int random_int = mt_random();
    const int max_random = (1<<16);
    float x = (float) (random_int%max_random) / (float) (max_random);
    float dm = log10(PopIIIUpperMassCutoff / PopIIILowerMassCutoff) / 
        (float) (IMF_TABLE_ENTRIES-1);

    /* (binary) search for the mass bin corresponding to the random
        number */

    int width = IMF_TABLE_ENTRIES/2;
    int bin_number = IMF_TABLE_ENTRIES/2;
    
    while (width > 1) {
        width /= 2;
        if (x > IMFData[bin_number])
        bin_number += width;
        else if (x < IMFData[bin_number])
        bin_number -= width;
        else
        break;
    }
  
    float FinalMass = PopIIILowerMassCutoff * POW(10.0, bin_number * dm);
    /* certain mass have no effect since they collapse to black holes */
    if (FinalMass < 10 || (FinalMass >40 && FinalMass < 140) || (FinalMass > 260)){
        printf("Mass Outside Supernova Range!\n");
        return SUCCESS;
    }
    if (FinalMass > 0.1*BaryonField[DensNum][index]*MassUnits){
        printf("Cell too small for PIII star!");
        return FAIL;
    }
    /* Now, calculate feedback parameters as in Star_CalculateFeedbackParameters.C */

    // parameters of supernovae //
    const float TypeIILowerMass = 11, TypeIIUpperMass = 20.;
    const float HNLowerMass = 20.1, HNUpperMass = 40.1;
    const float PISNLowerMass = 140, PISNUpperMass = 260;

    // From Nomoto et al. (2006)
    const float HypernovaMetals[] = {3.36, 3.53, 5.48, 7.03, 8.59}; // SolarMass
    const float HypernovaEnergy[] = {10, 10, 20, 25, 30}; // 1e51 erg 
    const float CoreCollapseMetals[] = {3.63, 4.41, 6.71, 8.95, 11.19}; // SolarMass
    const float CoreCollapseEnergy[] = {1, 1, 1, 1, 1}; // 1e51 erg

    const float SNExplosionMass[] = {19.99, 25, 30, 35, 40.01};  // SolarMass
    const float *SNExplosionMetals = (PopIIIUseHypernova ==TRUE) ? 
        HypernovaMetals : CoreCollapseMetals;
    const float *SNExplosionEnergy = (PopIIIUseHypernova ==TRUE) ? 
        HypernovaEnergy : CoreCollapseEnergy;
    float HeCoreMass = 0;
    float SNEnergy = 0;
    float EjectaMetal = 0;
    float EjectaMass = FinalMass;
    // if (BaryonField[DensNum][index]*DensityUnits*pow(LengthUnits*CellWidth[0][0],3) < 5.0*FinalMass){
    //     /* Should probably remove from larger area if host cell is too small,
    //         but should only matter at VERY high resolution: M_cell < 50 M_sun*/
    //     if (debug) fprintf(stdout, "Not enough mass in cell for Seed Supernova: %f < %f", 
    //             BaryonField[DensNum][index]*DensityUnits*pow(LengthUnits*CellWidth[0][0],3), FinalMass);
    //     return FAIL;
    //     }
    /* Reverse CIC out the star mass */
    int np = 1; float MassRem = -1*EjectaMass;
    float LeftEdge[3] = {CellLeftEdge[0][0], CellLeftEdge[1][0], CellLeftEdge[2][0]};
    float cloudSize = CellWidth[0][0];
    if (debug) fprintf (stdout, "Removing PIII mass from grid\n");
    FORTRAN_NAME(cic_deposit)(&position[0], &position[1], &position[2],
        &GridRank,&np,&MassRem, BaryonField[DensNum], LeftEdge, 
        &GridDimension[0], &GridDimension[1], &GridDimension[2], 
        &CellWidth[0][0], &cloudSize);
    if (debug) fprintf(stdout, "PIII Mass: %f\n", FinalMass);
    if (FinalMass > PISNLowerMass && FinalMass < PISNUpperMass){
        HeCoreMass = (13./24.)*(FinalMass-20);
        SNEnergy = (5.0+1.304*(HeCoreMass-64))*1e51;
        EjectaMetal = HeCoreMass;
    }
    if (FinalMass > TypeIILowerMass && FinalMass < TypeIIUpperMass){
        SNEnergy = 1e51;
        EjectaMetal = 0.1077+0.3383*(FinalMass-11);
    }
    if (FinalMass > HNLowerMass && FinalMass < HNUpperMass){
        int bin = search_lower_bound((float*)SNExplosionMass, FinalMass, 0, 5, 5);
     	float frac = (SNExplosionMass[bin+1] - FinalMass) / 
	            (SNExplosionMass[bin+1] - SNExplosionMass[bin]);
	    SNEnergy = 1e51 * (SNExplosionEnergy[bin] + 
			   frac * (SNExplosionEnergy[bin+1] - SNExplosionEnergy[bin]));
        EjectaMetal = (SNExplosionMetals[bin] + 
		     frac * (SNExplosionMetals[bin+1] - SNExplosionMetals[bin]));
    }
    /* Need mu_field for deposition routine */

    float mu_field [size];
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	  
	         index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
	         mu_field[index] = 0.0;
	         // calculate mu

	         if (MultiSpecies == 0) {
	            mu_field[index] = Mu;
	         } else {

	            if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
				         HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
	               ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
	            }

	         mu_field[index] = BaryonField[DeNum][index] + BaryonField[HINum][index] + BaryonField[HIINum][index] +
	               (BaryonField[HeINum][index] + BaryonField[HeIINum][index] + BaryonField[HeIIINum][index])/4.0;
	         if (MultiSpecies > 1) {
	            mu_field[index] += BaryonField[HMNum][index] + (BaryonField[H2INum][index] + BaryonField[H2IINum][index])/2.0;
	         }
	         if (MultiSpecies > 2) {
	            mu_field[index] += (BaryonField[DINum][index] + BaryonField[DIINum][index])/2.0 + (BaryonField[HDINum][index]/3.0);
	         }
	    
	      }
	   }
   }
   }
    /*  Add this to the grid using MechStars_DepositFeedback */
    float vp=0, up=0, wp=0;
    if (debug) fprintf(stdout, "Calling DepositFeedback!\n");
    this->MechStars_DepositFeedback(SNEnergy, EjectaMass, EjectaMetal, totalMetal,temperature,
                            &up, &vp, &wp, &position[0], &position[1], &position[2],
                            ip, jp, kp, size, mu_field, 0, 0, 0, 0, 1);
    
    return SUCCESS;
}
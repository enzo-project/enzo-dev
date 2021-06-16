/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A STAR PARTICLE TEST)
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::TestStarParticleInitializeGrid(float TestStarParticleStarMass,
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
					 int NumberOfTestStars, float clusterRadius, 
			     char *TestStarInitializationFilename)
{
  /* declarations */

  float CentralMass = 1.0;

  int i, dim, j, k, size, active_size, index, cindex, n;
  float TestInitialdt = *Initialdt;
  int DensNum, TENum, GENum, DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum, Vel1Num, Vel2Num, Vel3Num; 
  float *density_field = NULL, *electron_density_field = NULL,
          *hi_field = NULL, *hii_field = NULL, *hei_field = NULL,
          *heii_field = NULL, *heiii_field = NULL, *h2i_field = NULL,
          *h2ii_field = NULL, *hm_field = NULL, *z_field = NULL,
          *v1_field = NULL, *v2_field = NULL, *v3_field = NULL,
          *ge_field = NULL;
  int gz = NumberOfGhostZones;
  
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;


  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;
  double MassUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }


  if (TestStarInitializationFilename){
    /* Initialize a grid based on an input h5 file. Mostly lifted from Grid_PhotonTestInitializeGrid.C - AIW*/
    active_size = 1;
    int ActiveDims[MAX_DIMENSION];
    for (dim = 0; dim < GridRank; dim++) {
      ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
      active_size *= ActiveDims[dim];
    }



  /* create fields */
    NumberOfBaryonFields = 0;
    FieldType[DensNum = NumberOfBaryonFields++] = Density;
    FieldType[TENum = NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[GENum = NumberOfBaryonFields++] = InternalEnergy;
    int ivel = NumberOfBaryonFields;
    FieldType[Vel1Num = NumberOfBaryonFields++] = Velocity1;
    FieldType[Vel2Num = NumberOfBaryonFields++] = Velocity2;
    FieldType[Vel3Num = NumberOfBaryonFields++] = Velocity3;
    if (MultiSpecies) {
      FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
      FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
      FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
      FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
      FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
      FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
      if (MultiSpecies > 1) {
        FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
        FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
        FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
      }
      if (MultiSpecies > 2) {
        FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
        FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
        FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
      }
    }
     if (MetalCooling)
    FieldType[MetalNum    = NumberOfBaryonFields++] = Metallicity; 

    
    this->AllocateGrids();

    // Read Fields
    char *filename;
    hsize_t OutDims[MAX_DIMENSION];
    herr_t h5error;
    hid_t file_id;
    char *delim = "/";
    /* fields we load with this method.  The HDF5 file needs to \
    have one dataset for each field, named as follows.  This method is 
    intended for multispecies = 2 only at the moment -AIW */

    char *density = "Density";
    char *edensity = "Electron_Density";
    char *ge = "GasEnergy";
    char *hII_density = "HII_Density";
    char *hI_density = "HI_Density";
    char *heI_density = "HeI_Density";
    char *heII_density = "HeII_Density";
    char *heIII_density = "HeIII_Density";
    char *h2I_density = "H2I_Density";
    char *h2II_density = "H2II_Density";
    char *metal_density = "Metal_Density";
    char *hm_density = "HM_Density";
    char *velocity1 = "x-velocity";
    char *velocity2 = "y-velocity";
    char *velocity3 = "z-velocity";
  

    filename = strtok(TestStarInitializationFilename, delim);
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      OutDims[GridRank-dim-1] = GridEndIndex[dim] - GridStartIndex[dim] + 1;

    if (file_id == -1) ENZO_FAIL("Error opening field file.");
    
    density_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, density, file_id,
		       HDF5_REAL, BaryonField[DensNum], FALSE, NULL, NULL);

    electron_density_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, edensity, file_id,
                      HDF5_REAL, BaryonField[DeNum], FALSE, NULL, NULL);

    ge_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, ge, file_id,
                      HDF5_REAL, BaryonField[GENum], FALSE, NULL, NULL);

    hi_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, hI_density, file_id,
                      HDF5_REAL, BaryonField[HINum], FALSE, NULL, NULL);

    hii_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, hII_density, file_id,
                      HDF5_REAL, BaryonField[HIINum], FALSE, NULL, NULL);

    hei_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, heI_density, file_id,
                      HDF5_REAL, BaryonField[HeINum], FALSE, NULL, NULL);

    heii_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, heII_density, file_id,
                      HDF5_REAL, BaryonField[HeIINum], FALSE, NULL, NULL);
    
    heiii_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, heIII_density, file_id,
                      HDF5_REAL, BaryonField[HeIIINum], FALSE, NULL, NULL);

    h2i_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, h2I_density, file_id,
                      HDF5_REAL, BaryonField[H2INum], FALSE, NULL, NULL);

    h2ii_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, h2II_density, file_id,
                      HDF5_REAL, BaryonField[H2IINum], FALSE, NULL, NULL);

    hm_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, hm_density, file_id,
                      HDF5_REAL, BaryonField[HMNum], FALSE, NULL, NULL);

    z_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, metal_density, file_id,
                      HDF5_REAL, BaryonField[MetalNum], FALSE, NULL, NULL);

    v1_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, velocity1, file_id,
                      HDF5_REAL, BaryonField[Vel1Num], FALSE, NULL, NULL);

    v2_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, velocity2, file_id,
                      HDF5_REAL, BaryonField[Vel2Num], FALSE, NULL, NULL);

    v3_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, velocity3, file_id,
                      HDF5_REAL, BaryonField[Vel3Num], FALSE, NULL, NULL);

    h5error = H5Fclose(file_id);
    if (h5error == -1) ENZO_FAIL("Error closing initial conditions file.");
  printf("TestStarParticleInitialize\n");
  printf("Active size = %d, OutDim = %d, gz = %d, Grid dim = %ld\n", ActiveDims[0], OutDims[0], gz, 
                (GridDimension[0]));
  fflush(stdout);
  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {
        float hfrac = 0.75;
        float hefrac = 0.24;
        float h2frac = 0.01;
      	if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
	                j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
	                k >= GridStartIndex[2] && k <= GridEndIndex[2]) {
	        cindex = (i-GridLeftEdge[0]) + ActiveDims[0] *
	                  ((j-GridLeftEdge[1]) + (k-GridLeftEdge[2])*ActiveDims[1]);
          
          BaryonField[TENum][cindex] = BaryonField[GENum][cindex]; //0.5 * density_field[cindex] 
              // * (
              //   v1_field[cindex] * v1_field[cindex] 
              // + v2_field[cindex] * v2_field[cindex] 
              // + v3_field[cindex] * v3_field[cindex]
              //   ) 
              // + ge_field[cindex];
      //     BaryonField[DeNum][i] = BaryonField[HIINum][i] +
      //       0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];
      //           if (MultiSpecies > 1)
      //       BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] -
      //         BaryonField[HMNum][i];
      }
      }
  }
  printf("Central Mass: %f \n",CentralMass);

  /* Get Units. */


  /* Set Central Mass in simulation units */

  CentralMass = TestStarParticleStarMass*1.99e33* pow(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;


  /* Set number of particles for this grid and allocate space. */

  NumberOfParticles = NumberOfTestStars;
  NumberOfParticleAttributes = 4;
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles\n", NumberOfParticles);


  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    if (STARFEED_METHOD(POP3_STAR))
      ParticleType[i] = -1*PARTICLE_TYPE_SINGLE_STAR;
    else
      ParticleType[i] = PARTICLE_TYPE_STAR;
  }
  float p1;

  /* Set central particle. */
  for (i = 0; i <= NumberOfParticles; i++){
    for (dim = 0; dim < GridRank; dim++) {
      if (NumberOfParticles == 1){
        p1 = 0.5;
      }else{
        int rng = clusterRadius*200;
        p1 = float(rand() % rng)/100.0+(0.5-clusterRadius);
      }
        ParticlePosition[dim][i] = p1*
        (DomainLeftEdge[dim]+DomainRightEdge[dim]) + 0.5*CellWidth[0][0];
      ParticleVelocity[dim][i] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
    }
    ParticleMass[i] = CentralMass;
    ParticleAttribute[0][i] = Time+1e-7; //creation time:make sure it is non-zero
    if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][i] = 10.0 * Myr_s/TimeUnits;
    if (STARFEED_METHOD(MOM_STAR))
      if(StarMakerExplosionDelayTime >= 0.0)
        ParticleAttribute[1][i] = 1.0;
      else
        ParticleAttribute[1][i] =10.0 * Myr_s/TimeUnits;
    if (STARFEED_METHOD(MECHANICAL)) {
      if (StarParticleRadiativeFeedback){
        ParticleAttribute[1][i] = 25 * Myr_s/TimeUnits; // radiate for 25 Myr
      }
      else{
       ParticleAttribute[1][i] = 3.953;
      }
    }
  ParticleAttribute[2][i] = 0.0;  // Metal fraction
  ParticleAttribute[3][i] = 0.0;  // metalfSNIa
  }
  // delete [] density_field; 
  // delete [] electron_density_field;
  // delete [] hi_field;
  // delete [] hii_field; 
  // delete [] hei_field;
  // delete [] heii_field; 
  // delete [] heiii_field; 
  // delete [] h2i_field;
  // delete [] h2ii_field; 
  // delete [] hm_field; 
  // delete [] z_field;
  // delete [] v1_field; 
  // delete [] v2_field; 
  // delete [] v3_field;
  // delete [] ge_field;
  return SUCCESS;
}


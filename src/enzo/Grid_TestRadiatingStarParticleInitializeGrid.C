/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A RADIATING POP III STAR PARTICLE TEST)
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:  Danielle Skinner
/  date:       July, 2022
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

static int TestStarParticleParticleCount = 0;

int grid::TestRadiatingStarParticleInitializeGrid(float TestStarParticleStarMass, 
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
					 FLOAT TestStarParticleStarPosition[],
           float TestStarParticleDensity,
           float TestStarParticleEnergy,
           float TestStarParticleVelocity[])
{
  /* declarations */

  float CentralMass = 1.0, temperature;
  int dim, i, j, k, m, field, size, active_size, index, cindex;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum, MetalIaNum;

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kdissH2IINum, kphHMNum, RPresNum1, RPresNum2, RPresNum3;

  int ExtraField[2];

  float TestInitialdt = *Initialdt, *density_field = NULL, *HII_field = NULL, *HeII_field = NULL, 
    *HeIII_field = NULL, *Temperature_field = NULL;

   /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
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

  if (TestProblemData.UseMetallicityField){
    MetalNum = NumberOfBaryonFields;
    FieldType[NumberOfBaryonFields++] = Metallicity;
  }
  if (RadiativeTransfer && (MultiSpecies < 1)) {
    ENZO_FAIL("Grid_TestRadiatingStarParticleInitialize: Radiative Transfer but not MultiSpecies set");
  }

  //   Allocate fields for photo ionization and heating rates
  if (RadiativeTransfer)
    if (MultiSpecies) {
      FieldType[kphHINum    = NumberOfBaryonFields++] = kphHI;
      FieldType[gammaNum    = NumberOfBaryonFields++] = PhotoGamma;
      if (RadiativeTransferHydrogenOnly == FALSE) {
        FieldType[kphHeINum   = NumberOfBaryonFields++] = kphHeI;
        FieldType[kphHeIINum  = NumberOfBaryonFields++] = kphHeII;
      }
      if (MultiSpecies > 1) {
        FieldType[kdissH2INum     = NumberOfBaryonFields++] = kdissH2I;
        FieldType[kdissH2IINum    = NumberOfBaryonFields++] = kdissH2II;
        FieldType[kphHMNum        = NumberOfBaryonFields++] = kphHM;
      }
    } 

  if (RadiationPressure && RadiativeTransfer) {
    FieldType[RPresNum1 = NumberOfBaryonFields++] = RadPressure0;
    FieldType[RPresNum2 = NumberOfBaryonFields++] = RadPressure1;
    FieldType[RPresNum3 = NumberOfBaryonFields++] = RadPressure2;
  }

  NumberOfPhotonPackages = 0;
  PhotonPackages-> NextPackage= NULL;
  
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Get Units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, CriticalDensity = 1, BoxLength = 1, mu = 0.6, mu_data;
  double MassUnits = 1;
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Set number of particles for this grid and allocate space. */

  NumberOfParticles = 1;
  NumberOfParticleAttributes = 4;
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles\n", NumberOfParticles);
 
  /* Set up the baryon field. */ 
   /* compute size of fields */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  /* allocate fields */

  this->AllocateGrids();  

  /* Initialize radiation fields */ // need this

  if (this->InitializeRadiativeTransferFields() == FAIL) {
    ENZO_FAIL("\nError in InitializeRadiativeTransferFields.\n");
  }

  for (i = 0; i < size; i++) {
    BaryonField[0][i] = TestStarParticleDensity;
    BaryonField[1][i] = TestStarParticleEnergy;
  }

  /* set velocities */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      BaryonField[vel+dim][i] = TestStarParticleVelocity[dim];

  /* Set internal energy if necessary. */
 
  if (DualEnergyFormalism)
    for (i = 0; i < size; i++)
      BaryonField[2][i] = BaryonField[1][i];

  for (i = 0; i < size; i++){

    if(MultiSpecies > 0){

      BaryonField[HIINum][i] = TestProblemData.HII_Fraction * 
	TestProblemData.HydrogenFractionByMass * TestStarParticleDensity;
 
      BaryonField[HeIINum][i] =  TestProblemData.HeII_Fraction *
	TestStarParticleDensity * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeIIINum][i] = TestProblemData.HeIII_Fraction *
	TestStarParticleDensity * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeINum][i] =
	(1.0 - TestProblemData.HydrogenFractionByMass)*TestStarParticleDensity -
	BaryonField[HeIINum][i] - BaryonField[HeIIINum][i];

      if(MultiSpecies > 1){
	BaryonField[HMNum][i] = TestProblemData.HM_Fraction *
	  BaryonField[HIINum][i];

	BaryonField[H2INum][i] = TestProblemData.H2I_Fraction *
	  BaryonField[0][i] * TestProblemData.HydrogenFractionByMass;

	BaryonField[H2IINum][i] = TestProblemData.H2II_Fraction * 2.0 *
	  BaryonField[HIINum][i];

    BaryonField[kdissH2INum][i] = 0.0; 
	  BaryonField[kphHMNum][i] = 0.0;
	  BaryonField[kdissH2IINum][i] = 0.0;
      }

      BaryonField[HINum][i] = TestProblemData.HydrogenFractionByMass*BaryonField[0][i]
	- BaryonField[HIINum][i];
      if (MultiSpecies > 1)
	BaryonField[HINum][i] -= (BaryonField[HMNum][i] + BaryonField[H2IINum][i]
				  + BaryonField[H2INum][i]);

      BaryonField[DeNum][i] = BaryonField[HIINum][i] +
	0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];
      if (MultiSpecies > 1)
	BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] -
	  BaryonField[HMNum][i];

      if(MultiSpecies > 2){
	BaryonField[DINum ][i]  = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HINum][i];
	BaryonField[DIINum][i] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HIINum][i];
	BaryonField[HDINum][i] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * BaryonField[H2INum][i];
      }

    } // if(TestProblemData.MultiSpecies)

    // metallicity fields
    if(TestProblemData.UseMetallicityField){
      BaryonField[MetalNum][i] = TestProblemData.MetallicityField_Fraction* TestStarParticleDensity;

      }

if (MultiSpecies) {
	  mu_data =  
	    0.25*(BaryonField[HeINum][i]  + BaryonField[HeIINum][i] +
		  BaryonField[HeIIINum][i]                        ) +
	    BaryonField[HINum][i]   + BaryonField[HIINum][i]  +
	    BaryonField[DeNum][i];
	  if (MultiSpecies > 1)
	    mu_data += BaryonField[HMNum][i]   +
	      0.5*(BaryonField[H2INum][i]  + BaryonField[H2IINum][i]);
	  mu_data = BaryonField[0][i] / mu_data;
	} else
	  mu_data = mu;


      } // for (i = 0; i < size; i++)

  /* Set Central Mass in simulation units */

  CentralMass = TestStarParticleStarMass*SolarMass* pow(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;
  
  printf("Central Mass: %f \n",CentralMass);

  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    ParticleType[i] = -PopIII;
    ParticleAttribute[0][i] = Time + 1e-7;
    ParticleAttribute[1][i] = StarMakerExplosionDelayTime;
  }
  ParticleMass[0] = CentralMass;

  /* Set central particle. */ 
  for (dim = 0; dim < GridRank; dim++) {
    ParticlePosition[dim][0] = TestStarParticleStarPosition[dim]*
      (DomainLeftEdge[dim]+DomainRightEdge[dim]) + 0.5*CellWidth[0][0];
    ParticleVelocity[dim][0] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
  }

  if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;
  if (STARFEED_METHOD(MOM_STAR) || STARFEED_METHOD(MECH_STAR))
    if(StarMakerExplosionDelayTime >= 0.0)
      ParticleAttribute[1][0] = 1.0;
    else
      ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;
  
  ParticleAttribute[2][0] = 0.0;  // Metal fraction
  ParticleAttribute[3][0] = 0.0;  // metalfSNIa

  delete [] density_field;
  delete [] HII_field;
  delete [] HeII_field;
  delete [] HeIII_field;
  delete [] Temperature_field;

  return SUCCESS;
}


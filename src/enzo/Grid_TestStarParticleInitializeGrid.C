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

// Used to compute Bonner-Ebert density profile
double ph_BE(double r);
double ph_q(double r);
double ph_Ang(double a1, double a2, double R, double r);

// Returns random velocity from Maxwellian distribution
double ph_Maxwellian(double c_tilda, double vel_unit, double mu, double gamma);

static int PhotonTestParticleCount = 0;

int grid::TestStarParticleInitializeGrid(float TestStarParticleStarMass, 
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
					 FLOAT TestStarParticleStarPosition[],
           float TestStarParticleDensity,
           float TestStarParticleEnergy,
           float TestStarParticleVelocity[],
           float TestStarParticleBField[],
           float InitialTemperature,
           float PhotonTestInitialFractionHII, 
			     float PhotonTestInitialFractionHeII,
			     float PhotonTestInitialFractionHeIII, 
			     float PhotonTestInitialFractionHM,
			     float PhotonTestInitialFractionH2I, 
			     float PhotonTestInitialFractionH2II)
{
  /* declarations */

  float CentralMass = 1.0;
  int dim, i, j, k, m, field, size, active_size, index, cindex;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum, MetalIaNum, B1Num, B2Num, B3Num, PhiNum, CRNum;

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kdissH2IINum, kphHMNum, RPresNum1, RPresNum2, RPresNum3;

  int ExtraField[2];

  float TestInitialdt = *Initialdt, *density_field = NULL, *HII_field = NULL, *HeII_field = NULL, 
    *HeIII_field = NULL, *Temperature_field = NULL;
  
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

   /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1 || HydroMethod > 2)
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2 || HydroMethod > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  if ( CRModel ) {
    CRNum = NumberOfBaryonFields;
    FieldType[NumberOfBaryonFields++] = CRDensity;
  }

  if (WritePotential)
    FieldType[NumberOfBaryonFields++] = GravPotential;


  int colorfields = NumberOfBaryonFields;

  // Enzo's standard multispecies (primordial chemistry - H, D, He)
  if (TestProblemData.MultiSpecies) {
    FieldType[DeNum     = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum     = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum    = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum    = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum   = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum  = NumberOfBaryonFields++] = HeIIIDensity;
    if (TestProblemData.MultiSpecies > 1) {
      FieldType[HMNum   = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum  = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (TestProblemData.MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }

  //  Metal fields, including the standard 'metallicity' as well 
  // as two extra fields
  if (TestProblemData.UseMetallicityField) {
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity;

    if (StarMakerTypeIaSNe)
      FieldType[MetalIaNum = NumberOfBaryonFields++] = MetalSNIaDensity;

    if(TestProblemData.MultiMetals){
      FieldType[ExtraField[0] = NumberOfBaryonFields++] = ExtraType0;
      FieldType[ExtraField[1] = NumberOfBaryonFields++] = ExtraType1;
    }
  }

  if (RadiativeTransfer && (MultiSpecies < 1)) {
    ENZO_FAIL("Grid_TestStarParticleInitialize: Radiative Transfer but not MultiSpecies set");
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
  
  /* Get Units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, CriticalDensity = 1, BoxLength = 1, mu = 0.6, mu_data;
  double MassUnits = 1;
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Set Central Mass in simulation units */

  CentralMass = TestStarParticleStarMass*1.99e33* pow(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;

  printf("Central Mass: %f \n",CentralMass);

  /* Set number of particles for this grid and allocate space. */

  NumberOfParticles = 1;
  NumberOfParticleAttributes = 4;
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles\n", NumberOfParticles);

   this->AllocateGrids();

  /* Initialize radiation fields */ // need this

  if (this->InitializeRadiativeTransferFields() == FAIL) {
    ENZO_FAIL("\nError in InitializeRadiativeTransferFields.\n");
  }

  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    ParticleType[i] = -PopIII;
    ParticleAttribute[0][i] = Time + 1e-7;
  }
  ParticleMass[0] = CentralMass;

  /* Set central particle. */ 
  for (dim = 0; dim < GridRank; dim++) {
    ParticlePosition[dim][0] = TestStarParticleStarPosition[dim]*
      (DomainLeftEdge[dim]+DomainRightEdge[dim]) + 0.5*CellWidth[0][0];
    ParticleVelocity[dim][0] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
  }

  if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;
  if (STARFEED_METHOD(MOM_STAR))
    if(StarMakerExplosionDelayTime >= 0.0)
      ParticleAttribute[1][0] = 1.0;
    else
      ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;
  
  ParticleAttribute[2][0] = 0.0;  // Metal fraction
  ParticleAttribute[3][0] = 0.0;  // metalfSNIa

    if (this->FindNewStarParticles(GridLevel) == FAIL) {
    ENZO_FAIL("Error in grid::FindNewStarParticles.");
    }

  /* Reset particle type to be positive*/
    for (i = 0; i < NumberOfParticles; i++) {
      ParticleNumber[i] = i;
      ParticleType[i] = PopIII;
    }
    Star *cstar;
    for (cstar = Stars; cstar; cstar = cstar->NextStar)
      cstar->type = PopIII; 


  /* Set up the baryon field. */ // need thid
  /* compute size of fields */
  active_size = 1;
  int ActiveDims[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    active_size *= ActiveDims[dim];
  }

  /* Loop over the mesh. */
  float density, dens1, Velocity[MAX_DIMENSION],
    temperature, temp1, sigma, sigma1, colour, outer_radius;
  float HII_Fraction, HeII_Fraction, HeIII_Fraction, H2I_Fraction;
  FLOAT r, x, y = 0, z = 0;
  int n = 0;

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {

	/* Compute position */ // need this 

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

  if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
	    j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
	    k >= GridStartIndex[2] && k <= GridEndIndex[2]) {
	  cindex = (i-GridStartIndex[0]) + ActiveDims[0] *
	    ((j-GridStartIndex[1]) + (k-GridStartIndex[2])*ActiveDims[1]);
	  if (density_field != NULL)
	    density = max(density_field[cindex], 1e-6);
	  else
	    density = 1.0;
	  if (HII_field != NULL)
	    HII_Fraction = HII_field[cindex];
	  else
	    HII_Fraction = PhotonTestInitialFractionHII;
	  if (HeII_field != NULL)
	    HeII_Fraction = HeII_field[cindex];
	  else
	    HeII_Fraction = PhotonTestInitialFractionHeII;
	  if (HeIII_field != NULL)
	    HeIII_Fraction = HeIII_field[cindex];
	  else
	    HeIII_Fraction = PhotonTestInitialFractionHeIII;
	  if (Temperature_field != NULL)
	    temperature = temp1 = Temperature_field[cindex];
	  else
	    temperature = temp1 = InitialTemperature;
	}

  /* set density, total energy */ // check
    BaryonField[0][n] = TestStarParticleDensity;
    BaryonField[1][n] = TestStarParticleEnergy;

  /* set velocities */ 

  for (dim = 0; dim < GridRank; dim++)
	  BaryonField[vel+dim][n] = Velocity[dim] + TestStarParticleVelocity[dim];

  /* Set internal energy if necessary. */ //check

  if (DualEnergyFormalism)
    BaryonField[2][n] = BaryonField[1][n];

	/* If doing multi-species (HI, etc.), set these. */

	if (MultiSpecies > 0) {
	  BaryonField[HIINum][n] = TestProblemData.HII_Fraction * 
	    TestProblemData.HydrogenFractionByMass * TestStarParticleDensity;

    BaryonField[HeIINum][n] =  TestProblemData.HeII_Fraction *
	    TestStarParticleDensity * (1.0-TestProblemData.HydrogenFractionByMass);

    BaryonField[HeIIINum][n] = TestProblemData.HeIII_Fraction *
	    TestStarParticleDensity * (1.0-TestProblemData.HydrogenFractionByMass);

    BaryonField[HeINum][n] =
	    (1.0 - TestProblemData.HydrogenFractionByMass)*TestStarParticleDensity -
	    BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];
	}

  if (MultiSpecies > 1) {
	  BaryonField[HMNum][n] = TestProblemData.HM_Fraction *
	    BaryonField[HIINum][n];

	  BaryonField[H2INum][n] = TestProblemData.H2I_Fraction *
	    BaryonField[0][n] * TestProblemData.HydrogenFractionByMass;

	  BaryonField[H2IINum][n] = TestProblemData.H2II_Fraction * 2.0 *
	    BaryonField[HIINum][n];

	  BaryonField[kdissH2INum][n] = 0.0; 
	  BaryonField[kphHMNum][n] = 0.0;
	  BaryonField[kdissH2IINum][n] = 0.0;
	}

  // HI density is calculated by subtracting off the various ionized fractions
  // from the total
  BaryonField[HINum][n] = TestProblemData.HydrogenFractionByMass*BaryonField[0][n]
	  - BaryonField[HIINum][n];

	if (MultiSpecies > 1)
	  BaryonField[HINum][n] -= (BaryonField[HMNum][n] + BaryonField[H2IINum][n]
				  + BaryonField[H2INum][n]);

  BaryonField[DeNum][n] = BaryonField[HIINum][n] +
	  0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];

	if (MultiSpecies > 1)
	  BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] - 
	    BaryonField[HMNum][n];

	/* Set Deuterium species (assumed to be negligible). */

	if (MultiSpecies > 2) {
	  BaryonField[DINum ][n]  = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HINum][n];
	  BaryonField[DIINum][n] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HIINum][n];
	  BaryonField[HDINum][n] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * BaryonField[H2INum][n];
	}

	/* Set energy (thermal and then total if necessary). */

	if (MultiSpecies) {
	  mu_data =  
	    0.25*(BaryonField[HeINum][n]  + BaryonField[HeIINum][n] +
		  BaryonField[HeIIINum][n]                        ) +
	    BaryonField[HINum][n]   + BaryonField[HIINum][n]  +
	    BaryonField[DeNum][n];
	  if (MultiSpecies > 1)
	    mu_data += BaryonField[HMNum][n]   +
	      0.5*(BaryonField[H2INum][n]  + BaryonField[H2IINum][n]);
	  mu_data = BaryonField[0][n] / mu_data;
	} else
	  mu_data = mu;

  BaryonField[1][n] = temperature/TemperatureUnits/
	  ((Gamma-1.0)*mu_data);

      }

  delete [] density_field;
  delete [] HII_field;
  delete [] HeII_field;
  delete [] HeIII_field;
  delete [] Temperature_field;

  return SUCCESS;
}


/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COOLING TEST)
/
/  written by: John Wise
/  date:       April, 2009
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
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/********************* PROTOTYPES *********************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);

/*******************************************************/

int grid::CoolingTestInitializeGrid(float MinimumDensity,
				float MaximumDensity,
				float MinimumTemperature,
				float MaximumTemperature,
				float MinimumColour,
				float MaximumColour,
				int UseMetals,
				int UseElectronFraction)
{

  /* declarations */

  const float Z_SOLAR = 0.0204;
  int dim, i, j, k, field, size, index;
  int DensNum, TENum, GENum, Vel1Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum, ColourNum;

  /* create fields */
  NumberOfBaryonFields = 0;
  FieldType[DensNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[GENum = NumberOfBaryonFields++] = InternalEnergy;
  FieldType[Vel1Num = NumberOfBaryonFields++] = Velocity1;
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
  if (UseMetals)
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Set various units. */

  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.67e-8,
               pi = 3.14159, mh = 1.67e-24, kboltz = 1.381e-16;
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    MassUnits, VelocityUnits, mu = 0.6;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, &MassUnits, Time);

  /* Set up the baryon field. */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* allocate fields */

  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];

  /* Loop over the mesh. */

  float density, temperature, colour;
  float delta_density, delta_temp, delta_colour;
  float log_density0, log_temp0, log_colour0;
  
  log_density0 = log10(MinimumDensity);
  log_temp0 = log10(MinimumTemperature);
  log_colour0 = log10(MinimumColour);
  delta_density = log10(MaximumDensity / MinimumDensity) / 
    (float) (GridDimension[0]-1);
  delta_temp = log10(MaximumTemperature / MinimumTemperature) / 
    (float) (GridDimension[1]-1);
  delta_colour = log10(MaximumColour / MinimumColour) / 
    (float) (GridDimension[2]-1);

  index = 0;
  for (k = 0; k < GridDimension[2]; k++) {
    colour = pow(10.0, log_colour0 + k * delta_colour);

    if (UseElectronFraction) {
      TestProblemData.HI_Fraction = CoolData.HydrogenFractionByMass * (1.0 - colour);
      TestProblemData.HII_Fraction = CoolData.HydrogenFractionByMass * colour;
      TestProblemData.HeI_Fraction = 4*(1.0 - CoolData.HydrogenFractionByMass) * 0.5*(1.0 - colour);
      TestProblemData.HeII_Fraction = 4*(1.0 - CoolData.HydrogenFractionByMass) * 0.5*(1.0 - colour);
      TestProblemData.HeIII_Fraction = 4*(1.0 - CoolData.HydrogenFractionByMass) * colour;
    }

    for (j = 0; j < GridDimension[1]; j++) {
      temperature = pow(10.0, log_temp0 + j * delta_temp);
      for (i = 0; i < GridDimension[0]; i++, index++) {

	density = pow(10.0, log_density0 + i * delta_density);

	/* Set density. */

	BaryonField[DensNum][index] = density;

	/* If doing multi-species (HI, etc.), set these. */

	if (MultiSpecies > 0) {
	  BaryonField[HINum][index] = density * TestProblemData.HI_Fraction;
	  BaryonField[HIINum][index] = density * TestProblemData.HII_Fraction;
	  BaryonField[HeINum][index] = density * TestProblemData.HeI_Fraction;
	  BaryonField[HeIINum][index] = density * TestProblemData.HeII_Fraction;
	  BaryonField[HeIIINum][index] = density * TestProblemData.HeIII_Fraction;
	}
	if (MultiSpecies > 1) {
	  BaryonField[HMNum][index] = TestProblemData.HM_Fraction * 
	    BaryonField[HIINum][index] * pow(temperature, float(0.88));
	  BaryonField[H2IINum][index] = TestProblemData.H2II_Fraction *
	    2.0*BaryonField[HIINum][index] * pow(temperature,float(1.8));
	  BaryonField[H2INum][index] = density * TestProblemData.H2I_Fraction;
	  BaryonField[HINum][index] -= BaryonField[HMNum][index]
	    + BaryonField[H2IINum][index]
	    + BaryonField[H2INum][index];
	}
	
	BaryonField[DeNum][index] = BaryonField[HIINum][index] + 
	  0.25*BaryonField[HeIINum][index] + 0.5*BaryonField[HeIIINum][index];
	if (MultiSpecies > 1)
	  BaryonField[DeNum][index] += 0.5*BaryonField[H2IINum][index] - 
	    BaryonField[HMNum][index];
	
	/* Set Deuterium species (assumed to be negligible). */
	
	if (MultiSpecies > 2) {
	  BaryonField[DINum][index] = CoolData.DeuteriumToHydrogenRatio*
	    BaryonField[HINum][index];
	  BaryonField[DIINum][index] = CoolData.DeuteriumToHydrogenRatio*
	    BaryonField[HIINum][index];
	  BaryonField[HDINum][index] = CoolData.DeuteriumToHydrogenRatio*
	    BaryonField[H2INum][index];
	}
	
	
	/* If there is a colour field, set it. */
	
	if (UseMetals)
	  BaryonField[MetalNum][index] = density * Z_SOLAR * colour;
	
	/* Set Velocities. */
	
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[Vel1Num+dim][index] = 0.0;
	
	/* Set energy (thermal and then total if necessary). */
	
	BaryonField[TENum][index] = temperature/TemperatureUnits/
	  ((Gamma-1.0)*mu);
	
	if (DualEnergyFormalism)
	  BaryonField[GENum][index] = BaryonField[TENum][index];
	
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[TENum][index] += 0.5*pow(BaryonField[Vel1Num+dim][index], 2);

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k
	
  return SUCCESS;

}


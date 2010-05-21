/***********************************************************************
/
/  GRID CLASS (OUTPUT STAR PARTICLE INFORMATION)
/
/  written by: Greg Bryan
/  date:       September, 1999
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
#include "CosmologyParameters.h"
#include "Grid.h"
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int grid::OutputStarParticleInformation(FILE *StarFile)
{
 
 
  /* Declarations */
 
  int i, j, dim;
 
  /* Compute the volume factor. */
 
  float CellVolume = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    CellVolume *= CellWidth[dim][0];
 
  /* Set the Conversion factors for Density and X-rays.  If using comoving
     coordinates use solar masses and Mpc as the intrinsic units.
     Note: The X-ray units have been multiplied by 1.0e-20 to stop overflow.
     Note: The temperature field already has units of K. */
 
  float MassConversion = CellVolume;
  float TemperatureUnits=1, DensityUnits=1, LengthUnits=1, VelocityUnits=1, 
    TimeUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
    const double SolarMass = 1.989e33 /* , Mpc = 3.0824e24 */ ;
    MassConversion *= float(double(DensityUnits)*POW(double(LengthUnits), 3)/
			   SolarMass);
  }
 
  /* If using star particles, write out star particle info. */
 
  float32 buffer[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES+4];
  for (i = 0; i < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES+4; i++)
    buffer[i] = 0;
  if (debug) printf("NumberOfParticles = %"ISYM"\n", NumberOfParticles);
  if (StarParticleCreation > 0 && StarFile != NULL &&
      NumberOfParticleAttributes > 0) {
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleAttribute[0][i] > 0) {
	for (j = 0; j < GridRank; j++)
	  buffer[j] = float32(ParticlePosition[j][i]);
	buffer[3] = float32(ParticleMass[i]*MassConversion);
	for (j = 0; j < NumberOfParticleAttributes; j++)
	  buffer[j+4] = float32(ParticleAttribute[j][i]);
	buffer[4] *= TimeUnits;
	buffer[5] *= TimeUnits;
	if (fwrite( (void*) buffer,  sizeof(float32),

		    NumberOfParticleAttributes+4, StarFile) !=
	    NumberOfParticleAttributes+4)
	  perror("error in fwrite");
      }
  }
 
  return SUCCESS;
}

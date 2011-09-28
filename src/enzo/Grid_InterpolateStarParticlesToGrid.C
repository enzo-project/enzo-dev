/***********************************************************************
/
/  GRID CLASS (GRID STAR PARTICLES ONTO THE MESH)
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/  modified1:  

/  PURPOSE: Make gridded outputs for (traditional) star particles 
/           (type = PARTICLE_TYPE_STAR)
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 

  
int grid::InterpolateStarParticlesToGrid(int NumberOfSPFields)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0 || OutputGriddedStarParticle <= 0)
    return SUCCESS;
 
  /* initialize */
 
  int dim, i, j, k, size, field;
  int xindex, yindex, zindex, index;
  int StarDensNum, FormingStarDensNum, SFRDensNum, CreationTimeNum;
  int ActiveDim[MAX_DIMENSION];
  int NumberOfStarParticlesInGrid = 0;
  double xv1 = 0.0, xv2 = 0.0, minitial = 0.0, mform = 0.0;
  float CellWidthTemp = float(CellWidth[0][0]);
  const double Msun = 1.989e33, yr = 3.15557e7, kpc = 3.086e21;
  FLOAT dtForSFR; 

  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
 
  dtForSFR = StarMakerMinimumDynamicalTime * yr / TimeUnits;  // = 1.35e-5 in code unit

  /* Compute size (in floats) of the current grid. */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    size *= ActiveDim[dim];
  }

  /* Assign field numbers (0-4 is left out because it is used in Grid_InterpolateParticlesToGrid) */   

#define NumberOfInterpolatedFieldsForDM 10

  StarDensNum        = NumberOfInterpolatedFieldsForDM;
  FormingStarDensNum = NumberOfInterpolatedFieldsForDM + 1;
  SFRDensNum         = NumberOfInterpolatedFieldsForDM + 2;
  CreationTimeNum    = NumberOfInterpolatedFieldsForDM + 3;  

  /* Set up empty fields */

  if (MyProcessorNumber == ProcessorNumber)
    for (field = NumberOfInterpolatedFieldsForDM; 
	 field < NumberOfInterpolatedFieldsForDM+NumberOfSPFields; field++) {
      InterpolatedField[field] = new float[size];  
      if (field == CreationTimeNum) 
	for (i = 0; i < size; i++)
	  InterpolatedField[field][i] = 0.0;
      else
	for (i = 0; i < size; i++)
	  InterpolatedField[field][i] = 1.0e-10;
    }

  if (NumberOfParticles == 0)
    return SUCCESS;


  /* ------------------------------------------------------------------- */
  /*                       NOW GRID STAR PARTICLES                       */
  /* ------------------------------------------------------------------- */


  if (NumberOfParticles > 0) {

    for (i = 0; i < NumberOfParticles; i++) {

#define PARTICLE_IN_GRID_CHECK

#ifdef PARTICLE_IN_GRID_CHECK      

      xindex = (int)((ParticlePosition[0][i] - CellLeftEdge[0][0]) / CellWidthTemp);
      yindex = (int)((ParticlePosition[1][i] - CellLeftEdge[1][0]) / CellWidthTemp); 
      zindex = (int)((ParticlePosition[2][i] - CellLeftEdge[2][0]) / CellWidthTemp); 

      if (xindex < 0 || xindex > GridDimension[0] || 
	  yindex < 0 || yindex > GridDimension[1] || 
	  zindex < 0 || zindex > GridDimension[2])
	fprintf(stdout, "particle out of grid (C level); xind, yind, zind = %d, %d, %d\n",
		xindex, yindex, zindex); 
#endif

      /* Now store particle info */

      if (ParticleType[i] == PARTICLE_TYPE_STAR) {

	NumberOfStarParticlesInGrid++;

	// now different [xyz]index
	xindex = (int)((ParticlePosition[0][i] - GridLeftEdge[0]) / CellWidthTemp);
	yindex = (int)((ParticlePosition[1][i] - GridLeftEdge[1]) / CellWidthTemp); 
	zindex = (int)((ParticlePosition[2][i] - GridLeftEdge[2]) / CellWidthTemp); 
	
	index = xindex + ActiveDim[0] * (yindex + ActiveDim[1]*zindex);
	
	// (10) star particle density (in code density unit)
	InterpolatedField[StarDensNum][index] += ParticleMass[i];  
	
	if (OutputGriddedStarParticle > 1) { 

	  xv1 = (Time            - ParticleAttribute[0][i])/ParticleAttribute[1][i];
	  xv2 = (Time + dtForSFR - ParticleAttribute[0][i])/ParticleAttribute[1][i];
	  
	  /* For actively star-forming particles ( active if (t-t_cr)/t_dyn <= 12 ).
	     This assumes a particular star formation criteria in star_maker7.src */

	  if (xv1 <= 12.0) {
	    
	    minitial = ParticleMass[i] / (1.0 - StarMassEjectionFraction*(1.0 - (1.0 + xv1)*exp(-xv1)));	    
	    mform = minitial * ((1.0 + xv1)*exp(-xv1) - (1.0 + xv2)*exp(-xv2));
	    mform = max(min(mform, ParticleMass[i]), 0.0);
	    
	    // (11) forming stellar mass density (in code density unit)
	    InterpolatedField[FormingStarDensNum][index] += (float)((1.0 - StarMassEjectionFraction)*mform); 
	    
	    // (12) SFR density in current timestep 'dtFixed'
	    InterpolatedField[SFRDensNum][index] += (float)((1.0 - StarMassEjectionFraction)*mform);
	  }

	  // (13) average creation time 
	  InterpolatedField[CreationTimeNum][index] += ParticleAttribute[0][i]; 

	} 

      } //if PARTICLE_TYPE_STAR

      
    } //ENDFOR i in NumberOfParticles

    /* Divide certain fields to get the correct value */ 

    if (OutputGriddedStarParticle > 1 && NumberOfStarParticlesInGrid > 0) 
      for (i = 0; i < size; i++) {
	
	// (12) SFR density (code density unit  ->  Ms/yr/kpc^3)
	InterpolatedField[SFRDensNum][i] *= 
	  DensityUnits / Msun * pow(kpc, 3) / ( dtForSFR * TimeUnits / yr );
	
	// (13) (in code time unit)
	InterpolatedField[CreationTimeNum][i] /= NumberOfStarParticlesInGrid;    

      }

  } //NumberOfParticles > 0


//    printf("Time = %g, dtForSFR = %g, xv1 = %g, xv2 = %g, minitial = %g, mform = %g, ParticleAttribute[0][i-1] = %g\n", 
//	   Time, dtForSFR, xv1, xv2, minitial, mform, ParticleAttribute[0][i-1]);  

  return SUCCESS;

}

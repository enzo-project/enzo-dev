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
 
extern "C" void FORTRAN_NAME(particle_splitter)(int *nx, int *ny, int *nz,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
	     FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, int *ibuff,
             int *npart,
             FLOAT *xpold, FLOAT *ypold, FLOAT *zpold, float *upold, float *vpold, float *wpold,
	     float *mpold, float *tdpold, float *tcpold, float *metalfold, int *typeold,
	     int *nmax, int *npartnew, int *children, int *level,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
	     float *mp, float *tdp, float *tcp, float *metalf, int *type, 
             int *iterations, float *separation, int *ran1_init); 

  
int grid::InterpolateStarParticlesToGrid(void)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* initialize */
 
  int dim, i, j, k, size, field, NumberOfFields;
  int xindex, yindex, zindex, index;
  int StarDensNum, FormingStarDensNum, SFRNum, CreationTimeNum;
  int ActiveDim[MAX_DIMENSION];
  float CellWidthTemp = float(CellWidth[0][0]);
  int NumberOfStarParticlesInGrid = 0;
  float xv1, xv2, minitial, mform;

  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
 
  /* Compute size (in floats) of the current grid. */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    size *= ActiveDim[dim];
  }

  /* Assign number of fields */

  switch (OutputGriddedStarParticle) {
  case 1: NumberOfFields = 1; break;  // star particle density
  case 2: NumberOfFields = 4; break;  // + forming stellar mass density + SFR + average creation time
  default: 
    fprintf(stdout, "Unrecognized value for OutputGriddedStarParticle = %"ISYM"\n",
	    OutputGriddedStarParticle);
    fprintf(stdout, "Setting to 1.  Outputting particle density only.\n");
    OutputGriddedStarParticle = 1;
    NumberOfFields = 1;
    break;
  } 

  /* Assign field numbers (0-4 is left out because it is used in Grid_InterpolateParticlesToGrid) */   

#define NumberOfInterpolatedFieldsForDM 10

  StarDensNum        = NumberOfInterpolatedFieldsForDM;
  FormingStarDensNum = NumberOfInterpolatedFieldsForDM + 1;
  SFRNum             = NumberOfInterpolatedFieldsForDM + 2;
  CreationTimeNum    = NumberOfInterpolatedFieldsForDM + 3;  

  /* Set up empty fields */

  if (MyProcessorNumber == ProcessorNumber)
    for (field = NumberOfInterpolatedFieldsForDM; 
	 field < NumberOfInterpolatedFieldsForDM+NumberOfFields; field++) {
      InterpolatedField[field] = new float[size];  
      for (i = 0; i < size; i++)
	InterpolatedField[field][i] = 1.0e-10;
    }

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
	
	xv1 = (Time           - ParticleAttribute[0][i])/ParticleAttribute[1][i];
	xv2 = (Time + dtFixed - ParticleAttribute[0][i])/ParticleAttribute[1][i];

	// For actively star-forming particles ( active if (t-t_cr)/t_dyn <= 12 )
	if(xv1 <= 12) {

	  minitial = ParticleMass[i] /(1.0 - StarMassEjectionFraction*(1.0 - (1.0 + xv1)*exp(-xv1)));

	  mform = minitial * ((1.0 + xv1)*exp(-xv1) - (1.0 + xv2)*exp(-xv2));
	  mform = max(min(mform, ParticleMass[i]), 0.0);

	  // (11) forming stellar mass density (in code density unit)
	  InterpolatedField[FormingStarDensNum][index] += mform; 

	}

	// (13) average creation time
	InterpolatedField[CreationTimeNum][index] += ParticleAttribute[0][i]; 

      } //if PARTICLE_TYPE_STAR

      
    } //ENDFOR i in NumberOfParticles

    // (12) SFR in current timestep 'dtFixed' (in Ms/yr)
    InterpolatedField[SFRNum][index] *= 
      DensityUnits*LengthUnits*LengthUnits*LengthUnits/dtFixed/TimeUnits; 

    // (13) average creation time (in code time unit)
    InterpolatedField[CreationTimeNum][index] /= 
      NumberOfStarParticlesInGrid;    

  } //NumberOfParicles > 0


  fprintf(stdout, "InterpolatedField[10][30] = %g\n", InterpolatedField[10][30]);  //#####
  fprintf(stdout, "InterpolatedField[11][30] = %g\n", InterpolatedField[11][30]);  //#####
  fprintf(stdout, "InterpolatedField[12][30] = %g\n", InterpolatedField[12][30]);  //#####
  fprintf(stdout, "InterpolatedField[13][30] = %g\n", InterpolatedField[13][30]);  //#####


  return SUCCESS;

}

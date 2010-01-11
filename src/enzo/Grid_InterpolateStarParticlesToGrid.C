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
  int StarDensNum, FormingStarDensNum, SFRDensNum, CreationTimeNum;
  int ActiveDim[MAX_DIMENSION];

  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
 
  float CellWidthTemp = float(CellWidth[0][0]);

  /* Compute size (in floats) of the current grid. */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    size *= ActiveDim[dim];
  }

  /* Assign number of fields */

  switch (OutputGriddedStarParticle) {
  case 1: NumberOfFields = 1; break;  // particle density
  case 2: NumberOfFields = 4; break;  // + forming star particle density + SFR density, etc.
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
  SFRDensNum         = NumberOfInterpolatedFieldsForDM + 2;
  CreationTimeNum    = NumberOfInterpolatedFieldsForDM + 3;  

  /* Set up empty fields */

  if (MyProcessorNumber == ProcessorNumber)
    for (field = NumberOfInterpolatedFieldsForDM; 
	 field < NumberOfInterpolatedFieldsForDM+NumberOfFields; field++) {
      InterpolatedField[field] = new float[size];  
      for (i = 0; i < size; i++)
	InterpolatedField[field][i] = 1.0e-6;
    }

  /* ------------------------------------------------------------------- */
  /*                       NOW GRID STAR PARTICLES                       */
  /* ------------------------------------------------------------------- */

  fprintf(stdout,"grid::InterpolateSPsToGrid: PROC = %d, ID = %d, NumberOfParticles = %d\n", 
	  MyProcessorNumber, ID, NumberOfParticles);

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
	fprintf(stdout, "particle out of grid (C level); xind, yind, zind = %d, %d, %d, %d\n",
		xindex, yindex, zindex); 
#endif

      // now different [xyz]index
      if (ParticleType[i] == PARTICLE_TYPE_STAR) {

	xindex = (int)((ParticlePosition[0][i] - GridLeftEdge[0]) / CellWidthTemp);
	yindex = (int)((ParticlePosition[1][i] - GridLeftEdge[1]) / CellWidthTemp); 
	zindex = (int)((ParticlePosition[2][i] - GridLeftEdge[2]) / CellWidthTemp); 
	
	index = xindex + ActiveDim[0] * (yindex + ActiveDim[1]*zindex);
	
	InterpolatedField[StarDensNum][index] += ParticleMass[i];  // in code density unit
	//	InterpolatedField[FormingStarDensNum][index] += ParticleMass[i]; 
	//	InterpolatedField[SFRDensNum][index] += ParticleMass[i]; 

      }

    //    InterpolatedField[SFRDensNum][index] /= dt;
      
    } //ENDFOR i

  } //NumberOfParicles > 0

  fprintf(stdout,"grid::InterpolateSPsToGrid: PROC = %d, InterpolatedField[StarDensNum][30] = %g\n", 
	  MyProcessorNumber, InterpolatedField[StarDensNum][30]);


  return SUCCESS;

}

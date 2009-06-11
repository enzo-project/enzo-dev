/***********************************************************************
/
/  GRID CLASS (OUTPUT GRID/DM/STAR DATA STORED ON THIS GRID)
/
/  written by: Greg Bryan
/  date:       December, 1999
/  modified1:
/
/  PURPOSE:  Output selected information from this grid into one of
/     the output files.  Done only for data within the selected region
/     (and for grid points which are not covered by a finer mesh).
/     The value of WriteTime is currently ignored for the grid, but
/     is used to move particle positions to an approximation of where
/     they should be a time = WriteTime.
/
/  RETURNS:
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
#include "Grid.h"
#include "StarParticleData.h"
 
#define OUTPUT_GRID_VELOCITIES
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
/* This macro converts a float and writes it to the local buffer, which,
   when full is written to the file pointed to by fptr. */
 
 
int grid::OutputGridMovieData(FILE *Gridfptr, FILE *DMfptr, FILE *Starfptr,
			      FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[],
			      FLOAT WriteOutTime, int NumberOfPoints[3],
			      int NumberOfValuesPerPoint[3],
			      char *PointValueNames[3][20], float BaseRadius)
{
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  float CellVolume = 1.0, density;
  int i, j, k, dim, index, n, size = 1, FieldPosition[] = {0,0,0};
 
  for (dim = 0; dim < GridRank; dim++) {
    Left[dim] = max(RegionLeftEdge[dim], GridLeftEdge[dim]);
    Right[dim] = min(RegionRightEdge[dim], GridRightEdge[dim]);
    if (Left[dim] >= Right[dim])
      return SUCCESS;
    size *= GridDimension[dim];
    CellVolume *= CellWidth[dim][0];
  }
 
  /* If WriteOutTime is < 0, then set it to the current time. */
 
  FLOAT WriteTime = WriteOutTime;
  if (WriteOutTime < 0)
    WriteTime = Time;
 
  /* ------------------------------------------------------------- */
  /* 1) Output grid information. */
 
  if (NumberOfBaryonFields > 0 && Gridfptr != NULL) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, ibuf;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      ENZO_FAIL("");
    }
 
    /* Compute coefficient factors for linear interpolation in time.
       Note: interp = coef1*Old + coef2*New. */
 
    float coef1 = 0, coef2 = 1;
    if (Time != WriteTime && OldBaryonField[DensNum] != NULL) {
      if (Time <= OldTime) {
	fprintf(stderr, "OGMD: fields are at the same time or worse.\n");
	//	ENZO_FAIL("");
      } else {
	coef1 = max((Time - WriteTime)/
		    (Time - OldTime), 0.0);
	coef2  = (1.0 - coef1);
      }
    }
 
    /* Compute temperature. */
 
    float *temperature = new float[size];
    if (this->ComputeTemperatureField(temperature) == FAIL) {
      fprintf(stderr, "Error in grid->ComputeTemperatureField.\n");
      ENZO_FAIL("");
    }
 
    /* Find metallicity field and set flag. */
 
    int MetallicityField = FALSE, MetalNum;
    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
	!= -1)
      MetallicityField = TRUE;
    else
      MetalNum = 0;
 
    /* Compute start and end index of region to be output. */
 
    int StartIndex[] = {0,0,0}, EndIndex[] = {0,0,0};
    for (dim = 0; dim < GridRank; dim++) {
      StartIndex[dim] = nint((Left[dim] - CellLeftEdge[dim][0])/
			     CellWidth[dim][0]);
      EndIndex[dim] = nint((Right[dim] - CellLeftEdge[dim][0])/
			   CellWidth[dim][0]) - 1;
    }
 
    /* Set the number of values per point (and give names). */
 
    NumberOfValuesPerPoint[0] = GridRank + 3
#ifdef OUTPUT_GRID_VELOCITIES
                              + GridRank /* for velocities. */
#endif
                              + MetallicityField;
  if (PointValueNames[0][0] == NULL) {
      i = 0;
      PointValueNames[0][i++] = "x_position (in units of box size)";
      if (GridRank > 1)
	PointValueNames[0][i++] = "y_position (in units of box size)";
      if (GridRank > 1)
	PointValueNames[0][i++] = "z_position (in units of box size)";
      PointValueNames[0][i++] = "mass (in units of box mass)";
      PointValueNames[0][i++] = "cell size (in units of box size)";
      PointValueNames[0][i++] = "temperature (K)";
      if (MetallicityField)
	PointValueNames[0][i++] = "metallicity (fraction by mass)";
#ifdef OUTPUT_GRID_VELOCITIES
      for (dim = 0; dim < GridRank; dim++)
	PointValueNames[0][i++] = "velocity (code units)";
#endif
    }
 
    /* Allocate a 32 bit write buffer. */
 
    int BufferSize = (EndIndex[0]-StartIndex[0]+1)*NumberOfValuesPerPoint[0];
    float32 *buffer = new float32[BufferSize];
 
    /* Loop over grid and write out data through write buffer
       (ignore points which are covered by finger grids). */
 
    for (k = StartIndex[2]; k <= EndIndex[2]; k++)
      for (j = StartIndex[1]; j <= EndIndex[1]; j++) {
	index = (j+k*GridDimension[1])*GridDimension[0] + StartIndex[0];
	ibuf = 0;
	for (i = StartIndex[0]; i <= EndIndex[0]; i++, index++)
	  if (BaryonField[NumberOfBaryonFields][index] == 0) {
	    buffer[ibuf++] = float32(CellLeftEdge[0][i]+0.5*CellWidth[0][i]);
	    if (GridRank > 1)
	      buffer[ibuf++] = float32(CellLeftEdge[1][j]+0.5*CellWidth[1][j]);
	    if (GridRank > 2)
	      buffer[ibuf++] = float32(CellLeftEdge[2][k]+0.5*CellWidth[2][k]);
	    if (coef2 == 1)
	      buffer[ibuf++] = float32(BaryonField[DensNum][index]*CellVolume);
	    else
	      buffer[ibuf++] = float32((coef1*OldBaryonField[DensNum][index] +
					coef2*BaryonField[DensNum][index])
				       *CellVolume);
	    //	    buffer[ibuf++] = float32(CellWidth[0][0]);
	    buffer[ibuf++] = float32(BaseRadius *
			  POW(BaryonField[DensNum][index], float(-1.0/3.0)));
	    buffer[ibuf++] = float32(temperature[index]);
	    if (MetallicityField)
	      buffer[ibuf++] = float32(BaryonField[MetalNum][index]/
				       BaryonField[DensNum][index]);
#ifdef OUTPUT_GRID_VELOCITIES
	    for (dim = 0; dim < GridRank; dim++)
	      buffer[ibuf++] = float32(BaryonField[Vel1Num+dim][index]);
#endif
	    NumberOfPoints[0]++;  /* increment # of grid pts written */
	  }
	fwrite((void*) buffer, sizeof(float32), ibuf, Gridfptr);
      }
 
    delete [] temperature;
    delete [] buffer;
 
  } // end: if (NumberOfBaryonFields > 0)
 
  /* ------------------------------------------------------------- */
  /* 2) Output particle information. */
 
  FLOAT a, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(WriteTime, &a, &dadt)
	== FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      ENZO_FAIL("");
    }
 
  FLOAT TempPos[MAX_DIMENSION];
  float Coefficient = (WriteTime-Time)/a;
 
  /* Set number of values per point (and their names). */
 
  int NumberOfDMValues = GridRank+3, NumberOfStarValues = 0;
  NumberOfValuesPerPoint[1] = NumberOfDMValues;
  if (StarParticleCreation > 0)
    NumberOfStarValues = GridRank + 3 +
      ((NumberOfParticleAttributes > 2) ? 2 : 0);
  NumberOfValuesPerPoint[2] = NumberOfStarValues;
  for (j = 1; j < 3; j++)
    if (PointValueNames[j][0] == NULL) {
      i = 0;
      PointValueNames[j][i++] = "x_position (in units of box size)";
      if (GridRank > 1)
	PointValueNames[j][i++] = "y_position (in units of box size)";
      if (GridRank > 1)
	PointValueNames[j][i++] = "z_position (in units of box size)";
      PointValueNames[j][i++] = "mass (in units of box mass)";
      PointValueNames[j][i++] = "particle size (in units of box size)";
      PointValueNames[j][i++] = "particle index";
      if (NumberOfParticleAttributes > 2 && j == 2) {
	PointValueNames[j][i++] = "creation time (code units)";
	PointValueNames[j][i++] = "metallicity (fraction by mass)";
      }
    }
 
  /* Allocate a 32 bit write buffer. */
 
  int BufferSize = 128, idmbuf = 0, istarbuf = 0;
  float32 *dmbuffer = new float32[BufferSize*NumberOfDMValues], *float32_point;
  float32 *starbuffer = new float32[BufferSize*NumberOfStarValues];
 
  /* Loop over particles, out dm and star particles */
 
  for (n = 0; n < NumberOfParticles; n++) {
 
    /* Check to see if particle is in region. */
 
    for (dim = 0; dim < GridRank; dim++) {
      TempPos[dim] = ParticlePosition[dim][n] +
	             Coefficient*ParticleVelocity[dim][n];
      if (TempPos[dim] < RegionLeftEdge[dim] ||
	  TempPos[dim] > RegionRightEdge[dim])
	break;
    }
 
    if (dim == GridRank) {
 
      /* Get particle density for radius computation. */
 
      if (GravitatingMassFieldParticles != NULL) {
	for (dim = 0; dim < GridRank; dim++)
	  FieldPosition[dim] =
	    int((ParticlePosition[dim][n] -
		 GravitatingMassFieldParticlesLeftEdge[dim])
		/GravitatingMassFieldParticlesCellSize);
	index = (FieldPosition[1] +
		 FieldPosition[2]*GravitatingMassFieldParticlesDimension[1])*
	  GravitatingMassFieldParticlesDimension[0] + FieldPosition[0];
	density = GravitatingMassFieldParticles[index];
      } else
	density = 0;
      if (density == 0) density = ParticleMass[n];
 
      /* What we output depends on particle type */
 
      if (ParticleType[n] == PARTICLE_TYPE_DARK_MATTER ||
	  (ParticleType[n] == PARTICLE_TYPE_MUST_REFINE &&
	   ParticleMass[n] > 0)) {
 
	/* Dark Matter. */
 
	for (dim = 0; dim < GridRank; dim++)
	  dmbuffer[idmbuf++] = float32(TempPos[dim]);
	dmbuffer[idmbuf++] = float32(ParticleMass[n]*CellVolume);
	//	dmbuffer[idmbuf++] = float32(CellWidth[0][0]);
	dmbuffer[idmbuf++] = float32(BaseRadius*POW(density, float(-1.0/3.0)));
	float32_point = (float32 *) (ParticleNumber+n);
	dmbuffer[idmbuf++] = *float32_point;  /* treat int as float -- bad! */
	NumberOfPoints[1]++;  /* increment # of grid pts written */
	if (idmbuf >= BufferSize*NumberOfDMValues) {
	  if (DMfptr != NULL)
	    fwrite((void*) dmbuffer, sizeof(float32), idmbuf, DMfptr);
	  idmbuf = 0;
	}
 
      } // end: if (ParticleType == dm)
 
      /* Star particle. */
 
      if (ParticleType[n] == PARTICLE_TYPE_STAR) {
 
	for (dim = 0; dim < GridRank; dim++)
	  starbuffer[istarbuf++] = float32(TempPos[dim]);
	starbuffer[istarbuf++] = float32(ParticleMass[n]*CellVolume);
	//	starbuffer[istarbuf++] = float32(CellWidth[0][0]);
	starbuffer[istarbuf++] = float32(BaseRadius*
					 POW(density, float(-1.0/3.0)));
	float32_point = (float32 *) (ParticleNumber+n);
	starbuffer[istarbuf++] = *float32_point;  /* again - bad! */
	if (NumberOfParticleAttributes > 2) {
	  starbuffer[istarbuf++] = float32(ParticleAttribute[0][n]);
	  starbuffer[istarbuf++] = float32(ParticleAttribute[2][n]);
	}
	NumberOfPoints[2]++;  /* increment # of grid pts written */
	if (istarbuf >= BufferSize*NumberOfStarValues) {
	  if (Starfptr != NULL)
	    fwrite((void*) starbuffer, sizeof(float32), istarbuf, Starfptr);
	  istarbuf = 0;
	}
 
      } // end: if (ParticleType == star)
 
    } // end: if (dim == GridRank)
 
  } // end loop over particles
 
  /* Flush all buffers. */
 
  if (DMfptr != NULL && idmbuf > 0)
    fwrite((void*) dmbuffer, sizeof(float32), idmbuf, DMfptr);
  if (Starfptr != NULL && istarbuf > 0)
    fwrite((void*) starbuffer, sizeof(float32), istarbuf, Starfptr);
 
  /* Clean up. */
 
  delete [] starbuffer;
  delete [] dmbuffer;
  delete [] BaryonField[NumberOfBaryonFields];
  BaryonField[NumberOfBaryonFields] = NULL;
 
  return SUCCESS;
}
 

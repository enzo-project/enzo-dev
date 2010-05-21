/***********************************************************************
/
/  GRID CLASS (SET THE TRACER PARTICLE VELOCITIES FROM GAS GRID)
/
/  written by: Greg Bryan
/  date:       March, 2004
/  modified1:
/
/  PURPOSE: This routine is used to set the velocities of the tracer
/     particles, which are used to follow stream-lines in the flow.
/     Tracer particles are massless and can be used to output values
/     of the gas as they advect with the fluid.  This routine sets
/     their velocity equal to the gas velocity in the current cell.
/
/  NOTE:
/
************************************************************************/
 
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
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int InterpolateTracerValues(FLOAT *Position[MAX_DIMENSION],
			float *InterpolatedValue[], int *ParticleType,
			int NumberOfParticles,
			FLOAT LeftEdge[MAX_DIMENSION], float CellSize,
			int Dimension[MAX_DIMENSION], int Rank,
			float *Fields[],
			int NumberOfFieldsToInterpolate, int FillMode);
 
int grid::TracerParticleSetVelocity()
{
  if (NumberOfParticles == 0 || NumberOfBaryonFields == 0 ||
      MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* Set the left edge of the field and the cellsize. */
 
  int dim;
  FLOAT LeftEdge[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++)
    LeftEdge[dim] = CellLeftEdge[dim][0];
  float CellSize = CellWidth[0][0]; // assumes all dimensions are the same
 
  /* Set the fields to interpolate (Baryon Velocity) */
 
  int Vel1Num;
  float *VelocityFields[MAX_DIMENSION];
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields)) < 0) {
    ENZO_FAIL("Could not find baryon velocity.\n");
  }
  for (dim = 0; dim < GridRank; dim++)
    VelocityFields[dim] = BaryonField[Vel1Num+dim];
 
  /* Interpolate the velocities from the grid for the tracer particles. */
 
  InterpolateTracerValues(ParticlePosition, ParticleVelocity, ParticleType,
			  NumberOfParticles, LeftEdge, CellSize,
			  GridDimension, GridRank,
			  VelocityFields, GridRank, 0);
 
  return SUCCESS;
}
 
 
/* This routine actually carries out the interpolation of tracer particles
   and is used both for the interpolation of the velocity fields and
   also for the density and temperature fields (as required in
   Grid_OutputTracerParticleData).  In the first case, the InterpolatedValue
   arrays are the ParticleVelocities and we overwrite the current velocity
   corresponding to that tracer particle with the Baryon velocity field.
   In the second case, we fill up an array with the particle positions,
   and the interpolated fields.
      Position[] - positions of particles
      InterpolatedValue[] - values interpolated from field are put here
      ParticleType - type of particle (only apply this to tracer particles)
      NumberOfParticles - number of all particles (tracer and other)
      LeftEdge - Array corresponding to left edge of the interpolated field
      CellSize - scalar of cell-size of interpolated field (must be cubic)
      Dimension - dimension of interpolated field
      Rank     - number of dimensions
      Field    - array of fields to be interpolated from
      NumberOfFieldsToInterpolate - size of Field array
      FillMode - This integer flag (0 or 1) indicates how InterpolateValue
                 array is filled (see above). */
 
int InterpolateTracerValues(FLOAT *Position[MAX_DIMENSION],
			float *InterpolatedValue[], int *ParticleType,
			int NumberOfParticles,
			FLOAT LeftEdge[MAX_DIMENSION], float CellSize,
			int Dimension[MAX_DIMENSION], int Rank,
			float *Field[],
			int NumberOfFieldsToInterpolate, int FillMode)
{
 
  /* Local declarations */
 
  int i1, j1, k1, field, n, fillindex = 0;
  float xpos, ypos, zpos, dx, dy, dz, fact, value;
  FLOAT edge1, edge2, edge3;
  const float half = 0.50001; // a bit more than 0.5 to ensure correct rounding
 
  /* Precompute factors */
 
  fact = 1.0/CellSize;
  edge1 = float(Dimension[0]) - half;
  if (Rank > 1)
    edge2 = float(Dimension[1]) - half;
  if (Rank > 2)
    edge3 = float(Dimension[2]) - half;
 
  /* ----------------------------------------------------------- */
 
  /* 1D Interpolation */
 
  if (Rank == 1) {
 
    for (n = 0; n < NumberOfParticles; n++) {
 
      /* Only do tracer particles */
 
      if (ParticleType[n] == PARTICLE_TYPE_TRACER) {
 
	/* Compute the position of the central cell */
 
	xpos = min(max((Position[0][n] - LeftEdge[0])*fact, half), edge1);
 
	/* Convert this into an integer index */
	
	i1  = int(xpos - 0.5);
 
	/* Compute the weights */
 
	dx = float(i1) + 1.5 - xpos;
 
	/* If in fill mode 1, then fill in tracer particle position. */
 
	if (FillMode == 1) {
	  InterpolatedValue[0][fillindex++] = float(Position[0][n]);
	}
 
	/* Interpolate from field into InterpolatedValue */
 
	for (field = 0; field < NumberOfFieldsToInterpolate; field++) {
	  value = Field[field][i1  ]*dx +
	          Field[field][i1+1]*(1.0-dx);
	  if (FillMode == 0)
	    InterpolatedValue[field][n] = value;
	  else
	    InterpolatedValue[0][fillindex++] = value;
	}
 
      } // end: if particle is tracer type
    } // end: loop over particles
  } // end: if Rank == 1
 
  /* ----------------------------------------------------------- */
 
  /* 2D Interpolation */
 
  if (Rank == 2) {
 
    for (n = 0; n < NumberOfParticles; n++) {
 
      if (ParticleType[n] == PARTICLE_TYPE_TRACER) {
 
	/* Compute the position of the central cell */
 
	xpos = min(max((Position[0][n] - LeftEdge[0])*fact, half), edge1);
	ypos = min(max((Position[1][n] - LeftEdge[1])*fact, half), edge2);
 
	/* Convert this into an integer index */
	
	i1  = int(xpos - 0.5);
	j1  = int(ypos - 0.5);
 
	/* Compute the weights */
 
	dx = float(i1) + 1.5 - xpos;
	dy = float(j1) + 1.5 - ypos;
 
	/* If in fill mode 1, then fill in tracer particle position. */
 
	if (FillMode == 1) {
	  InterpolatedValue[0][fillindex++] = float(Position[0][n]);
	  InterpolatedValue[0][fillindex++] = float(Position[1][n]);
	}
 
	/* Interpolate from field into InterpolatedValue */
 
	int index1 =  j1   *Dimension[0] + i1,
	    index2 = (j1+1)*Dimension[0] + i1;
 
	for (field = 0; field < NumberOfFieldsToInterpolate; field++) {
	  value =
	    Field[field][index1  ]*     dx *     dy  +
	    Field[field][index1+1]*(1.0-dx)*     dy  +
	    Field[field][index2  ]*     dx *(1.0-dy) +
	    Field[field][index2+1]*(1.0-dx)*(1.0-dy);
	  if (FillMode == 0)
	    InterpolatedValue[field][n] = value;
	  else
	    InterpolatedValue[0][fillindex++] = value;
	}
 
      } // end: if particle is tracer type
    } // end: loop over particles
  } // end: if Rank == 2
 
  /* ----------------------------------------------------------- */
 
  /* 3D Interpolation */
 
  if (Rank == 3) {
 
    for (n = 0; n < NumberOfParticles; n++) {
 
      if (ParticleType[n] == PARTICLE_TYPE_TRACER) {
 
	/* Compute the position of the central cell */
 
	xpos = min(max((Position[0][n] - LeftEdge[0])*fact, half), edge1);
	ypos = min(max((Position[1][n] - LeftEdge[1])*fact, half), edge2);
	zpos = min(max((Position[2][n] - LeftEdge[2])*fact, half), edge3);
 
	/* Convert this into an integer index */
	
	i1  = int(xpos - 0.5);
	j1  = int(ypos - 0.5);
	k1  = int(zpos - 0.5);
 
	/* Compute the weights */
 
	dx = float(i1) + 1.5 - xpos;
	dy = float(j1) + 1.5 - ypos;
	dz = float(k1) + 1.5 - zpos;
 
	/* If in fill mode 1, then fill in tracer particle position. */
 
	if (FillMode == 1) {
	  InterpolatedValue[0][fillindex++] = float(Position[0][n]);
	  InterpolatedValue[0][fillindex++] = float(Position[1][n]);
	  InterpolatedValue[0][fillindex++] = float(Position[2][n]);
	}
 
	/* Interpolate from field into InterpolatedValue */
 
	int index1 = ( k1   *Dimension[1] + j1  )*Dimension[0] + i1,
	    index2 = ( k1   *Dimension[1] + j1+1)*Dimension[0] + i1,
	    index3 = ((k1+1)*Dimension[1] + j1  )*Dimension[0] + i1,
    	    index4 = ((k1+1)*Dimension[1] + j1+1)*Dimension[0] + i1;
 
	for (field = 0; field < NumberOfFieldsToInterpolate; field++) {
	  value =
	    Field[field][index1  ]*     dx *     dy *     dz  +
	    Field[field][index1+1]*(1.0-dx)*     dy *     dz  +
	    Field[field][index2  ]*     dx *(1.0-dy)*     dz  +
	    Field[field][index2+1]*(1.0-dx)*(1.0-dy)*     dz  +
	    Field[field][index3  ]*     dx *     dy *(1.0-dz) +
	    Field[field][index3+1]*(1.0-dx)*     dy *(1.0-dz) +
	    Field[field][index4  ]*     dx *(1.0-dy)*(1.0-dz) +
	    Field[field][index4+1]*(1.0-dx)*(1.0-dy)*(1.0-dz);
	  if (FillMode == 0)

	    InterpolatedValue[field][n] = value;
	  else
	    InterpolatedValue[0][fillindex++] = value;
	}
 
      } // end: if particle is tracer type
    } // end: loop over particles
  } // end: if Rank == 3
 
  return SUCCESS;
}

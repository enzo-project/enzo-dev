/***********************************************************************
/
/  GRID CLASS (OUTPUT BARYON DATA FOR ANY TRACER PARTICLES ON THIS GRID)
/
/  written by: Greg Bryan
/  date:       March 2004
/  modified1:  BWO, 23 May 2006:  tracer particle velocities too! 
/                     (controlled by #ifdef TP_VELOCITY)
/  modified2:  BWO, Dec. 2013: updated tracer particles so they work for 
/                     1D and 2D simulations, also converted velocity output
/                     to a runtime parameter.
/
/  PURPOSE:  Output position, density, temperature of any tracer
/     particles in this grid.  Data is written in binary format.
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
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int InterpolateTracerValues(FLOAT *Position[MAX_DIMENSION],
			    float *InterpolatedValue[], int *ParticleType,
			    int NumberOfParticles,
			    FLOAT LeftEdge[MAX_DIMENSION], float CellSize,
			    int Dimension[MAX_DIMENSION], int Rank,
			    float *Fields[],
			    int NumberOfFieldsToInterpolate, int FillMode);
 
 
int grid::TracerParticleOutputData(FILE *fptr, FLOAT WriteOutTime)
{
 
  if (ProcessorNumber != MyProcessorNumber || NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* Count number of tracer particles on this grid. */
 
  int dim, n, NumberOfTracerParticles = 0, DensNum, size = 1;
  for (n = 0; n < NumberOfParticles; n++)
    if (ParticleType[n] == PARTICLE_TYPE_TRACER)
      NumberOfTracerParticles++;
 
  if (NumberOfTracerParticles == 0)
    return SUCCESS;
 
  /* Set the left edge of the field and the cellsize. */
  
  FLOAT LeftEdge[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    LeftEdge[dim] = CellLeftEdge[dim][0];
    size *= GridDimension[dim];
  }
  float CellSize = float(CellWidth[0][0]);  // assumes all dims are the same
 
  /* Compute temperature field. */
 
  float *FieldsToInterpolate[2], *temperature = new float[size];
  if (this->ComputeTemperatureField(temperature) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");
  }
 
  /* Set up an array of pointers (of size 2) pointing to the fields
     to be interpolated from: currently density and temperature. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    ENZO_FAIL("failed to find density field.\n");
  }

  FieldsToInterpolate[0] = BaryonField[DensNum];
  FieldsToInterpolate[1] = temperature;
 
  /* Allocate space for tracer particle write buffer. */
 
  int BufferSize = NumberOfTracerParticles*(GridRank+2);  // 2 fields + gridrank particle position files
  float *TracerParticleData = new float[BufferSize];
 
  /* Interpolate the values of these fields at the location of the
     tracer particles. */
 
  if (InterpolateTracerValues(ParticlePosition, &TracerParticleData,
			      ParticleType,
			      NumberOfParticles, LeftEdge, CellSize,
			      GridDimension, GridRank,
			      FieldsToInterpolate, 2, 1) == FAIL) {
    ENZO_FAIL("Error in InterpolateTracerValies.\n");
  }
 
  float *TracerParticleVelocityData;
  int PartVelBufferSize;

  if(TracerParticleOutputVelocity){
    /* allocate space for tracer particle velocity buffers */
    PartVelBufferSize = NumberOfTracerParticles*GridRank;
    TracerParticleVelocityData = new float[PartVelBufferSize];

    /* get updated velocity data */
    if (this->TracerParticleSetVelocity() == FAIL) {
      ENZO_FAIL("Error in grid->TracerParticleSetVelocity.\n");
    }

    /* Go through particle list.  We've already allocated a buffer that's the
       right size for our tracer particle data - now all we need to do is 
       sort through the entire list of particles, check to see which ones are 
       tracer particles, and for THOSE ones, add that velocity value to the 
       array of values, which we will write out later. */

    int numberofcountedparticles = 0;  // keep track of the total ptcls copied
    int fillindex = 0;
    for (n = 0; n < NumberOfParticles; n++)
      if (ParticleType[n] == PARTICLE_TYPE_TRACER){
	
	for (dim = 0; dim < GridRank; dim++) {
	  // set velocity data buffer (somewhat long-winded name)
	  TracerParticleVelocityData[fillindex] = float(ParticleVelocity[dim][n]);
	  fillindex++;
	}
	numberofcountedparticles++;  // increment counter!
      }

    // double check to make sure we didn't screw up
    if((numberofcountedparticles != NumberOfTracerParticles) || (fillindex/GridRank != NumberOfTracerParticles)){
      ENZO_VFAIL("error in grid->TracerParticleOutputData:  Particle num wrong! %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",
		 numberofcountedparticles, NumberOfTracerParticles, fillindex, GridRank)
	}

  } // if(TracerParticleOutputVelocity)

  /* Write data -- First, write number of tracer particles and then the
     data itself. */
 
  fwrite((void*) &NumberOfTracerParticles, sizeof(int), 1, fptr);
 
  /* Write floatpoint data (this is (2+GridRank) x NumberOfTracerParticles). */
 
  if (fwrite((void*) TracerParticleData, sizeof(float), BufferSize, fptr) != BufferSize) {
    ENZO_VFAIL("P(%"ISYM") Error writing tracer particle data\n",
	       MyProcessorNumber)
      }

  delete [] TracerParticleData;
 
  if(TracerParticleOutputVelocity){

    /* write more floating point data (which is GridRank x NumberOfTracerParticles).  This
       happens to be particle velocity. */
    if (fwrite((void*) TracerParticleVelocityData, sizeof(float), PartVelBufferSize, fptr) != PartVelBufferSize) {
      ENZO_VFAIL("P(%"ISYM") Error writing tracer particle data (2)\n", 
		 MyProcessorNumber)
	}

    delete [] TracerParticleVelocityData;

  } // if(TracerParticleOutputVelocity)

  /* Write tracer particle index data (this is 1 x NumberOfTracerParticles). */
 
  int *int_buffer = new int[NumberOfTracerParticles];
  int index = 0;
  for (n = 0; n < NumberOfParticles; n++)
    if (ParticleType[n] == PARTICLE_TYPE_TRACER){
      int_buffer[index++] = ParticleNumber[n];
      //printf("particle ID %d (%d)\n",ParticleNumber[n], int_buffer[index-1]);
    }

  fwrite((void*) int_buffer, sizeof(int), NumberOfTracerParticles, fptr);
  delete [] int_buffer;
 
  /* Clean up. */
 
  delete [] temperature;
 
  return SUCCESS;
}
 

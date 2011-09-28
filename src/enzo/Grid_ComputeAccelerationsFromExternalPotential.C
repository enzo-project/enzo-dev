/***********************************************************************
/
/  GRID CLASS (COMPUTE PARTICLE AND GRID ACCELERATIONS 
/               FROM AN EXTERNAL POTENTIAL FIELD)
/
/  written by: Elizabeth Tasker
/  date:       February, 2008
/  modified1:
/
/  PURPOSE: This routine overlaps considerably with Grid_ComputeAccelerations.C
/           and probably should be combined. It is currently seperate due
/           to problems with re-setting AccelerationField and PotentialField
/           in the wrong place.
/
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
extern "C" void FORTRAN_NAME(comp_accel)(float *source, float *dest1, 
			    float *dest2, float *dest3, int *ndim, int *iflag,
                           int *sdim1, int *sdim2, int *sdim3, int *ddim1, 
			    int *ddim2, int *ddim3, 
                           int *start1, int *start2, int *start3, 
			    float *delx, float *dely, float *delz);

extern "C" void FORTRAN_NAME(cic_interp)(FLOAT *posx, FLOAT *posy, 
			FLOAT *posz, int *ndim, int *npositions, 
                        float *sumfield, float *field, FLOAT *leftedge, 
                        int *dim1, int *dim2, int *dim3, FLOAT *cellsize);


int grid::ComputeAccelerationsFromExternalPotential(int DifferenceType, 
						    float *ExternalPotential, 
						    float *Field[])
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;


 int DiffFlag = 1, size = 1, i, dim, Offset[MAX_DIMENSION] = {0,0,0};
  float CellSize[MAX_DIMENSION] = {1,1,1};

  if (DifferenceType == PARTICLES || DifferenceType == ZEUS_GRIDS)
    DiffFlag = 0;

  /* Compute adot/a at time = t+1/2dt (time-centered). */

  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }

 
  /* Set cell size. */

  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    CellSize[dim] = a * GravitatingMassFieldCellSize;
    Offset[dim] = nint((CellLeftEdge[dim][0] - 
			GravitatingMassFieldLeftEdge[dim])/ CellWidth[dim][0]);
  }

  float *Acceleration[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    Acceleration[dim] = new float[size];

    for (i = 0; i < size; i++)
      Acceleration[dim][i] = 0.0;

  }


  /* Potential field for ZEUS and Particles have values at Grid edge. 
     Calculate acceleration field for these types also at grid edge */


  FORTRAN_NAME(comp_accel)(ExternalPotential, Acceleration[0], 
			   Acceleration[1], Acceleration[2], &GridRank, 
			   &DiffFlag, GravitatingMassFieldDimension, 
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   GridDimension, GridDimension+1, GridDimension+2,
			   Offset, Offset+1, Offset+2, CellSize, CellSize+1, 
			   CellSize+2);


 
  if (DifferenceType != PARTICLES) {
    
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < size; i++)
	Field[dim][i] = Acceleration[dim][i];
  }

  /* Interpolate acceleration field to get ParticleAcceleration */

  else if (DifferenceType == PARTICLES) 
    {
      
      this->UpdateParticlePosition(0.5*dtFixed);

      FLOAT LeftEdge[MAX_DIMENSION];
      for (dim = 0; dim < GridRank; dim++)
	LeftEdge[dim] = CellLeftEdge[dim][0]- 0.5*CellWidth[dim][0];

      float *PartAccel[GridRank];
      for (dim = 0; dim < GridRank; dim++){		

	PartAccel[dim] = new float[NumberOfParticles];
	
	for (i = 0; i < NumberOfParticles; i++)
	  PartAccel[dim][i] = 0.0;
      

	/* Adjust the grid position as the acceleration is face-centered 
	   for particles. (Needed for interpolation)*/

	for (i = 0; i < GridRank; i++)
	  LeftEdge[i] = CellLeftEdge[i][0] - ((dim == i)? (0.5*CellWidth[i][0]) : 0);


	FORTRAN_NAME(cic_interp)(ParticlePosition[0], ParticlePosition[1], 
				 ParticlePosition[2], &GridRank,
				 &NumberOfParticles, PartAccel[dim], 
				 Acceleration[dim], LeftEdge, GridDimension, 
				 GridDimension+1, GridDimension+2, 
				 &GravitatingMassFieldCellSize);	
      }

      for (dim = 0; dim < GridRank; dim++)
	for (i = 0; i < NumberOfParticles; i++)
	  Field[dim][i] = PartAccel[dim][i];

      this->UpdateParticlePosition(-0.5*dtFixed);
     

      /* Cleanup */

      // delete AccelPoints;
      for (dim = 0; dim < GridRank; dim++){
	delete [] PartAccel[dim];
	PartAccel[dim] == NULL;
      }
    

    }


  /* Cleanup */
  for (dim = 0; dim < GridRank; dim++){
    delete [] Acceleration[dim];
    Acceleration[dim] = NULL;
  }

  return SUCCESS;
}

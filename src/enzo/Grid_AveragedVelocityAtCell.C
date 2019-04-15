/***********************************************************************
/
/  RETURN SPATIALLY AVERAGED VELOCITY AT A CELL
/
/  written by: John Wise
/  date:       March, 2006
/  modified1:
/
/  PURPOSE: Used in the star formation routines to avoid "runaway" 
/           particles.
/
************************************************************************/
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"

float* grid::AveragedVelocityAtCell(int index, int DensNum, int Vel1Num)
{

  const int offset[] = {1, GridDimension[0], GridDimension[0]*GridDimension[1]};

  int i, dim, indices[7];
  float weight;
  float *vel = new float[MAX_DIMENSION];
  float *rho = BaryonField[DensNum];
  float *allv[] = {BaryonField[Vel1Num], 
		   BaryonField[Vel1Num+1], 
		   BaryonField[Vel1Num+2]};

  indices[0] = index;  //(i,j,k)
  indices[1] = index - offset[0]; //(i-1,j,k)
  indices[2] = index + offset[0]; //(i+1,j,k)
  indices[3] = index - offset[1]; //(i,j-1,k)
  indices[4] = index + offset[1]; //(i,j+1,k)
  indices[5] = index - offset[2]; //(i,j,k-1)
  indices[6] = index + offset[2]; //(i,j,k+1)

  if (HydroMethod == Zeus_Hydro) {

    for (dim = 0; dim < GridRank; dim++) {
      weight = 0;
      vel[dim] = 0;
      for (i = 0; i < 7; i++) {
	weight += rho[indices[i]];
	vel[dim] += 0.5 * (allv[dim][indices[i]] + allv[dim][indices[i]+offset[dim]]);
      }
      vel[dim] /= weight;
    } // ENDFOR dim
      
  } else {

    for (dim = 0; dim < GridRank; dim++) {
      weight = 0;
      vel[dim] = 0;
      for (i = 0; i < 7; i++) {
	weight += rho[indices[i]];
	vel[dim] += allv[dim][indices[i]] * rho[indices[i]];
      }
      vel[dim] /= weight;
    } // ENDFOR dim

  } // ENDELSE !Zeus

  return vel;

}

/***********************************************************************
/
/  RETURN MINIMUM OF ANGULAR MOMENTUM IN A CONTROL VOLUME
/
/  written by: John Regan
/  date:       March 2019
/  modified1:
/
/  PURPOSE: Used to test for minimum of the angular momentum
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
#include "phys_constants.h"

float grid::FindAngularMomentumMinimum(FLOAT *cellpos, FLOAT radius, int DensNum, int Vel1Num,
				       int Vel2Num, int Vel3Num)
{

  int i = 0, j = 0, k = 0;
  FLOAT pos[3];
  FLOAT dx = CellWidth[0][0];
  float AngularMinimum = 1e20;
  int index = 0;
  float* velx = BaryonField[Vel1Num];
  float* vely = BaryonField[Vel2Num];
  float* velz = BaryonField[Vel3Num];
  FLOAT AngularMomentum[3] = {0.0, 0.0, 0.0};
  double rhocell = 0.0, mcell = 0.0;
  double CellVolume = dx*dx*dx;
  /* Need to find gravitational minimum of this grid in advance */
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i, j, k);

	pos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	pos[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	/* Use only values within the control volume */
	if(sqrt((cellpos[0] - pos[0])*(cellpos[0] - pos[0]) +
		(cellpos[1] - pos[1])*(cellpos[1] - pos[1]) +
		(cellpos[2] - pos[2])*(cellpos[2] - pos[2])) > radius){
	  continue;
	}
	rhocell = BaryonField[DensNum][index];
	mcell = rhocell*CellVolume;
	AngularMomentum[0] = mcell*(pos[1]*velz[index] -
				     pos[2]*vely[index]);
	AngularMomentum[1] = mcell*(pos[2]*velx[index]-
				     pos[0]*velz[index]);
	AngularMomentum[2] = mcell*(pos[0]*vely[index] -
				     pos[1]*velx[index]);
	
	FLOAT AngularMomentumMag = sqrt(AngularMomentum[0]*AngularMomentum[0] +
			       AngularMomentum[1]*AngularMomentum[1] +
			       AngularMomentum[2]*AngularMomentum[2]); 
	

	{
	  if(AngularMomentumMag < AngularMinimum) {
	    AngularMinimum = AngularMomentumMag;
	  }
	}
      }
    }
  }

  return AngularMinimum;
}

/***********************************************************************
/
/  RETURN MINIMUM OF GRAVITATIONAL POTENTIAL IN A CONTROL VOLUME
/
/  written by: John Regan
/  date:       January 2017
/  modified1:
/
/  PURPOSE: Used to test for minimum of the potential
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

float grid::FindMinimumPotential(FLOAT *cellpos, FLOAT radius, float *PotentialField)
{

  int i = 0, j = 0, k = 0;
  FLOAT pos[3];
  FLOAT dx = CellWidth[0][0];
  float GravitationalMinimum = 1e20;
  int index = 0;
 
  /* Need to find gravitational minimum of this grid in advance */
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i, j, k);
	if(PotentialField[index] != PotentialField[index])
	  PotentialField[index] = 0.0;
	pos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	pos[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	
	/* Use only values within the control volume */
	if(sqrt((cellpos[0] - pos[0])*(cellpos[0] - pos[0]) +
		(cellpos[1] - pos[1])*(cellpos[1] - pos[1]) +
		(cellpos[2] - pos[2])*(cellpos[2] - pos[2])) <= radius)
	{
	  if(PotentialField[index] < GravitationalMinimum)
	    GravitationalMinimum = PotentialField[index];
	}
      }
    }
  }

  return GravitationalMinimum;
}

/*
 * Calculate the Potential Field in this grid by direct summation
 * This is an N^2 problem - be careful how frequently you call it. 
 * Phi = -G m_j / (x - x_j)
 */
void grid::CalculatePotentialField(float *PotentialField, int DensNum, float DensityUnits,
				   float TimeUnits, float LengthUnits)
{
  int i = 0, j = 0, k = 0, ii = 0, jj = 0, kk = 0;
  FLOAT sep = 0;
  FLOAT dx = CellWidth[0][0], pos[3], pos_j[3];
  float MassUnits = DensityUnits*LengthUnits*LengthUnits*LengthUnits;
  //float G = 4*M_PI*GravConst*DensityUnits*TimeUnits*TimeUnits;
  float G = GravitationalConstant;
  float mass_j = 0;
  int index = 0, lindex = 0;
  float *density = BaryonField[DensNum];
  
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i, j, k);
	pos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	pos[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	

	for (kk = GridStartIndex[2]; kk <= GridEndIndex[2]; kk++) {
	  for (jj = GridStartIndex[1]; jj <= GridEndIndex[1]; jj++) {
	    for (ii = GridStartIndex[0]; ii <= GridEndIndex[0]; ii++) {
	      lindex = GRIDINDEX_NOGHOST(ii, jj, kk);
	      mass_j = density[lindex]*dx*dx*dx;
	      pos_j[0] = CellLeftEdge[0][ii] + 0.5*CellWidth[0][ii];
	      pos_j[1] = CellLeftEdge[1][jj] + 0.5*CellWidth[1][jj];
	      pos_j[2] = CellLeftEdge[2][kk] + 0.5*CellWidth[2][kk];
	      sep = sqrt((pos[0] - pos_j[0])*(pos[0] - pos_j[0]) +
			 (pos[1] - pos_j[1])*(pos[1] - pos_j[1]) +
			 (pos[2] - pos_j[2])*(pos[2] - pos_j[2]));
	      if(index == lindex) {
		sep = CellWidth[0][ii]/2.0;
	      }
	      PotentialField[index] += -G * mass_j/sep;
	    }
	  }
	}
      }
    }
  }
  return;
}

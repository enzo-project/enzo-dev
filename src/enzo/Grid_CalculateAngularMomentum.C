/***********************************************************************
/
/  GRID CLASS (COMPUTE THE ANGULAR MOMENTUM, GIVEN CENTER)
/
/  written by: Greg Bryan
/  date:       December, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "CosmologyParameters.h"
#include "Grid.h"
 
 
int grid::CalculateAngularMomentum(FLOAT Center[], float AngularMomentum[],
				   float MeanVelocity[], float DMVelocity[],
				   FLOAT CenterOfMass[], FLOAT DMCofM[])
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Declarations */
 
  float mass, TotalMass = 0;
  FLOAT xpos, ypos, zpos;
  int dim, i, j, k, index;
 
  /* Set a.m. to zero. */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    AngularMomentum[dim] = 0;
    MeanVelocity[dim] = 0;
    CenterOfMass[dim] = 0;
  }
 
  /* Only works for dim=3. */
 
  if (GridRank != 3)
    return SUCCESS;
 
  /* Compute Angular Momentum for the baryons. */
 
  if (NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }
 
    /* Loop over grid. */
 
    float vel[] = {0, 0, 0};
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      zpos = Center[2] - (CellLeftEdge[2][k] + 0.5*CellWidth[2][k]);
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	ypos = Center[1] - (CellLeftEdge[1][j] + 0.5*CellWidth[1][j]);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	  xpos = Center[0] - (CellLeftEdge[0][i] + 0.5*CellWidth[0][i]);
 
	  index = k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i;
 
	  mass = BaryonField[DensNum][index]*
	         CellWidth[0][i]*CellWidth[1][j]*CellWidth[2][k];
 
	  /* Get Velocity (average to cell-center for Zeus). */
 
	  int mult = 1;
	  for (dim = 0; dim < GridRank; dim++) {
	    vel[dim] = BaryonField[Vel1Num+dim][index];
	    if (HydroMethod == Zeus_Hydro)
	      vel[dim] = 0.5*(vel[dim] +
			      BaryonField[Vel1Num+dim][index+1*mult]);
	    mult *= GridDimension[dim];
	  }
 
	  for (dim = 0; dim < GridRank; dim++)
	    MeanVelocity[dim] += mass*vel[dim];
	  TotalMass += mass;
 
	  CenterOfMass[0] += mass * xpos;
	  CenterOfMass[1] += mass * ypos;
	  CenterOfMass[2] += mass * zpos;
 
	  AngularMomentum[0] += mass * (ypos*vel[2] - zpos*vel[1]);
	  AngularMomentum[1] += mass * (zpos*vel[0] - xpos*vel[2]);
	  AngularMomentum[2] += mass * (xpos*vel[1] - ypos*vel[0]);
 
/*	  if (AngularMomentum[0] != 0)
	    printf("%"ISYM" %10.5"GSYM" %10.5"GSYM" %10.5"GSYM"  %10.5"GSYM" %10.5"GSYM" %10.5"GSYM"  %10.5"GSYM" %10.5"GSYM" %10.5"GSYM"\n", index, AngularMomentum[0],
		   AngularMomentum[1], AngularMomentum[2], xpos, ypos, zpos,
		   BaryonField[Vel1Num][index], BaryonField[Vel2Num][index],
		   BaryonField[Vel3Num][index]);
*/
 
	}
      }
    }
 
    for (dim = 0; dim < GridRank; dim++) {
      MeanVelocity[dim] /= TotalMass;
      CenterOfMass[dim] /= TotalMass;
    }
 
  } // end: if (NumberOfBaryonFields > 0)
 
  /* Compute Particle mean velocity */
 
  float DMMass;
  for (dim = 0; dim < GridRank; dim++) {
    DMVelocity[dim] = DMCofM[dim] = DMMass = 0;
    for (i = 0; i < NumberOfParticles; i++) {
      DMVelocity[dim] += ParticleVelocity[dim][i]*ParticleMass[i];
      DMCofM[dim]     += ParticlePosition[dim][i]*ParticleMass[i];
      DMMass          += ParticleMass[i];
    }
    if (DMMass > 0) {

      DMVelocity[dim] /= DMMass;
      DMCofM[dim]     /= DMMass;
    }
  }
 
  return SUCCESS;
}

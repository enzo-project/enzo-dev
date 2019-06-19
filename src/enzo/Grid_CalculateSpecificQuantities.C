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
#include "phys_constants.h"
#define N 8
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int grid::CalculateSpecificQuantities(FLOAT *SinkParticlePos, FLOAT *CLEdge,
				      float *vgas, float msink,
				      float *vsink, int *numpoints)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1,
    PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  double SpecificAngularMomentum[3][N*N*N];
  double SpecificEnergy[N*N*N];
  FLOAT xpos = 0.0, ypos = 0.0, zpos = 0.0;
  int dim = 0, i = 0, j = 0, k = 0, index = 0;
  double Gcode = GravConst*DensityUnits*TimeUnits*TimeUnits;
  /* Set to zero. */
 for(index = 0; index < N*N*N; index++) {
   SpecificEnergy[index] = 0.0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
      SpecificAngularMomentum[dim][index] = 0.0;
    }
  }
 double A = 0.0, B = 0.0;
  /* Compute Angular Momentum for the baryons. */
 double jspec = 0.0, espec = 0.0;
 FLOAT deltax = CellWidth[0][0]/ (FLOAT)N;
 double xvel = 0.0, yvel = 0.0, zvel = 0.0;
  for(k = 0; k < N; k++) {
    zpos = SinkParticlePos[2] - (CLEdge[2] + k*deltax);
    zvel = vsink[2] - vgas[2];
    for(j = 0; j < N; j++) {
      ypos = SinkParticlePos[1] - (CLEdge[1] + j*deltax);
      yvel = vsink[1] - vgas[1];
      for(i = 0; i < N; i++) {
	xpos = SinkParticlePos[0] - (CLEdge[0] + i*deltax);
	xvel = vsink[0] - vgas[0];
	index = k*N*N + j*N + i;
	
	SpecificAngularMomentum[0][index] = (ypos*zvel - zpos*yvel);
	SpecificAngularMomentum[1][index] = (zpos*xvel - xpos*zvel);
	SpecificAngularMomentum[2][index] = (xpos*yvel - ypos*xvel);
	jspec = sqrt(SpecificAngularMomentum[0][index]*SpecificAngularMomentum[0][index] + 
		     SpecificAngularMomentum[1][index]*SpecificAngularMomentum[1][index] +
		     SpecificAngularMomentum[2][index]*SpecificAngularMomentum[2][index]);

	/* KE */
	SpecificEnergy[index] = (xvel*xvel + 
				 yvel*yvel + 
				 zvel*zvel);
	
	/* plus PE */
	double r = sqrt(xpos*xpos * ypos*ypos + zpos*zpos);
	SpecificEnergy[index] -= Gcode*msink/r;
	espec = SpecificEnergy[index];
	
	double rmin = 0.0;
	if(espec > 0.0) /* The point is not bound to the sink */
	  rmin = huge_number;
	else {
	  B = 2*jspec*espec/((Gcode*msink)*(Gcode*msink));
	  B = max(-1.0, B);
	  A = (1.0 - sqrt(1 + B));
	  rmin = -Gcode*msink*A/(2.0*espec);
	}
	if(rmin != huge_number){
	  // printf("espec = %e\t", espec);
	  // printf("rmin = %e\t CellWidth[0][0]/4.0 = %e\n", rmin, CellWidth[0][0]/4.0);
	  // printf("Gcode = %e\n", Gcode);
	  // printf("A = %e\n", A);
	  // printf("msink = %e\n", msink);
	  // printf("jspec = %e\n", jspec);
	  // printf("2*jspec*espec = %e\n", 2*jspec*espec);
	  // printf("((Gcode*msink)*(Gcode*msink)) = %e\n", ((Gcode*msink)*(Gcode*msink)));
	  // printf("2*jspec*espec/((Gcode*msink)*(Gcode*msink))) = %e\n",
	  // 	 2*jspec*espec/((Gcode*msink)*(Gcode*msink)));
	  // printf("B = %e\n", B);
	}
	if(rmin <= (double)CellWidth[0][0]/4.0) /* going to be counted */
	  *numpoints++;
      }
    }
  }    
    
  return SUCCESS;
}

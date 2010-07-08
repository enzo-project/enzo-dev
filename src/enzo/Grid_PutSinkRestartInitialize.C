/***********************************************************************
/
/  GRID CLASS (RE-INITIALIZE FOR PUTTING IN SINK PARTICLES)
/
/  written by: Elizabeth Harper-Clark
/  date:       March, 2010
/  modified1:
/
/  PURPOSE: based on SupernovaeRestartInitialize
/
/  RETURNS: FAIL or SUCCESS
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
 int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int grid::PutSinkRestartInitialize(int level, int *NumberOfCellsSet)
{
  /* declarations */
 
  int dim, i, j, k,l, n;
  FLOAT delx, dely, delz, radius2, DomainWidth[MAX_DIMENSION];

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, TimeUnits = 1.0,
    VelocityUnits = 1.0;
  if (UsePhysicalUnit)
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  double MassUnits = DensityUnits*pow(LengthUnits,3);
  printf("Mass Units = %g \n",MassUnits);
  printf("Time Units = %g \n",TimeUnits);


  for (dim = 0; dim < GridRank; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  if (level == 0) {  // set it up on level zero and make it mustrefine




    NumberOfParticles = 64;
    NumberOfStars = 64;
    printf("Adding Sink Particles. \n");
    NumberOfParticleAttributes = 3;
    if (StellarWindFeedback) NumberOfParticleAttributes = 6;
    this->AllocateNewParticles(NumberOfParticles);
    double mass_m = 3.415*1.989e33; //Mass of massive stars
    double mass_s = 0.01*1.989e33; //Mass of small stars
    mass_m /= MassUnits;
    mass_s /= MassUnits;
    double dx = CellWidth[0][0];
    double den_m = mass_m / pow(dx,3);
    double den_s = mass_s / pow(dx,3);
    double t_dyn_m = sqrt(3*M_PI/(6.672e-8*den_m*DensityUnits));
    double t_dyn_s = sqrt(3*M_PI/(6.672e-8*den_s*DensityUnits));
    t_dyn_m /= TimeUnits;
    t_dyn_s /= TimeUnits;
    double dxm = dx / pow(RefineBy, MaximumRefinementLevel);
    for (k=0; k<4; k++){
      for (j=0; j<4; j++){
	for (i=0; i<4; i++){
	  l = i+4*j+16*k;
	  printf("Creating particle %i \n",l);
	  ParticleMass[l] = den_m;
	  ParticleNumber[l] = l;
	  ParticleType[l] = PARTICLE_TYPE_MUST_REFINE;
	  ParticlePosition[0][l] = 0.125+0.25*i+0.5*dxm;
	  ParticlePosition[1][l] = 0.125+0.25*j+0.5*dxm;
	  ParticlePosition[2][l] = 0.125+0.25*k+0.5*dxm;
	  ParticleVelocity[0][l] = 0.0;
	  ParticleVelocity[1][l] = 0.0;
	  ParticleVelocity[2][l] = 0.0;
	  ParticleAcceleration[0] = NULL;
	  ParticleAcceleration[1] = NULL;
	  ParticleAcceleration[2] = NULL;
	  ParticleAttribute[0][l] = 0.0; // creation time             
	  ParticleAttribute[1][l] = t_dyn_m; // t_dyn
	  ParticleAttribute[2][l] = mass_m;
	  if (StellarWindFeedback) {
	    ParticleAttribute[3][l] = 1.0;  
	    ParticleAttribute[4][l] = 0.0;
	    ParticleAttribute[5][l] = 0.0;
	  }
	  this->ClearParticleAccelerations();
	  printf("Completed particle %i, position %g,%g,%g \n",l,ParticlePosition[0][l],ParticlePosition[1][l],ParticlePosition[2][l]);
	}
      }
    }
  }



















  /* Loop over grid and set quantities. */
 
//   n = 0;
//   for (k = 0; k < GridDimension[2]; k++) {
//     if (GridRank > 1) {
//       delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - EjectaCenter[2];
//       delz = min(delz, DomainWidth[2]-delz);
//     }
//     else
//       delz = 0;
 
//     for (j = 0; j < GridDimension[1]; j++) {
//       if (GridRank > 0) {
// 	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - EjectaCenter[1];
// 	dely = min(dely, DomainWidth[1]-dely);
//       }
//       else
// 	dely = 0;
 
//       for (i = 0; i < GridDimension[0]; i++, n++) {
// 	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - EjectaCenter[0];
// 	delx = min(delx, DomainWidth[0]-delx);
 
// 	/* Compute square of distance from cell to center. */
 
// 	radius2 = delx*delx + dely*dely + delz*delz;
 
// 	if (radius2 <= EjectaRadius*EjectaRadius*1.2*1.2) {
 
// 	  float r1 = sqrt(radius2)/EjectaRadius;
// 	  float ramp = min(max(1.0 - (r1 - 0.8)/0.4, 0.01), 1.0);
 
// 	  /* Cell is within ejecta, so set density, etc. */
 
// 	  BaryonField[DensNum][n] = EjectaDensity;
// 	  BaryonField[TENum][n]   = EjectaThermalEnergy*ramp;
// 	  if (GENum >= 0)
// 	    BaryonField[GENum][n]   = BaryonField[TENum][n];
// 	  for (dim = 0; dim < GridRank; dim++)
// 	    BaryonField[Vel1Num+dim][n] = 0;
// 	  (*NumberOfCellsSet)++;
 
// 	} // end: if (radius2 < EjectaRadius*EjectaRadius)

 
//       } // next i
//     } // next j
//   } // next k
 
 
  return SUCCESS;
}

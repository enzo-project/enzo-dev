/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR FREE EXPANSION BLAST WAVE) 
/
/  written by: John Wise
/  date:       August, 2009
/  modified1:  
/
/  PURPOSE: Sets the velocity and density in a freely expanding blast wave.
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

int grid::FreeExpansionInitializeGrid(int FreeExpansionFullBox,
				      float FreeExpansionDensity,
				      double FreeExpansionEnergy,
				      float FreeExpansionMaxVelocity,
				      float FreeExpansionMass,
				      float FreeExpansionRadius,
				      float DensityUnits, float VelocityUnits,
				      float LengthUnits, float TimeUnits)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  const float DensitySlope = 9.0;  // density decreases as (v/vcore)^n
				   // outside the core.
  const double Msun = 1.989e33, pc = 3.086e18, G = 6.673e-8;
  const double mh = 1.673e-24;

  int i, j, k, dim, index;
  float delx, dely, delz, r2, radius, router2, speed, density;
  FLOAT Center[] = {0,0,0};

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
			     Vel3Num, TENum);

  /* Compute some normalization factors, core and maximum velocity,
     time, and radius where the free expansion stage ends */

  double R_M, v_max, BlastTime, v_core, rho_0, r_max, normalization;
  double M_ej, E_ej;
  
  rho_0 = FreeExpansionDensity * DensityUnits;
  M_ej = FreeExpansionMass * Msun;
  E_ej = FreeExpansionEnergy;
  R_M = POW((3 * M_ej) / (4 * M_PI * rho_0), 1.0/3);
  r_max = FreeExpansionRadius * LengthUnits;

  if (FreeExpansionMaxVelocity == FLOAT_UNDEFINED)
    v_max = 1.165 * sqrt(2.0 * E_ej / M_ej) / (1.0 + POW(r_max / R_M, 3.0)) 
      / VelocityUnits;
  else
    v_max = FreeExpansionMaxVelocity * 1e5 / VelocityUnits;

  BlastTime = (r_max/LengthUnits) / v_max;
  v_core = sqrt( (10.0 * E_ej * (DensitySlope-5)) /
		 (3.0 * M_ej * (DensitySlope-3)) );
  normalization = (10.0 * (DensitySlope-5) * E_ej) / (4.0 * M_PI * DensitySlope) /
    POW(v_core, 5.0);
    

  if (FreeExpansionFullBox)
    for (dim = 0; dim < GridRank; dim++)
      Center[dim] = 0.5;

  // Absorb units
  v_core /= VelocityUnits;

  router2 = FreeExpansionRadius * FreeExpansionRadius;

  for (k = 0; k < GridDimension[2]; k++) {
    delz = (GridRank > 2) ? (CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - 
			     Center[2]) : 0;
    for (j = 0; j < GridDimension[1]; j++) {
      dely = (GridRank > 1) ? (CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - 
			       Center[1]) : 0;
      index = GRIDINDEX_NOGHOST(0,j,k);
      for (i = 0; i < GridDimension[0]; i++, index++) {
	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - Center[0];
	r2 = delx*delx + dely*dely + delz*delz;
	if (r2 < router2) {

	  radius = sqrt(r2);
	  speed = radius / BlastTime;
	  if (speed <= v_core)
	    density = normalization / POW(BlastTime*TimeUnits, 3.0) / DensityUnits;
	  else if (speed <= v_max)
	    density = normalization / POW(BlastTime*TimeUnits, 3.0) / 
	      POW(speed/v_core, DensitySlope) / DensityUnits;
	  density = max(density, FreeExpansionDensity);

	  BaryonField[DensNum][index] = density;
	  BaryonField[Vel1Num][index] = speed * delx / radius;
	  if (GridRank > 1)
	    BaryonField[Vel2Num][index] = speed * dely / radius;
	  if (GridRank > 2)
	    BaryonField[Vel3Num][index] = speed * delz / radius;
	  BaryonField[TENum][index] = 0.5 * speed * speed;
	  
	} // ENDIF r2 < router2
      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  return SUCCESS;
}

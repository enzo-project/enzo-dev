/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (SETS OUTFLOW BOUNDARY CONDITIONS FOR WAVE POOL)
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
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
#include "WavePoolGlobalData.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
// Set the Left BoundaryValue of the chosen wave direction (set by
//  WavePoolAngle) to the appropriate inflow boundary condition.
//  This is a sinusoidal wave
 
int ExternalBoundary::SetWavePoolBoundary(FLOAT time)
{
  /* declarations */
 
  int i, j, dim, index;
  int NumberOfZones[MAX_DIMENSION], Offset[MAX_DIMENSION];
  float Omega, Perturbation, deltime, distance, pos[MAX_DIMENSION];
  float WaveDensity, WaveTotalEnergy, WavePressure, WaveVelocity;
  float WaveAverageVelocity, WaveVelocity123[MAX_DIMENSION];
  const float TwoPi = 6.283185;
 
  /* Compute the velocity and time constant for the wave. */
 
  WaveAverageVelocity = sqrt(Gamma*WavePoolPressure/WavePoolDensity);
  Omega               = TwoPi*WaveAverageVelocity/WavePoolWavelength;
 
  /* Compute size of entire mesh. */
 
  int size = 1;
  for (dim = 0; dim < BoundaryRank; dim++)
    size = size*BoundaryDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* set the appropriate BoundaryValues on the left side */
 
  for (dim = 0; dim < BoundaryRank; dim++)
    if (BoundaryDimension[dim] != 1) {
 
      /* If the BoundaryValue fields are missing, create them. */
 
      for (int field = 0; field < NumberOfBaryonFields; field++)
	if (BoundaryValue[field][dim][0] == NULL)
	  BoundaryValue[field][dim][0] =
	    new float[size/BoundaryDimension[dim]];
 
      /* Compute quantities needed for boundary face loop (below). */
 
      int dim1, dim2;
      dim1 = (dim == 0) ? 1 : 0;
      dim2 = dim1 + 1;
      dim2 = (dim2 == dim) ? dim2+1 : dim2;
      for (i = 0; i < 3; i++) {
	NumberOfZones[i] = max(BoundaryDimension[i] - 2*NumberOfGhostZones,1);
	Offset[i]        = min(NumberOfGhostZones, BoundaryDimension[i]) - 1;
      }
      pos[dim] = 0.0;
 
      /* Loop over the boundary face. */
 
      for (i = 0; i < BoundaryDimension[dim1]; i++)
	for (j = 0; j < BoundaryDimension[dim2]; j++) {
 
	  /* Find the 3D vector from the corner to the current location. */
 
	  pos[dim1] = (float(i-Offset[dim1]))*
	    (DomainRightEdge[dim1]-DomainLeftEdge[dim1]) /
	      float(NumberOfZones[dim1]);
	  pos[dim2] = (float(j-Offset[dim2]))*
	    (DomainRightEdge[dim2]-DomainLeftEdge[dim2]) /
	      float(NumberOfZones[dim2]);
 
	  /* Compute the distance along the wave propogation vector
	     (cos(angle), sin(angle), 0). Convert to radians. */
	
	  distance = pos[0]*cos(WavePoolAngle*TwoPi/360.0) +
	             pos[1]*sin(WavePoolAngle*TwoPi/360.0);
 
	  /* Find the difference between the current time and the time at
	     which the wave will reach this point. */
 
	  deltime = time - distance/WaveAverageVelocity;
 
	  /* If this point is currently inside the wave-train, compute the
	     appropriate density, velocity and energy. */
 
	  if (deltime*WaveAverageVelocity <
	      WavePoolWavelength*float(WavePoolNumberOfWaves) &&
	      deltime > 0.0) {
	    Perturbation    =  WavePoolAmplitude *
	                       sin(Omega*deltime) / WavePoolDensity;
	    WaveDensity     =  WavePoolDensity  *    (1+Perturbation);
	    WavePressure    =  WavePoolPressure * POW(1+Perturbation, Gamma);
	    WaveVelocity    =  WaveAverageVelocity *      Perturbation;
	    WaveTotalEnergy =  WavePressure/((Gamma-1)*WaveDensity);
	    if (HydroMethod != Zeus_Hydro)
	      WaveTotalEnergy += 0.5*WaveVelocity*WaveVelocity;
	    WaveVelocity123[0] = WaveVelocity*cos(WavePoolAngle*TwoPi/360.0);
	    WaveVelocity123[1] = WaveVelocity*sin(WavePoolAngle*TwoPi/360.0);
	    WaveVelocity123[2] = 0.0;
	  } else {
 
	    /* If not, set the fields to their ambient values. */
 
	    WaveDensity     = WavePoolDensity;
	    WaveTotalEnergy = WavePoolPressure/((Gamma-1)*WavePoolDensity);
	    WaveVelocity123[0] = 0.0;
	    WaveVelocity123[1] = 0.0;
	    WaveVelocity123[2] = 0.0;
	  }
 
	  /* Set the field values. */
 
	  index = j*BoundaryDimension[dim1] + i;
 
	  *(BoundaryValue[DensNum][dim][0] + index) = WaveDensity;
	  *(BoundaryValue[TENum][dim][0]   + index) = WaveTotalEnergy;
	  *(BoundaryValue[Vel1Num][dim][0] + index) = WaveVelocity123[0];
	  if (BoundaryRank > 1)
	    *(BoundaryValue[Vel2Num][dim][0] + index) = WaveVelocity123[1];
	  if (BoundaryRank > 2)

	    *(BoundaryValue[Vel3Num][dim][0] + index) = WaveVelocity123[2];
 
	} // end loop over boundary slice
 
    } // end loop over boundary directions
 
  return SUCCESS;
 
}

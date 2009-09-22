/***********************************************************************
/
/  ADD H2 DISSOCIATION EMISSION FROM SHINING PARTICLES
/
/  written by: John Wise
/  date:       March, 2006
/  modified1:
/
/  PURPOSE:
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
#include "Hierarchy.h"
#include "CosmologyParameters.h"
#include "Star.h"

int FindField(int f, int farray[], int n);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::AddRadiationImpulse(int field, double Luminosity, double sigma, 
			      FLOAT BirthTime, FLOAT* pos)
{

  const double pc = 3.086e19, clight = 3e10;

  int i, j, k, index, dim, FieldNum;
  float kdiss_r2;
  FLOAT delx, dely, delz, DomainWidth[MAX_DIMENSION];
  FLOAT innerFront, outerFront;
  double r2, radius2;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if ((FieldNum = FindField(field, FieldType, NumberOfBaryonFields)) < 0)
    ENZO_FAIL("Failed to find radiation field.");

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, 
    TimeUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Absorb the unit conversions into the cross-section
  sigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);

  for (dim = 0; dim < GridRank; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

  // Dilution factor (prevent breaking in the rate solver near the star)
  float dilutionRadius = 10.0 * pc / (double) LengthUnits;
  float dilRadius2 = dilutionRadius * dilutionRadius;
  float LightTravelDist = TimeUnits * clight / LengthUnits;

  // Determine the inner and outer radiation fronts
  innerFront = 0;
  outerFront = (PhotonTime - BirthTime) * LightTravelDist;
      
  /* Loop over cells */

  index = 0;
  kdiss_r2 = (float) (Luminosity * sigma / (4.0 * M_PI));
  for (k = 0; k < GridDimension[2]; k++) {
    delz = fabs(CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - pos[2]);
    delz = min(delz, DomainWidth[2]-delz);
    for (j = 0; j < GridDimension[1]; j++) {
      dely = fabs(CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - pos[1]);
      dely = min(dely, DomainWidth[1]-dely);
      for (i = 0; i < GridDimension[0]; i++, index++) {
	delx = fabs(CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - pos[0]);
	delx = min(delx, DomainWidth[0]-delx);
	radius2 = delx*delx + dely*dely + delz*delz;
	if (radius2 > outerFront*outerFront || radius2 < innerFront*innerFront)
	  continue;
	
	radius2 = max(radius2, dilRadius2);
	BaryonField[FieldNum][index] += kdiss_r2 / radius2;

      } // END: i-direction
    } // END: j-direction
  } // END: k-direction

  return SUCCESS;

}

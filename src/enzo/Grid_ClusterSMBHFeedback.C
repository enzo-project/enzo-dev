/***********************************************************************
/
/  GRID: ADD RADIO-MODE JET-LIKE FEEDBACK BASED ON STATIC SMBH
/
/  written by: Yuan Li and Greg Bryan
/  date:       December, 2011
/  modified1: 
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
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
///#include ""
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt); // need this?

int grid::ClusterSMBHFeedback(int level)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Return if not on most-refined level. */

  if (level != MaximumRefinementLevel)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                                             Vel3Num, TENum) == FAIL)   ///this or thisgrid
     ENZO_FAIL("Error in IdentifyPhysicalQuantities.");


  /* Compute the jet launching region
     (assume jet launched from PointSourceGravityPosition) */

  FLOAT JetLeftCorner[MAX_DIMENSION], JetRightCorner[MAX_DIMENSION];
  FLOAT JetCenter[MAX_DIMENSION];

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
    TimeUnits = 1.0, VelocityUnits = 1.0;
  double MassUnits = 1.0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, 0.0) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }
  MassUnits = DensityUnits*pow(LengthUnits,3);
  VelocityUnits = LengthUnits/TimeUnits;   // need this??

  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.67e-8,
               pi = 3.14159, mh = 1.67e-24, kboltz = 1.381e-16;
  int i, j, k = 0, size = 1, dim,  cellindex;
  int jet_dim = 2;  // z-axis (should make parameter?)

  int JetLaunchOffset = 10; // 10 cellwidth
  int JetRadius = 6; // cellwidths
  float JetScaleRadius = 3.0; // cellwidths

  float JetMdot = 10.0; // Jet mass flow in SolarMass/year (need to convert units)
  float JetVelocity = 1000.0; // Jet Velocity in km/s (should make parameter)

  for (dim = 0; dim < GridRank; dim++) {
    JetCenter[dim] = PointSourceGravityPosition[dim];
    JetLeftCorner[dim] = JetCenter[dim];
    JetRightCorner[dim] = JetCenter[dim];
    if (dim != jet_dim) {
      JetLeftCorner[dim] -= JetRadius*CellWidth[dim][0];
      JetRightCorner[dim] += JetRadius*CellWidth[dim][0];
    }
  }

  JetLeftCorner[jet_dim] -= JetLaunchOffset*CellWidth[jet_dim][0];
  JetRightCorner[jet_dim] += JetLaunchOffset*CellWidth[jet_dim][0];
 
  /* Compute indices of jet launch region. */

  int JetStartIndex[MAX_DIMENSION], JetEndIndex[MAX_DIMENSION];

  for (dim = 0; dim < GridRank; dim++) {

    /* Compute start and end indices of jet */

    JetStartIndex[dim] = nint((JetLeftCorner[dim] - CellLeftEdge[dim][0])/CellWidth[dim][0]);
    JetEndIndex[dim] = nint((JetRightCorner[dim] - CellLeftEdge[dim][0])/CellWidth[dim][0]);

    /* If Jet is not on this grid, return. */

    if (JetStartIndex[dim] > GridDimension[dim]-1 || JetEndIndex[dim] < 0)
      return SUCCESS;

  } // end: loop over dim

  /* Compute mass and momentum to be put into cells in code units. */

  float JetNormalization = 1.0, density_normalization = 1.0, radius = 1.0;
  for (j = JetStartIndex[1]; j <= JetEndIndex[1]; j++) {
    for (i = JetStartIndex[0]; i <= JetEndIndex[0]; i++) {
      radius = sqrt(pow((CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - JetCenter[0]), 2) + 
		    pow((CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - JetCenter[1]), 2))
	/CellWidth[0][0]; // in cell widths
	
      JetNormalization += exp(-pow(radius/JetScaleRadius,2));
    }
  }
  JetMdot = (JetMdot*SolarMass/MassUnits)/3.1557e7;  // in code units
  JetVelocity = JetVelocity*1.0e5/VelocityUnits; //from km/s to code units
  density_normalization = (JetMdot/JetNormalization)*CellWidth[0][0]/JetVelocity;

  /* Clip edge of jet launching disk so we don't set cell off the edge of the grid. */


  for (dim = 0; dim < GridRank; dim++) {
    if (dim != jet_dim) {
      JetStartIndex[dim] = max(JetStartIndex[dim], 0);
      JetEndIndex[dim] = min(JetEndIndex[dim], GridDimension[dim]-1);
    }
  }

  /* Loop over launch disks and set cell values (this code assumes jet_dim = 2). */

  if (JetStartIndex[jet_dim] >= 0) {
    k = JetStartIndex[jet_dim];
    for (j = JetStartIndex[1]; j <= JetEndIndex[1]; j++) {
      for (i = JetStartIndex[0]; i <= JetEndIndex[0]; i++) {
	radius = sqrt(pow((CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - JetCenter[0]), 2) + 
		      pow((CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - JetCenter[1]), 2))
	  /CellWidth[0][0]; // in cell widths
	
	BaryonField[DensNum][GRIDINDEX(i,j,k)] += density_normalization*exp(-pow(radius/JetScaleRadius,2));
	BaryonField[Vel3Num][GRIDINDEX(i,j,k)] += JetVelocity;
	//	BaryonField[GENum][GRIDINDEX(i,j,k)] += XXX;
      }
    }
  }

    

  /* loop over cells to be modified, add jet mass, momentum, and energy. */

  return SUCCESS;

}
 

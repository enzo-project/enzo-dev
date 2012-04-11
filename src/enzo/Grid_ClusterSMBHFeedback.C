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
#include "phys_constants.h"
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);


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
               &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }
  MassUnits = DensityUnits*pow(LengthUnits,3);

  int i, j, k, dim = 0;
  int jet_dim = 2;  // z-axis (should make parameter?)

  float JetScaleRadius; // cellwidths
  float JetMdot; // Jet mass flow in SolarMass/year (need to convert units)-- gets value from parameter ClusterSMBHJetMdot
  float JetVelocity, FastJetVelocity; // Jet Velocity in km/s (should make parameter)-- gets value from parameter ClusterSMBHJetVelocity
  JetScaleRadius = ClusterSMBHJetRadius/2.0;  //JetScaleRadius is half the radius of the jet launch region in cellwidths
  for (dim = 0; dim < GridRank; dim++) {
    JetCenter[dim] = PointSourceGravityPosition[dim];
    JetLeftCorner[dim] = JetCenter[dim];
    JetRightCorner[dim] = JetCenter[dim];
    if (dim != jet_dim) {
      JetLeftCorner[dim] -= ClusterSMBHJetRadius*CellWidth[dim][0];
      JetRightCorner[dim] += ClusterSMBHJetRadius*CellWidth[dim][0];
    }
  }

  JetLeftCorner[jet_dim] -= ClusterSMBHJetLaunchOffset*CellWidth[jet_dim][0];
  JetRightCorner[jet_dim] += ClusterSMBHJetLaunchOffset*CellWidth[jet_dim][0];
 
  printf("JetLeftCorner = %g %g %g\n", JetLeftCorner[0],JetLeftCorner[1],JetLeftCorner[2]);
  printf("JetRightCorner = %g %g %g\n", JetRightCorner[0],JetRightCorner[1],JetRightCorner[2]);
  /* Compute indices of jet launch region. */

  int JetStartIndex[MAX_DIMENSION], JetEndIndex[MAX_DIMENSION];

  for (dim = 0; dim < GridRank; dim++) {

    /* Compute start and end indices of jet */

    JetStartIndex[dim] = nint((JetLeftCorner[dim] - CellLeftEdge[dim][0] - 0.5*CellWidth[dim][0])/CellWidth[dim][0]);
    JetEndIndex[dim] = nint((JetRightCorner[dim] - CellLeftEdge[dim][0] - 0.5*CellWidth[dim][0])/CellWidth[dim][0]);

    /* If Jet is not on this grid, return. */

    if (JetStartIndex[dim] > GridDimension[dim]-1 || JetEndIndex[dim] < 0)
      return SUCCESS;

  } // end: loop over dim
    printf("JetStartIndex = %d %d %d\n", JetStartIndex[0],JetStartIndex[1],JetStartIndex[2]);
    printf("JetEndIndex = %d %d %d\n", JetEndIndex[0],JetEndIndex[1],JetEndIndex[2]);

  /* Compute mass and momentum to be put into cells in code units. */

  float JetNormalization = 0.0, density_normalization, radius, Tramp;
  for (j = JetStartIndex[1]; j <= JetEndIndex[1]; j++) {
    for (i = JetStartIndex[0]; i <= JetEndIndex[0]; i++) {
      radius = sqrt(pow((CellLeftEdge[0][0] + (i+0.5)*CellWidth[0][0] - JetCenter[0]), 2) +
                    pow((CellLeftEdge[1][0] + (j+0.5)*CellWidth[1][0] - JetCenter[1]), 2) )/
               CellWidth[0][0];
      JetNormalization += exp(-pow(radius/JetScaleRadius,2)/2.0);   //add 2!!!!!!!!!!!!!!, print stqtement
    }
  }
  /* Convert to code units. */
  JetMdot = (ClusterSMBHJetMdot*SolarMass/3.1557e7)/(MassUnits/TimeUnits);  // from M_sun/yr to code units
  density_normalization = (JetMdot/JetNormalization)*dtFixed/pow(CellWidth[0][0], 3);
  printf("ClusterSMBHTramp= %g\n", ClusterSMBHTramp);
  Tramp = ClusterSMBHTramp*1.0e6*3.1557e7/TimeUnits;  // from Myr to code units 
  /* If Time is earlier than ClusterSMBHStartTime, return. */
  printf("Time and ClusterSMBHStartTime= %g %g\n", Time, ClusterSMBHStartTime);
  if (Time-ClusterSMBHStartTime < 0.0)
    return SUCCESS;
  JetVelocity = ClusterSMBHJetVelocity*1.0e5/VelocityUnits; //from km/s to code units
  JetVelocity *= min((Time-ClusterSMBHStartTime)/Tramp, 1.0);     //linear ramp
  FastJetVelocity = ClusterSMBHFastJetVelocity*1.0e5/VelocityUnits; //from km/s to code units
  FastJetVelocity *= min((Time-ClusterSMBHStartTime)/Tramp, 1.0);     //linear ramp
//  JetVelocity *= 0.5*tanh(5.0*((Time-ClusterSMBHStartTime)/Tramp-0.5)+1.0);     // tanh ramp
  printf("density_normalization= %g\n", density_normalization);
  printf("JetVelocity= %g\n", JetVelocity);

  /* Clip edge of jet launching disk so we don't set cell off the edge of the grid. */


  for (dim = 0; dim < GridRank; dim++) {
    if (dim != jet_dim) {
      JetStartIndex[dim] = max(JetStartIndex[dim], 0);
      JetEndIndex[dim] = min(JetEndIndex[dim], GridDimension[dim]-1);
    }
  }

  /* Loop over launch disks and set cell values (this code assumes jet_dim = 2). */
  float density_ratio, density_add, xpos, ypos, JetVelocity_z, JetVelocity_xy;
  for (j = JetStartIndex[1]; j <= JetEndIndex[1]; j++) {
    for (i = JetStartIndex[0]; i <= JetEndIndex[0]; i++) {
      ///index = GRIDINDEX_NOGHOST(i,j,k);  ///replace GRIDINDEX(i,j,k)
        xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - JetCenter[0];  //in the cell surface center
        ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - JetCenter[1];  //not in cellwidth
	radius = sqrt(pow(xpos,2) + pow(ypos, 2))/CellWidth[0][0];  //in cell width
	density_add = density_normalization*exp(-pow(radius/JetScaleRadius,2)/2.0);
	JetVelocity = (radius > ClusterSMBHFastJetRadius) ? JetVelocity : FastJetVelocity;
printf("JetRadius, FastJetRadius and JetVelocity= %g %g \n", ClusterSMBHJetRadius, ClusterSMBHFastJetRadius, JetVelocity);
      if (ClusterSMBHJetAngleRadius < 0.1) {   // if jet openning angle = 0, set ClusterSMBHJetAngleRadius=0
	JetVelocity_z = JetVelocity;
	JetVelocity_xy = 0;
	}
      else {
	JetVelocity_z = JetVelocity * ClusterSMBHJetAngleRadius / sqrt(pow(ClusterSMBHJetAngleRadius, 2) + pow(radius, 2));
	JetVelocity_xy = JetVelocity * radius / sqrt(pow(ClusterSMBHJetAngleRadius, 2) + pow(radius, 2));
	}
	/*this is the bottom jet: */
      if (JetStartIndex[jet_dim] >= 0) {   
        k = JetStartIndex[jet_dim];  //start from the lower(outer) boundary of the cell
	k1 = (HydroMethod == Zeus_Hydro) ? (k+1):k
	BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += density_add;
	density_ratio = density_add/ BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
    printf("density_add and density_ratio upper jet= %g %g \n", density_add, density_ratio);
	BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k1)] = density_ratio * JetVelocity_xy * (xpos/CellWidth[0][0]) / radius + (1.0-density_ratio)*BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k1)];
	BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k1)] = density_ratio * JetVelocity_xy * (ypos/CellWidth[0][0]) / radius + (1.0-density_ratio)*BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k1)];
	BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k1)] = - density_ratio*JetVelocity_z + (1.0-density_ratio)*BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k1)];
	//	BaryonField[GENum][GRIDINDEX_NOGHOST(i,j,k)] += XXX;
printf("lower jet BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)] = %g \n", BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)]);
      }
	/*this is the top jet: */
      if (JetEndIndex[jet_dim] <= GridDimension[jet_dim]-1) { 
        k = JetEndIndex[jet_dim];
        BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += density_add;
        density_ratio = density_add/ BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
    printf("density_add and density_ratio lower jet= %g %g\n", density_add, density_ratio);
        BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_xy * (xpos/CellWidth[0][0]) / radius + (1.0-density_ratio)*BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_xy * (ypos/CellWidth[0][0]) / radius + (1.0-density_ratio)*BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)];
	BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio*JetVelocity_z + (1.0-density_ratio)*BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)];
        //      BaryonField[GENum][GRIDINDEX_NOGHOST(i,j,k)] += XXX;
printf("upper jet BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)] = %g \n", BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)]);
      }
    }
  }

  /* loop over cells to be modified, add jet mass, momentum, and energy. */

  return SUCCESS;

}
 

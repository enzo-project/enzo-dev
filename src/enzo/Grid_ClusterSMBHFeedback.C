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
extern float ClusterSMBHColdGasMass;

int grid::ClusterSMBHFeedback(int level)
{
  /* Only use Zeus, not PPM*/
  if (HydroMethod != Zeus_Hydro) 
    ENZO_FAIL("Error in HydroMethod. Please use Zeus.");
  

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Return if not on most-refined level. */   
  if (ClusterSMBHCalculateGasMass != 4 && level != MaximumRefinementLevel)
    return SUCCESS;   /*jet is not on the most-refined level only for Bondi*/

  /* Return if using method 1 or 2 and Switch is off. */
  if (ClusterSMBHCalculateGasMass != 0 && ClusterSMBHFeedbackSwitch == FALSE)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                                             Vel3Num, TENum) == FAIL)   ///this or thisgrid
     ENZO_FAIL("Error in IdentifyPhysicalQuantities.");


  /* Compute the jet launching region and the disk region
     (assume jet launched from PointSourceGravityPosition) */

  FLOAT JetLeftCorner[MAX_DIMENSION], JetRightCorner[MAX_DIMENSION];
  FLOAT JetCenter[MAX_DIMENSION];

  FLOAT DiskLeftCorner[MAX_DIMENSION], DiskRightCorner[MAX_DIMENSION];
  FLOAT DiskCenter[MAX_DIMENSION];

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
    TimeUnits = 1.0, VelocityUnits = 1.0;
  double MassUnits = 1.0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }
  MassUnits = DensityUnits*POW(LengthUnits,3);
  /* If Time is earlier than ClusterSMBHStartTime, return. */
  if (Time-ClusterSMBHStartTime < 0.0)
    return SUCCESS;

  int i, j, k, dim = 0;
  int jet_dim;  // z-axis (should make parameter?)
  jet_dim = ClusterSMBHJetDim % 3;
  float JetScaleRadius; // cellwidths
  float JetMdot; // Jet mass flow in SolarMass/year (need to convert units)-- gets value from parameter ClusterSMBHJetMdot
  float JetVelocity, FastJetVelocity; // Jet Velocity in km/s (should make parameter)-- gets value from parameter ClusterSMBHJetVelocity
  JetScaleRadius = ClusterSMBHJetRadius/2.0;  //JetScaleRadius is half the radius of the jet launch region in cellwidths
  float DiskRadius; //ClusterSMBHDiskRadius = 0.5;  //ClusterSMBHDiskRadiu now a parameter
  DiskRadius = ClusterSMBHDiskRadius*kpc/LengthUnits; //from kpc to codeunits 

  for (dim = 0; dim < GridRank; dim++) {
    JetCenter[dim] = PointSourceGravityPosition[dim];
    JetLeftCorner[dim] = JetCenter[dim];
    JetRightCorner[dim] = JetCenter[dim];
    if (dim != jet_dim) {
      JetLeftCorner[dim] -= ClusterSMBHJetRadius*CellWidth[dim][0];
      JetRightCorner[dim] += ClusterSMBHJetRadius*CellWidth[dim][0];
    }
    DiskCenter[dim] = PointSourceGravityPosition[dim];
    DiskLeftCorner[dim] = PointSourceGravityPosition[dim]- DiskRadius;
    DiskRightCorner[dim] = PointSourceGravityPosition[dim] + DiskRadius;
  }

  JetLeftCorner[jet_dim] -= ClusterSMBHJetLaunchOffset*kpc/LengthUnits;
  JetRightCorner[jet_dim] += ClusterSMBHJetLaunchOffset*kpc/LengthUnits; //from kpc to codeunits


  /* Compute indices of jet launch region. */

  int JetStartIndex[MAX_DIMENSION], JetEndIndex[MAX_DIMENSION];
  bool JetOnGrid = true, DiskOnGrid = true;
  int DiskStartIndex[MAX_DIMENSION], DiskEndIndex[MAX_DIMENSION];

  for (dim = 0; dim < GridRank; dim++) {

    /* Compute start and end indices of jet */

    JetStartIndex[dim] = nint((JetLeftCorner[dim] - CellLeftEdge[dim][0] - 0.5*CellWidth[dim][0])/CellWidth[dim][0]);
    JetEndIndex[dim] = nint((JetRightCorner[dim] - CellLeftEdge[dim][0] - 0.5*CellWidth[dim][0])/CellWidth[dim][0]);

    /* If Jet is not on this grid, return. */

    if (JetStartIndex[dim] > GridDimension[dim]-1 || JetEndIndex[dim] < 0)
      JetOnGrid = false;

    /*For Bondi*/
    if (ClusterSMBHCalculateGasMass == 4 && level != MultiRefineRegionMaximumOuterLevel)
      JetOnGrid = false;
 
    /*When not Bondi*/
    if (ClusterSMBHCalculateGasMass != 4 && level != MaximumRefinementLevel)
      JetOnGrid = false;

    DiskStartIndex[dim] = nint((DiskLeftCorner[dim] - CellLeftEdge[dim][0] - 0.5*CellWidth[dim][0])/CellWidth
[dim][0]);
    DiskEndIndex[dim] = nint((DiskRightCorner[dim] - CellLeftEdge[dim][0] - 0.5*CellWidth[dim][0])/CellWidth[
dim][0]);
    DiskStartIndex[dim] = max(DiskStartIndex[dim], GridStartIndex[dim]);
    DiskEndIndex[dim] = min(DiskEndIndex[dim], GridEndIndex[dim]);

    /* If Disk is not on this grid, return. */
    if (DiskStartIndex[dim] > GridEndIndex[dim] || DiskEndIndex[dim] < GridStartIndex[dim])
      DiskOnGrid = false;
  } // end: loop over dim

  /* Compute mass and momentum to be put into cells in code units if Jet is on this grid. */
  JetMdot = (ClusterSMBHJetMdot*SolarMass/3.1557e7)/(MassUnits/TimeUnits);  // from M_sun/yr to code units
if (JetOnGrid == true){
  float JetNormalization = 0.0, density_normalization, radius, Tramp;
  if (jet_dim == 2){
    for (j = JetStartIndex[1]; j <= JetEndIndex[1]; j++) {
      for (i = JetStartIndex[0]; i <= JetEndIndex[0]; i++) {
        radius = sqrt(POW((CellLeftEdge[0][0] + (i+0.5)*CellWidth[0][0] - JetCenter[0]), 2) +
                      POW((CellLeftEdge[1][0] + (j+0.5)*CellWidth[1][0] - JetCenter[1]), 2) )/
                 CellWidth[0][0];
        JetNormalization += exp(-POW(radius/JetScaleRadius,2)/2.0);   //add 2!!!!!!!!!!!!!!, print stqtement
      }
    }
  }
  else if(jet_dim == 0){
    for (j = JetStartIndex[1]; j <= JetEndIndex[1]; j++) {
      for (k = JetStartIndex[2]; k <= JetEndIndex[2]; k++) {
        radius = sqrt(POW((CellLeftEdge[2][0] + (k+0.5)*CellWidth[2][0] - JetCenter[2]), 2) +
                      POW((CellLeftEdge[1][0] + (j+0.5)*CellWidth[1][0] - JetCenter[1]), 2) )/
                 CellWidth[0][0];
        JetNormalization += exp(-POW(radius/JetScaleRadius,2)/2.0);   //add 2!!!!!!!!!!!!!!, print stqtement
      }
    }
  }
  if (jet_dim == 1){
    for (k = JetStartIndex[2]; k <= JetEndIndex[2]; k++) {
      for (i = JetStartIndex[0]; i <= JetEndIndex[0]; i++) {
        radius = sqrt(POW((CellLeftEdge[0][0] + (i+0.5)*CellWidth[0][0] - JetCenter[0]), 2) +
                      POW((CellLeftEdge[2][0] + (k+0.5)*CellWidth[2][0] - JetCenter[2]), 2) )/
                 CellWidth[0][0];
        JetNormalization += exp(-POW(radius/JetScaleRadius,2)/2.0);   //add 2!!!!!!!!!!!!!!, print stqtement
      }
    }
  }

  /* Convert to code units. */
  density_normalization = (JetMdot/JetNormalization)*dtFixed/POW(CellWidth[0][0], 3);
  Tramp = ClusterSMBHTramp*1.0e6*3.1557e7/TimeUnits;  // from Myr to code units 

  JetVelocity = sqrt((ClusterSMBHJetEdot*1.0e44*ClusterSMBHKineticFraction*2)/(ClusterSMBHJetMdot*SolarMass/3.1557e7))/VelocityUnits;
  JetVelocity *= min((Time-ClusterSMBHStartTime)/Tramp, 1.0);     //linear ramp
  
  /* Clip edge of jet launching disk so we don't set cell off the edge of the grid. */


  for (dim = 0; dim < GridRank; dim++) {
    if (dim != jet_dim) {
      JetStartIndex[dim] = max(JetStartIndex[dim], 0);
      JetEndIndex[dim] = min(JetEndIndex[dim], GridDimension[dim]-1);
    }
  }

  float density_ratio, density_add, energy_add, xpos, ypos, zpos, JetVelocity_z, JetVelocity_zy, JetVelocity_x, JetVelocity_y, JetVelocity_xy, JetVelocity_zx;

if (jet_dim == 2){
  /* Loop over launch disks and set cell values (this code assumes jet_dim = 2). */
  /* loop over cells to be modified, add jet mass, momentum, and energy. */
  for (j = JetStartIndex[1]; j <= JetEndIndex[1]; j++) {
    for (i = JetStartIndex[0]; i <= JetEndIndex[0]; i++) {
      xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - JetCenter[0];  //in the cell surface center
      ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - JetCenter[1];  //not in cellwidth
      radius = sqrt(POW(xpos,2) + POW(ypos, 2))/CellWidth[0][0];  //in cell width
      density_add = density_normalization*exp(-POW(radius/JetScaleRadius,2)/2.0);
      energy_add = ((1.0-ClusterSMBHKineticFraction)/ClusterSMBHKineticFraction)*0.5*density_add*POW(JetVelocity, 2.0);
      //JetVelocity = (radius > ClusterSMBHFastJetRadius) ? SlowJetVelocity : FastJetVelocity;
      if (ClusterSMBHJetOpenAngleRadius < 0.00001) {   // if jet openning angle = 0, set ClusterSMBHJetOpenAngleRadius=0
	JetVelocity_z = JetVelocity*cos(ClusterSMBHJetAngleTheta*pi);
	JetVelocity_x = JetVelocity;  // mutiplied by sincos later
	JetVelocity_y = JetVelocity;  // mutiplied by sincos later
        if (ClusterSMBHJetPrecessionPeriod > 0.00001 || ClusterSMBHJetPrecessionPeriod < -0.00001)
           ClusterSMBHJetAnglePhi = Time*2.0/(ClusterSMBHJetPrecessionPeriod*1.0e6*3.1557e7/TimeUnits);  // ClusterSMBHJetPrecessionPeriod from Myr to codeunit; *2.0 instead of 2*pi because pi is used later
      }
      else {
	JetVelocity_z = JetVelocity * ClusterSMBHJetOpenAngleRadius / sqrt(POW(ClusterSMBHJetOpenAngleRadius, 2) + POW(radius, 2));
	JetVelocity_xy = JetVelocity * radius / sqrt(POW(ClusterSMBHJetOpenAngleRadius, 2) + POW(radius, 2));
        JetVelocity_x = JetVelocity_xy * (xpos/CellWidth[0][0]) / radius;
        JetVelocity_y = JetVelocity_xy * (ypos/CellWidth[0][0]) / radius; 
      }

      /*this is the bottom jet: */

      if (JetStartIndex[jet_dim] >= 0) {   
        k = JetStartIndex[jet_dim];  //start from the lower(outer) boundary of the cell
        BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += density_add;
        density_ratio = density_add/ BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
        //printf("density_add and density_ratio upper jet= %g %g \n", density_add, density_ratio);
        BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_x*sin(ClusterSMBHJetAngleTheta*pi)*cos(ClusterSMBHJetAnglePhi*pi+pi) + (1.0-density_ratio)*BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_y*sin(ClusterSMBHJetAngleTheta*pi)*sin(ClusterSMBHJetAnglePhi*pi+pi) + (1.0-density_ratio)*BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)];
	BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] = BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)]*(1.0-density_ratio) + energy_add/BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];

       /* If Using Zeus, shift the index for z-velocity */

      if (HydroMethod == Zeus_Hydro) { 
 	if (k+1 <= GridDimension[jet_dim]-1) { // update velocity if it is still on the grid
          BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k+1)] = - density_ratio*JetVelocity_z + (1.0-density_ratio)*BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k+1)];
	}   
      }   //end Zeus
      else {
	BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)] = - density_ratio*JetVelocity_z + (1.0-density_ratio)*BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)];
        }
      }  //end bottom jet

	/*this is the top jet: */
      if (JetEndIndex[jet_dim] <= GridDimension[jet_dim]-1) { 
        k = JetEndIndex[jet_dim];
        BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += density_add;
        density_ratio = density_add/ BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_x*sin(ClusterSMBHJetAngleTheta*pi)*cos(ClusterSMBHJetAnglePhi*pi) + (1.0-density_ratio)*BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_y*sin(ClusterSMBHJetAngleTheta*pi)*sin(ClusterSMBHJetAnglePhi*pi) + (1.0-density_ratio)*BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)];
	BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio*JetVelocity_z + (1.0-density_ratio)*BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)];
	BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] = BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)]*(1.0-density_ratio) + energy_add/BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
      } //end top jet
    }
  }
} //end if jet_dim == 2

if (jet_dim == 0){
  /* Loop over launch disks and set cell values (this code assumes jet_dim = 0). */
  /* loop over cells to be modified, add jet mass, momentum, and energy. */
  for (j = JetStartIndex[1]; j <= JetEndIndex[1]; j++) {
    for (k = JetStartIndex[2]; k <= JetEndIndex[2]; k++) {
      zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - JetCenter[2];  //in the cell surface center
      ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - JetCenter[1];  //not in cellwidth
      radius = sqrt(POW(zpos,2) + POW(ypos, 2))/CellWidth[0][0];  //in cell width
      density_add = density_normalization*exp(-POW(radius/JetScaleRadius,2)/2.0);
      energy_add = ((1.0-ClusterSMBHKineticFraction)/ClusterSMBHKineticFraction)*0.5*density_add*POW(JetVelocity, 2.0);
      if (ClusterSMBHJetOpenAngleRadius < 0.00001) {   // if jet openning angle = 0, set ClusterSMBHJetOpenAngleRadius=0
        JetVelocity_x = JetVelocity*cos(ClusterSMBHJetAngleTheta*pi);
        JetVelocity_z = JetVelocity;  // mutiplied by sincos later
        JetVelocity_y = JetVelocity;  // mutiplied by sincos later
        if (ClusterSMBHJetPrecessionPeriod > 0.00001)
           ClusterSMBHJetAnglePhi = Time*2.0/(ClusterSMBHJetPrecessionPeriod*1.0e6*3.1557e7/TimeUnits);  // ClusterSMBHJetPrecessionPeriod from Myr to codeunit; *2.0 instead of 2*pi because pi is used later
      }
      else {
        JetVelocity_x = JetVelocity * ClusterSMBHJetOpenAngleRadius / sqrt(POW(ClusterSMBHJetOpenAngleRadius, 2) + POW(radius, 2));
        JetVelocity_zy = JetVelocity * radius / sqrt(POW(ClusterSMBHJetOpenAngleRadius, 2) + POW(radius, 2));
        JetVelocity_z = JetVelocity_zy * (zpos/CellWidth[0][0]) / radius;
        JetVelocity_y = JetVelocity_zy * (ypos/CellWidth[0][0]) / radius;
      }

      /*this is the bottom jet: */

      if (JetStartIndex[jet_dim] >= 0) {
        i = JetStartIndex[jet_dim];  //start from the lower(outer) boundary of the cell
        BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += density_add;
        density_ratio = density_add/ BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_z*sin(ClusterSMBHJetAngleTheta*pi)*cos(ClusterSMBHJetAnglePhi*pi+pi) + (1.0-density_ratio)*BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_y*sin(ClusterSMBHJetAngleTheta*pi)*sin(ClusterSMBHJetAnglePhi*pi+pi) + (1.0-density_ratio)*BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] = BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)]*(1.0-density_ratio) + energy_add/BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];

       /* If Using Zeus, shift the index for z-velocity */

      if (HydroMethod == Zeus_Hydro) {
        if (i+1 <= GridDimension[jet_dim]-1) { // update velocity if it is still on the grid
          BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i+1,j,k)] = - density_ratio*JetVelocity_x + (1.0-density_ratio)*BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i+1,j,k)];
        }
      }   //end Zeus
      else {
        BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)] = - density_ratio*JetVelocity_x + (1.0-density_ratio)*BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)];
        }
      }  //end bottom jet

        /*this is the top jet: */
      if (JetEndIndex[jet_dim] <= GridDimension[jet_dim]-1) {
        i = JetEndIndex[jet_dim];
        BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += density_add;
        density_ratio = density_add/ BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_z*sin(ClusterSMBHJetAngleTheta*pi)*cos(ClusterSMBHJetAnglePhi*pi) + (1.0-density_ratio)*BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_y*sin(ClusterSMBHJetAngleTheta*pi)*sin(ClusterSMBHJetAnglePhi*pi) + (1.0-density_ratio)*BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio*JetVelocity_x + (1.0-density_ratio)*BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] = BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)]*(1.0-density_ratio) +
energy_add/BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
      } //end top jet
    }
  }
} //end if jet_dim == 0


if (jet_dim == 1){
  /* Loop over launch disks and set cell values (this code assumes jet_dim = 0). */
  /* loop over cells to be modified, add jet mass, momentum, and energy. */
  for (i = JetStartIndex[0]; i <= JetEndIndex[0]; i++) {
    for (k = JetStartIndex[2]; k <= JetEndIndex[2]; k++) {
      zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - JetCenter[2];  //in the cell surface center
      xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - JetCenter[0];  //not in cellwidth
      radius = sqrt(POW(zpos,2) + POW(xpos, 2))/CellWidth[0][0];  //in cell width
      density_add = density_normalization*exp(-POW(radius/JetScaleRadius,2)/2.0);
      energy_add = ((1.0-ClusterSMBHKineticFraction)/ClusterSMBHKineticFraction)*0.5*density_add*POW(JetVelocity, 2.0);
      if (ClusterSMBHJetOpenAngleRadius < 0.00001) {   // if jet openning angle = 0, set ClusterSMBHJetOpenAngleRadius=0
        JetVelocity_y = JetVelocity*cos(ClusterSMBHJetAngleTheta*pi);
        JetVelocity_z = JetVelocity;  // mutiplied by sincos later
        JetVelocity_x = JetVelocity;  // mutiplied by sincos later
        if (ClusterSMBHJetPrecessionPeriod > 0.00001)
           ClusterSMBHJetAnglePhi = Time*2.0/(ClusterSMBHJetPrecessionPeriod*1.0e6*3.1557e7/TimeUnits);  // ClusterSMBHJetPrecessionPeriod from Myr to codeunit; *2.0 instead of 2*pi because pi is used later
      }
      else {
        JetVelocity_y = JetVelocity * ClusterSMBHJetOpenAngleRadius / sqrt(POW(ClusterSMBHJetOpenAngleRadius, 2) + POW(radius, 2));
        JetVelocity_zx = JetVelocity * radius / sqrt(POW(ClusterSMBHJetOpenAngleRadius, 2) + POW(radius, 2));
        JetVelocity_z = JetVelocity_zx * (zpos/CellWidth[0][0]) / radius;
        JetVelocity_x = JetVelocity_zx * (ypos/CellWidth[0][0]) / radius;
      }

      /*this is the bottom jet: */

      if (JetStartIndex[jet_dim] >= 0) { 
        j = JetStartIndex[jet_dim];  //start from the lower(outer) boundary of the cell
        BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += density_add;
        density_ratio = density_add/ BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_z*sin(ClusterSMBHJetAngleTheta*pi)*cos(ClusterSMBHJetAnglePhi*pi+pi) + (1.0-density_ratio)*BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_x*sin(ClusterSMBHJetAngleTheta*pi)*sin(ClusterSMBHJetAnglePhi*pi+pi) + (1.0-density_ratio)*BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] = BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)]*(1.0-density_ratio) + energy_add/BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];

       /* If Using Zeus, shift the index for z-velocity */

      if (HydroMethod == Zeus_Hydro) { 
        if (j+1 <= GridDimension[jet_dim]-1) { // update velocity if it is still on the grid
          BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j+1,k)] = - density_ratio*JetVelocity_y + (1.0-density_ratio)*BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j+1,k)];
        }
      }   //end Zeus
      else {
        BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)] = - density_ratio*JetVelocity_y + (1.0-density_ratio)*BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)];
        }
      }  //end bottom jet

        /*this is the top jet: */
      if (JetEndIndex[jet_dim] <= GridDimension[jet_dim]-1) {
        j = JetEndIndex[jet_dim];
        BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += density_add;
        density_ratio = density_add/ BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_z*sin(ClusterSMBHJetAngleTheta*pi)*cos(ClusterSMBHJetAnglePhi*pi) + (1.0-density_ratio)*BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio * JetVelocity_x*sin(ClusterSMBHJetAngleTheta*pi)*sin(ClusterSMBHJetAnglePhi*pi) + (1.0-density_ratio)*BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)] = density_ratio*JetVelocity_y + (1.0-density_ratio)*BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)];
        BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] = BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)]*(1.0-density_ratio) +
energy_add/BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)];
      } //end top jet
    }
  }
} //end if jet_dim == 1
} // end JetOnGrid==true


  /* loop over cells of disk, remove mass. */
  /* Return if not on most-refined level. */
  if (level != MaximumRefinementLevel)
    return SUCCESS;

if (DiskOnGrid == true && ClusterSMBHCalculateGasMass != 0){
  float AccretionRate = JetMdot*2.0; // in codeunit  *2 because Mdot is Mdot of one jet. There are two jets!
  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  float ColdGasTemperature = 3.0e4;       //in K--parameter?
  if (ClusterSMBHCalculateGasMass == 4){
  ColdGasTemperature = 3.0e8;       //basically whatever--everything gets accreted
  }
  float *BaryonFieldTemperature = new float[size];  // i.e. temperature
  if (BaryonFieldTemperature == NULL)
    ENZO_FAIL("Unable to allocate Temperature field in Grid_ClusterSMBHEachGridGasMass.");
  this->ComputeTemperatureField(BaryonFieldTemperature);
  for (k = DiskStartIndex[2]; k <= DiskEndIndex[2]; k++) {
    for (j = DiskStartIndex[1]; j <= DiskEndIndex[1]; j++) {
      for (i = DiskStartIndex[0]; i <= DiskEndIndex[0]; i++) {
        if (BaryonFieldTemperature[GRIDINDEX_NOGHOST(i,j,k)] < ColdGasTemperature)
          {if (ClusterSMBHCalculateGasMass == 4){
             BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)]=1.0e-35;
             }
             else{
                BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] *= 1.0 - AccretionRate*dtFixed/ClusterSMBHColdGasMass; //take out part of the mass
             }
          }
      }
    }
  }
  delete [] BaryonFieldTemperature;
}

  return SUCCESS;

}
 

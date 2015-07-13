/***********************************************************************
/
/  GRID: ADD RADIO-MODE JET-LIKE FEEDBACK BASED ON STATIC SMBH
/
/  written by: Yuan Li and Greg Bryan
/  date:       May, 2012
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

//define global variable here
extern float ClusterSMBHColdGasMass;

int grid::ClusterSMBHEachGridGasMass(int level)
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


  /* Compute the cold disk gas region
     (assume the center of the disk is PointSourceGravityPosition) */

  FLOAT DiskLeftCorner[MAX_DIMENSION], DiskRightCorner[MAX_DIMENSION];
  FLOAT DiskCenter[MAX_DIMENSION];

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
    TimeUnits = 1.0, VelocityUnits = 1.0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  int dim = 0;
  float DiskRadius, ClusterSMBHDiskRadius = 0.5;  //ClusterSMBHDiskRadiu make parameter?
  DiskRadius = ClusterSMBHDiskRadius*kpc/LengthUnits; //from kpc to codeunits 
  for (dim = 0; dim < GridRank; dim++) {
    DiskCenter[dim] = PointSourceGravityPosition[dim];
    DiskLeftCorner[dim] = PointSourceGravityPosition[dim]- DiskRadius;
    DiskRightCorner[dim] = PointSourceGravityPosition[dim] + DiskRadius;
  }
//  printf("DiskLeftCorner = %g %g %g\n", DiskLeftCorner[0],DiskLeftCorner[1],DiskLeftCorner[2]);
//  printf("DiskRightCorner = %g %g %g\n", DiskRightCorner[0],DiskRightCorner[1],DiskRightCorner[2]);

  /* Compute indices of disk region. */

  int DiskStartIndex[MAX_DIMENSION], DiskEndIndex[MAX_DIMENSION];

  for (dim = 0; dim < GridRank; dim++) {

    /* Compute start and end indices of jet */

    DiskStartIndex[dim] = nint((DiskLeftCorner[dim] - CellLeftEdge[dim][0] - 0.5*CellWidth[dim][0])/CellWidth[dim][0]);
    DiskEndIndex[dim] = nint((DiskRightCorner[dim] - CellLeftEdge[dim][0] - 0.5*CellWidth[dim][0])/CellWidth[dim][0]);
    DiskStartIndex[dim] = max(DiskStartIndex[dim], GridStartIndex[dim]);
    DiskEndIndex[dim] = min(DiskEndIndex[dim], GridEndIndex[dim]);
    /* If Disk is not on this grid, return. */

    if (DiskStartIndex[dim] > GridEndIndex[dim] || DiskEndIndex[dim] < GridStartIndex[dim])
      return SUCCESS;

  } // end: loop over dim

//    printf("DiskStartIndex = %d %d %d\n", DiskStartIndex[0],DiskStartIndex[1],DiskStartIndex[2]);
//    printf("DiskEndIndex = %d %d %d\n", DiskEndIndex[0],DiskEndIndex[1],DiskEndIndex[2]);

  int i, j, k, size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  float ColdGasTemperature = 3.0e4;       //in K--parameter?
  float *BaryonFieldTemperature = new float[size];  // i.e. temperature
  if (BaryonFieldTemperature == NULL)
    ENZO_FAIL("Unable to allocate Temperature field in Grid_ClusterSMBHEachGridGasMass.");
  this->ComputeTemperatureField(BaryonFieldTemperature);
  for (k = DiskStartIndex[2]; k <= DiskEndIndex[2]; k++) {
    for (j = DiskStartIndex[1]; j <= DiskEndIndex[1]; j++) {
      for (i = DiskStartIndex[0]; i <= DiskEndIndex[0]; i++) {
//        printf("BaryonFieldTemperature[GRIDINDEX_NOGHOST(i,j,k) = %g \n", BaryonFieldTemperature[GRIDINDEX_NOGHOST(i,j,k)]);
        if (BaryonFieldTemperature[GRIDINDEX_NOGHOST(i,j,k)] < ColdGasTemperature)
          ClusterSMBHColdGasMass += BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)]*pow(CellWidth[0][0],3);   //Assuming it is refined to the highest refinement level (otherwise we should use the CellWidth at the exact position.)
//          printf("BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] and ClusterSMBHColdGasMass = %g %g \n", BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)], ClusterSMBHColdGasMass);
//take out part of the mass in ClusterSMBHFeedback?
      }
    }
  }
// printf("Each Grid ClusterSMBHColdGasMass = %g \n", ClusterSMBHColdGasMass);
  delete [] BaryonFieldTemperature;
  return SUCCESS;

}
 

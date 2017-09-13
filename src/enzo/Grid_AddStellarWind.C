/***********************************************************************
/
/  GRID: ADD Stellar Wind from a fixed point (a single star)
/
/  written by: Yuan Li
/  date:       Aug, 2017
/  modified1: 
/
/  PURPOSE: Simulate AGB Wind from a star like Mira
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
             float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::AddStellarWind()
{
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                                             Vel3Num, TENum) == FAIL)   ///this or thisgrid
     ENZO_FAIL("Error in IdentifyPhysicalQuantities.");


  int i, j, k;
  Mu = 0.6;

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
    TimeUnits = 1.0, VelocityUnits = 1.0;
  double MassUnits=1.0 ;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }



  FLOAT r, x, y, z = 0;

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++) {

        /* Compute position */

        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        if (GridRank > 1)
          y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        if (GridRank > 2)
          z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

          /* Find distance from center of star (0.125, 0.5, 0.5). */

        FLOAT xpos, ypos, zpos;
        xpos = x-StellarWindCenterPosition[0];
        ypos = y-StellarWindCenterPosition[1];
        zpos = z-StellarWindCenterPosition[2];
        r = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
        r = max(r, 0.1*CellWidth[0][0]);
        if (r<StellarWindRadius){
          BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] = StellarWindDensity*POW(r/StellarWindRadius, -2);
          BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)] = StellarWindSpeed * (xpos/r)/VelocityUnits;
          BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)] = StellarWindSpeed * (ypos/r)/VelocityUnits;
          BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)] = StellarWindSpeed * (zpos/r)/VelocityUnits;
          BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] = StellarWindTemperature/TemperatureUnits/ ((Gamma-1.0)*Mu);
          if (HydroMethod != Zeus_Hydro){
              BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] += 0.5*POW(BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)], 2.0);
            if(GridRank > 1)
              BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] += 0.5*POW(BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)], 2.0);
            if(GridRank > 2)
              BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] += 0.5*POW(BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)], 2.0);
          }
        }
  }
  return SUCCESS;

}
 

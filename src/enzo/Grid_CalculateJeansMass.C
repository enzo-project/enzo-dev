/***********************************************************************
/
/  Calculate the Jeans Mass on the grid
/
/  written by: John Regan
/  date:       January 2017
/  modified1:
/
/  PURPOSE: Calculate the Jeans Mass
/
************************************************************************/
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "phys_constants.h"

#define TAKEAVERAGE 1
float grid::CalculateJeansMass(int DensNum, float *T, float DensityUnits)
{
  int i = 0, j = 0, k = 0;
  int index = 0, count = 0, maxindex = 0;
  float avgdensity = -1, JeansMass = 0.0;
  float *density = BaryonField[DensNum], avgtemp = -1, maxdensity = -1;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i, j, k);
	avgdensity += density[index];
	count++;
	avgtemp += T[index];
	if(density[index] > maxdensity) {
	  maxdensity = density[index];
	  maxindex = index;
	}
      }
    }
  }
#if TAKEAVERAGE
  avgdensity = avgdensity/(float)count;
  avgtemp = avgtemp/(float)count;
  avgdensity *= DensityUnits/mh;   //in cm^-3
  //printf("%s: avgdensity = %g\t avgtemp = %g\n", __FUNCTION__,
  //	avgdensity, avgtemp);
  JeansMass = 3e4*sqrt(POW(avgtemp, 3.0)/(avgdensity*1e6)); //In Msolar
#else
  JeansMass = 3e4*sqrt(POW(T[maxindex], 3.0)/(maxdensity*1e6)); //In Msolar  
#endif
  //printf("%s: JeansMass = %g\n", __FUNCTION__, JeansMass);
  return JeansMass; 
}

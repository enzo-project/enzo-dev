/***********************************************************************
/
/  CALCULATE TABULATED LW BACKGROUND
/
/  written by: Michele Trenti and Britton Smith
/  date:       April, 2009
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
  
int RadiationFieldLymanWernerTable(float Redshift,float *J21)
{

  int index;
  int midpoint,highpoint;
  float slope;

  // find redshift index
  if (Redshift > RadiationData.LymanWerner_redshift[1]) {
    index = 0;
  }
  else if (Redshift <= RadiationData.LymanWerner_redshift[RadiationData.NumberOfLWRedshiftBins-2]) {
    index = RadiationData.NumberOfLWRedshiftBins - 2;
  }
  else {

    // bisection
    index = 0;
    highpoint = RadiationData.NumberOfLWRedshiftBins - 1;
    while (highpoint - index > 1) {
      midpoint = (highpoint + index) >> 1;
      if (Redshift <= RadiationData.LymanWerner_redshift[midpoint]) index = midpoint;
      else highpoint = midpoint;
    }
  }

  slope = (RadiationData.LymanWerner_J21[index+1] - RadiationData.LymanWerner_J21[index]) /
    (RadiationData.LymanWerner_redshift[index+1] - RadiationData.LymanWerner_redshift[index]);

  *J21 = (Redshift - RadiationData.LymanWerner_redshift[index]) * slope +
    RadiationData.LymanWerner_J21[index];

//   fprintf(stderr,"Tabulated J21 redshift: %"ESYM" >= %"ESYM" >= %"ESYM".\n",
// 	  RadiationData.LymanWerner_redshift[index],Redshift,
// 	  RadiationData.LymanWerner_redshift[index+1]);
//   fprintf(stderr,"Tabulated J21: %"ESYM" <= %"ESYM" <= %"ESYM".\n",
// 	  RadiationData.LymanWerner_J21[index],*J21,
// 	  RadiationData.LymanWerner_J21[index+1]);

  return SUCCESS;
}

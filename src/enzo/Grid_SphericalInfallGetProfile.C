/***********************************************************************
/
/  GRID CLASS (FIND AND SAVE PROFILE FOR THE SPHERICAL INFALL TEST)
/
/  written by: Greg Bryan
/  date:       December, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "CosmologyParameters.h"
#include "Grid.h"
#include "SphericalInfall.h"
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
 
int grid::SphericalInfallGetProfile(int level, int ReportLevel)
{
 
  /* Return if this grid is not on this processor
      (this routine is not multi-processor safe). */
 
   if (MyProcessorNumber != ProcessorNumber)
     return SUCCESS;
 
  static int SphericalInfallGetProfileNumber = 0;
  static int SphericalInfallGetProfileLowestLevel = 0;
 
  /* ReportLevel == 0 means no reporting. */
 
  if (ReportLevel == 0 || NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* For ReportLevel > 1, ignore level 0. */
 
  if (level == 0 && ReportLevel > 1)
    return SUCCESS;
 
  /* declarations */
 
  int i, j, k, n, dim, Index[MAX_DIMENSION];
  char *SphericalInfallReportName, ProfileName[MAX_LINE_LENGTH];
  FILE *fptr;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
 
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
 
  /* Compute angular momentum. */
 
  float AngularMomentum[MAX_DIMENSION], MeanVelocity[MAX_DIMENSION],
        DMVelocity[MAX_DIMENSION];
  FLOAT CenterOfMass[MAX_DIMENSION], DMCofM[MAX_DIMENSION];
  this->CalculateAngularMomentum(SphericalInfallCenter, AngularMomentum,
			     MeanVelocity, DMVelocity, CenterOfMass, DMCofM);
  float LMod = sqrt(AngularMomentum[0]*AngularMomentum[0] +
		    AngularMomentum[1]*AngularMomentum[1] +
		    AngularMomentum[2]*AngularMomentum[2]  );
 
  /* Calculate Div-v. */
 
  int index, offset1 = GridDimension[0], offset2 = offset1*GridDimension[1];
  float DivVel, MaxDivVel = 0;
  for (k = GridStartIndex[2]+1; k <= GridEndIndex[2]-1; k++)
    for (j = GridStartIndex[1]+1; j <= GridEndIndex[1]-1; j++) {
      index = (j+k*GridDimension[1])*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]+1; i <= GridEndIndex[0]-1; i++) {
	DivVel = BaryonField[Vel1Num][index+1] - BaryonField[Vel1Num][index-1];
	if (GridRank > 1)
	  DivVel += BaryonField[Vel2Num][index+offset1] -
	            BaryonField[Vel2Num][index-offset1];
	if (GridRank > 2)
	  DivVel += BaryonField[Vel3Num][index+offset2] -
	            BaryonField[Vel3Num][index-offset2];
	MaxDivVel = max(MaxDivVel, fabs(DivVel));
      }
    }
  MaxDivVel /= 2.0*a*CellWidth[0][0];
 
  /* Compute max/min vel x,y,z. */
 
  float MinVel[MAX_DIMENSION], MaxVel[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    MinVel[dim] = MaxVel[dim] = 0;
    for (i = 0; i < GridDimension[0]*GridDimension[1]*GridDimension[2]; i++) {
      MinVel[dim] = min(MinVel[dim], BaryonField[Vel1Num+dim][i]);
      MaxVel[dim] = max(MaxVel[dim], BaryonField[Vel1Num+dim][i]);
    }
  }
  float MaxDens = 0, MinEntropy = huge_number, MaxDensVel = 0;
  for (i = 0; i < GridDimension[0]*GridDimension[1]*GridDimension[2]; i++) {
    MaxDens = max(MaxDens, BaryonField[DensNum][i]);
    if (MaxDens == BaryonField[DensNum][i])
      MaxDensVel = BaryonField[Vel1Num][i];
    MinEntropy = min(MinEntropy, BaryonField[GENum][i]/
                                 POW(BaryonField[DensNum][i], Gamma-1));
  }
 
  /* Open output file. */
 
  sprintf(ProfileName, "%s.L%1.1"ISYM".%4.4"ISYM, "SphericalInfallProfile", level,
	  SphericalInfallGetProfileNumber);
 
  if (level >= SphericalInfallGetProfileLowestLevel) {
    SphericalInfallGetProfileLowestLevel = level;
    SphericalInfallGetProfileNumber++;
  }
 
  char *mode;
  if (ReportLevel == 1) {
    mode = "a";
    SphericalInfallReportName = "SphericalInfallReport";
  }
  else {
    mode = "w";
    SphericalInfallReportName = ProfileName;
  }
 
  if ((fptr = fopen(SphericalInfallReportName, mode)) == FAIL) {
    ENZO_VFAIL("Error opening %s.\n", SphericalInfallReportName)
  }
 
  fprintf(fptr, "# l %"ISYM" t = %"FSYM" L = %"GSYM" %"GSYM" %"GSYM" %"GSYM" d(max) = %"GSYM" S(min) = %"GSYM" v_dmax = %"GSYM" div_v_max = %"GSYM"  dt/dt_divv = %"GSYM"\n",
	  level, Time,
	  AngularMomentum[0], AngularMomentum[1], AngularMomentum[2], LMod,
	  MaxDens, MinEntropy, MaxDensVel, MaxDivVel, MaxDivVel*dtFixed);
//	  MaxVel[0], MaxVel[1], MaxVel[2], MinVel[0], MinVel[1], MinVel[2]);
 
  /* Find center index. */
 
  for (dim = 0; dim < GridRank; dim++) {
    Index[dim] = int( (SphericalInfallCenter[dim] - *CellLeftEdge[dim]) /
                      *CellWidth[dim] );
    if (Index[dim] < 0 || Index[dim] >= GridDimension[dim]) {
      fclose(fptr);
      return SUCCESS;
    }
  }
 
  /* Loop from center to right edge of grid. */
 
  if (ReportLevel > 1) {

    int Offset = GridDimension[0]*(Index[1] + Index[2]*GridDimension[1]);
    for (i = Index[0]; i < GridDimension[0]; i++) {
      fprintf(fptr, "%"FSYM" ", CellLeftEdge[0][i] + 0.5*CellWidth[0][i] -
	                   SphericalInfallCenter[0]);
      for (n = 0; n < NumberOfBaryonFields; n++)
	fprintf(fptr, "%e ", BaryonField[n][i + Offset]);
      fprintf(fptr, "\n");
    }
  }
 
  /* Close file. */
 
  fclose(fptr);
 
  return SUCCESS;
}

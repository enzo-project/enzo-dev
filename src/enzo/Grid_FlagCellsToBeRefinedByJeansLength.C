/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE JEAN'S CRITERION)
/
/  written by: Greg Bryan
/  date:       February, 1998
/  modified1:  Alexei Kritsuk, Feb. 2004 mods for isothermal EOS.
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "phys_constants.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "hydro_rk/EOS.h"
 
/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
 
int grid::FlagCellsToBeRefinedByJeansLength()
{
  /* declarations */
 
  int i, dim;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Compute the temperature field. */
 
  float *temperature = NULL;
  if (ProblemType != 60 && ProblemType != 61 && (EOSType == 0)) { //AK
    temperature = new float[size];
    if (this->ComputeTemperatureField(temperature) == FAIL) {
      fprintf(stderr, "Error in grid->ComputeTemperature.\n");
      return -1;
    }
    /* This is less efficient, but it avoids too many conditionals */
    if(JeansRefinementColdTemperature > 0.0){
      for (i = 0; i < size; i++) temperature[i] =
        JeansRefinementColdTemperature;
    }
  }
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* Get density units. */
 
  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
 
  /* Compute constant for Jean's length computation.
      l_j = sqrt((Gamma*pi*k*T) / (G rho mu m_h))  . */
 
  FLOAT JLSquared = (double(Gamma*pi*kboltz/GravConst)/
		     (double(DensityUnits)*double(Mu)*double(mh))) /
                (double(LengthUnits)*double(LengthUnits));
 
  if (ProblemType == 60 || ProblemType == 61)
    JLSquared = double(4.0*3.14159*3.14159)/GravitationalConstant; //AK

  if (EOSType > 0)
    {
      float cs,dpdrho,dpde, eint, h, rho, p;
      EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 1) ;
      JLSquared = cs*cs*M_PI/GravConst/DensityUnits*VelocityUnits*VelocityUnits/LengthUnits/LengthUnits; // TA
    }

  /* This is the safety factor to decrease the Jean's length by. */
 
  JLSquared /= POW(RefineByJeansLengthSafetyFactor, 2);
 
/* printf("jl: JL, dx, t, d = %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", sqrt(JLSquared), CellWidth[0][3],
	 temperature[(3 + 3*GridDimension[1])*GridDimension[0]+3],
	 BaryonField[DensNum][(3 + 3*GridDimension[1])*GridDimension[0]+3]);*/
 
  /* Loop over grid. */
 
#ifdef UNUSED
  int j, k, index;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (j + k*GridDimension[1])*GridDimension[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	if ((CellWidth[0][i])*(CellWidth[0][i]) >
	    JLSquared*temperature[index+i]/BaryonField[DensNum][index+i])
	  FlaggingField[index+i]++;
    }
#endif /* UNUSED */
 
  FLOAT CellWidthSquared = CellWidth[0][0]*CellWidth[0][0];
  for (i = 0; i < size; i++)
    {
      if (EOSType == 0) {
	if (CellWidthSquared > JLSquared*temperature[i]/BaryonField[DensNum][i])
	  FlaggingField[i]++; 
      }
      else // isothermal and ploytropic sound speed version
	if (CellWidthSquared > JLSquared/BaryonField[DensNum][i])
	  FlaggingField[i]++; 
    }
 
  /* clean up */
 
  if (ProblemType != 60 && ProblemType != 61) //AK
    delete [] temperature;
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }
 
  return NumberOfFlaggedCells;
 
}

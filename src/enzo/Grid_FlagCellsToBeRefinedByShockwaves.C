/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINE BY Shockwaves)
/
/  written by: Samuel W. Skillman
/  date:       August 2008
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
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
#include "phys_constants.h"
 
int FindField(int field, int farray[], int numfields);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int grid::FlagCellsToBeRefinedByShockwaves(int level)
{
  /* declarations */
 
  int i, j, k, index, dim, NumberOfFlaggedCells = 0;
  float CellVolume;
 
  /* Return if this grid is not on this processor. */
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* error check */
   if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

   /* Find fields: density, total energy, velocity1-3. */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return -1;
  }


   /* Find Mach field.  If no Mach field exists, quit and yell at user. */
  int MachField = FALSE, MachNum;
  if ((MachNum = FindField(Mach, FieldType, NumberOfBaryonFields)) != -1){
    MachField = TRUE;
  } else{
    fprintf(stderr,"FlagCellsToBeRefinedByShockwaves: no Mach field!\n");
    return -1;
  } 

  if (FindShocksOnlyOnOutput == 1){
    fprintf(stderr, 
            "FlagCellsToBeRefinedByShockwaves: Refusing to refine when shocks \n"
            "are only found during output. Please change FindShocksOnlyOnOutput\n"
            "to 0 or 2.\n");
    ENZO_FAIL("Error in FlagCellsToBeRefinedByShockwaves.");
  } 

  /* Compute cell volume */
  
  CellVolume = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    CellVolume *= CellWidth[dim][0];
  
  /* compute size */ 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, MassUnits = 1, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }


  float *temperature = new float[size]; 
  if (this->ComputeTemperatureField(temperature) == FAIL){
    fprintf(stderr, "Error in grid->ComputeTemperatureField.\n");
    return FAIL;
  }
  float Csound;
  for(i=0; i<size; i++){
    Csound = sqrt(Gamma*kboltz*temperature[i]/(Mu*mh));
    if( (BaryonField[MachNum][i]*Csound >= ShockwaveRefinementMinVelocity) &&
	(level < ShockwaveRefinementMaxLevel) && 
	(BaryonField[MachNum][i] >= ShockwaveRefinementMinMach)){
      FlaggingField[i]++;

//     if( (BaryonField[MachNum][i] >= ShockwaveRefinementMinMach) &&
// 	(level < ShockwaveRefinementMaxLevel)){
//       FlaggingField[i]++;
      
    }
  }
  /* Count number of flagged Cells. */ 
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }
  return NumberOfFlaggedCells;
}

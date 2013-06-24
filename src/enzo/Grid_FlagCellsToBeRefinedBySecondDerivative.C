/***********************************************************************
  /
  /  GRID CLASS (FLAG CELLS TO BE REFINED BY SECOND DERIVATIVE)
  /
  /  written by: Samuel Skillman 
  /  date:       October, 2012
  /  modified1:  
  /
  /  PURPOSE: 
  /    flag cells based on the second derivative, normalized by an
  /    average of the first derivatives.
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

int grid::FlagCellsToBeRefinedBySecondDerivative()
{
  /* declarations */

  int i, j, k, index, dim;
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* error check */

  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  /* Make sure quantities are defined at least to dim 3 */

  for (dim = GridRank; dim < 3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }

  /* loop over all zones */
 
  for (dim = 0; dim < 3; dim++) {
    Start[dim] = GridStartIndex[dim];
    End[dim]   = GridEndIndex[dim];
  }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* allocate a temporary slope field. */
 
  float *TopBuffer = new float[size];
  float *BottomBuffer = new float[size];
 
  /* some problems do not need to check all the fields, in particular
     the velocity components */
 
  int NumberOfFields = NumberOfBaryonFields;

  // Override NumberOfFields for some specific problems:

  if (ProblemType ==  6) NumberOfFields = 1; // Implosion (AK)
  if (ProblemType ==  7) NumberOfFields = 2; // SedovBlast (AK)
  if (ProblemType == 11) NumberOfFields = 2; // RadiatingShock (BWO)
  if (ProblemType == 30) NumberOfFields = 2; // Cosmology (BWO 23 May 2006)

  bool doField=false;
  float MinimumSecondDerivativeForRefinementThis;
  int Offset = 1;
  int Offsets[3];
  for (dim=0; dim<GridRank; dim++)
    Offsets[dim] = 1;
  
  for (dim = 0; dim<GridRank-1; dim++){
    Offsets[dim+1] = Offsets[dim]*GridDimension[dim]; 
  }
  for (int field = 0; field < NumberOfFields; field++) {

    doField = false;
    if (SecondDerivativeFlaggingFields[0]==INT_UNDEFINED){ 
      MinimumSecondDerivativeForRefinementThis=MinimumSecondDerivativeForRefinement[0];
      doField=true;
    } else {
      for (int g=0; g<MAX_FLAGGING_METHODS; g++){
        if (SecondDerivativeFlaggingFields[g]==FieldType[field]){
          MinimumSecondDerivativeForRefinementThis=MinimumSecondDerivativeForRefinement[g];
          doField=true;
        }
      }
    }

    if (!doField)
      continue;

    /* loop over active dimensions */
    for (i = 0; i < size; i++){
      TopBuffer[i] = 0.0;
      BottomBuffer[i] = 1.0; //Protects from 0.0/0.0 later on.
    }

    /* compute second derivative criteria */
    for (k = Start[2]; k <= End[2]; k++)
      for (j = Start[1]; j <= End[1]; j++)
        for (i = Start[0]; i <= End[0]; i++) {
          index = i + j*GridDimension[0] +
              k*GridDimension[1]*GridDimension[0];
          BottomBuffer[index] = 0.0;
          for (int dimk = 0; dimk < GridRank; dimk++){
            for (int diml = 0; diml < GridRank; diml++){
              /* zero slope */
              TopBuffer[index] += 
                  0.125*POW(BaryonField[field][index + Offsets[dimk] + Offsets[diml]] -
                            BaryonField[field][index + Offsets[dimk] - Offsets[diml]] -
                            BaryonField[field][index - Offsets[dimk] + Offsets[diml]] +
                            BaryonField[field][index - Offsets[dimk] - Offsets[diml]], 2.0);
              BottomBuffer[index] +=
                  POW(0.5*(fabs(BaryonField[field][index + Offsets[dimk]] -
                                BaryonField[field][index]) +
                           fabs(BaryonField[field][index] -
                                BaryonField[field][index - Offsets[dimk]])) +
                      SecondDerivativeEpsilon * (fabs(BaryonField[field][index + Offsets[dimk] + Offsets[diml]]) +
                                        fabs(BaryonField[field][index + Offsets[dimk] - Offsets[diml]]) +
                                        fabs(BaryonField[field][index - Offsets[dimk] + Offsets[diml]]) +
                                        fabs(BaryonField[field][index - Offsets[dimk] - Offsets[dimk]])) , 2.0);
            }
          }
        }

    /* flag field based on second derivative */
    for (i = 0; i < size; i++){
      TopBuffer[i] = sqrt(TopBuffer[i]/BottomBuffer[i]);
      FlaggingField[i] += (TopBuffer[i] > MinimumSecondDerivativeForRefinementThis) ? 1 : 0;
    }
  }  // end loop over field


  /* delete buffer */

  delete [] TopBuffer, BottomBuffer;

  /* Count number of flagged Cells. */

  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i] > 0)
      NumberOfFlaggedCells++;

  return NumberOfFlaggedCells;

}

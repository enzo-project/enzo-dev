#include <stdlib.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "phys_constants.h"

int grid::CheckField(int FieldNum){
 /* -----------------------------------------
  * CheckField
  * -----------------------------------------
  * A. Emerick - Jun 2016
  *
  * Loops through get to check passed field for
  * negative values. Meant to be used as a
  * debugging tool.
  * -----------------------------------------*/

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  int size = 1;
  for (int dim = 0; dim < this->GridRank; dim++){
    size *= this->GridDimension[dim];
  }

  for(int i = 0; i < size; i ++){
    if (this->BaryonField[FieldNum][i] < 0.0){
        return FAIL;
    }
  }
  return SUCCESS;
}

int grid::CheckDensity(void){

  if (NumberOfBaryonFields == 0) return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  return (this->CheckField(DensNum));
}


int grid::CheckOTRadiation(void){

  if (NumberOfBaryonFields == 0) return SUCCESS;

  int OTLWkdissH2INum = FindField(OTLWkdissH2I, this->FieldType, this->NumberOfBaryonFields);
  int PeNum = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);

  int err = SUCCESS;
  if (PeNum > 0) err *= this->CheckField(PeNum);
  if (OTLWkdissH2INum > 0) err *= this->CheckField(OTLWkdissH2INum);

  return err;

}


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

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  return (this->CheckField(DensNum));
}

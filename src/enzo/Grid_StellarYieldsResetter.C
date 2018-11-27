#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"

/* function prototypes */
 

int grid::StellarYieldsResetter(int level)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum,GENum,Vel1Num,Vel2Num,
                                       Vel3Num,TENum) == FAIL){
    ENZO_FAIL("Error in IdentifyPhysicalQuantities");
  }

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (int i = 0; i < size; i ++){

    // identify field and zero
    for (int j = 0; j < StellarYieldsNumberOfSpecies; j ++){
      int anum = StellarYieldsResetAtomicNumbers[j];
      if (anum <= 0) break;

      int field_num;
      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, anum);
      BaryonField[field_num][i] = tiny_number * BaryonField[DensNum][i];
    }

  }

  return SUCCESS;
}

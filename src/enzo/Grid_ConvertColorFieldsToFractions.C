/***********************************************************************
/
/  COMMONLY-USED FRACTION-CONVERSION FUNCTION
/
/  written by: Matthew Turk
/  date:       May, 2011
/
/  PURPOSE:
/   This is a very tiny function.  It has its own source file so it can
/   call ENZO_FAIL.
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

void grid::ConvertColorFieldsToFractions() {
  /* Convert the species densities into fractional densities (i.e. divide
     by total baryonic density).  At the end we will multiply by the new
     density so that species fractions are maintained. */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num,H2INum, H2IINum;
  int i,j,k,field,index;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
        Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
  for (field = 0; field < NumberOfBaryonFields; field++)
    /* Note how fragile this color field detection mechanism is */
    if (FieldType[field] >= ElectronDensity && FieldType[field] <= ExtraType1)
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
        for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
          index = (k*GridDimension[1] + j)*GridDimension[0] +
            GridStartIndex[0];
          for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
            BaryonField[field][index] /= BaryonField[DensNum][index];
        }
}

void grid::ConvertColorFieldsFromFractions() {
  /* Convert the species densities into fractional densities (i.e. divide
     by total baryonic density).  At the end we will multiply by the new
     density so that species fractions are maintained. */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num,H2INum, H2IINum;
  int i,j,k,field,index;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
        Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
  for (field = 0; field < NumberOfBaryonFields; field++)
    /* Note how fragile this color field detection mechanism is */
    if (FieldType[field] >= ElectronDensity && FieldType[field] <= ExtraType1)
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
        for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
          index = (k*GridDimension[1] + j)*GridDimension[0] +
            GridStartIndex[0];
          for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
            BaryonField[field][index] *= BaryonField[DensNum][index];
        }
}

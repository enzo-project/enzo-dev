/*********************************************************************
/
/  CLEAR EMISSIVITY FUNCTION
/
/  written by: Geoffrey So
/  date:       Jan 9, 2009
/
/  PURPOSE:
/    To clear the emissivity array (implimented as a BaryonField) when
/    starting in the next timestep in the root grid, so that when
/    CalcEmiss(...) is called, BaryonField[EmisNum] will have all values 
/    initalized to 0
/
/
/
/
*********************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

int FindField(int field, int farray[], int numfields);

#ifdef EMISSIVITY
int grid::ClearEmissivity(){
  if( MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  int EtaNum;

  // Copying what's done in Grid_IdentifyPhysicalQuantities to check
  if ((EtaNum = FindField(Emissivity0, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find EmissivityField.\n");
    return FAIL;
  }
  for(int i=0; i<size; i++) {
    if (BaryonField[EtaNum][i] > 1e-40)
      printf("before clear was %22.16e\n", BaryonField[EtaNum][i]);
    BaryonField[EtaNum][i] = 0;
  }
  return SUCCESS;
}
#endif

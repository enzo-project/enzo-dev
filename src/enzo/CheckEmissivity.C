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
int grid::CheckEmissivity(){
  if( MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  //int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  int EtaNum;

  // Copying what's done in Grid_IdentifyPhysicalQuantities to check
  if ((EtaNum = FindField(Emissivity0, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find Emissivity0.\n");
    return FAIL;
  }
  int emis_x=0, emis_y=0, emis_z=0, emis_index=0;
    for(int emis_z=0; emis_z< GridDimension[2]; emis_z++)
      for(int emis_y=0; emis_y< GridDimension[1]; emis_y++)
	for(int emis_x=0; emis_x< GridDimension[0]; emis_x++){
	  emis_index = emis_x + GridDimension[0]*(emis_y + GridDimension[1]*emis_z);
	  if(BaryonField[EtaNum][emis_index] != 0){
	    printf("Check Emiss %0.12"GSYM" at %"ISYM" %"ISYM" %"ISYM"\n", BaryonField[EtaNum][emis_index], emis_x, emis_y, emis_z);
	    printf("This is processor %"ISYM" in EvolveLevel\n", MyProcessorNumber);
	  }
	}
    //  printf("DONE \n");
    return SUCCESS;
}
#endif

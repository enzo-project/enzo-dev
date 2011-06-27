/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE SHOCK SPECIES FIELDS)
/
/  written by: Samuel Skillman
/  date:       May, 2008
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
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
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
int grid::IdentifyShockSpeciesFields(int &MachNum, int &PSTempNum, int &PSDenNum)
{
  MachNum = PSTempNum = PSDenNum = 0;

  // Basic: Mach, CR Protons
  if (ShockMethod) {
    MachNum  = FindField(Mach, FieldType, NumberOfBaryonFields);
    if (StorePreShockFields){
      PSTempNum= FindField(PreShockTemperature, FieldType, NumberOfBaryonFields);
      PSDenNum = FindField(PreShockDensity, FieldType, NumberOfBaryonFields);    
    }     
  }

  if ((MachNum < 0) || (PSTempNum < 0) || (PSDenNum < 0)) {
    fprintf(stderr,"Error identifying species for ShockMethod = %"ISYM" MachNum= %"ISYM" PSTempNum = %"ISYM" PSDenNum = %"ISYM" NBaryonFs = %"ISYM".\n",
	    ShockMethod,MachNum,PSTempNum,PSDenNum,NumberOfBaryonFields);
    return FAIL;
  }
  return SUCCESS;
}

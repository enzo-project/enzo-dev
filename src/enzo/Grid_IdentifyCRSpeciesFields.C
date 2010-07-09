/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE SPECIES FIELDS FOR SAM SKILLMAN'S COSMIC RAYS)
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
 
int grid::IdentifyCRSpeciesFields(int &MachNum,int &CRNum, 
				  int &PSTempNum, int &PSDenNum)
{
  MachNum = CRNum = PSTempNum = PSDenNum = 0;

  // Basic: Mach, CR Protons
  if (CRModel) {
    MachNum  = FindField(Mach, FieldType, NumberOfBaryonFields);
    CRNum    = FindField(CRDensity, FieldType, NumberOfBaryonFields);
    if (StorePreShockFields){
      PSTempNum= FindField(PreShockTemperature, FieldType, NumberOfBaryonFields);
      PSDenNum = FindField(PreShockDensity, FieldType, NumberOfBaryonFields);    
    }     
  }

  if ((MachNum < 0) || (CRNum < 0) || (PSTempNum < 0) || (PSDenNum < 0)) {
    fprintf(stderr,"Error identifying species for CRModel = %"ISYM" MachNum= %"ISYM" CRNum = %"ISYM" PSTempNum = %"ISYM" PSDenNum = %"ISYM" NBaryonFs = %"ISYM".\n",
	    CRModel,MachNum,CRNum,PSTempNum,PSDenNum,NumberOfBaryonFields);
    return FAIL;
  }
  return SUCCESS;
}

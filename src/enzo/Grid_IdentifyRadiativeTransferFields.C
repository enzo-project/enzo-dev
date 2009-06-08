/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE MULTI-SPECIES FIELDS)
/
/  written by: Tom Abel (adopted from Grid_IdentifySpeciesFields)
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: Return indeces of the relevant radiative transfer fields
/
/  NOTE: 
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

int FindField(int f, int farray[], int n);

int grid::IdentifyRadiativeTransferFields(int &kphHINum, int &gammaHINum,
					  int &kphHeINum, int &gammaHeINum,
					  int &kphHeIINum, int &gammaHeIINum,
					  int &kdissH2INum) 
{
  
  kphHINum = gammaHINum = kphHeINum = gammaHeINum = kphHeIINum = gammaHeIINum = kdissH2INum = 0;

  if (RadiativeTransfer) {

    kphHINum    = FindField(kphHI, FieldType, NumberOfBaryonFields); 
    gammaHINum  = FindField(gammaHI, FieldType, NumberOfBaryonFields);	  
    kphHeINum   = FindField(kphHeI, FieldType, NumberOfBaryonFields);	  
    gammaHeINum = FindField(gammaHeI, FieldType, NumberOfBaryonFields);	  
    kphHeIINum  = FindField(kphHeII, FieldType, NumberOfBaryonFields);	  
    gammaHeIINum= FindField(gammaHeII, FieldType,  NumberOfBaryonFields);     
    if (MultiSpecies > 1) 
      kdissH2INum = FindField(kdissH2I, FieldType, NumberOfBaryonFields);     


    if (kphHINum<0 || gammaHINum<0 || kphHeINum<0 || gammaHeINum<0 ||
	kphHeIINum<0 || gammaHeIINum<0 || kdissH2INum<0) {
      fprintf(stderr, "Grid_IdentifyRadiativeTransferFields: failed\n");
      fprintf(stderr, "kphHI=%"ISYM", gammaHI=%"ISYM", kphHeI=%"ISYM", gammaHeI=%"ISYM", "
	      "kphHeII=%"ISYM", gammaHeII=%"ISYM", kdissH2I=%"ISYM"\n",
	      kphHINum, gammaHINum, kphHeINum, gammaHeINum, kphHeIINum, 
	      gammaHeIINum, kdissH2INum);
      return FAIL;
    }
  }

  return SUCCESS;
}

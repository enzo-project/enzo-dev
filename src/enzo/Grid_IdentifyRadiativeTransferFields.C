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
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

int FindField(int f, int farray[], int n);

int grid::IdentifyRadiativeTransferFields(int &kphHINum, int &gammaNum,
					  int &kphHeINum, 
					  int &kphHeIINum,
					  int &kdissH2INum,
					  int &kphHMNum,
					  int &kdissH2IINum) 
{
  
  kphHINum = gammaNum = kphHeINum = kphHeIINum = kdissH2INum = kphHMNum = kdissH2IINum = 0;

  if ((RadiativeTransfer) || (RadiativeTransferFLD)) {

    kphHINum    = FindField(kphHI, FieldType, NumberOfBaryonFields); 
    gammaNum    = FindField(PhotoGamma, FieldType, NumberOfBaryonFields);	  
    if (RadiativeTransferHydrogenOnly == FALSE) {
      kphHeINum   = FindField(kphHeI, FieldType, NumberOfBaryonFields);	  
      kphHeIINum  = FindField(kphHeII, FieldType, NumberOfBaryonFields);	  
    }
    if (MultiSpecies > 1) {
      kdissH2INum = FindField(kdissH2I, FieldType, NumberOfBaryonFields);     
      kphHMNum = FindField(kphHM, FieldType, NumberOfBaryonFields); 
      kdissH2IINum = FindField(kdissH2II, FieldType, NumberOfBaryonFields);   
    }
    if (kphHINum<0 || gammaNum<0 || kphHeINum<0 || kphHeIINum<0 || 
	kdissH2INum<0 || kphHMNum<0 || kdissH2IINum<0) {
      ENZO_VFAIL("Could not identify a RadiativeTransferField!  kphHI=%"ISYM", gamma=%"ISYM", "
		 "kphHeI=%"ISYM", kphHeII=%"ISYM", kdissH2I=%"ISYM", kphHM=%"ISYM", "
		 "kdissH2II=%"ISYM"\n", 
		 kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);
    }
  }

  return SUCCESS;
}

int grid::IdentifyRadiativeTransferFieldsMD(int &kdissHDINum,
					    int &kphCINum,
					    int &kphOINum,
					    int &kdissCONum,
					    int &kdissOHNum,
					    int &kdissH2ONum) 
{
  
  kdissHDINum = kphCINum = kphOINum = kdissCONum = kdissOHNum = kdissH2ONum = 0;

  if ((RadiativeTransfer) || (RadiativeTransferFLD)) {

    if (MultiSpecies > 2 && RadiativeTransferHDIDiss) {
      kdissHDINum = FindField(kdissHDI, FieldType, NumberOfBaryonFields);     
    }
    if (MetalChemistry && RadiativeTransferMetalIon) {
      kphCINum  = FindField(kphCI , FieldType, NumberOfBaryonFields);     
      kphOINum  = FindField(kphOI , FieldType, NumberOfBaryonFields);     
    }
    if (MetalChemistry && RadiativeTransferMetalDiss) {
      kdissCONum  = FindField(kdissCO , FieldType, NumberOfBaryonFields);     
      kdissOHNum  = FindField(kdissOH , FieldType, NumberOfBaryonFields);     
      kdissH2ONum = FindField(kdissH2O, FieldType, NumberOfBaryonFields);     
    }
    if (kdissHDINum<0 || kphCINum<0 || kphOINum<0 || kdissCONum<0 || kdissOHNum<0 || kdissH2ONum<0) {
      ENZO_VFAIL("Could not identify a RadiativeTransferFieldMD!  kdissHDI=%"ISYM", "
		 "kphCI=%"ISYM", kphOI=%"ISYM", " 
		 "kdissCO=%"ISYM", kdissOH=%"ISYM", kdissH2O=%"ISYM"\n", 
                 kdissHDINum, kphCINum, kphOINum, kdissCONum, kdissOHNum, kdissH2ONum);
    }
  }

  return SUCCESS;
}

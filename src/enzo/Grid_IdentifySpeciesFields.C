/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE MULTI-SPECIES FIELDS)
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
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
 
 
 
int grid::IdentifySpeciesFields(int &DeNum, int &HINum, int &HIINum,
				int &HeINum, int &HeIINum, int &HeIIINum,
                                int &HMNum, int &H2INum, int &H2IINum,
                                int &DINum, int &DIINum, int &HDINum)
{
 
  DeNum = HINum = HIINum = HeINum = HeIINum = HeIIINum = 0;
  HMNum = H2INum = H2IINum = DINum = DIINum = HDINum = 0;
 
  /* Find Fields for the 6-species model. */
 
  DeNum    = FindField(ElectronDensity, FieldType, NumberOfBaryonFields);
  HINum    = FindField(HIDensity, FieldType, NumberOfBaryonFields);
  HIINum   = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
  HeINum   = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
  HeIINum  = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
  HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);
 
  /* Error if any not found. */
 
  if (DeNum < 0 || HINum < 0 || HIINum < 0 || HeINum < 0 || HeIINum < 0 ||
      HeIIINum < 0) {
    fprintf(stderr, "De=%"ISYM", HI=%"ISYM", HII=%"ISYM", HeI=%"ISYM", HeII=%"ISYM", HeIII=%"ISYM"\n",
	    DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum);
    return FAIL;
  }
 
  /* Find Fields for the 9-species model. */
 
  if (MultiSpecies > 1) {
    HMNum    = FindField(HMDensity, FieldType, NumberOfBaryonFields);
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
    H2IINum  = FindField(H2IIDensity, FieldType, NumberOfBaryonFields);
 
    if (HMNum < 0 || H2INum < 0 || H2IINum < 0) {
      fprintf(stderr, "H2 related field missing.\n");
      return FAIL;
    }
 
  }
 
  /* Find Fields for the 12-species model. */
 
  if (MultiSpecies > 2) {
    DINum   = FindField(DIDensity, FieldType, NumberOfBaryonFields);
    DIINum  = FindField(DIIDensity, FieldType, NumberOfBaryonFields);
    HDINum  = FindField(HDIDensity, FieldType, NumberOfBaryonFields);
 
    if (DINum < 0 || DIINum < 0 || HDINum < 0) {
      fprintf(stderr, "HD related field missing.\n");
      return FAIL;
    }
 
  }
 
  return SUCCESS;
}

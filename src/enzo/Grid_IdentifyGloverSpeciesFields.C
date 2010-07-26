/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE SPECIES FIELDS FOR SIMON GLOVER'S COOLING)
/
/  written by: Britton Smith
/  date:       September, 2007
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
 
int grid::IdentifyGloverSpeciesFields(int &HIINum,int &HINum,int &H2INum,
				      int &DINum,int &DIINum,int &HDINum,
				      int &HeINum,int &HeIINum,int &HeIIINum,
				      int &CINum,int &CIINum,int &OINum,
				      int &OIINum,int &SiINum,int &SiIINum,
				      int &SiIIINum,int &CHINum,int &CH2INum,
				      int &CH3IINum,int &C2INum,int &COINum,
				      int &HCOIINum,int &OHINum,int &H2OINum,
				      int &O2INum)
{
 
  HIINum = HINum = H2INum = DINum = DIINum = HDINum = HeINum = HeIINum = 0;
  HeIIINum = CINum = CIINum = OINum = OIINum = SiINum = SiIINum = SiIIINum = 0;
  CHINum = CH2INum = CH3IINum = C2INum = COINum = HCOIINum = OHINum = H2OINum = O2INum = 0;

  // Primordial: H, D, and He
  if (GloverChemistryModel == 1) {
    HIINum   = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
    HINum    = FindField(HIDensity, FieldType, NumberOfBaryonFields);
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
    DINum    = FindField(DIDensity, FieldType, NumberOfBaryonFields);
    DIINum   = FindField(DIIDensity, FieldType, NumberOfBaryonFields);
    HDINum   = FindField(HDIDensity, FieldType, NumberOfBaryonFields);
    HeINum   = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
    HeIINum  = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
    HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);
  }

  // Low Z: H, D, He, C, O, and Si
  else if (GloverChemistryModel == 2) {
    HIINum   = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
    HINum    = FindField(HIDensity, FieldType, NumberOfBaryonFields);
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
    DINum    = FindField(DIDensity, FieldType, NumberOfBaryonFields);
    DIINum   = FindField(DIIDensity, FieldType, NumberOfBaryonFields);
    HDINum   = FindField(HDIDensity, FieldType, NumberOfBaryonFields);
    HeINum   = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
    HeIINum  = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
    HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);
    CINum    = FindField(CIDensity, FieldType, NumberOfBaryonFields);
    CIINum   = FindField(CIIDensity, FieldType, NumberOfBaryonFields);
    OINum    = FindField(OIDensity, FieldType, NumberOfBaryonFields);
    OIINum   = FindField(OIIDensity, FieldType, NumberOfBaryonFields);
    SiINum   = FindField(SiIDensity, FieldType, NumberOfBaryonFields);
    SiIINum  = FindField(SiIIDensity, FieldType, NumberOfBaryonFields);
    SiIIINum = FindField(SiIIIDensity, FieldType, NumberOfBaryonFields);
  }

  // Molecular: H, D, He, C, O, Si, plus molecules.
  else if (GloverChemistryModel == 3) {
    HIINum   = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
    HINum    = FindField(HIDensity, FieldType, NumberOfBaryonFields);
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
    DINum    = FindField(DIDensity, FieldType, NumberOfBaryonFields);
    DIINum   = FindField(DIIDensity, FieldType, NumberOfBaryonFields);
    HDINum   = FindField(HDIDensity, FieldType, NumberOfBaryonFields);
    HeINum   = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
    HeIINum  = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
    HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);
    CINum    = FindField(CIDensity, FieldType, NumberOfBaryonFields);
    CIINum   = FindField(CIIDensity, FieldType, NumberOfBaryonFields);
    OINum    = FindField(OIDensity, FieldType, NumberOfBaryonFields);
    OIINum   = FindField(OIIDensity, FieldType, NumberOfBaryonFields);
    SiINum   = FindField(SiIDensity, FieldType, NumberOfBaryonFields);
    SiIINum  = FindField(SiIIDensity, FieldType, NumberOfBaryonFields);
    SiIIINum = FindField(SiIIIDensity, FieldType, NumberOfBaryonFields);
    CHINum   = FindField(CHIDensity, FieldType, NumberOfBaryonFields);
    CH2INum  = FindField(CH2IDensity, FieldType, NumberOfBaryonFields);
    CH3IINum = FindField(CH3IIDensity, FieldType, NumberOfBaryonFields);
    C2INum   = FindField(C2IDensity, FieldType, NumberOfBaryonFields);
    COINum   = FindField(COIDensity, FieldType, NumberOfBaryonFields);
    HCOIINum = FindField(HCOIIDensity, FieldType, NumberOfBaryonFields);
    OHINum   = FindField(OHIDensity, FieldType, NumberOfBaryonFields);
    H2OINum  = FindField(H2OIDensity, FieldType, NumberOfBaryonFields);
    O2INum   = FindField(O2IDensity, FieldType, NumberOfBaryonFields);
  }

  // Simple: H
  else if (GloverChemistryModel == 4) {
    HIINum   = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
    HINum    = FindField(HIDensity, FieldType, NumberOfBaryonFields);
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
  }

  // Slightly less simple: H, CO
  else if (GloverChemistryModel == 5) {
    HIINum   = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
    HINum    = FindField(HIDensity, FieldType, NumberOfBaryonFields);
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
    COINum   = FindField(COIDensity, FieldType, NumberOfBaryonFields);
  }

  // GMC: H, He, C, O, plus some molecules
  else if (GloverChemistryModel == 7) {
    HIINum   = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
    HINum    = FindField(HIDensity, FieldType, NumberOfBaryonFields);
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
    DINum    = FindField(DIDensity, FieldType, NumberOfBaryonFields);
    DIINum   = FindField(DIIDensity, FieldType, NumberOfBaryonFields);
    HDINum   = FindField(HDIDensity, FieldType, NumberOfBaryonFields);
    HeINum   = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
    HeIINum  = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
    HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);
    CINum    = FindField(CIDensity, FieldType, NumberOfBaryonFields);
    CIINum   = FindField(CIIDensity, FieldType, NumberOfBaryonFields);
    OINum    = FindField(OIDensity, FieldType, NumberOfBaryonFields);
    OIINum   = FindField(OIIDensity, FieldType, NumberOfBaryonFields);
    CHINum   = FindField(CHIDensity, FieldType, NumberOfBaryonFields);
    CH2INum  = FindField(CH2IDensity, FieldType, NumberOfBaryonFields);
    CH3IINum = FindField(CH3IIDensity, FieldType, NumberOfBaryonFields);
    C2INum   = FindField(C2IDensity, FieldType, NumberOfBaryonFields);
    COINum   = FindField(COIDensity, FieldType, NumberOfBaryonFields);
    HCOIINum = FindField(HCOIIDensity, FieldType, NumberOfBaryonFields);
    OHINum   = FindField(OHIDensity, FieldType, NumberOfBaryonFields);
    H2OINum  = FindField(H2OIDensity, FieldType, NumberOfBaryonFields);
    O2INum   = FindField(O2IDensity, FieldType, NumberOfBaryonFields);
  }


  if ((HIINum < 0) || (HINum < 0) || (H2INum < 0) || (DINum < 0) || (DIINum < 0) || 
      (HDINum < 0) || (HeINum < 0) || (HeIINum < 0) || (HeIIINum < 0) || (CINum < 0) || 
      (CIINum < 0) || (OINum < 0) || (OIINum < 0) || (SiINum < 0) || (SiIINum < 0) || 
      (SiIIINum < 0) || (CHINum < 0) || (CH2INum < 0) || (CH3IINum < 0) || (C2INum < 0) || 
      (COINum < 0) || (HCOIINum < 0) || (OHINum < 0) || (H2OINum < 0) || (O2INum < 0)) {


    ENZO_VFAIL("Error identifying species for GloverChemistryModel = %"ISYM".\n",
	    GloverChemistryModel)

  }

  return SUCCESS;
}

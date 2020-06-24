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
    ENZO_VFAIL("De=%"ISYM", HI=%"ISYM", HII=%"ISYM", HeI=%"ISYM", HeII=%"ISYM", HeIII=%"ISYM"\n",
	    DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum)
  }
 
  /* Find Fields for the 9-species model. */
 
  if (MultiSpecies > 1) {
    HMNum    = FindField(HMDensity, FieldType, NumberOfBaryonFields);
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
    H2IINum  = FindField(H2IIDensity, FieldType, NumberOfBaryonFields);
 
    if (HMNum < 0 || H2INum < 0 || H2IINum < 0) {
      ENZO_FAIL("H2 related field missing.\n");
    }
 
  }
 
  /* Find Fields for the 12-species model. */
 
  if (MultiSpecies > 2) {
    DINum   = FindField(DIDensity, FieldType, NumberOfBaryonFields);
    DIINum  = FindField(DIIDensity, FieldType, NumberOfBaryonFields);
    HDINum  = FindField(HDIDensity, FieldType, NumberOfBaryonFields);
 
    if (DINum < 0 || DIINum < 0 || HDINum < 0) {
      ENZO_FAIL("HD related field missing.\n");

    }
 
  }
 
  return SUCCESS;
}

#ifdef GRACKLE_MD
int grid::IdentifySpeciesFieldsMD( int &HeHIINum , int &DMNum   , int &HDIINum
                                 , int &CINum    , int &CIINum  , int &CONum     , int &CO2Num
                                 , int &OINum    , int &OHNum   , int &H2ONum    , int &O2Num
                                 , int &SiINum   , int &SiOINum , int &SiO2INum
                                 , int &CHNum    , int &CH2Num  , int &COIINum   , int &OIINum
                                 , int &OHIINum  , int &H2OIINum, int &H3OIINum  , int &O2IINum
                                 , int &MgNum    , int &AlNum   , int &SNum      , int &FeNum
                                 , int &SiMNum   , int &FeMNum  , int &Mg2SiO4Num
                                 , int &MgSiO3Num, int &Fe3O4Num, int &ACNum
                                 , int &SiO2DNum , int &MgONum  , int &FeSNum    , int &Al2O3Num)
{

    HeHIINum = DMNum     = HDIINum
  = CINum    =  CIINum   = CONum      = CO2Num    = OINum    = OHNum
  = H2ONum   =  O2Num    = SiINum     = SiOINum   = SiO2INum
  = CHNum    =  CH2Num   = COIINum    = OIINum    = OHIINum  = H2OIINum = H3OIINum = O2IINum
  = MgNum    =  AlNum    = SNum       = FeNum
  = SiMNum   =  FeMNum   = Mg2SiO4Num = MgSiO3Num = Fe3O4Num
  = ACNum    =  SiO2DNum = MgONum     = FeSNum    = Al2O3Num
  = 0;

  if (MultiSpecies > 3) {
    HeHIINum   = FindField(HeHIIDensity  , FieldType, NumberOfBaryonFields);
    DMNum      = FindField(DMDensity     , FieldType, NumberOfBaryonFields);
    HDIINum    = FindField(HDIIDensity   , FieldType, NumberOfBaryonFields);
    if(HeHIINum  < 0 || DMNum   < 0 || HDIINum     < 0) {
      ENZO_FAIL("D, He related field missing.\n");
    }
  }

  if (MetalChemistry > 0) {
    CINum      = FindField(CIDensity     , FieldType, NumberOfBaryonFields);
    CIINum     = FindField(CIIDensity    , FieldType, NumberOfBaryonFields);
    CONum      = FindField(COIDensity    , FieldType, NumberOfBaryonFields);
    CO2Num     = FindField(CO2IDensity   , FieldType, NumberOfBaryonFields);
    OINum      = FindField(OIDensity     , FieldType, NumberOfBaryonFields);
    OHNum      = FindField(OHIDensity    , FieldType, NumberOfBaryonFields);
    if(CINum  < 0 || CIINum   < 0 || CONum      < 0 || CO2Num    < 0 || OINum    < 0 || OHNum < 0 ) {
      ENZO_FAIL("C related field missing.\n");
    }

    H2ONum     = FindField(H2OIDensity   , FieldType, NumberOfBaryonFields);
    O2Num      = FindField(O2IDensity    , FieldType, NumberOfBaryonFields);
    SiINum     = FindField(SiIDensity    , FieldType, NumberOfBaryonFields);
    SiOINum    = FindField(SiOIDensity   , FieldType, NumberOfBaryonFields);
    SiO2INum   = FindField(SiO2IDensity  , FieldType, NumberOfBaryonFields);
    if(H2ONum < 0 || O2Num    < 0 || SiINum     < 0 || SiOINum   < 0 || SiO2INum < 0 ) {
      ENZO_FAIL("Si related field missing.\n");
    }

    CHNum      = FindField(CHIDensity    , FieldType, NumberOfBaryonFields);
    CH2Num     = FindField(CH2IDensity   , FieldType, NumberOfBaryonFields);
    COIINum    = FindField(COIIDensity   , FieldType, NumberOfBaryonFields);
    OIINum     = FindField(OIIDensity    , FieldType, NumberOfBaryonFields);
    OHIINum    = FindField(OHIIDensity   , FieldType, NumberOfBaryonFields);
    H2OIINum   = FindField(H2OIIDensity  , FieldType, NumberOfBaryonFields);
    H3OIINum   = FindField(H3OIIDensity  , FieldType, NumberOfBaryonFields);
    O2IINum    = FindField(O2IIDensity   , FieldType, NumberOfBaryonFields);
    if(CHNum  < 0 || CH2Num   < 0 || COIINum    < 0 || OIINum    < 0 || OHIINum  < 0 || H2OIINum < 0 || H3OIINum < 0 || O2IINum < 0 ) {
      ENZO_FAIL("C related field missing.\n");
    }
  }

  if (GrainGrowth) {
    MgNum      = FindField(MgDensity     , FieldType, NumberOfBaryonFields);
    AlNum      = FindField(AlDensity     , FieldType, NumberOfBaryonFields);
    SNum       = FindField(SDensity      , FieldType, NumberOfBaryonFields);
    FeNum      = FindField(FeDensity     , FieldType, NumberOfBaryonFields);
    if(MgNum  < 0 || AlNum    < 0 || SNum       < 0 || FeNum < 0 ) {
      ENZO_FAIL("Metal related field missing.\n");
    }

    SiMNum     = FindField(SiMDensity    , FieldType, NumberOfBaryonFields);
    FeMNum     = FindField(FeMDensity    , FieldType, NumberOfBaryonFields);
    Mg2SiO4Num = FindField(Mg2SiO4Density, FieldType, NumberOfBaryonFields);
    MgSiO3Num  = FindField(MgSiO3Density , FieldType, NumberOfBaryonFields);
    Fe3O4Num   = FindField(Fe3O4Density  , FieldType, NumberOfBaryonFields);
    if(SiMNum < 0 || FeMNum   < 0 || Mg2SiO4Num < 0 || MgSiO3Num < 0 || Fe3O4Num < 0 ) {
      ENZO_FAIL("Dust related field missing.\n");
    }

    ACNum      = FindField(ACDensity     , FieldType, NumberOfBaryonFields);
    SiO2DNum   = FindField(SiO2DDensity  , FieldType, NumberOfBaryonFields);
    MgONum     = FindField(MgODensity    , FieldType, NumberOfBaryonFields);
    FeSNum     = FindField(FeSDensity    , FieldType, NumberOfBaryonFields);
    Al2O3Num   = FindField(Al2O3Density  , FieldType, NumberOfBaryonFields);
    if(ACNum  < 0 || SiO2DNum < 0 || MgONum     < 0 || FeSNum    < 0 || Al2O3Num < 0 ) {
      ENZO_FAIL("Dust related field missing.\n");
    }
  }

  return SUCCESS;

}
#endif

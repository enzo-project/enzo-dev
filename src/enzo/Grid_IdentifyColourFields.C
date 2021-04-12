/***********************************************************************
/
/  GRID CLASS (IDENTIFY COLOUR FIELDS THAT ARE FREQUENTLY USED)
/
/  written by: Ji-hoon Kim
/  date:       October, 2009
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
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);

int grid::IdentifyColourFields(int &SNColourNum, int &MetalNum, 
			       int &MetalIaNum, int &MetalIINum, int &MBHColourNum,
			       int &Galaxy1ColourNum, int &Galaxy2ColourNum)
{
 
  SNColourNum = MetalNum = MetalIaNum = MBHColourNum = Galaxy1ColourNum = 
    MetalIINum = Galaxy2ColourNum = 0;
 
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  MetalIaNum = FindField(MetalSNIaDensity, FieldType, NumberOfBaryonFields);
  MetalIINum = FindField(MetalSNIIDensity, FieldType, NumberOfBaryonFields);
  MBHColourNum = FindField(MBHColour, FieldType, NumberOfBaryonFields);
  Galaxy1ColourNum = FindField(Galaxy1Colour, FieldType, NumberOfBaryonFields);
  Galaxy2ColourNum = FindField(Galaxy2Colour, FieldType, NumberOfBaryonFields);

  /*
  if ((SNColourNum < 0) && (MetalNum < 0) && (MBHColourNum < 0) && 
      (Galaxy1ColourNum < 0) && (Galaxy2ColourNum < 0)) {
"No colour field identified;    ENZO_FAIL( while this could happen, check if it was expected.\n");

  }
  */

  return SUCCESS;
}

int grid::IdentifyExtraTypeFields(
      int &ExtraType0Num, int &ExtraType1Num, int &ExtraType2Num, int &ExtraType3Num,
      int &ExtraType4Num, int &ExtraType5Num, int &ExtraType6Num, int &ExtraType7Num,
      int &ExtraType8Num, int &ExtraType9Num, int &ExtraType10Num,int &ExtraType11Num
        )
{
    ExtraType0Num = ExtraType1Num = ExtraType2Num = ExtraType3Num
  = ExtraType4Num = ExtraType5Num = ExtraType6Num = ExtraType7Num
  = ExtraType8Num = ExtraType9Num = ExtraType10Num= ExtraType11Num
  = 0;

  ExtraType0Num = FindField(ExtraType0, FieldType, NumberOfBaryonFields);
  ExtraType1Num = FindField(ExtraType1, FieldType, NumberOfBaryonFields);
  ExtraType2Num = FindField(ExtraType2, FieldType, NumberOfBaryonFields);
  ExtraType3Num = FindField(ExtraType3, FieldType, NumberOfBaryonFields);
  ExtraType4Num = FindField(ExtraType4, FieldType, NumberOfBaryonFields);
  ExtraType5Num = FindField(ExtraType5, FieldType, NumberOfBaryonFields);
  ExtraType6Num = FindField(ExtraType6, FieldType, NumberOfBaryonFields);
  ExtraType7Num = FindField(ExtraType7, FieldType, NumberOfBaryonFields);
  ExtraType8Num = FindField(ExtraType8, FieldType, NumberOfBaryonFields);
  ExtraType9Num = FindField(ExtraType9, FieldType, NumberOfBaryonFields);
  ExtraType10Num= FindField(ExtraType10,FieldType, NumberOfBaryonFields);
  ExtraType11Num= FindField(ExtraType11,FieldType, NumberOfBaryonFields);

  return SUCCESS;
}

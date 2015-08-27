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


/***********************************************************************
/
/  GRID CLASS (RETURNS THE NUMBER OF COLOURS)
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);

int grid::SetNumberOfColours(void)
{

  int _nc = NColor; // alter only a local variable and return it.
  int MetalNum, SNColourNum;

  // If already counted, return.
  if (_nc != INT_UNDEFINED)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  _nc = 0;

  // Metallicity fields

  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1) {
    _nc++;
    if (MultiMetals || TestProblemData.MultiMetals)
      _nc += 2;
  }

  // Supernova "colour"
  
  if ((SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields)) != -1)
    _nc++;
  
  /* Treat these colour (i.e. metal) fields as species fields in the
     MUSCL solvers. */

  NColor = 0;
  NSpecies += _nc;

  return SUCCESS;

}

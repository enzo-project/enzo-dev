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

  // If already counted, return.
  if (_nc != INT_UNDEFINED)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  _nc = 0;

  /* Count colours */  

  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum; 

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    return FAIL;
  }
  
  if (MetalNum != -1) {
    _nc++;
    if (MultiMetals || TestProblemData.MultiMetals) {
      _nc += 2;
    }
  }

  if (SNColourNum      != -1) _nc++;
  if (MBHColourNum     != -1) _nc++;
  if (Galaxy1ColourNum != -1) _nc++;
  if (Galaxy2ColourNum != -1) _nc++;

  /*   //#####
  int MetalNum, SNColourNum;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1) {
    _nc++;
    if (MultiMetals || TestProblemData.MultiMetals)
      _nc += 2;
  }

  if ((SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields)) != -1)
    _nc++;
  */

  
  /* Treat these colour (i.e. metal) fields as species fields in the
     MUSCL solvers. */

  NColor = 0;
  NSpecies += _nc;

  return SUCCESS;

}

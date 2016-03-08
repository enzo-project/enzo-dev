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
#include "ErrorExceptions.h"
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

  int _nc = NColor;  // alter only a local variable and return it.

  if (_nc != INT_UNDEFINED)   // If already counted, return.
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  _nc = 0;

  /* Count colours */  

  int SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, Galaxy1ColourNum, 
    Galaxy2ColourNum; 

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    return FAIL;
  }
  
  /*
  fprintf(stdout, "grid:SetNumberOfColours: %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" \n", 
  	  SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum); 
  */

  if (MetalNum != -1) {
    _nc++;
    if (StarMakerTypeIaSNe) _nc++;
    if (StarMakerTypeIISNeMetalField) _nc++;
    if (MultiMetals || TestProblemData.MultiMetals) {
      _nc += 2;
    }
  }

  if (SNColourNum      != -1) _nc++;
  /*   //##### These fields are currently not being used and only causing interpolation problems
  if (MBHColourNum     != -1) _nc++;
  if (Galaxy1ColourNum != -1) _nc++;
  if (Galaxy2ColourNum != -1) _nc++;
  */


  /* Treat these colour (i.e. metal) fields as species fields in the
     MUSCL solvers. */

  NColor = 0;  
  NSpecies += _nc;

  
  if (NSpecies != 0 && MultiSpecies == 0)
    NoMultiSpeciesButColors = 1;

  if (debug) 
    fprintf(stdout, "grid:SetNumberOfColours: NEQ_HYDRO = %"ISYM", NSpecies = %"ISYM", NColor = %"ISYM"\n", 
  	  NEQ_HYDRO, NSpecies, NColor); 


  return SUCCESS;

}

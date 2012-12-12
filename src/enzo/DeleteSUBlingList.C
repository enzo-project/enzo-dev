/*****************************************************************************
 *                                                                           *
 * Copyright 2005 Rick Wagner                                                *
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  DELETE SUBLINGLIST
/
/  written by: Rick Wagner
/  date:       May, 2005
/
/
/  PURPOSE: Try to forget about all of  dirty little subgrids that
/           touched us.
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
 
 
int DeleteSUBlingList( int NumberOfGrids,
		      LevelHierarchyEntry **SUBlingList)
{

  if( FluxCorrection != TRUE )
    return SUCCESS;
  LevelHierarchyEntry *LastEntry, *NextEntry;
 
  /* Add all the SUBgrids to the list */
  for (int grid = 0; grid < NumberOfGrids; grid++){
    LastEntry = SUBlingList[grid];
    while( LastEntry ){
      NextEntry = LastEntry->NextGridThisLevel;
      delete LastEntry;
      LastEntry = NextEntry;
    }
  }
 
  delete [] SUBlingList;
 
  return SUCCESS;
 
}

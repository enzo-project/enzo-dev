 
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
/  CREATE SUBLINGLIST
/
/  written by: Rick Wagner
/  date:       May, 2005
/
/
/  PURPOSE: Create a list, for each grid, of all grids on the next finer level
/           that share a face ONLY.
/           Note that this configuration does NOT recognize grids that are strict subgrids
/           AT ALL, so if you're running with fewer than 2 processors in each direction
/           AND periodic boundary conditions, there will be grids at the domain edge that may
/           cause conservation errors (since they'll both share a face AND be subgrids.)
/           In that case, the user will need to do some work to a.) get those grids included
/           in this list AND b.) avoid correcting that grid multiple times.  dcc.
/
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
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
 
#ifdef FLUX_FIX
int CreateSUBlingList(TopGridData *MetaData,
		      HierarchyEntry *Grids[],
		      int NumberOfGrids,
		      LevelHierarchyEntry ***SUBlingList)
{
 
  int grid, othergrid;
  /* Create a SUBling list of the subgrids */
 
  *SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];
  LevelHierarchyEntry *NextEntry, *LastEntry = NULL;
  HierarchyEntry *NextGrid;
 
  /* Add all the SUBgrids to the list */
 
  for (grid = 0; grid < NumberOfGrids; grid++){
 
    (*SUBlingList)[grid] = NULL;
    for (othergrid = 0; othergrid < NumberOfGrids; othergrid++){
      if( othergrid != grid ){
 
	if (Grids[grid]->GridData->
	    CheckForSharedFace(Grids[othergrid]->GridData,
				    MetaData->LeftFaceBoundaryCondition,
				    MetaData->RightFaceBoundaryCondition
				    ) == TRUE){
 
	  NextGrid = Grids[othergrid]->NextGridNextLevel;
	  while( NextGrid ){
 
	    if (Grids[grid]->GridData->
		CheckForSharedFace(NextGrid->GridData,
					MetaData->LeftFaceBoundaryCondition,
					MetaData->RightFaceBoundaryCondition
					) == TRUE){
 
	      if( LastEntry != NULL ){
		LastEntry->NextGridThisLevel = new LevelHierarchyEntry;
		LastEntry = LastEntry->NextGridThisLevel;
	      }else{
		(*SUBlingList)[grid] = new LevelHierarchyEntry;
		LastEntry = (*SUBlingList)[grid];
	      }
	      LastEntry->GridHierarchyEntry = NextGrid;
	      LastEntry->GridData = NextGrid->GridData;
	      LastEntry->NextGridThisLevel = NULL;
	
	    }//match
 
	    NextGrid = NextGrid->NextGridThisLevel;
	  } /* while( NextGrid ) -- loop over subgrids */
	}/* if grids overlap */
 
      }
    } /* loop over other grids this level */
  } /* loop over all grids this level */
 
 
  return SUCCESS;
 
}
#endif

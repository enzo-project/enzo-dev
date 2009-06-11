/***********************************************************************
/
/  GRID CLASS (INHERIT BASIC GRID PROPERTIES FROM THE PARENT GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
//
//  Assign basic values to a grid (allocate fields)
//

#include <stdio.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
void grid::InheritProperties(grid *ParentGrid)
{
 
  /*  Set rank and current time */
 
  GridRank = ParentGrid->GridRank;
  Time     = ParentGrid->Time;
 
  /*  Baryons only: set up field quantities and allocate fields
       (we assume here the grid is uniform in each dimension) */
 
  NumberOfBaryonFields = ParentGrid->NumberOfBaryonFields;
 
//  printf("InheritProperties: NOBF %"ISYM"\n", NumberOfBaryonFields);
 
  for (int field = 0; field < NumberOfBaryonFields; field++)
    FieldType[field]      = ParentGrid->FieldType[field];
 
  CourantSafetyNumber    = ParentGrid->CourantSafetyNumber;
  PPMFlatteningParameter = ParentGrid->PPMFlatteningParameter;
  PPMDiffusionParameter  = ParentGrid->PPMDiffusionParameter;
  PPMSteepeningParameter = ParentGrid->PPMSteepeningParameter;
 
  /* For gravity, this must be a subgrid, so set the boundary accordingly. */
 
  GravityBoundaryType = SubGridIsolated;
 
}

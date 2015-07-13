/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (PREPARE OBJECT FROM FRIENDLY GRID TOPGRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

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
 
// Set the parameters of the external boundary from a grid passed in the
// argument line.  This routine must be a friend of class grid.
//
 
/* function prototypes */
 
void WriteListOfInts(FILE *fptr, int N, int nums[]);
 
int ExternalBoundary::Prepare(grid *TopGrid)
{
  BoundaryRank = TopGrid->GridRank;
  NumberOfBaryonFields = TopGrid->NumberOfBaryonFields;
  if (ParallelRootGridIO == TRUE)
    TopGrid->NumberOfBaryonFields = 0; /* bad kludge! */
 
  if (debug) {
    printf("ExtBndry: BoundaryRank = %"ISYM"\n", BoundaryRank);
    printf("ExtBndry: GridDimension = ");
    WriteListOfInts(stdout, BoundaryRank, TopGrid->GridDimension);
    printf("ExtBndry: NumberOfBaryonFields = %"ISYM"\n", NumberOfBaryonFields);
  }
 
  int size = 1;
  for (int i = 0; i < BoundaryRank; i++)
    size = size*TopGrid->GridDimension[i];
 
  for (int dim = 0; dim < BoundaryRank; dim++) {
    BoundaryDimension[dim] = TopGrid->GridDimension[dim];
    if (BoundaryDimension[dim] != 1)
      for (int field = 0; field < NumberOfBaryonFields; field++) {
	BoundaryFieldType[field] = TopGrid->FieldType[field];
	BoundaryType[field][dim][0]   = NULL;
	BoundaryType[field][dim][1]   = NULL;
	BoundaryValue[field][dim][0]  = NULL;
	BoundaryValue[field][dim][1]  = NULL;
      }
  }
 
  if(UseMHDCT)
    {
      for(int field=0;field<3; field++)
	for(int dim=0;dim<3; dim++){
	  MagneticBoundaryDims[field][dim]=TopGrid->GridDimension[dim];
	  
	  if(field == dim) MagneticBoundaryDims[field][dim]++;
	}

      for(int dim=0;dim<3;dim++)
	for(int field = 0; field<3; field++){
	  MagneticBoundaryValue[field][dim][0] = NULL;
	  MagneticBoundaryValue[field][dim][1] = NULL;
	}
    }

  return SUCCESS;
 
}

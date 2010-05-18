#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "Hierarchy.h"

void grid::FlagGridArray( HierarchyEntry ***GridArray, int *dx,
			  float *cell_width, HierarchyEntry *my_HE ){

  int index_start[MAX_DIMENSION];
  int index_end[MAX_DIMENSION];	 
  int i, j, k;
  
  for(i = 0; i < MAX_DIMENSION; i++){
    index_start[i] = (int)((GridLeftEdge[i] - DomainLeftEdge[i] + 
			    0.5 *cell_width[i] ) / cell_width[i]);
    index_end[i] = (int)((GridRightEdge[i] - DomainLeftEdge[i] + 
			  0.5 * cell_width[i] ) / cell_width[i]);
  }
  
  for (k = index_start[2]; k < index_end[2]; k++)
    for (j = index_start[1]; j < index_end[1]; j++)
      for (i = index_start[0]; i < index_end[0]; i++)
	(*GridArray)[i + j*dx[0] +
		     k*dx[0]*dx[1]] = my_HE;
  
}

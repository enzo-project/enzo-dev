/***********************************************************************
/
/  GRID CLASS (WRITE OUT GRID-TO-TASK MAP)
/
/  written by: Robert Harkness
/  date:       January 2006
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI 
#include <mpi.h>
#endif
#include <string.h>
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
void my_exit(int status);
 
// function prototypes
 
int grid::WriteTaskMap(FILE *fptr, char *base_name, int grid_id)
{
 
  int i, j, k;
  int file_status;
  int output_cube;

/*
  int crap1 = MyProcessorNumber;
  int crap2 = this->ProcessorNumber;
  int crap3 = grid_id;

  fprintf(stderr, "CRAP1 %"ISYM"  CRAP2 %"ISYM"  CRAP3 %"ISYM"\n", crap1, crap2, crap3);
*/

  if ( this->ProcessorNumber == MyProcessorNumber )
    fprintf(fptr, "Grid %8"ISYM"  PN %8"ISYM"\n", grid_id, MyProcessorNumber);

  return SUCCESS;
 
}

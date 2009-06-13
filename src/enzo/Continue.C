/***********************************************************************
/
/  ENZO Test for continue
/
/  written by: Robert Harkness
/  date:       February, 2007
/
************************************************************************/

#ifdef USE_MPI
#include <mpi.h>
#endif

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
#include "CosmologyParameters.h"
#include "version.def"


int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
void my_exit(int status);

void ContinueExecution(void)
{

  FILE *con;
  int flag, zero;

  zero = 0;
  flag = 1;

  if (MyProcessorNumber == 0) {

    con = fopen("ContinueFlag", "r");

    if (con != NULL) {
      fscanf(con, "Continue = %"ISYM, &flag);
      fclose(con);
    }

    fprintf(stderr, "Continuation Flag = %"ISYM"\n", flag);

  }

  CommunicationBroadcastValue(&flag, zero);

  if ( flag == FALSE ) {

    my_exit(zero);

  }

}

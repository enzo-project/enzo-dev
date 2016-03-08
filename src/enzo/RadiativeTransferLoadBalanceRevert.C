/***********************************************************************
/
/  COMMUNICATION ROUTINE: CLEANUP GRIDS AFTER EVOLVE PHOTON
/
/  written by: John Wise
/  date:       September, 2010
/  modified1:
/
/  PURPOSE: Cleanup "fake" replicated grids after ray tracing with RT 
/           load balancing on.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "communication.h"
#include "CommunicationUtilities.h"

double ReturnWallTime(void);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

#define MIN_LEVEL 1

int RadiativeTransferLoadBalanceRevert(HierarchyEntry **Grids[], int *NumberOfGrids)
{

  if (NumberOfProcessors == 1)
    return SUCCESS;

  int i, index, level, temp_proc, ori_proc, nph, GridsMoved, TotalNumberOfGrids;

  /* Send RadiationPresent fields */

  TotalNumberOfGrids = 0;
  for (level = MIN_LEVEL; level < MAX_DEPTH_OF_HIERARCHY; level++)
    TotalNumberOfGrids += NumberOfGrids[level];
  
  index = 0;
  int *RadiationPresent = new int[TotalNumberOfGrids];
  for (level = MIN_LEVEL; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (i = 0; i < NumberOfGrids[level]; i++) {
      temp_proc = Grids[level][i]->GridData->ReturnProcessorNumber();
      if (MyProcessorNumber == temp_proc)
	RadiationPresent[index] = 
	  Grids[level][i]->GridData->RadiationPresent();
      else
	RadiationPresent[index] = 0;
      index++;
    } // ENDFOR grids

  CommunicationAllSumValues(RadiationPresent, TotalNumberOfGrids);

  index = 0;
  for (level = MIN_LEVEL; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (i = 0; i < NumberOfGrids[level]; i++)
      Grids[level][i]->GridData->SetRadiation(RadiationPresent[index++]);

  delete [] RadiationPresent;

  /* Send updated baryon fields (species and energy in particular)
     back to the original processor */

  /* Now we know where the grids are going, transfer them. */

  for (level = MIN_LEVEL; level < MAX_DEPTH_OF_HIERARCHY; level++) {

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  GridsMoved = 0;
  for (i = 0; i < NumberOfGrids[level]; i++) {
    ori_proc = Grids[level][i]->GridData->ReturnOriginalProcessorNumber();
    temp_proc = Grids[level][i]->GridData->ReturnProcessorNumber();
    if (ori_proc != temp_proc) {
      Grids[level][i]->GridData->CommunicationMoveGrid(ori_proc, FALSE, FALSE);
      GridsMoved++;
    }
  }

  /* Send grids */

  CommunicationDirection = COMMUNICATION_SEND;

  for (i = 0; i < NumberOfGrids[level]; i++) {
    ori_proc = Grids[level][i]->GridData->ReturnOriginalProcessorNumber();
    temp_proc = Grids[level][i]->GridData->ReturnProcessorNumber();
    if (ori_proc != temp_proc) {
      if (RandomForcing)  //AK
	Grids[level][i]->GridData->AppendForcingToBaryonFields();
      Grids[level][i]->GridData->CommunicationMoveGrid(ori_proc, FALSE, FALSE);
    }
  }

  /* Receive grids */

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("CommunicationReceiveHandler() failed!\n");

  CommunicationBarrier();

//  if (MyProcessorNumber == ROOT_PROCESSOR && GridsMoved > 0)
//    printf("RevertPhotonLoadBalance[%d]: Number of grids moved = %"ISYM" out of %"ISYM"\n",
//	   level, GridsMoved, NumberOfGrids[level]);

  } // ENDFOR level

  /* Delete baryon fields in replicated grid */

  for (level = MIN_LEVEL; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (i = 0; i < NumberOfGrids[level]; i++) {
      ori_proc = Grids[level][i]->GridData->ReturnOriginalProcessorNumber();
      temp_proc = Grids[level][i]->GridData->ReturnProcessorNumber();
      if (MyProcessorNumber == temp_proc && ori_proc != temp_proc) {
	Grids[level][i]->GridData->DeleteAllFields();
	Grids[level][i]->GridData->DeleteSubgridMarker();
      }
    } // ENDFOR grids
  } // ENDFOR Level

  /* Set processor number back to its original number */

  for (level = MIN_LEVEL; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    for (i = 0; i < NumberOfGrids[level]; i++) {
      ori_proc = Grids[level][i]->GridData->ReturnOriginalProcessorNumber();
      temp_proc = Grids[level][i]->GridData->ReturnProcessorNumber();
      if (ori_proc != temp_proc)
	Grids[level][i]->GridData->SetProcessorNumber(ori_proc);
    } // ENDFOR grids

  } // ENDFOR level

  return SUCCESS;

}

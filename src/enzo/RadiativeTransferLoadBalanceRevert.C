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

int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

int RadiativeTransferLoadBalanceRevert(HierarchyEntry **Grids[], int *NumberOfGrids)
{

  int i, level, temp_proc, ori_proc, nph;

  /* Delete baryon fields in replicated grid */

  for (level = 1; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (i = 0; i < NumberOfGrids[level]; i++) {
      ori_proc = Grids[level][i]->GridData->ReturnOriginalProcessorNumber();
      temp_proc = Grids[level][i]->GridData->ReturnProcessorNumber();
      if (MyProcessorNumber == temp_proc && ori_proc != temp_proc) {
	Grids[level][i]->GridData->DeleteAllFields();
	Grids[level][i]->GridData->DeleteSubgridMarker();
      }
    } // ENDFOR grids
  } // ENDFOR Level

  /* Send photon packages from replicated grid to actual grid */

  PhotonPackageEntry *PP;

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (level = 1; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (i = 0; i < NumberOfGrids[level]; i++) {
      ori_proc = Grids[level][i]->GridData->ReturnOriginalProcessorNumber();
      temp_proc = Grids[level][i]->GridData->ReturnProcessorNumber();
      if (ori_proc != temp_proc) {
	PP = Grids[level][i]->GridData->ReturnPhotonPackagePointer();
	nph = Grids[level][i]->GridData->ReturnNumberOfPhotonPackages();
	if (PP->NextPackage != NULL)
	  Grids[level][i]->GridData->CommunicationSendPhotonPackages
	    (Grids[level][i]->GridData, ori_proc, nph, nph, &PP->NextPackage);
      } // ENDIF
    } // ENDFOR grids
  } // ENDFOR level

  /* Send photons */
  
  CommunicationDirection = COMMUNICATION_SEND;

  for (level = 1; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (i = 0; i < NumberOfGrids[level]; i++) {
      ori_proc = Grids[level][i]->GridData->ReturnOriginalProcessorNumber();
      temp_proc = Grids[level][i]->GridData->ReturnProcessorNumber();
      if (ori_proc != temp_proc) {
	PP = Grids[level][i]->GridData->ReturnPhotonPackagePointer();
	nph = Grids[level][i]->GridData->ReturnNumberOfPhotonPackages();
	if (PP->NextPackage != NULL)
	  Grids[level][i]->GridData->CommunicationSendPhotonPackages
	    (Grids[level][i]->GridData, ori_proc, nph, nph, &PP->NextPackage);
      } // ENDIF
    } // ENDFOR grids
  } // ENDFOR level

  /* Receive photons */

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("CommunicationReceiveHandler() failed!\n");

  /* Return processor number back to original */

  for (level = 1; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (i = 0; i < NumberOfGrids[level]; i++) {
      ori_proc = Grids[level][i]->GridData->ReturnOriginalProcessorNumber();
      temp_proc = Grids[level][i]->GridData->ReturnProcessorNumber();
      if (ori_proc != temp_proc)
	Grids[level][i]->GridData->SetProcessorNumber(ori_proc);
    } // ENDFOR grids
  } // ENDFOR level

  return SUCCESS;

}

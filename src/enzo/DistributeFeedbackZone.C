/***********************************************************************
/
/  Copy data from a 'fake' feedback zone grid back to 
/  the real grids
/
/  written by: Nathan Goldbaum
/  date:       June 2012
/
************************************************************************/

#ifdef USE_MPI
#endif

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "phys_constants.h"
#include "CommunicationUtilities.h"
#include "communication.h"

int CommunicationBufferPurge(void);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

int DistributeFeedbackZone(grid* FeedbackZone, HierarchyEntry** Grids, 
			   int NumberOfGrids, int SendField)
{
  int i,j;

  FLOAT ZeroVector[] = {0, 0, 0};

  /* Post receives */
  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (i = 0; i < NumberOfGrids; i++) 
    if (Grids[i]->GridData->CopyActiveZonesFromGrid(FeedbackZone,ZeroVector,SendField) == FAIL)
      ENZO_FAIL("FeedbackZone copy failed!\n");
    
  /* Send data */
    
  CommunicationDirection = COMMUNICATION_SEND;

  for (i = 0; i < NumberOfGrids; i++) 
    if (Grids[i]->GridData->CopyActiveZonesFromGrid(FeedbackZone,ZeroVector,SendField) == FAIL)
      ENZO_FAIL("FeedbackZone copy failed!\n");
  
  /* Receive data */
  
  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("CommunicationReceiveHandler() failed!\n");

#ifdef USE_MPI
  CommunicationBufferPurge();
#endif

  return SUCCESS;
}

/***********************************************************************
/
/  GRID CLASS (MOVE A GRID FROM ONE PROCESSOR TO ANOTHER)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  January, 2012 (John Wise) move fields and particles+ 
/              separately
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif

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
#include "communication.h"
 
/* function prototypes */
 
 
 
int grid::CommunicationMoveGrid2(int ToProcessor, int MoveParticles)
{

  int dim;
  int Zero[] = {0, 0, 0};
  FLOAT FZero[] = {0.0, 0.0, 0.0};

  //CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  if ((MyProcessorNumber == ProcessorNumber ||
       MyProcessorNumber == ToProcessor) &&
      ProcessorNumber != ToProcessor) {

    /* Copy particles. */

    if (NumberOfParticles > 0 && MoveParticles == TRUE)
      this->CommunicationSendParticles(this, ToProcessor, 0,
				       NumberOfParticles, 0);

    /* Copy stars */

    if (NumberOfStars > 0 && MoveParticles == TRUE)
      this->CommunicationSendStars(this, ToProcessor);

    /* Copy photon packages */

#ifdef TRANSFER
    PhotonPackageEntry *PP = PhotonPackages->NextPackage;
    if (PP != NULL)
      this->CommunicationSendPhotonPackages(this, ToProcessor, 
					    NumberOfPhotonPackages, 
					    NumberOfPhotonPackages, &PP);
#endif /* TRANSFER */    

    /* Delete fields on old grid. */
 
    if (MyProcessorNumber == ProcessorNumber && ProcessorNumber != ToProcessor &&
	(CommunicationDirection == COMMUNICATION_SEND ||
	 CommunicationDirection == COMMUNICATION_SEND_RECEIVE)) {
      if (MoveParticles) this->DeleteParticles();
    }
    
  } // ENDIF right processor
 
  /* Update processor number. */
  
  if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
    ProcessorNumber = ToProcessor;
 
  return SUCCESS;
}
 

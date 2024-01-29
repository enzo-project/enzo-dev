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
 
 
 
int grid::CommunicationMoveGrid1(int ToProcessor)
{

  int dim;
  int Zero[] = {0, 0, 0};
  FLOAT FZero[] = {0.0, 0.0, 0.0};

  //CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  if ((MyProcessorNumber == ProcessorNumber ||
       MyProcessorNumber == ToProcessor) &&
      ProcessorNumber != ToProcessor) {

    /* Copy baryons. */
 
    if (NumberOfBaryonFields > 0) {
#ifdef USE_MPI
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = this;
	CommunicationReceiveCallType[CommunicationReceiveIndex] = 16;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  CommunicationReceiveArgumentInt[dim][CommunicationReceiveIndex] =
	    GridDimension[dim];
      }
#endif
      this->CommunicationSendRegion(this, ToProcessor, ALL_FIELDS,
				    NEW_ONLY, Zero, GridDimension);
    }
 
    /* Delete fields on old grid. */
 
    if (MyProcessorNumber == ProcessorNumber && ProcessorNumber != ToProcessor &&
	(CommunicationDirection == COMMUNICATION_SEND ||
	 CommunicationDirection == COMMUNICATION_SEND_RECEIVE)) {
      this->DeleteAllButParticles();
    }
    
  } // ENDIF right processor
 
  /* Update processor number. */
  
//  if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
//    ProcessorNumber = ToProcessor;
 
  return SUCCESS;
}
 

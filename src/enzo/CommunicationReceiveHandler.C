/***********************************************************************
/
/  EVOLVE LEVEL ROUTINES (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: This routine processes the receives stored in the 
/           CommunicationReceive stack.  Each receive is tagged with a 
/           type which indicates which method to call 
/           (and a record of the arguments).
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "TopGridData.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "communication.h"
 
#ifdef USE_MPI
static MPI_Arg ListOfIndices[MAX_RECEIVE_BUFFERS];
static MPI_Status ListOfStatuses[MAX_RECEIVE_BUFFERS];
#endif /* USE_MPI */

double ReturnWallTime(void);

int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[],
				int NumberOfSubgrids[],
				int FluxFlag, TopGridData* MetaData)
{

#ifdef USE_MPI

  int NoErrorSoFar = TRUE;
  int Zero[] = {0, 0, 0};

  /* Set the communication mode. */

  CommunicationDirection = COMMUNICATION_RECEIVE;

//  printf("P(%"ISYM") in CRH with %"ISYM" requests\n", MyProcessorNumber,
//  	 CommunicationReceiveIndex);

  MPI_Arg NumberOfCompleteRequests, TotalReceives;
  int ReceivesCompletedToDate = 0, index, errcode, SUBling, level,
    igrid, isubgrid, dim, FromStart, FromNumber, ToStart, ToNumber;
  int GridDimension[MAX_DIMENSION];
  FLOAT EdgeOffset[MAX_DIMENSION];
  grid *grid_one, *grid_two;
  TotalReceives = CommunicationReceiveIndex;
#ifdef TRANSFER
  PhotonPackageEntry *PP;
#endif

  /* Define a temporary flux holder for the refined fluxes. */

  fluxes SubgridFluxesRefined;
  InitializeFluxes(&SubgridFluxesRefined);

  while (ReceivesCompletedToDate < TotalReceives) {

    /* Call the MPI wait handler. */

    float time1 = ReturnWallTime();

    MPI_Waitsome(TotalReceives, CommunicationReceiveMPI_Request,
		 &NumberOfCompleteRequests, ListOfIndices, ListOfStatuses);
//    printf("MPI: %"ISYM" %"ISYM" %"ISYM"\n", TotalReceives, 
//	   ReceivesCompletedToDate, NumberOfCompleteRequests);

    CommunicationTime += ReturnWallTime() - time1;

    /* Error check */

#if 0
    if (NumberOfCompleteRequests == MPI_UNDEFINED) {
      ENZO_FAIL("Error in MPI_Waitsome\n");
    }
#endif

    /* Should loop over newly received completions and check error msgs now. */
    for (index = 0; index < NumberOfCompleteRequests; index++)
      if (ListOfStatuses[index].MPI_ERROR != 0) {
	if (NoErrorSoFar) {
	  fprintf(stderr, "MPI Error on processor %"ISYM". "
		  "Error number %"ISYM" on request %"ISYM"\n",
		  MyProcessorNumber, ListOfStatuses[index].MPI_ERROR, index);
	  NoErrorSoFar = FALSE;
	}
	fprintf(stdout, "P(%"ISYM") index %"ISYM" -- mpi error %"ISYM"\n", 
		MyProcessorNumber, index, ListOfStatuses[index].MPI_ERROR);
	fprintf(stdout, "%"ISYM": Type = %"ISYM", Grid1 = %x, Request = %"ISYM", "
		"DependsOn = %"ISYM"\n", index, 
		CommunicationReceiveCallType[index],
		CommunicationReceiveGridOne[index],
		CommunicationReceiveMPI_Request[index],
		CommunicationReceiveDependsOn[index]);
      }

    /* Loop over the receive handles, looking for completed (i.e. null)
       requests associated with unprocessed (i.e. non-null) grids. 
       It's insufficient to just loop over newly completed receives because
       there may be some completed receives which were not processed due
       to dependence issues. */

    for (index = 0; index < TotalReceives; index++) {
      if (CommunicationReceiveGridOne[index] != NULL &&
	  CommunicationReceiveMPI_Request[index] == MPI_REQUEST_NULL) {

	// fprintf(stdout, "::MPI:: %d %d %d %d %d\n", index, 
 	// 	CommunicationReceiveCallType[index],
 	// 	CommunicationReceiveGridOne[index],
 	// 	CommunicationReceiveMPI_Request[index],
 	// 	CommunicationReceiveDependsOn[index]);

	/* If this depends on an un-processed receive, then skip it. */

	if (CommunicationReceiveDependsOn[index] != 
	    COMMUNICATION_NO_DEPENDENCE)
	  if (CommunicationReceiveGridOne[CommunicationReceiveDependsOn[index]]
	      != NULL)
	    continue;

	grid_one = CommunicationReceiveGridOne[index];
	grid_two = CommunicationReceiveGridTwo[index];
	CommunicationReceiveIndex = index;

	/* Copy out the argument if needed */

	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  EdgeOffset[dim] = CommunicationReceiveArgument[dim][index];
	  
	/* Handle the buffers received, calling the appropriate method. */

	switch (CommunicationReceiveCallType[index]) {

	case 1:
	  errcode = grid_one->InterpolateBoundaryFromParent(grid_two);
	  break;

	case 2:
	  errcode = grid_one->CopyZonesFromGrid(grid_two, EdgeOffset);
	  break;

	case 3:
	  errcode = grid_one->DepositParticlePositions(grid_two,
			CommunicationReceiveArgument[0][index],
			CommunicationReceiveArgumentInt[0][index]);
	  break;

	case 4:
	  errcode = grid_one->CopyParentToGravitatingFieldBoundary(grid_two);
	  break;

	case 5:
	  errcode = grid_one->DepositBaryons(grid_two,
			           CommunicationReceiveArgument[0][index]);
	  break;

	case 6:
	  errcode = grid_one->AddOverlappingParticleMassField(grid_two,
							      EdgeOffset);
	  break;

	case 7:
	  errcode = grid_one->PreparePotentialField(grid_two);
	  break;

	case 8:
	  errcode = grid_one->CopyOverlappingMassField(grid_two, EdgeOffset);
	  break;

	case 9:
	  errcode = grid_one->CopyPotentialField(grid_two, EdgeOffset);
	  break;

	case 10:
	  errcode = grid_one->InterpolateAccelerations(grid_two);
	  break;
      
	case 11:  /* Note this one involves two calls. */

	  /* Project subgrid's refined fluxes to the level of this grid. */

	  if (grid_one->GetProjectedBoundaryFluxes(grid_two, 
					       SubgridFluxesRefined) == FAIL) {
	    ENZO_FAIL("Error in grid->GetProjectedBoundaryFluxes.\n");
	  }
	
	  /* Correct this grid for the refined fluxes (step #19)
	     (this also deletes the fields in SubgridFluxesRefined). */

	  // For SUBlings, the subgrid number is flagged by setting it to negative

	  igrid = CommunicationReceiveArgumentInt[0][index];
	  isubgrid = CommunicationReceiveArgumentInt[1][index];
	  SUBling = CommunicationReceiveArgumentInt[2][index];
	  if ((errcode = grid_two->CorrectForRefinedFluxes
	      (SubgridFluxesEstimate[igrid][isubgrid], &SubgridFluxesRefined, 
	       SubgridFluxesEstimate[igrid][NumberOfSubgrids[igrid] - 1],
	       SUBling, MetaData)) == FAIL) {
	    ENZO_FAIL("Error in grid->CorrectForRefinedFluxes.\n");
	  }
	  break;

	case 12:
	  errcode = grid_one->ProjectSolutionToParentGrid(*grid_two);
	  break;

	case 13:
	  errcode = grid_one->SetParticleMassFlaggingField();
	  break;

	case 14:
	  FromStart = CommunicationReceiveArgumentInt[0][index];
	  FromNumber = CommunicationReceiveArgumentInt[1][index];
	  ToStart = CommunicationReceiveArgumentInt[2][index];
	  errcode = grid_one->CommunicationSendParticles
	    (grid_two, MyProcessorNumber, FromStart, FromNumber, ToStart);
	  break;

#ifdef TRANSFER
	case 15:
	  ToNumber = CommunicationReceiveArgumentInt[0][index];
	  FromNumber = CommunicationReceiveArgumentInt[1][index];
	  PP = grid_one->ReturnPhotonPackagePointer();
	  errcode = grid_one->CommunicationSendPhotonPackages
	    (grid_two, MyProcessorNumber, ToNumber, FromNumber,
	     &PP);
	  break;
#endif /* TRANSFER */

	case 16:
	  for (dim = 0; dim < MAX_DIMENSION; dim++)
	    GridDimension[dim] = CommunicationReceiveArgumentInt[dim][index];
	  errcode = grid_one->CommunicationSendRegion
	    (grid_two, MyProcessorNumber, ALL_FIELDS, NEW_ONLY, Zero,
	     GridDimension);
	  break;

	case 17:
	  errcode = grid_one->InterpolateParticlesToGrid(NULL);
	  break;

	case 18:
	  errcode = grid_one->CommunicationSendStars(grid_two, 
						     MyProcessorNumber);
	  break;

#ifdef TRANSFER
	case 19:
	  level = CommunicationReceiveArgumentInt[0][index];
	  errcode = grid_one->SetSubgridMarkerFromParent(grid_two, level);
	  break;

	case 20:
	  errcode = grid_one->CommunicationSendSubgridMarker(grid_two, 
							     MyProcessorNumber);
	  break;
#endif

	default:
	  ENZO_VFAIL("Unrecognized call type %"ISYM"\n", 
		  CommunicationReceiveCallType[index])

	} // end: switch on call type

	/* Report error if there has been one in any of the above calls. */

	if (errcode == FAIL) {
	  ENZO_VFAIL("Error in CommunicationReceiveHandler, method %"ISYM"\n",
		  CommunicationReceiveCallType[index])

	}

	/* Mark this receive complete. */

	CommunicationReceiveGridOne[index] = NULL;
	//	MPI_Request_free(CommunicationReceiveMPI_Request+index);
	ReceivesCompletedToDate++;

      } // end: if statement to check if receive should be processed

    } // end: loop over all receives

  } // end: while loop waiting for all receives to be processed

  CommunicationReceiveIndex = 0;

  /* Reset the communication mode. */

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
  //  printf("P(%d) out of CRH\n", MyProcessorNumber);

#endif /* USE_MPI */

  return SUCCESS;

}

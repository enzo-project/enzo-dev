/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE GRIDS BY RAY TRACING WORK
/
/  written by: John Wise
/  date:       September, 2010
/  modified1:
/
/  PURPOSE:
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

#define LOAD_BALANCE_RATIO 1.05
#define MIN_LEVEL 1

/* function prototypes */

int LoadBalanceHilbertCurve(grid *GridPointers[], int NumberOfGrids, 
			    float *InputGridWork, int *NewProcessorNumber);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
double ReturnWallTime(void);
int FindField(int field, int farray[], int numfields);
void fpcol(Eflt32 *x, int n, int m, FILE *log_fptr);
void icol(int *x, int n, int m, FILE *log_fptr);

Eint32 _compare(const void * a, const void * b)
{
  return ( *(float*)a - *(float*)b );
}

int CommunicationLoadBalancePhotonGrids(HierarchyEntry **Grids[], int *NumberOfGrids,
					int FirstTimeAfterRestart)
{

  if (NumberOfProcessors == 1)
    return SUCCESS;

#ifdef TRANSFER // nee TRANSFER

  if (RadiativeTransfer == FALSE)
    return SUCCESS;

  if (RadiativeTransferLoadBalance == FALSE)
    return SUCCESS;

  double tt0, tt1;

  /* Initialize */

  float *NonZeroComputeTime, *ComputeTime;
  int *NewProcessorNumber, *NonZeroList;
  grid **NonZeroGrids;

  int i, index, index2, lvl, dim, proc, GridsMoved, TotalNumberOfGrids;
  int NumberOfBaryonFields, Nonzero;
  int FieldTypes[MAX_NUMBER_OF_BARYON_FIELDS];

  GridsMoved = 0;

  TotalNumberOfGrids = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    TotalNumberOfGrids += NumberOfGrids[lvl];

  /* Allocate memory */

  ComputeTime = new float[TotalNumberOfGrids];
  NewProcessorNumber = new int[TotalNumberOfGrids];

  for (i = 0; i < TotalNumberOfGrids; i++)
    ComputeTime[i] = 0.0;
  
  /* Compute work for each grid. */

  NumberOfBaryonFields = Grids[0][0]->GridData->
    ReturnNumberOfBaryonFields();
  Grids[0][0]->GridData->ReturnFieldType(FieldTypes);
  int RaySegNum = FindField(RaySegments, FieldTypes, NumberOfBaryonFields);

  index = 0;
  Nonzero = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    for (i = 0; i < NumberOfGrids[lvl]; i++, index++) {

      proc = Grids[lvl][i]->GridData->ReturnProcessorNumber();
      Grids[lvl][i]->GridData->SetOriginalProcessorNumber(proc);
      if (MyProcessorNumber == proc) {
//	if (FirstTimeAfterRestart)  // Possible to have no ray segment data
//	  ComputeTime[index] = Grids[lvl][i]->GridData->
//	    CountRadiationCells();
//	else
	ComputeTime[index] = Grids[lvl][i]->GridData->
	  ReturnTotalNumberOfRaySegments(RaySegNum);
	if (ComputeTime[index] > 0) Nonzero++;
      } else
	ComputeTime[index] = 0.0;

      NewProcessorNumber[index] = proc;

  } // ENDFOR grids

  /* Get total compute time over all processors */

  CommunicationSumValues(ComputeTime, TotalNumberOfGrids);
  CommunicationSumValues(&Nonzero, 1);

  if (MyProcessorNumber == ROOT_PROCESSOR) {

  /* Find grids with non-zero amount of work and only load balance
     those grids.  Construct a contiguous grid list with work. */

  NonZeroComputeTime = new float[Nonzero];
  NonZeroList = new int[Nonzero];
  for (i = 0, index = 0; i < TotalNumberOfGrids; i++)
    if (ComputeTime[i] > 0) {
      NonZeroComputeTime[index] = ComputeTime[i];
      NonZeroList[index++] = i;
    } // ENDIF

  NonZeroGrids = new grid*[Nonzero];
  index = 0;
  index2 = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    for (i = 0; i < NumberOfGrids[lvl]; i++, index++)
      if (ComputeTime[index] > 0)
	NonZeroGrids[index2++] = Grids[lvl][i]->GridData;

  /* Transfer grids from heavily-loaded processors. */

  int *NonZeroProcs = new int[Nonzero];
  LoadBalanceHilbertCurve(NonZeroGrids, Nonzero, NonZeroComputeTime,
			  NonZeroProcs);

  // Change processor number in list.
  for (i = 0; i < Nonzero; i++)
    NewProcessorNumber[NonZeroList[i]] = NonZeroProcs[i];

  delete [] NonZeroComputeTime;
  delete [] NonZeroList;
  delete [] NonZeroGrids;
  delete [] NonZeroProcs;

  } // ENDIF ROOT_PROCESSOR

#ifdef USE_MPI
  MPI_Arg Root = ROOT_PROCESSOR;
  MPI_Arg Count = TotalNumberOfGrids;
  MPI_Bcast ((void*) NewProcessorNumber, Count, IntDataType,
	     Root, MPI_COMM_WORLD);
#endif

  /* Now we know where the grids are going, transfer them. */

  index2 = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++) {

  GridsMoved = 0;
  tt0 = ReturnWallTime();

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  index = index2;
  for (i = 0; i < NumberOfGrids[lvl]; i++, index++) 
    if (Grids[lvl][i]->GridData->ReturnProcessorNumber() !=
	NewProcessorNumber[index]) {
      Grids[lvl][i]->GridData->
	CommunicationMoveGrid(NewProcessorNumber[index], FALSE, FALSE, TRUE);
      GridsMoved++;
    }

  /* Send grids */

  CommunicationDirection = COMMUNICATION_SEND;

  index = index2;
  for (i = 0; i < NumberOfGrids[lvl]; i++, index++)
    if (Grids[lvl][i]->GridData->ReturnProcessorNumber() !=
	NewProcessorNumber[index]) {
      if (RandomForcing)  //AK
	Grids[lvl][i]->GridData->AppendForcingToBaryonFields();
      Grids[lvl][i]->GridData->
	CommunicationMoveGrid(NewProcessorNumber[index], FALSE, FALSE, TRUE);
    }

  /* Receive grids */

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("CommunicationReceiveHandler() failed!\n");

  /* Update processor numbers */

  index = index2;
  for (i = 0; i < NumberOfGrids[lvl]; i++, index++) {
    Grids[lvl][i]->GridData->SetProcessorNumber(NewProcessorNumber[index]);
    if (RandomForcing)  //AK
      Grids[lvl][i]->GridData->RemoveForcingFromBaryonFields();
  }

  /* Post-process SubgridMarker to obtain grid pointers */

  index = index2;
  for (i = 0; i < NumberOfGrids[lvl]; i++, index++)
    Grids[lvl][i]->GridData->SubgridMarkerPostParallel(Grids, NumberOfGrids);

  index2 += NumberOfGrids[lvl];
  CommunicationBarrier();

  if (MyProcessorNumber == ROOT_PROCESSOR && GridsMoved > 0) {
    tt1 = ReturnWallTime();
    printf("PhotonLoadBalance[%"ISYM"]: Number of grids moved = %"ISYM" out of %"ISYM" "
	   "(%lg seconds elapsed)\n", lvl, GridsMoved, NumberOfGrids[lvl], tt1-tt0);
  }

  } // ENDFOR level

  /* Cleanup */

  delete [] ComputeTime;
  delete [] NewProcessorNumber;

#endif /* TRANSFER */

  return SUCCESS;
}

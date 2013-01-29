/***********************************************************************
/
/  GRID CLASS (DEPOSIT BARYON FIELD IN TO TARGET GRAVITATING MASS FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1999
/  modified1:  Robert Harkness
/  date:       March, 2004
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
 
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
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */
 
extern "C" void FORTRAN_NAME(dep_grid_cic)(
                               float *source, float *dest, float *temp,
			       float *velx, float *vely, float *velz,
			       float *dt, float *rfield, int *ndim,
                                   hydro_method *ihydro,
			       float *delx, float *dely, float *delz,
			       int *sdim1, int *sdim2, int *sdim3,
			       int *sstart1, int *sstart2, int *sstart3,
			       int *send1, int *send2, int *send3,
			       float *offset1, float *offset2, float *offset3,
			       int *ddim1, int *ddim2, int *ddim3,
			       int *refine1, int *refine2, int *refine3);
 
int RK2SecondStepBaryonDeposit = 0;

/* InterpolateBoundaryFromParent function */
 
int grid::DepositBaryons(grid *TargetGrid, FLOAT DepositTime)
{
 
  /* If this doesn't concern us, return. */
 
  if (TargetGrid->CommunicationMethodShouldExit(this) ||
      NumberOfBaryonFields == 0)
    return SUCCESS;
 
  TargetGrid->DebugCheck("DepositBaryons_target");
  this->DebugCheck("DepositBaryons_this");
 
  /* Declarations. */
 
  float GridStart[MAX_DIMENSION] = {0,0,0};
  int GridOffset[MAX_DIMENSION] = {0,0,0}, Refinement[MAX_DIMENSION],
      RegionDim[MAX_DIMENSION] = {1,1,1}, GridOffsetEnd[MAX_DIMENSION],
      i, j, k, index, gmindex, dim, size = 1;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* Error check: subgrid covering field must exist on entry. */
 
  if (MyProcessorNumber == ProcessorNumber &&
      BaryonField[NumberOfBaryonFields] == NULL) {
    ENZO_FAIL("subgrid covering field missing\n");
  }
 
  /* Compute refinement factors. */
 
  TargetGrid->ComputeRefinementFactors(this, Refinement);
 
  /* This routine will create a temporary patch with cell width equal to
     the target grid.  The current grid then deposits into this patch.
     Compute the TargetOffset (in grid units) and TargetStartIndex and
     the region dim (in Target units). */
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* GridOffset is the number of TargetGrid cells from the edge of
       the TargetGrid mass field and the region to be deposited by this
       grid.  It must not extended beyond the active region of TargetGrid
       (if we are depositing in self). */
 
    GridOffset[dim] = nint((GridLeftEdge[dim] -
			    TargetGrid->GravitatingMassFieldLeftEdge[dim])/

			   TargetGrid->GravitatingMassFieldCellSize) - 1; 

    if (TargetGrid == this)
      GridOffset[dim] = max(GridOffset[dim],
	nint((TargetGrid->GridLeftEdge[dim] -
	      TargetGrid->GravitatingMassFieldLeftEdge[dim])/
	     TargetGrid->GravitatingMassFieldCellSize) );
 
    /* GridStart is the distance (float) in target grid cell units between
       the exact value of GridOffset and it's integer version. */
 
    GridStart[dim] = (GridLeftEdge[dim] -
		      TargetGrid->GravitatingMassFieldLeftEdge[dim])/
                      TargetGrid->GravitatingMassFieldCellSize -
                      float(GridOffset[dim]);
 
    /* RegionDim is the size, in TargetGrid cell units, of the region (patch)
       to be deposited. It must not extend beyond the edge of the active
       region of TargetGrid. */
 
    GridOffsetEnd[dim] = nint((GridRightEdge[dim] -
			    TargetGrid->GravitatingMassFieldLeftEdge[dim])/
			   TargetGrid->GravitatingMassFieldCellSize);

    if (TargetGrid == this)
      GridOffsetEnd[dim] = min(GridOffsetEnd[dim],
	nint((TargetGrid->GridRightEdge[dim] -
	      TargetGrid->GravitatingMassFieldLeftEdge[dim])/
	     TargetGrid->GravitatingMassFieldCellSize)-1 );

    RegionDim[dim] = GridOffsetEnd[dim] - GridOffset[dim] + 1;
		
    size *= RegionDim[dim];
 
    if (TargetGrid != this && GridOffset[dim] < 0) {
      fprintf(stderr, "GridOffsetEnd[%"ISYM"] = %"ISYM" \n", dim, GridOffsetEnd[dim]);
      fprintf(stderr, "GridOffset[%"ISYM"] = %"ISYM" \n", dim, GridOffset[dim]);
      ENZO_VFAIL("GridOffset[%"ISYM"] = %"GSYM" < 0.\n", dim,GridOffset[dim])
    }
 
    if (RegionDim[dim] < 2) {
      fprintf(stderr, "GridStart[%"ISYM"] = %"ISYM" \n", dim, GridStart[dim]);
      fprintf(stderr, "GridOffsetEnd[%"ISYM"] = %"ISYM"\n", dim, GridOffsetEnd[dim]);
      fprintf(stderr, "GridOffset[%"ISYM"] = %"ISYM"\n", dim, GridOffset[dim]);
      ENZO_VFAIL("RegionDim[%"ISYM"] = %"ISYM" < 2!\n", dim, RegionDim[dim])
    }
 
  }
 
  /* Prepare the density field. */
 
  float *dens_field = NULL;
#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    dens_field = CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
#endif /* USE_MPI */
    dens_field = new float[size];
 
  if (ProcessorNumber == MyProcessorNumber) {
 
    /* Compute the dt to advance from current time to DepositTime. */
 
    FLOAT a = 1, dadt;
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(0.5*(Time+DepositTime), &a, &dadt)
	  == FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
    float dt = (DepositTime - Time)/a;
    if (HydroMethod < 3) dt = 0;
    //    dt = 0;

    /* Set up a float version of cell size to pass to fortran. */
 
    float dxfloat[MAX_DIMENSION] = {0,0,0};
    for (dim = 0; dim < GridRank; dim++)
      dxfloat[dim] = float(CellWidth[dim][0]);
 
    /* Allocate a density and velocity mesh for this grid. */
 
    float *vel_field = new float[size*4];
 
    /* Generate the density field advanced by dt using smoothed
       velocity field. */
 
//  fprintf(stderr, "Grid_DepositBaryons - call dep_grid_cic\n");

    float *input_density = BaryonField[DensNum];
    float *input_velx    = BaryonField[Vel1Num];
    float *input_vely    = BaryonField[Vel2Num];
    float *input_velz    = BaryonField[Vel3Num];
  
    float *av_dens = NULL;
    // supply zero velocity field and the time averaged density to the deposit routine for second step
    if (RK2SecondStepBaryonDeposit == 1) { 
      int current_size = 1;
      dt = 0;
      for (dim = 0; dim < GridRank; dim++)
	current_size *= GridDimension[dim];

      if (OldBaryonField[DensNum] != NULL) {
	av_dens = new float[current_size];
	for (int i=0; i<current_size; i++) 
	  av_dens[i] = (BaryonField[DensNum][i] + OldBaryonField[DensNum][i]) / 2.;
      } else 	
	av_dens = BaryonField[DensNum];

      input_density = av_dens;
    }

    //    printf("DepositBaryons, %i\n", RK2SecondStepBaryonDeposit);
    
    FORTRAN_NAME(dep_grid_cic)(input_density, dens_field, vel_field,
			       input_velx, input_vely, input_velz,
			       &dt,
			       BaryonField[NumberOfBaryonFields], &GridRank,
			       &HydroMethod,
			       dxfloat, dxfloat+1, dxfloat+2,
			       GridDimension, GridDimension+1, GridDimension+2,
			       GridStartIndex, GridStartIndex+1, GridStartIndex+2,
			       GridEndIndex, GridEndIndex+1, GridEndIndex+2,
			       GridStart, GridStart+1, GridStart+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Refinement, Refinement+1, Refinement+2);
 
    delete [] vel_field;
    if ( RK2SecondStepBaryonDeposit ) 
      if (OldBaryonField[DensNum] != NULL) 
	delete [] av_dens;
    
 
  } // end: if (ProcessorNumber == MyProcessorNumber)
 
  /* If necessary, copy data from this processor to target grid's processor.
     Note: this really needs to be put into it's own transfer routine. */
 
  if (ProcessorNumber != TargetGrid->ProcessorNumber) {
 
#ifdef USE_MPI
 
    MPI_Status status;
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg Count;
    MPI_Arg Source;

    Count = size;
    Source = ProcessorNumber;
 
    double time1 = MPI_Wtime();
 
    /* If posting a receive, then record details of call. */

    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
      CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = TargetGrid;
      CommunicationReceiveCallType[CommunicationReceiveIndex] = 5;
      CommunicationReceiveArgument[0][CommunicationReceiveIndex] = DepositTime;
    }

    /* Send Mode */

    if (MyProcessorNumber == ProcessorNumber)
      CommunicationBufferedSend(dens_field, size, DataType, 
				TargetGrid->ProcessorNumber, MPI_SENDREGION_TAG, 
				MPI_COMM_WORLD, BUFFER_IN_PLACE);

    /* Send/Recv Mode */

    if (MyProcessorNumber == TargetGrid->ProcessorNumber &&
	CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
      MPI_Recv(dens_field, size, DataType, ProcessorNumber, 
	       MPI_SENDREGION_TAG, MPI_COMM_WORLD, &status);

    /* Post receive call */

    if (MyProcessorNumber == TargetGrid->ProcessorNumber &&
	CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      MPI_Irecv(dens_field, size, DataType, ProcessorNumber, 
	        MPI_SENDREGION_TAG, MPI_COMM_WORLD, 
	        CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
      CommunicationReceiveBuffer[CommunicationReceiveIndex] = dens_field;
      CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
	CommunicationReceiveCurrentDependsOn;
      CommunicationReceiveIndex++;
    }
 
    double time3 = MPI_Wtime();
 
    CommunicationTime += time3 - time1;
 
#endif /* USE_MPI */
 
  } // end: if (ProcessorNumber != TargetGrid->ProcessorNumber)

  if (MyProcessorNumber == ProcessorNumber) {
    delete [] BaryonField[NumberOfBaryonFields];
    BaryonField[NumberOfBaryonFields] = NULL;
  }
 
  /* Return if this is not our concern. */
 
  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||

      MyProcessorNumber != TargetGrid->ProcessorNumber)
    return SUCCESS;
 
  /* Add dens_field to GravitatingMassField in target grid. */
 
  index = 0;
  for (k = 0; k < RegionDim[2]; k++)
    for (j = 0; j < RegionDim[1]; j++) {
      gmindex = (j+GridOffset[1] +
               (k+GridOffset[2])*TargetGrid->GravitatingMassFieldDimension[1])*
	      TargetGrid->GravitatingMassFieldDimension[0] + GridOffset[0];
      for (i = 0; i < RegionDim[0]; i++, gmindex++, index++)
	TargetGrid->GravitatingMassField[gmindex] += dens_field[index];
    }
 
  /* Clean up */
 
  delete [] dens_field;
 
  return SUCCESS;
 
}

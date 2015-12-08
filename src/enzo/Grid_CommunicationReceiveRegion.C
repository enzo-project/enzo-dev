/***********************************************************************
/
/  GRID CLASS (RECEIVES FROM 'FAKE' GRID TO REAL GRID)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  Robert Harkness
/  date:       January, 2004
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
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
#include "CommunicationUtilities.h"

// function prototypes
 
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */
 
 
int grid::CommunicationReceiveRegion(grid *FromGrid, int FromProcessor,
				     int SendField, int NewOrOld,
				     int RegionStart[], int RegionDim[],
				     int IncludeBoundary)
{
#ifdef USE_MPI 
  MPI_Request  RequestHandle;
  MPI_Status Status;
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
  MPI_Arg Source;
  MPI_Arg ZeroTag;

  int i, index, field, dim, Zero[] = {0, 0, 0};

  if (CommunicationShouldExit(FromProcessor, ProcessorNumber))
    return SUCCESS;

//  if (MyProcessorNumber != ProcessorNumber &&
//      MyProcessorNumber != FromProcessor)
//    return SUCCESS;

  // Compute size of region to transfer
  int NumberOfFields, RegionSize, TransferSize;
 
int SendAllBaryonFields = FALSE;

  switch( SendField ){
  case ALL_FIELDS:
  case JUST_BARYONS:
  case BARYONS_ELECTRIC:
  case BARYONS_MAGNETIC:
    SendAllBaryonFields = TRUE;
    break;
  default:
    SendAllBaryonFields = FALSE;
    break;
  }

  NumberOfFields = 0;
  if( SendAllBaryonFields == TRUE ) 
    NumberOfFields = NumberOfBaryonFields;

  //if SendField >= 0 then send only that field.
  if( SendField >= 0 )
    NumberOfFields = 1;

  if( NewOrOld == NEW_AND_OLD )
    NumberOfFields *= 2;

  if (SendField == INTERPOLATED_FIELDS) {
    switch (OutputSmoothedDarkMatter) {
    case 1: NumberOfFields = 1; break;  // density
    case 2: NumberOfFields = 5; break;  // + rms velocity + 3-velocity
    }
  }
  RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];

  TransferSize = RegionSize*NumberOfFields;

  int FromDim[MAX_DIMENSION], FromOffset[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    FromOffset[dim] = (dim < GridRank && IncludeBoundary == FALSE &&
		       SendField != INTERPOLATED_FIELDS) ?
      NumberOfGhostZones : 0;
    FromDim[dim] = RegionDim[dim] + 2*FromOffset[dim];
  }

  //MHD if for the Magnetic Field
  //MHDe is for the ElectricField
  
  int MHDRegionDim[3][3], MHDRegionSize[3]={1,1,1}, MHDFromDim[3][3];
  int MHDeRegionDim[3][3], MHDeRegionSize[3]={1,1,1}, MHDeFromDim[3][3];
  int ThisIsAFaceProjection=FALSE;
  int MHD_SendBFlag[3]={FALSE,FALSE,FALSE};
  int MHD_SendEFlag[3]={FALSE,FALSE,FALSE};
  //This is used for determining the proper Electric dimensions.
  int MHD_BoundaryOnly[3]={FALSE,FALSE,FALSE};

  //This complicated Field Specific Send business is only important in MHD_ProjectFace.
  if( UseMHDCT ){
    switch(SendField){
    case ALL_FIELDS:
      for(dim=0;dim<3;dim++){
	MHD_SendBFlag[dim]=TRUE;
	MHD_SendEFlag[dim]=TRUE;  
	MHD_BoundaryOnly[dim]=FALSE;
      }
      break;

    case BARYONS_MAGNETIC:
      for(dim=0;dim<3;dim++){
	MHD_SendBFlag[dim]=TRUE;
        MHD_SendEFlag[dim]=FALSE;
        MHD_BoundaryOnly[dim]=FALSE;
      }
      break;
    case BARYONS_ELECTRIC:
      for(dim=0;dim<3;dim++){
	MHD_SendBFlag[dim]=FALSE;
        MHD_SendEFlag[dim]=TRUE;
        MHD_BoundaryOnly[dim]=FALSE;
      }
      break;

    case JUST_BARYONS:
      for(dim=0;dim<3;dim++){
	MHD_SendBFlag[dim]=FALSE;
	MHD_SendEFlag[dim]=FALSE;
      }
      break;


    case ELECTRIC_FIELD:
      for(dim=0;dim<3;dim++){
	MHD_SendEFlag[dim] = ( MHD_ProjectThisFace[dim] == TRUE ) ? FALSE : TRUE;
	MHD_BoundaryOnly[dim] = MHD_ProjectThisFace[dim]; 

      }//dim
      break;
      
    case MAGNETIC_FIELD:

      for(dim=0;dim<3;dim++){

	if( MHD_ProjectThisFace[dim] == TRUE ){
	  MHD_BoundaryOnly[dim] = TRUE;
	  MHD_SendBFlag[dim] = TRUE;
	}else{
	  MHD_BoundaryOnly[dim] = FALSE;
	  MHD_SendBFlag[dim]=FALSE;
	}

      }//dim

      break;
    default:
      break;
    }

    // For the cell centered magnetic field:
    TransferSize += ((SendField == ALL_FIELDS)? 3*RegionSize : 0 )*
      ((NewOrOld == NEW_AND_OLD)? 2 : 1);
    
    //
    // Face Centered Magnetic Field 
    //

    for(field =0; field<3;field++){

      if( MHD_SendBFlag[field] == TRUE ){
	
	//Calculate size of transfer region
	
	for(dim=0;dim<3;dim++){
	  
	  //Only expand MHD region (in standard way) if 
	  //NOT a face projection.
	  
	  MHDRegionDim[field][dim] = RegionDim[dim]
	    +( (MHD_BoundaryOnly[dim]==TRUE)? 0: ( (field==dim) ? 1:0) );
	
	  MHDRegionSize[field] *= MHDRegionDim[field][dim];
	  MHDFromDim[field][dim] = FromDim[dim] +( (field==dim) ? 1:0 ) ;
	}//dim
      	
	//increase allocation for Face Centered field
   
	TransferSize += MHDRegionSize[field] *
	  ((NewOrOld == NEW_AND_OLD)? 2 : 1);

      }//SendBFlag
    }//field
    
    //
    // Transfer size for electric field
    //

    for(field=0;field<3;field++){

      if(MHD_SendEFlag[field]==TRUE){
	
	for( dim=0; dim<3;dim++){
	  MHDeRegionDim[field][dim] = RegionDim[dim] + 
	    ( (field==dim) ? 0 : ( (MHD_BoundaryOnly[dim]==TRUE) ? 0:1 ) );
	  MHDeRegionSize[field] *= MHDeRegionDim[field][dim];
	  MHDeFromDim[field][dim] = FromDim[dim] +( (field==dim) ? 0 : 1 );
	}
	
	TransferSize += MHDeRegionSize[field];
	
      }//Efield
    }//field
    
  }//UseMHDCT

  // Allocate buffer
 
  float *buffer = NULL;
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
    buffer = new float[TransferSize];

  if (MyProcessorNumber == FromProcessor) {
 
    index = 0;
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < FromGrid->NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS || SendAllBaryonFields == TRUE) {
	  FORTRAN_NAME(copy3d)(FromGrid->BaryonField[field], &buffer[index],
			       FromDim, FromDim+1, FromDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       FromOffset, FromOffset+1, FromOffset+2);
	  index += RegionSize;
	}
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < FromGrid->NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS || SendAllBaryonFields == TRUE ){
	  FORTRAN_NAME(copy3d)(FromGrid->OldBaryonField[field], &buffer[index],
			       FromDim, FromDim+1, FromDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       FromOffset, FromOffset+1, FromOffset+2);
	  index += RegionSize;
      }

    if (SendField == INTERPOLATED_FIELDS) {
      for (field = 0; field < NumberOfFields; field++) {
	FORTRAN_NAME(copy3d)(FromGrid->InterpolatedField[field], &buffer[index],
			     FromDim, FromDim+1, FromDim+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     FromOffset, FromOffset+1, FromOffset+2);
	index += RegionSize;
      }
    }

    if( UseMHDCT ){
      
      /* send Face B */
      if( NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY )
	for( field=0;field<3;field++)
	  if( MHD_SendBFlag[field]==TRUE){

	    FORTRAN_NAME(copy3d)(FromGrid->MagneticField[field], &buffer[index],
				 MHDFromDim[field], MHDFromDim[field]+1, MHDFromDim[field]+2,
				 MHDRegionDim[field], MHDRegionDim[field]+1, MHDRegionDim[field]+2,
				 Zero, Zero+1, Zero+2,
				 FromOffset, FromOffset+1, FromOffset+2);
	    index += MHDRegionSize[field];
	  }
      
      if( NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY )
	  for( field=0;field<3;field++)
	  if( MHD_SendBFlag[field]==TRUE){
	    FORTRAN_NAME(copy3d)(FromGrid->OldMagneticField[field], &buffer[index],
				 MHDFromDim[field], MHDFromDim[field]+1, MHDFromDim[field]+2,
				 MHDRegionDim[field], MHDRegionDim[field]+1, MHDRegionDim[field]+2,
				 Zero, Zero+1, Zero+2,
				 FromOffset, FromOffset+1, FromOffset+2);
	    index += MHDRegionSize[field];
	  }
 
      for( field=0; field<3; field++)
	if( MHD_SendEFlag[field]==TRUE){	
	  FORTRAN_NAME(copy3d)(FromGrid->ElectricField[field], &buffer[index],
		       MHDeFromDim[field], MHDeFromDim[field]+1, MHDeFromDim[field]+2,
		       MHDeRegionDim[field], MHDeRegionDim[field]+1, MHDeRegionDim[field]+2, 
		       Zero, Zero+1, Zero+2,
		       FromOffset, FromOffset+1, FromOffset+2);
	  index += MHDeRegionSize[field];
	}//field
    }//UseMHDCT

  } // ENDIF FromProcessor
 
  /* Send buffer */
 
  // Only send if processor numbers are not identical
 
  if (ProcessorNumber != FromProcessor) {
 
#ifdef MPI_INSTRUMENTATION
    starttime=MPI_Wtime();
#endif /* MPI_INSTRUMENTATION */
 
    if (MyProcessorNumber == FromProcessor) {
#ifdef MPI_INSTRUMENTATION
      if (traceMPI) 
	fprintf(tracePtr, "CRR RF: Sending %"ISYM" floats from %"ISYM" to %"ISYM"\n", 
		TransferSize, FromProcessor, ProcessorNumber);
#endif
      CommunicationBufferedSend(buffer, TransferSize, DataType, ProcessorNumber, 
				0, MPI_COMM_WORLD, BUFFER_IN_PLACE);
    } // ENDIF from processor
    
    if (MyProcessorNumber == ProcessorNumber) {

      /* Post the receive message without waiting for the message to
	 be received.  When the data arrives, this will be called again
	 in (the real) receive mode. */

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	MPI_Irecv(buffer, TransferSize, DataType, FromProcessor, 0, 
		  MPI_COMM_WORLD, 
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
	CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;      }

      /* If in send-receive mode, then wait for the message now. */

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
	MPI_Recv(buffer, TransferSize, DataType, FromProcessor, 0, 
		 MPI_COMM_WORLD, &Status);
      }

    } // ENDIF grid processor
 
#ifdef MPI_INSTRUMENTATION
    endtime=MPI_Wtime();
    timer[5] += endtime-starttime;
    counter[5] ++;
    timer[6] += double(TransferSize);
    RecvComm += endtime-starttime;
    CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
  } // ENDIF different processors
 
  /* If this is the to processor, unpack fields */
 
  int GridSize = GridDimension[0]*GridDimension[1]*GridDimension[2];
  int ActiveDims[MAX_DIMENSION], ActiveSize = 1;
  for (dim = 0; dim < GridRank; dim++) {
    ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    ActiveSize *= ActiveDims[dim];
  }
 
  if (MyProcessorNumber == ProcessorNumber &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {
 
    index = 0;
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS || SendAllBaryonFields == TRUE ){
	  if (BaryonField[field] == NULL) {
	    BaryonField[field] = new float[GridSize];
	    for (i = 0; i < GridSize; i++)
	      BaryonField[field][i] = 0;
          }
	  FORTRAN_NAME(copy3d)(&buffer[index], BaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionStart, RegionStart+1, RegionStart+2,
			       Zero, Zero+1, Zero+2);
 
	  index += RegionSize;
	}
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS || SendAllBaryonFields == TRUE) {
	  if (OldBaryonField[field] == NULL) {
	    OldBaryonField[field] = new float[GridSize];
	    for (i = 0; i < GridSize; i++)
	      BaryonField[field][i] = 0;
          }
	  FORTRAN_NAME(copy3d)(&buffer[index], OldBaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionStart, RegionStart+1, RegionStart+2,
			       Zero, Zero+1, Zero+2);
 
	  index += RegionSize;
	}

    if (SendField == INTERPOLATED_FIELDS)
      for (field = 0; field < NumberOfFields; field++) {
	if (InterpolatedField[field] == NULL) {
	  InterpolatedField[field] = new float[ActiveSize];
	  for (i = 0; i < ActiveSize; i++)
	    InterpolatedField[field][i] = 0;
	}
	FORTRAN_NAME(copy3d)(&buffer[index], InterpolatedField[field],
			     RegionDim, RegionDim+1, RegionDim+2,
			     ActiveDims, ActiveDims+1, ActiveDims+2,
			     RegionStart, RegionStart+1, RegionStart+2,
			     Zero, Zero+1, Zero+2);
 
	index += RegionSize;
      }

    if( UseMHDCT ){     
      
      /* unpack face B */
      if( NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY )
	for(field = 0; field<3; field++)
	  if( MHD_SendBFlag[field]==TRUE){
	    if(MagneticField[field] == NULL){
	      MagneticField[field] = new float[MagneticSize[field] ];
	      for(i=0;i<MagneticSize[field];i++) MagneticField[field][i] = 0.0;
	    }//allocate Bf
	    
	    FORTRAN_NAME(copy3d)(&buffer[index], MagneticField[field],
				 MHDRegionDim[field], MHDRegionDim[field]+1, MHDRegionDim[field]+2,
				 MagneticDims[field], MagneticDims[field]+1, MagneticDims[field]+2,
				 RegionStart, RegionStart+1, RegionStart+2,
				 Zero, Zero+1, Zero+2);
	    index += MHDRegionSize[field];
	  }//unpack new bf
    
      if( NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY )
	for(field = 0; field<3; field++)
	  if( MHD_SendBFlag[field]==TRUE ){
	    if(OldMagneticField[field] == NULL){
	      OldMagneticField[field] = new float[MagneticSize[field] ];
	      for(i=0;i<MagneticSize[field];i++) OldMagneticField[field][i] = 0.0;
	    }//allocate Bf
	    
	    FORTRAN_NAME(copy3d)(&buffer[index], OldMagneticField[field],
				 MHDRegionDim[field], MHDRegionDim[field]+1, MHDRegionDim[field]+2,
				 MagneticDims[field], MagneticDims[field]+1, MagneticDims[field]+2,
				 RegionStart, RegionStart+1, RegionStart+2,
				 Zero, Zero+1, Zero+2);
	    index += MHDRegionSize[field];
	  }//unpack old bf
      
      for(field=0; field<3;field++)
	if(MHD_SendEFlag[field]==TRUE){

	  if( ElectricField[field] == NULL){
	    ElectricField[field] = new float[ ElectricSize[field] ];
	    for(i=0;i<ElectricSize[field];i++)
	      ElectricField[field][i] = 0.0;
	  }

	  FORTRAN_NAME(copy3d)(&buffer[index], ElectricField[field],
		       MHDeRegionDim[field], MHDeRegionDim[field]+1, MHDeRegionDim[field]+2,
		       ElectricDims[field],ElectricDims[field]+1,ElectricDims[field]+2,
		       RegionStart, RegionStart + 1, RegionStart + 2,
		       Zero, Zero+1, Zero+2);
	  index += MHDeRegionSize[field];
	}
            
    }//UseMHDCT
 
    /* Clean up */
 
    delete [] buffer;

  } // ENDIF unpack
 
#endif /* USE_MPI */ 

  return SUCCESS;
}

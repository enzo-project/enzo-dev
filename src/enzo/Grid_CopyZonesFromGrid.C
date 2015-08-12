

/***********************************************************************
/
/  GRID CLASS (COPY OVERLAPPING ZONES FROM GRID IN ARGUMENT TO THIS GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
// This routine copies zones which overlap from the grid in the argument
//   to the current grid.  We use only the active region of the OtherGrid,
//   but copy into the entire region (including boundaries) of this grid.
//
// The input argument EdgeOffset is the amount the corner of this grid is
//   considered to have moved for grid comparison and copying purposes.
//   See Grid_CheckForOverlappingZones for more details.

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
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

extern "C" void FORTRAN_NAME(copy3drel)(float *source, float *dest,
                                   int *dim1, int *dim2, int *dim3,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dstart3);

int FindField(int field, int farray[], int numfields);

 
int grid::CopyZonesFromGrid(grid *OtherGrid, FLOAT EdgeOffset[MAX_DIMENSION])
{
 
  /* Return if this doesn't involve us. */
 

 
  if (ProcessorNumber != MyProcessorNumber &&
      OtherGrid->ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  this->DebugCheck("CopyZonesFromGrid (before)");
 
  /* declarations */
 
  int dim;

  bool shiftPos, shiftNeg; float delta; FLOAT L;

  if (ShearingBoundaryDirection!=-1){
    L=(DomainRightEdge[ShearingBoundaryDirection]-DomainLeftEdge[ShearingBoundaryDirection]);

   
    bool noMove=false;

    delta=L*AngularVelocity*VelocityGradient;
  
    //printf("L: %"GSYM" Delta: %"GSYM" %"GSYM" (%"GSYM" %"GSYM")\n", L, delta, delta, AngularVelocity, VelocityGradient);

    if (fabs(EdgeOffset[ShearingBoundaryDirection]-FLOAT(1.0)*L)<=
	CellWidth[ShearingBoundaryDirection][0]*0.1) 
      shiftPos=true;
    else shiftPos=false;
    if (fabs(EdgeOffset[ShearingBoundaryDirection]+FLOAT(1.0)*L)<=
	CellWidth[ShearingBoundaryDirection][0]*0.1) 
      shiftNeg=true;
    else shiftNeg=false;
  }

  bool isShearing= (ShearingBoundaryDirection!=-1 && (shiftNeg || shiftPos));

  /* Compute the left and right edges of this grid (including ghost zones). */
 
  FLOAT GridLeft[MAX_DIMENSION]; FLOAT GridRight[MAX_DIMENSION];

  for (dim = 0; dim < GridRank; dim++) {
    GridLeft[dim]  = CellLeftEdge[dim][0] + EdgeOffset[dim];
    GridRight[dim] = CellLeftEdge[dim][GridDimension[dim]-1] +
      CellWidth[dim][GridDimension[dim]-1]    +
      EdgeOffset[dim];
  }
 
  /* Do a quick check to see if there is any overlap. */
 
  for (dim = 0; dim < GridRank; dim++)
    if (GridLeft[dim]  >= OtherGrid->GridRightEdge[dim] ||
        GridRight[dim] <= OtherGrid->GridLeftEdge[dim]   )   
      return SUCCESS;
 

  /* There is some overlap, so copy overlapping region */
 
  //FLOAT Left, Right;
  FLOAT Left[MAX_DIMENSION];
  FLOAT Right[MAX_DIMENSION];
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
  int StartOther[MAX_DIMENSION], Dim[MAX_DIMENSION];
  int OtherDim[MAX_DIMENSION];
  
 
  /* compute start and stop indicies of overlapping region for both this
     grid and the Other grid. */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Start[dim]      = 0;
    End[dim]        = 0;
    StartOther[dim] = 0;
    OtherDim[dim]   = 1;
    Dim[dim]        = 1;
  }


 
  //bool isShearing=false;

  for (dim = 0; dim < GridRank; dim++){
    if (GridDimension[dim] > 1) {
 
      //printf("Dim %d", dim);
      /* Compute left and right positions in problem space.
	 note: include buffer zones of this grid but not the other grid. */
 
      Left[dim]  = max(GridLeft[dim], OtherGrid->GridLeftEdge[dim]);
      Right[dim] = min(GridRight[dim], OtherGrid->GridRightEdge[dim]);
 
      /* Convert this to index positions in this grid */
 
      Start[dim] = nint((Left[dim]  - GridLeft[dim]) / CellWidth[dim][0]);
      End[dim]   = nint((Right[dim] - GridLeft[dim]) / CellWidth[dim][0]) - 1;

      if (ShearingVelocityDirection==dim && isShearing){

       
	Start[dim]=(int) ceil ((Left[dim]  - GridLeft[dim]) / CellWidth[dim][0]);	
	End[dim] = (int) ceil ((Right[dim] - GridLeft[dim]) / CellWidth[dim][0]) - 1;

	if (Start[dim] >= GridDimension[dim] || End[dim] >= GridDimension [dim] ) return SUCCESS;

      }
      
      if (End[dim] - Start[dim] < 0)
	return SUCCESS;
    

      Dim[dim] = End[dim] - Start[dim] + 1;  

      /* Compute index positions in the other grid */
 
      StartOther[dim] = nint((Left[dim] - OtherGrid->CellLeftEdge[dim][0])/
			     CellWidth[dim][0]);
 
      if (isShearing && ShearingVelocityDirection==dim){
	StartOther[dim] = (int) floor((Left[dim] - OtherGrid->CellLeftEdge[dim][0]) / 
				      CellWidth[dim][0]);
      }

   

      /* Copy dimensions into temporary space */
 
      OtherDim[dim] = OtherGrid->GridDimension[dim];
    }}


  //Shearing Boundary Variables
  float rho, vx, vy, vz, v2, b2, bx, by, bz=0.0;
  int thisindex, otherindex=0;
  FLOAT a,b;  FLOAT val1=-9999;FLOAT val2=-9999;
  
  if (isShearing){
    a=(CellLeftEdge[ShearingVelocityDirection][Start[ShearingVelocityDirection]] + 
       EdgeOffset[ShearingVelocityDirection]-
       OtherGrid->CellLeftEdge[ShearingVelocityDirection][StartOther[ShearingVelocityDirection]])/
      CellWidth[ShearingVelocityDirection][Start[ShearingVelocityDirection]];
    
    b=1.0-a;
  }
  
  
  /* Calculate dimensions */
  
  // for (dim = 0; dim < MAX_DIMENSION; dim++)

 
  
  //need extra cell if you want to do shearing boundaries 
  //and that needs to be communicated from other grids possibly

  int ShearingCommunicationDims[MAX_DIMENSION];  
  
  for (dim = 0; dim <  MAX_DIMENSION; dim++){
    ShearingCommunicationDims[dim] = Dim[dim];
    if (isShearing && dim==ShearingVelocityDirection)
      ShearingCommunicationDims[dim] = Dim[dim] + 1;
  }
  
  /* If posting a receive, then record details of call. */

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE &&
      MyProcessorNumber == ProcessorNumber) {
    CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
    CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = OtherGrid;
    CommunicationReceiveCallType[CommunicationReceiveIndex] = 2;
    for (dim = 0; dim < GridRank; dim++)
      CommunicationReceiveArgument[dim][CommunicationReceiveIndex] = 
	EdgeOffset[dim];
  }
#endif /* USE_MPI */

  /* Copy data from other processor if needed (modify OtherDim and
     StartOther to reflect the fact that we are only coping part of
     the grid. */
 
  if (traceMPI) 
    fprintf(tracePtr, "CopyZones SendRegion from %"ISYM" to %"ISYM"\n", 
	    ProcessorNumber, OtherGrid->ProcessorNumber);
 
  
  if (ProcessorNumber != OtherGrid->ProcessorNumber) {
    OtherGrid->CommunicationSendRegion(OtherGrid, ProcessorNumber,
				       ALL_FIELDS, NEW_ONLY, StartOther, 
				       ShearingCommunicationDims);
    
    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	CommunicationDirection == COMMUNICATION_SEND)
      return SUCCESS;    
    
    for (dim = 0; dim < GridRank; dim++) {
      OtherDim[dim]=ShearingCommunicationDims[dim];
      StartOther[dim] = 0;
    }
  }


  /* Return if this is not our concern. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;



  // If DualEnergyFormalism is turned off, subtract v and b from Etotal to get Eint

  int iden, ivx, ivy, ivz, ietot, ieint, iBx, iBy, iBz;
  if (isShearing){

    iden=FindField(Density, FieldType, NumberOfBaryonFields);
    ivx=FindField(Velocity1, FieldType, NumberOfBaryonFields);
    ivy=FindField(Velocity2, FieldType, NumberOfBaryonFields);
    ivz=FindField(Velocity3, FieldType, NumberOfBaryonFields);
    ietot=FindField(TotalEnergy, FieldType, NumberOfBaryonFields);
    if (DualEnergyFormalism) ieint=FindField(InternalEnergy, FieldType, NumberOfBaryonFields);
    
    if (UseMHD){
      iBx=FindField(Bfield1, FieldType, NumberOfBaryonFields);
      iBy=FindField(Bfield2, FieldType, NumberOfBaryonFields);
      if (GridRank==3) iBz=FindField(Bfield3, FieldType, NumberOfBaryonFields);
      }
  }
    


  //   PrintToScreenBoundaries(BaryonField[ieint], "Eint before a copy");
  //   PrintToScreenBoundaries(BaryonField[ietot], "Etot before a copy");

  /* Copy zones */

 

  int addDim[3] = {1, OtherDim[0], OtherDim[0]*OtherDim[1]};
  int velocityTypes[3]={Velocity1, Velocity2, Velocity3};
  int Zero[3] = {0,0,0};

  if (!isShearing)
    for (int field = 0; field < NumberOfBaryonFields; field++)
      FORTRAN_NAME(copy3drel)(OtherGrid->BaryonField[field], BaryonField[field],
			      Dim, Dim+1, Dim+2,
			      OtherDim, OtherDim+1, OtherDim+2,
			      GridDimension, GridDimension+1, GridDimension+2,
			      StartOther, StartOther+1, StartOther+2,
			      Start, Start+1, Start+2);

  if (isShearing) {
    for (int field = 0; field < NumberOfBaryonFields; field++)
      for (int k = 0; k < Dim[2]; k++)
	for (int j = 0; j < Dim[1]; j++) {
	  thisindex = (0 + Start[0]) + (j + Start[1])*GridDimension[0] +
	    (k + Start[2])*GridDimension[0]*GridDimension[1];
	  otherindex = (0 + StartOther[0]) + (j + StartOther[1])*OtherDim[0] +
	    (k + StartOther[2])*OtherDim[0]*OtherDim[1];
	  for (int i = 0; i < Dim[0]; i++, thisindex++, otherindex++){

	    int otherindexB=otherindex+ addDim[ShearingVelocityDirection];
	    
	    val1=OtherGrid->BaryonField[field][otherindex];
	    val2=OtherGrid->BaryonField[field][otherindexB];
	    
	    if (DualEnergyFormalism==0 && FieldType[field]==TotalEnergy) {
	      for (int loop=0; loop<=1; loop++){
		int iLoop=otherindex;
		if (loop==1) iLoop= otherindexB;

		float vx, vy, vz, v2, rho;
		v2=0.0;
		rho= OtherGrid->BaryonField[iden][iLoop];
		vx= OtherGrid->BaryonField[ivx][iLoop];
		vy= OtherGrid->BaryonField[ivy][iLoop];  
		if (GridRank==3) vz=OtherGrid->BaryonField[ivz][iLoop];  
		else vz=0.0;
		v2=vx*vx+vy*vy+vz*vz;
		
		float bx, by, bz, b2;
		b2=0.0;
		if (UseMHD) {
		  bx= OtherGrid->BaryonField[iBx][iLoop];
		  by= OtherGrid->BaryonField[iBy][iLoop];  
		  if (GridRank==3) bz= OtherGrid->BaryonField[iBz][iLoop];  
		  else
		    bz=0.0;
		  b2=bx*bx+by*by+bz*bz;
		}
		
		if (loop==0) val1=val1- 0.5*v2 - 0.5*b2/rho;
		else if (loop==1) val2=val2- 0.5*v2 - 0.5*b2/rho;

	
	      }  
	    }
	    
	    BaryonField[field][thisindex] = (float) (b*val1+a*val2);
		     

	    if (FieldType[field]==velocityTypes[ShearingVelocityDirection]){ 
	      if (shiftNeg){
		BaryonField[field][thisindex] -=delta;	  
	      }
	      else if (shiftPos){
		BaryonField[field][thisindex] +=delta;
	      }
	    }
	  } // ENDFOR i
	} // ENDFOR j
  } // ENDIF isShearing	      

   //Update the energys due to sheared boundaries
  if (isShearing){
    
    for (int k = 0; k < Dim[2]; k++)
      for (int j = 0; j < Dim[1]; j++) {
	thisindex = (0 + Start[0]) + (j + Start[1])*GridDimension[0] +
	  (k + Start[2])*GridDimension[0]*GridDimension[1];
	for (int i = 0; i < Dim[0]; i++, thisindex++){
	  float vx, vy, vz, v2, rho, bx, by, bz, b2;
	  rho= BaryonField[iden][thisindex];  
	  vx= BaryonField[ivx][thisindex];
	  vy= BaryonField[ivy][thisindex];  
	  if (GridRank==3) vz=BaryonField[ivz][thisindex];  
	  else vz=0.0;
	  v2=vx*vx+vy*vy+vz*vz;
	  
	  b2=0.0;
	  if (UseMHD) {
	    bx= BaryonField[iBx][thisindex];
	    by= BaryonField[iBy][thisindex];  
	    if (GridRank==3) bz= BaryonField[iBz][thisindex];  
	    else
	      bz=0.0;
	    b2=bx*bx+by*by+bz*bz;
	  }
	  


	  if  (DualEnergyFormalism==0){
	   
		BaryonField[ietot][thisindex] = BaryonField[ietot][thisindex] + 0.5*v2 +0.5*b2/rho;
	  }
	  else{
	
		BaryonField[ietot][thisindex] = BaryonField[ieint][thisindex] + 0.5*v2 +0.5*b2/rho;
	  }
	  
	}}}
 
  
  int i,j,k,field;
  if( UseMHDCT ){
    
    //
    // Face Centered Magnetic Field
    //
    
    // set up some Magnetic Field properties
    // The DimToAdd business controlls whether or not to add the first active
    // face centered field.  (i.e., Bzf(z=0) on the z boundary.)
    // Only the right edge, face centered field (i.e. Bzf(z=1) on the z face) is coppied
    // for boundary calls.  This ensures correct periodicity. dcc.
    
    int MHDDim[3][3], MHDOtherDim[3][3], MHDShift[3][3], DimToAdd=0;
    
    for(field=0;field<3;field++)
      for(dim=0;dim<3;dim++){
	MHDShift[field][dim]=0;
	DimToAdd = (End[dim] == NumberOfGhostZones-1 && field == dim ) ? 0 : 1;
	MHDDim[field][dim] = Dim[dim] + ( (field == dim) ? 1 : 0 );
	MHDOtherDim[field][dim] = OtherDim[dim] + ((field == dim) ?1:0);
      }
    
    int othersize[3]={1,1,1};
    for (field =0; field<3; field++){
      
      if( MagneticField[field] == NULL )
	ENZO_VFAIL("Severe Error: Grid_CopyZonesFromGrid.  MagneticField[%d] == NULL..\n", field);
      
      othersize[field] = MHDOtherDim[field][0]*MHDOtherDim[field][1]*MHDOtherDim[field][2];
      for( k=0; k<MHDDim[field][2]; k++)
	for( j=0; j<MHDDim[field][1]; j++)
	  for( i=0; i<MHDDim[field][0]; i++){
	    thisindex = ( (i + Start[0]+MHDShift[field][0])
			  +(j+ Start[1]+MHDShift[field][1])*MagneticDims[field][0]
			  +(k+ Start[2]+MHDShift[field][2])*MagneticDims[field][1]*MagneticDims[field][0] );
	    
	    otherindex= ( (i + StartOther[0]+MHDShift[field][0] )
			  +(j+ StartOther[1]+MHDShift[field][1] )*(MHDOtherDim[field][0])
			  +(k+ StartOther[2]+MHDShift[field][2] )*(MHDOtherDim[field][0]*MHDOtherDim[field][1]));
	    
	    //MagneticField[field][thisindex] = ((dccCounter < 0 )? 0 : OtherGrid->MagneticField[field][otherindex]);
	    MagneticField[field][thisindex] = OtherGrid->MagneticField[field][otherindex];
	    }//i
    }//field
    
    
    /* Clean up if we have transfered data. */
  }//UseMHDCT

  /* Clean up if we have transfered data. */
  
  if (MyProcessorNumber != OtherGrid->ProcessorNumber)
    OtherGrid->DeleteAllFields();

 
 
  this->DebugCheck("CopyZonesFromGrid (after)");
 //  PrintToScreenBoundaries(BaryonField[ieint], "Eint after a copy");
//   PrintToScreenBoundaries(BaryonField[ietot], "Etot after a copy");
  //  printf("***Labels copy %d \n", FieldType[ivy]);

   
  return SUCCESS;
  
}
  

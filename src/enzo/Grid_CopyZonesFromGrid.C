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


int FindField(int field, int farray[], int numfields);

 
int grid::CopyZonesFromGrid(grid *OtherGrid, FLOAT EdgeOffset[MAX_DIMENSION])
{
 
  /* Return if this doesn't involve us. */
 

 
  if (ProcessorNumber != MyProcessorNumber &&
      OtherGrid->ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 

//  printf("CopyZonesFromGrid: %"ISYM"\n", NumberOfBaryonFields);
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;


 
  this->DebugCheck("CopyZonesFromGrid (before)");
 


 
  /* declarations */
 
  int dim;

  int StartSave[3];
  
  bool shiftPos, shiftNeg; float delta;
  if (ShearingBoundaryDirection!=-1){
    FLOAT L=(DomainRightEdge[ShearingBoundaryDirection]-DomainLeftEdge[ShearingBoundaryDirection]);

    int dim;
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
  }

  //PrintToScreenBoundaries(BaryonField[3], "Vz Before\n");
 
  for (dim = 0; dim < GridRank; dim++)
    if (GridDimension[dim] > 1) {
 
      /* Compute left and right positions in problem space.
	 note: include buffer zones of this grid but not the other grid. */
 
      Left[dim]  = max(GridLeft[dim], OtherGrid->GridLeftEdge[dim]);
      Right[dim] = min(GridRight[dim], OtherGrid->GridRightEdge[dim]);
 
      /* Convert this to index positions in this grid */
 
      Start[dim] = nint((Left[dim]  - GridLeft[dim]) / CellWidth[dim][0]);
      End[dim]   = nint((Right[dim] - GridLeft[dim]) / CellWidth[dim][0]) - 1;

      if (ShearingVelocityDirection==dim && isShearing){


	Start[dim] = (int) floor((Left[dim]  - GridLeft[dim]) / CellWidth[dim][0]);
	End[dim]   = (int) floor((Right[dim] - GridLeft[dim]) / CellWidth[dim][0]) - 1;
      }
 
      if (End[dim] - Start[dim] < 0)
	return SUCCESS;

    

      /* Compute index positions in the other grid */
 
      StartOther[dim] = nint((Left[dim] - OtherGrid->CellLeftEdge[dim][0])/
			     CellWidth[dim][0]);
 

      if (isShearing && ShearingVelocityDirection==dim)
	StartOther[dim] = (int) floor((Left[dim] - OtherGrid->CellLeftEdge[dim][0]) / 
				      CellWidth[dim][0]);
      

      StartSave[dim]=StartOther[dim];

      /* Copy dimensions into temporary space */
 
      OtherDim[dim] = OtherGrid->GridDimension[dim];
    }

  
 
  /* Calculate dimensions */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Dim[dim] = End[dim] - Start[dim] + 1;

  //need extra cell if you want to do shearing boundaries 
  //and that needs to be communicated from other grids possibly

  int ShearingCommunicationDims[3];  

  //We want two more cells on the other grid for the interpolation

  if (isShearing){
    int x=Dim[ShearingVelocityDirection];
    x=x+1;
    //printf("*** %d (%d)\n", x, OtherDim[ShearingVelocityDirection] );
   


    for (dim = 0; dim < MAX_DIMENSION; dim++)
      if (dim==ShearingVelocityDirection){
	ShearingCommunicationDims[dim]=x;
	//if (x>OtherDim[dim]) printf("BAAAAAAAD!!!!!\n");
      }
      else
	ShearingCommunicationDims[dim] = Dim[dim];

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
    if (isShearing)
      OtherGrid->CommunicationSendRegion(OtherGrid, ProcessorNumber,
					 ALL_FIELDS, NEW_ONLY, StartOther, 
					 ShearingCommunicationDims);
    else    
      OtherGrid->CommunicationSendRegion(OtherGrid, ProcessorNumber,
					 ALL_FIELDS, NEW_ONLY, StartOther, Dim);
    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	CommunicationDirection == COMMUNICATION_SEND)
      return SUCCESS;    
    for (dim = 0; dim < GridRank; dim++) {
      if (isShearing) {
	OtherDim[dim]=ShearingCommunicationDims[dim];
	StartOther[dim] = 0;
      }
      else{
	OtherDim[dim] = Dim[dim];
	StartOther[dim] = 0;
      }
    }
  }

  

  /* Return if this is not our concern. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;


  /* Copy zones */
 

  //Shearing Boundary Variables
  float rho, vx, vy, vz, v2, b2, bx, by, bz;
  int thisindex, otherindex;
  FLOAT a,b;  FLOAT val1, val2;
  
  if (isShearing){
    a=GridLeft[ShearingVelocityDirection]/CellWidth[ShearingVelocityDirection][0]-
      floor(GridLeft[ShearingVelocityDirection]/CellWidth[ShearingVelocityDirection][0]);
    b=1.0-a;
  }

  int addDim[3] = {1, OtherDim[0], OtherDim[0]*OtherDim[1]};
  int velocityTypes[3]={Velocity1, Velocity2, Velocity3};



  
 
  for (int field = 0; field < NumberOfBaryonFields; field++)
    for (int k = 0; k < Dim[2]; k++)
      for (int j = 0; j < Dim[1]; j++) {
	thisindex = (0 + Start[0]) + (j + Start[1])*GridDimension[0] +
                    (k + Start[2])*GridDimension[0]*GridDimension[1];
	otherindex = (0 + StartOther[0]) + (j + StartOther[1])*OtherDim[0] +
                     (k + StartOther[2])*OtherDim[0]*OtherDim[1];
	for (int i = 0; i < Dim[0]; i++, thisindex++, otherindex++){

		  
	  if (!isShearing) {  
	    BaryonField[field][thisindex] = OtherGrid->BaryonField[field][otherindex];
	  }
	  else {

	  
	    val1=OtherGrid->BaryonField[field][otherindex];//This guaranteed to be in the active zone of the other grid
	    val2=OtherGrid->BaryonField[field][otherindex+ addDim[ShearingVelocityDirection]];
	    BaryonField[field][thisindex] = (float) (b*val1+a*val2);


		     

	    if (FieldType[field]==velocityTypes[ShearingVelocityDirection]){ 
	      if (shiftNeg){
		BaryonField[field][thisindex] -=delta;	  
	      }
	      else if (shiftPos){
		BaryonField[field][thisindex] +=delta;
	      }
	      
	    }
	  }
	  
	  
	}}
 
  //Update the energys due to sheared boundaries



  if (isShearing){

    int iden=FindField(Density, FieldType, NumberOfBaryonFields);
    int ivx=FindField(Velocity1, FieldType, NumberOfBaryonFields);
    int ivy=FindField(Velocity2, FieldType, NumberOfBaryonFields);
    int ivz=FindField(Velocity3, FieldType, NumberOfBaryonFields);
    int ietot=FindField(TotalEnergy, FieldType, NumberOfBaryonFields);
    int ieint=FindField(InternalEnergy, FieldType, NumberOfBaryonFields);
    
    int iBx, iBy, iBz;
    if (useMHD){
      iBx=FindField(Velocity1, FieldType, NumberOfBaryonFields);
      iBy=FindField(Velocity2, FieldType, NumberOfBaryonFields);
      if (GridRank==3) iBz=FindField(Velocity3, FieldType, NumberOfBaryonFields);
      
    }
    
    for (int k = 0; k < Dim[2]; k++)
      for (int j = 0; j < Dim[1]; j++) {
	thisindex = (0 + Start[0]) + (j + Start[1])*GridDimension[0] +
	  (k + Start[2])*GridDimension[0]*GridDimension[1];
	for (int i = 0; i < Dim[0]; i++, thisindex++){
	  rho= BaryonField[iden][thisindex];  
	  
	    
	  vx= BaryonField[ivx][thisindex];
	  vy= BaryonField[ivy][thisindex];  
	  if (GridRank==3) vz=BaryonField[ivz][thisindex];  
	  else vz=0.0;
	  v2=vx*vx+vy*vy+vz*vz;
	  
	  b2=0.0;
	  if (useMHD) {
	    bx= BaryonField[iBx][thisindex];
	    by= BaryonField[iBy][thisindex];  
	    if (GridRank==3) bz= BaryonField[iBz][thisindex];  
	    else
	      bz=0.0;
	    b2=bx*bx+by*by+bz*bz;
	  }
	  
	  BaryonField[ietot][thisindex] = BaryonField[ieint][thisindex]  + 0.5*v2 +0.5*b2/rho;
	  
	}}}
  
  /* Clean up if we have transfered data. */
  
  if (MyProcessorNumber != OtherGrid->ProcessorNumber)
    OtherGrid->DeleteAllFields();



  this->DebugCheck("CopyZonesFromGrid (after)");
  

   
  return SUCCESS;
  
}
  

/***********************************************************************
/
/  GRID CLASS Project MHD face
/
/  written by: David Collins
/  date:       January, 2005
/  modified1:
/
/  PURPOSE:  This does projection on Electric or Magnetic fields.
/            Only acts on the face of a grid.
/            The reason this is in a seperate routine is because the BaryonFields
/            are all cell centered, so there will be no need to project to non-parents.
/
/
/  NOTE: 
/
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

/* function prototypes */
int grid::MHD_ProjectFace(grid &ParentGrid, 
				      boundary_type LeftFaceBoundaryCondition[],
				      boundary_type RightFaceBoundaryCondition[])
{

  //
  // checks that there is something do be done on this processor, with this pair of grids, now.
  //

   if( UseMHDCT != TRUE)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber && 
      MyProcessorNumber != ParentGrid.ProcessorNumber){
    return SUCCESS;
  }  

  // What this actually says: Unwrapping all this Negative Logic
  // If SEND: Only run if a.) the proceesors are different and b.) I have the Subgrid. ("this")
  // if Receive: run if the processors are the same OR I have the Parent Grid.

  if (CommunicationDirection == COMMUNICATION_SEND &&
      (MyProcessorNumber == ParentGrid.ProcessorNumber || 
       ProcessorNumber == ParentGrid.ProcessorNumber)){
    return SUCCESS;
  }

  if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE &&
      MyProcessorNumber != ParentGrid.ProcessorNumber &&
      ProcessorNumber != ParentGrid.ProcessorNumber){
    return SUCCESS;
  }
  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  //
  //Variables
  //

  int i, j, k, dim, field, One = 1, Zero = 0, skipi, skipj, skipk;
  int SkipShift=FALSE, ThisIsAFaceProjection;
  int shift[3];
  int SendField;
  int ParentStartIndex[MAX_DIMENSION], ParentDim[MAX_DIMENSION],
    ParentEndIndex[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION], Dim[MAX_DIMENSION];
  int SubgridStartIndex[MAX_DIMENSION];
  int i1, j1, k1, pindex, gindex;
  int ParentRegionDim[MAX_DIMENSION];
  float RelativeVolume = 1.0;

  FLOAT MovedSubRight[3], MovedSubLeft[3], OverlapLeft, OverlapRight;
  FLOAT Some;  

  //Initialize everything to zero
  
  //MHD_ProjectThisFace is a global variable that gets used in the Communication.
  //ThisIsAFaceProjection is a local variable for proper indexing.

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    ParentDim[dim]        = 1;
    ParentStartIndex[dim] = 0;
    ParentEndIndex[dim]   = 0;
    Dim[dim]              = 1;
    SubgridStartIndex[dim]= 0;
    MHD_ProjectThisFace[dim]     = FALSE;
    ThisIsAFaceProjection = FALSE;
  }
  
  // clearly labeld: 
  
  ParentGrid.ComputeRefinementFactors(this, Refinement);
  
  //'Some' eliminates floting-point logical errors, since comparing equality of two floats
  //is dangerous.  it's 0.1 times the minimum cell size.  (where 0.1 is arbitrary,less than 1)
  

  Some = (CellWidth[0][0] < CellWidth[1][0]) ? CellWidth[0][0]: CellWidth[1][0];
  Some = (Some < CellWidth[2][0] ) ? Some: CellWidth[2][0];
  Some *= 0.1;


  //
  //
  //Shift the subgrid by one domain size in each direction.  This is to check corner alignment
  //of parent and subgrid electric fields
  //
  //

  
  for(shift[2]=-1;shift[2]<=1;shift[2]++)
    for(shift[1]=-1;shift[1]<=1;shift[1]++)
      for(shift[0]=-1;shift[0]<=1;shift[0]++){


	SkipShift=FALSE;
	//
	//If any of the  shifted directions isn't periodic, don't do the shift.
	//If you know a slicker way to do the exit, lemme know.
	//

	for(dim=0;dim<GridRank;dim++)
	  if( LeftFaceBoundaryCondition[dim] != periodic && shift[dim] != 0) 
	    SkipShift = TRUE;

	if( SkipShift == TRUE ) 
	  continue;

	//
	// Now do the actual shifting.
	//

	for(dim=0;dim<GridRank;dim++){
	  
	  MovedSubRight[dim] = GridRightEdge[dim] 
	    + shift[dim]*(DomainRightEdge[dim]-DomainLeftEdge[dim]);
	  MovedSubLeft[dim] = GridLeftEdge[dim] 
	    + shift[dim]*(DomainRightEdge[dim]-DomainLeftEdge[dim]);
	  
	}//move 
	
	//
	// Set up index arrays.  These will be modified based on facial overlap.
	//

	for(dim=0;dim<3;dim++){

	  //This is used for indexing of the parent grid for zeroing purposes.
	  //For communication, it will be re-set later.
	  ParentDim[dim] = ParentGrid.GridDimension[dim];
	  MHD_ProjectThisFace[dim]=FALSE;	
	  OverlapLeft = max(ParentGrid.GridLeftEdge[dim], MovedSubLeft[dim] );
	  OverlapRight = min(ParentGrid.GridRightEdge[dim],MovedSubRight[dim]);

	  ParentStartIndex[dim]=
	    nint( (OverlapLeft-ParentGrid.GridLeftEdge[dim])/
		  ParentGrid.CellWidth[dim][0] ) +
	    ParentGrid.GridStartIndex[dim];

	  SubgridStartIndex[dim] = 
	    nint( (OverlapLeft-MovedSubLeft[dim])/this->CellWidth[dim][0] )+
	    this->GridStartIndex[dim];

	  Dim[dim] = nint( ( OverlapRight - OverlapLeft )/ this->CellWidth[dim][0] );

	  if( (MyProcessorNumber == ProcessorNumber) &&
	      (ProcessorNumber != ParentGrid.ProcessorNumber ) ){
	    ParentStartIndex[dim]=0;
	    ParentDim[dim]=Dim[dim]/Refinement[dim];
	  }

	  ParentEndIndex[dim] = ParentStartIndex[dim] + Dim[dim]/Refinement[dim] - 1;
	  
	}//index setup
	ThisIsAFaceProjection = FALSE;	
	//
	// Compare the actual faces for overlap and whatnot.
	//

	for(dim=0;dim<3;dim++){


	  //if the subgrids are totally disjoint, move on to the next shift.

	  if( MovedSubRight[dim] < ParentGrid.GridLeftEdge[dim] - Some ){
	    SkipShift=TRUE;
	    break;
	  }
	  if( (MovedSubLeft[dim] > ParentGrid.GridRightEdge[dim] + Some ) ){
	    SkipShift=TRUE;
	    break;
	  }
	    
	}//Disjoint

	if( SkipShift == TRUE )
	  continue;
	
	for(dim=0;dim<3;dim++){

	  //Check Left Overlap.

	  //There should be flags to only check these for certain shifts, that might save
	  //us some time.  

	  if( fabs( MovedSubLeft[dim] - ParentGrid.GridRightEdge[dim] ) < Some ){

	    //These 'face' flags indicate that only variables stored in the plane
	    //are to be projected, not the whole cell.
	    
	    MHD_ProjectThisFace[dim]=TRUE;
	    ThisIsAFaceProjection   =TRUE;
	    ParentStartIndex[dim] = ParentGrid.GridEndIndex[dim] + 1;
	    Dim[dim]=1;
	    
	    if( (MyProcessorNumber == ProcessorNumber ) &&
		(MyProcessorNumber != ParentGrid.ProcessorNumber) ){
	      ParentStartIndex[dim]=0;
	      ParentDim[dim]=1;
	    }
	    
	    ParentEndIndex[dim]= ParentStartIndex[dim];
	    
	    //Sub Left
	  }else
	    if( fabs( MovedSubRight[dim] - ParentGrid.GridLeftEdge[dim] ) < Some ){
	      
	    MHD_ProjectThisFace[dim]=TRUE;
	    ThisIsAFaceProjection   =TRUE;
	    SubgridStartIndex[dim] = GridEndIndex[dim]+1;
	    ParentStartIndex[dim] = ParentGrid.GridStartIndex[dim];

	    if( (MyProcessorNumber == ProcessorNumber) &&
		(ProcessorNumber != ParentGrid.ProcessorNumber ) ){
	      ParentStartIndex[dim]=0;
	      ParentDim[dim]=1;
	    }
	    ParentEndIndex[dim]=ParentStartIndex[dim];
	    Dim[dim]=1;
	  }

	}//Index Modification for face


	//If all three boundaries line up, then there isn't actually any overlap,
	//( just a corner, and I don't have corner centerd variables) Move on.

	if( MHD_ProjectThisFace[0]==TRUE && 
	    MHD_ProjectThisFace[1]==TRUE && 
	    MHD_ProjectThisFace[2]==TRUE )
	  continue;


	//More than one face-plane alignment means an edge.
	//We only have elctric fields on edges.
	if( MHD_ProjectB == TRUE ){
	  int numberoffaces=0;
	  for( dim=0;dim<3;dim++ )
	    if( MHD_ProjectThisFace[dim] == TRUE ) 
	      numberoffaces++;
	  if( numberoffaces > 1 )
	    continue;
	}

	//I think this skip shift here is redundant.  Whatever.
	if( SkipShift == TRUE )
	  continue;

	//
	//
	//  The Grail
	//
	//

	if( MHD_ProjectE == TRUE )
	  SendField = ELECTRIC_FIELD;
	if( MHD_ProjectB == TRUE )
	  SendField = MAGNETIC_FIELD;

	RelativeVolume = 1.0;
	for(dim=0;dim<3;dim++)
	  RelativeVolume /= float(Refinement[dim]);

	//If this is the subgrid processor, fill the parent (or dummy parent) with the 
	//projected field.
	
	if(ProcessorNumber == MyProcessorNumber){
	  if( MHD_ProjectE == TRUE ){


	    int MHDeDim[3][3], MHDeParentDim[3][3], MHDeParentSize[3]={1,1,1};
	    int MHDeParentEndIndex[3][3];
	    float RefineInv;
	    
	    for(field=0;field<3;field++){
	      
	      //Only project those fields that lie in the plane.
	      
	      if( MHD_ProjectThisFace[field] == TRUE) 
		continue;
	      
	      for(dim=0;dim<3;dim++){
		
		//This controls the extents of the actual projection.
		//note that we ONLY want to project into common faces, so the standard expansion
		//of the electric dataset (Ex(nx,ny+1,nz+1), etc)isn't valid here.
		
		MHDeDim[field][dim]=Dim[dim]+
		  ((field==dim)?0: ( MHD_ProjectThisFace[dim]==TRUE) ? 0:1);
		
		//This is only used for zeroing the parent electric field.
		MHDeParentEndIndex[field][dim] = ParentEndIndex[dim]+
		  ((field==dim)?0: ( MHD_ProjectThisFace[dim]==TRUE) ? 0:1);
		
		//This is the data structure dimension, and as such needs to follow the electric
		//dimenion convention.
		//If child and parent processor are different, there's twice 
		//as much space alocated as
		//there absolutely needs to (since we only communicate the plane in this routine)
		//If you feel really strongly, and want to code the extra logic to take that into
		//account, you MUST also change it in Grid_CommunicationReceiveRegion.

		MHDeParentDim[field][dim]=ParentDim[dim]+((field==dim)?0:1);
		MHDeParentSize[field]*=MHDeParentDim[field][dim];

	      }//dim for field dimensions
	    
	      // Allocate if necessary

	      if(ParentGrid.ProcessorNumber != MyProcessorNumber ){
		if(ParentGrid.ElectricField[field] != NULL ){
		  delete ParentGrid.ElectricField[field];
		}
		ParentGrid.ElectricField[field]=new float[ MHDeParentSize[field] ];
		for(i=0;i<MHDeParentSize[field];i++)
		  ParentGrid.ElectricField[field][i]=0.0;
	      }
	      
	      for(k=ParentStartIndex[2];k<=MHDeParentEndIndex[field][2]; k++)
		for(j=ParentStartIndex[1];j<=MHDeParentEndIndex[field][1]; j++)
		  for(i=ParentStartIndex[0];i<=MHDeParentEndIndex[field][0]; i++){
		    pindex=i+MHDeParentDim[field][0]*(j+MHDeParentDim[field][1]*k);
		    
		    ParentGrid.ElectricField[field][pindex]=0.0;
		  }
	      
	    }//field
	    
	    //Now do the actual projection
	    //Since the Parent and Subgrid electric fields are co-located along one axis,
	    //we skip the interspacing points when doing the projection.

	    for(field=0;field<3;field++){
	      
	      //Electric fields along an axis DO NOT lie in that face,
	      //so DO NOT project from the subgrid to the parent grid if this is only a 
	      //boundary overlap.
	      
	      if( MHD_ProjectThisFace[field] == TRUE ) continue;
	      
	      
	      RefineInv=1.0/Refinement[field];
	      for(k=0;k<MHDeDim[field][2];k+=((field==2)?1:Refinement[2]) ){
		k1=k/Refinement[2];
		for(j=0;j<MHDeDim[field][1];j+=((field==1)?1:Refinement[1])){
		  j1=j/Refinement[1];
		  
		  pindex= 0+ParentStartIndex[0]
		    +(j1+ParentStartIndex[1])*MHDeParentDim[field][0]
		    +(k1+ParentStartIndex[2])*MHDeParentDim[field][1]*MHDeParentDim[field][0];
		  
		  gindex = 0 + SubgridStartIndex[0]
		    +(j+SubgridStartIndex[1])*ElectricDims[field][0]
		    +(k+SubgridStartIndex[2])*ElectricDims[field][1]*ElectricDims[field][0];
		  
		  
		  //Note that we use AvgElectricField on the subgrid, but ElectricField on the 
		  //Parent.  This is because Parent.ElectricField must reflect the time structure
		  //of the subgrid advance.
		  
		  for(i=0;i<MHDeDim[field][0];i+=((field==0)?1:Refinement[0])){
		    i1=i/Refinement[0];
		    ParentGrid.ElectricField[field][pindex+i1] += 
		      AvgElectricField[field][gindex+i]*RefineInv;
		    
		  }//i
		  
		}//j
	      }//k  

	    }//field
	    //Note that Electric Fields are NOT located on their faces: the OTHER TWO fields are
	    
	  }//MHD_ProjectE

	  //This should be installed, for generality.

	  if(MHD_ProjectB == TRUE){
	    int MHDDim[3][3], MHDParentDim[3][3], MHDParentSize[3]={1,1,1},
	      MHDParentEndIndex[3][3];
	      
	      for(field=0;field<3;field++){

		//If this is a face projection, I ONLY want the field centered there.
		if( ThisIsAFaceProjection == TRUE &&
		    MHD_ProjectThisFace[field] != TRUE ) {
		  continue;
		}

		for(dim=0;dim<3;dim++){

		  //this is the variable that controlls the extents of the projection.
		  //If this is a face projection, we only need a single plane.
		  
		  MHDDim[field][dim] = Dim[dim]+
		    ( (field==dim) ? ( ( MHD_ProjectThisFace[field]==TRUE)? 0:1 ): 0 );

		  //this is only used for zeroing the parent field. 
		  MHDParentEndIndex[field][dim] = ParentEndIndex[dim]+
		    ( (field==dim) ? ( ( MHD_ProjectThisFace[field]==TRUE)? 0:1 ): 0 );
		  
		  MHDParentDim[field][dim] = ParentDim[dim]+ (( field==dim ) ? 1:0) ;
		  
		  MHDParentSize[field] *= MHDParentDim[field][dim];
		}
		
		//allocate new (if different procs)
		if( ParentGrid.ProcessorNumber != MyProcessorNumber) {
		  delete ParentGrid.MagneticField[field];
		  ParentGrid.MagneticField[field] = new float[MHDParentSize[field]];
		}
		
		//zero parent.  (always necessary)
		for (k = ParentStartIndex[2]; k <= MHDParentEndIndex[field][2]; k++)
		  for (j = ParentStartIndex[1]; j <= MHDParentEndIndex[field][1]; j++) 
		    for (i = ParentStartIndex[0]; i <= 	MHDParentEndIndex[field][0]; i++){
		      
		      pindex = i+(k*MHDParentDim[field][1] + j)*MHDParentDim[field][0];
		      ParentGrid.MagneticField[field][pindex] = 0.0;
		      
		    }
		
		skipi = skipj = skipk = 1;
		float weight = RelativeVolume*Refinement[field];

		if (field == 0) skipi = Refinement[0];
		if (field == 1) skipj = Refinement[1];
		if (field == 2) skipk = Refinement[2];
		
		for (k = 0; k < MHDDim[field][2]; k += skipk) {
		  k1 = k/Refinement[2];
		  for (j = 0; j < MHDDim[field][1]; j += skipj) {
		    j1 = j/Refinement[1];
		    
		    pindex = (0  + ParentStartIndex[0])                            + 
		      (j1 + ParentStartIndex[1])*MHDParentDim[field][0]     +
		      (k1 + ParentStartIndex[2])*MHDParentDim[field][0]*MHDParentDim[field][1];
		    
		    gindex = 0+SubgridStartIndex[0]                                      + 
		      (j+SubgridStartIndex[1])*MagneticDims[field][0]              +
		      (k+SubgridStartIndex[2])*MagneticDims[field][0]*MagneticDims[field][1];
		    
		    for (i = 0; i < MHDDim[field][0]; i += skipi) { 
		      i1 = i/Refinement[0];
		      ParentGrid.MagneticField[field][pindex+i1]
			+=MagneticField[field][gindex+i]*weight;

		    }
		  }
		}
		
	      }//field
	      
	  }//MHD_ProjectB
	
	
	}//Processor=Mine
      	
	
	/* If necessary, copy the projected field from the 'fake' ParentGrid to
	   the real one. */
	if (ProcessorNumber != ParentGrid.ProcessorNumber) {
	  for (dim = 0; dim < MAX_DIMENSION; dim++)
	    ParentRegionDim[dim] = ParentEndIndex[dim] - ParentStartIndex[dim] + 1;
	  ParentGrid.CommunicationReceiveRegion(&ParentGrid, ProcessorNumber,
			SendField, NEW_ONLY, ParentStartIndex, ParentRegionDim, TRUE);
	}
	//Clean up:
      
	//set this to false, so nothing gets contaminated.
	for(dim=0;dim<3;dim++)
	  MHD_ProjectThisFace[dim]=FALSE;
	
	
	if( MyProcessorNumber != ParentGrid.ProcessorNumber &&
	    ProcessorNumber != ParentGrid.ProcessorNumber ){
	  
	  for(field=0;field<3;field++){
	    
	    if( ParentGrid.ElectricField[field] != NULL){
	      delete ParentGrid.ElectricField[field];
	      ParentGrid.ElectricField[field]=NULL;
	    }
	    
	    if(ParentGrid.MagneticField[field] != NULL){
	      delete ParentGrid.MagneticField[field] ;
	      ParentGrid.MagneticField[field]=NULL;
	    }
	    
	  }//field
	}
      }//shift[] loops.
  
  this->DebugCheck("MHD ProjectFace, after");
  
  return SUCCESS;
  
}

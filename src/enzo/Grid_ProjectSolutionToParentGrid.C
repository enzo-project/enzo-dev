/***********************************************************************
/
/  GRID CLASS (PROJECT SOLUTION IN CURRENT GRID TO PARENT GRID
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  NOTE: This routine assumes that the parent and current grids have the
/        same baryon fields.
/
************************************************************************/

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
#include "fortran.def"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
extern "C" void FORTRAN_NAME(mult3d)(float *source, float *dest,
				     int *sdim1, int *sdim2, int *sdim3,
				     int *ddim1, int *ddim2, int *ddim3,
				     int *sstart1, int *sstart2, int *sstart3,
				     int *dstart1, int *dstart2, int *dstart3);
extern "C" void FORTRAN_NAME(div3d)(float *source, float *dest,
				    int *sdim1, int *sdim2, int *sdim3,
				    int *ddim1, int *ddim2, int *ddim3,
				    int *sstart1, int *sstart2, int *sstart3,
				    int *dstart1, int *dstart2, int *dstart3,
                                    int *rstart1, int *rstart2, int *rstart3,
                                    int *rend1, int *rend2, int *rend3);
int MakeFieldConservative(field_type field);
int grid::ProjectSolutionToParentGrid(grid &ParentGrid)
{
  /* Return if this doesn't involve us. */
 
  if (ParentGrid.CommunicationMethodShouldExit(this) ||
      NumberOfBaryonFields == 0)
    return SUCCESS;
 
  this->DebugCheck("ProjectSolutionToParentGrid (before)");
 
  /* declarations */
 
  int i, j, k, dim, field, One = 1, Zero = 0, skipi, skipj, skipk;
  int ParentStartIndex[MAX_DIMENSION], ParentDim[MAX_DIMENSION],
      ParentEndIndex[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION], Dim[MAX_DIMENSION];
 
  /* compute size of current grid fields */
 
  int ParentSize = 1, Size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    Size *= GridDimension[dim];
    ParentSize *= ParentGrid.GridDimension[dim];
    ParentDim[dim] = ParentGrid.GridDimension[dim];
    Dim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
  }
 
  /* set values that are needed for triply-nested loops. */
 
  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    ParentDim[dim]        = 1;
    ParentStartIndex[dim] = 0;
    ParentEndIndex[dim]   = 0;
    Dim[dim]              = 1;
  }
 
  /* compute refinement factor */
 
  ParentGrid.ComputeRefinementFactors(this, Refinement);
 
  /* compute the offset (in parent grid index units) from the edge of the
     parent grid to the beginning of the active region in this grid. */
 
  for (dim = 0; dim < GridRank; dim++) {
    if (GridLeftEdge[dim] >= ParentGrid.GridRightEdge[dim] ||
	GridRightEdge[dim] <= ParentGrid.GridLeftEdge[dim])
      return SUCCESS;
    ParentStartIndex[dim] = nint((GridLeftEdge[dim] -
				  ParentGrid.GridLeftEdge[dim])/
				 (ParentGrid.CellWidth[dim][0]))
                          + ParentGrid.GridStartIndex[dim];
    ParentEndIndex[dim] = ParentStartIndex[dim] + Dim[dim]/Refinement[dim] - 1;
  }
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifyPhysicalQuantities.\n");
  }
 
  /* Compute the ratio Volume[ThisGridCell]/Volume[ParentCell]. */
 
  float RelativeVolume = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    RelativeVolume /= float(Refinement[dim]);
 
  /* Multiply all fields by the density to get conserved quantities. */
 
  if (ProcessorNumber == MyProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++)
      if (MakeFieldConservative(FieldType[field]))
      FORTRAN_NAME(mult3d)(BaryonField[DensNum], BaryonField[field],
			   &Size, &One, &One, &Size, &One, &One,
			   &Zero, &Zero, &Zero, &Zero, &Zero, &Zero);
 
  /* Allocate Parent temps if it is in another processor. */
 
  if (ParentGrid.ProcessorNumber != MyProcessorNumber) {
    int ParentSize = 1;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      ParentStartIndex[dim] = 0;
      ParentDim[dim]        = Dim[dim]/Refinement[dim];
      ParentEndIndex[dim]   = ParentDim[dim] - 1;
      ParentSize *= ParentDim[dim];
    }
    for (field = 0; field < NumberOfBaryonFields; field++) {
      delete ParentGrid.BaryonField[field];
      ParentGrid.BaryonField[field] = new float[ParentSize];
    }
  }
 
  /* For each field, zero the appropriate parental zones. */
 
  int i1, j1, k1, pindex, gindex;
  if (ProcessorNumber == MyProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++){
      // Do not zero out fields that won't be projected to parent.
      if (FieldTypeNoInterpolate(FieldType[field]) == TRUE){
	continue;
      }
      for (k = ParentStartIndex[2]; k <= ParentEndIndex[2]; k++)
	for (j = ParentStartIndex[1]; j <= ParentEndIndex[1]; j++) {
	  pindex = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, pindex++)
	    ParentGrid.BaryonField[field][pindex] = 0.0;
	}
    }
  /* For each field, accumulate it's conserved quantities in the parent
     grid. */
 
  if (ProcessorNumber == MyProcessorNumber){
    for (field = 0; field < NumberOfBaryonFields; field++) {
      if (FieldTypeNoInterpolate(FieldType[field]) == TRUE)
	continue;
      skipi = skipj = skipk = 1;
      float weight = RelativeVolume;
      if (HydroMethod == Zeus_Hydro) {
	if (FieldType[field] == Velocity1) skipi = Refinement[0];
	if (FieldType[field] == Velocity2) skipj = Refinement[1];
	if (FieldType[field] == Velocity3) skipk = Refinement[2];
      }
      if (skipi*skipj*skipk != 1)
	weight *= float(skipi*skipj*skipk);
      for (k = 0; k < Dim[2]; k += skipk) {
	k1 = k/Refinement[2];
	for (j = 0; j < Dim[1]; j += skipj) {
	  j1 = j/Refinement[1];
 
	  pindex = (0  + ParentStartIndex[0])                            +
	           (j1 + ParentStartIndex[1])*ParentDim[0]               +
	           (k1 + ParentStartIndex[2])*ParentDim[0]*ParentDim[1];
 
	  gindex = 0+GridStartIndex[0]                                      +
		  (j+GridStartIndex[1])*GridDimension[0]                    +
		  (k+GridStartIndex[2])*GridDimension[0]*GridDimension[1];
 
	  for (i = 0; i < Dim[0]; i += skipi) {
	    i1 = i/Refinement[0];
	    ParentGrid.BaryonField[field][pindex+i1] +=
	      BaryonField[field][gindex+i]*weight;
	  }
	}
      }
    }

    if(UseMHDCT == TRUE ){

      if(MHD_ProjectE == TRUE ){
       
        int MHDeDim[3][3], MHDeParentDim[3][3], MHDeParentSize[3]={1,1,1};
        float RefineInv;
        for(field=0;field<3;field++){
          for(dim=0;dim<3;dim++){
            MHDeDim[field][dim]=Dim[dim]+((field==dim)?0:1);
            MHDeParentDim[field][dim]=ParentDim[dim]+((field==dim)?0:1);
            MHDeParentSize[field]*=MHDeParentDim[field][dim];
          }
	  
          if(ParentGrid.ProcessorNumber != MyProcessorNumber ){
            if(ParentGrid.ElectricField[field] != NULL ){
//             fprintf(stderr,"ProjectSolution: ElectricField not null where it should be.\n");
//             fprintf(stderr,"Find out why, and where.\n");
              delete [] ParentGrid.ElectricField[field];
            }
	    ParentGrid.ElectricField[field]=new float[ MHDeParentSize[field] ];
	    
          }
	  
          for(k=ParentStartIndex[2];k<=ParentEndIndex[2]+((field==2)?0:1); k++)
            for(j=ParentStartIndex[1];j<=ParentEndIndex[1]+((field==1)?0:1);j++)
              for(i=ParentStartIndex[0];i<=ParentEndIndex[0]+((field==0)?0:1);i++){

                pindex=i+MHDeParentDim[field][0]*(j+MHDeParentDim[field][1]*k);		
                if( pindex >= MHDeParentSize[field] )
		  ENZO_FAIL("ProjectSolutionToParentGrid: Memory Violation.\n");
                ParentGrid.ElectricField[field][pindex]=0.0;

              }

        }//field
	
	//Now do the actual projection
	//Since the Parent and Subgrid electric fields are co-located along one axis,
	//we skip the interspacing points when doing the projection.
	
	
	for(field=0;field<3;field++){
	  RefineInv=1.0/Refinement[field];
	  for(k=0;k<MHDeDim[field][2];k+=((field==2)?1:Refinement[2]) ){
	    k1=k/Refinement[2];
	    for(j=0;j<MHDeDim[field][1];j+=((field==1)?1:Refinement[1])){
	      j1=j/Refinement[1];
	
	      pindex= 0+ParentStartIndex[0]
		+(j1+ParentStartIndex[1])*MHDeParentDim[field][0]
		+(k1+ParentStartIndex[2])*MHDeParentDim[field][1]*MHDeParentDim[field][0];
	      
	      gindex = 0 + GridStartIndex[0]
		+(j+GridStartIndex[1])*ElectricDims[field][0]
		+(k+GridStartIndex[2])*ElectricDims[field][1]*ElectricDims[field][0];
                

	      //Note that we use AvgElectricField on the subgrid, but ElectricField on the 
	      //Parent.  This is because Parent.ElectricField must reflect the time structure
	      //of the subgrid advance.
              
	      for(i=0;i<MHDeDim[field][0];i+=((field==0)?1:Refinement[0])){
		i1=i/Refinement[0];
                if( pindex+i1 >= MHDeParentSize[field] )
		  ENZO_FAIL("ProjectSolutionToParentGrid: Memory Violation 2\n");
		ParentGrid.ElectricField[field][pindex+i1] += 
		  AvgElectricField[field][gindex+i]*RefineInv;
	      }//i
	      
	    }//j
	  }//k  
	}//field
        
      }//MHD_ProjectE

      if(MHD_ProjectB == TRUE){

	fprintf(stderr, "PROJB my proc %"ISYM" parent proc %"ISYM"\n", MyProcessorNumber, ParentGrid.ReturnProcessorNumber());
	int MHDDim[3][3], MHDParentDim[3][3], MHDParentSize[3]={1,1,1};
	
	for(field=0;field<3;field++){
	  for(dim=0;dim<3;dim++){
	    MHDDim[field][dim] = Dim[dim]+MHDAdd[field][dim];
	    MHDParentDim[field][dim] = ParentDim[dim]+MHDAdd[field][dim];
	    MHDParentSize[field] *= MHDParentDim[field][dim];
	  }
	  if( ParentGrid.ProcessorNumber != MyProcessorNumber) {
	    fprintf(stderr,"allocating magnetic field\n");
	    delete ParentGrid.MagneticField[field];
	    ParentGrid.MagneticField[field] = new float[MHDParentSize[field]];
	    
	  }
	  
	  
	}//field
	
	for (field = 0; field < 3; field++)
	  for (k = ParentStartIndex[2]; k <= ParentEndIndex[2]+MHDAdd[field][2]; k++)
	    for (j = ParentStartIndex[1]; j <= ParentEndIndex[1]+MHDAdd[field][1]; j++) 
	      for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]+MHDAdd[field][0]; i++){

		pindex = i+(k*MHDParentDim[field][1] + j)*MHDParentDim[field][0];
		ParentGrid.MagneticField[field][pindex] = 0.0;
	      }
	
	for (field = 0; field < 3; field++) {
	  skipi = skipj = skipk = 1;
	  float weight = RelativeVolume;
	  
	  if (field == 0) skipi = Refinement[0];
	  if (field == 1) skipj = Refinement[1];
	  if (field == 2) skipk = Refinement[2];
	  
	  weight *= float(skipi*skipj*skipk);
	  
	  for (k = 0; k < MHDDim[field][2]; k += skipk) {
	    k1 = k/Refinement[2];
	    for (j = 0; j < MHDDim[field][1]; j += skipj) {
	      j1 = j/Refinement[1];
	      
	      pindex = (0  + ParentStartIndex[0])                            + 
		(j1 + ParentStartIndex[1])*MHDParentDim[field][0]     +
		(k1 + ParentStartIndex[2])*MHDParentDim[field][0]*MHDParentDim[field][1];
	      
	      gindex = 0+GridStartIndex[0]                                      + 
		(j+GridStartIndex[1])*MagneticDims[field][0]              +
		(k+GridStartIndex[2])*MagneticDims[field][0]*MagneticDims[field][1];
	      
	      for (i = 0; i < MHDDim[field][0]; i += skipi) { 
		i1 = i/Refinement[0];
		ParentGrid.MagneticField[field][pindex+i1]
		  +=MagneticField[field][gindex+i]*weight;
		
		
	      }
	    }
	  }
	  
	}//field
      }//Proj B
      
    }//UseMHDCT
   
  } // if (ProcessorNumber == MyProcessorNumber)
    
 
  /* If necessary, copy the projected field from the 'fake' ParentGrid to
     the real one. */
 
  if (ProcessorNumber != ParentGrid.ProcessorNumber) {

    /* If posting a receive, then record details of call. */
  int FieldToSend = JUST_BARYONS; 

  if( UseMHDCT == TRUE ){

    if( MHD_ProjectB == TRUE ){
      FieldToSend = BARYONS_MAGNETIC;
    }
    if( MHD_ProjectE == TRUE ){
      FieldToSend =  BARYONS_ELECTRIC;
    }

  }

#ifdef USE_MPI
    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
      CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = &ParentGrid;
      CommunicationReceiveCallType[CommunicationReceiveIndex] = 12;
    }
#endif /* USE_MPI */

    int ParentRegionDim[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ParentRegionDim[dim] = ParentEndIndex[dim] - ParentStartIndex[dim] + 1;
    ParentGrid.CommunicationReceiveRegion(&ParentGrid, ProcessorNumber,
	      FieldToSend, NEW_ONLY, ParentStartIndex, ParentRegionDim, TRUE);

    /* Return if only posting the receive, not actually getting the data. */

    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE)
      return SUCCESS;

  }
 
  /* Divide all fields by mass to return to original quantity. */
 
  for (field = 0; field < NumberOfBaryonFields; field++)
    if ( MakeFieldConservative(FieldType[field]) ) {
      if (ProcessorNumber == MyProcessorNumber)
	FORTRAN_NAME(div3d)(BaryonField[DensNum], BaryonField[field],
			    &Size, &One, &One, &Size, &One, &One,
			    &Zero, &Zero, &Zero, &Zero, &Zero, &Zero,
			    &Zero, &Zero, &Zero, &Size, &Zero, &Zero);
      if (ParentGrid.ProcessorNumber == MyProcessorNumber)
	FORTRAN_NAME(div3d)(ParentGrid.BaryonField[DensNum],
			    ParentGrid.BaryonField[field],
			    ParentGrid.GridDimension,
			      ParentGrid.GridDimension+1,
			      ParentGrid.GridDimension+2,
			    ParentGrid.GridDimension,
			      ParentGrid.GridDimension+1,
			      ParentGrid.GridDimension+2,
			    &Zero, &Zero, &Zero, &Zero, &Zero, &Zero,
			    ParentStartIndex, ParentStartIndex+1,
                             ParentStartIndex+2,
			    ParentEndIndex, ParentEndIndex+1, ParentEndIndex+2);
			
    }
 
  /* If appropriate, restore consistency between total and internal
     energy in projected regions. */
 
  if (ParentGrid.ProcessorNumber == MyProcessorNumber)
   if (DualEnergyFormalism)
    for (k = ParentStartIndex[2]; k <= ParentEndIndex[2]; k++)
      for (j = ParentStartIndex[1]; j <= ParentEndIndex[1]; j++) {
 
	i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];
 
	for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++)
	  ParentGrid.BaryonField[TENum][i1] =
	    ParentGrid.BaryonField[GENum][i1] + 0.5*
	    ParentGrid.BaryonField[Vel1Num][i1] *
	    ParentGrid.BaryonField[Vel1Num][i1];
 
	i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];
 
	if (GridRank > 1)
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++)
	    ParentGrid.BaryonField[TENum][i1] += 0.5*
	      ParentGrid.BaryonField[Vel2Num][i1] *
	      ParentGrid.BaryonField[Vel2Num][i1];
 
	i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];
 
	if (GridRank > 2)
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++)
	    ParentGrid.BaryonField[TENum][i1] += 0.5*
	      ParentGrid.BaryonField[Vel3Num][i1] *
	      ParentGrid.BaryonField[Vel3Num][i1];

	if (HydroMethod == MHD_RK) {
	  float B2;
	  i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++) {
	    B2 = pow(ParentGrid.BaryonField[B1Num][i1],2) + 
      	         pow(ParentGrid.BaryonField[B2Num][i1],2) +
	         pow(ParentGrid.BaryonField[B3Num][i1],2);
	    ParentGrid.BaryonField[TENum][i1] += 
	      0.5 * B2 / ParentGrid.BaryonField[DensNum][i1];  
	  }
	}

        if(UseMHDCT==TRUE){
	  i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];
          for(i = ParentStartIndex[0];i <= ParentEndIndex[0];i++,i1++){
            ParentGrid.BaryonField[TENum][i1]+=0.5*(pow(ParentGrid.CenteredB[0][i1],2)+pow(ParentGrid.CenteredB[1][i1],2)
						    +pow(ParentGrid.CenteredB[2][i1],2))/ParentGrid.BaryonField[DensNum][i1];
          }
	}
		
      } // end: loop over faces
 
  /* Clean up the fake ParentGrid. */
 
  if (ParentGrid.ProcessorNumber != MyProcessorNumber)

    for (field = 0; field < NumberOfBaryonFields; field++) {
      delete ParentGrid.BaryonField[field];
      ParentGrid.BaryonField[field] = NULL;
    }
 
  ParentGrid.DebugCheck("ProjectSolutionToParentGrid (Parent, after)");
 
  return SUCCESS;
}

/***********************************************************************
/
/  GRID CLASS (INTERPOLATE BOUNDARIES FROM PARENT GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/    This function interpolates boundary values from the parent grid
/    (specified in the argument) to the current grid.  The interpolation
/    used should be monotonic and conservative (and preferably third-
/    order accurate).  The values are also interpolated linearly in time,
/    using the OldBaryonField and BaryonField of the parent.  Make sure
/    that both of these fields have intact boundary values.
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
 
int FindField(int f, int farray[], int n);

extern "C" void FORTRAN_NAME(mhd_interpolate)
  (float *parentx, float *parenty, float *parentz,
   int parentdim[], int refine[],
   float * childx, float * childy, float *childz,
   int childdim[], int childstart[], int refinedim[],
   float * otherbx, float * otherby, float * otherbz,
   int otherdim[], int otherstart[],
   float* DyBx, float* DzBx, float* DyzBx,
   float* DxBy, float* DzBy, float* DxzBy,
   float* DxBz, float* DyBz, float* DxyBz,
   int* DBxFlag,int* DByFlag,int* DBzFlag,
   FLOAT *cellsizex, FLOAT *cellsizey, FLOAT *cellsizez,
   int * face, int * step, int * counter);

extern "C" void FORTRAN_NAME(interpolate)
                             (int *rank, float *pfield, int pdim[],
			      int pis[], int pie[], int r[],
			      float *field, int dim[], int is[], float *work,
			      interpolation_type *imethod, int *posflag,
			      int *ierror);
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
extern "C" void FORTRAN_NAME(combine3d)(
               float *source1, float *weight1, float *source2, float *weight2,
	       float *dest, int *sdim1, int *sdim2, int *sdim3,
	       int *ddim1, int *ddim2, int *ddim3,
	       int *sstart1, int *sstart2, int *sstart3,
	       int *dstart1, int *dstart2, int *dstart3,
	       int *ivel_flag, int *irefine);
 
/* InterpolateBoundaryFromParent function */
int MakeFieldConservative(field_type field); 
int grid::InterpolateBoundaryFromParent(grid *ParentGrid)
{
 
  /* Return if this doesn't involve us. */
 
  if (this->CommunicationMethodShouldExit(ParentGrid))
    return SUCCESS;
 
  /* declarations */
 
  int StartIndex[MAX_DIMENSION], EndIndex[MAX_DIMENSION];
  int ParentStartIndex[MAX_DIMENSION], ParentTempEndIndex[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION], Offset[MAX_DIMENSION];
  int ZeroVector[MAX_DIMENSION], ParentTempStartIndex[MAX_DIMENSION];
  int ParentTempDim[MAX_DIMENSION], TempDim[MAX_DIMENSION];
  int ParentDim[MAX_DIMENSION];
  int ParentTempSize, WorkSize, TempSize, One = 1, Zero = 0;
  int i, j, k, dim, field, fieldindex, tempindex, interp_error;
  float *TemporaryField, *TemporaryDensityField, *Work,
        *ParentTemp[MAX_NUMBER_OF_BARYON_FIELDS], *FieldPointer;
  interpolation_type FieldInterpolationMethod;

  if (NumberOfBaryonFields > 0) {

    interp_error = FALSE;
 
    /* Compute refinement factors and set zero. */
 
    ParentGrid->ComputeRefinementFactors(this, Refinement);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ZeroVector[dim] = 0;
 
    /* Find density field */
 
    int densfield = -1;

    if(AccelerationHack != TRUE) {  //this code is also used to set the acceleration field.
      if ((densfield=FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("No density field!\n");
      }
    }
 
    /* Set up array of flags if we are using SecondOrderB interpolation
       method.  These flags indicate if the quantity must always by > 0.
       For Zeus interpolation, they are overloaded as cell vs. face flags. */

    // Which interpolation method for AccelerationHack?
 
    int SecondOrderBFlag[MAX_NUMBER_OF_BARYON_FIELDS];

    for (int i=0; i<MAX_NUMBER_OF_BARYON_FIELDS; i++) {
      SecondOrderBFlag[i] = 0;
    }

    if (InterpolationMethod == SecondOrderB)
      for (field = 0; field < NumberOfBaryonFields; field++) {
        if (FieldType[field] == TotalEnergy || FieldType[field] == Pressure ||
            FieldType[field] == InternalEnergy)
          SecondOrderBFlag[field] = 2;   // enforce monotonicity
        else
          SecondOrderBFlag[field] = 2;   // enforce only positivity
        if (FieldType[field] >= Velocity1 && FieldType[field] <= Velocity3)
          SecondOrderBFlag[field] = 2;   //  no positivity for velocity
      }
    if (HydroMethod == Zeus_Hydro)
      for (field = 0; field < NumberOfBaryonFields; field++)
        if ( (FieldType[field] >= Velocity1 && FieldType[field] <= Velocity3 ) ||
             (FieldType[field] >= Acceleration0 && FieldType[field] <= Acceleration2 ) ){
          SecondOrderBFlag[field] = FieldType[field] - Velocity1 + 1;
          }else{
          SecondOrderBFlag[field] = 0;
          }

 
    /* Compute coefficient factors for linear interpolation in time.
       Note: interp = coef1*Old + coef2*New. */
 
    float coef1 = 0, coef2 = 1;
    if (Time != ParentGrid->Time) {
      if (ParentGrid->Time <= ParentGrid->OldTime) {
	ENZO_FAIL("ParentGrid fields are at the same time or worse.\n");
      }
      coef1 = max((ParentGrid->Time -                Time)/
                  (ParentGrid->Time - ParentGrid->OldTime), 0.0);
      coef2  = (1.0 - coef1);
    }
 
    /* Compute the start and end indicies (in this grid units) of this grid
       within it's parent */
    /* StartIndex = cells from left edge of parent active region
	                      to left edge of grid total region
		  + boundary cells of parent grid (current grid units). */
    /*   EndIndex = cells from left  edge of parent active region
	                      to right edge of grid total region
		  + boundary cells of parent grid (current grid units). */
    /* Ee adjust StartIndex and EndIndex to start and end at a parental
       cell edge - this means interpolating a larger region than is
       necessary. */
 
    ParentTempSize = 1;
    TempSize       = 1;
    WorkSize       = 1;
    for (dim = 0; dim < GridRank; dim++) {
 
      StartIndex[dim] = nint((CellLeftEdge[dim][0] -
			      ParentGrid->CellLeftEdge[dim][0])/
			     CellWidth[dim][0]);
      EndIndex[dim] = nint((CellLeftEdge[dim][GridDimension[dim]-1] +
			    CellWidth[dim][GridDimension[dim]-1]
			    - ParentGrid->CellLeftEdge[dim][0])/
			   CellWidth[dim][0])
		      - 1;
 
      /* Record the offset between the real StartIndex and the one adjusted
	 to conform with the parent. */
 
      Offset[dim] = StartIndex[dim] -
	            int(StartIndex[dim]/Refinement[dim])*Refinement[dim];
 
      /* Adjust Start/end index to conform with parental cell edges. */
 
      StartIndex[dim] = int(StartIndex[dim]/Refinement[dim])*Refinement[dim];
      EndIndex[dim] = (int(EndIndex[dim]/Refinement[dim])+1)*Refinement[dim]-1;
 
      /* Compute values for the ParentTemp fields. */
 
      ParentStartIndex[dim]     = StartIndex[dim] / Refinement[dim];
      ParentTempDim[dim]        = (EndIndex[dim] - StartIndex[dim] + 1) /
	                          Refinement[dim];
      ParentDim[dim]            = ParentGrid->GridDimension[dim];
 
      /* Set the start and and end indicies to pass to interpolate. */
 
      ParentTempStartIndex[dim] = Refinement[dim];
      ParentTempEndIndex[dim]   = Refinement[dim]*(ParentTempDim[dim]+1) - 1;
 
      /* Add an extra cell to each side of the ParentTemp field since the
         interpolator will require it (and error check). */
 
      ParentStartIndex[dim]     -= 1;
      ParentTempDim[dim]        += 2;
      if (ParentStartIndex[dim] < 0 ||
          ParentStartIndex[dim]+ParentTempDim[dim] >
          ParentGrid->GridDimension[dim]) {
        ENZO_VFAIL("Parent grid not big enough for interpolation!  ParentStartIndex[%"ISYM"] = %"ISYM"  ParentTempDim = %"ISYM"ParentGrid->GridDimension = %"ISYM"\n",dim, ParentStartIndex[dim], ParentTempDim[dim], ParentGrid->GridDimension[dim])
      }
 
      /* Compute the dimensions of the current grid temporary field. */
 
      TempDim[dim]            = EndIndex[dim] - StartIndex[dim] + 1;
 
      /* Compute size of current grid (in floats) for the temp fields. */
 
      ParentTempSize *= ParentTempDim[dim];
      TempSize       *= TempDim[dim];
      WorkSize       *= (TempDim[dim]/Refinement[dim] + 1);
    }
 
    /* Fill out the rest if dim < MAX_DIMENSION. */
 
    for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
      ParentTempDim[dim]      = 1;
      ParentStartIndex[dim]   = 0;
      ParentTempEndIndex[dim] = 0;
      ParentDim[dim]          = 1;
      TempDim[dim]            = 1;
      Offset[dim]             = 0;
    }
 
    /* If posting a receive, then record details of call. */

#ifdef USE_MPI
    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
      CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = ParentGrid;
      CommunicationReceiveCallType[CommunicationReceiveIndex] = 1;
    }
#endif /* USE_MPI */

    /* Copy data from other processor if needed (modify ParentDim and
       ParentStartIndex to reflect the fact that we are only coping part of
       the grid. */
 
    if (ProcessorNumber != ParentGrid->ProcessorNumber) {
      int NewOrOld = (Time == ParentGrid->Time) ? NEW_ONLY : NEW_AND_OLD;
      ParentGrid->CommunicationSendRegion(ParentGrid, ProcessorNumber,
		   ALL_FIELDS, NewOrOld, ParentStartIndex, ParentTempDim);
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	  CommunicationDirection == COMMUNICATION_SEND)
	return SUCCESS;
      for (dim = 0; dim < GridRank; dim++) {
	ParentDim[dim] = ParentTempDim[dim];
	ParentStartIndex[dim] = 0;
      }
    }
 
    /* Return if this is not our concern. */
 
    if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;
 
    /* Allocate temporary space. */
 
    TemporaryField        = new float[TempSize]();
    TemporaryDensityField = new float[TempSize]();
    Work                  = new float[WorkSize]();
    for (field = 0; field < NumberOfBaryonFields; field++)
      ParentTemp[field]     = new float[ParentTempSize]();
 
    /* Copy just the required section from the parent fields to the temp
       space, doing the linear interpolation in time as we do it. */
 
    int VelocityShiftFlag;

    for (field = 0; field < NumberOfBaryonFields; field++) {
      VelocityShiftFlag = 0;

/*      if (HydroMethod == Zeus_Hydro &&
	  FieldType[field] >= Velocity1 && FieldType[field] <= Velocity3)
	  VelocityShiftFlag = FieldType[field] - Velocity1 + 1; */

      float *ParentOld = ParentGrid->OldBaryonField[field];
      if (Time == ParentGrid->Time)
	ParentOld = ParentGrid->BaryonField[field]; // not used

      if (ParentOld != NULL && ParentGrid->BaryonField[field] != NULL)
	FORTRAN_NAME(combine3d)(
				ParentOld, &coef1, ParentGrid->BaryonField[field], &coef2,
				ParentTemp[field], ParentDim, ParentDim+1, ParentDim+2,
				ParentTempDim, ParentTempDim+1, ParentTempDim+2,
				&Zero, &Zero, &Zero,
				ParentStartIndex, ParentStartIndex+1, ParentStartIndex+2,
				&VelocityShiftFlag, Refinement);
    }
 
    /* Multiply ParentTemp fields by their own density to get conserved
       quantities. */
 
    if (ConservativeInterpolation)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (MakeFieldConservative(FieldType[field])) {
	  FORTRAN_NAME(mult3d)(ParentTemp[densfield], ParentTemp[field],
			       &ParentTempSize, &One, &One,
			       &ParentTempSize, &One, &One,
			       &Zero, &Zero, &Zero, &Zero, &Zero, &Zero);
	}
    
    /* Do the interpolation for the density field. */
 
    if (HydroMethod == Zeus_Hydro)
      InterpolationMethod = (SecondOrderBFlag[densfield] == 0) ?
	SecondOrderA : SecondOrderC;

    if( AccelerationHack != TRUE ) {
      FORTRAN_NAME(interpolate)(&GridRank,
			      ParentTemp[densfield], ParentTempDim,
			      ParentTempStartIndex, ParentTempEndIndex,
                                 Refinement,
			      TemporaryDensityField, TempDim, ZeroVector, Work,
			      &InterpolationMethod,
			      &SecondOrderBFlag[densfield], &interp_error);
      if (interp_error) {
	printf("P%"ISYM": Error interpolating density.\n"
		   "ParentGrid ID = %"ISYM"\n"
		   "\t LeftEdge  = %"PSYM" %"PSYM" %"PSYM"\n"
		   "\t RightEdge = %"PSYM" %"PSYM" %"PSYM"\n"
		   "ThisGrid ID = %"ISYM"\n"
		   "\t LeftEdge  = %"PSYM" %"PSYM" %"PSYM"\n"
		   "\t RightEdge = %"PSYM" %"PSYM" %"PSYM"\n",
		   MyProcessorNumber, ParentGrid->ID, 
		   ParentGrid->GridLeftEdge[0], ParentGrid->GridLeftEdge[1], 
		   ParentGrid->GridLeftEdge[2], ParentGrid->GridRightEdge[0], 
		   ParentGrid->GridRightEdge[1], ParentGrid->GridRightEdge[2],
		   this->ID, 
		   this->GridLeftEdge[0], this->GridLeftEdge[1], 
		   this->GridLeftEdge[2], this->GridRightEdge[0], 
	       this->GridRightEdge[1], this->GridRightEdge[2]);
	ENZO_FAIL("");
      }
    } // ENDIF !AccelerationHack

    /* Loop over all the fields. */
 
    for (field = 0; field < NumberOfBaryonFields; field++) {
 
      if (HydroMethod == Zeus_Hydro)
        InterpolationMethod = (SecondOrderBFlag[field] == 0) ?
            SecondOrderA : SecondOrderC;
      
      // Set FieldInterpolationMethod to be FirstOrderA for 
      // fields that shouldn't be interpolated.'
      FieldInterpolationMethod = InterpolationMethod;
      if (FieldTypeNoInterpolate(FieldType[field]) == TRUE)
        FieldInterpolationMethod = FirstOrderA; 
 
      /* Interpolating from the ParentTemp field to a Temporary field.  This
	 is done for the entire current grid, not just it's boundaries.
	 (skip density since we did it already) */

      if (FieldType[field] != Density && FieldType[field] != DebugField) {
	//      if (FieldType[field] != Density) {
	FORTRAN_NAME(interpolate)(&GridRank,
				  ParentTemp[field], ParentTempDim,
				  ParentTempStartIndex, ParentTempEndIndex,
                                     Refinement,
				  TemporaryField, TempDim, ZeroVector, Work,
				  &FieldInterpolationMethod,
				  &SecondOrderBFlag[field], &interp_error);
	if (interp_error) {
	  printf("P%"ISYM": Error interpolating field %"ISYM" (%s).\n"
		     "ParentGrid ID = %"ISYM"\n"
		     "\t LeftEdge  = %"PSYM" %"PSYM" %"PSYM"\n"
		     "\t RightEdge = %"PSYM" %"PSYM" %"PSYM"\n"
		     "ThisGrid ID = %"ISYM"\n"
		     "\t LeftEdge  = %"PSYM" %"PSYM" %"PSYM"\n"
		 "\t RightEdge = %"PSYM" %"PSYM" %"PSYM"\n",
		     MyProcessorNumber, field, DataLabel[field], ParentGrid->ID, 
		     ParentGrid->GridLeftEdge[0], ParentGrid->GridLeftEdge[1], 
		     ParentGrid->GridLeftEdge[2], ParentGrid->GridRightEdge[0], 
		     ParentGrid->GridRightEdge[1], ParentGrid->GridRightEdge[2],
		     this->ID, 
		     this->GridLeftEdge[0], this->GridLeftEdge[1], 
		     this->GridLeftEdge[2], this->GridRightEdge[0], 
		 this->GridRightEdge[1], this->GridRightEdge[2]);
	  ENZO_FAIL("");
	}
      }
 
      /* Divide by density field to convert from conserved to physical
         variables (skipping density). */
 
      if (ConservativeInterpolation)
	if (MakeFieldConservative(FieldType[field]))
	  FORTRAN_NAME(div3d)(TemporaryDensityField, TemporaryField,
			      &TempSize, &One, &One,
			      &TempSize, &One, &One,
			      &Zero, &Zero, &Zero, &Zero, &Zero, &Zero,
			      &Zero, &Zero, &Zero, &TempSize, &Zero, &Zero);
 
      /* Set FieldPointer to either the correct field (density or the one we
	 just interpolated to). */
 
      if (FieldType[field] == Density)
	FieldPointer = TemporaryDensityField;
      else 
	FieldPointer = TemporaryField;
 
      /* Copy needed portion of temp field to current grid. */
 
      /* a) j/k slices. */
 
      for (k = 0; k < GridDimension[2]; k++)
	for (j = 0; j < GridDimension[1]; j++) {
	  tempindex = ((k + Offset[2])*TempDim[1] + (j + Offset[1]))*TempDim[0]
	            +  (0 + Offset[0]);
	  fieldindex = (k*GridDimension[1] + j)*GridDimension[0];
	  for (i = 0; i < GridStartIndex[0]; i++)
	    BaryonField[field][fieldindex+i] = FieldPointer[tempindex+i];
	  for (i = GridEndIndex[0]+1; i < GridDimension[0]; i++)
	    BaryonField[field][fieldindex+i] = FieldPointer[tempindex+i];
	}
 
      /* b) k/i slices. */
 
      for (j = 0; j < GridStartIndex[1]; j++)
	for (k = 0; k < GridDimension[2]; k++) {
	  tempindex = ((k + Offset[2])*TempDim[1] + (j + Offset[1]))*TempDim[0]
	            +  (0 + Offset[0]);
	  fieldindex = (k*GridDimension[1] + j)*GridDimension[0];
	  for (i = 0; i < GridDimension[0]; i++, fieldindex++, tempindex++)
	    BaryonField[field][fieldindex] = FieldPointer[tempindex];
	}
      for (j = GridEndIndex[1]+1; j < GridDimension[1]; j++)
	for (k = 0; k < GridDimension[2]; k++) {
	  tempindex = ((k + Offset[2])*TempDim[1] + (j + Offset[1]))*TempDim[0]
	            +  (0 + Offset[0]);
	  fieldindex = (k*GridDimension[1] + j)*GridDimension[0];
	  for (i = 0; i < GridDimension[0]; i++, fieldindex++, tempindex++)
	    BaryonField[field][fieldindex] = FieldPointer[tempindex];
	}
 
      /* c) i/j slices. */
 
      for (k = 0; k < GridStartIndex[2]; k++)
	for (j = 0; j < GridDimension[1]; j++) {
	  tempindex = ((k + Offset[2])*TempDim[1] + (j + Offset[1]))*TempDim[0]
	            +  (0 + Offset[0]);
	  fieldindex = (k*GridDimension[1] + j)*GridDimension[0];
	  for (i = 0; i < GridDimension[0]; i++, fieldindex++, tempindex++)
	    BaryonField[field][fieldindex] = FieldPointer[tempindex];
	}
      for (k = GridEndIndex[2]+1; k < GridDimension[2]; k++)
	for (j = 0; j < GridDimension[1]; j++) {
	  tempindex = ((k + Offset[2])*TempDim[1] + (j + Offset[1]))*TempDim[0]
	            +  (0 + Offset[0]);
	  fieldindex = (k*GridDimension[1] + j)*GridDimension[0];
	  for (i = 0; i < GridDimension[0]; i++, fieldindex++, tempindex++)
	    BaryonField[field][fieldindex] = FieldPointer[tempindex];
	}

    } // end loop over fields
  
    delete [] Work;
    delete [] TemporaryField;
    delete [] TemporaryDensityField;
    for (field = 0; field < NumberOfBaryonFields; field++)
      delete [] ParentTemp[field];
 
    /* If using the dual energy formalism, then modify the total energy field
       to maintain consistency between the total and internal energy fields.
       This is necessary because the interpolation introduces small
       descrepancies between the two fields which are normally kept in sync. */

#ifdef SAB
    if (AccelerationHack != TRUE)
#endif 
    if (DualEnergyFormalism)
      if (this->RestoreEnergyConsistency(ONLY_BOUNDARY) == FAIL) {
	ENZO_FAIL("Error in grid->RestoreEnergyConsisitency.\n");
      }


     if(UseMHDCT == TRUE) {
      
       float  *MHDChildTemp[3], *dummy = new float;;
      
       int MHDParentTempDims[3][3], MHDChildTempDims[3][3], MHDParentDims[3][3];

       int MHDParentTempSize[3]={1,1,1}, MHDChildTempSize[3]={1,1,1}, One[3] = {1,1,1};
      
       *dummy = 1;

       //These get saved for the call to CID
       ParentDx =ParentGrid->CellWidth[0][0];
       ParentDy =ParentGrid->CellWidth[1][0];
       ParentDz =ParentGrid->CellWidth[2][0];
      
      
       //
       // the field loop
       //
      
       for( field=0;field<3;field++ ){
	
	 for(dim=0;dim<3;dim++){

	   //Also saved for CID
	   MHDParentTempPermanent[dim] = ParentTempDim[dim];
	   MHDRefinementFactors[dim] = Refinement[dim];

	   //Dimension of parent.  If it's off processor, it's only the dimension of the relevant 
	   //region.
	   MHDParentDims[field][dim]=ParentDim[dim]+MHDAdd[field][dim]; 
	  
	   //Dimension of relevant region.
	   MHDParentTempDims[field][dim] = ParentTempDim[dim]+MHDAdd[field][dim];

	  
	   MHDChildTempDims[field][dim] = TempDim[dim] + MHDAdd[field][dim];
	   MHDParentTempSize[field] *= MHDParentTempDims[field][dim];
	   MHDChildTempSize[field] *= MHDChildTempDims[field][dim];
	  
	 }
	
	 if(MHDParentTemp[field] != NULL ) 
	   delete MHDParentTemp[field];

	 MHDParentTemp[field] = new float[ MHDParentTempSize[field] ];
	 MHDChildTemp[field] = new float[ MHDChildTempSize[field] ];
	 
	 for(i=0;i<MHDParentTempSize[field];i++) MHDParentTemp[field][i] = -10000000.0;
	
	 float *ParentOld = ParentGrid->OldMagneticField[field];

	 if (Time == ParentGrid->Time)
	   ParentOld = ParentGrid->MagneticField[field];


	 //Linear combination in time of the parent grid.  


	 FORTRAN_NAME(combine3d)(
				 ParentOld, &coef1, ParentGrid->MagneticField[field], &coef2,
				 MHDParentTemp[field], 
				 MHDParentDims[field],
				 MHDParentDims[field]+1,
				 MHDParentDims[field]+2,
				 MHDParentTempDims[field],
				 MHDParentTempDims[field]+1,
				 MHDParentTempDims[field]+2,
				 &Zero,&Zero,&Zero,
				 ParentStartIndex, ParentStartIndex+1,
				 ParentStartIndex+2, &Zero, &Zero);

	 if(ParentGrid->ProcessorNumber != ProcessorNumber && 0==1)
	   fprintf(stderr,"Parent StartIndex %"ISYM" %"ISYM" %"ISYM" \n"
		   ,ParentStartIndex[0],ParentStartIndex[1],ParentStartIndex[2]);
	
	 if(ParentGrid->ProcessorNumber != ProcessorNumber && 1==0)
	   for(k=0;k<MHDParentTempDims[field][2];k++)
	     for(j=0;j<MHDParentTempDims[field][1];j++)
	       for(i=0;i<MHDParentTempDims[field][0];i++){
		 int index =  i + MHDParentTempDims[field][0]*(j+MHDParentTempDims[field][1]*k);
		 fprintf(stderr,"MHD: (%"ISYM",%"ISYM",%"ISYM",%"ISYM") %"FSYM"*%"FSYM" + %"FSYM"*%"FSYM" = %"FSYM"\n",
			 field,i,j,k,coef1,ParentOld[index],coef2,
			 ParentGrid->MagneticField[field][index], 
			 MHDParentTemp[field][index]);
	       }//i,j,k
	
      } //  for( field=0;field<3;field++ ) 

       
       MHD_ProlongAllocate(TempDim);
       int Step = 1;
       int mhd_interpolate_face = 0;

       //FILL
       FORTRAN_NAME(mhd_interpolate)(MHDParentTemp[0], MHDParentTemp[1], MHDParentTemp[2],
				     ParentTempDim,    Refinement,
				     MHDChildTemp[0], MHDChildTemp[1], MHDChildTemp[2],
				     TempDim, ZeroVector, TempDim,
				     dummy, dummy, dummy, ZeroVector, ZeroVector,
				     DyBx, DzBx, DyzBx,
				     DxBy, DzBy, DxzBy,
				     DxBz, DyBz, DxyBz,
				     DBxFlag,DByFlag,DBzFlag,
				     ParentGrid->CellWidth[0],
				     ParentGrid->CellWidth[1],
				     ParentGrid->CellWidth[2],
				     &mhd_interpolate_face, &Step, &Step);
       
      Step = 2;

      FORTRAN_NAME(mhd_interpolate)(MHDParentTemp[0], MHDParentTemp[1], MHDParentTemp[2],
                                    ParentTempDim,    Refinement,
                                    MHDChildTemp[0], MHDChildTemp[1], MHDChildTemp[2],
                                    TempDim, ZeroVector, TempDim,
                                    dummy, dummy, dummy, ZeroVector, ZeroVector,
				    DyBx, DzBx, DyzBx,
				    DxBy, DzBy, DxzBy,
				    DxBz, DyBz, DxyBz,
                                    DBxFlag,DByFlag,DBzFlag,
                                    ParentGrid->CellWidth[0],
                                    ParentGrid->CellWidth[1],
                                    ParentGrid->CellWidth[2],
				    &mhd_interpolate_face, &Step, &Step);

      for(field=0; field<3; field++){
	for(k=0;k<MHDChildTempDims[field][2];k++)
	  for(j=0;j<MHDChildTempDims[field][1];j++)
	    for(i=0;i<MHDChildTempDims[field][0];i++){
	      
	      tempindex = ((i )
			   +(j)*MHDChildTempDims[field][0]
			   +(k)*MHDChildTempDims[field][0]*MHDChildTempDims[field][1]);
	      
	      if(  MHDChildTemp[field][tempindex] !=  MHDChildTemp[field][tempindex] ){

		fprintf(stderr,"Error: Bad Child Temp. %"ISYM" (%"ISYM",%"ISYM",%"ISYM")\n",
			field,i,j,k);
	      }

	    }//i,j,k
      
      }//  for( field=0;field<3;field++ )


      for(field=0; field<3; field++){
	
	// i faces
	for( k=0; k<MagneticDims[field][2]; k++)
	  for(j=0; j<MagneticDims[field][1]; j++){
	    tempindex = ((0 + Offset[0])
			 +(j+ Offset[1])*MHDChildTempDims[field][0]
			 +(k+ Offset[2])*MHDChildTempDims[field][0]*MHDChildTempDims[field][1]);
	    fieldindex = (k*MagneticDims[field][1] + j)*MagneticDims[field][0];
	    for( i=0; i<MHDStartIndex[field][0]; i++){
	      MagneticField[field][fieldindex+i] = MHDChildTemp[field][tempindex+i];
	      
	      if(fieldindex >= MagneticSize[field] || tempindex >= MHDChildTempSize[field] ) {
              ENZO_FAIL("InterpolateBoundaryFromParent Bounds Check.\n");
	      }
	    }
	    for( i=MHDEndIndex[field][0]+1;i<MagneticDims[field][0]; i++){
	      MagneticField[field][fieldindex+i] = MHDChildTemp[field][tempindex+i];

	      if(fieldindex >= MagneticSize[field] || tempindex >= MHDChildTempSize[field] ) {
		ENZO_FAIL("InterpolateBoundaryFromParent Bounds Check.\n");
	      }
	      
	    }
	    
	  }//i faces
	
	// Left j face
	
	for( k=0; k<MagneticDims[field][2]; k++)
	  for(j=0; j<MHDStartIndex[field][1]; j++)
	    for(i=0; i<MagneticDims[field][0]; i++){
	      tempindex = ((i + Offset[0])
			   +(j+ Offset[1])*MHDChildTempDims[field][0]
			   +(k+ Offset[2])*MHDChildTempDims[field][0]*MHDChildTempDims[field][1]);
	      fieldindex = (k*MagneticDims[field][1] + j)*MagneticDims[field][0] +i;
	      MagneticField[field][fieldindex] = MHDChildTemp[field][tempindex];

	      if(fieldindex >= MagneticSize[field] || tempindex >= MHDChildTempSize[field] ) {
		ENZO_FAIL("InterpolateBoundaryFromParent Bounds Check.\n");
	      }
	      
	    }//left j face (i loop)
	
	// Right j face
	for( k=0; k<MagneticDims[field][2]; k++)
	  for(j=MHDEndIndex[field][1]+1; j<MagneticDims[field][1]; j++)
	    for(i=0; i<MagneticDims[field][0]; i++){
	      tempindex = ((i + Offset[0])
			   +(j+ Offset[1])*MHDChildTempDims[field][0]
			   +(k+ Offset[2])*MHDChildTempDims[field][0]*MHDChildTempDims[field][1]);
	      fieldindex = (k*MagneticDims[field][1] + j)*MagneticDims[field][0] +i;
	      MagneticField[field][fieldindex] = MHDChildTemp[field][tempindex];

	      if(fieldindex >= MagneticSize[field] || tempindex >= MHDChildTempSize[field] ) {
		ENZO_FAIL("InterpolateBoundaryFromParent Bounds Check.\n");
	      }
	      
	    }//right j face (i loop)
	
	// Left k face
	
	for( k=0; k<MHDStartIndex[field][2]; k++)
	  for( j=0; j<MagneticDims[field][1]; j++)
	    for( i=0; i<MagneticDims[field][0]; i++){
	      tempindex = ((i + Offset[0])
			   +(j+ Offset[1])*MHDChildTempDims[field][0]
			   +(k+ Offset[2])*MHDChildTempDims[field][0]*MHDChildTempDims[field][1]);
	      fieldindex = (k*MagneticDims[field][1] + j)*MagneticDims[field][0] +i;
	      MagneticField[field][fieldindex] = MHDChildTemp[field][tempindex];

	      if(fieldindex >= MagneticSize[field] || tempindex >= MHDChildTempSize[field] ) {
		ENZO_FAIL("InterpolateBoundaryFromParent Bounds Check.\n");
	      }
	      
	    }// left k face ( i loop)
	
	// right k face
	
	for( k=MHDEndIndex[field][2]+1; k<MagneticDims[field][2]; k++)
	  for( j=0; j<MagneticDims[field][1]; j++)
	    for( i=0; i<MagneticDims[field][0]; i++){
	      tempindex = ((i + Offset[0])
			   +(j+ Offset[1])*MHDChildTempDims[field][0]
			   +(k+ Offset[2])*MHDChildTempDims[field][0]*MHDChildTempDims[field][1]);
	      fieldindex = (k*MagneticDims[field][1] + j)*MagneticDims[field][0] +i;
	      
	      MagneticField[field][fieldindex] = MHDChildTemp[field][tempindex];
	      if(fieldindex >= MagneticSize[field] || tempindex >= MHDChildTempSize[field] ) {
		ENZO_FAIL("InterpolateBoundaryFromParent Bounds Check.\n");
	      }	      
	    }// right k face ( i loop)
      }//field

      for (field = 0; field < 3; field++){
	delete MHDParentTemp[field];
	MHDParentTemp[field] = NULL;
	delete MHDChildTemp[field];
	
	if(ParentGrid->ProcessorNumber != MyProcessorNumber ){
	  if( ParentGrid->MagneticField[field] != NULL ){
	    delete ParentGrid->MagneticField[field];
	    ParentGrid->MagneticField[field] = NULL;
	  }
	  if( ParentGrid->OldMagneticField[field] != NULL ){
	    delete ParentGrid->OldMagneticField[field];
	    ParentGrid->OldMagneticField[field] = NULL;
	  }
	  if( ParentGrid->ElectricField[field] != NULL ){
	    delete ParentGrid->ElectricField[field];
	    ParentGrid->ElectricField[field] = NULL;
	  }
	}//ParentProc
      }//field
      this->MHD_ProlongFree();
    }//UseMHDCT
  } // end: if (NumberOfBaryonFields > 0)
 
  this->DebugCheck("InterpolateBoundaryFromParent (after)");

  /* Clean up if we have transfered data. */

  if (MyProcessorNumber != ParentGrid->ProcessorNumber)

    ParentGrid->DeleteAllFields();
 
  return SUCCESS;
 
}

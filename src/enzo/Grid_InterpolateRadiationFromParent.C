/***********************************************************************
/
/  GRID CLASS (INTERPOLATE RADIATION FIELD FROM PARENT GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  July, 2009 by John Wise to only consider radiation
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
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
extern "C" void FORTRAN_NAME(interpolate)
                             (int *rank, float *pfield, int pdim[],
			      int pis[], int pie[], int r[],
			      float *field, int dim[], int is[], float *work,
			      interpolation_type *imethod, int *posflag,
			      int *ierror);
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
				     int *sdim1, int *sdim2, int *sdim3,
				     int *ddim1, int *ddim2, int *ddim3,
				     int *sstart1, int *sstart2, int *sstart3,
				     int *dstart1, int *dstart2, int *dstart3);
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
 
int grid::InterpolateRadiationFromParent(grid *ParentGrid, int Field)
{
 
  /* Return if this doesn't involve us. */

  if (RadiativeTransferFLD == FALSE)
    return SUCCESS;
 
  if (ProcessorNumber != MyProcessorNumber &&
      ParentGrid->ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* declarations */

  int StartIndex[MAX_DIMENSION], EndIndex[MAX_DIMENSION];
  int ParentStartIndex[MAX_DIMENSION], ParentTempEndIndex[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION], Offset[MAX_DIMENSION];
  int ZeroVector[MAX_DIMENSION], ParentTempStartIndex[MAX_DIMENSION];
  int ParentTempDim[MAX_DIMENSION], TempDim[MAX_DIMENSION];
  int ParentDim[MAX_DIMENSION];
  int ParentTempSize, WorkSize, TempSize, GridSize, One = 1, Zero = 0;
  int dim, interp_error;
  float *TemporaryField, *Work, *ParentTemp;

  interp_error = FALSE;

  /* Compute refinement factors and set zero. */
 
  ParentGrid->ComputeRefinementFactors(this, Refinement);
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    ZeroVector[dim] = 0;
 
  /* Find field */
 
  int FieldNum;
  if ((FieldNum = FindField(Field, FieldType, NumberOfBaryonFields)) < 0)
    ENZO_FAIL("Cannot find field to interpolate!");
 
  /* Set up array of flags if we are using SecondOrderB interpolation
     method.  These flags indicate if the quantity must always by > 0.
     For Zeus interpolation, they are overloaded as cell vs. face flags. */
 
  int SecondOrderBFlag;
  SecondOrderBFlag = (InterpolationMethod == SecondOrderB) ? 2 : 0;
 
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
 
  ParentTempSize = TempSize = WorkSize = GridSize = 1;
  for (dim = 0; dim < GridRank; dim++) {
 
    StartIndex[dim] = nint((CellLeftEdge[dim][0] -
			    ParentGrid->CellLeftEdge[dim][0])/
			   CellWidth[dim][0]);
    EndIndex[dim] = nint((CellLeftEdge[dim][GridDimension[dim]-1] +
			  CellWidth[dim][GridDimension[dim]-1]
			  - ParentGrid->CellLeftEdge[dim][0])/
			 CellWidth[dim][0]) - 1;
 
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
      ENZO_VFAIL("Parent grid not big enough for interpolation!  ParentStartIndex[%"ISYM"] = %"ISYM"  ParentTempDim = %"ISYM"\n", dim, ParentStartIndex[dim], ParentTempDim[dim])
    }
 
    /* Compute the dimensions of the current grid temporary field. */
 
    TempDim[dim]            = EndIndex[dim] - StartIndex[dim] + 1;
 
    /* Compute size of current grid (in floats) for the temp fields. */
 
    ParentTempSize *= ParentTempDim[dim];
    TempSize       *= TempDim[dim];
    WorkSize       *= (TempDim[dim]/Refinement[dim] + 1);
    GridSize       *= GridDimension[dim];
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
 
  /* Copy data from other processor if needed (modify ParentDim and
     ParentStartIndex to reflect the fact that we are only coping part of
     the grid. */
 
  if (ProcessorNumber != ParentGrid->ProcessorNumber) {
    ParentGrid->CommunicationSendRegion(ParentGrid, ProcessorNumber,
			FieldNum, NEW_ONLY, ParentStartIndex, ParentTempDim);
    for (dim = 0; dim < GridRank; dim++) {
      ParentDim[dim] = ParentTempDim[dim];
      ParentStartIndex[dim] = 0;
    }
  }
 
  /* Return if this is not our concern. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Allocate temporary space. */
 
  TemporaryField = new float[TempSize];
  Work           = new float[WorkSize];
  ParentTemp     = new float[ParentTempSize];
 
  /* Copy just the required section from the parent fields to the
     temp space. */

  FORTRAN_NAME(copy3d)(ParentGrid->BaryonField[FieldNum], ParentTemp,
		       ParentDim, ParentDim+1, ParentDim+2,
		       ParentTempDim, ParentTempDim+1, ParentTempDim+2,
		       &Zero, &Zero, &Zero,
		       ParentStartIndex, ParentStartIndex+1,
		       ParentStartIndex+2);
 
  /* Do the interpolation for the emissivity fields. */
 
  FORTRAN_NAME(interpolate)(&GridRank,
			    ParentTemp, ParentTempDim,
			    ParentTempStartIndex, ParentTempEndIndex,
			    Refinement,
			    TemporaryField, TempDim, ZeroVector, Work,
			    &InterpolationMethod,
			    &SecondOrderBFlag, &interp_error);
 
    if (interp_error) {
      ENZO_VFAIL("P%d: Error interpolating emissivity.\n"
		 "ParentGrid ID = %d\n"
		 "\t LeftEdge  = %"PSYM" %"PSYM" %"PSYM"\n"
		 "\t RightEdge = %"PSYM" %"PSYM" %"PSYM"\n"
		 "ThisGrid ID = %d\n"
		 "\t LeftEdge  = %"PSYM" %"PSYM" %"PSYM"\n"
		 "\t RightEdge = %"PSYM" %"PSYM" %"PSYM"\n",
		 MyProcessorNumber, ParentGrid->ID, 
		 ParentGrid->GridLeftEdge[0], ParentGrid->GridLeftEdge[1], 
		 ParentGrid->GridLeftEdge[2], ParentGrid->GridRightEdge[0], 
		 ParentGrid->GridRightEdge[1], ParentGrid->GridRightEdge[2],
		 this->ID, 
		 this->GridLeftEdge[0], this->GridLeftEdge[1], 
		 this->GridLeftEdge[2], this->GridRightEdge[0], 
		 this->GridRightEdge[1], this->GridRightEdge[2])
    }

 
  /* Copy needed portion of temp field to current grid. */
 
  if (BaryonField[Field] == NULL)
    BaryonField[Field] = new float[GridSize];
  if (BaryonField[Field] == NULL)
    ENZO_FAIL("malloc error (out of memory?)");

  FORTRAN_NAME(copy3d)(TemporaryField, BaryonField[FieldNum],
		       TempDim, TempDim+1, TempDim+2,
		       GridDimension, GridDimension+1, GridDimension+2,
		       &Zero, &Zero, &Zero,
		       Offset, Offset+1, Offset+2);
 
  delete [] Work;
  delete [] TemporaryField;
  delete [] ParentTemp;
 
  /* Clean up if we have transfered data. */
 
  if (MyProcessorNumber != ParentGrid->ProcessorNumber) {

    delete [] BaryonField[FieldNum];
    BaryonField[FieldNum] = NULL;
  }
 
  return SUCCESS;
}

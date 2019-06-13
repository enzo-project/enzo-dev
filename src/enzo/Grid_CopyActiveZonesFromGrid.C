/***********************************************************************
/
/  GRID CLASS (COPY OVERLAPPING ACTIVE ZONES FROM GRID IN ARGUMENT TO THIS GRID)
/
/  written by: Nathan Goldbaum
/  date:       July 2012 (adapted from Grid_CopyZonesFromGrid.C)
/  modified1:  Stephen Skory, Sept 2012
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
// This routine copies active zones (no ghost zones) which overlap
//   from the grid in the argument to the current grid.  We use only
//   the active region of the OtherGrid, but copy into the entire
//   region (including boundaries) of this grid.
//
// The input argument EdgeOffset is the amount the corner of this grid is
//   considered to have moved for grid comparison and copying purposes.
//   See Grid_CheckForOverlappingZones for more details.

#ifdef USE_MPI
#endif /* USE_MPI */
 
#include "preincludes.h"
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

int check_overlap(FLOAT a[MAX_DIMENSION], FLOAT b[MAX_DIMENSION],
                  FLOAT s[MAX_DIMENSION], FLOAT t[MAX_DIMENSION],
                  int dims, FLOAT period[MAX_DIMENSION],
                  int *shift, int skipto);

int grid::CopyActiveZonesFromGrid(grid *OtherGrid, FLOAT EdgeOffset[MAX_DIMENSION],
    int SendField)
{
 
  /* Return if this doesn't involve us. */
  if (ProcessorNumber != MyProcessorNumber &&
      OtherGrid->ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
    }
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  if (this->GetCellWidth(0,0) != OtherGrid->GetCellWidth(0,0))
    return SUCCESS;

  /* Compute the left and right edges of this grid (including ghost zones). */
 
  FLOAT GridLeft[MAX_DIMENSION]; FLOAT GridRight[MAX_DIMENSION];
  FLOAT ActiveLeft[MAX_DIMENSION]; FLOAT ActiveRight[MAX_DIMENSION];
  FLOAT period[MAX_DIMENSION];
  int *shift;

  int adjust;

  shift = new int[MAX_DIMENSION];

  int dim;

  for (dim = 0; dim < GridRank; dim++) {
    period[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    shift[dim] = -1;
  }

  /* Check to see if there's overlap. Grids can overlap with more than one
     shift vector, so we need to keep checking until we've got them all.
     The loop breaks when check_overlap returns -1.
   */
  int overlap = 0;
  while (true) 
    {
   
      for (dim = 0; dim < GridRank; dim++) {
        ActiveLeft[dim]  = GridLeftEdge[dim]  + EdgeOffset[dim];
        ActiveRight[dim] = GridRightEdge[dim] + EdgeOffset[dim];
        GridLeft[dim]  = CellLeftEdge[dim][0] + EdgeOffset[dim];
        GridRight[dim] = CellLeftEdge[dim][GridDimension[dim]-1] +
          CellWidth[dim][GridDimension[dim]-1]    +
          EdgeOffset[dim];
      }

      overlap = check_overlap(ActiveLeft, ActiveRight,
        OtherGrid->GridLeftEdge, OtherGrid->GridRightEdge,
        GridRank, period, shift, overlap);
      if (overlap == -1) {
        break;
      }
      /* If we're here, then there is overlap, and we should fix our edges
         so the math works out. */
      FLOAT shift_temp;
      for (dim = 0; dim < GridRank; dim++) {
        shift_temp = period[dim] * FLOAT(shift[dim]);
        ActiveLeft[dim]  = GridLeftEdge[dim]  + shift_temp;
        ActiveRight[dim] = GridRightEdge[dim] + shift_temp;
        GridLeft[dim]  = CellLeftEdge[dim][0] + shift_temp;
        GridRight[dim] = CellLeftEdge[dim][GridDimension[dim]-1] +
          CellWidth[dim][GridDimension[dim]-1] + shift_temp;
      }

      /* There is some overlap, so copy overlapping region */
     
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
     
      for (dim = 0; dim < GridRank; dim++){
        if (GridDimension[dim] > 1) {
        shift_temp = period[dim] * FLOAT(shift[dim]);
          /* Compute left and right positions in problem space.
         note: include buffer zones of this grid but not the other grid. */
     
          Left[dim]  = max(ActiveLeft[dim], OtherGrid->GridLeftEdge[dim]);
          Right[dim] = min(ActiveRight[dim], OtherGrid->GridRightEdge[dim]);
     
          /* Convert this to index positions in this grid */
     
          Start[dim] = nint((Left[dim]  - GridLeft[dim]) / CellWidth[dim][0]);
          End[dim]   = nint((Right[dim] - GridLeft[dim]) / CellWidth[dim][0]) - 1;
    
          if (SendField == GRAVITATING_MASS_FIELD) {
            adjust = (GravitatingMassFieldDimension[dim] - GridDimension[dim]) / 2;
            Start[dim] += adjust;
            End[dim] += adjust;
          }
    
          if (End[dim] - Start[dim] < 0) {
            delete [] shift;
            return SUCCESS;
          }

          Dim[dim] = End[dim] - Start[dim] + 1;  
    
          /* Compute index positions in the other grid */
          StartOther[dim] = nint((Left[dim] - OtherGrid->CellLeftEdge[dim][0])/
                   CellWidth[dim][0]);
          
          if (SendField == GRAVITATING_MASS_FIELD) {
            adjust = (OtherGrid->GravitatingMassFieldDimension[dim] - 
              OtherGrid->GridDimension[dim]) / 2;
            StartOther[dim] += adjust;
          }
     
          /* Copy dimensions into temporary space */
          
          if (SendField == ALL_FIELDS) {
            OtherDim[dim] = OtherGrid->GridDimension[dim];
          } else if (SendField == GRAVITATING_MASS_FIELD) {
            OtherDim[dim] = OtherGrid->GravitatingMassFieldDimension[dim];
          }
        } // if GD > 1
      } // dim
      /* Calculate dimensions */
      
      /* If posting a receive, then record details of call. */
    
    #ifdef USE_MPI
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE &&
          MyProcessorNumber == ProcessorNumber) {
        CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
        CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = OtherGrid;
        CommunicationReceiveCallType[CommunicationReceiveIndex] = 21;
        for (dim = 0; dim < GridRank; dim++)
          CommunicationReceiveArgument[dim][CommunicationReceiveIndex] = 
            EdgeOffset[dim];
        }
        CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] =
            SendField;
    #endif /* USE_MPI */
      /* Copy data from other processor if needed (modify OtherDim and
         StartOther to reflect the fact that we are only copying part of
         the grid. */
     
      if (traceMPI) 
        fprintf(tracePtr, "CopyZones SendRegion from %"ISYM" to %"ISYM"\n", 
            ProcessorNumber, OtherGrid->ProcessorNumber);
      
      if (ProcessorNumber != OtherGrid->ProcessorNumber) {
        OtherGrid->CommunicationSendRegion(OtherGrid, ProcessorNumber,
                           SendField, NEW_ONLY, StartOther, Dim);
        
        if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	    CommunicationDirection == COMMUNICATION_SEND) {
          delete [] shift;
	  return SUCCESS;
	}

        for (dim = 0; dim < GridRank; dim++) {
          OtherDim[dim]=Dim[dim];
          StartOther[dim] = 0;
        }
      }
    
      /* Return if this is not our concern. */
     
      if (ProcessorNumber != MyProcessorNumber) {
        delete [] shift;
	return SUCCESS;
        }
    
      if (SendField == ALL_FIELDS) {
          for (int field = 0; field < NumberOfBaryonFields; field++) {
            FORTRAN_NAME(copy3drel)(OtherGrid->BaryonField[field], BaryonField[field],
                        Dim, Dim+1, Dim+2,
                        OtherDim, OtherDim+1, OtherDim+2,
                        GridDimension, GridDimension+1, GridDimension+2,
                        StartOther, StartOther+1, StartOther+2,
                        Start, Start+1, Start+2);
           }
      } else if (SendField == GRAVITATING_MASS_FIELD) {
        FORTRAN_NAME(copy3drel)(OtherGrid->GravitatingMassField,
                    GravitatingMassField,
                    Dim, Dim+1, Dim+2,
                    OtherDim, OtherDim+1, OtherDim+2,
                    GravitatingMassFieldDimension,
                    GravitatingMassFieldDimension+1,
                    GravitatingMassFieldDimension+2,
                    StartOther, StartOther+1, StartOther+2,
                    Start, Start+1, Start+2);
      }
      
      /* Clean up if we have transfered data. */

  } // end while(true)
  if (MyProcessorNumber != OtherGrid->ProcessorNumber) {
    OtherGrid->DeleteAllFields();
    }

  delete [] shift;

  this->DebugCheck("CopyActiveZonesFromGrid (after)");

  return SUCCESS;
}
  

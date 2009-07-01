/***********************************************************************
/
/  INTERPOLATE PARTICLE QUANTITIES TO GRID
/
/  written by: John Wise
/  date:       June, 2009
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <string.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "communication.h"

#include "FOF_allvars.h"
#include "FOF_nrutil.h"
#include "FOF_proto.h"

int FindField(int field, int farray[], int numfields);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */

int grid::InterpolateParticlesToGrid(FOFData *D)
{

  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
      CommunicationDirection == COMMUNICATION_RECEIVE)
    if (NumberOfProcessors == 1 || MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

  if (CommunicationDirection == COMMUNICATION_SEND && D == NULL)
    ENZO_FAIL("D (FOFData) cannot be null in send-mode!");

  int ActiveDim[MAX_DIMENSION];
  int min_slab, max_slab, NumberOfFields = 0;
  int i, j, k, n, field, dim, proc, size;
  int DensNum, VrmsNum, Vel1Num, Vel2Num, Vel3Num;

#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
  MPI_Arg Source;
  MPI_Arg Tag;
  float *buffer;
#endif /* USE_MPI */

  /* Calculate the overlapping slabs */

  if (NumberOfProcessors == 1) {
    min_slab = 0;
    max_slab = 0;
  } else {
    min_slab = (int) (NumberOfProcessors * (GridLeftEdge[0] - DomainLeftEdge[0]) /
		      (DomainRightEdge[0] - DomainLeftEdge[0]));
    max_slab = (int) (NumberOfProcessors * (GridRightEdge[0] - DomainLeftEdge[0]) /
		      (DomainRightEdge[0] - DomainLeftEdge[0]));
    min_slab = max(min_slab, 0);
    max_slab = min(max_slab, NumberOfProcessors-1);
  }

  /* Post receive mode :: allocate memory and post receives */

  switch (OutputSmoothedDarkMatter) {
  case 1: NumberOfFields = 1; break;  // density
  case 2: NumberOfFields = 5; break;  // + rms velocity + 3-velocity
  default: 
    fprintf(stdout, "Unrecognized value for OutputSmoothedDarkMatter = %"ISYM"\n",
	    OutputSmoothedDarkMatter);
    fprintf(stdout, "Setting to 1.  Outputting smoothed density only.\n");
    OutputSmoothedDarkMatter = 1;
    NumberOfFields = 1;
    break;
  } // ENDSWITCH

  /* Assign field numbers */

  DensNum = 0;
  VrmsNum = 1;
  Vel1Num = 2;
  Vel2Num = 3;
  Vel3Num = 4;

  size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    size *= ActiveDim[dim];
  }

  /************************* POST-RECEIVE MODE *************************/

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {

    Count = size;
    for (proc = min_slab; proc <= max_slab; proc++)
      if (proc != MyProcessorNumber) {
	
	Source = proc;
	for (field = 0; field < NumberOfFields; field++) {

	  buffer = new float[size];
	  Tag = MPI_SENDPARTFIELD_TAG+field;
	  MPI_Irecv(buffer, Count, DataType, Source, Tag, MPI_COMM_WORLD,
		    CommunicationReceiveMPI_Request+CommunicationReceiveIndex);

	  CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	  CommunicationReceiveGridTwo[CommunicationReceiveIndex] = NULL;
	  CommunicationReceiveCallType[CommunicationReceiveIndex] = 17;
	  CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer;
	  CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = field;
	  CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
	    CommunicationReceiveCurrentDependsOn;
	  CommunicationReceiveIndex++;

	} // ENDFOR field

      } // ENDIF different processor
  } // ENDIF post_receive
#endif /* USE_MPI */    

  /************************* SEND MODE *************************/
  // Interpolate from particles to grid here.  Then send if needed.

  if (CommunicationDirection == COMMUNICATION_SEND) {

    float *r2list = NULL;
    int *ngblist = NULL;


    int slab, ind, ik, SlabStartIndex, SlabEndIndex, index;
    FLOAT SlabLeftEdge, SlabRightEdge;
    double CellPos[MAX_DIMENSION];
    double r, h, h2, hinv, hinv3, u, delv, weight;
    double *wk = NULL;

    FLOAT a, dadt, CurrentRedshift = 0.0;
    float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
      DensityUnits, LengthConversion;
    float UnitConversion[MAX_NUMBER_OF_BARYON_FIELDS];


    // Exit if not overlapping.  But we still need to allocate memory if
    // this is the host processor.

    if (MyProcessorNumber < min_slab || MyProcessorNumber > max_slab) {
      if (MyProcessorNumber == ProcessorNumber)
	for (field = 0; field < NumberOfFields; field++) {
	  InterpolatedField[field] = new float[size];
	  for (i = 0; i < size; i++)
	    InterpolatedField[field][i] = 0.0;
	}
      return SUCCESS;
    }


    // Get units
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	     &TimeUnits, &VelocityUnits, Time);
    if (ComovingCoordinates) {
      CosmologyComputeExpansionFactor(Time, &a, &dadt);
      CurrentRedshift = (1 + InitialRedshift)/a - 1;
    }
    
    // 1e10 Msun / (comoving kpc)^3 -> DensityUnits
    // 6.77e-22 = 1e10 Msun / kpc^3 in amu
    UnitConversion[DensNum] = 6.76779e-22 * pow(1+CurrentRedshift, 3) /
      DensityUnits;

    // proper km/s -> VelocityUnits
    for (field = VrmsNum; field <= Vel3Num; field++)
      UnitConversion[field] = 1e5 / VelocityUnits;

    // Determine where the slab starts / ends in the grid.
    slab = MyProcessorNumber;
    SlabLeftEdge = slab * (DomainRightEdge[0] - DomainLeftEdge[0]) /
      NumberOfProcessors;
    SlabRightEdge = (slab+1) * (DomainRightEdge[0] - DomainLeftEdge[0]) /
      NumberOfProcessors;
    
    SlabStartIndex = (int) ((SlabLeftEdge - GridLeftEdge[0]) / CellWidth[0][0]);
    SlabStartIndex = max(SlabStartIndex, 0);

    SlabEndIndex = (int) ((SlabRightEdge - GridLeftEdge[0]) / CellWidth[0][0]) - 1;
    SlabEndIndex = min(SlabEndIndex, ActiveDim[0]-1);

    // Allocate and zero memory
    for (field = 0; field < NumberOfFields; field++) {
      InterpolatedField[field] = new float[size];
      for (i = 0; i < size; i++)
	InterpolatedField[field][i] = 0.0;
    }

    wk = new double[D->DesDensityNgb];

    // Tree is built in comoving kpc, so CellPos is.
    LengthConversion = (LengthUnits / 3.086e21) * (1.0 + CurrentRedshift);

    for (k = 0; k < ActiveDim[2]; k++) {
      CellPos[2] = LengthConversion * 
	(GridLeftEdge[2] + (k+0.5) * CellWidth[2][0]);
      for (j = 0; j < ActiveDim[1]; j++) {
	CellPos[1] = LengthConversion * 
	  (GridLeftEdge[1] + (j+0.5) * CellWidth[1][0]);
	index = SlabStartIndex + ActiveDim[0] * (j + ActiveDim[1]*k);
	for (i = SlabStartIndex; i <= SlabEndIndex; i++, index++) {

	  CellPos[0] = LengthConversion * 
	    (GridLeftEdge[0] + (i+0.5) * CellWidth[0][0]);
	  h2 = ngb_treefind(D->P, CellPos, D->DesDensityNgb, 0, &ngblist, &r2list);
	  h = sqrt(h2);
	  hinv = 1.0 / h;
	  hinv3 = hinv*hinv*hinv;
	  weight = 0;

	  for (n = 0; n < D->DesDensityNgb; n++) {
	    r = sqrt(r2list[n]);
	    ind = ngblist[n];
	    if (r < h) {
	      u = r*hinv;
	      ik = (int) (u * KERNEL_TABLE);
	      wk[n] = D->P[ind].Mass * 
		( D->Kernel[ik] + (D->Kernel[ik+1] - D->Kernel[ik]) *
		  (u - D->KernelRad[ik]) * KERNEL_TABLE);
	      weight += D->P[ind].Mass;

	      // Density
	      InterpolatedField[DensNum][index] += wk[n];

	      // Particle velocities and velocity dispersions (later)
	      if (OutputSmoothedDarkMatter > 1)
		for (dim = 0; dim < 3; dim++)
		  InterpolatedField[Vel1Num+dim][index] += 
		    wk[n] * D->P[ind].Vel[dim];
	    } // ENDIF r < h
	    else
	      wk[n] = 0.0;

	  } // ENDFOR n

	  // After we've calculated the bulk velocity, we can compute
	  // the rms velocity.
	  if (OutputSmoothedDarkMatter > 1)
	    for (n = 0; n < D->DesDensityNgb; n++) {
	      ind = ngblist[n];
	      for (dim = 0; dim < 3; dim++) {
		delv = InterpolatedField[Vel1Num+dim][index] - D->P[ind].Vel[dim];
		InterpolatedField[VrmsNum][index] += wk[n] * delv * delv;
	      } // ENDFOR dim
	    } // ENDFOR n

	  // Normalize and convert to enzo units
	  InterpolatedField[DensNum][index] *= hinv3;

	  for (field = VrmsNum; field <= Vel3Num && field < NumberOfFields; field++)
	    InterpolatedField[field][index] /= weight;

	  if (OutputSmoothedDarkMatter > 1)
	    InterpolatedField[VrmsNum][index] =
	      sqrt(InterpolatedField[VrmsNum][index]);

	  for (field = 0; field < NumberOfFields; field++)
	    InterpolatedField[field][index] *= UnitConversion[field];

	} // ENDFOR i
      } // ENDFOR j
    } // ENDFOR k

    delete [] wk;

#ifdef USE_MPI
    if (MyProcessorNumber != ProcessorNumber) {

      for (field = 0; field < NumberOfFields; field++) {
	Tag = MPI_SENDPARTFIELD_TAG + field;
	CommunicationBufferedSend(InterpolatedField[field], size, DataType,
				  ProcessorNumber, Tag, MPI_COMM_WORLD,
				  size * sizeof(float));
	delete [] InterpolatedField[field];
	InterpolatedField[field] = NULL;
      } // ENDFOR field
      
    } // ENDIF different processors
#endif /* USE_MPI */

  } // ENDIF send

  /************************* RECEIVE MODE *************************/  
  // Sum the data that's received.

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_RECEIVE) {

    buffer = CommunicationReceiveBuffer[CommunicationReceiveIndex];
    field = CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex];
    for (i = 0; i < size; i++)
      InterpolatedField[field][i] += buffer[i];
    delete [] buffer;

  } // ENDIF receive
#endif /* USE_MPI */

  return SUCCESS;

}

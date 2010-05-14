/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

/* function prototypes */

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

#define PHOTON_WEIGHT 5.1

int CommunicationLoadBalancePhotonGrids(HierarchyEntry *GridHierarchyPointer[],
					int NumberOfGrids)
{

  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;

#ifdef TRANSFER // nee TRANSFER

  if (RadiativeTransfer == FALSE)
    return SUCCESS;

  /* Initialize */

  int i, dim, GridMemory, NumberOfCells, CellsTotal, Particles;
  float AxialRatio, GridVolume;
  float rad_factor;
  float *ComputeTime = new float[NumberOfGrids];
  float *ProcessorComputeTime = new float[NumberOfProcessors];

  int *GridRadiation = new int[NumberOfGrids];
  FLOAT r_photon[MAX_DIMENSION];
//  for (dim = 0; dim < MAX_DIMENSION; dim++)
//    r_photon[dim] = new FLOAT[NumberOfGrids];
  long Nside = (long) ceil(sqrt(NumberOfProcessors / 12.0));
  long Npix  = 12 * Nside * Nside;

  for (i = 0; i < NumberOfProcessors; i++)
    ProcessorComputeTime[i] = 0;
  
  /* Compute work for each grid. */

  for (i = 0; i < NumberOfGrids; i++) {

    if (GridHierarchyPointer[i]->GridData->RadiationPresent() == TRUE)
      continue;

    GridHierarchyPointer[i]->GridData->
      CollectGridInformation(GridMemory, GridVolume, NumberOfCells, 
			     AxialRatio, CellsTotal, Particles);
    //    ComputeTime[i] = GridMemory; // roughly speaking
    ComputeTime[i] = float(NumberOfCells);

  /* For grids with radiation and if there exists some radiation
     source(s), add an additional ComputeTime, find the closest
     radiation source, and determine the processor by grouping the
     grids by HEALPIX to reduce the communication of photons. */

    if (GlobalRadiationSources->NextSource != NULL) {

      FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION], Center[MAX_DIMENSION];
      int Rank, Dims[MAX_DIMENSION];
      
      GridHierarchyPointer[i]->GridData->ReturnGridInfo(&Rank,Dims,Left,Right);
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	Center[dim] = 0.5 * (Left[dim] + Right[dim]);

      /* Find the closest radiation source and record the r-vector
	 from the source to the grid center */

      RadiationSourceEntry *RS;
      RS = GlobalRadiationSources->NextSource;
      float radius2, delx, dely, delz, min_radius2 = huge_number;
      int ToProcessor, FromProcessor;
      long ipix;

      while (RS != NULL) {
	radius2 = 0;
	delx = RS->Position[0] - Center[0];
	dely = RS->Position[1] - Center[1];
	delz = RS->Position[2] - Center[2];
	radius2 = delx*delx + dely*dely + delz*delz;

	if (radius2 < min_radius2) {
	  min_radius2 = radius2;
	  r_photon[0] = (FLOAT) delx;
	  r_photon[1] = (FLOAT) dely;
	  r_photon[2] = (FLOAT) delz;
	}

	RS = RS->NextSource;

      } // ENDWHILE RADIATION SOURCE

      /* Compute the HEALPIX pixel where this grid exists, and move it
	 to the corresponding processor if not there */

      if (vec2pix_nest(Nside, r_photon, &ipix) == FAIL) {
	fprintf(stderr, "Error in vec2pix_nest.\n");
	ENZO_FAIL("");
      }

      // Group pixels together according to processor number
      // i.e. for 12 pixels and np=4 :: 0 0 0 1 1 1 2 2 2 3 3 3
      ToProcessor = (int)
	floor(float(ipix) / ((12.0*Nside*Nside) / float(NumberOfProcessors)));
      FromProcessor = GridHierarchyPointer[i]->GridData->ReturnProcessorNumber();

      GridHierarchyPointer[i]->GridData->CommunicationMoveGrid(ToProcessor);

      ProcessorComputeTime[ToProcessor] += ComputeTime[i];

//      if (debug) {
//	fprintf(stdout, "LB[G%"ISYM"]: nside = %"ISYM", ncells = %"ISYM", Center = %"GSYM" %"GSYM" %"GSYM"\n",
//		i, (int)Nside, NumberOfCells, Center[0], Center[1], Center[2]);
//	fprintf(stdout, "LB[G%"ISYM"]: r_vec = %"GSYM" %"GSYM" %"GSYM", ipix = %"ISYM", ToProc = %"ISYM"\n",
//		i, delx, dely, delz, (int)ipix, ToProcessor);
//      }
//      if (FromProcessor != ToProcessor) {
//	if (debug)
//	  printf("LB[G%"ISYM"]: Moving grid P%"ISYM"->P%"ISYM"\n", i,FromProcessor,ToProcessor);
//	GridHierarchyPointer[i]->GridData->CommunicationMoveGrid(ToProcessor);
//      }

    }

  } // ENDFOR grids

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("LoadBalancePhotons (grids=%"ISYM"): ", NumberOfGrids);
    WriteListOfFloats(stdout, NumberOfProcessors, ProcessorComputeTime);
  }

  delete [] ComputeTime;
  delete [] ProcessorComputeTime;

#endif /* TRANSFER */

  return SUCCESS;
}

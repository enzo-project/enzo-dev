/***********************************************************************
/
/  INLINE HALO FINDER HELPER FUNCTION :: INITIALIZATION
/
/  written by: John Wise
/  date:       June, 2009
/
/  PURPOSE:    Moves particle data from enzo's memory to the halo 
/              finder's memory.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
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
#include "CosmologyParameters.h"

#include "FOF_allvars.h"
#include "FOF_nrutil.h"
#include "FOF_proto.h"

/************************************************************************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/************************************************************************/

void FOF_Initialize(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], 
		    FOFData &D)
{

  /* Get enzo units */

  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits, MassUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, &MassUnits, MetaData->Time);
  
  // Mpc/h -> kpc
  D.BoxSize = 1e3 * ComovingBoxSize / HubbleConstantNow;
  
  // Time = a = 1/(1+z)
  FLOAT CurrentRedshift, a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(MetaData->Time, &a, &dadt);
  D.Time = a;

  // Sometimes MassUnits is infinite (in cgs) when using single
  // precision, so calculate it in double precision.

  double EnzoMassUnits = (double) DensityUnits * pow(LengthUnits, 3.0);

  // Copy other cosmology parameters
  D.Omega = OmegaMatterNow;
  D.OmegaLambda = OmegaLambdaNow;

  /****************** MOVE PARTICLES TO P-GROUPFINDER ******************/

  int i, j, proc, level, slab;
  int GridNum, Index, NumberOfLocalParticles, ptype_size;
  LevelHierarchyEntry *Temp;
  FOF_particle_data *Plocal;
  MPI_Arg *Nslab_local, *NtoLeft_local, *NtoRight_local;

#ifdef USE_MPI
  MPI_Datatype IntType = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
#endif
  ptype_size = sizeof(FOF_particle_data);

  // Count local particles
  NumberOfLocalParticles = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber())
	NumberOfLocalParticles += Temp->GridData->ReturnNumberOfParticles();

  // Allocate memory
  Plocal = new FOF_particle_data[NumberOfLocalParticles];

  // Move particles
  Index = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level], GridNum = 0; Temp; 
	 Temp = Temp->NextGridThisLevel, GridNum++)
      Temp->GridData->MoveParticlesFOF(level, GridNum, Plocal, Index, D,
				       VelocityUnits, EnzoMassUnits, COPY_OUT);

  /************* Count particles in each slab and shadow ****************/

  Nslab_local = new MPI_Arg[NumberOfProcessors];
  NtoLeft_local = new MPI_Arg[NumberOfProcessors];
  NtoRight_local = new MPI_Arg[NumberOfProcessors];
  D.Nslab = new int[NumberOfProcessors];
  D.NtoLeft = new int[NumberOfProcessors];
  D.NtoRight = new int[NumberOfProcessors];
  D.Nshadow = new int[NumberOfProcessors];
  D.Noffset = new int[NumberOfProcessors];

  for (proc = 0; proc < NumberOfProcessors; proc++) {
    Nslab_local[proc] = 0;
    NtoLeft_local[proc] = 0;
    NtoRight_local[proc] = 0;
  } 

  for (i = 0; i < NumberOfLocalParticles; i++) {

    slab = Plocal[i].slab;
    Nslab_local[slab]++;

    // Left and right "shadows"
    if (Plocal[i].Pos[0] < (slab * D.BoxSize / NumberOfProcessors + 
			    D.SearchRadius))
      NtoLeft_local[slab]++;

    if (Plocal[i].Pos[0] > ((slab+1) * D.BoxSize / NumberOfProcessors -
			    D.SearchRadius))
      NtoRight_local[slab]++;

  } // ENDFOR particles

#ifdef USE_MPI
  // Get counts over all processors
  MPI_Allreduce(Nslab_local, D.Nslab, NumberOfProcessors, IntType,
		MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(NtoLeft_local, D.NtoLeft, NumberOfProcessors, IntType,
		MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(NtoRight_local, D.NtoRight, NumberOfProcessors, IntType,
		MPI_SUM, MPI_COMM_WORLD);
#endif

  D.NumPart = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++) {
    D.NumPart += D.Nslab[proc];

    if (proc < NumberOfProcessors-1)
      D.Nshadow[i] += D.NtoLeft[i+1];
    else
      D.Nshadow[i] += D.NtoLeft[0];

    if (proc > 0)
      D.Nshadow[i] += D.NtoRight[i-1];
    else
      D.Nshadow[i] += D.NtoRight[NumberOfProcessors-1];

  } // ENDFOR processors

  for (proc = 0; proc < NumberOfProcessors; proc++)
    for (i = 0, D.Noffset[i] = 0; j < i; j++)
      D.Noffset[i] += D.Nslab[j];

  /* Now we need to sort the local particles by slab, so we can send
     them to their correct processor (==slab) */

  qsort(Plocal, NumberOfLocalParticles, ptype_size, compare_slab);

  /************************ COMMUNICATION ***************************
    We now send particles to the processor which hosts their slab.
    Only after each processor has their slab, we can construct the 
    tree and find the halos
  *******************************************************************/

  int TotalLocal, TotalRecv;
  MPI_Arg *disp_local, *disp_recv, *Nslab_recv;

  disp_local = new MPI_Arg[NumberOfProcessors];
  disp_recv  = new MPI_Arg[NumberOfProcessors];
  Nslab_recv = new MPI_Arg[NumberOfProcessors];

  allocate_memory(D);
  TotalLocal = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++) {
    disp_local[proc] = TotalLocal;
    TotalLocal += Nslab_local[proc];
  }

  // First gather the number of local particles on each processor
#ifdef USE_MPI
  MPI_Alltoall(Nslab_local, 1, IntType, Nslab_recv, 1, IntType, MPI_COMM_WORLD);
#endif

  TotalRecv = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++) {
    disp_recv[proc] = TotalRecv;
    TotalRecv += Nslab_recv[proc];
  }
  D.Nlocal = TotalRecv;

  // Because we're sending the structures as bytes, multiply
  // everything by the size of the particle structure
  for (proc = 0; proc < NumberOfProcessors; proc++) {
    Nslab_local[proc] *= ptype_size;
    Nslab_recv[proc]  *= ptype_size;
    disp_local[proc]  *= ptype_size;
    disp_recv[proc]   *= ptype_size;
  }

  // Now we can do the big collective call
#ifdef USE_MPI
  MPI_Alltoallv(Plocal, Nslab_local, disp_local, MPI_BYTE,
		D.P+1, Nslab_recv, disp_recv, MPI_BYTE,
		MPI_COMM_WORLD);
#endif

  delete [] Plocal;
  delete [] Nslab_local;
  delete [] NtoLeft_local;
  delete [] NtoRight_local;
  delete [] Nslab_recv;
  delete [] disp_local;
  delete [] disp_recv;

  return;

}

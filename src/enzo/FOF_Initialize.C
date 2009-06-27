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
#include <unistd.h>
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
int SysMkdir (char *startDir, char *directory);

/************************************************************************/

void FOF_Initialize(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], 
		    FOFData &D)
{

  /* Check if the FOF directory exists */

  int unixresult;
  char FOF_dir[MAX_LINE_LENGTH];

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    strcpy(FOF_dir, MetaData->GlobalDir);
    strcat(FOF_dir, "/FOF");
    if (access(FOF_dir, F_OK) == -1)
      SysMkdir("", FOF_dir);
  }

  /* Get enzo units */

  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits, MassUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, &MassUnits, MetaData->Time);
  
  // Mpc/h -> kpc
  D.BoxSize = 1e3 * ComovingBoxSize / HubbleConstantNow;
  
  // Time = a = 1/(1+z).  In enzo, the scale factor is in units of (1+z0).
  FLOAT CurrentRedshift, a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(MetaData->Time, &a, &dadt);
  D.Time = a / (1 + InitialRedshift);

  // Critical density in units of Msun / kpc^3
  D.RhoCritical0 = 1.4775867e31 * 
    ((3 * pow(100 * HubbleConstantNow / 3.086e19, 2)) / (8 * M_PI * GRAVITY));
  //D.RhoCritical /= pow(D.Time, 3);

  // Sometimes MassUnits is infinite (in cgs) when using single
  // precision, so calculate it in double precision.

  int TopGridDims3 = MetaData->TopGridDims[0] * MetaData->TopGridDims[1] * 
    MetaData->TopGridDims[2];

  double EnzoMassUnits = (double) DensityUnits * pow(LengthUnits, 3.0) /
    TopGridDims3;

  // Copy other cosmology parameters
  D.Omega = OmegaMatterNow;
  D.OmegaLambda = OmegaLambdaNow;

  // Get total number of particles to set search radius and softening
  // length for potential computation.  Be careful about nested grid
  // simulations.

  int i, level, FinestStaticLevel = 0;
  LevelHierarchyEntry *Temp;
  
  D.NumPart = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      D.NumPart += Temp->GridData->ReturnNumberOfParticles();

  // Search for the finest nested static grid
  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    FinestStaticLevel = max(StaticRefineRegionLevel[i]+1, FinestStaticLevel);

  // TODO: corrections to SearchRadius and Epsilon for nested grid sims.

  D.SearchRadius = D.LinkLength * D.BoxSize / pow(TopGridDims3, 1.0/3) / 
    pow(RefineBy, FinestStaticLevel);

  // softening length for potential computation
  D.Epsilon = 0.05 * D.BoxSize / pow(TopGridDims3, 1.0/3) / 
    pow(RefineBy, FinestStaticLevel);
  
  if (debug) {
    fprintf(stdout, "Inline halo finder starting...\n");
    fprintf(stdout, "FOF: Comoving linking length: %g kpc\n", D.SearchRadius);
  }


  /****************** MOVE PARTICLES TO P-GROUPFINDER ******************/

  int j, proc, slab;
  int GridNum, Index, NumberOfLocalParticles, ptype_size;
  FOF_particle_data *Plocal;
  MPI_Arg *Nslab_local, *NtoLeft_local, *NtoRight_local;

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
	 Temp = Temp->NextGridThisLevel, GridNum++) {
      Temp->GridData->MoveParticlesFOF(level, GridNum, Plocal, Index, D,
				       VelocityUnits, EnzoMassUnits, COPY_OUT);
      //Temp->GridData->SetNumberOfParticles(0);
    }

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
    D.Nslab[proc] = 0;
    D.NtoLeft[proc] = 0;
    D.NtoRight[proc] = 0;
    D.Nshadow[proc] = 0;
    D.Noffset[proc] = 0;
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
  MPI_Allreduce(Nslab_local, D.Nslab, NumberOfProcessors, IntDataType,
		MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(NtoLeft_local, D.NtoLeft, NumberOfProcessors, IntDataType,
		MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(NtoRight_local, D.NtoRight, NumberOfProcessors, IntDataType,
		MPI_SUM, MPI_COMM_WORLD);
#endif

  for (proc = 0; proc < NumberOfProcessors; proc++) {

    if (proc < NumberOfProcessors-1)
      D.Nshadow[proc] += D.NtoLeft[proc+1];
    else
      D.Nshadow[proc] += D.NtoLeft[0];

    if (proc > 0)
      D.Nshadow[proc] += D.NtoRight[proc-1];
    else
      D.Nshadow[proc] += D.NtoRight[NumberOfProcessors-1];

  } // ENDFOR processors

  for (proc = 0; proc < NumberOfProcessors; proc++)
    for (j = 0, D.Noffset[proc] = 0; j < proc; j++)
      D.Noffset[proc] += D.Nslab[j];

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
  MPI_Alltoall(Nslab_local, 1, IntDataType, Nslab_recv, 1, IntDataType, MPI_COMM_WORLD);
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

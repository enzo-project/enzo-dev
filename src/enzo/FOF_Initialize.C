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
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int SysMkdir (char *startDir, char *directory);

/************************************************************************/

void FOF_Initialize(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], 
		    FOFData &D, bool SmoothData)
{

  /* Check if the FOF directory exists */

  char FOF_dir[MAX_LINE_LENGTH];

  if (MyProcessorNumber == ROOT_PROCESSOR && !SmoothData) {
    strcpy(FOF_dir, MetaData->GlobalDir);
    strcat(FOF_dir, "/FOF");
    if (access(FOF_dir, F_OK) == -1)
      SysMkdir("", FOF_dir);
  }

  /* Get enzo units */

  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits;
  double MassUnits=1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, &MassUnits, MetaData->Time);
  
  // Time = a = 1/(1+z).  In enzo, the scale factor is in units of (1+z0).
  FLOAT a = 1, dadt;
  if (ComovingCoordinates) {
    // Mpc/h -> kpc
    D.BoxSize = 1e3 * ComovingBoxSize / HubbleConstantNow;
  
    CosmologyComputeExpansionFactor(MetaData->Time, &a, &dadt);
    D.Time = a / (1 + InitialRedshift);
  }
  else {
    D.BoxSize = LengthUnits / 3.086e21;
    D.Time = 1.0;
  }

  // Critical density in units of Msun / kpc^3
  D.RhoCritical0 = 1.4775867e31 * 
    ((3 * pow(100 * HubbleConstantNow / 3.086e19, 2)) / (8 * M_PI * GRAVITY));
  //D.RhoCritical /= pow(D.Time, 3);

  // Sometimes MassUnits is infinite (in cgs) when using single
  // precision, so calculate it in double precision.

  int TopGridDims3 = MetaData->TopGridDims[0] * MetaData->TopGridDims[1] * 
    MetaData->TopGridDims[2];

  MassUnits = (double) DensityUnits * pow(LengthUnits, 3.0) /
    TopGridDims3;

  // Copy other cosmology parameters
  D.Omega = OmegaMatterNow;
  D.OmegaLambda = OmegaLambdaNow;

  // Get total number of particles to set search radius and softening
  // length for potential computation.  Be careful about nested grid
  // simulations.

  float StaticRegionCellWidth[MAX_STATIC_REGIONS+1];
  int i, level, FinestStaticLevel = 0;
  LevelHierarchyEntry *Temp;
  
  D.NumPart = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      D.NumPart += Temp->GridData->ReturnNumberOfParticles();

  // Search for the finest nested static grid
  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    FinestStaticLevel = max(StaticRefineRegionLevel[i]+1, FinestStaticLevel);

  D.SearchRadius = D.LinkLength * D.BoxSize / pow(TopGridDims3, 1.0/3) / 
    pow(RefineBy, FinestStaticLevel);

  // softening length for potential computation
  D.Epsilon = 0.05 * D.BoxSize / pow(TopGridDims3, 1.0/3) / 
    pow(RefineBy, FinestStaticLevel);

  // Pre-compute cell widths for each static region (for adaptive smoothing)
  StaticRegionCellWidth[0] = D.BoxSize / MetaData->TopGridDims[0];
  for (i = 0; i < MAX_STATIC_REGIONS; i++) {
    if (StaticRefineRegionRightEdge[i][0] > 0)
      StaticRegionCellWidth[i+1] = D.BoxSize / MetaData->TopGridDims[0] /
	pow(RefineBy, StaticRefineRegionLevel[i]+1);
    else
      StaticRegionCellWidth[i+1] = 0;
  }

  if (debug && !SmoothData) {
    fprintf(stdout, "Inline halo finder starting...\n");
    fprintf(stdout, "FOF: Comoving linking length: %g kpc\n", D.SearchRadius);
  }


  /****************** MOVE PARTICLES TO P-GROUPFINDER ******************/

  bool inside;
  int dim, j, proc, slab, region;
  int Index, NumberOfLocalParticles, ptype_size;
  double sr;
  FOF_particle_data *Plocal;
  PINT *Nslab_local, *NtoLeft_local, *NtoRight_local;
  MPI_Arg *MPI_Nslab_local, *MPI_Nslab_recv, *MPI_disp_local, *MPI_disp_recv;

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
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
      Temp->GridData->MoveParticlesFOF(level, Plocal, Index, D,
				       VelocityUnits, MassUnits, COPY_OUT);
      Temp->GridData->SetNumberOfParticles(0);
    }

  /************* Count particles in each slab and shadow ****************/

  D.Nslab = new PINT[NumberOfProcessors];
  D.NtoLeft = new PINT[NumberOfProcessors];
  D.NtoRight = new PINT[NumberOfProcessors];
  D.Nshadow = new PINT[NumberOfProcessors];
  D.Noffset = new PINT[NumberOfProcessors];

  if (NumberOfProcessors == 1) {

    D.Nlocal = NumberOfLocalParticles;
    D.Nslab[0] = NumberOfLocalParticles;
    D.NtoLeft[0] = 0;
    D.NtoRight[0] = 0;
    D.Nshadow[0] = 0;
    D.Noffset[0] = 0;

    allocate_memory(D);
    memcpy(D.P+1, Plocal, sizeof(FOF_particle_data)*NumberOfLocalParticles);
    delete [] Plocal;

  } // ENDIF serial

  else {

    Nslab_local = new PINT[NumberOfProcessors];
    NtoLeft_local = new PINT[NumberOfProcessors];
    NtoRight_local = new PINT[NumberOfProcessors];

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

      /* If we're smoothing the data, have a varying search radius,
	 depending on which static grid the particle is in. */

      if (SmoothData) {

	sr = D.LinkLength * StaticRegionCellWidth[0];
	for (region = MAX_STATIC_REGIONS-1; region >= 0; region--) {
	  if (StaticRefineRegionLevel[region] < 0)
	    continue;
	  inside = true;
	  for (dim = 0; dim < 3; dim++) {
	    inside &= 
	      Plocal[i].Pos[dim] >= StaticRefineRegionLeftEdge[region][dim] &&
	      Plocal[i].Pos[dim] <= StaticRefineRegionRightEdge[region][dim];
	    if (!inside)
	      break;
	  } // ENDFOR dim

	  if (inside) {
	    sr = D.LinkLength * StaticRegionCellWidth[region+1];
	    break;
	  }

	} // ENDFOR region

      } // ENDIF SmoothData
      else
	sr = D.SearchRadius;

      // Left and right "shadows"
      if (Plocal[i].Pos[0] < (slab*(D.BoxSize / NumberOfProcessors) + sr))
	NtoLeft_local[slab]++;

      if (Plocal[i].Pos[0] > ((slab+1)*(D.BoxSize / NumberOfProcessors) - sr))
	NtoRight_local[slab]++;

    } // ENDFOR particles

#ifdef USE_MPI
    // Get counts over all processors
    MPI_Allreduce(Nslab_local, D.Nslab, NumberOfProcessors, PINTDataType,
		  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(NtoLeft_local, D.NtoLeft, NumberOfProcessors, PINTDataType,
		  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(NtoRight_local, D.NtoRight, NumberOfProcessors, PINTDataType,
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
    PINT *disp_local, *disp_recv, *Nslab_recv;

    disp_local = new PINT[NumberOfProcessors];
    disp_recv  = new PINT[NumberOfProcessors];
    Nslab_recv = new PINT[NumberOfProcessors];

    MPI_disp_local = new MPI_Arg[NumberOfProcessors];
    MPI_disp_recv  = new MPI_Arg[NumberOfProcessors];
    MPI_Nslab_local = new MPI_Arg[NumberOfProcessors];
    MPI_Nslab_recv = new MPI_Arg[NumberOfProcessors];

    allocate_memory(D);
    TotalLocal = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++) {
      disp_local[proc] = TotalLocal;
      TotalLocal += Nslab_local[proc];
    }

    // First gather the number of local particles on each processor
#ifdef USE_MPI
    MPI_Alltoall(Nslab_local, 1, PINTDataType, Nslab_recv, 1, PINTDataType,
		 MPI_COMM_WORLD);
#endif

    TotalRecv = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++) {
      disp_recv[proc] = TotalRecv;
      TotalRecv += Nslab_recv[proc];
    }
    D.Nlocal = TotalRecv;

    // Because we're sending the structures as bytes, multiply
    // everything by the size of the particle structure (careful about
    // data types)
    for (proc = 0; proc < NumberOfProcessors; proc++) {
      MPI_Nslab_local[proc] = ptype_size * Nslab_local[proc];
      MPI_Nslab_recv[proc]  = ptype_size * Nslab_recv[proc];
      MPI_disp_local[proc]  = ptype_size * disp_local[proc];
      MPI_disp_recv[proc]   = ptype_size * disp_recv[proc];
    }

    // Now we can do the big collective call
#ifdef USE_MPI
    MPI_Alltoallv(Plocal, MPI_Nslab_local, MPI_disp_local, MPI_BYTE,
		  D.P+1, MPI_Nslab_recv, MPI_disp_recv, MPI_BYTE,
		  MPI_COMM_WORLD);
#endif

    delete [] Plocal;
    delete [] Nslab_local;
    delete [] NtoLeft_local;
    delete [] NtoRight_local;
    delete [] Nslab_recv;
    delete [] disp_local;
    delete [] disp_recv;
    delete [] MPI_disp_local;
    delete [] MPI_disp_recv;
    delete [] MPI_Nslab_local;
    delete [] MPI_Nslab_recv;

  } // ENDELSE multi-processor

}

/***********************************************************************
/
/  COMMUNICATION ROUTINE: MERGE STAR PARTICLE ON THE MAXIMUM LEVEL
/
/  written by: Peng Wang
/  date:       Januaray, 2009
/  modified1:
/
/  PURPOSE:
/    This routine merges star particles on the maximum level within the merging
/    radius.
/
/  DESCRIPTION:
/    All the particles are communicated to every processors. Every processor
/    do the same merger calculation and then delete and merge particle locally.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdio.h>
#include <map>
#include <string>
#include <math.h>
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

int CommunicationBarrier();
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
void ParticleMergeSmallToBig(ParticleEntry *List, const int &Size, 
			     const float &MergeMass, const FLOAT &MergeDistance, 
			     int *Flag, int &GroupSize);
void ParticleMergeSmallGroup(ParticleEntry *List, const int &Size, 
			     const float &MergeMass, const FLOAT &MergeDistance, 
			     int *Flag, int &GroupSize);
void MergeToNewList(ParticleEntry *List, const int &Size, int *Flag, 
		    const int &NewSize, ParticleEntry *NewList);
int CheckMergeFlagList(ParticleEntry *List, const int &Size, int *Flag, const int &GroupSize,
		       const double &MassUnits);
int CommunicationAllSumValues(int *Values, int Number);
int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],
				       int NumberOfGrids);
float ReturnCPUTime();
double ReturnWallTime();

#ifdef USE_MPI
static MPI_Datatype MPI_ParticleEntry;
#endif

int CommunicationMergeStarParticle(HierarchyEntry *Grids[],				   
				   int NumberOfGrids)
{
  //printf("CommunicationMergeStarParticle running......................\n");
#ifdef USE_MPI
  double time1 = ReturnWallTime();

  /* count particles on this processor */

  int i, n, dim, grid;
  Eint32 ParticlesToSend = 0;

  for (grid = 0; grid < NumberOfGrids; grid++) {
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber)
      ParticlesToSend += Grids[grid]->GridData->ReturnNumberOfParticles();
  }

  /* collect particle information on this processor */

  ParticleEntry *SendList = NULL;

  if (ParticlesToSend > 0) {
    
    SendList = new ParticleEntry[ParticlesToSend];

    int c = 0;
    for (grid = 0; grid < NumberOfGrids; grid++) {
      if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
	int np = Grids[grid]->GridData->ReturnParticleEntry(&SendList[c]);
	c += np;
      }
    }

  }

  //for (i = 0; i < count; i++)
  //printf("P(%"ISYM"): %"ISYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" \n", MyProcessorNumber,
  //     ParticlePartialList[i].Number, 
  //     ParticlePartialList[i].Mass, ParticlePartialList[i].Position[0], 
  //     ParticlePartialList[i].Velocity[0], ParticlePartialList[i].Attribute[0]);
	
  static int FirstTimeCalled = TRUE;

  /* define a MPI type for ParticleEntry */

  if (FirstTimeCalled) {
    MPI_Type_contiguous(sizeof(ParticleEntry), MPI_BYTE, &MPI_ParticleEntry);
    MPI_Type_commit(&MPI_ParticleEntry);
    FirstTimeCalled = FALSE;
  }

  /* communicate to get particle counts on every processory */

  Eint32 *SendListCount = new Eint32[NumberOfProcessors];

  MPI_Allgather(&ParticlesToSend, 1, MPI_INT, SendListCount, 1, MPI_INT, MPI_COMM_WORLD);

  int NumberOfSharedParticles = 0;
  for (i = 0; i < NumberOfProcessors; i++)
    NumberOfSharedParticles += SendListCount[i];

  if (NumberOfSharedParticles == 0) {
    delete [] SendListCount;
    return SUCCESS;
  }
  
  /* calculate memory displacement */

  Eint32 *SendListDisplacements = new Eint32[NumberOfProcessors];
  SendListDisplacements[0] = 0;
  for (i = 1; i < NumberOfProcessors; i++)
    SendListDisplacements[i] = SendListDisplacements[i-1] + SendListCount[i-1];

  /* gather all particles to every processors */

  ParticleEntry *SharedList = new ParticleEntry[NumberOfSharedParticles];

  MPI_Allgatherv(SendList, ParticlesToSend, MPI_ParticleEntry,
		 SharedList, SendListCount, SendListDisplacements, MPI_ParticleEntry,
		 MPI_COMM_WORLD);

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, 
    TimeUnits = 1.0, VelocityUnits = 1.0;
  double MassUnits = 1.0;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
	   &VelocityUnits, &MassUnits, 1.0);

  double Msun = 1.989e33;
  if(UsePhysicalUnit) MassUnits = DensityUnits*pow(LengthUnits,3)/Msun;

  /*
  if (debug) {
    for (i = 0; i < NumberOfSharedParticles; i++) {
      int found = 0;
      for (grid = 0; grid < NumberOfGrids; grid++)
	found += Grids[grid]->GridData->CheckGridBoundaries(SharedList[i].Position);
      if (!found) {
	float v = sqrt(pow(SharedList[i].Velocity[0],2) + pow(SharedList[i].Velocity[1],2) +
		       pow(SharedList[i].Velocity[2],2))*VelocityUnits/1e5;
	printf("Particle %"ISYM" is not in any grid!\n", SharedList[i].Number);
	printf("m=%"GSYM", v=%"GSYM", tc=%"GSYM", dm=%"GSYM"\n", SharedList[i].Mass*MassUnits,
	       v, SharedList[i].Attribute[0], SharedList[i].Attribute[2]);
	return FAIL;
      }
    }
  }
  */
  //  printf("\n CMSP: SinkMergeDistance = %"GSYM", SinkMergeMass = %"GSYM", MassUnits = %"GSYM"\n \n",
  //	 SinkMergeDistance, SinkMergeMass, MassUnits); 


  int *MergeFlagList = new int[NumberOfSharedParticles]; // group ID of every merged particle
  for (i = 0; i < NumberOfSharedParticles; i++)
    MergeFlagList[i] = -1; // default is not merged
  int NumberOfGroups = 0;

  /* first, merge small particles to big ones */
  //printf("Merge small particles to big ones \n");
  ParticleMergeSmallToBig(SharedList, NumberOfSharedParticles, 
			  SinkMergeMass/MassUnits, SinkMergeDistance, 
			  MergeFlagList, NumberOfGroups);

  /* second, group small particles using FOF and merge */
  //printf("Merge small particles together \n");
  ParticleMergeSmallGroup(SharedList, NumberOfSharedParticles, 
			  SinkMergeMass/MassUnits, SinkMergeDistance,
			  MergeFlagList, NumberOfGroups);

  /* delete merged old particles */
  //printf("Delete old particles \n");
  for (grid = 0; grid < NumberOfGrids; grid++) {
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      Grids[grid]->GridData->RemoveMergedParticles(SharedList, NumberOfSharedParticles, 
						   MergeFlagList);
      Grids[grid]->GridData->CleanUpMovedParticles();
    }
  }

  /* create a list of merged new particles */

  ParticleEntry *NewList = new ParticleEntry[NumberOfGroups];
  MergeToNewList(SharedList, NumberOfSharedParticles, MergeFlagList, NumberOfGroups, NewList);

  /* add new merged particles to the grids */
  
  int *PartialAdded = new int[NumberOfGroups]; // flag whether a group is merged
  int *TotalAdded = new int[NumberOfGroups];
  for (i = 0; i < NumberOfGroups; i++)
    PartialAdded[i] = 0;

  /* Adjust for periodic boundary conditions */

  float DomainWidth[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
  for (i = 0; i < NumberOfGroups; i++)
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      if (NewList[i].Position[dim] < DomainLeftEdge[dim])
	NewList[i].Position[dim] += DomainWidth[dim];
      else if (NewList[i].Position[dim] > DomainRightEdge[dim])
	NewList[i].Position[dim] -= DomainWidth[dim];
    }

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber)
      Grids[grid]->GridData->AddParticlesFromList(NewList, NumberOfGroups, PartialAdded);

  /* communicate to check whether all the particles are added */

  MPI_Allreduce(PartialAdded, TotalAdded, NumberOfGroups, IntDataType, MPI_SUM, MPI_COMM_WORLD);

  int total = 0;
  for (i = 0; i < NumberOfGroups; i++)
    total += TotalAdded[i];

  int Rank, closest, Dims[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  bool found;
  float r2, r2min, dx[MAX_DIMENSION];
  if (total != NumberOfGroups) {
    if (debug)
      printf("CommunicationMergeParticle: total %"ISYM" != NumberOfGroups %"ISYM"\n", 
	     total, NumberOfGroups);
    for (i = 0; i < NumberOfGroups; i++) {
      if (!TotalAdded[i]) {
	float v = sqrt(pow(NewList[i].Velocity[0],2) + pow(NewList[i].Velocity[1],2) +
		       pow(NewList[i].Velocity[2],2));
	if (debug) {
	  printf("Group %"ISYM" failed to be added.\n", NewList[i].Number);
	  printf("m=%"GSYM", x=(%"GSYM", %"GSYM", %"GSYM"), v=%"GSYM"\n", 
		 NewList[i].Mass*MassUnits, NewList[i].Position[0],
		 NewList[i].Position[1], NewList[i].Position[2], v*VelocityUnits/1e5);
	}
	found = false;
	for (grid = 0; grid < NumberOfGrids; grid++)
	  found |= Grids[grid]->GridData->CheckGridBoundaries(NewList[i].Position);
	if (!found) {
	  if (debug)
	    printf("Inserting merged particle into nearest grid.  RebuildHierarchy will place it\n"
		   "into the correct grid in the next call.\n");
	  closest = 0;
	  r2min = huge_number;
	  for (grid = 0; grid < NumberOfGrids; grid++) {
	    Grids[grid]->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
	    r2 = 0;
	    for (dim = 0; dim < MAX_DIMENSION; dim++) {
	      dx[dim] = 0.5*(Left[dim] + Right[dim]) - NewList[i].Position[dim];
	      r2 += dx[dim] * dx[dim];
	    }
	    if (r2 < r2min) {
	      r2min = r2;
	      closest = grid;
	    }
	  } // ENDFOR grids
	  if (debug) {
	    Grids[closest]->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
	    printf("Adding particle to grid %"ISYM"\n"
		   "left edge  = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n"
		   "right edge = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",		   
		   closest, Left[0], Left[1], Left[2],
		   Right[0], Right[1], Right[2]);
	  }
	  Grids[closest]->GridData->AddOneParticleFromList(NewList, i);
	} // ENDIF !found

	  /* if this is a small problem, ignore it */
//	  if (NewList[i].Mass*MassUnits > 0.1) {
//	    CheckMergeFlagList(SharedList, NumberOfSharedParticles, MergeFlagList, NumberOfGroups, MassUnits);      
//	    return FAIL;
//	  }
      } // ENDIF !TotalAdded
      
      //CheckMergeFlagList(SharedList, NumberOfSharedParticles, MergeFlagList, NumberOfGroups, MassUnits);      
    }
  }

  CommunicationSyncNumberOfParticles(Grids, NumberOfGrids);

  delete [] PartialAdded;
  delete [] TotalAdded;
  delete [] SendList;
  delete [] SendListCount;
  delete [] SendListDisplacements;
  delete [] SharedList;
  delete [] MergeFlagList;
  delete [] NewList;

  //  PerformanceTimers[34] += ReturnWallTime() - time1;
#endif /* MPI */

  return SUCCESS;
}

int CheckMergeFlagList(ParticleEntry *List, const int &Size, int *Flag, const int &GroupSize,
		       const double &MassUnits)
{
  int i, group;
  for (group = 0; group < GroupSize; group++) {
    printf("group %"ISYM" ", group);
    for (i = 0; i < Size; i++) {
      if (Flag[i] == group)
	printf("%"GSYM" ", List[i].Mass*MassUnits);
    }
    printf("\n");
  }

  return 1;
}

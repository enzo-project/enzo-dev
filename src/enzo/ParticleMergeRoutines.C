/*******************************************************************
/
/ FUNCTION: ParticleMergeSmallToBig
/ 
/ written by: Peng Wang
/ date: January, 2009
/
/ PURPOSE: merge small particles to its nearest big ones
/ 
/ INPUT: 
/   List: particle list
/   Size: list size
/   MergeMass: threshold mass for small particles
/   MergeDistance: maximum distance for merging
/
/ OUTPUT:
/   Flag: group id of every input particles. Those won't
/         merge has flag=-1.
/   GroupSize: number of groups
/
**********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"

void ParticleMergeSmallToBig(ParticleEntry *List, const int &Size, 
			     const float &MergeMass, const FLOAT &MergeDistance, 
			     int *Flag, int &GroupSize)
{
  int *Center = new int[Size];
  for (int i = 0; i < Size; i++)
    Center[i] = -1;

  FLOAT MergeDist2 = pow(MergeDistance, 2);
  PINT MergeID, MergeIndex;
  
  for (int i = 0; i < Size; i++) {
    
    if (List[i].Mass > MergeMass) continue;
    
    /* Find the nearest big particle */

    FLOAT MinDist2 = 1e20;
    for (int j = 0; j < Size; j++) {
      if (List[j].Mass <= MergeMass) continue;

      FLOAT dist2 = pow(List[i].Position[0] - List[j].Position[0],2) +
	pow(List[i].Position[1] - List[j].Position[1],2) +
	pow(List[i].Position[2] - List[j].Position[2],2);
      if (dist2 < MinDist2) {
	MinDist2 = dist2;
	MergeID = List[j].Number;
	MergeIndex = j;
      }
    }

    if (MinDist2 < MergeDist2) {

      /* if a group for this MergeID already exist, add to it */

      for (int j = 0; j < i; j++) 
	if (Center[j] == MergeID) {
	  Flag[i] = Flag[j];
	  Center[i] = Center[j];
	}

      /* if no group exists, create a new group */

      if (Center[i] == -1) {
	Center[i] = MergeID; // remember the big one
	Flag[i] = GroupSize; // flag the small one
	Flag[MergeIndex] = GroupSize; // flag the big one
	GroupSize++;
      }

    }
  }

  delete [] Center;

}

/****************************************************************************
/
/ FUNCTION: ParticleMergeSmallGroup
/
/ written by: Peng Wang
/ date: January, 2009
/
/ PURPOSE: merge remaining small particles using Friend-Of-Friend algorithm
/ 
/ INPUT: 
/   List: particle list
/   Size: list size
/   MergeMass: threshold mass for small particles
/   MergeDistance: maximum distance for merging
/
/ OUTPUT:
/   Flag: group id of every input particles. Those won't
/         merge has flag=-1.
/   GroupSize: number of groups
/
********************************************************************************/

/* declarations */

int fof(FLOAT *x, FLOAT *y, FLOAT *z, const int &np, const FLOAT &l,
        int *group, int &ng);

void ParticleMergeSmallGroup(ParticleEntry *List, const int &Size, 
			     const float &MergeMass, const FLOAT &MergeDistance, 
			     int *Flag, int &GroupSize)
{
  //printf("\n Particle Merge Small group MergeDistance = %g\n \n",MergeDistance);
  /* count the number of remaining small particles */

  int NumberOfRemainingSmallParticles = 0;
  for (int i = 0; i < Size; i++)
    if (List[i].Mass < MergeMass && Flag[i] == -1)
      NumberOfRemainingSmallParticles++;

  /* nothing needs to be done, return */

  if (NumberOfRemainingSmallParticles <= 1)
    return;

  /* find out small particles that did not merge to big ones */

  int *IndexArray = new int[NumberOfRemainingSmallParticles];
  int n = 0;
  for (int i = 0; i < Size; i++)
    if (List[i].Mass < MergeMass && Flag[i] == -1)
      IndexArray[n++] = i;

  FLOAT *xp = new FLOAT[NumberOfRemainingSmallParticles];
  FLOAT *yp = new FLOAT[NumberOfRemainingSmallParticles];
  FLOAT *zp = new FLOAT[NumberOfRemainingSmallParticles];

  for (int i = 0; i < NumberOfRemainingSmallParticles; i++) {
    xp[i] = List[IndexArray[i]].Position[0];
    yp[i] = List[IndexArray[i]].Position[1];
    zp[i] = List[IndexArray[i]].Position[2];
  }

  int *TempFlag = new int[NumberOfRemainingSmallParticles];  
  for (int i = 0; i < NumberOfRemainingSmallParticles; i++)
    TempFlag[i] = -1;

  /* call the FOF routine */

  fof(xp, yp, zp, NumberOfRemainingSmallParticles, MergeDistance, TempFlag, GroupSize);

  /* copy the results to the output */

  for (int i = 0; i < NumberOfRemainingSmallParticles; i++)
    if (TempFlag[i] >= 0) 
      Flag[IndexArray[i]] = TempFlag[i];

  delete [] IndexArray;
  delete [] TempFlag;
  delete [] xp;
  delete [] yp;
  delete [] zp;

}

/**************************************************************************
/
/ FUNCTION: MergeToNewList
/ 
/ written by: Peng Wang
/ date: January, 2009
/
/ PURPOSE: Merge flagged particles, creating a list of merged new particles 
/
/ INPUT: List, Size, Flag, NewSize
/   
/ OUTPUT: NewList
/
****************************************************************************/

void MergeToNewList(ParticleEntry *List, const int &Size, int *Flag, 
		    const int &NewSize, ParticleEntry *NewList)
{
  //printf("\n Particle Merge to new list \n \n");
  for (int group = 0; group < NewSize; group++) {

    float mass = 0.0, maxmass = -1.0;
    FLOAT position[3] = {0.0, 0.0, 0.0};
    float velocity[3] = {0.0, 0.0, 0.0};
    float dm = 0.0;
    int index; 
    
    /* transverse the List to calculate the properties of the ith new particle */
    
    for (int j = 0; j < Size; j++) {
      if (Flag[j] == group) {
	dm += List[j].Attribute[2]; // mass increase used for stellar wind feedback
	for (int dim = 0; dim < 3; dim++) {
	  position[dim] = (mass*position[dim] + List[j].Mass*List[j].Position[dim])/
	    (mass + List[j].Mass); // new position is the center of mass
	  velocity[dim] = (mass*velocity[dim] + List[j].Mass*List[j].Velocity[dim])/
	    (mass + List[j].Mass); // momentum conservation
	}
	mass += List[j].Mass;

	if (List[j].Mass > maxmass) { // record the maximum mass one inside a group 
	  maxmass = List[j].Mass;
	  index = j;
	}
      }
    }

    NewList[group].Number = List[index].Number; // ID & attribute is the one with maximum mass
    NewList[group].Type = List[index].Type;
    NewList[group].Mass = mass;
    for (int dim = 0; dim < 3; dim++) {
      NewList[group].Position[dim] = position[dim];
      NewList[group].Velocity[dim] = velocity[dim];
    }
    for (int n = 0; n < NumberOfParticleAttributes; n++)
      NewList[group].Attribute[n] = List[index].Attribute[n];
    NewList[group].Attribute[2] = dm; // reset the second one to what we calculated above

  }
}

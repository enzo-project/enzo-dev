#include <stdlib.h>
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
#include "LevelHierarchy.h"

int loop_count;
void DeleteSourceClusteringTree(SuperSourceEntry * &leaf);
int ReassignSuperSources(LevelHierarchyEntry *LevelArray[]);
void PrintSourceClusteringTree(SuperSourceEntry *leaf, FILE *fptr);

Eint32 compare_x (const void *a, const void *b)
{
  SuperSourceData *ia = (SuperSourceData*) a;
  SuperSourceData *ib = (SuperSourceData*) b;
  if ( ia->Position[0] - ib->Position[0] < 0)
    return -1;
  else if ( ia->Position[0] - ib->Position[0] > 0)
    return 1;
  return 0;
}
Eint32 compare_y (const void *a, const void *b)
{
  SuperSourceData *ia = (SuperSourceData*) a;
  SuperSourceData *ib = (SuperSourceData*) b;
  if ( ia->Position[1] - ib->Position[1] < 0)
    return -1;
  else if ( ia->Position[1] - ib->Position[1] > 0)
    return 1;
  return 0;
}
Eint32 compare_z (const void *a, const void *b)
{
  SuperSourceData *ia = (SuperSourceData*) a;
  SuperSourceData *ib = (SuperSourceData*) b;
  if ( ia->Position[2] - ib->Position[2] < 0)
    return -1;
  else if ( ia->Position[2] - ib->Position[2] > 0)
    return 1;
  return 0;
}

int CreateSourceClusteringTree(int nShine, SuperSourceData *SourceList,
			       LevelHierarchyEntry *LevelArray[])
{

  if (GlobalRadiationSources == NULL)
    return SUCCESS;

  const int LymanWernerBin = 3;

  int i, j, num[2], dim, sort_dim, median, nleft, nright;
  bool top_level = false;
  SuperSourceEntry *new_leaf = NULL;
  SuperSourceData *temp = NULL; // workspace

  /* Create work arrays and delete old tree (if it exists) if this is
     the first time (indicated by SourceList == NULL) */

  if (SourceList == NULL) {
    nShine = 0;
    loop_count = 0;
    top_level = true;
    RadiationSourceEntry *RadSource = GlobalRadiationSources->NextSource;
    while (RadSource != NULL) {
      nShine++;
      RadSource = RadSource->NextSource;
    }
    if (nShine <= 1) 
	return SUCCESS;

    SourceList = new SuperSourceData[nShine];
    RadSource = GlobalRadiationSources->NextSource;
    for (i = 0; i < nShine; i++) {
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	SourceList[i].Position[dim] = RadSource->Position[dim];
      SourceList[i].Luminosity = RadSource->Luminosity;
      if (RadSource->EnergyBins > LymanWernerBin)
	SourceList[i].LWLuminosity = RadSource->Luminosity *
	  RadSource->SED[LymanWernerBin];
      else
	SourceList[i].LWLuminosity = RadSource->LWLuminosity;
      SourceList[i].Source = RadSource;
      RadSource = RadSource->NextSource;
    }

    // Copy clustering tree from previous timestep
    // TODO: Rebuild only branches that have changed.
    if (ReassignSuperSources(LevelArray) == FAIL) {
      ENZO_FAIL("Error in ReassignSuperSources.\n");
    }
    if (OldSourceClusteringTree != NULL)
      DeleteSourceClusteringTree(OldSourceClusteringTree);
    OldSourceClusteringTree = SourceClusteringTree;
    SourceClusteringTree = NULL;

  } // ENDIF SourceList == NULL (first time)  

  /* Calculate "center of light" first and assign it to the tree. */

  FLOAT center[MAX_DIMENSION];
  double weight = 0.0;
  float lw_lum = 0.0;
  
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    center[dim] = 0.0;

  for (i = 0; i < nShine; i++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      center[dim] += SourceList[i].Position[dim] * SourceList[i].Luminosity;
    weight += SourceList[i].Luminosity;
    lw_lum += SourceList[i].LWLuminosity;
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    center[dim] /= weight;

  /* Now calculate the maximum separation of the center and particles */
  
  FLOAT dx;
  float radius2, max_separation = -1e20;
  for (i = 0; i < nShine; i++) {
    radius2 = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      dx = center[dim] - SourceList[i].Position[dim];
      radius2 += dx*dx;
    }
    max_separation = max(max_separation, radius2);
  }
  max_separation = sqrt(max_separation);

  /* Create new leaf and insert into tree */

  new_leaf = new SuperSourceEntry;
  new_leaf->ParentSource = NULL;
  for (i = 0; i < MAX_LEAF; i++)
    new_leaf->ChildSource[i] = NULL;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    new_leaf->Position[dim] = center[dim];
  new_leaf->ClusteringRadius = max_separation;
  new_leaf->LeafID = loop_count;
  new_leaf->LWLuminosity = lw_lum;

  if (SourceClusteringTree == NULL) // top-grid (first time through)
    SourceClusteringTree = new_leaf;
  else
    for (i = 0; i < MAX_LEAF; i++)
      if (SourceClusteringTree->ChildSource[i] == NULL) {
	new_leaf->ParentSource = SourceClusteringTree;
	if (SourceClusteringTree->ClusteringRadius < new_leaf->ClusteringRadius)
	  new_leaf->ClusteringRadius = 0.9 * SourceClusteringTree->ClusteringRadius;
	SourceClusteringTree->ChildSource[i] = new_leaf;
	SourceClusteringTree = new_leaf;
	break;
      }

  /* Assign this leaf to the particle.  It will be overwritten by
     finer levels later since we want the particle to know its finest
     leaf then it can work its way up the hierarchy of super
     sources. */

  for (i = 0; i < nShine; i++)
    SourceList[i].Source->SuperSource = SourceClusteringTree;

  sort_dim = loop_count % MAX_DIMENSION;
  switch (sort_dim) {
  case 0:
    qsort(SourceList, nShine, sizeof(SuperSourceData), compare_x);
    break;
  case 1:
    qsort(SourceList, nShine, sizeof(SuperSourceData), compare_y);
    break;
  case 2:
    qsort(SourceList, nShine, sizeof(SuperSourceData), compare_z);
    break;
  default:
    ENZO_VFAIL("sort_dim = %"ISYM" ?!  This should never be greater than 2.\n",
	    sort_dim)
  } // ENDSWITCH
  loop_count++;

//  printf("%"ISYM" (%"ISYM", %"ISYM") :: %"FSYM" %"FSYM" %"FSYM"\n", 
//	 loop_count-1, sort_dim, nShine,
//	 center[0], center[1], center[2]);
//  for (i = 0; i < nShine; i++)
//    printf("==> %"FSYM" %"FSYM" %"FSYM"\n", SourceList[i].Position[0], 
//	   SourceList[i].Position[1], SourceList[i].Position[2]);

  FLOAT leftdiff, rightdiff;

  if (nShine == 3) {
    leftdiff = SourceList[1].Position[sort_dim] - SourceList[0].Position[sort_dim];
    rightdiff = SourceList[2].Position[sort_dim] - SourceList[1].Position[sort_dim];
    if (rightdiff > leftdiff) {
      median = 1;
      nleft = 2;
      nright = 1;
    } else {
      median = 2;
      nleft = 1;
      nright = 2;
    }
  } else {
    median = nShine/2;
    nleft = (nShine+1)/2;
    nright = nShine-nleft;
  }
  /* Divide into children if there are more than one source */
  
  if (nShine > 2) {

    // Left leaf (copy to temp, make tree, copy back)
    if (nleft > 1) {
      temp = new SuperSourceData[nleft];
      for (i = 0; i < nleft; i++)
	temp[i] = SourceList[i];
      CreateSourceClusteringTree(nleft, temp, NULL);
      for (i = 0; i < nleft; i++)
	SourceList[i] = temp[i];
      SourceClusteringTree = SourceClusteringTree->ParentSource;
      delete [] temp;
    } 

    // Right leaf
    if (nright > 1) {
      temp = new SuperSourceData[nright];
      for (i = 0; i < nright; i++)
	temp[i] = SourceList[nleft+i];
      CreateSourceClusteringTree(nright, temp, NULL);
      for (i = 0; i < nright; i++)
	SourceList[nleft+i] = temp[i];
      SourceClusteringTree = SourceClusteringTree->ParentSource;
      delete [] temp;
    } 
  } // ENDIF nShine > 2

  else {

    sort_dim = SourceClusteringTree->LeafID % MAX_DIMENSION;
    if (nShine > 1) {
      if (SourceList[0].Position[sort_dim] <
	  SourceList[1].Position[sort_dim]) {
	num[0] = 0;
	num[1] = 1;
      } else {
	num[0] = 1;
	num[1] = 0;
      }
    } else {
      num[0] = 0;
    }
    for (i = 0; i < nShine; i++) {
      new_leaf = new SuperSourceEntry;
      new_leaf->ParentSource = NULL;
      for (j = 0; j < MAX_LEAF; j++)
	new_leaf->ChildSource[j] = NULL;
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	new_leaf->Position[dim] = SourceList[i].Position[dim];
      new_leaf->ClusteringRadius = 0;
      new_leaf->LeafID = INT_UNDEFINED;
      new_leaf->LWLuminosity = SourceList[i].LWLuminosity;
      new_leaf->ParentSource = SourceClusteringTree;
      SourceClusteringTree->ChildSource[num[i]] = new_leaf;
    } // ENDFOR i

  }

  if (top_level)
    delete [] SourceList;

  return SUCCESS;

}

void DeleteSourceClusteringTree(SuperSourceEntry * &leaf)
{

  int i;
  for (i = 0; i < MAX_LEAF; i++)
    if (leaf->ChildSource[i] != NULL)
      DeleteSourceClusteringTree(leaf->ChildSource[i]);
  delete leaf;
  leaf = NULL;

}

void PrintSourceClusteringTree(SuperSourceEntry *leaf, FILE *fptr)
{

  int i;

  if (leaf == NULL)
    return;

  fprintf(fptr, 
	 "Source clustering[P%"ISYM"]: leaf %"ISYM", SRC = %x, parent = %x,\n"
	 "                        children = %x %x\n"
	 "                        pos = %"FSYM" %"FSYM" %"FSYM", cradius = %"GSYM"\n",
	 MyProcessorNumber, leaf->LeafID, leaf, leaf->ParentSource,
	 leaf->ChildSource[0], leaf->ChildSource[1], leaf->Position[0],
	 leaf->Position[1], leaf->Position[2], leaf->ClusteringRadius);

  for (i = 0; i < MAX_LEAF; i++)
    if (leaf->ChildSource[i] != NULL)
      PrintSourceClusteringTree(leaf->ChildSource[i], fptr);

}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode)
{

  NewNode->PreviousPackage = Node;
  NewNode->NextPackage = Node->NextPackage;
  if (Node->NextPackage != NULL)
    Node->NextPackage->PreviousPackage = NewNode;
  Node->NextPackage = NewNode;
  return;

}

PhotonPackageEntry *PopPhoton(PhotonPackageEntry * &Node)
{

  if (Node->PreviousPackage != NULL)
    Node->PreviousPackage->NextPackage = Node->NextPackage;
  if (Node->NextPackage != NULL)
    Node->NextPackage->PreviousPackage = Node->PreviousPackage;
  return Node;


}

PhotonPackageEntry *LinkedListToArray(PhotonPackageEntry *Node, int n)
{
  int dim, count = 0;
  PhotonPackageEntry *result = new PhotonPackageEntry[n];
  PhotonPackageEntry *tmp = Node;
  while (tmp != NULL) {
    result[count].NextPackage = tmp->NextPackage;
    result[count].PreviousPackage = tmp->PreviousPackage;
    result[count].CurrentSource = tmp->CurrentSource;
    //result[count].CurrentSource->ParentSource = tmp->CurrentSource;
    result[count].Photons = tmp->Photons;
    result[count].Type = tmp->Type;
    result[count].Energy = tmp->Energy;
    result[count].CrossSection = tmp->CrossSection;
    result[count].EmissionTimeInterval = tmp->EmissionTimeInterval;
    result[count].EmissionTime = tmp->EmissionTime;
    result[count].CurrentTime = tmp->CurrentTime;
    result[count].Radius = tmp->Radius;
    result[count].ColumnDensity = tmp->ColumnDensity;
    result[count].ipix = tmp->ipix;
    result[count].level = tmp->level;
    result[count].SourcePositionDiff = tmp->SourcePositionDiff;
    for (dim = 0; dim < 3; dim++)
      result[count].SourcePosition[dim] = tmp->SourcePosition[dim];
    count++;
    tmp = tmp->NextPackage;
  }
  return result;
}

#ifdef UNUSED
/**********************************************************************/
/*                    QUICKSORT FOR LINKED LISTS                      */
/**********************************************************************/

// Insert after a known node and return the node before inserted node
PhotonPackageEntry *QSinsert(PhotonPackageEntry *Node, 
			     PhotonPackageEntry *NewNode)
{
  NewNode->PreviousPackage = Node;
  NewNode->NextPackage = Node->NextPackage;
  if (Node->NextPackage != NULL)
    Node->NextPackage->PreviousPackage = NewNode;
  Node->NextPackage = NewNode;
  return Node;
}

void QSdestroylist(PhotonPackageEntry **Node)
{
  PhotonPackageEntry *tmp;
  if (!Node || !*Node) return;
  while (*Node) {
    tmp = (*Node)->NextPackage;
    delete *Node;
    *Node = tmp;    
  }
  return;
}

// Given a node in list, prints entire list after the node
void QSprintlist(PhotonPackageEntry *Node)
{
  PhotonPackageEntry *t = Node;
  while (t) {
    printf("Photon %x :: level %"ISYM", ipix %"ISYM"\n", t, t->level, t->ipix);
    t = t->NextPackage;
  }
  return;
}

// Find previous node of curr node
PhotonPackageEntry *QSfindPrev(PhotonPackageEntry *list, 
			       PhotonPackageEntry *curr)
{
  for(;list && list->NextPackage; list=list->NextPackage)
    if (curr == list->NextPackage)
      break;
  // Didn't find any element because we were searching for the
  // beginning node
  if (list->NextPackage != curr)
    return curr;
  return list;
}

/* A complete swap algorithm which cares of 
several scenarios while swapping two nodes in 
a linked list which doesn't have any special nodes
scenarios considered while swapping

1) two nodes which are far away

2) two nodes which are far away, one is node is at beginning of the
   list

3) two node which are neighbors

4) two nodes which are neibhors, one node is at beginning of the list

*/
#endif

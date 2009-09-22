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

void InsertStarAfter(Star * &Node, Star * &NewNode)
{
  if (Node == NULL)
    Node = NewNode;
  else {
    if (NewNode == Node->NextStar) {
      printf("Node already in list?!\n");
      exit(1);
    }
    NewNode->PrevStar = Node;
    NewNode->NextStar = Node->NextStar;
    if (Node->NextStar != NULL)
      Node->NextStar->PrevStar = NewNode;
    Node->NextStar = NewNode;
  }
  return;
}

Star *PopStar(Star * &Node)
{
  Star *result = Node;
  if (Node->PrevStar != NULL)
    Node->PrevStar->NextStar = Node->NextStar;
  if (Node->NextStar != NULL)
    Node->NextStar->PrevStar = Node->PrevStar;
  Node = Node->NextStar;
  result->NextStar = NULL;
  result->PrevStar = NULL;
  return result;
}

void DeleteStar(Star * &Node)
{
  Star *Orphan = PopStar(Node);
  //Node = Node->NextStar;
  if (Orphan != NULL) delete Orphan;
  return;
}

void DeleteStarList(Star * &Node)
{
  Star *tmp = Node;
  while (tmp)  // delete all linked stars
    DeleteStar(tmp);
  Node = NULL;
  return;
}

Star *StarListToArray(Star *Node, int n)
{
  int dim, count = 0;
  Star *result = new Star[n];
  Star *tmp = Node;
  while (tmp != NULL) {
    result[count++] = *tmp;
    tmp = tmp->NextStar;
  }
  return result;
}

/* Since InsertStarAfter puts the node after the head node.  We insert
   the nodes in a fashion to preserve the order of the array. */

Star* StarBufferToList(StarBuffer buffer)
{
  Star *result = NULL;
  result = new Star(buffer);
  return result;
}

Star* StarBufferToList(StarBuffer *buffer, int n) 
{
  int i;
  Star *result = NULL, *NewNode = NULL;
  if (n > 0) {
    NewNode = new Star(buffer, 0);
    InsertStarAfter(result, NewNode);
  }
  for (i = n-1; i > 0; i--) {
    NewNode = new Star(buffer, i);
    InsertStarAfter(result, NewNode);
  }
  return result;
}

#define DEBUG 0
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "PhotonPackage.h"

/* Find super source leaf in the binary tree, given an ID. */

void PrintSourceClusteringTree(SuperSourceEntry *leaf, FILE *fptr);

int FindSuperSource(PhotonPackageEntry **PP, int &LeafID, 
		    int SearchNewTree = TRUE)
{

  // If no super source (top-level), assign CurrentSource = NULL
  if (LeafID < 0) {
    (*PP)->CurrentSource = NULL;
    return SUCCESS;
  }

  int i;
  SuperSourceEntry *temp, *last = NULL;
  temp = (SearchNewTree) ? SourceClusteringTree : OldSourceClusteringTree;

  while (1) {

    if (temp == NULL) {
      if (DEBUG) {
	FILE *fptr = fopen("sources.dat", "w");
	PrintSourceClusteringTree(SourceClusteringTree, fptr);
	fclose(fptr);
	fprintf(stdout, "WARNING: NULL leaf in clustering tree. Source was deleted? "
		"LeafID = %"ISYM", temp = %p, last = %p\n", LeafID, temp, last);
	if (last != NULL)
	  fprintf(stdout, "\t last->leafID = %d\n", last->LeafID);
      }
      // Will search by position in grid::ReassignSuperSources
      LeafID = INT_UNDEFINED;
      (*PP)->CurrentSource = NULL;
      break;
      //ENZO_FAIL("");
    }

    // Found it!
    if (LeafID == temp->LeafID) {
      (*PP)->CurrentSource = temp;
      break;
    }

    /* Tree is numbered like so...
            0
       1         6
     2   5     7    10
    3 4       8 9 11
     */

    // Check which branch to transverse (general for any MAX_LEAF)
    for (i = 0; i < MAX_LEAF-1; i++) {
      last = temp;
      if (i == MAX_LEAF-2 && temp->ChildSource[i+1] == NULL) {
	// next child doesn't exist, it has to be this child
	temp = temp->ChildSource[i];
      } else {

	// Check if the ID is between this and the next child
	if (LeafID >= temp->ChildSource[i]->LeafID &&
	    LeafID < temp->ChildSource[i+1]->LeafID) {
	  temp = temp->ChildSource[i];  // this child

	// IF next leaf is the last leaf, this is the only choice left
	} else if (i == MAX_LEAF-2) {
	  temp = temp->ChildSource[i+1];
	}

      } // ENDELSE
    } // ENDFOR

  } // ENDWHILE
  
  return SUCCESS;

}

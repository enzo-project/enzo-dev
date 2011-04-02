#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "PhotonPackage.h"

/* Find super source leaf in the binary tree, given a photon package.
   Search by position and photon radius. */

/* Travel down the tree until a child leaf is found with a clustering
   radius smaller than the merging radius of the photon. */

/* ONLY VALID FOR MAX_LEAF=2! (i.e. binary tree) */

int FindSuperSourceByPosition(PhotonPackageEntry **PP)
{

  int i, dim_search, loop_count = 0, found = FALSE;
  SuperSourceEntry *temp = OldSourceClusteringTree;
  float SearchRadius;  // min clustering search radius
  SearchRadius = 0.5*(*PP)->Radius / RadiativeTransferPhotonMergeRadius;

  while (1) {

    if (temp == NULL) {
      ENZO_VFAIL("FindSuperSourceByPosition: NULL leaf in clustering tree.  "
	      "This shouldn't happen. LeafID = %"ISYM"\n", temp->LeafID)
    }

    //dim_search = loop_count % MAX_DIMENSION;
    dim_search = temp->LeafID % MAX_DIMENSION;

    // Left leaf
    if ((*PP)->SourcePosition[dim_search] < temp->Position[dim_search]) {
      if (temp->ChildSource[0] == NULL)
	break;
      if (temp->ChildSource[0]->ClusteringRadius < SearchRadius)
	break;
      temp = temp->ChildSource[0];
    } // ENDIF left leaf

    // Right leaf
    else {
      if (temp->ChildSource[1] == NULL)
	break;
      if (temp->ChildSource[1]->ClusteringRadius < SearchRadius)

	break;
      temp = temp->ChildSource[1];
    } // ENDELSE right left
    
    loop_count++;

  } // ENDWHILE

  (*PP)->CurrentSource = temp;
  
  return SUCCESS;

}

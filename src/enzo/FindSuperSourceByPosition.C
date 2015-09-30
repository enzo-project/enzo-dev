#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <xmmintrin.h>
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

/**********************************************************************/
// Overloaded function with a position instead of a photon package as
// the input.  Using the new clustering tree.

int FindSuperSourceByPosition(FLOAT *pos, SuperSourceEntry **result,
			      int DEBUG)
{

  int i, dim, dim_search, loop_count = 0;
  SuperSourceEntry *temp = SourceClusteringTree;
  float merge_inv = 0.5 / RadiativeTransferPhotonMergeRadius;
  FLOAT SearchRadius;  // min clustering search radius
  FLOAT dx;
  //SearchRadius = 0.5*(*PP)->Radius / RadiativeTransferPhotonMergeRadius;

  while (1) {

    if (temp == NULL) {
      ENZO_VFAIL("FindSuperSourceByPosition: NULL leaf in clustering tree.  "
	      "This shouldn't happen. LeafID = %"ISYM"\n", temp->LeafID)
    }

    //dim_search = loop_count % MAX_DIMENSION;
    if (temp->LeafID >= 0)
      dim_search = temp->LeafID % MAX_DIMENSION;
    else
      dim_search = temp->ParentSource->LeafID % MAX_DIMENSION;

    if (DEBUG)
      printf("leaf %d: dim = %d, pos = %f %f %f, leafpos = %f %f %f\n",
	     temp->LeafID, dim_search, pos[0], pos[1], pos[2],
	     temp->Position[0], temp->Position[1],
	     temp->Position[2]);

    // Left leaf
    if (pos[dim_search] < temp->Position[dim_search]) {
      if (DEBUG) printf("--> Picking left leaf\n");
      if (temp->ChildSource[0] == NULL)
	break;
      // Compute distance from position to the first child
      SearchRadius = 0.0;
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	dx = pos[dim] - temp->ChildSource[0]->Position[dim];
	SearchRadius += dx*dx;
      }
      SearchRadius = merge_inv * sqrt(SearchRadius);
      if (DEBUG)
	printf("->> cluster_radius = %g, search = %g\n",
	       temp->ChildSource[0]->ClusteringRadius, SearchRadius);
      if (temp->ChildSource[0]->ClusteringRadius < SearchRadius &&
	  temp->ChildSource[0]->ClusteringRadius > 0)
	break;
      temp = temp->ChildSource[0];
    } // ENDIF left leaf

    // Right leaf
    else {
      if (DEBUG) printf("--> Picking right leaf\n");
      if (temp->ChildSource[1] == NULL)
	break;
      // Compute distance from position to the second child
      SearchRadius = 0.0;
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	dx = pos[dim] - temp->ChildSource[1]->Position[dim];
	SearchRadius += dx*dx;
      }
      SearchRadius = merge_inv * sqrt(SearchRadius);
      if (DEBUG)
	printf("->> cluster_radius = %g, search = %g\n",
	       temp->ChildSource[1]->ClusteringRadius, SearchRadius);
      if (temp->ChildSource[1]->ClusteringRadius < SearchRadius &&
	  temp->ChildSource[1]->ClusteringRadius > 0)
	break;
      temp = temp->ChildSource[1];
    } // ENDELSE right left
    
    loop_count++;

  } // ENDWHILE

  *result = temp;
  
  return SUCCESS;

}

/* SSE intrinsic approximate inverse sqrt.  IEEE precision isn't
   required to choose the correct leafs. */

inline void vrsqrt(Eflt32* __x, Eflt32* __outrsqrt)
{
  __m128 x = _mm_set_ss(*__x);
  __m128 recip = _mm_rsqrt_ss(x);
  _mm_store_ss(__outrsqrt, recip);
// _m128* precip = (__m128 *)__outrsqrt;
// *precip = _mm_mul_ss(_mm_set_ss(0.5f), _mm_add_ss(recip, _mm_rcp_ss(_mm_mul_ss(x, recip))));
}

float CalculateLWFromTree(const FLOAT pos[], 
			  const float angle, 
			  const SuperSourceEntry *Leaf, 
			  const float min_radius, 
			  float result0)
{

  int dim;
  FLOAT dx, radius2;
  float tan_angle, result;
  Eflt32 temp, radius_inv;

  if (Leaf == NULL) 
    return result0;

  result = result0;
  radius2 = 0.0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    dx = Leaf->Position[dim] - pos[dim];
    radius2 += dx*dx;
  }

  temp = (Eflt32)radius2;
  temp = max(min_radius, temp);
  vrsqrt(&temp, &radius_inv);
  //radius_inv = 1.0 / sqrtf((float)radius2);
  tan_angle = Leaf->ClusteringRadius * radius_inv;

//  int pid = (Leaf->ParentSource == NULL) ? -1 : Leaf->ParentSource->LeafID;
//  printf("Leaf->ID = %d (%d), cradius = %g, radius = %g, tan_angle = %g, result0 = %g\n",
//	 Leaf->LeafID, pid, Leaf->ClusteringRadius, sqrt(radius2), tan_angle, result0);

  // Larger than opening angle -> go to children
  if (tan_angle > angle) {
    result = CalculateLWFromTree(pos, angle, Leaf->ChildSource[0], min_radius, result);
    result = CalculateLWFromTree(pos, angle, Leaf->ChildSource[1], min_radius, result);
  }

  // Smaller than opening angle -> use this in the calculation
  else {
    result += Leaf->LWLuminosity * radius_inv * radius_inv;
  }

  //printf("\t after[%d] -- result = %g\n", Leaf->LeafID, result);

  return result;

}

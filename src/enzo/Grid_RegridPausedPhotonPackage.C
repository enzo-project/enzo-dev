#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (RE-GRID PAUSED PHOTON PACKAGES AND MOVE TO OTHER GRID 
/              IF NECESSARY)
/
/  written by: John Wise
/  date:       May, 2011
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "SortCompareFunctions.h"

#define MAX_HEALPIX_LEVEL 13

int grid::RegridPausedPhotonPackage(PhotonPackageEntry** PP, grid* ParentGrid,
				    grid** MoveToGrid, int &DeltaLevel,
				    int &DeleteMe, const float *DomainWidth,
				    const float LightSpeed)
{

  if ((*PP) == NULL) {
    ENZO_VFAIL("Called grid::RegridPausedPhotonPackage with an invalid pointer.\n"
	    "\t %p %p %p %p\n",
	    (*PP), (*PP)->PreviousPackage, (*PP)->PreviousPackage->NextPackage,
	    PhotonPackages)
  }

  /* Reassign source position and recalculate HEALPix pixel number
     based on the normal vector between the ray and new super
     source.  Assign super source to the parent. */

  int dim, index;
  int int_pos[MAX_DIMENSION];
  float length;
  FLOAT new_pos[MAX_DIMENSION];
  FLOAT original_vec[MAX_DIMENSION], vec[MAX_DIMENSION], 
    new_vec[MAX_DIMENSION];
  float dx_inv = 1.0 / this->CellWidth[0][0];
  float dx2_inv = dx_inv * dx_inv;
  const float ln2_inv = 1.0/M_LN2;

  // Calculate original unit directional vector
  if (pix2vec_nest((long) (1 << (*PP)->level), (*PP)->ipix, original_vec) == FAIL) {
    ENZO_VFAIL("grid::RegridPausedPhotonPackages(1) -- "
	       "pix2vec_nest %"ISYM" %"ISYM" %"GSYM"\n",
	       (long) (1 << (*PP)->level), (*PP)->ipix, (*PP)->Photons)
  }

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    new_pos[dim] = (*PP)->SourcePosition[dim] + original_vec[dim]*(*PP)->Radius;

  length = 0.0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    // Unnormalized vector from super source to this package
    vec[dim] = ((*PP)->SourcePosition[dim] + (*PP)->Radius * original_vec[dim])
      - (*PP)->CurrentSource->Position[dim];
    length += vec[dim] * vec[dim];
	
    // Change source to super source.
    (*PP)->SourcePosition[dim] = (*PP)->CurrentSource->Position[dim];
  } // ENDFOR dim
  (*PP)->SourcePositionDiff = 0.0;

//  printf("before %p: lvl %"ISYM" pix %"ISYM" :: r=%"GSYM", "
//	 "x=%"GSYM" %"GSYM" %"GSYM"\n", 
//	 (*PP), (*PP)->level, (*PP)->ipix, (*PP)->Radius,
//	 new_pos[0], new_pos[1], new_pos[2]);
  
  length = sqrt(length);
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    vec[dim] /= length;

  // With the new radius, calculate new HEALPix level
  (*PP)->Radius = length;
  (*PP)->level = (int) (0.5*ln2_inv * 
		     logf(3 * M_1_PI * ((*PP)->Radius*(*PP)->Radius * dx2_inv)));
  (*PP)->level = min(max((*PP)->level, 0), MAX_HEALPIX_LEVEL);

  // Adjust CurrentTime to equal (Radius / c)
  (*PP)->CurrentTime = length / LightSpeed;

  // Calculate new pixel number with the super source
  if (vec2pix_nest( (long) (1 << (*PP)->level), vec, &((*PP)->ipix) ) == FAIL) {
    ENZO_VFAIL("grid::RegridPausedPhotonPackages -- "
	       "vec2pix_nest %"ISYM" %"ISYM" %"GSYM"\n",
	       (long) (1 << (*PP)->level), (*PP)->ipix, (*PP)->Photons)
  }
  
  // Calculate new unit directional vector for merged package
  if (pix2vec_nest((long) (1 << (*PP)->level), (*PP)->ipix, new_vec) == FAIL) {
    ENZO_VFAIL("grid::RegridPausedPhotonPackages(2) -- "
	       "pix2vec_nest %"ISYM" %"ISYM" %"GSYM"\n",
	       (long) (1 << (*PP)->level), (*PP)->ipix, (*PP)->Photons)
  }

  // Calculate new photon package and see if it needs to be moved.
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    new_pos[dim] = (*PP)->SourcePosition[dim] + new_vec[dim]*(*PP)->Radius;
    int_pos[dim] = (int) ((new_pos[dim] - this->GridLeftEdge[dim]) * dx_inv + 
			  GridStartIndex[dim]);
  }
//  printf("after %p:  lvl %"ISYM" pix %"ISYM" :: r=%"GSYM", "
//	 "x=%"GSYM" %"GSYM" %"GSYM"\n", 
//	 (*PP), (*PP)->level, (*PP)->ipix, (*PP)->Radius,
//	 new_pos[0], new_pos[1], new_pos[2]);
    
  index = GRIDINDEX_NOGHOST(int_pos[0], int_pos[1], int_pos[2]);
  if (SubgridMarker[index] != this)
    this->FindPhotonNewGrid(index, new_pos, new_vec, *PP, *MoveToGrid,
			    DeltaLevel, DomainWidth, DeleteMe, 
			    ParentGrid);

  return SUCCESS;

}

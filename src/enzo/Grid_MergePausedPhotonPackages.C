#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (MERGE PAUSED PHOTON PACKAGES AND PUT BACK INTO MAIN LIST)
/
/  written by: John Wise
/  date:       September, 2008
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

PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);
PhotonPackageEntry *LinkedListToArray(PhotonPackageEntry *Node, int n);

Eint32 compare_pix (const void *a, const void *b)
{
  PhotonPackageEntry *ia = (PhotonPackageEntry*) a;
  PhotonPackageEntry *ib = (PhotonPackageEntry*) b;
  if ( ia->ipix < ib->ipix)
    return -1;
  else if ( ia->ipix > ib->ipix)
    return 1;
  return 0;
}
Eint32 compare_lvl (const void *a, const void *b)
{
  PhotonPackageEntry *ia = (PhotonPackageEntry*) a;
  PhotonPackageEntry *ib = (PhotonPackageEntry*) b;
  if ( ia->level < ib->level)
    return -1;
  else if ( ia->level > ib->level)
    return 1;
  return 0;
}
Eint32 compare_ss (const void *a, const void *b)
{
  PhotonPackageEntry *ia = (PhotonPackageEntry*) a;
  PhotonPackageEntry *ib = (PhotonPackageEntry*) b;
  // Sort by source, then level, then pixel number, then photon type
  if ( ia->CurrentSource < ib->CurrentSource )
    return -1;
  else if ( ia->CurrentSource > ib->CurrentSource )
    return 1;
  else {
    if ( ia->level < ib->level)
      return -1;
    else if ( ia->level > ib->level)
      return 1;
    else {
      if ( ia->ipix < ib->ipix)
	return -1;
      else if ( ia->ipix > ib->ipix)
	return 1;
      else {
	if ( ia->Type < ib->Type)
	  return -1;
	else if ( ia->Type > ib->Type)
	  return 1;
	else
	  return 0;
      }
    }
  }
}

int grid::MergePausedPhotonPackages() {

  if (PausedPhotonPackages->NextPackage == NULL)
    return 0;

  /* Reassign source position and recalculate HEALPix pixel number
     based on the normal vector between the ray and new super
     source.  Assign super source to the parent. */
  
  int i, dim, nphotons, count;
  FLOAT length, length2, length_inv;
  FLOAT original_vec[MAX_DIMENSION], vec[MAX_DIMENSION];
  PhotonPackageEntry *PP = PausedPhotonPackages->NextPackage;
  float dx2 = this->CellWidth[0][0] * this->CellWidth[0][0];
  const float ln2_inv = 1.0/M_LN2;
  
  nphotons = 0;
  while (PP != NULL) {

    // Calculate original unit directional vector
    if (pix2vec_nest((long) (1 << PP->level), PP->ipix, original_vec) == FAIL) {
      ENZO_VFAIL("grid::MergePausedPhotonPackages -- pix2vec_nest %"ISYM" %"ISYM" %"GSYM"\n",
	      (long) (1 << PP->level), PP->ipix, PP->Photons)
    }

    length2 = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      // Unnormalized vector from super source to this package
      vec[dim] = (PP->SourcePosition[dim] + PP->Radius * original_vec[dim])
	- PP->CurrentSource->Position[dim];
      length2 += vec[dim] * vec[dim];
	
      // Change source to super source.
      PP->SourcePosition[dim] = PP->CurrentSource->Position[dim];
    } // ENDFOR dim
    PP->SourcePositionDiff = 0.0;

    //      printf("before %x: lvl %"ISYM" pix %"ISYM" :: r=%"GSYM", u=%"GSYM" %"GSYM" %"GSYM"\n", 
    //	     PP, PP->level, PP->ipix, PP->Radius,
    //	     original_vec[0], original_vec[1], original_vec[2]);

    length = sqrt(length2);
    length_inv = 1.0/length;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      vec[dim] *= length_inv;

    // With the new radius, calculate new HEALPix level
    PP->Radius = length2;
    PP->level = (int) (0.5*ln2_inv * 
		       logf(3 * M_1_PI * (PP->Radius*PP->Radius/dx2) * 
			    RadiativeTransferRaysPerCell));
    PP->level = min(max(PP->level, 0), MAX_HEALPIX_LEVEL);

    // Calculate new pixel number with the super source
    if (vec2pix_nest( (long) (1 << PP->level), vec, &(PP->ipix) ) == FAIL) {
      ENZO_VFAIL("grid::MergePausedPhotonPackages -- vec2pix_nest %"ISYM" %"ISYM" %"GSYM"\n",
	      (long) (1 << PP->level), PP->ipix, PP->Photons)
    }
    //      printf("after %x:  lvl %"ISYM" pix %"ISYM" :: r=%"GSYM", u=%"GSYM" %"GSYM" %"GSYM"\n", 
    //	     PP, PP->level, PP->ipix, PP->Radius,
    //	     vec[0], vec[1], vec[2]);

    nphotons++;
    PP = PP->NextPackage;

  } // ENDWHILE photons

  if (nphotons < 2)
    return SUCCESS;

  /* It's easier to sort an array with qsort rather than a linked
     list, so let's put the photons in a temp. array.  After sorting
     we should reassign the links. */

  int min_level = 99999, max_level = -99999;
  PhotonPackageEntry *TempPP;
  int ilvl, level, cumulative_sum1, cumulative_sum2;
  int photons_this_level, photons_this_source;
  int count_ss, count_lvl;

  PP = PausedPhotonPackages->NextPackage;
  TempPP = LinkedListToArray(PP, nphotons);

  /* Sort by super source, then on level for each source, then on
     pixel number, last on photon type. */

  if (DEBUG) {
    printf("========== BEFORE SORTING ==========\n");
    for (i = 0; i < nphotons; i++)
      printf("photon %"ISYM": type %"ISYM", lvl %"ISYM", pix %"ISYM", r=%"GSYM", L=%"GSYM", CSRC=%x\n", i, TempPP[i].Type, TempPP[i].level,
	     TempPP[i].ipix, TempPP[i].Radius, TempPP[i].Photons, 
	     TempPP[i].CurrentSource);
  }

  qsort(TempPP, nphotons, sizeof(PhotonPackageEntry), compare_ss);
  //std::sort(TempPP, TempPP+nphotons, cmp_ss());

  if (DEBUG) {
    printf("========== AFTER ALL SORTING ==========\n");
    for (i = 0; i < nphotons; i++)
      printf("photon %"ISYM": type %"ISYM", lvl %"ISYM", pix %"ISYM", r=%"GSYM", L=%"GSYM", CSRC=%x\n", i, TempPP[i].Type, TempPP[i].level,
	     TempPP[i].ipix, TempPP[i].Radius, TempPP[i].Photons, 
	     TempPP[i].CurrentSource);
  }

  /* Now that the list is sorted, we can easily merge pixels with the
     same pixel number, level, and type */

  if (DEBUG)
    printf("========== AFTER MERGING ==========\n");

  PhotonPackageEntry *NewPack = NULL;
  int match, merges = 0;
  float weight;
  for (i = 0; i < nphotons; i++) {

//    if (TempPP[i].Photons <= tiny_number)
//      continue;

    if (i > 0)
      match = ((TempPP[i].level == TempPP[i-1].level) &&
	       (TempPP[i].ipix  == TempPP[i-1].ipix)  &&
	       (TempPP[i].Type  == TempPP[i-1].Type) &&
	       (TempPP[i].CurrentSource == TempPP[i-1].CurrentSource));
    else
      match = FALSE;

    if (match) {
      // Add photons to the current package
      weight = TempPP[i].Photons;
      NewPack->Photons += weight;
      NewPack->EmissionTimeInterval += TempPP[i].EmissionTimeInterval * weight;
      //NewPack->Radius += TempPP[i].Radius * weight;
      NewPack->ColumnDensity += TempPP[i].ColumnDensity * weight;
      this->NumberOfPhotonPackages--;
    } else { // ENDIF match

      // First put the previous package (after correcting several
      // values for weighted averages) in the linked list if not NULL
      if (NewPack != NULL) {
	//NewPack->Radius /= NewPack->Photons;
	NewPack->EmissionTimeInterval /= NewPack->Photons;
	NewPack->ColumnDensity /= NewPack->Photons;
	if (DEBUG)
	  printf("photon %"ISYM": type %"ISYM", lvl %"ISYM", pix %"ISYM", r=%"GSYM", L=%"GSYM", CSRC=%x\n", merges, 
		 NewPack->Type, NewPack->level,
		 NewPack->ipix, NewPack->Radius, NewPack->Photons, 
		 NewPack->CurrentSource);
	InsertPhotonAfter(this->PhotonPackages, NewPack);
      }

      // Create a new package
      NewPack = new PhotonPackageEntry;
      weight = TempPP[i].Photons;
      NewPack->Photons = TempPP[i].Photons;
      NewPack->Type = TempPP[i].Type;
      NewPack->Energy = TempPP[i].Energy;
      NewPack->CrossSection = TempPP[i].CrossSection;
      NewPack->EmissionTimeInterval = TempPP[i].EmissionTimeInterval * weight;
      NewPack->EmissionTime = TempPP[i].EmissionTime;
      NewPack->CurrentTime = TempPP[i].CurrentTime;
      NewPack->Radius = TempPP[i].Radius;
      NewPack->ColumnDensity = TempPP[i].ColumnDensity * weight;
      NewPack->ipix = TempPP[i].ipix;
      NewPack->level = TempPP[i].level;
      NewPack->SourcePositionDiff = 0.0;
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	NewPack->SourcePosition[dim] = TempPP[i].SourcePosition[dim];
      
      // go up one level up in the source tree
      NewPack->CurrentSource = TempPP[i].CurrentSource->ParentSource;

      merges++;
    } // ENDELSE match
  } // ENDFOR photons
  
  // Insert the last ray into the linked list
  if (NewPack != NULL) {
    //NewPack->Radius /= NewPack->Photons;
    NewPack->EmissionTimeInterval /= NewPack->Photons;
    NewPack->ColumnDensity /= NewPack->Photons;
    if (DEBUG)
      printf("photon %"ISYM": type %"ISYM", lvl %"ISYM", pix %"ISYM", r=%"GSYM", L=%"GSYM", CSRC=%x\n", merges, 
	     NewPack->Type, NewPack->level,
	     NewPack->ipix, NewPack->Radius, NewPack->Photons, 
	     NewPack->CurrentSource);
    InsertPhotonAfter(this->PhotonPackages, NewPack);
  }

  /* Delete all paused packages and cleanup temporary arrays */

  PP = PausedPhotonPackages->NextPackage;
  while (PP != NULL) {
    PP = DeletePhotonPackage(PP);
    PP = PP->NextPackage;
  }
  PausedPhotonPackages->NextPackage = NULL;
  PausedPhotonPackages->PreviousPackage = NULL;
  delete [] TempPP;

  if (DEBUG)
    printf("P%d: MergePausedPhotonPackages: %"ISYM" => %"ISYM" photons\n", 
	   MyProcessorNumber, nphotons, merges);

  return SUCCESS;
}

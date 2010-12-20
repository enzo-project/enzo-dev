/***********************************************************************
/
/  MERGES ALL CLUSTERED NEW STAR PARTICLES
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  NOTES:
/  (March 2009) For now, we will retain the old and inefficient O(N^2)
/  method of finding new star particle neighbors and merging.  In a
/  later revision, we should use FOF or HOP or something else to merge
/  star particles.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
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

Star *PopStar(Star * &Node);
void InsertStarAfter(Star * &Node, Star * &NewNode);
void DeleteStar(Star * &Node);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int StarParticleMergeNew(LevelHierarchyEntry *LevelArray[], Star *&AllStars)
{

  Star *ThisStar, *OtherStar, *LastStar, *MoveStar;
  LevelHierarchyEntry *Temp;
  float rmerge2, rmerge2o, dx, dx2;
  FLOAT TimeNow;
  int dim, level;
  const float pc = 3.086e18;

  /* Get the time at the finest level */
  
  for (level = MAX_DEPTH_OF_HIERARCHY-1; level >= 0; level--)
    if (LevelArray[level] != NULL) {
      TimeNow = LevelArray[level]->GridData->ReturnTime();
      break;
    }

  /* Set the units. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
	   &VelocityUnits, TimeNow);

  /* Merge all yet-to-be born stars within r_merge */

  rmerge2o = powf(StarClusterCombineRadius * pc / LengthUnits, 2.0f);

  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar) {
    if (ThisStar->IsActive() || ThisStar->MarkedToDelete())
      continue;

    // Merge stars with separations less than a
    // StarClusterCombineRadius or a cell width
    dx = TopGridDx[0] * POW(RefineBy, -ThisStar->ReturnLevel());
    dx2 = dx*dx;
    rmerge2 = max(rmerge2o, dx2);

    for (OtherStar = ThisStar->NextStar; OtherStar;
	 OtherStar = OtherStar->NextStar) {
      if (ThisStar->ReturnID() == OtherStar->ReturnID()) {
	if (debug) {
	  printf("%"ISYM" -- merging duplicate particle??\n", ThisStar->ReturnID());
	  printf("ThisStar:\n");
	  ThisStar->PrintInfo();
	  printf("OtherStar:\n");
	  OtherStar->PrintInfo();
	}
	ENZO_FAIL("Merging Duplicate Particle!?\n");
      }
      if (ThisStar->Mergable(*OtherStar))
	if (ThisStar->Separation2(*OtherStar) <= rmerge2) {
	  ThisStar->Merge(OtherStar);
	  OtherStar->MarkForDeletion();
//	  printf("Merging stars %"ISYM" and %"ISYM"\n", ThisStar->ReturnID(),
//		 OtherStar->ReturnID());
	} // ENDIF radius2 < rmerge2
    } // ENDFOR OtherStar
  } // ENDFOR ThisStar

  /* Delete all marked star particles and their associated normal
     particles */
  
  ThisStar = AllStars;
  AllStars = NULL;
  LastStar = NULL;
  while (ThisStar) {
    MoveStar = PopStar(ThisStar);  // ThisStar becomes the next star in PopStar()
    if (MoveStar->MarkedToDelete()) {
      MoveStar->DeleteCopyInGrid();
      MoveStar->DisableParticle(LevelArray); // convert to a massless particle
      DeleteStar(MoveStar);
    } else {
      // Re-insert at the end of the list to keep the ordering the
      // same
      if (LastStar == NULL)
	InsertStarAfter(AllStars, MoveStar);
      else
	InsertStarAfter(LastStar, MoveStar);
      LastStar = MoveStar;
    }
  } // ENDWHILE

  return SUCCESS;

}

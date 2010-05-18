/***********************************************************************
/
/  MERGES MBH PARTICLES THAT ARE CLOSE ENOUGH
/
/  written by: Ji-hoon Kim (Heavily adopted from StarParticleMergeNew.C)
/  date:       July, 2009
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

void DeleteStar(Star * &Node);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int StarParticleMergeMBH(LevelHierarchyEntry *LevelArray[], Star *&AllStars)
{

  Star *ThisStar, *OtherStar, *PrevStar;
  LevelHierarchyEntry *Temp;
  float rmerge2;
  double vcirc2;
  FLOAT TimeNow;
  int dim, level;
  const float pc = 3.086e18;
  const double Grav = 6.673e-8, Msun = 1.989e33;

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

  /* Merge multiple MBH particles within r_merge */

  rmerge2 = powf(MBHCombineRadius * pc / LengthUnits, 2.0f);

  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar) {
    if (ThisStar->ReturnType() != MBH || ThisStar->IsActive() || ThisStar->MarkedToDelete())
      continue;
    for (OtherStar = ThisStar->NextStar; OtherStar;
	 OtherStar = OtherStar->NextStar) {
      if (ThisStar->ReturnID() == OtherStar->ReturnID()) {
	printf("%"ISYM" -- merging duplicate particle??\n", ThisStar->ReturnID());
	ENZO_FAIL("");
      }


      /* To merge two MBH particles you should satisfy 3 conditions 
         (1) the same types 
         (2) separation less than MBHCombineRadius 
         (3) relative velocity less than the circular velocity */

      vcirc2 = Grav * (ThisStar->ReturnMass() + OtherStar->ReturnMass()) * Msun / 
	ThisStar->Separation(OtherStar); // reduced mass dynamics

      if (ThisStar->MergableMBH(OtherStar) == 1 &&
	  ThisStar->Separation2(OtherStar) < rmerge2 &&
	  ThisStar->RelativeVelocity2(OtherStar) < vcirc2) {
	ThisStar->Merge(OtherStar);
	OtherStar->MarkForDeletion();
//	  printf("Merging stars %"ISYM" and %"ISYM"\n", ThisStar->ReturnID(),
//		 OtherStar->ReturnID());
      } 

    } // ENDFOR OtherStar
  } // ENDFOR ThisStar

  /* Delete all marked star particles and their associated normal
     particles */
  
  ThisStar = AllStars;
  while (ThisStar)
    if (ThisStar->MarkedToDelete()) {
      ThisStar->DeleteCopyInGrid();
      ThisStar->DisableParticle(LevelArray);
      DeleteStar(ThisStar); // ThisStar becomes the next star in DeleteStar()
    } else
      ThisStar = ThisStar->NextStar;

  return SUCCESS;

}

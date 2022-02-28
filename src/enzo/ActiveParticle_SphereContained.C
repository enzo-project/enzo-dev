/***********************************************************************
/
/  FIND ACCRETION SPHERE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  November, 2011 (JHW) -- porting to active particles
/
/  PURPOSE: When we remove baryons from the grid to add to the star
/           particle, look for a sphere that contains twice its mass.
/           Stepping outward by a cell width.
/
************************************************************************/
#ifdef USE_MPI
#endif /* USE_MPI */

#include "preincludes.h"
#include "performance.h"
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
#include "CommunicationUtilities.h"
#include "ActiveParticle.h"

int ActiveParticleType::SphereContained(LevelHierarchyEntry *LevelArray[], int level, 
					float Radius)
{

  LevelHierarchyEntry *Temp;
  int i, dim, direction, cornersContained, Rank, result;
  bool inside;
  Eint32 cornerDone[8];
  int Dims[MAX_DIMENSION];
  FLOAT corners[MAX_DIMENSION][8];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  /**************************************************************

     Compute corners of cube that contains a sphere of r=Radius

   **************************************************************/

  for (i = 0; i < 8; i++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {

      // If the bit is true, forward.  If not, reverse.
      direction = (i >> dim & 1) ? 1 : -1;
      corners[dim][i] = pos[dim] + direction * Radius;
    }
    cornerDone[i] = 0;
  }

  /* Check if the influenced sphere is contained within the grids on
     this level */

  Temp = LevelArray[level];
  for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
    // SG. This does almost nothing.
    // if (Temp->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      Temp->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);

      for (i = 0; i < 8; i++) {
	if (cornerDone[i]) continue;  // Skip if already locally found
	inside = true;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  inside &= (corners[dim][i] >= LeftEdge[dim] &&
		     corners[dim][i] <= RightEdge[dim]);
	if (inside)
	  cornerDone[i] = 1;
      } // ENDFOR corners

    //} // ENDIF MyProcessorNumber == ProcessorNumber
  } // ENDFOR grids

  cornersContained = 0;
  for (i = 0; i < 8; i++)
    cornersContained += cornerDone[i];
    //fprintf(stderr, "cornersContained = %d\n", __FUNCTION__, cornersContained);
  result = (cornersContained == 8);
  return result;

}

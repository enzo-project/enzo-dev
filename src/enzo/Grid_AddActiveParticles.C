/***********************************************************************
/
/  GRID CLASS (ADD ACTIVE PARTICLE DATA TO GRID)
/
/  written by: John Wise
/  date:       December, 2011
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#include "preincludes.h"

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
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"

int grid::AddActiveParticles(ActiveParticleList<ActiveParticleType> &NewParticles, int start, int end)
{

  if (NewParticles.size() == 0)
    return SUCCESS;

  this->NumberOfActiveParticles += (end - start);
  for (int i=start; i < end; i++)
  {
    NewParticles[i]->SetGridID(this->ID);
    NewParticles[i]->AssignCurrentGrid(this);
    this->ActiveParticles.copy_and_insert(*NewParticles[i]);
  }

#define NO_DEBUG_APs
#ifdef DEBUG_APs
  int dim, inside;
  FLOAT *pos;
  float TotalMass = 0;

  if (NumberOfActiveParticles != this->ActiveParticles.size()) {
    printf("Active particle count mismatch!\n");
    printf("NumberOfActiveParticles = %"ISYM"\n", NumberOfActiveParticles);
    printf("ActiveParticles.size() = %"ISYM"\n", ActiveParticles.size());
    ENZO_FAIL("")
  }

  for (int i = start; i < end; i++) {
    pos = NewParticles[i]->ReturnPosition();
    TotalMass += NewParticles[i]->ReturnMass();
    inside = this->PointInGrid(pos);
    if (inside == FALSE) {
      fprintf(stderr,"pos[0]: %"PSYM", pos[1]: %"PSYM", pos[2]: %"PSYM"\n",pos[0],pos[1],pos[2]);
      fprintf(stderr,"mass: %"FSYM"\n");
      fprintf(stderr,"GridLeftEdge[0]: %"PSYM", GridLeftEdge[1]: %"PSYM", GridLeftEdge[2]: %"PSYM"\n",
	      GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
      fprintf(stderr,"GridRightEdge[0]: %"PSYM", GridRightEdge[1]: %"PSYM", GridRightEdge[2]: %"PSYM"\n",
	      GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);      
      ENZO_FAIL("ActiveParticle outside!\n");
    }
  }
  fprintf(stdout,"AddActiveParticles: Total Mass added to grid = %"FSYM"\n",TotalMass);
#endif /* DEBUG */  

  return SUCCESS;

}

 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

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

//variables  for RandomForcing

int CallProblemSpecificRoutines(TopGridData * MetaData, HierarchyEntry *ThisGrid,
				int GridNum, float *norm, float TopGridTimeStep, 
				int level, int LevelCycleCount[])
{

  /* Add RandomForcing fields to velocities after the copying of current
     fields to old. I also update the total energy accordingly here.
     It makes no sense to force on the very first time step. */
 
  if (MetaData->CycleNumber > 0)
    ThisGrid->GridData->AddRandomForcing(norm, TopGridTimeStep);

  //dcc cut stop Forcing

  if (ProblemType == 24)
    ThisGrid->GridData->SphericalInfallGetProfile(level, 1);
  if (ProblemType == 30)
    ThisGrid->GridData->AnalyzeTrackPeaks(level, 0);
  if (ProblemType == 27)
    if (ThisGrid->GridData->ReturnProcessorNumber()==MyProcessorNumber){
      float AM[3], MeanVelocity[3], DMVelocity[3];
      FLOAT Center[] = {0,0,0}, CenterOfMass[3], DMCofM[3];
      ThisGrid->GridData->CalculateAngularMomentum
	(Center, AM, MeanVelocity, DMVelocity, CenterOfMass, DMCofM);
      fprintf(stdout, 
	      "level = %"ISYM" %"ISYM" %"ISYM"  "
	      "Vel %"FSYM" %"FSYM" %"FSYM"  "
	      "DMVel %"FSYM" %"FSYM" %"FSYM"  "
	      "CofM %"PSYM" %"PSYM" %"PSYM"  "
	      "DMCofM %"FSYM" %"FSYM" %"FSYM"\n",
	      level, LevelCycleCount[level], GridNum, MeanVelocity[0],
	      MeanVelocity[1], MeanVelocity[2],
	      DMVelocity[0], DMVelocity[1], DMVelocity[2],
	      -CenterOfMass[0], -CenterOfMass[1], -CenterOfMass[2],
	      DMCofM[0], DMCofM[1], DMCofM[2]);
    }

  /* Solve analytical free-fall */
  if (ProblemType == 63) {
    ThisGrid->GridData->SolveOneZoneFreefall();
  }

  return SUCCESS;
}

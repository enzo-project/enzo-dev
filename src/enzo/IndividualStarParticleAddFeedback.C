#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
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


int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);


int IndividualStarParticleAddFeedback(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
                                      int level, Star* &AllStars, bool* &AddedFeedback){




  Star *cstar;
  LevelHierarchyEntry *Temp;

  bool FeedbackOnGrid;

  FLOAT Time;
  float dxThisLevel;
  FLOAT *pos;
  float *vel;

  if (AllStars == NULL)
    return SUCCESS;

  LCAPERF_START("IndividualStarParticleAddFeedback");

  Temp        = LevelArray[level];
  Time        = Temp->GridData->ReturnTime();

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);

  /* make sure feedback will happen somewhere on this grid */
  int Rank, Dims[MAX_DIMENSION];
  float CellWidth;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  Temp->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
  CellWidth = (RightEdge[0] - LeftEdge[0]) / (Dims[0] - 2*NumberOfGhostZones);

  int count = 0;

  for (cstar = AllStars; cstar; cstar = cstar->NextStar, count++){

    AddedFeedback[count] = false;

    printf("cstar type = %"ISYM"\n", cstar->ReturnType());
    printf("cstar flag = %"ISYM"\n", cstar->ReturnFeedbackFlag());


    if( ABS(cstar->ReturnType()) != IndividualStar &&
        ABS(cstar->ReturnType()) != IndividualStarWD &&
        ABS(cstar->ReturnType()) != IndividualStarRemnant){
      continue; // what are we doing here in the first place?
    }

    /* Check feedback flag and see if particle is doing anything intersting */
    if( cstar->ReturnFeedbackFlag() < INDIVIDUAL_STAR_STELLAR_WIND ||
        cstar->ReturnFeedbackFlag() > INDIVIDUAL_STAR_WIND_AND_SN){
      continue; // skip to next star
    }

    //
    // for winds and supernova, feedback is applied on a n x n cell stencil
    // check and see if any of the affected cells are on this grid level
    //
    printf("cstar level and grid level = %"ISYM" %"ISYM"\n",cstar->ReturnLevel(),level);
    if( cstar->ReturnLevel() != level){
      continue;
    }

    /* make sure star's feedback will land somewhere on this grid - skip if not*/
    int ncell = (int) ((IndividualStarFeedbackStencilSize+1)/2.0 - 1);

    pos = cstar->ReturnPosition();
    vel = cstar->ReturnVelocity();

    printf("position %"FSYM" %"FSYM " %"FSYM"\n",pos[0], pos[1], pos[2]);
    printf("ncell cell width %"ISYM"\n",ncell);

    if( (pos[0] - (ncell + 0.5)*CellWidth > RightEdge[0]) ||
        (pos[0] + (ncell + 0.5)*CellWidth < LeftEdge[0])  ||
        (pos[1] - (ncell + 0.5)*CellWidth > RightEdge[1]) ||
        (pos[1] + (ncell + 0.5)*CellWidth < LeftEdge[1])  ||
        (pos[2] - (ncell + 0.5)*CellWidth > RightEdge[2]) ||
        (pos[2] + (ncell + 0.5)*CellWidth < LeftEdge[2])){
      // particle feedback zone is not on grid at all. skip
      continue;
    }
    printf("passing position check\n");

    /* Add stellar wind feedback if present */
    double particle_mass = cstar->ReturnMass();
    int feedback_mode = -9999;
    // AJE TO DO NEED TO CHECK MASS UNITS HERE
    if( cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_STELLAR_WIND ||
        cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_WIND_AND_SN){
      feedback_mode = -1;
    }
    if( cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_SNII ||
        cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_WIND_AND_SN){
      feedback_mode = 1;
    }
    if( cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_SNIA){
      feedback_mode = 2;
      if(Temp->GridData->isLocal()){
        Temp->GridData->IndividualStarAddFeedbackGeneral(pos[0], pos[1], pos[2],
                                                         vel[0], vel[1], vel[2],
                                                         cstar->ReturnBirthMass(), cstar->ReturnLifetime(),
                                                         Time - cstar->ReturnBirthTime(),
                                                         cstar->ReturnMetallicity(), &particle_mass, 2);
      }
      AddedFeedback[count] = true;
    }

    if (feedback_mode > -9999){

      for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){

        for (Temp = LevelArray[l]; Temp; Temp = Temp ->NextGridThisLevel){

          if(Temp->GridData->isLocal()){
            Temp->GridData->IndividualStarAddFeedbackGeneral(pos[0], pos[1], pos[2],
                                                           vel[0], vel[1], vel[2],
                                                           cstar->ReturnBirthMass(), cstar->ReturnLifetime(),
                                                           Time - cstar->ReturnBirthTime(),
                                                           cstar->ReturnMetallicity(), &particle_mass, 2);
            cstar->SetNewMass(particle_mass); // CHECK MASS UNITS
            AddedFeedback[count] = true;
          }
        }
      }
    }

  } // end stars loop

}



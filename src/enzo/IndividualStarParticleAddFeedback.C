/******************************************************************
/
/ HANDLES FEEDBACK FOR INDIVIDUAL STAR PARTICLES
/
/ written by: Andrew Emerick
/ date:       August 2016
/ modified:
/
/ PURPOSE: Routine applies feedback from individual star particles
/          (stellar winds or supernovae) using the CIC stencil
/          interpolation methods from Simpson et. al. 2015. This
/          routine exists following the pre-existing Star Particle
/          class formalism in order to allow consistent feedback
/          across grids when particles are near grid boundaries.
/          This is a substantial improvement over previous method
/          which "shifted" feedback zone to avoid this.
/
/
/ OUTSTANDIG ISSUES: Type Ia supernovae rely on random number generator
/                    to go off... need to do this consistently (somehow)
/                    across all processors. Do in Set feedback flag
/
/
**********************************************************************/

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



int IndividualStarParticleAddFeedback(TopGridData *MetaData,
                                      LevelHierarchyEntry *LevelArray[],
                                      int level, Star* &AllStars,
                                      bool* &AddedFeedback){

  Star *cstar;
  LevelHierarchyEntry *Temp;

  FLOAT Time;
  FLOAT *pos;
  float dxThisLevel;
  float *vel;

  if (AllStars == NULL)
    return SUCCESS;

  LCAPERF_START("IndividualStarParticleAddFeedback");

  Temp        = LevelArray[level];
  Time        = Temp->GridData->ReturnTime();


  /* figure out some grid properties */
  int Rank, Dims[MAX_DIMENSION];
  float CellWidth;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  Temp->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
  CellWidth = (RightEdge[0] - LeftEdge[0]) / (Dims[0] - 2*NumberOfGhostZones);

  int count = 0;

  /* Loop over all stars, checking properties before feedback */
  for (cstar = AllStars; cstar; cstar = cstar->NextStar, count++){

    AddedFeedback[count] = false;

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

    if( cstar->ReturnLevel() != level){
      continue; // only apply feedback on level of star
    }

    /* make sure star's feedback will land somewhere on this grid - skip if not*/
    int ncell = (int) ((IndividualStarFeedbackStencilSize+1)/2.0 - 1);

    pos = cstar->ReturnPosition();
    vel = cstar->ReturnVelocity();

    if( (pos[0] - (ncell + 0.5)*CellWidth > RightEdge[0]) ||
        (pos[0] + (ncell + 0.5)*CellWidth < LeftEdge[0])  ||
        (pos[1] - (ncell + 0.5)*CellWidth > RightEdge[1]) ||
        (pos[1] + (ncell + 0.5)*CellWidth < LeftEdge[1])  ||
        (pos[2] - (ncell + 0.5)*CellWidth > RightEdge[2]) ||
        (pos[2] + (ncell + 0.5)*CellWidth < LeftEdge[2])){
      // particle feedback zone is not on grid at all. skip
      // this check is performed also in actual feedback routines, but
      // redundancy is O.K. here
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

  LCAPERF_END("IndividualStarParticleAddFeedback");
  return SUCCESS;
}



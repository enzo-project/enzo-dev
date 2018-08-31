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
#include "EnzoTiming.h"
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


int IsParticleFeedbackInGrid(float *pos, int ncell, LevelHierarchyEntry *Temp);

int IndividualStarParticleAddFeedback(TopGridData *MetaData,
                                      LevelHierarchyEntry *LevelArray[],
                                      int level, Star* &AllStars,
                                      bool* &AddedFeedback){

  /* stars and hierarchy */
  Star *cstar;
  LevelHierarchyEntry *Temp;

  /* pos and vel for star */
  FLOAT *pos;
  float *vel;

  if (AllStars == NULL)
    return SUCCESS;

  TIMER_START("IndividualStarParticleAddFeedback");

  int count = 0;

  /* Loop over all stars, checking properties before doing feedback */
  for (cstar = AllStars; cstar; cstar = cstar->NextStar, count++){

    AddedFeedback[count] = false;

    if( ABS(cstar->ReturnType()) != IndividualStar &&
        ABS(cstar->ReturnType()) != IndividualStarWD &&
        ABS(cstar->ReturnType()) != IndividualStarRemnant &&
        ABS(cstar->ReturnType()) != IndividualStarUnresolved ){
      continue; // This probably should never need to be checked
    }

    /* Check feedback flag - skip if particle isn't doing anything interesting */
    if( cstar->ReturnFeedbackFlag() < INDIVIDUAL_STAR_STELLAR_WIND ||
        cstar->ReturnFeedbackFlag() > INDIVIDUAL_STAR_WIND_AND_SN){
      continue; // skip to next star
    }

    if( cstar->ReturnLevel() > level){
      continue; // only apply feedback on level of star
    }

    if(cstar->ReturnMass() < 0.0){
        cstar->PrintInfo();
        ENZO_FAIL("Particle Mass initially negative in IndividualStarParticleAddFeedback");
    }


    /* feedback is done in a cic interpolation. This is number of cells
       on eiher side of central cell (i.e. 3x3 CIC -> ncell = 1) */
//    int ncell = (int) ((IndividualStarFeedbackStencilSize+1)/2.0 - 1);

    int ncell = (int) (IndividualStarFeedbackStencilSize) + 1;

    pos = cstar->ReturnPosition();
    vel = cstar->ReturnVelocity();

    double particle_mass;

    //
    // Check Stellar Winds - Apply if particle on grid and grid local
    //
    if( cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_STELLAR_WIND ||
        cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_WIND_AND_SN){
      for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){

        for (Temp = LevelArray[l]; Temp; Temp = Temp ->NextGridThisLevel){

          if(Temp->GridData->isLocal() && IsParticleFeedbackInGrid(pos, ncell, Temp) ){
            // refresh mass every time to prevent double+ counting loss
            particle_mass = cstar->ReturnMass();

            Temp->GridData->IndividualStarAddFeedbackSphere(cstar, pos[0], pos[1], pos[2],
                                                            vel[0], vel[1], vel[2],
                                                            cstar->ReturnBirthMass(), cstar->ReturnLifetime(),
                                                            Temp->GridData->ReturnTime() - cstar->ReturnBirthTime(),
                                                            cstar->ReturnMetallicity(), &particle_mass, -1);

/*
            Temp->GridData->IndividualStarAddFeedbackGeneral(pos[0], pos[1], pos[2],
                                                           vel[0], vel[1], vel[2],
                                                           cstar->ReturnBirthMass(), cstar->ReturnLifetime(),
                                                           Temp->GridData->ReturnTime() - cstar->ReturnBirthTime(),
                                                           cstar->ReturnMetallicity(), &particle_mass, -1);
                                                           // < 0 in last arg signifies stellar winds
*/
            AddedFeedback[count] = TRUE;
          }
        }
      }

      if (AddedFeedback[count]){ // only if this particle did something
//        cstar->PrintInfo();
        float old_mass = cstar->ReturnMass();
        cstar->SetNewMass(particle_mass); // update mass (only once)
        cstar->AddToWindMassEjected(old_mass - particle_mass);
//        cstar->PrintInfo();
      }
    }

    //
    // Check Core Collapse Supernova
    if( cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_SNII ||
        cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_WIND_AND_SN){
      for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){

        for (Temp = LevelArray[l]; Temp; Temp = Temp ->NextGridThisLevel){

          if(Temp->GridData->isLocal() && IsParticleFeedbackInGrid(pos, ncell, Temp)){
            // refresh mass every time to prevent double+ counting
            particle_mass = cstar->ReturnMass();
            Temp->GridData->IndividualStarAddFeedbackSphere(cstar, pos[0], pos[1], pos[2],
                                                            vel[0], vel[1], vel[2],
                                                            cstar->ReturnBirthMass(), cstar->ReturnLifetime(),
                                                            Temp->GridData->ReturnTime() - cstar->ReturnBirthTime(),
                                                            cstar->ReturnMetallicity(), &particle_mass, 1);

/*
            Temp->GridData->IndividualStarAddFeedbackGeneral(pos[0], pos[1], pos[2],
                                                           vel[0], vel[1], vel[2],
                                                           cstar->ReturnBirthMass(), cstar->ReturnLifetime(),
                                                           Temp->GridData->ReturnTime() - cstar->ReturnBirthTime(),
                                                           cstar->ReturnMetallicity(), &particle_mass, 1);
                                                           // 1 in last arg signifies Core collapse SN
*/
            AddedFeedback[count] = TRUE;
          }
        }
      }

      if (AddedFeedback[count]){
        float old_mass = cstar->ReturnMass();
        AddedFeedback[count] = true;
        cstar->SetFeedbackFlag(INDIVIDUAL_STAR_SN_COMPLETE);
        cstar->SetNewMass(particle_mass); // update mass (only once)
        cstar->AddToSNMassEjected(old_mass - particle_mass);
      }
    }

    //
    // Check Type Ia Supernova
    //
    if( cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_SNIA){
      for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){
        for (Temp = LevelArray[l]; Temp; Temp = Temp ->NextGridThisLevel){

          if(Temp->GridData->isLocal() && IsParticleFeedbackInGrid(pos, ncell, Temp)){
            particle_mass = cstar->ReturnMass();

            Temp->GridData->IndividualStarAddFeedbackSphere(cstar, pos[0], pos[1], pos[2],
                                                            vel[0], vel[1], vel[2],
                                                            cstar->ReturnBirthMass(), cstar->ReturnLifetime(),
                                                            Temp->GridData->ReturnTime() - cstar->ReturnBirthTime(),
                                                            cstar->ReturnMetallicity(), &particle_mass, 2);


/*
            Temp->GridData->IndividualStarAddFeedbackGeneral(pos[0], pos[1], pos[2],
                                                           vel[0], vel[1], vel[2],
                                                           cstar->ReturnBirthMass(), cstar->ReturnLifetime(),
                                                           Temp->GridData->ReturnTime() - cstar->ReturnBirthTime(),
                                                           cstar->ReturnMetallicity(), &particle_mass, 2);
                                                           // 2 in last arg signifies Type 1a
*/
          }
        }
      }

      AddedFeedback[count] = true;
      cstar->SetFeedbackFlag(INDIVIDUAL_STAR_SN_COMPLETE);
      cstar->AddToSNMassEjected(cstar->ReturnMass()); // not the actual mass ejection from SNIA !!!
      cstar->SetNewMass(0.0); // now a massless tracer
    }

    if(cstar->ReturnMass() < 0.0){
        cstar->PrintInfo();
        ENZO_FAIL("Particle Mass going negative in IndividualStarParticleAddFeedback");
    }

  } // end stars loop


  // debugging loop to ensure validity of mass ejection
  if (TRUE) {
    for (cstar = AllStars; cstar; cstar = cstar->NextStar){
      cstar->CheckMassEjectionValidity();
    }
  }


  TIMER_STOP("IndividualStarParticleAddFeedback");
  return SUCCESS;
}


int IsParticleFeedbackInGrid(float *pos, int ncell, LevelHierarchyEntry *Temp){
/* Check and see if particle feedback zone overlaps with any portion
   of grid before calling feedback functions. This is checked as well
   in feedback functions but reduces overhead somewhat by having
   a redundant check here
*/

  int Rank, Dims[MAX_DIMENSION];
  float CellWidth;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  Temp->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
  CellWidth = (RightEdge[0] - LeftEdge[0]) / Dims[0];

  float fudge = 0.1;

  if( (pos[0] - (ncell + fudge)*CellWidth > RightEdge[0]) ||
      (pos[0] + (ncell + fudge)*CellWidth < LeftEdge[0])  ||
      (pos[1] - (ncell + fudge)*CellWidth > RightEdge[1]) ||
      (pos[1] + (ncell + fudge)*CellWidth < LeftEdge[1])  ||
      (pos[2] - (ncell + fudge)*CellWidth > RightEdge[2]) ||
      (pos[2] + (ncell + fudge)*CellWidth < LeftEdge[2])){
      // particle feedback zone is not on grid at all. skip
      // this check is performed also in actual feedback routines, but
      // redundancy is O.K. here
      return FALSE;
  }

  return TRUE;
}



/***********************************************************************
/
/  FROM THE ATTRIBUTES, DETERMINE WHETHER THE STAR DIES
/
/  written by: John Wise
/  date:       April, 2009
/  modified1:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
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

#include "StarParticleData.h" // AJE

#define NO_DEATH 0
#define KILL_STAR 1
#define KILL_ALL 2

int Star::HitEndpoint(FLOAT Time)
{

  const float TypeIILowerMass = 11, TypeIIUpperMass = 40.1;
  const float PISNLowerMass = 140, PISNUpperMass = 260;

  /* First check if the star's past its lifetime and then check other
     constrains based on its star type */

  int result = NO_DEATH;
  if ((Time > this->BirthTime + this->LifeTime) && this->type >=0)
    result = KILL_STAR;
  else
    return result;

  switch (this->type) {

  case PopIII:
    // If a Pop III star is going supernova, only kill it after it has
    // applied its feedback sphere
    if ((this->Mass >= PISNLowerMass && this->Mass <= PISNUpperMass) ||
        ((this->Mass >= TypeIILowerMass && this->Mass <= TypeIIUpperMass) &&
            PopIIISupernovaExplosions == TRUE)) {

      // Needs to be non-zero (multiply by a small number to retain
      // memory of mass)
      if (this->FeedbackFlag == DEATH) {
	this->Mass *= tiny_number;  


	// Set lifetime so the time of death is exactly now.
	this->LifeTime = Time - this->BirthTime;

	//this->FeedbackFlag = NO_FEEDBACK;
	result = KILL_STAR;
	//result = NO_DEATH;
      } else {
	result = NO_DEATH;
      }

    // Check mass: Don't want to kill tracer SN particles formed
    // (above) in the previous timesteps.

    } else if (this->Mass > 1e-9) {
      // Turn particle into a black hole (either radiative or tracer)
      if (PopIIIBlackHoles) {
	this->type = BlackHole;
	this->LifeTime = huge_number;
	this->FeedbackFlag = NO_FEEDBACK;
	result = NO_DEATH;
      } else {
	this->type = PARTICLE_TYPE_DARK_MATTER;
	result = KILL_STAR;
      }
    } else // SN tracers (must refine)
      result = NO_DEATH;

    if (debug)
      printf("HitEndpoint[%"ISYM"]: type = %"ISYM", mass = %"GOUTSYM", result = %"ISYM", feedback = %"ISYM", Time = %"PSYM"/%"PSYM"\n",
	     this->Identifier, this->type, this->Mass, result, this->FeedbackFlag, Time,
	     this->BirthTime+this->LifeTime);

    break;
    
  case PopII:
    break;

  case BlackHole:
    break;
    
  case MBH:
    break;

  case PopIII_CF:
    break;

  case NormalStar:
    break;

  case IndividualStar:
    {
      result = NO_DEATH; // never kill these stars

      float mproj;
      mproj       = this->BirthMass;

      /* check mass */
      if(mproj > IndividualStarSNIIMassCutoff && mproj < IndividualStarDirectCollapseThreshold){

          if(this->FeedbackFlag == INDIVIDUAL_STAR_SN_COMPLETE){
             //(this->FeedbackFlag == INDIVIDUAL_STAR_SNII ||
             //this->FeedbackFlag == INDIVIDUAL_STAR_WIND_AND_SN)){
              this->LifeTime = huge_number*this->LifeTime; // make to a Hubble time
              this->type = IndividualStarRemnant;
              this->FeedbackFlag = NO_FEEDBACK;
          } else { // havent blown up yet
              this->FeedbackFlag = INDIVIDUAL_STAR_SNII;
          }

      } else if (mproj > IndividualStarDirectCollapseThreshold){
        this->FeedbackFlag = NO_FEEDBACK;
        this->LifeTime     = huge_number * this->LifeTime;
        this->type         = IndividualStarRemnant;
      }

      if(mproj >= IndividualStarWDMinimumMass && mproj <= IndividualStarWDMaximumMass){

          if(this->FeedbackFlag == INDIVIDUAL_STAR_STELLAR_WIND){

              this->type = -IndividualStarWD;           // set negative
//              this->LifeTime = huge_number ; // very large number

              float wd_mass;
              if(   mproj < 4.0){ wd_mass = 0.134 * mproj + 0.331;}
              else if (mproj >=4.0){ wd_mass = 0.047 * mproj + 0.679;}
              this->Mass = wd_mass;

              this->FeedbackFlag = NO_FEEDBACK;
          } else{

            // otherwise, we haven't done the last phase of AGB wind yet
            // keep this particle around for another timestep
            this->FeedbackFlag = INDIVIDUAL_STAR_STELLAR_WIND;

          }

      } // end WD check


      break;
    }
  case IndividualStarWD:
    result = NO_DEATH;

    if ( this->FeedbackFlag == INDIVIDUAL_STAR_SN_COMPLETE){
         //this->FeedbackFlag == INDIVIDUAL_STAR_SNIA ){
      this->Mass         = 0.0;
      this->LifeTime     = huge_number * this->LifeTime;
      this->FeedbackFlag = NO_FEEDBACK;
    } else{
      this->FeedbackFlag = INDIVIDUAL_STAR_SNIA;
    }

    break;

  case IndividualStarRemnant:
  case IndividualStarUnresolved:
    result = NO_DEATH; //printf("Individual Star remnamt in hit endpoint\n");
    this->FeedbackFlag = NO_FEEDBACK;
    this->LifeTime     = huge_number * this->LifeTime;
    break;

  } // ENDSWITCH

  return result;
}

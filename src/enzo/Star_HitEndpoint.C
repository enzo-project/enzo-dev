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
#include "phys_constants.h"

#define NO_DEATH 0
#define KILL_STAR 1
#define KILL_ALL 2

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int Star::HitEndpoint(FLOAT Time)
{

//const float TypeIILowerMass = 11, TypeIIUpperMass = 40.1;
  const float TypeIILowerMass  = 8, TypeIIUpperMass  = 40.1;
  const float FaintSNLowerMass = 8, FaintSNUpperMass = 140;
  const float PISNLowerMass = 140, PISNUpperMass = 260;
  const float tlife[] = 
           {  2.54502e+07,  9.27383e+06,  6.69605e+06,  5.58739e+06
           ,  4.21361e+06,  2.43907e+06,  2.28077e+06,  2.12901e+06  };
             // lifetime of PopII stars [year]
             // approximate Pop III stellar lifetime (Schaerer 2002)
  float tdyn, old_mass, frac;

  /* Set the units. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  /* First check if the star's past its lifetime and then check other
     constrains based on its star type */

  int result = NO_DEATH;
  if ((Time > this->BirthTime + this->LifeTime) && this->type >=0)
    result = KILL_STAR;
  else if (this->type != PopII)
    return result;

  switch (this->type) {

  case PopIII:
    // If a Pop III star is going supernova, only kill it after it has
    // applied its feedback sphere
    if(((MetalChemistry == 0 &&
         ((this->Mass >= PISNLowerMass && this->Mass <= PISNUpperMass) ||
          (this->Mass >= TypeIILowerMass && this->Mass <= TypeIIUpperMass))) ||
        (MetalChemistry >  0 &&
         (((this->Metallicity <  PopIIIMetalCriticalFraction) && (this->Mass >= FaintSNLowerMass && this->Mass <= FaintSNUpperMass)) ||
          ((this->Metallicity >= PopIIIMetalCriticalFraction) && (this->Mass >= TypeIILowerMass  && this->Mass <= TypeIIUpperMass ))))) &&
            PopIIISupernovaExplosions == TRUE)
         {

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
    if (MetalChemistry) {
      old_mass = (float)(this->Mass);

      tdyn = (Time - this->BirthTime) * (TimeUnits/yr_s); // year
      if( tdyn > tlife[7] && this->Mass > 260.0 ) {
        this->Mass -= (this->Mass > 301.0 ? 300.0 : this->Mass - 1.0);
      } else if( tdyn > tlife[6] && this->Mass > 185.0 ) {
        this->Mass -= (this->Mass > 201.0 ? 200.0 : this->Mass - 1.0);
      } else if( tdyn > tlife[5] && this->Mass > 140.0 ) {
        this->Mass -= (this->Mass > 171.0 ? 170.0 : this->Mass - 1.0);
      } else if( tdyn > tlife[4] && this->Mass > 40.0 ) {
        this->Mass -= (this->Mass >  91.0 ?  90.0 : this->Mass - 1.0);
      } else if( tdyn > tlife[3] && this->Mass > 27.5 ) {
        this->Mass -= (this->Mass >  31.0 ?  30.0 : this->Mass - 1.0);
      } else if( tdyn > tlife[2] && this->Mass > 22.5 ) {
        this->Mass -= (this->Mass >  26.0 ?  25.0 : this->Mass - 1.0);
      } else if( tdyn > tlife[1] && this->Mass > 16.5 ) {
        this->Mass -= (this->Mass >  21.0 ?  20.0 : this->Mass - 1.0);
      } else if( tdyn > tlife[0] && this->Mass > 8.0 ) {
        this->Mass -= (this->Mass >  14.0 ?  13.0 : this->Mass - 1.0);
      }
    }
#ifdef UNUSED
    /* Momentum conservation */
    frac = old_mass / this->Mass;
    this->vel[0] *= frac;
    this->vel[1] *= frac;
    this->vel[2] *= frac;
#endif
    break;

  case BlackHole:
    break;
    
  case MBH:
    break;

  case PopIII_CF:
    break;

  } // ENDSWITCH

  return result;
}


/***********************************************************************
/
/  SUBTRACT MASS FROM STAR AFTER FEEDBACK
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/  modified1: 
/
/  PURPOSE: Currently this file is not being used; now the functionality 
/           of this file has moved to
/
/           MBH_THERMAL -> no need to do this because EjectaDensity = 0
/           MBH_JETS    -> now done in Grid_AddFeedbackSphere.C
/
************************************************************************/
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
#include "CosmologyParameters.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

int FindField(int field, int farray[], int numfields);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int Star::RemoveMassFromStarAfterFeedback(float &Radius, double &EjectaDensity, 
					  float DensityUnits, float LengthUnits,
					  int &CellsModified)
{

  double old_mass;
  double Msun = 1.989e33;

  /* Check if the star type is correct */

  if ((ABS(this->type) != MBH) || 
      (this->FeedbackFlag != MBH_THERMAL && this->FeedbackFlag != MBH_JETS))
    return SUCCESS;

  /* Check if there really was a feedback; this is crucial for MBH_JETS
     because if the jet mass doesn't exceed the threshold it won't do any feedback */

  if (EjectaDensity != 0 ||
      CellsModified == 0 ||
      this->CurrentGrid == NULL)
    return SUCCESS;

  float BubbleVolume = (4.0 * M_PI / 3.0) * Radius * Radius * Radius;

  old_mass = this->Mass;
//  printf("star::RMFSAF: before: old_mass = %lf\n", old_mass);  
//  printf("star::RMFSAF: vel = %"FSYM" %"FSYM" %"FSYM"\n", 
//	 this->vel[0], this->vel[1], this->vel[2]);

  /* Now let's start working!  Because we injected the mass in Grid_AddFeedbackSphere, 
     we should subtract the mass from the star */

  switch (this->FeedbackFlag) {

  case MBH_THERMAL:  

    //this actually would not do anything because EjectaDensity = 0 for MBH_THERMAL
    //unless one changes the current scheme - Ji-hoon Kim, Jan.2010

    this->Mass -= EjectaDensity * DensityUnits * BubbleVolume * pow(LengthUnits,3.0) / Msun;  
    break;

  case MBH_JETS:

    //remove the mass; and now that we injected the NotEjectedMass, set it to zero 

    this->Mass -= this->NotEjectedMass;  
    this->NotEjectedMass = 0.0;
    break;

  } //ENDSWITCH

  /* Because the mass is subtracted out with zero net momentum, 
     increase the velocity accordingly */

  this->vel[0] *= old_mass / this->Mass; 
  this->vel[1] *= old_mass / this->Mass;
  this->vel[2] *= old_mass / this->Mass; 

//  printf("star::RMFSAF: after : this->Mass = %lf\n", this->Mass);  
//  printf("star::RMFSAF: vel = %"FSYM" %"FSYM" %"FSYM"\n", 
//	 this->vel[0], this->vel[1], this->vel[2]);

  return SUCCESS;

}

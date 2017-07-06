/***********************************************************************
/
/  CHECK FOR TIME ACTION
/
/  written by: Greg Bryan
/  date:       September, 2000
/  modified1:
/
/  PURPOSE:
/    This routine checks to see if the time has arrived for a given
/      "time-action" to occur (i.e. some action that is applied to
/       the data at a given time/redshift).
/
************************************************************************/

#include <stdlib.h> 
#include <math.h>
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
#include "CosmologyParameters.h"


/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int CommunicationBroadcastValues(float *Values, int Number, int BroadcastProcessor);

int IsParticleFeedbackInGrid(float *pos, int ncell, LevelHierarchyEntry *Temp);

int search_lower_bound(float *arr, float value, int low, int high, int total);

int CheckForTimeAction(LevelHierarchyEntry *LevelArray[],
		       TopGridData &MetaData)
{
 
  /* Declarations. */
 
  int i, level;
 
  /* Check for time against list of actions. */
 
  for (i = 0; i < MAX_TIME_ACTIONS; i++){
    if (MetaData.Time >= TimeActionTime[i] && TimeActionTime[i] > 0) {
 
      if (debug)
	printf("Applying TimeAction %"ISYM" at t=%"GOUTSYM"\n", i, TimeActionTime[i]);


      /* Following Miao's time action SN injection */
      if (TimeActionType[i] == 3){ // individual star SN

        float DensityUnits, LengthUnits,TemperatureUnits, TimeUnits,VelocityUnits;
        double MassUnits;
        if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                  &TimeUnits, &VelocityUnits, MetaData.Time) == FAIL) {
          ENZO_FAIL("Error in GetUnits in CheckForTimeAction.C.");
        }
        MassUnits = DensityUnits * LengthUnits * LengthUnits * LengthUnits;

        float Myr = 3.1556E13;
        float Mpc = 3.086E24;
        float pi  = 3.14156926;

        // in Myr
        float sntime = IndividualStarICSupernovaTime*Myr/TimeUnits;


        if (IndividualStarICSupernovaFromFile){ //  interpolate along tabulated rates
           // interpolate along icsupernovarate - icsupernovatime arrays
           if(MetaData.Time > sntime){
             TimeActionParameter[i] = -1;
             TimeActionTime[i]      = -1;
           } else {
             float time_now = MetaData.Time * TimeUnits / Myr;

             int index = search_lower_bound(ICSupernovaTimeArray, time_now, 0,
                                            ICSupernovaNumberOfPoints, ICSupernovaNumberOfPoints);
             float slope = (ICSupernovaSNRArray[index+1] - ICSupernovaSNRArray[index])/
                             (ICSupernovaTimeArray[index+1] - ICSupernovaTimeArray[index]);

             TimeActionParameter[i] = slope*(time_now - ICSupernovaTimeArray[index])
                                         + ICSupernovaSNRArray[index];
             TimeActionParameter[i] = 1.0 / TimeActionParameter[i];

             TimeActionTime[i] += TimeActionParameter[i] * 3.1556E7 / TimeUnits;

           }


        } else if (IndividualStarICSupernovaTime <= 0){ // go until stars form

          if (MetaData.NumberOfParticles > 0){ // turn off SN
            TimeActionParameter[i] = -1;
            TimeActionTime[i]      = -1;

          } else {
            TimeActionParameter[i] = 1.0 / IndividualStarICSupernovaRate;
            TimeActionTime[i]     += TimeActionParameter[i] * 3.1556E7 / TimeUnits;
          }

        } else{

          // otherwise, time sets a cutoff to turn off SN
          if (MetaData.Time > 0.5*sntime){
            TimeActionParameter[i] = IndividualStarICSupernovaRate * (sntime - MetaData.Time)/(0.5*sntime);
          } else{
            TimeActionParameter[i] = IndividualStarICSupernovaRate;
          }

          TimeActionParameter[i] = 1.0 / TimeActionParameter[i]; // convert to time from rate

          /* set next explosion time if there is one */
          if (MetaData.Time >= sntime){
              TimeActionTime[i] = -1;
          } else{
              TimeActionTime[i] += TimeActionParameter[i] * 3.1556E7 / TimeUnits;
              TimeActionTime[i] = min(TimeActionTime[i], sntime);
          }

        } // end supernova time check


        // compute the SN location
        // for simplicity, this will be randomly chosen in a uniform disk
        float SNPosition[3];
        for(int i = 0; i < 3; i++)
          SNPosition[i] = IndividualStarICSupernovaPos[i];

        float rmax = IndividualStarICSupernovaR * Mpc / LengthUnits;
        float zmax = IndividualStarICSupernovaZ * Mpc / LengthUnits;

        if (MyProcessorNumber == ROOT_PROCESSOR){
          float theta, z;

          theta = (rand()*1.0/RAND_MAX) * 2.0 * pi;

          if (IndividualStarICSupernovaMethod == 1){
            // choose points uniformly in disk / cylinder
            float r;

            r     = sqrt( rand()*1.0/RAND_MAX ) * rmax;
            z     = (rand()*1.0/RAND_MAX) * 2.0 *zmax - zmax;

            SNPosition[0] += r * cos(theta);
            SNPosition[1] += r * sin(theta);
            SNPosition[2] += z;
          } else if (IndividualStarICSupernovaMethod == 2){
            // choose points uniformly on sphere
            float x, y, u;

            u = ((rand()*1.0/RAND_MAX)*2.0 - 1.0);
            x = sqrt(1.0 - u*u) * cos(theta) * rmax;
            u = ((rand()*1.0/RAND_MAX)*2.0 - 1.0);
            y = sqrt(1.0 - u*u) * sin(theta) * rmax;

            z = ((rand()*1.0/RAND_MAX)*2.0 - 1.0) * rmax;

            SNPosition[0] += x;
            SNPosition[1] += y;
            SNPosition[2] += z;
          } // end SN method check for position

        } // end ROOT check


        CommunicationBroadcastValues(SNPosition, 3, ROOT_PROCESSOR);

        if (MyProcessorNumber == ROOT_PROCESSOR)
          printf("SNPosition = %f     %f     %f\n",SNPosition[0],SNPosition[1],SNPosition[2]);


        float m_eject   = 10.0 * 1.989E33 / MassUnits; // need to divide by dx before sending
        float E_thermal = IndividualStarSupernovaEnergy * 1.0E51 / (MassUnits*VelocityUnits*VelocityUnits); 

        for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++){
          LevelHierarchyEntry *Temp = LevelArray[level];
          while (Temp != NULL){

            // if explosion is on this grid (and grid is local), add feedback
            float dx = float(Temp->GridData->ReturnCellWidth());
            float m_eject_dens = m_eject / (dx*dx*dx);
            float e_eject      = E_thermal / (dx*dx*dx);


            if (IndividualStarICSupernovaInjectionMethod == 1){ // spherical bubble injection
              if (Temp->GridData->isLocal() && IsParticleFeedbackInGrid(SNPosition,
                                                                        IndividualStarFeedbackStencilSize,
                                                                        Temp)){

                  if ( Temp->GridData->IndividualStarInjectSphericalFeedback(NULL,
                                             SNPosition[0], SNPosition[1], SNPosition[2],
                                             m_eject_dens, e_eject, NULL, 0) == FAIL){
                    ENZO_FAIL("Error in grid->IndividualStarInjectSphericalFeedback called from CheckForTimeAction\n");
                  }

              } // end if on grid

            } else if (IndividualStarICSupernovaInjectionMethod == 2){ // Simpson et. al. thermal + kinetic

              if(Temp->GridData->isLocal() && IsParticleFeedbackInGrid(SNPosition,
                                                                           IndividualStarFeedbackStencilSize*2,
                                                                           Temp)){
                float f_k = 0.25; // kinetic fraction
                float e_eject_kin   = f_k * e_eject;
                float e_eject_therm = (1.0 - f_k) * e_eject;
                float vx = 0.0, vy = 0.0, vz = 0.0;   // really should set to local gas velocity
                float p_feedback = 0.0; // no explicit momentum

                float metal_mass[3];
                metal_mass[0] = 0.0186; // total metal
                metal_mass[1] = 0.72  ; // H
                metal_mass[2] = (1.0 - metal_mass[0] - metal_mass[1]); // He
                for (int i = 0; i < 3; i++) metal_mass[i] *= m_eject_dens;

                if ( Temp->GridData->IndividualStarInjectFeedbackToGrid(SNPosition[0],SNPosition[1],SNPosition[2],
                                                 vx, vy, vz, m_eject_dens, e_eject_therm, e_eject_kin,
                                                 p_feedback, metal_mass, 0, 0) == FAIL){
                  ENZO_FAIL("Error in grid->IndividualStarInjectFeedbackToGrid called from CheckTimeAction \n");
                }

              } // if on grid

            } // feedback injection method

            Temp = Temp->NextGridThisLevel;
          } // end while
        } //end loop over hierarchy


        for (level = MaximumRefinementLevel; level > 0; level--){
          LevelHierarchyEntry *Temp = LevelArray[level];
          while (Temp != NULL) {
            if (Temp->GridData->ProjectSolutionToParentGrid(*Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL){
              fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid\n");
              return FAIL;
            }
            Temp = Temp->NextGridThisLevel;
          }


        }

      } else {
        /* Done, turn it off (-1 in redshift indicates off). */

        TimeActionTime[i] = 0;
        TimeActionRedshift[i] = -1;

        /* Now apply to all grids. */
        for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      	  LevelHierarchyEntry *Temp = LevelArray[level];
  	  while (Temp != NULL) {
	    if (Temp->GridData->ApplyTimeAction(TimeActionType[i],
	    				   TimeActionParameter[i]) == FAIL) {
	      ENZO_FAIL("Errot in grid->ApplyTimeAction\n");

    	    }
	    Temp = Temp->NextGridThisLevel;
	  }
        }
      } // end time action type check for doing action
    } // end if time action is NOW
  } // end loop over time action types

  return SUCCESS;
}


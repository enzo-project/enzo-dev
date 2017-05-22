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

int CheckForTimeAction(LevelHierarchyEntry *LevelArray[],
		       TopGridData &MetaData)
{
 
  /* Declarations. */
 
  int i, level;
 
  /* Check for time against list of actions. */
 
  for (i = 0; i < MAX_TIME_ACTIONS; i++)
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


        // compute the SN location
        // for simplicity, this will be randomly chosen in a uniform disk
        float SNPosition[3] = {0.5, 0.5, 0.5};
        float rmax = IndividualStarICSupernovaR * Mpc / LengthUnits;
        float zmax = IndividualStarICSupernovaZ * Mpc / LengthUnits;

        if (MyProcessorNumber == ROOT_PROCESSOR){
          float r, theta;
          float z;

          r     = sqrt( rand()*1.0/RAND_MAX ) * rmax;
          theta = (rand()*1.0/RAND_MAX)  * 2.0 * pi;
          z     = (rand()*1.0/RAND_MAX) * 2.0 *zmax - zmax;

          SNPosition[0] += r * cos(theta);
          SNPosition[1] += r * sin(theta);
          SNPosition[2] += z;
        }


        CommunicationBroadcastValues(SNPosition, 3, ROOT_PROCESSOR);

        if (MyProcessorNumber == ROOT_PROCESSOR)
          printf("SNPosition = %f     %f     %f\n",SNPosition[0],SNPosition[1],SNPosition[2]);


        float m_eject   = 10.0 * 1.989E33 / MassUnits; // need to divide by dx before sending
        float E_thermal = IndividualStarSupernovaEnergy * 1.0E51 / (MassUnits*VelocityUnits*VelocityUnits); 

        for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++){
          LevelHierarchyEntry *Temp = LevelArray[level];
          while (Temp != NULL){

            // if explosion is on this grid (and grid is local), add feedback
            if (Temp->GridData->isLocal() && IsParticleFeedbackInGrid(SNPosition,
                                                                      IndividualStarFeedbackStencilSize,
                                                                      Temp)){

              float dx = float(Temp->GridData->ReturnCellWidth());
              float m_eject_dens = m_eject / (dx*dx*dx);
              float e_eject      = E_thermal / (dx*dx*dx);

              if ( Temp->GridData->IndividualStarInjectSphericalFeedback(NULL,
                                             SNPosition[0], SNPosition[1], SNPosition[2],
                                             m_eject_dens, e_eject, NULL, 0) == FAIL){
                ENZO_FAIL("Error in grid->function()");
              }

            }

            Temp = Temp->NextGridThisLevel;
          }
        }

        
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

        printf("finished applying explosion\n");
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
    } 



    }
 
  return SUCCESS;
}

/***********************************************************************
/
/  RECALIBRATE MBH FEEDBACK THERMAL RADIUS WHEN REQUESTED
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/  modified1: 
/
/  PURPOSE: This routine is implemented to apply MBH thermal feedback 
/           always to a constant mass, not to a constant volume.
/           Invoked when MBHFeedbackThermalRadius < 0. Enlarge Radius 
/           so that the thermal feedback affects the sphere enclosing 
/           a constant mass: n_ISM * 4pi/3 * (-MBHFeedbackThemralRadius)^3
/           assuming n_ISM = 1 Ms/pc^3 = 40.4 /cm^3
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

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
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int RecalibrateMBHFeedbackThermalRadius(FLOAT star_pos[], LevelHierarchyEntry *LevelArray[], 
					int level, float &Radius, 
					double &EjectaDensity, double &EjectaMetalDensity,
					double &EjectaThermalEnergy)
{

  const double pc = 3.086e18, Msun = 1.989e33, pMass = 1.673e-24, 
    gravConst = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13;

  float AvgVelocity[MAX_DIMENSION], MassEnclosed = 0, Metallicity = 0, ColdGasMass = 0;
  float OneOverRSquaredSum, initialRadius; 
  int i, l, dim, FirstLoop = TRUE, MBHFeedbackThermalRadiusTooSmall;
  LevelHierarchyEntry *Temp, *Temp2;
  FLOAT Time;

  /* This routine is invoked only when MBHFeedbackThermalRadius < 0 */

  if (MBHFeedback != 1 || MBHFeedbackThermalRadius >= 0) 
    return SUCCESS;

//  printf("RecalibrateMFTR: MyProcNum = %d, star_pos[] = %g, %g, %g\n", 
//	 MyProcessorNumber, star_pos[0], star_pos[1], star_pos[2]); 

  /* Set the units. */

  Temp = LevelArray[level];
  Time = Temp->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  /* Get cell width */

  int Rank, Dims[MAX_DIMENSION];
  float CellWidth;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  LevelArray[level]->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
  CellWidth = (RightEdge[0] - LeftEdge[0]) / (Dims[0] - 2*NumberOfGhostZones);

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    AvgVelocity[dim] = 0.0;

  /***********************************************************************
                           RECALIBRATE THE RADIUS
  ***********************************************************************/

  initialRadius = Radius;

  MBHFeedbackThermalRadiusTooSmall = (MBHFeedbackThermalRadius < 0);

  while (MBHFeedbackThermalRadiusTooSmall) { 
    Radius += CellWidth;
    MassEnclosed = 0;
    Metallicity = 0;
    ColdGasMass = 0;
    OneOverRSquaredSum = 0;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      AvgVelocity[dim] = 0.0;
    for (l = 0; l < MAX_DEPTH_OF_HIERARCHY; l++) {
      Temp = LevelArray[l];
      while (Temp != NULL) {

	/* Zero under subgrid field */

	if (FirstLoop == 1) {
	  Temp->GridData->
	    ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = LevelArray[l+1];
	  while (Temp2 != NULL) {
	    Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						     ZERO_UNDER_SUBGRID_FIELD);
	    Temp2 = Temp2->NextGridThisLevel;
	  }
	}

	/* Sum enclosed mass in this grid */

	if (Temp->GridData->GetEnclosedMass(star_pos, Radius, MassEnclosed, 
					    Metallicity, ColdGasMass, 
					    AvgVelocity, OneOverRSquaredSum) == FAIL) {
	  ENZO_FAIL("Error in GetEnclosedMass.");
	}

	Temp = Temp->NextGridThisLevel;

      } // END: Grids

      if (l == MAX_DEPTH_OF_HIERARCHY-1)
	FirstLoop = 0;

    } // END: level

    CommunicationAllSumValues(&Metallicity, 1);
    CommunicationAllSumValues(&MassEnclosed, 1);
    CommunicationAllSumValues(&ColdGasMass, 1);
    CommunicationAllSumValues(&OneOverRSquaredSum, 1);
    CommunicationAllSumValues(AvgVelocity, 3);

    if (MassEnclosed == 0) {
      ENZO_FAIL("RecalibrateMFTR: MassEnclosed = 0, something is wrong!\n");
    }

    Metallicity /= MassEnclosed;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      AvgVelocity[dim] /= MassEnclosed;

    MBHFeedbackThermalRadiusTooSmall = 
      (MassEnclosed < 4*M_PI/3.0*pow(-MBHFeedbackThermalRadius, 3)); 

//    fprintf(stdout, "RecalibrateMFTR: MassEnclosed = %g, MassEnclosedTarget = %g, Radius = %g\n", 
//	    MassEnclosed, 4*M_PI/3.0*pow(-MBHFeedbackThermalRadius, 3), Radius); 
    
  }  // ENDWHILE (too little mass)


  /* Reduce EjectaDensity */

  EjectaDensity *= pow(initialRadius/Radius, 3); 
  EjectaMetalDensity *= pow(initialRadius/Radius, 3); 

  /* Find EjectaThermalEnergy in a new expanded Radius; 
     check Grid_AddFeedbackSphere for different options */

  EjectaThermalEnergy *= pow(initialRadius/Radius, 3); 

#ifdef USE_ONE_OVER_RSQUARED
  EjectaThermalEnergy *= pow(initialRadius/Radius, 3)/OneOverRSquaredSum; 
#endif

#ifdef CONSTANT_SPECIFIC
  EjectaThermalEnergy *= 4*M_PI/3.0*pow(-MBHFeedbackThermalRadius, 3)/MassEnclosed;
#endif

//  fprintf(stdout, "RecalibrateMFTR: OneOverRSquaredSum = %g\n", OneOverRSquaredSum); 
//  fprintf(stdout, "RecalibrateMFTR: Radius = %g -> %g, EjectaThermalEnergy = %g -> %g, MassEnclosed = %g Msun\n", 
//	  initialRadius, Radius, EjectaThermalEnergy/pow(initialRadius/Radius,3), 
//	  EjectaThermalEnergy, MassEnclosed); 

  return SUCCESS;

}

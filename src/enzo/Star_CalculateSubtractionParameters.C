/***********************************************************************
/
/  CALCULATE PARAMETERS FOR SUBTRACTING MASS
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/
/  modified1: calculates parameters for subtracting mass from cells 
/             for MBH particles 
/
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

void Star::CalculateSubtractionParameters(LevelHierarchyEntry *LevelArray[], float &Radius, 
					  float RootCellWidth,
					  double &EjectaDensity,
					  float DensityUnits, float LengthUnits, 
					  float TemperatureUnits, float TimeUnits,
					  float VelocityUnits, float dtForThisStar)
{

  const double pc = 3.086e18, Msun = 1.989e33, Grav = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13, 
    k_b = 1.38e-16, m_h = 1.673e-24, c = 3.0e10, sigma_T = 6.65e-25, h=0.70;

  float StarLevelCellWidth, mdot, AccretedMass;
  float MassEnclosed = 0, Metallicity = 0, ColdGasMass = 0, AvgVelocity[MAX_DIMENSION];

  int igrid[MAX_DIMENSION], dim, l, index;
  int size=1, FirstLoop = TRUE;
  LevelHierarchyEntry *Temp, *Temp2;

  if (this->type != MBH || MBHAccretion <= 0) {
    fprintf(stderr, "You're not supposed to be coming in here; something's wrong!\n");
    return;
  }

  Radius = 0.0;
  EjectaDensity = 0.0;
  AccretedMass = 0.0;
  StarLevelCellWidth = RootCellWidth / powf(float(RefineBy), float(this->level));

  /* Find mdot and Radius */

  mdot = isnan(this->last_accretion_rate) ? 0.0 : this->last_accretion_rate;  
  
  Radius = MBHAccretionRadius * pc / LengthUnits;
  Radius = max(Radius, 2*StarLevelCellWidth);

  /* use negative EjectaDensity so the mass can be subtracted */

  EjectaDensity = -mdot * Msun * dtForThisStar * TimeUnits /   
    (4*M_PI/3.0 * pow(Radius*LengthUnits, 3)) / DensityUnits; 

//  fprintf(stdout, "star::CSP: mdot = %g, EjectaDensity=%lf, Radius=%g\n", mdot, EjectaDensity, Radius); 
//  fprintf(stdout, "star::CSP: star - old_mass = %lf  ->  new_mass = %lf\n", 
//          Mass, Mass + mdot * dtForThisStar * TimeUnits); 


  /* Below method currently not working because of CommunicationAllSumValues;
     if you want to use this approach, you may want to accordingly change 
     Grid_SubtractAccretedMassFromSphere.C as well */ 

#ifdef UNUSED   
  /* Find enclosed mass within Radius */

  MassEnclosed = 0;
  Metallicity = 0;
  ColdGasMass = 0;
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
      
      if (Temp->GridData->GetEnclosedMass(this, Radius, MassEnclosed, 
					  Metallicity, ColdGasMass, 
					  AvgVelocity) == FAIL) {
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
  CommunicationAllSumValues(AvgVelocity, 3);

  Metallicity /= MassEnclosed;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    AvgVelocity[dim] /= MassEnclosed;

//  fprintf(stdout, "star::CSP-2: MassEnclosed=%g\n", MassEnclosed); 
  
  /* Find ejecta characteristics (i.e. characteristics for the regions after mass subtraction) */

  AccretedMass = mdot * dtForThisStar * TimeUnits; //in Msun

  EjectaDensity = (float) 
    (double(Msun * (MassEnclosed - AccretedMass)) / 
     double(4*M_PI/3.0 * pow(Radius*LengthUnits, 3)) /
     DensityUnits);

  if (EjectaDensity < 0) {
    fprintf(stderr, "Your parameter (most likely MBHAccretionRadius) is set wrongly.\n");
    ENZO_FAIL("");
  }
#endif

  return;
}





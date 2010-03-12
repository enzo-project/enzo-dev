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
  float *temperature, density, old_mass, c_s, mu, number_density;

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

#ifdef ACCURATE_TEMPERATURE_AND_MU
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  /* Find fields: density, total energy, velocity1-3. */

  if (CurrentGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					      Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }

  /* Find Multi-species fields. */

  if (MultiSpecies)
    if (CurrentGrid->
	IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
			      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) 
	== FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      ENZO_FAIL("");
    }

  /* Find temperature */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= CurrentGrid->GridDimension[dim];
    igrid[dim] = (int) (pos[dim] - CurrentGrid->GridLeftEdge[dim]) /
      CurrentGrid->CellWidth[0][0];
  }

  temperature = new float[size];
  if (CurrentGrid->ComputeTemperatureField(temperature) == FAIL) {
    fprintf(stderr, "Error in ComputeTemperatureField.\n");
    ENZO_FAIL("");
  }

  /* Calculate mu inside cell */

  index = 
    ((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
     igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
    igrid[0] + CurrentGrid->GridStartIndex[0];
  density = CurrentGrid->BaryonField[DensNum][index];

  if (MultiSpecies == 0) {
    number_density = density * DensityUnits / (DEFAULT_MU * m_h);
    mu = DEFAULT_MU;
  } else {
    number_density = 
      CurrentGrid->BaryonField[HINum][index] + 
      CurrentGrid->BaryonField[HIINum][index] +
      CurrentGrid->BaryonField[DeNum][index] +
      0.25 * (CurrentGrid->BaryonField[HeINum][index] +
	      CurrentGrid->BaryonField[HeIINum][index] +
	      CurrentGrid->BaryonField[HeIIINum][index]);
    if (MultiSpecies > 1)
      number_density += 
	CurrentGrid->BaryonField[HMNum][index] +
	0.5 * (CurrentGrid->BaryonField[H2INum][index] +
	       CurrentGrid->BaryonField[H2IINum][index]);
    mu = density / number_density;
  }

  /* If requested, fix the temperature */

  if (this->type == MBH && (MBHAccretion == 2 || MBHAccretion == 12)) {
    temperature[index] = MBHAccretionFixedTemperature;   
  }

  /* Calculate c_s */

  c_s = sqrt(Gamma * k_b * temperature[index] / (mu * m_h));

#endif


  /* Find mdot and Radius;
     for negative MBHAccretionRadius, then use Bondi accretion radius instead */

  /* When calculating Bondi radius, let's not bother to get accurate temperature and mu 
     just use the MBHAccretionFixedTemperature and default mu */

  c_s = sqrt(Gamma * k_b * MBHAccretionFixedTemperature / (DEFAULT_MU * m_h));
  old_mass = (float)(this->Mass);
  mdot = isnan(this->last_accretion_rate) ? 0.0 : this->last_accretion_rate;  
  
  if (MBHAccretionRadius > 0) 
    Radius = MBHAccretionRadius * pc / LengthUnits;
  else {
    float safetyfactor = 2.0;
    Radius = safetyfactor * 2.0 * Grav * old_mass * Msun / (c_s * c_s) / LengthUnits;
  }

  Radius = max(Radius, 2*StarLevelCellWidth);

  /* use negative EjectaDensity so the mass can be subtracted */

  EjectaDensity = -mdot * Msun * dtForThisStar * TimeUnits /   
    (4*M_PI/3.0 * pow(Radius*LengthUnits, 3)) / DensityUnits; 

//  fprintf(stdout, "star::CSP: mdot = %g, EjectaDensity=%lf, Radius=%g\n", 
//	  mdot, EjectaDensity, Radius);  
//  fprintf(stdout, "star::CSP: star - old_mass = %lf  ->  new_mass = %lf\n", 
//          Mass - mdot * dtForThisStar * TimeUnits, Mass); 


#ifdef UNUSED   
  /* Below method currently not working because of CommunicationAllSumValues;
     if you want to use this approach, you may want to accordingly change 
     Grid_SubtractAccretedMassFromSphere.C as well */ 

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


#ifdef ACCURATE_TEMPERATURE_AND_MU
  delete [] temperature;
#endif 

  return;
}





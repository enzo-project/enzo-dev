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

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);   

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int CalculateSubtractionParameters(LevelHierarchyEntry *LevelArray[], int level, FLOAT star_pos[],
				   double star_mass, float star_last_accretion_rate, star_type star_type,
				   grid *star_CurrentGrid,
				   float StarLevelCellWidth, float dtForThisStar,
				   float &Radius, double &Subtraction)
{

  const double pc = 3.086e18, Msun = 1.989e33, Grav = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13, 
    k_b = 1.38e-16, m_h = 1.673e-24, c = 3.0e10, sigma_T = 6.65e-25, h=0.70;

  float mdot, AccretedMass, SafetyFactor;
  float MassEnclosed = 0, Metallicity = 0, ColdGasMass = 0, OneOverRSquaredSum, AvgVelocity[MAX_DIMENSION];
  float *temperature, density, old_mass, mu, number_density;
  
  int igrid[MAX_DIMENSION], dim, l, index, c_s;
  int size=1, FirstLoop = TRUE;
  LevelHierarchyEntry *Temp, *Temp2;
  FLOAT Time;

  /* This routine is invoked only when there is a MBH accretion */

  if (star_type != MBH || MBHAccretion <= 0 || MBHAccretionRadius == 0) {
    fprintf(stderr, "You're not supposed to be coming in here; something's wrong!\n");
    return SUCCESS;
  }

  /* Set the units. */

  Temp = LevelArray[level];
  Time = Temp->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  Radius = 0.0;
  Subtraction = 0.0;
  AccretedMass = 0.0;



#ifdef UNUSED
  /* If star_CurrentGrid is on the current processor, 
     find mu, temperature[], and c_s so that we can eventually calculate Bondi radius.
     One has to broadcast c_s to other processors */
  
  if (star_CurrentGrid->ReturnProcessorNumber() ==  MyProcessorNumber) {
    
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
    
    /* Find fields: density, total energy, velocity1-3. */
    
    if (star_CurrentGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
						     Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      ENZO_FAIL("");
    }

    /* Find Multi-species fields. */
    
    if (MultiSpecies)
      if (star_CurrentGrid->
	  IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
				HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) 
	  == FAIL) {
	fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
	ENZO_FAIL("");
      }
    
    /* Find temperature */
    
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      size *= star_CurrentGrid->ReturnGridDimension()[dim];
      igrid[dim] = (int) ((star_pos[dim] - star_CurrentGrid->ReturnGridLeftEdge()[dim]) /
			  StarLevelCellWidth);
    }
    
    temperature = new float[size];
    if (star_CurrentGrid->ComputeTemperatureField(temperature) == FAIL) {
      fprintf(stderr, "Error in ComputeTemperatureField.\n");
      ENZO_FAIL("");
    }

    /* Calculate mu inside cell */
    
    index = 
      ((igrid[2] + star_CurrentGrid->ReturnGridStartIndex()[2]) * star_CurrentGrid->ReturnGridDimension()[1] + 
       igrid[1] + star_CurrentGrid->ReturnGridStartIndex()[1]) * star_CurrentGrid->ReturnGridDimension()[0] + 
      igrid[0] + star_CurrentGrid->ReturnGridStartIndex()[0];
    density = star_CurrentGrid->ReturnBaryonField(DensNum)[index];
    
    if (MultiSpecies == 0) {
      number_density = density * DensityUnits / (Mu * m_h);
      mu = Mu;
    } else {
      number_density = 
	star_CurrentGrid->ReturnBaryonField(HINum)[index] + 
	star_CurrentGrid->ReturnBaryonField(HIINum)[index] +
	star_CurrentGrid->ReturnBaryonField(DeNum)[index] +
	0.25 * (star_CurrentGrid->ReturnBaryonField(HeINum)[index] +
		star_CurrentGrid->ReturnBaryonField(HeIINum)[index] +
		star_CurrentGrid->ReturnBaryonField(HeIIINum)[index]);
      if (MultiSpecies > 1)
	number_density += 
	  star_CurrentGrid->ReturnBaryonField(HMNum)[index] +
	  0.5 * (star_CurrentGrid->ReturnBaryonField(H2INum)[index] +
		 star_CurrentGrid->ReturnBaryonField(H2IINum)[index]);
      mu = density / number_density;
    }

    /* If requested, fix the temperature */
    
    if (star_type == MBH && (MBHAccretion == 2 || MBHAccretion == 12)) {
      temperature[index] = MBHAccretionFixedTemperature;   
    }

    /* Calculate c_s, Here I boldly typecast c_s into an integer 
       because only integer can be broadcasted at the moment;
       anyway, this shouldn't make a huge difference overall. */
    
    c_s = (int)(sqrt(Gamma * k_b * temperature[index] / (mu * m_h)));

    fprintf(stdout, "index = %d, temp = %g, mu = %g, density = %g, number_density = %g, c_s = %d\n",
	    index, temperature[index], mu, density, number_density, c_s);  

    delete [] temperature;

  }  // END OF if (star_CurrentGrid->ReturnProcessorNumber() ==  MyProcessorNumber) 


  /* Now broadcast c_s to other processors */

  CommunicationBroadcastValue(&c_s, star_CurrentGrid->ReturnProcessorNumber());

  fprintf(stdout, "MyProc = %d, star_CurrentGrid_ProcNum = %d, c_s = %d\n", 
	  MyProcessorNumber, star_CurrentGrid->ReturnProcessorNumber(), c_s);
#endif // UNUSED



#define APPROXIMATE_TEMPERATURE_AND_MU
#ifdef APPROXIMATE_TEMPERATURE_AND_MU
  /* When calculating Bondi radius, let's not bother to get accurate temperature and mu 
     just use the MBHAccretionFixedTemperature and default mu */

    c_s = (int)(sqrt(Gamma * k_b * MBHAccretionFixedTemperature / (Mu * m_h)));
#endif



  /* Find mdot and Radius (from c_s);
     for negative MBHAccretionRadius, then use Bondi accretion radius instead */

  old_mass = (float)(star_mass);
  mdot = isnan(star_last_accretion_rate) ? 0.0 : star_last_accretion_rate;  
  AccretedMass = mdot * dtForThisStar * TimeUnits; //in Msun

  
  if (MBHAccretionRadius > 0) 
    Radius = MBHAccretionRadius * pc / LengthUnits;
  else {
    SafetyFactor = -MBHAccretionRadius;
    Radius = SafetyFactor * 2.0 * Grav * old_mass * Msun / (c_s * c_s) / LengthUnits;
  }

  Radius = min(max(Radius, 4*StarLevelCellWidth), 100*StarLevelCellWidth);


  while (MassEnclosed <= 0) {

    /* Find enclosed mass within Radius */
    
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
      ENZO_FAIL("CSP: MassEnclosed = 0; something is wrong!\n");
    }

    Metallicity /= MassEnclosed;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      AvgVelocity[dim] /= MassEnclosed;

//  fprintf(stdout, "CSP: MassEnclosed, OneOverRSquaredSum = %g %g \n", 
//	  MassEnclosed, OneOverRSquaredSum); 
  
  }

  /* here Subtraction is the mass ratio of 
     accreted mass (thus, to-be-subtracted mass) to total mass in the sphere */

  Subtraction = (double)(AccretedMass) / (double)(MassEnclosed); 

  if (Subtraction < 0) {
    ENZO_FAIL("CSP: parameters (most likely MBHAccretionRadius) set wrongly.\n");
  }


#ifdef SUBTRACTION_UNIFORM
  /* find negative Subtraction so the mass can be uniformly subtracted out of the sphere;
     here Subtraction is the density that should be subtracted.
     if one wants to use this, change Grid_SAMFS.C accordingly. */

  Subtraction = -AccretedMass * Msun /   
    (4*M_PI/3.0 * pow(Radius*LengthUnits, 3)) / DensityUnits; 
#endif



//  fprintf(stdout, "CSP: mdot = %g, Subtraction=%g, Radius=%g\n", mdot, Subtraction, Radius);   
//  fprintf(stdout, "CSP: star - old_mass = %lf  ->  new_mass = %lf\n", star_mass - AccretedMass, star_mass);  

  return SUCCESS;
}





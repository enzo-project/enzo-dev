/***********************************************************************
/
/  CALCULATE FEEDBACK SPHERE PARAMETERS
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: Ji-hoon Kim
/             July, 2009
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
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "StarParticleData.h"

void Star::CalculateFeedbackParameters(float &Radius, 
				       float RootCellWidth,
				       float SNe_dt, double &EjectaDensity,
				       double &EjectaThermalEnergy,
				       double &EjectaMetalDensity,
				       float DensityUnits, float LengthUnits, 
				       float TemperatureUnits, float TimeUnits,
				       float VelocityUnits)
{

  // Parameters for the Stroemgen sphere in Whalen et al. (2004)
  const float	BirthRadius	  = 50;		// pc
  const float	WhalenTemperature = 20000;	// K
  const float	WhalenDensity	  = 1;	// cm^-3
  const float	WhalenMaxVelocity = 35;		// km/s

  const double pc = 3.086e18, Msun = 1.989e33, Grav = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13, 
    k_b = 1.38e-16, m_h = 1.673e-24, c = 3.0e10, sigma_T = 6.65e-25;

  float StarLevelCellWidth;
  double EjectaVolume, SNEnergy, HeliumCoreMass, Delta_SF;

  int igrid[MAX_DIMENSION], dim, index;
  int size=1;
  float c_s, mu, number_density, old_mass, delta_mass, mdot, mdot_Edd, v_rel, dvel;
  float *temperature, density;

  Radius = 0.0;
  EjectaDensity = 0.0;
  EjectaThermalEnergy = 0.0;
  EjectaMetalDensity = 0.0;
  StarLevelCellWidth = RootCellWidth / powf(float(RefineBy), float(this->level));

  switch (this->FeedbackFlag) {
  case SUPERNOVA:  // pair-instability SNe
    Radius = PopIIISupernovaRadius * pc / LengthUnits;
    Radius = max(Radius, 3.5*StarLevelCellWidth);
    EjectaVolume = 4.0/3.0 * 3.14159 * pow(PopIIISupernovaRadius*pc, 3);
    EjectaDensity = Mass * Msun / EjectaVolume / DensityUnits;
    HeliumCoreMass = (13./24.) * (Mass - 20);
    SNEnergy = (5.0 + 1.304 * (HeliumCoreMass - 64)) * 1e51;
    EjectaMetalDensity = HeliumCoreMass * Msun / EjectaVolume / 
      DensityUnits;
    EjectaThermalEnergy = SNEnergy / (Mass * Msun) / VelocityUnits /
      VelocityUnits;

    // Exaggerate influence radius because the blastwave will enter
    // into some of the surrounding parent grids within the next
    // timestep if we inject the energy into a small radius.
    Radius *= 8.0;
    break;

  case STROEMGREN:
    Radius = BirthRadius * pc / LengthUnits;
    Radius = max(Radius, 1.5*StarLevelCellWidth);
    EjectaDensity = WhalenDensity * m_h / DensityUnits;
    EjectaThermalEnergy =
      WhalenTemperature / (TemperatureUnits * (Gamma-1.0) * 0.6);
    break;

  case FORMATION:
    Radius = 0;
    break;

  case CONT_SUPERNOVA:
    // Inject energy into a sphere
    Radius = StarClusterSNRadius * pc / LengthUnits;
    Radius = max(Radius, 2*StarLevelCellWidth);

    // Release SNe energy constantly over 16 Myr (t = 4-20 Myr), which is defined in Star_SetFeedbackFlag.C.
    Delta_SF = Mass * SNe_dt * TimeUnits / (16.0*Myr);
    EjectaVolume = 4.0/3.0 * 3.14159 * pow(Radius*LengthUnits, 3);
    EjectaDensity = Delta_SF * Msun / EjectaVolume / DensityUnits; 
    EjectaMetalDensity = EjectaDensity * StarMetalYield;
    EjectaThermalEnergy = StarClusterSNEnergy / Msun /
      (VelocityUnits * VelocityUnits);
    break;

  case MBH_THERMAL:
    if (this->type != MBH || this->CurrentGrid ==  NULL) break;

    /* Using the code snippets adopted from Star_CalculateMassAccretion.C (case LOCAL_ACCRETION)
       estimate the Bondi accretion rate.  - Ji-hoon Kim */

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

    if (CurrentGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
						Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      ENZO_FAIL("");
    }

    if (MultiSpecies)
      if (CurrentGrid->
	  IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
				HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) 
	  == FAIL) {
	fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
	ENZO_FAIL("");
      }

   
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

    // Calculate gas density inside cell
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
    c_s = sqrt(Gamma * k_b * temperature[index] / (mu * m_h));
    old_mass = this->Mass; //Msun

    // Calculate gas relative velocity (cm/s)
    v_rel = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      delta_vel[dim] = vel[dim] - CurrentGrid->BaryonField[Vel1Num+dim][index];
      v_rel += delta_vel[dim] * delta_vel[dim];
    }
    v_rel = sqrt(v_rel) * VelocityUnits;

    // Calculate Bondi accretion rate in Msun/s 
    mdot = 4.0 * PI * Grav*Grav * (old_mass * old_mass * Msun) * 
      (density * DensityUnits) / pow(c_s * c_s + v_rel * v_rel, 1.5);

    /* end of the code snippets from Star_CalculateMassAccretion.C */

    // Calculate Eddington accretion rate in Msun/s; the Eddington limit for feedback
    mdot_Edd = 4.0 * PI * Grav * old_mass * m_h /
      MBHFeedbackRadiativeEfficiency / sigma_T / c; 

    // Inject energy into a sphere
    Radius = MBHFeedbackRadius * pc / LengthUnits;
    Radius = max(Radius, 2*StarLevelCellWidth);

    // Release MBH-AGN thermal energy constantly. Here no mass is released.
    EjectaVolume = 4.0/3.0 * 3.14159 * pow(Radius*LengthUnits, 3);
    EjectaDensity = 0.0;
    EjectaMetalDensity = 0.0; 

    /* Now calculate the feedback parameter based on mdot estimated above.  - Ji-hoon Kim 
       For CONT_SUPERNOVA, the unit of EjectaThermalEnergy was ergs/g, 
       but here for MBH_THERMAL, the unit of EjectaThermalEnergy is ergs/cm^3.
       This is because EjectaDensity = 0 in this case; see Grid_AddFeedbackSphere.C  - Ji-hoon Kim */

    EjectaThermalEnergy = MBHFeedbackThermalCoupling * MBHFeedbackRadiativeEfficiency * 
      min(mdot, mdot_Edd) * Msun * c * c * CurrentGrid->dtFixed * TimeUnits / EjectaVolume / 
      DensityUnits / (VelocityUnits * VelocityUnits) ; //Eq.(34) in Springel (2005) 

#define SEDOV_TEST
#ifdef SEDOV_TEST
    EjectaThermalEnergy = 1.0e50 / EjectaVolume / DensityUnits / (VelocityUnits * VelocityUnits);  
    
    // For the continuous energy injection case (variation of Sedov test)
    /*
    EjectaThermalEnergy = 1.0e52 * CurrentGrid->dtFixed * TimeUnits / 9e14
      / EjectaVolume / DensityUnits / (VelocityUnits * VelocityUnits);  
    */
#endif

    break;

  } // ENDSWITCH FeedbackFlag
  
  return;
}





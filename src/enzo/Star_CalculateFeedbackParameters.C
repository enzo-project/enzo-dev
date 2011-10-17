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

int search_lower_bound(float *arr, float value, int low, int high, 
		       int total);

void Star::CalculateFeedbackParameters(float &Radius, 
				       float RootCellWidth,
				       float SNe_dt, double &EjectaDensity,
				       double &EjectaThermalEnergy,
				       double &EjectaMetalDensity,
				       float DensityUnits, float LengthUnits, 
				       float TemperatureUnits, float TimeUnits,
				       float VelocityUnits, float dtForThisStar,
				       FLOAT Time)
{

  // Parameters for the Stroemgen sphere in Whalen et al. (2004)
  const float	BirthRadius	  = 50;		// pc
  const float	WhalenTemperature = 20000;	// K
  const float	WhalenDensity	  = 1;	        // cm^-3
  const float	WhalenMaxVelocity = 35;		// km/s

  const double pc = 3.086e18, Msun = 1.989e33, Grav = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13, 
    k_b = 1.38e-16, m_h = 1.673e-24, c = 3.0e10, sigma_T = 6.65e-25, h=0.70;

  const float TypeIILowerMass = 11, TypeIIUpperMass = 40;
  const float PISNLowerMass = 140, PISNUpperMass = 260;

  // From Nomoto et al. (2006)
  const float HypernovaMass[] = {19.99, 25, 30, 35, 40.01};  // Msun
  const float HypernovaMetals[] = {3.36, 3.53, 5.48, 7.03, 8.59}; // Msun
  const float HypernovaEnergy[] = {10, 10, 20, 25, 30};  // 1e51 erg

  float StarLevelCellWidth, tdyn, frac;
  double EjectaVolume, SNEnergy, HeliumCoreMass, Delta_SF, MetalMass;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum;

  int igrid[MAX_DIMENSION], dim, index, bin;
  int size=1;
  float mdot;

  Radius = 0.0;
  EjectaDensity = 0.0;
  EjectaThermalEnergy = 0.0;
  EjectaMetalDensity = 0.0;
  StarLevelCellWidth = RootCellWidth / powf(float(RefineBy), float(this->level));

  switch (this->FeedbackFlag) {
  case SUPERNOVA:  // Single thermal bubble of SN feedback
    Radius = PopIIISupernovaRadius * pc / LengthUnits;
    Radius = max(Radius, 3.5*StarLevelCellWidth);
    EjectaVolume = 4.0/3.0 * 3.14159 * pow(PopIIISupernovaRadius*pc, 3);
    EjectaDensity = Mass * Msun / EjectaVolume / DensityUnits;

    // pair-instability SNe
    if (this->Mass >= PISNLowerMass && this->Mass <= PISNUpperMass) {
      HeliumCoreMass = (13./24.) * (Mass - 20);
      SNEnergy = (5.0 + 1.304 * (HeliumCoreMass - 64)) * 1e51;
      EjectaMetalDensity = HeliumCoreMass * Msun / EjectaVolume / 
	DensityUnits;
    } 
    
    // Type II SNe
    else if (this->Mass >= TypeIILowerMass && this->Mass <= TypeIIUpperMass) {
      if (this->Mass < 20.0) { // Normal Type II
	SNEnergy = 1e51;
	MetalMass = 0.1077 + 0.3383 * (this->Mass - 11.0);  // Fit to Nomoto+06
      } else { // Hypernova (should we add the "failed" SNe?)
	bin = search_lower_bound((float*)HypernovaMass, this->Mass, 0, 5, 5);
	frac = (HypernovaMass[bin+1] - this->Mass) / 
	  (HypernovaMass[bin+1] - HypernovaMass[bin]);
	SNEnergy = 1e51 * (HypernovaEnergy[bin] + 
			   frac * (HypernovaEnergy[bin+1] - HypernovaEnergy[bin]));
	MetalMass = (HypernovaMetals[bin] + 
		     frac * (HypernovaMetals[bin+1] - HypernovaMetals[bin]));
      }
      EjectaMetalDensity = MetalMass * Msun / EjectaVolume / DensityUnits;
    }
    EjectaThermalEnergy = SNEnergy / (Mass * Msun) / VelocityUnits /
      VelocityUnits;

    // Exaggerate influence radius because the blastwave will enter
    // into some of the surrounding parent grids within the next
    // timestep if we inject the energy into a small radius.
    Radius *= 1.0;
    
#define DEBUG
#ifdef DEBUG
    if (MyProcessorNumber == ROOT_PROCESSOR)
      if (Radius > 10*PopIIISupernovaRadius * pc / LengthUnits) {
	printf("WARNING: Inserting PISN into a large radius!  %g pc\n",
	       Radius * LengthUnits / pc);
      }
#endif    

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
    //Delta_SF = StarMassEjectionFraction * Mass * SNe_dt * TimeUnits / (16.0*Myr);
    if (StarClusterUnresolvedModel) {
      // lifetime = 5*tdyn
      tdyn = 0.2 * LifeTime;
      frac = (Time-BirthTime) / tdyn;
      Delta_SF = dtForThisStar * StarMassEjectionFraction * (Mass/tdyn) *
	frac * exp(-frac);
//      if (debug)
//	printf("Star %d: Delta_SF = %g, Mass = %g, frac = %g (%f %f)\n", 
//	       Identifier, Delta_SF, Mass, frac, Time, BirthTime);
    } else {
      Delta_SF = StarMassEjectionFraction * Mass * dtForThisStar * 
	TimeUnits / (16.0*Myr);
    }
    EjectaVolume = 4.0/3.0 * 3.14159 * pow(Radius*LengthUnits, 3);   
    EjectaDensity = Delta_SF * Msun / EjectaVolume / DensityUnits;   
    EjectaMetalDensity = EjectaDensity * StarMetalYield;
    EjectaThermalEnergy = StarClusterSNEnergy / Msun /   
      (VelocityUnits * VelocityUnits);
    break;

  case MBH_THERMAL:
    if (this->type != MBH) 
      ENZO_FAIL("Applying MBH_THERMAL feedback to non-MBH particle!");

    /* find mdot */
    mdot = isnan(this->last_accretion_rate) ? 0.0 : this->last_accretion_rate;  
    
    /* Inject energy into a sphere */
    Radius = MBHFeedbackThermalRadius * pc / LengthUnits;
    Radius = max(Radius, 2*StarLevelCellWidth);

    /* Only EjectaVolume is in physical units; all others are in code units. */
    EjectaVolume = 4.0/3.0 * PI * pow(Radius*LengthUnits, 3);  
    EjectaDensity = mdot * Msun * dtForThisStar * TimeUnits * MBHFeedbackMassEjectionFraction /
      EjectaVolume / DensityUnits; 
    EjectaMetalDensity = EjectaDensity * MBHFeedbackMetalYield; 

    /* When injected energy is uniform throughout the volume;
       The unit of EjectaThermalEnergy is ergs/cm3 = (cm^2/s^2) * (g/cm3).
       This value will be recalibrated in RecalibrateMFTR */
    EjectaThermalEnergy = MBHFeedbackEnergyCoupling * MBHFeedbackRadiativeEfficiency * 
      mdot * Msun * c * c * dtForThisStar * TimeUnits / 
      EjectaVolume / DensityUnits / (VelocityUnits * VelocityUnits); 

#ifdef CONSTANT_SPECIFIC
    /* When injected energy is proportional to the cell mass;
       The unit of EjectaThermalEnergy is ergs/g = cm^2/s^2. */
    EjectaThermalEnergy = MBHFeedbackEnergyCoupling * MBHFeedbackRadiativeEfficiency * 
      mdot * Msun * c * c * dtForThisStar * TimeUnits / 
      (4.0/3.0 * PI * pow(-MBHFeedbackThermalRadius, 3) * Msun) / (VelocityUnits * VelocityUnits);
#endif    

#ifdef SEDOV_TEST
    // For Sedov test, here the unit of EjectaThermalEnergy is ergs/cm^3.  
//    EjectaDensity = 0.0;
//    EjectaThermalEnergy = 1.0e50 /
//      EjectaVolume / DensityUnits / (VelocityUnits * VelocityUnits);  
    
    // For Ostriker & McKee test (the continuous energy injection case, variation of Sedov test)
    EjectaDensity = 0.0;
    EjectaThermalEnergy = 1.0e40 * dtForThisStar * TimeUnits /
      EjectaVolume / DensityUnits / (VelocityUnits * VelocityUnits);  
#endif

    if (isnan(EjectaThermalEnergy)) EjectaThermalEnergy = 0.0;

    break;

  case MBH_JETS:
    if (this->type != MBH) 
      ENZO_FAIL("Applying MBH_JETS feedback to non-MBH particle!");

    /* find mdot */
    mdot = isnan(this->last_accretion_rate) ? 0.0 : this->last_accretion_rate;  
    
    /* Inject energy into a sphere */
    Radius = MBHFeedbackThermalRadius * pc / LengthUnits;
    Radius = max(Radius, 2*StarLevelCellWidth);

    /* Release MBH-AGN thermal energy constantly. 
       Only EjectaVolume is in physical units; all others are in code units. */
    EjectaVolume = 4.0/3.0 * PI * pow(Radius*LengthUnits, 3);  
    EjectaDensity = mdot * Msun * dtForThisStar * TimeUnits * MBHFeedbackMassEjectionFraction /
      EjectaVolume / DensityUnits; 
    EjectaMetalDensity = EjectaDensity * MBHFeedbackMetalYield; 

    /* Now calculate the feedback parameter based on mdot estimated above.  
       The unit of EjectaThermalEnergy is ergs/g = cm^2/s^2. */
    EjectaThermalEnergy = MBHFeedbackEnergyCoupling * MBHFeedbackRadiativeEfficiency * 
      mdot * Msun * c * c * dtForThisStar * TimeUnits / 
      (EjectaDensity * DensityUnits) / EjectaVolume / (VelocityUnits * VelocityUnits);
    if (isnan(EjectaThermalEnergy)) EjectaThermalEnergy = 0.0;

    
    break;

  } // ENDSWITCH FeedbackFlag
  
//    fprintf(stdout, "star::CFP:  EjectaThermalEnergy = %g, EjectaDensity = %g, 
//                Radius = %g, mdot = %g, dtForThisStar = %g\n", 
//    	    EjectaThermalEnergy, EjectaDensity, Radius, mdot, dtForThisStar);  

  return;
}





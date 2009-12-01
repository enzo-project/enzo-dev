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
				       float VelocityUnits, float dtForThisStar)
{

  // Parameters for the Stroemgen sphere in Whalen et al. (2004)
  const float	BirthRadius	  = 50;		// pc
  const float	WhalenTemperature = 20000;	// K
  const float	WhalenDensity	  = 1;	        // cm^-3
  const float	WhalenMaxVelocity = 35;		// km/s

  const double pc = 3.086e18, Msun = 1.989e33, Grav = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13, 
    k_b = 1.38e-16, m_h = 1.673e-24, c = 3.0e10, sigma_T = 6.65e-25, h=0.70;

  float StarLevelCellWidth;
  double EjectaVolume, SNEnergy, HeliumCoreMass, Delta_SF;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum;

  int igrid[MAX_DIMENSION], dim, index;
  int size=1;
  float c_s, mu, number_density, old_mass, delta_mass, mdot, mdot_Edd, mdot_UpperLimit, v_rel, dvel;
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
  case MBH_JETS:
    if (this->type != MBH || this->CurrentGrid ==  NULL) break;

    /**********************************************
                     CALCULATE mdot
    **********************************************/

    /* [1] Method 1
       Use DeltaMass calculated in the previous timestep in Star_CalculateMassAccretion.C
       This turned out to be unsuccessful, however. */

    // mdot = this->DeltaMass / CurrentGrid->dtFixed / TimeUnits; // in Msun/sec

    /* [2] Method 2
       Use code snippets from Star_CalculateMassAccretion.C. (many comments omitted) 
       This is redundant, but for now, works great!  */

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
    old_mass = this->Mass;

    // Calculate gas relative velocity (cm/s)
    v_rel = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      delta_vel[dim] = vel[dim] - CurrentGrid->BaryonField[Vel1Num+dim][index];
      v_rel += delta_vel[dim] * delta_vel[dim];
    }
    v_rel = sqrt(v_rel) * VelocityUnits;

    // Calculate accretion rate in Msun/s
    mdot = 4.0 * PI * Grav*Grav * (old_mass * old_mass * Msun) * 
      (density * DensityUnits) / pow(c_s * c_s + v_rel * v_rel, 1.5);

    // Don't take out too much mass suddenly; mdot should leave at least 75% of the gas in the grids.
    mdot_UpperLimit = 0.25 * density * DensityUnits * 
      pow(CurrentGrid->CellWidth[0][0]*LengthUnits, 3.0) / Msun / 
      (CurrentGrid->dtFixed) / TimeUnits;
    mdot = min(mdot, mdot_UpperLimit);

    // No accretion if the BH is in some low-density and cold cell.
    if (density < tiny_number || temperature[index] < 10 || isnan(mdot) || MBHAccretion != 1)
      mdot = 0.0;

    if (this->type == MBH) { 
      mdot *= MBHAccretingMassRatio;

      mdot_Edd = 4.0 * PI * Grav * old_mass * m_h /
	MBHFeedbackRadiativeEfficiency / sigma_T / c; 

      mdot = min(mdot, mdot_Edd); 
    }

    /* End of the code snippets from Star_CalculateMassAccretion.C */



    // Inject energy into a sphere
    Radius = MBHFeedbackThermalRadius * pc / LengthUnits;
    Radius = max(Radius, 2*StarLevelCellWidth);

    /* Release MBH-AGN thermal energy constantly. 
       Only EjectaVolume is in physical units; all others are in code units. */
    EjectaVolume = 4.0/3.0 * PI * pow(Radius*LengthUnits, 3);  
    EjectaDensity = mdot * Msun * dtForThisStar * TimeUnits * MBHFeedbackMassEjectionFraction /
      EjectaVolume / DensityUnits; 
    EjectaMetalDensity = EjectaDensity * MBHFeedbackMetalYield; //very fiducial

    /* Now calculate the feedback parameter based on mdot estimated above.  
       The unit of EjectaThermalEnergy is ergs/g = cm^2/s^2.  - Ji-hoon Kim  Aug.2009 */
    EjectaThermalEnergy = MBHFeedbackThermalCoupling * MBHFeedbackRadiativeEfficiency * 
      mdot * Msun * c * c * dtForThisStar * TimeUnits / 
      (EjectaDensity * DensityUnits) / EjectaVolume / (VelocityUnits * VelocityUnits) ; //Eq.(34) in Springel (2005) 
    

#define NOT_SEDOV_TEST
#ifdef SEDOV_TEST
    /* For Sedov test, here the unit of EjectaThermalEnergy is ergs/cm^3.  
       This is because EjectaDensity = 0 in this case; see Grid_AddFeedbackSphere.C  */
    /*
    EjectaDensity = 0.0;
    EjectaThermalEnergy = 1.0e50 /
      EjectaVolume / DensityUnits / (VelocityUnits * VelocityUnits);  
    */
    
    // For Ostriker & McKee test (the continuous energy injection case, variation of Sedov test)
    EjectaDensity = 0.0;
    EjectaThermalEnergy = 1.0e40 * dtForThisStar * TimeUnits /
      EjectaVolume / DensityUnits / (VelocityUnits * VelocityUnits);  
#endif

    //    fprintf(stdout, "star::CFP:  EjectaThermalEnergy = %g, EjectaDensity = %g, 
    //            Radius = %g, mdot = %g, dtForThisStar = %g\n", 
    //	    EjectaThermalEnergy, EjectaDensity, Radius, mdot, dtForThisStar); 
    break;

  } // ENDSWITCH FeedbackFlag
  
  return;
}





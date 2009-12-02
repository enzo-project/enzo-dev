/***********************************************************************
/
/  CALCULATE MASS ACCRETION ONTO STAR PARTICLE
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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int Star::CalculateMassAccretion(void)
{

  if ((this->type != BlackHole && ABS(this->type) != MBH) || 
      (this->CurrentGrid == NULL)) 
    return SUCCESS;

  const double Grav = 6.673e-8, k_b = 1.38e-16, m_h = 1.673e-24;
  const double Msun = 1.989e33, yr = 3.1557e7, sigma_T = 6.65e-25, c = 3.0e10;
  const int AccretionType = LOCAL_ACCRETION;
  FLOAT time = CurrentGrid->OldTime;

  if (time <= 0)
    time = CurrentGrid->Time - CurrentGrid->dtFixed;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, time);

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

  int igrid[MAX_DIMENSION], dim, index, size = 1;
  float c_s, mu, number_density, old_mass, delta_mass, mdot, mdot_UpperLimit, mdot_Edd, v_rel, dvel;
  float *temperature, density;

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

  /* Reset the accretion rate (DeltaMass) */

  this->ResetAccretion();

  if (AccretionType == LOCAL_ACCRETION) {

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
      /* For MBH, MBHAccretingMassRatio is implemented; when we resolve Bondi radius, 
	 the local density used to calculate mdot can be higher than what was supposed 
	 to be used, resulting in the overestimation of mdot. MBHAccretingMassRatio 
	 (=< 1) can be used to fix this.  -Ji-hoon Kim, Sep.2009 */
      mdot *= MBHAccretingMassRatio;

      // Calculate Eddington accretion rate in Msun/s; the Eddington limit for feedback
      mdot_Edd = 4.0 * PI * Grav * old_mass * m_h /
	MBHFeedbackRadiativeEfficiency / sigma_T / c; 

      mdot = min(mdot, mdot_Edd); 
//    fprintf(stdout, "mdot_UpperLimit=%g, mdot_Edd=%g, mdot=%g\n", mdot_UpperLimit, mdot_Edd, mdot); 
    }

    //this->DeltaMass += mdot * (CurrentGrid->dtFixed * TimeUnits);

    this->naccretions = 1;
    if (this->accretion_rate == NULL) {
      this->accretion_rate = new float[naccretions];
      this->accretion_time = new FLOAT[naccretions];
    }
    this->accretion_rate[0] = mdot;
    this->accretion_time[0] = time;

    if (mdot > 0.0)
      fprintf(stdout, "BH Accretion[%"ISYM"]: time = %"FSYM", mdot = %"GSYM" Msun/yr, "
	      "M_BH = %"GSYM" Msun, rho = %"GSYM" g/cm3, T = %"GSYM" K, v_rel = %"GSYM" cm/s\n",
	      Identifier, time, mdot*yr, Mass, density*DensityUnits,
	      temperature[index], v_rel);

  } // ENDIF LOCAL_ACCRETION  


  
  if (AccretionType == BONDI_ACCRETION ||
      AccretionType == RADIAL_ACCRETION) {
    fprintf(stderr, "AccretionType = %"ISYM" not implemented yet.\n", 
	    AccretionType);
    ENZO_FAIL("");
  }

  delete [] temperature;

  return SUCCESS;
}

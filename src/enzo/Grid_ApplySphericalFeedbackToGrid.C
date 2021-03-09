/***********************************************************************
/
/  Algorithm for applying thermal feedback to the temporary grid of an active particle
/
/  written by: John Regan
/  date:       December, 2020
/
/  note: Based on methods originally implemented by Stephen Skory
************************************************************************/
#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "phys_constants.h"
#include "ActiveParticle_SmartStar.h"
#define MAX_TEMPERATURE 1e8

int grid::ApplySphericalFeedbackToGrid(ActiveParticleType** ThisParticle,
				       float EjectaDensity, 
				       float EjectaThermalEnergy,
				       float EjectaMetalDensity) {
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1,
    PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  double MassUnits = 1.0;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, this->ReturnTime()) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
         Vel3Num, TENum) == FAIL) {
     ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
   }
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (MultiSpecies) 
    if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				    HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				    DIINum, HDINum) == FAIL) {
        ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }
  int SNColourNum, MetalNum, Metal2Num, MBHColourNum, Galaxy1ColourNum, 
    Galaxy2ColourNum, MetalIaNum, MetalIINum;
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, Metal2Num, MetalIaNum, 
				 MetalIINum, MBHColourNum, Galaxy1ColourNum, 
				 Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  MetalNum = max(Metal2Num, SNColourNum);
  MetallicityField = (MetalNum > 0) ? TRUE : FALSE;
  ActiveParticleType_SmartStar *SS = static_cast<ActiveParticleType_SmartStar*>(* ThisParticle);
  FLOAT radius = max(64*this->CellWidth[0][0], SS->InfluenceRadius);
  //printf("%s: radius (in cellwidths) = %f\n", __FUNCTION__, radius/this->CellWidth[0][0]);
  float MetalRadius = 1.0;
  FLOAT MetalRadius2 = radius * radius * MetalRadius * MetalRadius;
  float dx = float(this->CellWidth[0][0]);
  FLOAT *pos = SS->ReturnPosition();
  FLOAT outerRadius2 = POW(1.2*radius, 2.0);
  float maxGE = MAX_TEMPERATURE / (TemperatureUnits * (Gamma-1.0) * 0.6);
  float delta_fz = 0.0;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      int index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	
	  FLOAT radius2 = POW(CellLeftEdge[0][i] + 0.5*dx - pos[0],2.0) +
	    POW(CellLeftEdge[1][j] + 0.5*dx - pos[1],2.0) +
	    POW(CellLeftEdge[2][k] + 0.5*dx - pos[2],2.0);

	  if (radius2 < outerRadius2) {
	    float r1 = sqrt(radius2) / radius;
	    float norm = 0.98;
	    float ramp = norm*(0.5 - 0.5 * tanh(10.0*(r1-1.0)));
	    /* 1/1.2^3 factor to dilute the density since we're
	       depositing a uniform ejecta in a sphere of 1.2*radius
	       without a ramp.  The ramp is only applied to the
	       energy*density factor. */
	    float factor = 0.578704;
	    
	    float OldDensity = this->BaryonField[DensNum][index];
	    BaryonField[DensNum][index] += factor*EjectaDensity;
	    /* Get specific energy */
	    if (GENum >= 0 && DualEnergyFormalism) {

	      /* When injected energy is uniform throughout the volume;
		 EjectaThermalEnergy in EnergyUnits/VolumeUnits */
	      float oldGE =  this->BaryonField[GENum][index];
	      float newGE = (OldDensity * this->BaryonField[GENum][index] +
		       ramp * factor * EjectaThermalEnergy * EjectaDensity)
		/ BaryonField[DensNum][index] ;

	      newGE = min(newGE, maxGE);  
	      //printf("%s: Energy Before = %e\t Energy injected = %e\t Increase = %e\n", __FUNCTION__, 
	      //	     oldGE,ramp * factor * EjectaThermalEnergy / Density, (newGE - oldGE)/oldGE);
	      fflush(stdout);
	      
	      this->BaryonField[GENum][index] = newGE;
	      this->BaryonField[TENum][index] = newGE;

	      for (int dim = 0; dim < GridRank; dim++)
		this->BaryonField[TENum][index] += 
		  0.5 * this->BaryonField[Vel1Num+dim][index] * 
		  this->BaryonField[Vel1Num+dim][index];

	      //printf("%s: Increase in GE energy is %e\n", __FUNCTION__, (newGE - oldGE)/oldGE);
	      
	    } else {

	      float newGE = (OldDensity * this->BaryonField[TENum][index] +
		       ramp * factor * EjectaDensity * EjectaThermalEnergy) /
		BaryonField[DensNum][index];

	      newGE = min(newGE, maxGE);  
	      this->BaryonField[TENum][index] = newGE;

	    } //end if(GENum >= 0 && DualEnergyFormalism)

	    /* Update species and colour fields */
	    if (MetallicityField == TRUE && radius2 <= MetalRadius2)
	      delta_fz = EjectaMetalDensity / OldDensity;
	    else
	      delta_fz = 0.0;
	    float increase = BaryonField[DensNum][index] / OldDensity - delta_fz;

	    if (MultiSpecies) {
	      BaryonField[DeNum][index] *= increase;
	      BaryonField[HINum][index] *= increase;
	      BaryonField[HIINum][index] *= increase;
	      BaryonField[HeINum][index] *= increase;
	      BaryonField[HeIINum][index] *= increase;
	      BaryonField[HeIIINum][index] *= increase;
	    }
	    if (MultiSpecies > 1) {
	      BaryonField[HMNum][index] *= increase;
	      BaryonField[H2INum][index] *= increase;
	      BaryonField[H2IINum][index] *= increase;
	    }
	    if (MultiSpecies > 2) {
	      BaryonField[DINum][index] *= increase;
	      BaryonField[DIINum][index] *= increase;
	      BaryonField[HDINum][index] *= increase;
	    }

	    if (MetallicityField == TRUE)
	      BaryonField[MetalNum][index] += EjectaMetalDensity;

	    /* MBHColour injected */
	    if (MBHColourNum > 0)
	      BaryonField[MBHColourNum][index] += factor*EjectaDensity;

	  } // END if inside radius
	}  // END i-direction
      }  // END j-direction
    }  // END k-direction
  return SUCCESS;
}

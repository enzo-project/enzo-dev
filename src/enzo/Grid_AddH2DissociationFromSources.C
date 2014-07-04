/***********************************************************************
/
/  ADD H2 DISSOCIATION EMISSION FROM SHINING PARTICLES
/
/  written by: John Wise
/  date:       March, 2006
/  modified1:
/
/  PURPOSE:
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
#include "CosmologyParameters.h"
#include "Star.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::AddH2DissociationFromSources(Star *AllStars)
{

  Star *cstar;
  FLOAT DomainWidth[MAX_DIMENSION];
  FLOAT *ddr2[MAX_DIMENSION];
  FLOAT innerFront, outerFront, innerFront2, outerFront2;
  double Luminosity[MAX_ENERGY_BINS];
  float energies[MAX_ENERGY_BINS], kdiss_r2;
  int ipart, dim, i, j, k, index, indixe, nbins;
  int ActiveDims[MAX_DIMENSION];
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  const double pc = 3.086e18, clight = 3e10;
  double H2Luminosity, H2ISigma = 3.71e-18;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  this->DebugCheck((char*) "Grid_AddH2Dissociation");

  /* Find Multi-species fields. */

  if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				  HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				  DIINum, HDINum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }

  /* Get photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, kphHeIINum, 
				  kdissH2INum);

  /* For now, initialize H2 photo-dissociation field. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  // Now done in Grid_InitializeRadiativeTransferFields.C
//  for (i = 0; i < size; i++)
//    BaryonField[kdissH2INum][i] = 0;
  
  /* Exit if no star particles and not Photon Test */

  if (AllStars == NULL && ProblemType != 50)
    return SUCCESS;

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, 
    TimeUnits, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Absorb the unit conversions into the cross-section
  H2ISigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);

  for (dim = 0; dim < GridRank; dim++) {
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    ddr2[dim] = new FLOAT[ActiveDims[dim]];
  }

  // Dilution factor (prevent breaking in the rate solver near the star)
  float dilutionRadius = 4.848e-6 * pc / (double) LengthUnits;  // 1 AU
  float dilRadius2 = dilutionRadius * dilutionRadius;
  float LightTravelDist = TimeUnits * clight / LengthUnits;

  // Convert from #/s to RT units
  double LConv = (double) TimeUnits / pow(LengthUnits,3);

  /* Loop over radiation sources or star particles in the grid */

  if (ProblemType == 50) {

    RadiationSourceEntry *RS;
    for (RS = GlobalRadiationSources->NextSource; RS; RS = RS->NextSource) {

      if (PhotonTime < RS->CreationTime && 
	  PhotonTime > RS->CreationTime + RS->LifeTime)
	continue;

      H2Luminosity = RS->SED[3] * RS->Luminosity / LConv;

      /* Pre-calculate distances from cells to source */

      for (dim = 0; dim < GridRank; dim++)
	for (i = 0, index = GridStartIndex[dim]; i < ActiveDims[dim]; 
	     i++, index++) {

	  // Calculate dr_i first, then square it
	  ddr2[dim][i] = 
	    fabs(CellLeftEdge[dim][index] + 0.5*CellWidth[dim][index] - 
		 RS->Position[dim]);
	  ddr2[dim][i] = min(ddr2[dim][i], DomainWidth[dim]-ddr2[dim][i]);
	  ddr2[dim][i] = ddr2[dim][i] * ddr2[dim][i];
	}

      /* Loop over cells */

      double radius2, radius2_yz;

      kdiss_r2 = (float) (H2Luminosity * H2ISigma / (4.0 * M_PI));
      for (k = 0; k < ActiveDims[2]; k++) {
	for (j = 0; j < ActiveDims[1]; j++) {
	  radius2_yz = ddr2[1][j] + ddr2[2][k];
	  index = GRIDINDEX(0, j, k);
	  for (i = 0; i < ActiveDims[0]; i++, index++) {
	    radius2 = radius2_yz + ddr2[0][i];
	    //if (radius2 < outerFront2 && radius2 > innerFront2) {
	    //radius2 = max(radius2, dilRadius2);
	    if (radius2 < dilRadius2)
	      BaryonField[kdissH2INum][index] += kdiss_r2 / dilRadius2;
	    else
	      BaryonField[kdissH2INum][index] += kdiss_r2 / radius2;
	    //} // ENDIF
	  } // END: i-direction
	} // END: j-direction
      } // END: k-direction
      
    } // ENDFOR sources

  } // ENDIF ProblemType == 50

  else {

    for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

      // Skip if not 'living'
      if (!(cstar->FeedbackFlag == NO_FEEDBACK ||
	    cstar->FeedbackFlag == CONT_SUPERNOVA)) 
	continue;
      
      /* Determine H2 emission rate */

      if (cstar->ComputePhotonRates(TimeUnits, nbins, energies, Luminosity) == FAIL) {
	ENZO_FAIL("Error in ComputePhotonRates.\n");
      }
      H2Luminosity = Luminosity[3];

      /* Pre-calculate distances from cells to source */

      for (dim = 0; dim < GridRank; dim++)
	for (i = 0, index = GridStartIndex[dim]; i < ActiveDims[dim]; 
	     i++, index++) {
	  
	  // Calculate dr_i first, then square it
	  ddr2[dim][i] = 
	    fabs(CellLeftEdge[dim][index] + 0.5*CellWidth[dim][index] -
		 cstar->pos[dim]);
	  ddr2[dim][i] = min(ddr2[dim][i], DomainWidth[dim]-ddr2[dim][i]);
	  ddr2[dim][i] = ddr2[dim][i] * ddr2[dim][i];
	}

      /* Loop over cells */

      double radius2, radius2_yz;

      kdiss_r2 = (float) (H2Luminosity * H2ISigma / (4.0 * M_PI));
      for (k = 0; k < ActiveDims[2]; k++) {
	for (j = 0; j < ActiveDims[1]; j++) {
	  radius2_yz = ddr2[1][j] + ddr2[2][k];
	  index = GRIDINDEX(0, j, k);
	  for (i = 0; i < ActiveDims[0]; i++, index++) {
	    radius2 = radius2_yz + ddr2[0][i];
	    //if (radius2 < outerFront2 && radius2 > innerFront2) {
	    //radius2 = max(radius2, dilRadius2);
	    if (radius2 < dilRadius2)

	      BaryonField[kdissH2INum][index] += kdiss_r2 / dilRadius2;
	    else
	      BaryonField[kdissH2INum][index] += kdiss_r2 / radius2;
	    //} // ENDIF
	  } // END: i-direction
	} // END: j-direction
      } // END: k-direction
    } // ENDFOR stars

  } // ENDELSE ProblemType == 50

  for (dim = 0; dim < GridRank; dim++)
    delete [] ddr2[dim];

  return SUCCESS;

}

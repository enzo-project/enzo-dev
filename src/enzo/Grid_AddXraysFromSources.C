/***********************************************************************
/
/  ADD X-RAYS EMISSION FROM SHINING PARTICLES
/
/  written by: John Wise
/  date:       April, 2012
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
FLOAT FindCrossSection(int type, float energy);

int grid::AddXraysFromSources(Star *AllStars)
{

  const float EnergyThresholds[] = {13.6, 24.6, 54.4};
  const float PopulationFractions[] = {1.0, 0.25, 0.25};
  const float erg_eV = 1.602176e-12;
  const double clight = 2.99792e10;
  const double pc = 3.086e18;

  Star *cstar;
  FLOAT DomainWidth[MAX_DIMENSION];
  FLOAT *ddr2[MAX_DIMENSION];
  double Luminosity[MAX_ENERGY_BINS];
  float energies[MAX_ENERGY_BINS], kph_r2, gamma_r2, dkph, dE;
  int ipart, dim, a, i, j, k, bin, index, indixe, nbins;
  int ActiveDims[MAX_DIMENSION];
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  double XrayLuminosity, LConv, CrossSectionConv;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Exit if no star particles and not Photon Test */

  if (AllStars == NULL && ProblemType != 50)
    return SUCCESS;

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
  const int kphNum[] = {kphHINum, kphHeINum, kphHeIINum};

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, 
    TimeUnits, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Absorb the unit conversions into the cross-section
  CrossSectionConv = (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);

  // Convert from #/s to RT units
  LConv = (double) TimeUnits / pow(LengthUnits,3);

  for (dim = 0; dim < GridRank; dim++) {
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    ddr2[dim] = new FLOAT[ActiveDims[dim]];
  }

  /* Loop over radiation sources or star particles in the grid */

  FLOAT CrossSections[3]; // HI, HeI, HeII
  float xx, heat_factor, ion2_factor[3];
  float nSecondaryHII = 1, nSecondaryHeII = 1;

  ion2_factor[2] = 1.0;
  if (RadiationXRaySecondaryIon == FALSE) {
    heat_factor = 1.0;
    ion2_factor[0] = ion2_factor[1] = 1.0;
  }

  if (ProblemType == 50) {

    RadiationSourceEntry *RS;
    for (RS = GlobalRadiationSources->NextSource; RS; RS = RS->NextSource) {

      if (PhotonTime < RS->CreationTime && 
	  PhotonTime > RS->CreationTime + RS->LifeTime)
	continue;

      /* Loop over energy bins and consider X-rays to be E >= 100 eV */

      for (bin = 0; bin < RS->EnergyBins; bin++) {

	if (RS->Energy[bin] < 100) continue;

	XrayLuminosity = RS->SED[bin] * RS->Luminosity / LConv;
	for (i = 0; i < 3; i++)
	  CrossSections[i] = FindCrossSection(i, RS->Energy[bin]) * CrossSectionConv;

	if (RadiationXRaySecondaryIon == TRUE) {
	  nSecondaryHII = RS->Energy[bin] / 13.6;
	  nSecondaryHeII = RS->Energy[bin] / 24.6;
	}

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

	/* Loop over absorbers then cells */

	double radius2_yz;

	for (a = 0; a < 3; a++) {
	  kph_r2 = (float) (XrayLuminosity * CrossSections[a] / (4.0 * M_PI));
	  dE = (RS->SED[bin] - EnergyThresholds[a]);
	  for (k = 0; k < ActiveDims[2]; k++) {
	    for (j = 0; j < ActiveDims[1]; j++) {
	      radius2_yz = ddr2[1][j] + ddr2[2][k];
	      index = GRIDINDEX(0, j, k);
	      for (i = 0; i < ActiveDims[0]; i++, index++) {
		dkph = kph_r2 / (radius2_yz + ddr2[0][i]);
		if (RadiationXRaySecondaryIon == TRUE) {
		  xx = max(BaryonField[HIINum][index] / 
			   (BaryonField[HINum][index] + BaryonField[HIINum][index]),
			   1e-4);
		  heat_factor = 0.9971 * (1.0 - powf(1.0 - powf(xx, 0.2663f), 1.3163));
		  ion2_factor[0] = 0.3908 * nSecondaryHII * 
		    powf(1 - powf(xx, 0.4092f), 1.7592f);
		  ion2_factor[1] = 0.0554 * nSecondaryHeII * 
		    powf(1 - powf(xx, 0.4614f), 1.6660f);
		}
		BaryonField[kphNum[a]][index] += dkph * ion2_factor[a];
		BaryonField[gammaNum][index] += dkph * dE * heat_factor;
	      } // END: i-direction
	    } // END: j-direction
	  } // END: k-direction
	} // END: absorbers (a)

      } // ENDFOR bin
      
    } // ENDFOR sources

  } // ENDIF ProblemType == 50

  else {

#ifdef UNUSED
    for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

      // Skip if not 'living'
      if (!(cstar->FeedbackFlag == NO_FEEDBACK ||
	    cstar->FeedbackFlag == CONT_SUPERNOVA)) 
	continue;
      
      /* Determine H2 emission rate */

      if (cstar->ComputePhotonRates(nbins, energies, Luminosity) == FAIL) {
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
#else /* UNUSED */
    ENZO_FAIL("Optically thin X-rays not implemented for stars yet.");
#endif

  } // ENDELSE ProblemType == 50

  for (dim = 0; dim < GridRank; dim++)
    delete [] ddr2[dim];

  return SUCCESS;

}

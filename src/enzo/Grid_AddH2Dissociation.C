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

int grid::AddH2Dissociation(Star *AllStars)
{

  Star *cstar;
  RadiationSourceEntry *RS;
  FLOAT delx, dely, delz, DomainWidth[MAX_DIMENSION];
  FLOAT innerFront, outerFront;
  double r2, radius2, Luminosity[4];
  float energies[4], kdiss_r2;
  int ipart, dim, i, j, k, index, indixe;
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
    fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
    ENZO_FAIL("");
  }

  /* Get photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaHINum, gammaHeINum, gammaHeIINum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				      gammaHeINum, kphHeIINum, gammaHeIINum, 
				      kdissH2INum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyRadiativeTransferFields.\n");
    ENZO_FAIL("");
  }

  /* For now, initialize H2 photo-dissociation field. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (i = 0; i < size; i++)
    BaryonField[kdissH2INum][i] = 0;
  
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

  for (dim = 0; dim < GridRank; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

  // Dilution factor (prevent breaking in the rate solver near the star)
  float dilutionRadius = 10.0 * pc / (double) LengthUnits;
  float dilRadius2 = dilutionRadius * dilutionRadius;
  float LightTravelDist = TimeUnits * clight / LengthUnits;

  // Convert from #/s to RT units
  double LConv = (double) TimeUnits / pow(LengthUnits,3);

  /* Loop over radiation sources or star particles in the grid */

  if (ProblemType == 50) {

    for (RS = GlobalRadiationSources->NextSource; RS; RS = RS->NextSource) {

      if (PhotonTime < RS->CreationTime && 
	  PhotonTime > RS->CreationTime + RS->LifeTime)
	continue;

      H2Luminosity = RS->SED[3] * RS->Luminosity / LConv;

      // Determine the inner and outer radiation fronts
      innerFront = 0;
      outerFront = (PhotonTime - RS->CreationTime) * LightTravelDist;
      
      /* Loop over cells */

      index = 0;
      kdiss_r2 = (float) (H2Luminosity * H2ISigma / (4.0 * M_PI));
      for (k = 0; k < GridDimension[2]; k++) {
	delz = fabs(CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - RS->Position[2]);
	delz = min(delz, DomainWidth[2]-delz);
	for (j = 0; j < GridDimension[1]; j++) {
	  dely = fabs(CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - RS->Position[1]);
	  dely = min(dely, DomainWidth[1]-dely);
	  for (i = 0; i < GridDimension[0]; i++, index++) {
	    delx = fabs(CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - RS->Position[0]);
	    delx = min(delx, DomainWidth[0]-delx);
	    radius2 = delx*delx + dely*dely + delz*delz;
	    if (radius2 > outerFront*outerFront || radius2 < innerFront*innerFront)
	      continue;

	    radius2 = max(radius2, dilRadius2);
	    BaryonField[kdissH2INum][index] += kdiss_r2 / radius2;

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

      if (cstar->ComputePhotonRates(energies, Luminosity) == FAIL) {
	fprintf(stderr, "Error in ComputePhotonRates.\n");
	ENZO_FAIL("");
      }
      H2Luminosity = Luminosity[3];

      // Determine the inner and outer radiation fronts
      outerFront = (PhotonTime - cstar->BirthTime) * LightTravelDist;

      if (PhotonTime > (cstar->BirthTime + cstar->LifeTime))
	innerFront = (PhotonTime - cstar->BirthTime - cstar->LifeTime)
	  * LightTravelDist;
      else
	innerFront = 0;

      /* Loop over cells */

      index = 0;
      kdiss_r2 = (float) (H2Luminosity * H2ISigma / (4.0 * M_PI));
      for (k = 0; k < GridDimension[2]; k++) {
	delz = fabs(CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2]);
	delz = min(delz, DomainWidth[2]-delz);
	for (j = 0; j < GridDimension[1]; j++) {
	  dely = fabs(CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1]);
	  dely = min(dely, DomainWidth[1]-dely);
	  for (i = 0; i < GridDimension[0]; i++, index++) {
	    delx = fabs(CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0]);
	    delx = min(delx, DomainWidth[0]-delx);
	    radius2 = delx*delx + dely*dely + delz*delz;
	    if (radius2 > outerFront*outerFront || radius2 < innerFront*innerFront)
	      continue;

	    radius2 = max(radius2, dilRadius2);
	    BaryonField[kdissH2INum][index] += kdiss_r2 / radius2;

	  } // END: i-direction
	} // END: j-direction
      } // END: k-direction

    } // END: stars

  } // ENDELSE ProblemType == 50

  return SUCCESS;

}

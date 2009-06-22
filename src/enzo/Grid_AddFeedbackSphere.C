/***********************************************************************
/
/  GRID: ADD SPHERICAL STAR PARTICLE FEEDBACK TO CELLS
/
/  written by: John Wise
/  date:       September, 2005
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

#ifdef CONFIG_BFLOAT_4
#define TOLERANCE 1e-06
#endif
#ifdef CONFIG_BFLOAT_8
#define TOLERANCE 1e-12
#endif
#ifdef CONFIG_BFLOAT_16
#define TOLERANCE 1e-15
#endif

#define MAX_TEMPERATURE 1e8

int FindField(int field, int farray[], int numfields);

int grid::AddFeedbackSphere(Star *cstar, int level, float radius, float VelocityUnits, 
			    float TemperatureUnits, double EjectaDensity, 
			    double EjectaMetalDensity, double EjectaThermalEnergy, 
			    int &CellsModified)
{

  const float WhalenMaxVelocity = 35;		// km/s

  int dim, i, j, k, index;
  int sx, sy, sz;
  FLOAT delx, dely, delz, radius2, Radius, DomainWidth[MAX_DIMENSION];
  float coef, speed, maxVelocity;
  float OldDensity;
  float r1, norm, ramp, factor, newGE, increase, fh;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Check if sphere overlaps with this grid */

  for (dim = 0; dim < GridRank; dim++)
    if (cstar->pos[dim] - radius > GridRightEdge[dim] ||
	cstar->pos[dim] + radius < GridLeftEdge[dim])
      return SUCCESS;

  for (dim = 0; dim < GridRank; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

  /* Find metallicity field and set flag. */

  int ZNum, ZField;
  int MetallicityField = FALSE, MetalNum;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) 
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;

  /* Find SN colour field */

  int UseColour = FALSE, SNColourNum;
  if ((SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields)) 
      != -1)
    UseColour = TRUE;
  else
    SNColourNum = 0;

  ZNum = max(MetalNum, SNColourNum);
  ZField = max(MetallicityField, UseColour);

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				  HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				  DIINum, HDINum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
    ENZO_FAIL("");
  }

  /***********************************************************************
                                 SUPERNOVAE
  ************************************************************************/

  // Assume that the remnant is still in the free expansion stage and
  // hasn't had any radiative losses.  In this case, the ejecta will
  // be at 3/4 the radius of the shock front (see Ostriker & McKee
  // 1988 or Tenorio-Tagle 1996).

  const float MetalRadius = 0.75;
  float ionizedFraction = 0.999;  // Assume supernova is ionized
  float maxGE, MetalRadius2, PrimordialDensity, metallicity, fhz, fhez;

  // Correct for exaggerated influence radius for pair-instability supernovae
  if (cstar->FeedbackFlag == SUPERNOVA)
    radius /= 8.0;

  // Correct if the volume with 27 cells is larger than the energy bubble volume
  float BoxVolume = 27 * CellWidth[0][0] * CellWidth[0][0] * CellWidth[0][0];
  float BubbleVolume = (4.0 * M_PI / 3.0) * radius * radius * radius;
  //printf("BoxVolume = %lg, BubbleVolume = %lg\n", BoxVolume, BubbleVolume);
  if (BoxVolume > BubbleVolume) {
    //printf("Reducing ejecta density by %g\n", BubbleVolume / BoxVolume);
    EjectaDensity *= BubbleVolume / BoxVolume;
    EjectaThermalEnergy *= BubbleVolume / BoxVolume;
  }
  if (cstar->level > level) {
//    printf("Reducing ejecta density and energy by 10%% on "
//	   "level %"ISYM" to avoid crashing.\n", level);
    EjectaDensity *= 0.1;
    EjectaThermalEnergy *= 0.1;
  }

  // Correct for smaller enrichment radius
  EjectaMetalDensity *= pow(MetalRadius, -3.0);
  PrimordialDensity = EjectaDensity - EjectaMetalDensity;
  fh = CoolData.HydrogenFractionByMass;
  MetalRadius2 = radius * radius * MetalRadius * MetalRadius;

  if (cstar->FeedbackFlag == SUPERNOVA || cstar->FeedbackFlag == CONT_SUPERNOVA) {

//    printf("SN: pos = %"FSYM" %"FSYM" %"FSYM"\n", 
//	   cstar->pos[0], cstar->pos[1], cstar->pos[2]);
    maxGE = MAX_TEMPERATURE / (TemperatureUnits * (Gamma-1.0) * 0.6);

    for (k = 0; k < GridDimension[2]; k++) {

      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
      sz = sign(delz);
      delz = fabs(delz);
      delz = min(delz, DomainWidth[2]-delz);

      for (j = 0; j < GridDimension[1]; j++) {

	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
	sy = sign(dely);
	dely = fabs(dely);
	dely = min(dely, DomainWidth[1]-dely);

	index = (k*GridDimension[1] + j)*GridDimension[0];
	for (i = 0; i < GridDimension[0]; i++, index++) {

	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	  sx = sign(delx);
	  delx = fabs(delx);
	  delx = min(delx, DomainWidth[0]-delx);

	  radius2 = delx*delx + dely*dely + delz*delz;
	  if (radius2 <= 1.2*1.2*radius*radius) {

	    r1 = sqrt(radius2) / radius;
	    norm = 0.98;
	    ramp = norm*(0.5 - 0.5 * tanh(10.0*(r1-1.0)));
	    //	     ramp = min(max(1.0 - (r1 - 0.8)/0.4, 0.01), 1.0);

	    // 1/1.2^3 factor to dilute the density since we're
	    // depositing a uniform ejecta in a sphere of 1.2*radius
	    // without a ramp.  The ramp is only applied to the
	    // energy*density factor.
	    factor = 0.578704;

	    OldDensity = BaryonField[DensNum][index];
	    BaryonField[DensNum][index] += factor*EjectaDensity;

	    /* Add total energies of spheres together, then divide by
	       density to get specific energy */
	    
	    newGE = (OldDensity * BaryonField[GENum][index] +
		     ramp * factor*EjectaDensity * EjectaThermalEnergy) /
	      BaryonField[DensNum][index];
	    newGE = min(newGE, maxGE);
//	    newGE = ramp * EjectaThermalEnergy;
//	    printf("AddSN: rho = %"GSYM"=>%"GSYM", GE = %"GSYM"=>%"GSYM", drho = %"GSYM", dE = %"GSYM"\n",
//		   OldDensity, BaryonField[DensNum][index], 
//		   BaryonField[GENum][index], newGE, EjectaDensity,
//		   EjectaThermalEnergy);

	    if (GENum >= 0 && DualEnergyFormalism) {
	      BaryonField[GENum][index] = newGE;
	      BaryonField[TENum][index] = newGE;
	    } else
	      BaryonField[TENum][index] = newGE;

	    for (dim = 0; dim < GridRank; dim++)
	      BaryonField[TENum][index] += 
		0.5 * BaryonField[Vel1Num+dim][index] * 
		BaryonField[Vel1Num+dim][index];

	    //increase = BaryonField[DensNum][index] / OldDensity;
	    if (ZField == TRUE) {
	      if (radius2 <= MetalRadius2) {
		metallicity = (BaryonField[ZNum][index] + EjectaMetalDensity) /
		  BaryonField[DensNum][index];
	      } else {
		metallicity = BaryonField[ZNum][index] / BaryonField[DensNum][index];
	      }
	    } else
	      metallicity = 0;

	    fhz = fh * (1-metallicity);
	    fhez = (1-fh) * (1-metallicity);

	    if (MultiSpecies) {
	      BaryonField[DeNum][index] = 
		BaryonField[DensNum][index] * ionizedFraction;
	      BaryonField[HINum][index] = 
		BaryonField[DensNum][index] * fhz * (1-ionizedFraction);
	      BaryonField[HIINum][index] =
		BaryonField[DensNum][index] * fhz * ionizedFraction;
	      BaryonField[HeINum][index] =
		0.5*BaryonField[DensNum][index] * fhez * (1-ionizedFraction);
	      BaryonField[HeIINum][index] =
		0.5*BaryonField[DensNum][index] * fhez * (1-ionizedFraction);
	      BaryonField[HeIIINum][index] =
		BaryonField[DensNum][index] * fhez * ionizedFraction;
	    }
	    if (MultiSpecies > 1) {
	      BaryonField[HMNum][index] = tiny_number * BaryonField[DensNum][index];
	      BaryonField[H2INum][index] = 
		tiny_number * BaryonField[DensNum][index];
	      BaryonField[H2IINum][index] = 
		tiny_number * BaryonField[DensNum][index];
	    }
	    if (MultiSpecies > 2) {
	      BaryonField[DINum][index] = BaryonField[DensNum][index] * fh *
		CoolData.DeuteriumToHydrogenRatio * (1-ionizedFraction);
	      BaryonField[DIINum][index] = BaryonField[DensNum][index] * fh *
		CoolData.DeuteriumToHydrogenRatio * ionizedFraction;
	      BaryonField[HDINum][index] = 
		tiny_number * BaryonField[DensNum][index];
	    }

	    if (ZField == TRUE)
	      BaryonField[ZNum][index] = metallicity * BaryonField[DensNum][index];

	    CellsModified++;

	  } // END if inside radius
	}  // END i-direction
      }  // END j-direction
    }  // END k-direction

  }  // END Supernova

  /***********************************************************************
          STROEMGREN SPHERE FROM A POP III STAR (WHALEN ET AL. 2004)
  ************************************************************************/

  if (cstar->FeedbackFlag == STROEMGREN) {
    
    maxVelocity = 1e5*WhalenMaxVelocity / VelocityUnits;
    coef = maxVelocity / ( 0.8*0.8 + 2*0.8 );
    index = 0;

    for (k = 0; k < GridDimension[2]; k++) {

      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
      sz = sign(delz);
      delz = fabs(delz);
      delz = min(delz, DomainWidth[2]-delz);

      for (j = 0; j < GridDimension[1]; j++) {

	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
	sy = sign(dely);
	dely = fabs(dely);
	dely = min(dely, DomainWidth[1]-dely);

	for (i = 0; i < GridDimension[0]; i++, index++) {

	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	  sx = sign(delx);
	  delx = fabs(delx);
	  delx = min(delx, DomainWidth[0]-delx);

	  radius2 = delx*delx + dely*dely + delz*delz;
	  if (radius2 <= 1.2*1.2*radius*radius) {

	    float r1 = sqrt(radius2) / radius;
	    float ramp = min(max(1.0 - (r1 - 0.8)/0.4, 0.01), 1.0);
	    
	    BaryonField[DensNum][index] = EjectaDensity * ramp;
	  
	    // Estimate of velocity profile
	    speed = coef * (r1*r1 + 2*r1) * ramp;
	    Radius = sqrt(radius2);
	    
	    if (Radius > TOLERANCE) {
	      BaryonField[Vel1Num][index] = sx*speed*delx/Radius + cstar->vel[0];
	      BaryonField[Vel2Num][index] = sy*speed*dely/Radius + cstar->vel[1];
	      BaryonField[Vel3Num][index] = sz*speed*delz/Radius + cstar->vel[2];
	    } else {
	      BaryonField[Vel1Num][index] = cstar->vel[0];
	      BaryonField[Vel2Num][index] = cstar->vel[1];
	      BaryonField[Vel3Num][index] = cstar->vel[2];
	    }

	    BaryonField[TENum][index] = EjectaThermalEnergy*ramp;
	    if (GENum >= 0)
	      BaryonField[GENum][index] = EjectaThermalEnergy*ramp;
	    
	    for (dim = 0; dim < GridRank; dim++)
	      BaryonField[TENum][index] += 
		0.5 * BaryonField[Vel1Num+dim][index] *
		BaryonField[Vel1Num+dim][index];
	    
	    CellsModified++;

	  }  // END if inside radius

	}  // END i-direction
      }  // END j-direction
    }  // END k-direction
  }  // END Pop III sphere

  /***********************************************************************
          BARYON REMOVAL AFTER ACCRETION OR STAR PARTICLE CREATION
         -- Enforce a minimum temperature for cold gas accretion --
  ************************************************************************/

  float MinimumTemperature = 1.0, AdditionalEnergy, GasEnergy;

  if (cstar->FeedbackFlag == FORMATION) {

    index = 0;
    if (cstar->type == PopII)
      MinimumTemperature = 1e4;
    
    for (k = 0; k < GridDimension[2]; k++) {

      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
      sz = sign(delz);
      delz = fabs(delz);
      delz = min(delz, DomainWidth[2]-delz);

      for (j = 0; j < GridDimension[1]; j++) {

	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
	sy = sign(dely);
	dely = fabs(dely);
	dely = min(dely, DomainWidth[1]-dely);

	for (i = 0; i < GridDimension[0]; i++, index++) {

	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	  sx = sign(delx);
	  delx = fabs(delx);
	  delx = min(delx, DomainWidth[0]-delx);

	  radius2 = delx*delx + dely*dely + delz*delz;
	  if (radius2 <= radius*radius) {

	    factor = EjectaDensity / BaryonField[DensNum][index];
	    BaryonField[DensNum][index] *= factor;

	    if (MultiSpecies) {
	      BaryonField[DeNum][index] *= factor;
	      BaryonField[HINum][index] *= factor;
	      BaryonField[HIINum][index] *= factor;
	      BaryonField[HeINum][index] *= factor;
	      BaryonField[HeIINum][index] *= factor;
	      BaryonField[HeIIINum][index] *= factor;
	    }
	    if (MultiSpecies > 1) {
	      BaryonField[HMNum][index] *= factor;
	      BaryonField[H2INum][index] *= factor;
	      BaryonField[H2IINum][index] *= factor;
	    }
	    if (MultiSpecies > 2) {
	      BaryonField[DINum][index] *= factor;
	      BaryonField[DIINum][index] *= factor;
	      BaryonField[HIINum][index] *= factor;
	      BaryonField[HDINum][index] *= factor;
	    }

	    // For cold gas accretion, set a minimum temperature of
	    // 1e4 K since it has been accreted onto the star

	    if (DualEnergyFormalism) {
	      GasEnergy = BaryonField[GENum][index];
	    } else {
	      GasEnergy = BaryonField[TENum][index];
	      for (dim = 0; dim < GridRank; dim++)
		GasEnergy -= 0.5*BaryonField[Vel1Num+dim][index] * 
		  BaryonField[Vel1Num+dim][index];
	    }
	    AdditionalEnergy = 
	      MinimumTemperature / (TemperatureUnits * (Gamma-1.0) * 0.6) - 
	      GasEnergy;
	    AdditionalEnergy = max(AdditionalEnergy, 0.0);
	    BaryonField[TENum][index] += AdditionalEnergy;
	    if (DualEnergyFormalism)
	      BaryonField[GENum][index] += AdditionalEnergy;

	    CellsModified++;

	  }  // END if inside radius

	}  // END i-direction
      }  // END j-direction
    }  // END k-direction
  }  // END star birth

  if (cstar->FeedbackFlag == COLOR_FIELD) {

    int ColorField = FindField(ForbiddenRefinement, FieldType, NumberOfBaryonFields); 
    index = 0;

    for (k = 0; k < GridDimension[2]; k++) {

      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
      sz = sign(delz);
      delz = fabs(delz);
      delz = min(delz, DomainWidth[2]-delz);

      for (j = 0; j < GridDimension[1]; j++) {

	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
	sy = sign(dely);
	dely = fabs(dely);
	dely = min(dely, DomainWidth[1]-dely);

	for (i = 0; i < GridDimension[0]; i++, index++) {

	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	  sx = sign(delx);
	  delx = fabs(delx);
	  delx = min(delx, DomainWidth[0]-delx);

	  radius2 = delx*delx + dely*dely + delz*delz;
	  if (radius2 <= radius*radius) {

	    BaryonField[ColorField][index] =
            BaryonField[DensNum][index];

	    CellsModified++;

	  }  // END if inside radius

	}  // END i-direction
      }  // END j-direction
    }  // END k-direction
  }

  /* Now it's done, unmark. */

  //cstar->FeedbackFlag = NO_FEEDBACK;

  return SUCCESS;

}
